"""
Lab Formulation Calculator -- Streamlit app (v2).
Modes:
  - Molar ratios -> batch   (forward: design a charge)
  - Grams -> ratios & theory (inverse: analyze/modify an existing formula)
Plus: wt%% columns, and adjust-to-yield reformulation in grams mode.

Local test: pip install streamlit pandas requests
            streamlit run formulation_app.py
"""
import json
import time
from datetime import date
from pathlib import Path
from urllib.parse import quote

import pandas as pd
import streamlit as st

from formulation_core import (solve_scale, batch_summary, charge_sheet_text,
                              grams_to_ratios, weight_percents,
                              adjust_to_yield, cold_blend_solids,
                              dilute_to_solids)

st.set_page_config(page_title="Formulation Calculator", page_icon="⚗️",
                   layout="wide")

CACHE_FILE = Path(__file__).parent / "mw_cache.json"
PUBCHEM_PROP = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
                "{}/property/MolecularWeight,IUPACName,MolecularFormula/JSON")

GROUP_HELP = ("'anhydride' = f×acid eq, first ring-opening releases no H2O "
              "(MA: f=2). 'capper' = addition capper, consumes an acid eq, "
              "no H2O (DCPD: f=1). 'inert' = solvent/non-reactive.")


def load_mw_cache():
    if CACHE_FILE.exists():
        try:
            return json.loads(CACHE_FILE.read_text())
        except json.JSONDecodeError:
            return {}
    return {}


def fetch_mw(identifier: str):
    cache = load_mw_cache()
    key = identifier.lower().strip()
    if key in cache:
        c = cache[key]
        return c["mw"], c["iupac"], c["formula"]
    import requests
    try:
        r = requests.get(PUBCHEM_PROP.format(quote(identifier)), timeout=15)
        if r.status_code == 200:
            p = r.json()["PropertyTable"]["Properties"][0]
            mw = float(p["MolecularWeight"])
            cache[key] = {"mw": mw, "iupac": p.get("IUPACName", ""),
                          "formula": p.get("MolecularFormula", "")}
            CACHE_FILE.write_text(json.dumps(cache, indent=1))
            time.sleep(0.22)
            return mw, cache[key]["iupac"], cache[key]["formula"]
    except Exception:
        pass
    return None, None, None


DEFAULT_ROWS = pd.DataFrame([
    {"Component": "", "CAS / name for lookup": "", "MW (g/mol)": None,
     "Assay %": 100.0, "Functionality": 2.0, "Group": "acid",
     "Amount": 1.0},
])

DEFAULT_BLEND = pd.DataFrame([
    {"Ingredient": "Eastek 1400", "% solids": 30.0, "Parts (by weight)": 50.0},
    {"Ingredient": "Water", "% solids": 0.0, "Parts (by weight)": 50.0},
])

if "table" not in st.session_state:
    st.session_state.table = DEFAULT_ROWS.copy()
if "blend_table" not in st.session_state:
    st.session_state.blend_table = DEFAULT_BLEND.copy()
if "changelog" not in st.session_state:
    st.session_state.changelog = []

st.title("⚗️ Lab Formulation Calculator")

# ---------------- save / load formulas ----------------
with st.expander("💾 Save / load formula"):
    sc1, sc2 = st.columns(2)
    with sc1:
        fname = st.text_input("Formula name", "formula")
        payload = json.dumps({
            "name": fname,
            "saved": str(date.today()),
            "table": st.session_state.table.to_dict(orient="records"),
            "blend_table": st.session_state.blend_table.to_dict(
                orient="records"),
            "changelog": st.session_state.changelog,
        }, indent=1)
        st.download_button("⬇️ Save current formula (.json)", payload,
                           file_name=f"{fname}.json", mime="application/json")
    with sc2:
        up = st.file_uploader("Load a saved formula", type=["json"])
        if up is not None and st.button("Load it"):
            try:
                data = json.loads(up.read().decode())
                st.session_state.table = pd.DataFrame(data["table"])
                if "blend_table" in data:
                    st.session_state.blend_table = pd.DataFrame(
                        data["blend_table"])
                st.session_state.changelog = data.get("changelog", [])
                st.success(f"Loaded '{data.get('name','formula')}' "
                           f"(saved {data.get('saved','?')}). ")
                st.rerun()
            except Exception as e:
                st.error(f"Couldn't load file: {e}")

mode = st.radio("Mode", ["Molar ratios → batch", "Grams → ratios & theory",
                         "Cold blend (wt% + solids)"],
                horizontal=True,
                help="Forward: design a charge from ratios. Inverse: analyze "
                     "an existing formula in grams. Cold blend: combine "
                     "dispersions/solvents by weight parts with per-"
                     "ingredient % solids — no MW or chemistry needed.")

# ===================== COLD BLEND MODE =====================
if mode.startswith("Cold blend"):
    from formulation_core import (blend_summary, blend_solve_diluent,
                                  blend_batch_grams, blend_sheet_text,
                                  US_FLOZ_ML)
    st.subheader("1 · Ingredients (weight parts)")
    st.caption("Parts are relative by weight — 50/50, 2:1:1, real grams, "
               "wt% — any consistent basis works. % solids: dispersions per "
               "their TDS (Eastek 1400 = 30), solvents/water = 0.")
    bt = st.data_editor(
        st.session_state.blend_table, num_rows="dynamic",
        use_container_width=True,
        column_config={
            "Ingredient": st.column_config.TextColumn(required=True),
            "% solids": st.column_config.NumberColumn(min_value=0.0,
                                                      max_value=100.0,
                                                      format="%.1f"),
            "Parts (by weight)": st.column_config.NumberColumn(
                min_value=0.0, format="%.3f"),
        }, key="blend_editor")
    st.session_state.blend_table = bt

    bvalid = bt[(bt["Ingredient"].astype(str).str.strip() != "") &
                (bt["Parts (by weight)"] > 0)]
    if len(bvalid):
        names = bvalid["Ingredient"].tolist()
        solids = [float(x) for x in bvalid["% solids"]]
        parts = [float(x) for x in bvalid["Parts (by weight)"]]
        wt, blend_pct, _ = blend_summary(parts, solids)

        st.subheader("2 · Blend")
        bdf = pd.DataFrame({
            "Ingredient": names, "% solids": solids,
            "Parts": [round(p, 3) for p in parts],
            "wt%": [round(100 * w, 2) for w in wt],
            "Solids contribution (per 100 g)":
                [round(w * s, 2) for w, s in zip(wt, solids)],
        })
        st.dataframe(bdf, use_container_width=True, hide_index=True)
        st.metric("Blend solids", f"{blend_pct:.2f}%")

        # ---- solve diluent to target solids ----
        st.subheader("3 · Hit a target % solids")
        t1, t2, t3 = st.columns(3)
        with t1:
            dil = st.selectbox("Adjust which ingredient", names)
        with t2:
            tgt = st.number_input("Target % solids", min_value=0.01,
                                  max_value=99.9, value=10.0, format="%.2f")
        di = names.index(dil)
        with t3:
            st.write("")
            if st.button("Solve & apply"):
                new_parts, warn = blend_solve_diluent(parts, solids, di, tgt)
                if warn:
                    st.warning(warn)
                else:
                    st.session_state.changelog.append(
                        f"[{date.today()}] Blend solids {blend_pct:.2f}% → "
                        f"{tgt:.2f}%: {dil} {parts[di]:.3f} → "
                        f"{new_parts[di]:.3f} parts")
                    tbl = st.session_state.blend_table.copy()
                    for i2 in tbl.index:
                        nm = str(tbl.at[i2, "Ingredient"]).strip()
                        if nm in names:
                            tbl.at[i2, "Parts (by weight)"] = \
                                new_parts[names.index(nm)]
                    st.session_state.blend_table = tbl
                    st.rerun()

        # ---- batch scaling ----
        st.subheader("4 · Scale to batch size")
        s1, s2 = st.columns(2)
        with s1:
            density = st.number_input(
                "Blend density (g/mL) — for volume sizes", min_value=0.5,
                max_value=2.5, value=1.00, format="%.3f",
                help="Waterborne blends are usually 1.00–1.05; check or "
                     "measure for accuracy. Mass-based sizes ignore this.")
        with s2:
            custom_g = st.number_input("Custom batch (g, 0 = off)",
                                       min_value=0.0, value=0.0,
                                       format="%.1f")
        sizes = [("4 oz", 4 * US_FLOZ_ML * density),
                 ("Pint", 16 * US_FLOZ_ML * density),
                 ("Quart", 32 * US_FLOZ_ML * density),
                 ("Gallon", 128 * US_FLOZ_ML * density),
                 ("100 g", 100.0), ("1 kg", 1000.0)]
        if custom_g > 0:
            sizes.append((f"{custom_g:.0f} g", custom_g))
        scale_df = pd.DataFrame({"Ingredient": names, "% solids": solids})
        for label, grams_total in sizes:
            g = blend_batch_grams(parts, grams_total)
            scale_df[f"{label} (g)"] = [round(x, 2) for x in g]
        st.dataframe(scale_df, use_container_width=True, hide_index=True)

        # ---- charge sheet ----
        st.subheader("5 · Charge sheet")
        pick = st.selectbox("Sheet for batch size",
                            [lab for lab, _ in sizes])
        gtot = dict(sizes)[pick]
        h1, h2, h3 = st.columns(3)
        meta = {"project": h1.text_input("Project code", "", key="bp"),
                "exp_id": h2.text_input("Batch ID", "", key="be"),
                "chemist": h3.text_input("Chemist", "", key="bc"),
                "date": str(date.today())}
        sheet = blend_sheet_text(names, solids,
                                 blend_batch_grams(parts, gtot),
                                 pick, blend_pct, meta)
        st.code(sheet, language=None)
        st.download_button("Download charge sheet (.txt)", sheet,
                           file_name=f"blend_{meta['exp_id'] or 'batch'}.txt")

        if st.session_state.changelog:
            with st.expander(
                    f"📜 Change log ({len(st.session_state.changelog)})"):
                for line in reversed(st.session_state.changelog):
                    st.text(line)

    st.divider()
    st.caption("Blend solids = Σ(wt% × ingredient solids). Volume sizes use "
               "the density you enter — measure it for tight work. Planning "
               "tool, not a CoA.")
    st.stop()

grams_mode = mode.startswith("Grams")
amount_label = "As-is grams" if grams_mode else "Molar ratio"

st.subheader("1 · Components")
st.caption(f"**Amount column = {amount_label}** in this mode. MW auto-fills "
           "from PubChem; enter MW manually for UVCBs (effective MW).")

edited = st.data_editor(
    st.session_state.table,
    num_rows="dynamic",
    use_container_width=True,
    column_config={
        "Component": st.column_config.TextColumn(required=True),
        "CAS / name for lookup": st.column_config.TextColumn(),
        "MW (g/mol)": st.column_config.NumberColumn(min_value=0.0,
                                                    format="%.2f"),
        "Assay %": st.column_config.NumberColumn(min_value=0.1,
                                                 max_value=100.0,
                                                 format="%.1f"),
        "Functionality": st.column_config.NumberColumn(min_value=0.0,
                                                       format="%.2f"),
        "Group": st.column_config.SelectboxColumn(
            options=["acid", "amine", "hydroxyl", "anhydride", "capper",
                     "inert"], help=GROUP_HELP),
        "Amount": st.column_config.NumberColumn(
            label=amount_label, min_value=0.0, format="%.4f"),
    },
    key="editor",
)

if st.button("🔎 Fill missing MWs from PubChem"):
    misses = []
    for i, row in edited.iterrows():
        ident = (row["CAS / name for lookup"] or row["Component"] or "").strip()
        if ident and (pd.isna(row["MW (g/mol)"]) or not row["MW (g/mol)"]):
            mw, _, formula = fetch_mw(ident)
            if mw:
                edited.at[i, "MW (g/mol)"] = mw
                st.toast(f"{ident}: {mw:.2f} g/mol ({formula})")
            else:
                misses.append(ident)
    st.session_state.table = edited
    if misses:
        st.warning(f"No PubChem match for: {', '.join(misses)} — enter MW "
                   "manually (typical for UVCBs).")
    st.rerun()

valid = edited.dropna(subset=["MW (g/mol)"])
valid = valid[(valid["Component"].astype(str).str.strip() != "") &
              (valid["Amount"] > 0)]
names = valid["Component"].tolist()

# ---------------- controls ----------------
st.subheader("2 · " + ("Reaction" if grams_mode else "Scale anchor & reaction"))
if grams_mode:
    c1, c2 = st.columns(2)
    anchor_type = anchor_value = anchor_comp = None
else:
    c1, c2, c3, c4 = st.columns(4)
    with c3:
        anchor_type = st.selectbox("Scale by", [
            "Total charge mass (g)", "Fixed component mass (g)",
            "Fixed component moles"])
    with c4:
        anchor_value = st.number_input("Target value", min_value=0.0,
                                       value=500.0, format="%.3f")
    anchor_comp = None
    if anchor_type != "Total charge mass (g)" and names:
        anchor_comp = st.selectbox("Anchor component", names)
with c1:
    reaction = st.selectbox("Condensation reaction", [
        "none", "amidation (acid+amine)", "esterification (acid+OH)"])
with c2:
    extent = 1.0
    if reaction != "none":
        extent = st.slider("Conversion of limiting group (p)", 0.50, 1.00,
                           1.00, 0.005)

run = st.button("Calculate batch", type="primary", use_container_width=True)

if run and len(valid):
    comps = [{
        "name": r["Component"], "cas": r["CAS / name for lookup"],
        "mw": float(r["MW (g/mol)"]), "assay": float(r["Assay %"]) / 100.0,
        "functionality": float(r["Functionality"]),
        "group": r["Group"],
    } for _, r in valid.iterrows()]
    rxn = {"none": "none", "amidation (acid+amine)": "amidation",
           "esterification (acid+OH)": "esterification"}[reaction]

    try:
        if grams_mode:
            for c, (_, r) in zip(comps, valid.iterrows()):
                c["g_asis_in"] = float(r["Amount"])
            comps, norm_ratios = grams_to_ratios(comps)
            scale = 1.0
        else:
            for c, (_, r) in zip(comps, valid.iterrows()):
                c["ratio"] = float(r["Amount"])
            atype = {"Total charge mass (g)": "total_asis_mass",
                     "Fixed component mass (g)": "component_asis_mass",
                     "Fixed component moles": "component_moles"}[anchor_type]
            scale = solve_scale(comps, atype, anchor_value, anchor_comp)
            norm_ratios = None
        out = batch_summary(comps, scale, reaction=rxn, extent=extent)
    except Exception as e:
        st.error(f"Calculation failed: {e}")
        st.stop()

    st.session_state.last = {"out": out, "rxn": rxn, "extent": extent,
                             "grams_mode": grams_mode,
                             "norm_ratios": norm_ratios,
                             "comps": comps}

# ---------------- results (persist across reruns for the adjust panel) ----
if "last" in st.session_state:
    L = st.session_state.last
    out, rxn = L["out"], L["rxn"]
    wp = weight_percents(out["rows"])

    st.subheader("3 · Charge table")
    rows_disp = []
    for i, r in enumerate(out["rows"]):
        d = {"Component": r["name"], "CAS": r.get("cas", ""),
             "Moles": round(r["moles"], 4), "Eq": round(r["eq"], 4),
             "Active g": round(r["g_active"], 2),
             "As-is g": round(r["g_asis"], 2),
             "wt% (as-is)": round(wp[i][0], 2),
             "wt% (active)": round(wp[i][1], 2)}
        if L.get("norm_ratios"):
            d["Molar ratio"] = round(L["norm_ratios"][i], 4)
        rows_disp.append(d)
    st.dataframe(pd.DataFrame(rows_disp), use_container_width=True,
                 hide_index=True)

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Total charge (as-is)", f"{out['charge_mass']:.1f} g")
    if rxn != "none":
        m2.metric("Predicted condensate", f"{out['cond']['water_g']:.2f} g H₂O")
        m3.metric("Theoretical final mass", f"{out['final_mass']:.1f} g")
        m4.metric("Theoretical solids", f"{out['solids_pct']:.1f}%")
    else:
        cb_s, cb_t, cb_pct = cold_blend_solids(out["rows"])
        m2.metric("Cold blend solids", f"{cb_pct:.1f}%",
                  help="Assay-corrected active mass of all non-inert "
                       "components ÷ total as-is mass")
        m3.metric("Solids mass", f"{cb_s:.1f} g")

    # ---------------- dilute to % solids (cold blend, grams mode) ----------
    if L["grams_mode"] and rxn == "none":
        inert_names = [r["name"] for r in out["rows"] if r["group"] == "inert"]
        st.subheader("Dilute / concentrate to target % solids")
        if not inert_names:
            st.caption("Mark a diluent component as 'inert' to enable this.")
        else:
            d1, d2, d3 = st.columns(3)
            with d1:
                diluent = st.selectbox("Diluent (inert)", inert_names)
            with d2:
                target_pct = st.number_input("Target % solids", min_value=0.1,
                                             max_value=100.0,
                                             value=round(cb_pct, 1),
                                             format="%.1f")
            comp_names_all = [r["name"] for r in out["rows"]]
            di = comp_names_all.index(diluent)
            with d3:
                st.write("")
                if st.button("Apply dilution"):
                    amounts = [r["g_asis"] for r in out["rows"]]
                    new_amts, warn = dilute_to_solids(amounts, out["rows"],
                                                      di, target_pct)
                    if warn:
                        st.warning(warn)
                    else:
                        st.session_state.changelog.append(
                            f"[{date.today()}] Solids {cb_pct:.1f}% → "
                            f"{target_pct:.1f}%: {diluent} "
                            f"{amounts[di]:.2f} → {new_amts[di]:.2f} g")
                        tbl = st.session_state.table.copy()
                        for i2 in tbl.index:
                            rname = str(tbl.at[i2, "Component"]).strip()
                            if rname in comp_names_all:
                                tbl.at[i2, "Amount"] = \
                                    new_amts[comp_names_all.index(rname)]
                        st.session_state.table = tbl
                        st.success("Diluent updated — press **Calculate "
                                   "batch** to recompute.")
                        st.rerun()

    ev, car = out["end_values"], out["carothers"]
    if rxn != "none":
        st.subheader("4 · Stoichiometric feedback")
        f1, f2, f3 = st.columns(3)
        if ev:
            f1.metric("Theoretical acid value", f"{ev['acid_value']:.1f}")
            f2.metric("Theoretical amine value", f"{ev['amine_value']:.1f}")
            f3.metric("Theoretical OH value", f"{ev['hydroxyl_value']:.1f}")
        if car:
            xn_s = "∞" if car["Xn"] == float("inf") else f"{car['Xn']:.1f}"
            mn_s = "∞" if car["Mn"] == float("inf") else f"{car['Mn']:.0f} g/mol"
            st.info(f"r = {car['r']:.3f} ({car['excess_pct']:.1f}% excess "
                    f"{car['excess_group']}) · at p = {L['extent']:.3f}: "
                    f"Xn ≈ {xn_s}, Mn ≈ {mn_s}")
            if car["r"] > 0.995 and L["extent"] > 0.98:
                st.warning("⚠️ Near-perfect stoichiometry at high conversion "
                           "— Mn unbounded (gel risk for thermoplastic targets).")

    # ---------------- adjust to yield (grams mode) ----------------
    if L["grams_mode"]:
        st.subheader("5 · Modify at constant yield")
        st.caption("Change one ingredient; total batch mass stays fixed. "
                   "Choose what absorbs the difference.")
        comp_names = [r["name"] for r in out["rows"]]
        a1, a2, a3 = st.columns(3)
        with a1:
            target = st.selectbox("Ingredient to change", comp_names)
        ti = comp_names.index(target)
        with a2:
            new_g = st.number_input(
                "New as-is grams", min_value=0.0,
                value=float(out["rows"][ti]["g_asis"]), format="%.2f")
        with a3:
            how = st.radio("Absorb difference in",
                           ["All other ingredients (pro-rata)",
                            "Selected ingredients only"])
        if how.startswith("Selected"):
            absorbers = st.multiselect(
                "Absorber ingredients",
                [n for n in comp_names if n != target])
            absorber_idx = [comp_names.index(n) for n in absorbers]
        else:
            absorber_idx = [i for i in range(len(comp_names)) if i != ti]

        if st.button("Apply adjustment & recalculate"):
            amounts = [r["g_asis"] for r in out["rows"]]
            new_amts, warn = adjust_to_yield(amounts, ti, new_g, absorber_idx)
            if warn:
                st.warning(warn)
            absorber_names = [comp_names[i] for i in absorber_idx]
            av_note = ""
            if out["end_values"]:
                av_note = (f"; prior AV {out['end_values']['acid_value']:.1f}"
                           f", OHV {out['end_values']['hydroxyl_value']:.1f}")
            st.session_state.changelog.append(
                f"[{date.today()}] {target}: "
                f"{amounts[ti]:.2f} → {new_g:.2f} g; absorbed by "
                f"{', '.join(absorber_names)}{av_note}")
            tbl = st.session_state.table.copy()
            j = 0
            for i2 in tbl.index:
                rname = str(tbl.at[i2, "Component"]).strip()
                if rname in comp_names:
                    tbl.at[i2, "Amount"] = new_amts[comp_names.index(rname)]
            st.session_state.table = tbl
            st.success("Amounts updated in the table — press "
                       "**Calculate batch** to recompute everything.")
            st.rerun()

        if st.session_state.changelog:
            with st.expander(
                    f"📜 Change log ({len(st.session_state.changelog)})"):
                for line in reversed(st.session_state.changelog):
                    st.text(line)
                if st.button("Clear change log"):
                    st.session_state.changelog = []
                    st.rerun()

    # ---------------- charge sheet ----------------
    st.subheader("6 · Charge sheet" if L["grams_mode"] else "5 · Charge sheet")
    h1, h2, h3 = st.columns(3)
    meta = {"project": h1.text_input("Project code", ""),
            "exp_id": h2.text_input("Experiment ID", ""),
            "chemist": h3.text_input("Chemist", ""),
            "date": str(date.today()),
            "description": st.text_input("Batch description", "")}
    sheet = charge_sheet_text(out, meta)
    st.code(sheet, language=None)
    st.download_button("Download charge sheet (.txt)", sheet,
                       file_name=f"charge_{meta['exp_id'] or 'batch'}.txt")

st.divider()
st.caption("Verify MWs and assays against CoAs. Theoretical values assume "
           "ideal step growth; water prediction credits cappers and anhydride "
           "first-openings as water-free. Planning tool, not a CoA.")
