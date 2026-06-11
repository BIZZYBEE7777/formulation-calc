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
    import blend_ui
    blend_ui.render()
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

# ---------------- materials library ----------------
LIB_FILE = Path(__file__).parent / "materials.csv"
if LIB_FILE.exists():
    try:
        _lib = pd.read_csv(LIB_FILE, dtype=str).fillna("")
        lib_names = _lib["Name"].tolist()
        lc1, lc2 = st.columns([3, 1])
        with lc1:
            pick_mat = st.selectbox("📚 Add from materials library",
                                    ["—"] + lib_names)
        with lc2:
            st.write("")
            if st.button("Add to table") and pick_mat != "—":
                row = _lib[_lib["Name"] == pick_mat].iloc[0]
                new_row = {
                    "Component": row.get("Name", ""),
                    "CAS / name for lookup": row.get("CAS", ""),
                    "MW (g/mol)": float(row.get("MW") or 0) or None,
                    "Assay %": float(row.get("Assay %") or 100),
                    "Functionality": float(row.get("Functionality") or 2),
                    "Group": (row.get("Group") or "acid").strip().lower(),
                    "Amount": 1.0,
                }
                st.session_state.table = pd.concat(
                    [st.session_state.table, pd.DataFrame([new_row])],
                    ignore_index=True)
                st.rerun()
    except Exception as e:
        st.warning(f"materials.csv found but unreadable: {e}")
else:
    st.caption("📚 Tip: add a `materials.csv` to the repo (columns: Name, "
               "CAS, MW, Functionality, Assay %, Group, Notes) and your "
               "lab's raw materials become a pick list here.")

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
        "none", "amidation (acid+amine)",
        "amidation + imidazoline (acid+amine)",
        "esterification (acid+OH)"])
with c2:
    extent = 1.0
    if reaction != "none":
        extent = st.slider("Conversion of limiting group (p)", 0.50, 1.00,
                           1.00, 0.005)
cyc_extent = 0.0
if "imidazoline" in reaction:
    cyc_extent = st.slider(
        "Cyclization extent (amide → imidazoline)", 0.0, 1.0, 0.0, 0.01,
        help="Fraction of eligible amide pairs ring-closed, releasing a "
             "SECOND water each (one ring max per polyamine molecule). "
             "Receiver water vs. prediction is your live measure of this.")

# ---------------- imidazoline designer ----------------
if "imidazoline" in reaction:
    with st.expander("🧬 Imidazoline designer — ratio · conversion · "
                     "cyclization → % imidazoline (and back)"):
        st.caption("For a MONOACID (TOFA-type) + DETA-type polyamine. "
                   "Sequential acylation, one ring per polyamine. The ratio "
                   "sets the ceiling and species mix; **cyclization — a "
                   "process variable (temperature, time, vacuum) — sets the "
                   "%**. Verify against receiver water and the ~1605 cm⁻¹ "
                   "C=N band.")
        from formulation_core import imidazoline_species, solve_imidazoline
        i1, i2, i3, i4 = st.columns(4)
        with i1:
            R_im = st.number_input("Acid : amine molar ratio", 0.05, 2.0,
                                   1.0, 0.05, format="%.2f")
        with i2:
            p_im = st.number_input("Amidation conversion p", 0.50, 1.0, 1.0,
                                   0.01, format="%.2f")
        with i3:
            mw_a = st.number_input("Acid MW (TOFA ≈ 285)", 50.0, 1500.0,
                                   285.0, format="%.1f")
        with i4:
            mw_n = st.number_input("Amine MW (DETA 103.2)", 50.0, 500.0,
                                   103.17, format="%.2f")

        tab_fwd, tab_inv = st.tabs(["Forward: c → % imid",
                                    "Inverse: % imid → required c"])
        with tab_fwd:
            c_im = st.slider("Cyclization extent c", 0.0, 1.0, 0.70, 0.01)
            res = imidazoline_species(R_im, p_im, c_im, mw_a, mw_n)
            f1, f2, f3 = st.columns(3)
            f1.metric("% imid of acylated product",
                      f"{res['pct_imid_of_acylated']:.1f}%",
                      help="imidazoline ÷ (imidazoline + amide species)")
            f2.metric("Rings per polyamine",
                      f"{res['pct_imid_per_amine']:.1f}%")
            f3.metric("wt% ring species in product",
                      f"{res['wt_pct_imid_species']:.1f}%")
            sp = res["species_per_amine"]
            st.dataframe(pd.DataFrame(
                {"Species": list(sp.keys()),
                 "mol per mol amine": [round(v, 4) for v in sp.values()]}),
                hide_index=True, use_container_width=True)
            if res["free_acid_per_amine"] > 1e-6:
                st.caption(f"Unreacted acid riding in product: "
                           f"{res['free_acid_per_amine']:.3f} mol/mol amine.")
        with tab_inv:
            v1, v2, v3 = st.columns(3)
            with v1:
                tgt_im = st.number_input("Target % imidazoline", 1.0, 100.0,
                                         70.0, format="%.1f")
            with v2:
                basis = st.selectbox("Basis", ["acylated", "amine", "weight"],
                                     format_func={
                                         "acylated": "of acylated product",
                                         "amine": "rings per polyamine",
                                         "weight": "wt% ring species"}.get)
            with v3:
                st.write("")
                if st.button("Solve required cyclization"):
                    c_req, res2, warn_im = solve_imidazoline(
                        tgt_im, basis, R_im, p_im, mw_a, mw_n)
                    if warn_im:
                        st.warning(warn_im)
                    else:
                        st.success(f"**Required cyclization c = "
                                   f"{c_req:.3f}** at R = {R_im:g}, "
                                   f"p = {p_im:g}")
                        st.info(f"At that c: {res2['pct_imid_of_acylated']:.1f}% "
                                f"of acylated · {res2['pct_imid_per_amine']:.1f}% "
                                f"per amine · {res2['wt_pct_imid_species']:.1f} "
                                f"wt% — set the main cyclization slider to "
                                f"{c_req:.2f} to see water/receiver "
                                f"predictions for the full charge.")

run = st.button("Calculate batch", type="primary", use_container_width=True)

if run and len(valid):
    comps = [{
        "name": r["Component"], "cas": r["CAS / name for lookup"],
        "mw": float(r["MW (g/mol)"]), "assay": float(r["Assay %"]) / 100.0,
        "functionality": float(r["Functionality"]),
        "group": r["Group"],
    } for _, r in valid.iterrows()]
    rxn = {"none": "none", "amidation (acid+amine)": "amidation",
           "amidation + imidazoline (acid+amine)": "amidation",
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
        out = batch_summary(comps, scale, reaction=rxn, extent=extent,
                            cyc_extent=cyc_extent)
    except Exception as e:
        st.error(f"Calculation failed: {e}")
        st.stop()

    st.session_state.last = {"out": out, "rxn": rxn, "extent": extent,
                             "cyc": cyc_extent, "grams_mode": grams_mode,
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
        cyc_note = (f" (amide {out['cond']['water_amide_g']:.1f} + "
                    f"ring {out['cond']['water_cyc_g']:.1f})"
                    if out['cond'].get('water_cyc_g', 0) > 0 else "")
        m2.metric("Predicted condensate",
                  f"{out['cond']['water_g']:.2f} g H₂O", cyc_note or None)
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

    # ---------------- cook tracker ----------------
    if rxn != "none":
        with st.expander("🧪 Cook tracker — where am I in the reaction?"):
            st.caption("Enter a mid-cook titration; get implied conversion, "
                       "predicted water to this point, and projected finals.")
            k1, k2, k3 = st.columns(3)
            with k1:
                mkey = st.selectbox("Measured value", [
                    "acid_value", "amine_value", "hydroxyl_value"],
                    format_func=lambda s: s.replace("_", " ").title())
            with k2:
                mval = st.number_input("Titrated value (mg KOH/g)",
                                       min_value=0.0, value=40.0,
                                       format="%.1f")
            with k3:
                st.write("")
                go_track = st.button("Locate cook")
            if go_track:
                from formulation_core import p_from_measured
                p_now, summ, twarn = p_from_measured(
                    L["comps"], rxn, mkey, mval, cyc_extent=L.get("cyc", 0.0))
                if twarn:
                    st.warning(twarn)
                else:
                    frac = out["cond"]["water_g"]
                    w_now = summ["cond"]["water_g"] / summ["charge_mass"]                         * out["charge_mass"]
                    t1c, t2c, t3c = st.columns(3)
                    t1c.metric("Implied conversion p", f"{p_now:.3f}")
                    t2c.metric("Water evolved by now",
                               f"{w_now:.1f} g",
                               help="Compare to the receiver. Receiver HIGH "
                                    "vs this → more cyclization than the "
                                    "slider assumes; LOW → less.")
                    t3c.metric("Water remaining to p=1",
                               f"{max(frac - w_now, 0):.1f} g")
                    fin = out["end_values"]
                    st.info(f"Projected finals at p = {L['extent']:.3f}: "
                            f"AV {fin.get('acid_value', 0):.1f} · AmV "
                            f"{fin.get('amine_value', 0):.1f} · OHV "
                            f"{fin.get('hydroxyl_value', 0):.1f}")

    # ---------------- spec targeting (ratio mode) ----------------
    if rxn != "none" and not L["grams_mode"]:
        with st.expander("🎯 Hit a target spec — solve a component ratio"):
            st.caption("Pick the component to vary; everything else holds. "
                       "Solved at the current p and cyclization settings.")
            comp_names_t = [c["name"] for c in L["comps"]]
            g1, g2, g3, g4 = st.columns(4)
            with g1:
                vary = st.selectbox("Vary component", comp_names_t)
            with g2:
                tkey = st.selectbox("Target", [
                    "acid_value", "amine_value", "hydroxyl_value"],
                    format_func=lambda s: s.replace("_", " ").title(),
                    key="tkey")
            with g3:
                tval = st.number_input("Target value", min_value=0.0,
                                       value=280.0, format="%.1f")
            with g4:
                st.write("")
                go_solve = st.button("Solve ratio")
            if go_solve:
                from formulation_core import solve_ratio_for_target
                r_new, summ2, swarn = solve_ratio_for_target(
                    L["comps"], vary, rxn, tkey, tval,
                    extent=L["extent"], cyc_extent=L.get("cyc", 0.0))
                if swarn:
                    st.warning(swarn)
                else:
                    old_r = next(c["ratio"] for c in L["comps"]
                                 if c["name"] == vary)
                    st.success(f"**{vary}: molar ratio {old_r:.4f} → "
                               f"{r_new:.4f}** hits "
                               f"{tkey.replace('_', ' ')} = {tval:.1f}")
                    e2 = summ2["end_values"]
                    st.info(f"Resulting theory: AV {e2['acid_value']:.1f} · "
                            f"AmV {e2['amine_value']:.1f} · OHV "
                            f"{e2['hydroxyl_value']:.1f}")
                    st.session_state.changelog.append(
                        f"[{date.today()}] Spec solve: {vary} ratio "
                        f"{old_r:.4f} → {r_new:.4f} for "
                        f"{tkey.replace('_', ' ')} {tval:.1f}")
                    tbl = st.session_state.table.copy()
                    for i2 in tbl.index:
                        if str(tbl.at[i2, "Component"]).strip() == vary:
                            tbl.at[i2, "Amount"] = r_new
                    st.session_state.table = tbl
                    st.caption("Ratio written to the table — press "
                               "**Calculate batch** to refresh all numbers.")

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
