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
                              adjust_to_yield)

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

if "table" not in st.session_state:
    st.session_state.table = DEFAULT_ROWS.copy()

st.title("⚗️ Lab Formulation Calculator")

mode = st.radio("Mode", ["Molar ratios → batch", "Grams → ratios & theory"],
                horizontal=True,
                help="Forward: design a charge from ratios. Inverse: enter an "
                     "existing formula in grams to get molar ratios, wt%, and "
                     "all theoretical values — then modify it at constant yield.")
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
