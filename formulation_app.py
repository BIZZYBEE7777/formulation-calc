"""
Lab Formulation Calculator -- Streamlit app.
Deploys identically to the TSCA tool: repo with this file, formulation_core.py,
requirements.txt -> share.streamlit.io (or Railway).

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

from formulation_core import (solve_scale, batch_summary, charge_sheet_text)

st.set_page_config(page_title="Formulation Calculator", page_icon="⚗️",
                   layout="wide")

CACHE_FILE = Path(__file__).parent / "mw_cache.json"
PUBCHEM_PROP = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
                "{}/property/MolecularWeight,IUPACName,MolecularFormula/JSON")


# ---------------- PubChem MW lookup (cached) ----------------
def load_mw_cache():
    if CACHE_FILE.exists():
        try:
            return json.loads(CACHE_FILE.read_text())
        except json.JSONDecodeError:
            return {}
    return {}


def fetch_mw(identifier: str):
    """Return (mw, iupac_name, formula) or (None, None, None)."""
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
            iupac = p.get("IUPACName", "")
            formula = p.get("MolecularFormula", "")
            cache[key] = {"mw": mw, "iupac": iupac, "formula": formula}
            CACHE_FILE.write_text(json.dumps(cache, indent=1))
            time.sleep(0.22)
            return mw, iupac, formula
    except Exception:
        pass
    return None, None, None


# ---------------- session state ----------------
DEFAULT_ROWS = pd.DataFrame([
    {"Component": "", "CAS / name for lookup": "", "MW (g/mol)": None,
     "Assay %": 100.0, "Functionality": 2.0, "Group": "acid",
     "Molar ratio": 1.0},
])

if "table" not in st.session_state:
    st.session_state.table = DEFAULT_ROWS.copy()

st.title("⚗️ Lab Formulation Calculator")
st.caption("Molar ratios in → scalable gram charges, condensate prediction, "
           "theoretical end-group values, and stoichiometric feedback out. "
           "MW auto-fills from PubChem by CAS or name; enter MW manually for "
           "UVCBs (use effective MW or equivalent weight × functionality).")

# ---------------- component table ----------------
st.subheader("1 · Components")
edited = st.data_editor(
    st.session_state.table,
    num_rows="dynamic",
    use_container_width=True,
    column_config={
        "Component": st.column_config.TextColumn(required=True),
        "CAS / name for lookup": st.column_config.TextColumn(
            help="CAS number or chemical name; used for the PubChem MW lookup"),
        "MW (g/mol)": st.column_config.NumberColumn(
            min_value=0.0, format="%.2f",
            help="Auto-filled by lookup, or enter manually (UVCBs)"),
        "Assay %": st.column_config.NumberColumn(
            min_value=0.1, max_value=100.0, format="%.1f",
            help="Purity / solids of the as-supplied material"),
        "Functionality": st.column_config.NumberColumn(
            min_value=0.0, format="%.2f",
            help="Reactive groups per molecule (diacid = 2, DETA ≈ 2 primary "
                 "NH2 for amidation, glycerol = 3...)"),
        "Group": st.column_config.SelectboxColumn(
            options=["acid", "amine", "hydroxyl", "inert"],
            help="'inert' = solvent/non-reactive; excluded from stoichiometry"),
        "Molar ratio": st.column_config.NumberColumn(
            min_value=0.0, format="%.4f"),
    },
    key="editor",
)

if st.button("🔎 Fill missing MWs from PubChem"):
    misses = []
    for i, row in edited.iterrows():
        ident = (row["CAS / name for lookup"] or row["Component"] or "").strip()
        if ident and (pd.isna(row["MW (g/mol)"]) or not row["MW (g/mol)"]):
            mw, iupac, formula = fetch_mw(ident)
            if mw:
                edited.at[i, "MW (g/mol)"] = mw
                st.toast(f"{ident}: {mw:.2f} g/mol ({formula})")
            else:
                misses.append(ident)
    st.session_state.table = edited
    if misses:
        st.warning(f"No PubChem match for: {', '.join(misses)} — enter MW "
                   "manually (typical for UVCBs like dimer acids; use "
                   "effective MW).")
    st.rerun()

# ---------------- batch + reaction controls ----------------
st.subheader("2 · Scale anchor & reaction")
valid = edited.dropna(subset=["MW (g/mol)"])
valid = valid[(valid["Component"].astype(str).str.strip() != "") &
              (valid["Molar ratio"] > 0)]
names = valid["Component"].tolist()

c1, c2, c3, c4 = st.columns(4)
with c1:
    anchor_type = st.selectbox("Scale by", [
        "Total charge mass (g)", "Fixed component mass (g)",
        "Fixed component moles"])
with c2:
    anchor_value = st.number_input("Target value", min_value=0.0,
                                   value=500.0, format="%.3f")
with c3:
    anchor_comp = st.selectbox("Anchor component", names) if \
        anchor_type != "Total charge mass (g)" and names else None
with c4:
    reaction = st.selectbox("Condensation reaction", [
        "none", "amidation (acid+amine)", "esterification (acid+OH)"])

extent = 1.0
if reaction != "none":
    extent = st.slider("Conversion of limiting group (p)", 0.50, 1.00, 1.00,
                       0.005, help="1.00 = run to completion; lower to model "
                       "a target conversion / staged cook")

# ---------------- compute ----------------
if len(valid) and st.button("Calculate batch", type="primary",
                            use_container_width=True):
    comps = [{
        "name": r["Component"], "cas": r["CAS / name for lookup"],
        "mw": float(r["MW (g/mol)"]), "assay": float(r["Assay %"]) / 100.0,
        "functionality": float(r["Functionality"]),
        "group": r["Group"], "ratio": float(r["Molar ratio"]),
    } for _, r in valid.iterrows()]

    atype = {"Total charge mass (g)": "total_asis_mass",
             "Fixed component mass (g)": "component_asis_mass",
             "Fixed component moles": "component_moles"}[anchor_type]
    rxn = {"none": "none", "amidation (acid+amine)": "amidation",
           "esterification (acid+OH)": "esterification"}[reaction]

    try:
        scale = solve_scale(comps, atype, anchor_value, anchor_comp)
        out = batch_summary(comps, scale, reaction=rxn, extent=extent)
    except Exception as e:
        st.error(f"Calculation failed: {e}")
        st.stop()

    st.subheader("3 · Charge table")
    df = pd.DataFrame([{
        "Component": r["name"], "CAS": r.get("cas", ""),
        "Moles": round(r["moles"], 4), "Equivalents": round(r["eq"], 4),
        "Active g": round(r["g_active"], 2),
        "As-is g to weigh": round(r["g_asis"], 2),
    } for r in out["rows"]])
    st.dataframe(df, use_container_width=True, hide_index=True)

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Total charge (as-is)", f"{out['charge_mass']:.1f} g")
    if rxn != "none":
        m2.metric("Predicted condensate", f"{out['cond']['water_g']:.2f} g H₂O")
        m3.metric("Theoretical final mass", f"{out['final_mass']:.1f} g")
        m4.metric("Theoretical solids", f"{out['solids_pct']:.1f}%")

    # ---- feedback panel ----
    ev, car, t = out["end_values"], out["carothers"], out["totals"]
    if rxn != "none":
        st.subheader("4 · Stoichiometric feedback")
        f1, f2, f3 = st.columns(3)
        if ev:
            f1.metric("Theoretical acid value", f"{ev['acid_value']:.1f}",
                      help="mg KOH / g resin — compare to titration")
            f2.metric("Theoretical amine value", f"{ev['amine_value']:.1f}")
            f3.metric("Theoretical OH value", f"{ev['hydroxyl_value']:.1f}")
        if car:
            st.markdown(
                f"**Stoichiometric balance:** r = `{car['r']:.3f}` "
                f"({car['excess_pct']:.1f}% excess **{car['excess_group']}**)")
            xn_s = "∞" if car["Xn"] == float("inf") else f"{car['Xn']:.1f}"
            mn_s = "∞" if car["Mn"] == float("inf") else f"{car['Mn']:.0f} g/mol"
            st.info(f"At p = {extent:.3f}: expected **Xn ≈ {xn_s}**, "
                    f"**Mn ≈ {mn_s}**, chains predominantly "
                    f"**{car['excess_group']}-terminated**.")
            if car["r"] > 0.995 and extent > 0.98:
                st.warning("⚠️ Near-perfect stoichiometry at high conversion — "
                           "Mn is unbounded. If this is a thermoset/gel target, "
                           "fine; if you want a stable liquid resin, build in "
                           "an end-group excess or cap conversion.")

    # ---- charge sheet ----
    st.subheader("5 · Charge sheet")
    h1, h2, h3 = st.columns(3)
    meta = {
        "project": h1.text_input("Project code", ""),
        "exp_id": h2.text_input("Experiment ID", ""),
        "chemist": h3.text_input("Chemist", ""),
        "date": str(date.today()),
        "description": st.text_input("Batch description", ""),
    }
    sheet = charge_sheet_text(out, meta)
    st.code(sheet, language=None)
    st.download_button("Download charge sheet (.txt)", sheet,
                       file_name=f"charge_{meta['exp_id'] or 'batch'}.txt")

st.divider()
st.caption("Verify MWs and assays against CoAs. Theoretical values assume "
           "ideal step-growth behavior; real cooks deviate (side reactions, "
           "volatilized amine, incomplete condensate recovery). This is a "
           "planning tool, not a CoA.")
