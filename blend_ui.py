"""Cold blend UI — wt% + per-ingredient solids. Two input methods:
   1) Solve from constraints (fixed wt%, solids-ratio group, balance)
   2) Known parts (enter all weights directly)
Shared downstream: blend table, batch scaling, charge sheet, change log."""
from datetime import date

import pandas as pd
import streamlit as st

from formulation_core import (blend_summary, blend_solve_diluent,
                              blend_batch_grams, blend_sheet_text,
                              solve_blend, US_FLOZ_ML)

DEFAULT_SPECS = pd.DataFrame([
    {"Ingredient": "Dispersion A", "% solids": 30.0,
     "Role": "ratio (solids carrier)", "Value": 2.0},
    {"Ingredient": "Dispersion B", "% solids": 30.0,
     "Role": "ratio (solids carrier)", "Value": 1.0},
    {"Ingredient": "Acetone", "% solids": 0.0,
     "Role": "fixed wt%", "Value": 10.0},
    {"Ingredient": "Water", "% solids": 0.0,
     "Role": "balance (fills to 100%)", "Value": 0.0},
])

DEFAULT_PARTS = pd.DataFrame([
    {"Ingredient": "Eastek 1400", "% solids": 30.0, "Parts (by weight)": 50.0},
    {"Ingredient": "Water", "% solids": 0.0, "Parts (by weight)": 50.0},
])

ROLE_MAP = {"fixed wt%": "fixed", "ratio (solids carrier)": "ratio",
            "balance (fills to 100%)": "balance"}

# Typical specific-gravity RANGES by resin family (g/mL). Reference only —
# prefills the density field; the user measures/overrides. NOT a spec.
SG_LIBRARY = {
    "Long-oil alkyd": (0.93, 0.97),
    "Medium-oil alkyd": (0.95, 0.99),
    "Short-oil alkyd": (0.98, 1.04),
    "Imidazoline": (0.94, 0.98),
    "Polyamide resin": (0.95, 0.99),
    "UPR (in styrene)": (1.05, 1.15),
    "Vinyl ester (in styrene)": (1.04, 1.10),
    "Acrylic polyol (solvent)": (1.02, 1.08),
    "Urethane prepolymer": (1.05, 1.15),
    "Waterborne dispersion": (1.02, 1.06),
    "Water": (1.00, 1.00),
}


def render():
    if "blend_specs" not in st.session_state:
        st.session_state.blend_specs = DEFAULT_SPECS.copy()
    if "blend_table" not in st.session_state:
        st.session_state.blend_table = DEFAULT_PARTS.copy()
    if "blend_result" not in st.session_state:
        st.session_state.blend_result = None
    if "blend_ver" not in st.session_state:
        st.session_state.blend_ver = 0

    method = st.radio("Blend input", ["Solve from constraints", "Known parts"],
                      horizontal=True,
                      help="Solve: tell it the rules (acetone = 10% of final, "
                           "dispersions 2:1, water fills the rest, hit 10% "
                           "solids) and it computes the recipe. Known parts: "
                           "you already have the weights.")

    names = solids = parts = None

    # ---------------- method 1: constraint solver ----------------
    if method == "Solve from constraints":
        st.subheader("1 · Ingredients & rules")
        st.caption("Roles — **fixed wt%**: ingredient is exactly this % of "
                   "the final blend (Value = wt%). **ratio (solids "
                   "carrier)**: held at this relative weight ratio vs other "
                   "carriers, scaled to hit the solids target (Value = ratio "
                   "number, e.g. 2 and 1). **balance**: exactly one, fills "
                   "to 100% (Value ignored) — usually water.")
        spec_df = st.data_editor(
            st.session_state.blend_specs, num_rows="dynamic",
            width="stretch",
            column_config={
                "Ingredient": st.column_config.TextColumn(required=True),
                "% solids": st.column_config.NumberColumn(
                    min_value=0.0, max_value=100.0, format="%.1f"),
                "Role": st.column_config.SelectboxColumn(
                    options=list(ROLE_MAP.keys())),
                "Value": st.column_config.NumberColumn(min_value=0.0,
                                                       format="%.3f"),
            }, key=f"spec_editor_v{st.session_state.blend_ver}")

        tcol, bcol = st.columns([1, 2])
        with tcol:
            target = st.number_input("Target % solids of final blend",
                                     min_value=0.1, max_value=99.0,
                                     value=10.0, format="%.2f")
        with bcol:
            st.write("")
            st.write("")
            if st.button("Solve recipe", type="primary"):
                sv = spec_df[spec_df["Ingredient"].astype(str)
                             .str.strip() != ""]
                specs = [{"name": r["Ingredient"],
                          "solids_pct": float(r["% solids"] or 0),
                          "role": ROLE_MAP.get(r["Role"], "fixed"),
                          "value": float(r["Value"] or 0)}
                         for _, r in sv.iterrows()]
                solved, warn = solve_blend(specs, target)
                if warn:
                    st.session_state.blend_result = None
                    st.error(warn)
                else:
                    st.session_state.blend_result = {
                        "names": [s["name"] for s in specs],
                        "solids": [s["solids_pct"] for s in specs],
                        "parts": solved, "target": target}
                    st.session_state.changelog.append(
                        f"[{date.today()}] Solved blend to {target:.2f}% "
                        f"solids: " + ", ".join(
                            f"{s['name']} {p:.2f}"
                            for s, p in zip(specs, solved)) + " (per 100 g)")
        if st.session_state.blend_result:
            R = st.session_state.blend_result
            names, solids, parts = R["names"], R["solids"], R["parts"]

    # ---------------- method 2: known parts ----------------
    else:
        st.subheader("1 · Ingredients (weight parts)")
        st.caption("Any consistent weight basis — 50/50, 2:1:1, real grams, "
                   "wt%. % solids per TDS; solvents/water = 0.")
        bt = st.data_editor(
            st.session_state.blend_table, num_rows="dynamic",
            width="stretch",
            column_config={
                "Ingredient": st.column_config.TextColumn(required=True),
                "% solids": st.column_config.NumberColumn(
                    min_value=0.0, max_value=100.0, format="%.1f"),
                "Parts (by weight)": st.column_config.NumberColumn(
                    min_value=0.0, format="%.3f"),
            }, key=f"blend_editor_v{st.session_state.blend_ver}")
        bv = bt[(bt["Ingredient"].astype(str).str.strip() != "") &
                (bt["Parts (by weight)"].fillna(0) > 0)]
        if len(bv):
            names = bv["Ingredient"].tolist()
            solids = [float(x) for x in bv["% solids"].fillna(0)]
            parts = [float(x) for x in bv["Parts (by weight)"]]

    if not names:
        st.info("Enter ingredients above" +
                (" and press **Solve recipe**." if method.startswith("Solve")
                 else "."))
        return

    # ---------------- shared downstream ----------------
    wt, blend_pct, _ = blend_summary(parts, solids)
    st.subheader("2 · Blend")
    bdf = pd.DataFrame({
        "Ingredient": names, "% solids": solids,
        "Parts / g per 100 g": [round(p, 3) for p in parts],
        "wt%": [round(100 * w, 2) for w in wt],
        "Solids contribution (g/100 g)":
            [round(w * s, 2) for w, s in zip(wt, solids)],
    })
    st.dataframe(bdf, width="stretch", hide_index=True)
    st.metric("Blend solids", f"{blend_pct:.2f}%")

    # simple re-target panel (known-parts path only; solver has its own)
    if method == "Known parts":
        st.subheader("3 · Hit a target % solids (adjust one ingredient)")
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
                        f"[{date.today()}] Blend {blend_pct:.2f}% → "
                        f"{tgt:.2f}% solids: {dil} {parts[di]:.3f} → "
                        f"{new_parts[di]:.3f} parts")
                    tbl = st.session_state.blend_table.copy()
                    for i2 in tbl.index:
                        nm = str(tbl.at[i2, "Ingredient"]).strip()
                        if nm in names:
                            tbl.at[i2, "Parts (by weight)"] = \
                                new_parts[names.index(nm)]
                    st.session_state.blend_table = tbl
                    st.session_state.blend_ver += 1
                    st.rerun()

    # ---------------- batch scaling ----------------
    st.subheader("4 · Scale to batch size" if method == "Known parts"
                 else "3 · Scale to batch size")
    sg_pick = st.selectbox(
        "Typical SG prefill (optional — measure yours)",
        ["— (use 1.00)"] + list(SG_LIBRARY.keys()),
        help="Per-family TYPICAL specific-gravity ranges. Sets the density "
             "default below; override with your measured value. Reference "
             "only, not a spec.")
    if sg_pick != "— (use 1.00)":
        lo, hi = SG_LIBRARY[sg_pick]
        sg_default = round((lo + hi) / 2.0, 3)
        st.caption(f"{sg_pick}: typical SG {lo:.2f}–{hi:.2f} g/mL "
                   f"(prefilled {sg_default:.3f} — measure yours).")
    else:
        sg_default = 1.00
    s1, s2 = st.columns(2)
    with s1:
        density = st.number_input(
            "Blend density (g/mL) — for volume sizes", min_value=0.5,
            max_value=2.5, value=sg_default, format="%.3f",
            help="Sets volume-size scaling. Prefilled from the typical-SG "
                 "picker above; measure for tight work. Mass sizes ignore "
                 "this.")
    with s2:
        custom_g = st.number_input("Custom batch (g, 0 = off)",
                                   min_value=0.0, value=0.0, format="%.1f")
    sizes = [("4 oz", 4 * US_FLOZ_ML * density),
             ("Pint", 16 * US_FLOZ_ML * density),
             ("Quart", 32 * US_FLOZ_ML * density),
             ("Gallon", 128 * US_FLOZ_ML * density),
             ("100 g", 100.0), ("1 kg", 1000.0)]
    if custom_g > 0:
        sizes.append((f"{custom_g:.0f} g", custom_g))
    scale_df = pd.DataFrame({"Ingredient": names, "% solids": solids})
    for label, gtot in sizes:
        scale_df[f"{label} (g)"] = [round(x, 2)
                                    for x in blend_batch_grams(parts, gtot)]
    st.dataframe(scale_df, width="stretch", hide_index=True)

    # ---------------- charge sheet ----------------
    st.subheader("Charge sheet")
    pick = st.selectbox("Sheet for batch size", [lab for lab, _ in sizes])
    gtot = dict(sizes)[pick]
    h1, h2, h3 = st.columns(3)
    meta = {"project": h1.text_input("Project code", "", key="bp"),
            "exp_id": h2.text_input("Batch ID", "", key="be"),
            "chemist": h3.text_input("Chemist", "", key="bc"),
            "date": str(date.today())}
    sheet = blend_sheet_text(names, solids, blend_batch_grams(parts, gtot),
                             pick, blend_pct, meta)
    st.code(sheet, language=None)
    from pdf_util import text_to_pdf
    dl1, dl2 = st.columns(2)
    dl1.download_button("Download charge sheet (.txt)", sheet,
                        file_name=f"blend_{meta['exp_id'] or 'batch'}.txt")
    dl2.download_button("⬇️ Download charge sheet (.pdf)",
                        text_to_pdf(sheet, "Cold blend charge sheet"),
                        file_name=f"blend_{meta['exp_id'] or 'batch'}.pdf",
                        mime="application/pdf")

    if st.session_state.changelog:
        with st.expander(f"📜 Change log ({len(st.session_state.changelog)})"):
            for line in reversed(st.session_state.changelog):
                st.text(line)

    st.divider()
    st.caption("Blend solids = Σ(wt% × ingredient solids). Volume sizes use "
               "the entered density. Planning tool, not a CoA.")
