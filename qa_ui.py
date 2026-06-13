"""QA / QC tools — bench/plant calculators, kept deliberately simple for techs.

Two independent tools, each reachable as its OWN sidebar entry so a tech sees
only the controls for the job in front of them:
  - '% Solids adjustment' : exact mass-balance solids raise/dilute (Adjust +
    Compare carriers), lab- and plant-scale.
  - 'pH adjustment'       : dose scaler from a bench titration, with automatic
    titrant/sample concentration + dilution correction.

render(section) renders one tool ('solids' / 'ph') or both ('all').
Shares st.session_state.changelog with the rest of the app.
"""
from datetime import date

import pandas as pd
import streamlit as st

from formulation_core import (solve_solids_adjustment,
                              compare_solids_carriers, qa_solids_sheet_text)

DEFAULT_STREAMS = pd.DataFrame([
    {"Add stream": "60% resin", "% solids": 60.0, "Mode": "solve",
     "Amount": 1.0},
])
DEFAULT_CANDS = pd.DataFrame([
    {"Carrier": "60% resin", "% solids": 60.0},
    {"Carrier": "45% dispersion", "% solids": 45.0},
])
MODE_OPTS = ["fixed", "solve"]

PRESET_FMTS = ["% w/w", "% w/v", "Molarity (mol/L)", "Normality (eq/L)"]


def _titrant_inputs(prefix):
    """Render a titrant concentration block; return eq per gram of solution
    (or None if inputs incomplete). Preset path = MW+eq baked in; direct path
    = eq/kg or eq/L."""
    from formulation_core import (TITRANT_PRESETS, titrant_eq_per_g,
                                  ph_format_needs)
    src = st.radio("Concentration source",
                   ["Preset titrant", "Direct equivalents"],
                   key=f"{prefix}_src", horizontal=True)
    eqpg = None
    if src == "Preset titrant":
        names = list(TITRANT_PRESETS.keys()) + ["Other (enter MW + eq)"]
        choice = st.selectbox("Titrant", names, key=f"{prefix}_preset")
        if choice == "Other (enter MW + eq)":
            o1, o2 = st.columns(2)
            mw = o1.number_input("MW (g/mol)", min_value=0.01, value=40.0,
                                 format="%.2f", key=f"{prefix}_mw")
            eqm = o2.number_input("Equivalents per mole", min_value=1,
                                  value=1, step=1, key=f"{prefix}_eqm")
        else:
            mw, eqm = TITRANT_PRESETS[choice]
        fmt = st.selectbox("Concentration format", PRESET_FMTS,
                           key=f"{prefix}_fmt")
        needs_d, _ = ph_format_needs(fmt)
        c1, c2 = st.columns(2)
        val = c1.number_input(f"Value ({fmt})", min_value=0.0, value=50.0,
                              format="%.4f", key=f"{prefix}_val")
        dens = 1.0
        if needs_d:
            dens = c2.number_input("Solution density (g/mL)", min_value=0.1,
                                   value=1.000, format="%.3f",
                                   key=f"{prefix}_dens",
                                   help="Required for volume-based formats "
                                        "(% w/v, Molarity, Normality).")
        try:
            eqpg = titrant_eq_per_g(fmt, val, eqm, mw, dens)
        except Exception:
            eqpg = None
    else:
        d1, d2 = st.columns(2)
        unit = d1.selectbox("Unit", ["eq/kg", "eq/L"], key=f"{prefix}_du")
        v = d2.number_input(f"Value ({unit})", min_value=0.0, value=12.5,
                            format="%.4f", key=f"{prefix}_dv")
        if unit == "eq/L":
            dd = st.number_input("Solution density (g/mL)", min_value=0.1,
                                 value=1.000, format="%.3f", key=f"{prefix}_dd")
            eqpg = titrant_eq_per_g("eq/L", v, density=dd)
        else:
            eqpg = titrant_eq_per_g("eq/kg", v)
    if eqpg:
        st.caption(f"→ {eqpg:.5g} eq per gram of solution")
    return eqpg


def _solids_section():
    st.subheader("% Solids adjustment")
    st.caption("Exact mass balance, unit-agnostic (g at the bench, lb/kg at "
               "the kettle — output matches input). Raising low solids and "
               "diluting high solids both supported.")

    # ---------------- batch ----------------
    b1, b2, b3, b4 = st.columns(4)
    with b1:
        batch_mass = st.number_input("Batch size", min_value=0.0,
                                     value=1000.0, format="%.2f")
    with b2:
        unit = st.selectbox("Unit", ["g", "kg", "lb"], index=2)
    with b3:
        cur_pct = st.number_input("Current % solids", min_value=0.0,
                                  max_value=100.0, value=28.0, format="%.2f")
    with b4:
        tgt_pct = st.number_input("Target % solids", min_value=0.0,
                                  max_value=100.0, value=35.0, format="%.2f")

    # =================== ADJUST (solve) ===================
    st.markdown("#### Adjust — solve the addition")
    st.caption("**fixed** = an addition you must make in a known amount. "
               "**solve** = the tool computes it. ONE solve stream → a single "
               "carrier. MULTIPLE solve streams → carriers added **together "
               "in fixed proportion** (their Amount = relative ratio weight, "
               "e.g. 2 and 1 for 2:1); the tool solves the common scale.")
    streams_df = st.data_editor(
        st.session_state.qa_streams, num_rows="dynamic",
        width="stretch",
        column_config={
            "Add stream": st.column_config.TextColumn(required=True),
            "% solids": st.column_config.NumberColumn(
                min_value=0.0, max_value=100.0, format="%.1f",
                help="Solids content of THIS stream (water/solvent = 0)"),
            "Mode": st.column_config.SelectboxColumn(
                options=MODE_OPTS,
                help="'solve' for the carrier(s); 'fixed' for known adds"),
            "Amount": st.column_config.NumberColumn(
                min_value=0.0, format="%.3f",
                help="Fixed stream: as-added amount in the batch unit. "
                     "Solve stream: relative ratio weight when combining "
                     "multiple carriers (single solve ignores this)."),
        }, key=f"qa_streams_v{st.session_state.qa_ver}")

    sv = streams_df[streams_df["Add stream"].astype(str).str.strip() != ""]
    streams = [{"name": r["Add stream"],
                "solids_pct": float(r["% solids"] or 0),
                "mode": str(r["Mode"]).strip() or "fixed",
                "amount": float(r["Amount"] or 0)}
               for _, r in sv.iterrows()]

    if st.button("Solve adjustment", type="primary"):
        if not streams:
            st.info("Add at least one stream (one marked 'solve').")
        else:
            res = solve_solids_adjustment(batch_mass, cur_pct, tgt_pct,
                                          streams)
            if res.get("warning"):
                st.session_state.qa_result = None
                st.warning(res["warning"])
            else:
                st.session_state.qa_result = {
                    "batch_mass": batch_mass, "unit": unit, "cur_pct": cur_pct,
                    "tgt_pct": tgt_pct, "streams": streams, "res": res}
                names = ", ".join(res["solved_names"])
                st.session_state.changelog.append(
                    f"[{date.today()}] QA solids {cur_pct:.2f}%→{tgt_pct:.2f}% "
                    f"on {batch_mass:.1f} {unit}: add "
                    f"{res['solved_total']:.2f} {unit} ({names}) → "
                    f"{res['final_mass']:.1f} {unit} @ "
                    f"{res['final_solids_pct']:.2f}%")

    R = st.session_state.get("qa_result")
    if R:
        res = R["res"]
        u = R["unit"]
        st.subheader("Result")
        m1, m2, m3 = st.columns(3)
        m1.metric(f"Total carrier added", f"{res['solved_total']:.2f} {u}",
                  help="Sum across solve streams" if
                  len(res["solved_names"]) > 1 else None)
        m2.metric("Final batch mass", f"{res['final_mass']:.2f} {u}")
        m3.metric("Final % solids", f"{res['final_solids_pct']:.2f}%")
        add_df = pd.DataFrame({
            "Add stream": [s["name"] for s in R["streams"]],
            "% solids": [s["solids_pct"] for s in R["streams"]],
            "Mode": [s["mode"] for s in R["streams"]],
            f"Add ({u})": [round(m, 3) for m in res["stream_masses"]],
        })
        st.dataframe(add_df, width="stretch", hide_index=True)
        st.caption(f"Start {R['batch_mass']:.2f} {u} @ {R['cur_pct']:.2f}% → "
                   f"add the above → {res['final_mass']:.2f} {u} @ "
                   f"{res['final_solids_pct']:.2f}%. Exact mass balance — "
                   f"confirm with a measured solids check, not a yield model.")

        st.subheader("Charge sheet")
        h1, h2, h3 = st.columns(3)
        meta = {"project": h1.text_input("Project code", "", key="qa_p"),
                "exp_id": h2.text_input("Batch ID", "", key="qa_e"),
                "chemist": h3.text_input("Chemist", "", key="qa_c"),
                "date": str(date.today())}
        sheet = qa_solids_sheet_text(
            R["batch_mass"], u, R["cur_pct"], R["tgt_pct"],
            [s["name"] for s in R["streams"]],
            [s["solids_pct"] for s in R["streams"]],
            res["stream_masses"], res["final_mass"],
            res["final_solids_pct"], meta)
        st.code(sheet, language=None)
        from pdf_util import text_to_pdf
        dl1, dl2 = st.columns(2)
        dl1.download_button("Download charge sheet (.txt)", sheet,
                            file_name=f"qa_solids_{meta['exp_id'] or 'batch'}"
                                      f".txt")
        dl2.download_button("⬇️ Download charge sheet (.pdf)",
                            text_to_pdf(sheet, "QA % solids charge sheet"),
                            file_name=f"qa_solids_{meta['exp_id'] or 'batch'}"
                                      f".pdf", mime="application/pdf")

    # =================== COMPARE candidates ===================
    st.markdown("#### Compare carriers — raise low solids")
    st.caption("List candidate solids carriers; the tool shows how much of "
               "**each one alone** reaches the target, so you can pick. Uses "
               "the batch size / current / target set above. Carriers at or "
               "below the target are flagged (can't raise solids).")
    cands_df = st.data_editor(
        st.session_state.qa_cands, num_rows="dynamic",
        width="stretch",
        column_config={
            "Carrier": st.column_config.TextColumn(required=True),
            "% solids": st.column_config.NumberColumn(
                min_value=0.0, max_value=100.0, format="%.1f"),
        }, key=f"qa_cands_v{st.session_state.qa_ver}")
    cv = cands_df[cands_df["Carrier"].astype(str).str.strip() != ""]
    candidates = [{"name": r["Carrier"],
                   "solids_pct": float(r["% solids"] or 0)}
                  for _, r in cv.iterrows()]

    if st.button("Compare carriers"):
        if not candidates:
            st.info("List at least one candidate carrier.")
        else:
            rows = compare_solids_carriers(batch_mass, cur_pct, tgt_pct,
                                           candidates)
            disp = []
            for row in rows:
                if row["feasible"]:
                    disp.append({
                        "Carrier": row["name"],
                        "% solids": row["solids_pct"],
                        f"Add ({unit})": round(row["amount"], 2),
                        f"Final mass ({unit})": round(row["final_mass"], 1),
                        "Final % solids": round(row["final_pct"], 2),
                        "Feasible": "✓"})
                else:
                    disp.append({
                        "Carrier": row["name"],
                        "% solids": row["solids_pct"],
                        f"Add ({unit})": None,
                        f"Final mass ({unit})": None,
                        "Final % solids": None,
                        "Feasible": "✗ at/below target"})
            st.dataframe(pd.DataFrame(disp), width="stretch",
                         hide_index=True)
            feas = [r for r in rows if r["feasible"]]
            if feas:
                lo = min(feas, key=lambda r: r["amount"])
                st.caption(f"Least material to reach {tgt_pct:.2f}%: "
                           f"**{lo['name']}** at {lo['amount']:.2f} {unit}. "
                           f"(Lower-solids carriers need more mass and add "
                           f"more volume.)")

    st.caption("% solids = exact mass balance (confirm with a measured solids "
               "check). Planning / QC tool, not a CoA.")


def _ph_section():
    st.subheader("pH adjustment — dose from a bench titration")
    st.caption("Titrate a weighed sample to your target endpoint; that "
               "captures the batch's own buffering. The tool scales it to the "
               "batch and **corrects for any change in titrant or sample "
               "concentration** — the step that gets missed. It does not "
               "predict pH from volume.")

    st.markdown("**1 · Calibration titration** (on the bench sample)")
    cc1, cc2, cc3 = st.columns(3)
    with cc1:
        samp_mass = st.number_input("Sample mass titrated (g)",
                                    min_value=0.0, value=200.0, format="%.3f")
    with cc2:
        prod_pct = st.number_input("% product in sample", min_value=0.0,
                                   max_value=100.0, value=100.0, format="%.2f",
                                   help="100 = neat. If the titrated sample "
                                        "was diluted, enter the % that is neat "
                                        "product — the dose is referenced to "
                                        "neat product, not the diluted mass.")
    with cc3:
        titr_g = st.number_input("Titrant used (g)", min_value=0.0,
                                 value=4.0, format="%.4f",
                                 help="Grams of titrant to reach the target "
                                      "endpoint on this sample.")
    st.markdown("Calibration titrant:")
    eqpg_cal = _titrant_inputs("ph_cal")

    st.markdown("**2 · Dosing titrant** (what you'll add to the batch)")
    same = st.checkbox("Same titrant & strength as calibration", value=True,
                       key="ph_same")
    if same:
        eqpg_dose = eqpg_cal
    else:
        eqpg_dose = _titrant_inputs("ph_dose")
    diluted = st.checkbox("Dosing titrant is cut from stock before use",
                          value=False, key="ph_dil")
    if diluted:
        x1, x2 = st.columns(2)
        ps = x1.number_input("Parts stock (by mass)", min_value=0.0,
                             value=1.0, format="%.3f", key="ph_ps")
        pd_ = x2.number_input("Parts diluent (by mass)", min_value=0.0,
                             value=1.0, format="%.3f", key="ph_pd")
        if eqpg_dose and (ps + pd_) > 0:
            frac = ps / (ps + pd_)
            eqpg_dose = eqpg_dose * frac
            st.caption(f"Cut to {100*frac:.1f}% stock by mass → "
                       f"{eqpg_dose:.5g} eq/g working strength "
                       f"(C₁V₁=C₂V₂, mass basis).")

    st.markdown("**3 · Batch to adjust**")
    bb1, bb2 = st.columns(2)
    with bb1:
        ph_batch = st.number_input("Batch size (neat product)",
                                   min_value=0.0, value=2000.0, format="%.2f",
                                   key="ph_batch")
    with bb2:
        ph_unit = st.selectbox("Batch unit", ["g", "kg", "lb"], index=2,
                               key="ph_unit")
    factor = {"g": 1.0, "kg": 1000.0, "lb": 453.59237}[ph_unit]

    if st.button("Compute dose", type="primary", key="ph_go"):
        from formulation_core import ph_dose_from_titration
        if not eqpg_cal or not eqpg_dose:
            st.warning("Complete both titrant concentration blocks first.")
        else:
            r = ph_dose_from_titration(samp_mass, prod_pct, titr_g, eqpg_cal,
                                       ph_batch * factor, eqpg_dose)
            if r.get("warning"):
                st.warning(r["warning"])
            else:
                # dose returned in grams; convert to batch unit for display
                corr_u = r["corrected_dose_g"] / factor
                man_u = r["manual_dose_g"] / factor
                st.subheader("Dose")
                p1, p2, p3 = st.columns(3)
                p1.metric("Corrected dose", f"{corr_u:.3f} {ph_unit}",
                          help="Concentration- and neat-product-corrected.")
                p2.metric("Naive manual dose", f"{man_u:.3f} {ph_unit}",
                          help="(g titrant / g sample) × batch — no "
                               "correction.")
                err = r["pct_error"]
                if abs(err) < 0.05:
                    p3.metric("Manual error", "≈ 0%")
                    st.success("No concentration/dilution change — manual and "
                               "corrected agree.")
                else:
                    word = "OVER-dosing" if err > 0 else "UNDER-dosing"
                    p3.metric("Manual error", f"{err:+.1f}%",
                              delta=f"{word}", delta_color="inverse")
                    st.warning(f"The uncorrected manual figure would be "
                               f"**{word.lower()} by {abs(err):.1f}%** — the "
                               f"titrant or sample concentration changed. Use "
                               f"the corrected dose.")
                st.caption(f"Neat product in sample: {r['w_neat']:.2f} g · "
                           f"specific demand {r['specific_demand']:.5g} eq per "
                           f"g product · batch needs {r['E_batch']:.4g} eq. "
                           f"Assumes the same titrant chemistry; if the "
                           f"batch's composition differs materially from the "
                           f"sample, re-titrate a fresh sample.")
                st.session_state.changelog.append(
                    f"[{date.today()}] QA pH dose on {ph_batch:.1f} {ph_unit}: "
                    f"corrected {corr_u:.3f} {ph_unit} (manual {man_u:.3f}, "
                    f"err {err:+.1f}%)")

    st.caption("pH dose scales your bench titration with exact concentration "
               "math; it does not predict pH from volume, and a material "
               "change in batch composition means re-titrate. Planning / QC "
               "tool, not a CoA.")


def _changelog():
    if st.session_state.get("changelog"):
        with st.expander(f"📜 Change log ({len(st.session_state.changelog)})"):
            for line in reversed(st.session_state.changelog):
                st.text(line)


def render(section="all"):
    """Render one QA tool ('solids' or 'ph') or both ('all')."""
    if "qa_ver" not in st.session_state:
        st.session_state.qa_ver = 0
    if "qa_streams" not in st.session_state:
        st.session_state.qa_streams = DEFAULT_STREAMS.copy()
    if "qa_cands" not in st.session_state:
        st.session_state.qa_cands = DEFAULT_CANDS.copy()
    if "changelog" not in st.session_state:
        st.session_state.changelog = []

    if section in ("all", "solids"):
        _solids_section()
    if section == "all":
        st.divider()
    if section in ("all", "ph"):
        _ph_section()
    _changelog()
