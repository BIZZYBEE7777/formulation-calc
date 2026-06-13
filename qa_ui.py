"""QA / QC tools tab.
Loop 1: % solids adjustment (exact mass balance), lab- and plant-scale.
  - Adjust: fixed additions + one-or-more 'solve' carriers. One solve = a
    single carrier; multiple solves = carriers combined in fixed ratio.
  - Compare: each candidate carrier evaluated individually to raise low solids.
Loop 2 (next): pH dose scaler — sample titration ratio x batch weight with
  automatic titrant/sample concentration correction + format conversion.
Mirrors blend_ui.render(); shares st.session_state.changelog.
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


def render():
    if "qa_ver" not in st.session_state:
        st.session_state.qa_ver = 0
    if "qa_streams" not in st.session_state:
        st.session_state.qa_streams = DEFAULT_STREAMS.copy()
    if "qa_cands" not in st.session_state:
        st.session_state.qa_cands = DEFAULT_CANDS.copy()
    if "changelog" not in st.session_state:
        st.session_state.changelog = []

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
        use_container_width=True,
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
        st.dataframe(add_df, use_container_width=True, hide_index=True)
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
        st.download_button("Download charge sheet (.txt)", sheet,
                           file_name=f"qa_solids_{meta['exp_id'] or 'batch'}"
                                     f".txt")

    # =================== COMPARE candidates ===================
    st.markdown("#### Compare carriers — raise low solids")
    st.caption("List candidate solids carriers; the tool shows how much of "
               "**each one alone** reaches the target, so you can pick. Uses "
               "the batch size / current / target set above. Carriers at or "
               "below the target are flagged (can't raise solids).")
    cands_df = st.data_editor(
        st.session_state.qa_cands, num_rows="dynamic",
        use_container_width=True,
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
            st.dataframe(pd.DataFrame(disp), use_container_width=True,
                         hide_index=True)
            feas = [r for r in rows if r["feasible"]]
            if feas:
                lo = min(feas, key=lambda r: r["amount"])
                st.caption(f"Least material to reach {tgt_pct:.2f}%: "
                           f"**{lo['name']}** at {lo['amount']:.2f} {unit}. "
                           f"(Lower-solids carriers need more mass and add "
                           f"more volume.)")

    if st.session_state.changelog:
        with st.expander(f"📜 Change log ({len(st.session_state.changelog)})"):
            for line in reversed(st.session_state.changelog):
                st.text(line)

    st.divider()
    st.info("**pH adjustment — next build loop.** Sample-titration dose "
            "scaler: (g titrant / g neat product) × batch weight, with "
            "**automatic correction when titrant or sample concentration "
            "changes** (format conversion % w/w · % w/v · M · N · eq/L, plus "
            "titrant cut C₁V₁=C₂V₂ and sample-dilution basis). Shows the "
            "manual number beside the corrected one so the forgotten "
            "concentration adjustment is impossible to miss. It will not "
            "predict pH from volume, and a material change in batch "
            "composition means re-titrate a fresh sample.")
    st.caption("Planning/QC tool, not a CoA. % solids adjustment is exact "
               "mass balance; confirm with a measured solids check.")
