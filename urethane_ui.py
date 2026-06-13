"""2K urethane tab — resin (polyol) + isocyanate hardener mix ratio to a
target NCO:OH index. A convenience FRONT-DOOR over the same step-growth engine
(isocyanate + OH is an addition bond, no condensate); it adds no new chemistry,
just the stoichiometry every 2K formulation needs.

Mirrors blend_ui.render() / qa_ui.render(); shares st.session_state.changelog.
For a full multi-component urethane charge (water/byproduct accounting, %NCO of
the prepolymer, Carothers Mn), use the main 'Molar ratios -> batch' mode with
the 'urethane (NCO+OH)' reaction and 'isocyanate' group — this tab is the quick
two-component ratio answer.
"""
from datetime import date

import streamlit as st

from formulation_core import (solve_2k_index_ratio, nco_eq_per_g, oh_eq_per_g)


def render():
    if "changelog" not in st.session_state:
        st.session_state.changelog = []

    st.subheader("2K urethane — mix ratio to a target NCO:OH index")
    st.caption("Enter the resin's hydroxyl value and the hardener's %NCO (both "
               "off the CoA); the tool solves the weight ratio to hit a target "
               "index. **Index 100 = stoichiometric**, 105-110 = typical "
               "NCO-rich. Same engine as every family — NCO+OH is just an "
               "addition bond; this is the two-component shortcut.")

    c1, c2 = st.columns(2)
    with c1:
        st.markdown("**Resin (polyol)**")
        resin_name = st.text_input("Name", "Polyol resin", key="u_rn")
        oh_src = st.radio("Hydroxyl spec", ["Hydroxyl value (mg KOH/g)",
                                            "OH equivalent weight (g/eq)"],
                          key="u_ohsrc")
        if oh_src.startswith("Hydroxyl value"):
            resin_ohv = st.number_input("OHV (mg KOH/g)", min_value=0.0,
                                        value=100.0, format="%.2f", key="u_ohv")
        else:
            oh_ew = st.number_input("OH eq. weight (g/eq)", min_value=0.01,
                                    value=561.0, format="%.2f", key="u_ohew")
            resin_ohv = 56100.0 / oh_ew if oh_ew else 0.0
            st.caption(f"= OHV {resin_ohv:.2f} mg KOH/g")
    with c2:
        st.markdown("**Isocyanate hardener**")
        hard_name = st.text_input("Name", "Polyiso hardener", key="u_hn")
        nco_src = st.radio("NCO spec", ["%NCO (by weight)",
                                        "NCO equivalent weight (g/eq)"],
                           key="u_ncosrc")
        if nco_src.startswith("%NCO"):
            hard_pct_nco = st.number_input("%NCO", min_value=0.0, value=20.0,
                                           format="%.2f", key="u_pnco")
        else:
            nco_ew = st.number_input("NCO eq. weight (g/eq)", min_value=0.01,
                                     value=210.0, format="%.2f", key="u_ncoew")
            hard_pct_nco = 4201.7 / nco_ew if nco_ew else 0.0
            st.caption(f"= {hard_pct_nco:.2f} %NCO")

    t1, t2 = st.columns(2)
    with t1:
        target_index = st.number_input(
            "Target NCO:OH index", min_value=1.0, max_value=300.0, value=105.0,
            format="%.1f", help="100 = 1:1 NCO:OH. >100 leaves residual NCO "
                                "(moisture-cure / crosslink density); <100 "
                                "leaves residual OH.")
    with t2:
        batch = st.number_input("Batch size (g, 0 = per 100 g resin)",
                                min_value=0.0, value=0.0, format="%.1f",
                                key="u_batch")

    if st.button("Solve mix ratio", type="primary"):
        res = solve_2k_index_ratio(resin_ohv, hard_pct_nco, target_index)
        if res.get("warning"):
            st.warning(res["warning"])
        else:
            st.session_state.u_result = {
                "res": res, "resin_name": resin_name, "hard_name": hard_name,
                "resin_ohv": resin_ohv, "hard_pct_nco": hard_pct_nco,
                "target_index": target_index, "batch": batch}
            st.session_state.changelog.append(
                f"[{date.today()}] 2K urethane: {resin_name} OHV "
                f"{resin_ohv:.1f} + {hard_name} {hard_pct_nco:.2f}%NCO at "
                f"index {target_index:.1f} -> {res['mix_ratio']:.4f} g "
                f"hardener per g resin ({res['wt_pct_hardener']:.1f} wt% "
                f"hardener)")

    R = st.session_state.get("u_result")
    if R:
        res = R["res"]
        st.subheader("Mix ratio")
        m1, m2, m3 = st.columns(3)
        m1.metric("Hardener : resin (by weight)",
                  f"{res['mix_ratio']:.4f} : 1")
        m2.metric("wt% hardener", f"{res['wt_pct_hardener']:.2f}%")
        m3.metric("Resulting index", f"{res['index']:.1f}")

        basis = R["batch"] if R["batch"] > 0 else 100.0
        scale = basis / res["total"]
        rg, hg = res["resin_parts"] * scale, res["hardener_parts"] * scale
        unit_note = (f"for a {R['batch']:.0f} g batch" if R["batch"] > 0
                     else "per 100 g resin basis")
        import pandas as pd
        st.dataframe(pd.DataFrame({
            "Component": [R["resin_name"], R["hard_name"], "TOTAL"],
            "Spec": [f"OHV {R['resin_ohv']:.1f}",
                     f"{R['hard_pct_nco']:.2f} %NCO", ""],
            "Grams": [round(rg, 2), round(hg, 2), round(rg + hg, 2)],
            "wt%": [round(100 * rg / (rg + hg), 2),
                    round(100 * hg / (rg + hg), 2), 100.0],
        }), use_container_width=True, hide_index=True)
        st.caption(f"Amounts shown {unit_note}. OH = {res['eq_oh']*scale:.4f} "
                   f"eq, NCO = {res['eq_nco']*scale:.4f} eq "
                   f"(index {res['index']:.1f}).")

    st.divider()
    st.caption("Stoichiometry only: models the primary NCO+OH reaction and "
               "ignores NCO-moisture (urea/CO2), NCO-amine, and "
               "allophanate/biuret branching. The index is a charge ratio, not "
               "a pot life or cure prediction. Confirm against a measured %NCO "
               "or a draw-down, not a yield model. Planning tool, not a CoA.")
