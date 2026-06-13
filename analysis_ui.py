"""Post-engine family analysis tools (the pattern going forward).

These are NOT additions to the shared step-growth engine — they are render
functions that take the engine's `out` summary and apply family-specific
science on top of it. Each states its assumptions/source in the UI and
distinguishes theory from empirical fit. `formulation_app.py` calls the
relevant ones per active family.

  render_oil_length(out)            — alkyd (charge bookkeeping)
  render_gel_point(out, extent)     — alkyd + UPR (Carothers theory)
  render_unsaturation(out)          — UPR (C=C bookkeeping)
  render_viscosity_ohv()            — alkyd/ester (EMPIRICAL fit to user data)
  render_spec_strip(out, basic_n)   — all families (live value vs targets)
"""
import pandas as pd
import streamlit as st

from formulation_core import (oil_length, gel_point_carothers,
                              unsaturation_stats, fit_log_viscosity_ohv,
                              predict_ohv_from_viscosity, total_amine_value)

# ---- General (all-family) live spec strip: value vs user target range ----
# (key, short label, value format, unit)
_SPEC_DEFS = [
    ("acid_value", "AV", "{:.1f}", "mg KOH/g"),
    ("hydroxyl_value", "OHV", "{:.1f}", "mg KOH/g"),
    ("amine_value", "AmV", "{:.1f}", "mg KOH/g"),
    ("total_amine_value", "TAV", "{:.1f}", "mg KOH/g"),
    ("pct_nco", "%NCO", "{:.2f}", "%"),
]
_LABEL2KEY = {lab: k for k, lab, _, _ in _SPEC_DEFS}
DEFAULT_SPEC_TARGETS = pd.DataFrame(
    [{"Spec": lab, "Track": False, "Low": 0.0, "High": 0.0}
     for _, lab, _, _ in _SPEC_DEFS])


def _spec_status(v, lo, hi):
    """🟢 in [lo,hi], 🟡 within a 10%-of-range margin outside, 🔴 beyond,
    ⚪ if no value."""
    if v is None:
        return "⚪"
    if hi < lo:
        lo, hi = hi, lo
    if lo <= v <= hi:
        return "🟢"
    width = hi - lo
    margin = 0.10 * width if width > 0 else 0.05 * max(abs(hi), 1.0)
    return "🟡" if (lo - margin) <= v <= (hi + margin) else "🔴"


def render_spec_strip(out, basic_n):
    """Live AV/OHV/AmV/TAV/%NCO vs user-set target ranges, green/amber/red.
    Targets persist in the formula JSON (pin button). All-family — a planning
    aid, not a CoA."""
    ev = out.get("end_values", {})
    vals = {"acid_value": ev.get("acid_value"),
            "hydroxyl_value": ev.get("hydroxyl_value"),
            "amine_value": ev.get("amine_value"),
            "total_amine_value": total_amine_value(out, basic_n),
            "pct_nco": ev.get("pct_nco")}
    if "spec_targets" not in st.session_state:
        st.session_state.spec_targets = DEFAULT_SPEC_TARGETS.copy()
    if "spec_tgt_ver" not in st.session_state:
        st.session_state.spec_tgt_ver = 0

    with st.expander("🎯 Spec targets — tick a spec and set its range"):
        ed = st.data_editor(
            st.session_state.spec_targets, width="stretch", hide_index=True,
            column_config={
                "Spec": st.column_config.TextColumn(disabled=True),
                "Track": st.column_config.CheckboxColumn(
                    help="Show this spec in the status strip"),
                "Low": st.column_config.NumberColumn(format="%.2f"),
                "High": st.column_config.NumberColumn(format="%.2f"),
            }, key=f"spec_tgt_v{st.session_state.spec_tgt_ver}")
        if st.button("📌 Pin targets to formula", key="spec_pin",
                     help="Save these ranges in the formula JSON. The strip "
                          "already uses the table live; pin before you Save."):
            st.session_state.spec_targets = ed.copy()
            st.session_state.spec_tgt_ver += 1
            st.success("Spec targets pinned — they'll save with the formula.")
            st.rerun()

    tracked = ed[ed["Track"] == True]  # noqa: E712
    if tracked.empty:
        return
    fmt = {lab: (f, u) for _, lab, f, u in _SPEC_DEFS}
    cols = st.columns(len(tracked))
    for col, (_, row) in zip(cols, tracked.iterrows()):
        lab = str(row["Spec"])
        lo, hi = float(row["Low"]), float(row["High"])
        v = vals.get(_LABEL2KEY.get(lab))
        f, u = fmt.get(lab, ("{:.1f}", ""))
        vtxt = f.format(v) if v is not None else "—"
        col.markdown(f"{_spec_status(v, lo, hi)} **{lab}**  \n{vtxt} {u}  \n"
                     f"<small>target {lo:g}–{hi:g}</small>",
                     unsafe_allow_html=True)

DEFAULT_VISC = pd.DataFrame({
    "OHV": pd.Series([None, None, None], dtype="float"),
    "Viscosity": pd.Series([None, None, None], dtype="float"),
    "Temp (°C)": pd.Series([None, None, None], dtype="float"),
})


def render_oil_length(out):
    st.markdown("**Oil length**")
    oil_names = st.multiselect(
        "Tag the oil / fatty-acid component(s)",
        [r["name"] for r in out["rows"]], key="oil_tag",
        help="Oil length = active wt of the tagged rows ÷ theoretical resin "
             "solids × 100. NB: classic oil length is on a triglyceride-oil "
             "basis; this is the charged-active-weight ratio of whatever you "
             "tag.")
    if oil_names:
        ol = oil_length(out, oil_names)
        bucket = "short" if ol < 45 else "medium" if ol < 60 else "long"
        st.metric("Oil length", f"{ol:.1f}%",
                  help="<45 short · 45–60 medium · >60 long oil")
        st.caption(f"Tagged active ÷ resin solids → {bucket} oil "
                   f"(by charged active weight).")


def render_gel_point(out, extent):
    g = gel_point_carothers(out)
    if g["p_gel"] is not None:
        st.metric("Carothers gel point", f"p_gel ≈ {g['p_gel']:.3f}",
                  delta=f"f_avg {g['f_avg']:.2f}", delta_color="off",
                  help="Branched step-growth gels at p ≈ 2/f_avg (Carothers; "
                       "over-predicts vs the Flory statistical treatment). "
                       "Assumes equal reactivity & balanced stoichiometry.")
        if extent >= g["p_gel"]:
            st.warning(f"⚠️ Conversion p = {extent:.3f} is at or past the "
                       f"Carothers gel point ({g['p_gel']:.3f}) — expect "
                       f"gelation.")
    else:
        st.caption(f"Average functionality f_avg = {g['f_avg']:.2f} ≤ 2 → "
                   f"linear, no Carothers gel point.")


def render_unsaturation(out):
    u = unsaturation_stats(out)
    if u["cc_moles"] > 0:
        uu1, uu2 = st.columns(2)
        eqw = ("∞" if u["cc_eq_weight"] == float("inf")
               else f"{u['cc_eq_weight']:.0f} g/mol")
        uu1.metric("C=C equivalent weight", eqw,
                   help="Resin solids per mole of reactive C=C. Lower = more "
                        "unsaturation = higher crosslink-density potential.")
        uu2.metric("Unsaturation", f"{u['mmol_per_g']:.2f} mmol C=C/g",
                   help="= mol C=C per kg resin. A backbone count for peroxide "
                        "/ UV radical cure — NOT a cure, gel time or exotherm "
                        "prediction (those are kinetic / measured).")
    else:
        st.caption("Set the **C=C / mol** column on the unsaturated monomer(s) "
                   "(maleic/fumaric ≈ 1, acrylate per double bond) to see C=C "
                   "equivalent weight and unsaturation per kg.")


def render_viscosity_ohv():
    """EMPIRICAL viscosity↔OHV calibration on the user's own cook data."""
    st.markdown("**Viscosity → OHV (empirical calibration)**")
    st.caption("Fits ln(viscosity) = a·OHV + b to YOUR cook history at one "
               "fixed temperature, then predicts OHV from a new cap-viscosity "
               "reading. **Empirical fit — your data, your resin; not "
               "transferable to another formula or temperature without "
               "recalibrating.**")

    if "visc_cal" not in st.session_state:
        st.session_state.visc_cal = DEFAULT_VISC.copy()
    if "visc_cal_ver" not in st.session_state:
        st.session_state.visc_cal_ver = 0

    ed = st.data_editor(
        st.session_state.visc_cal, num_rows="dynamic", width="stretch",
        column_config={
            "OHV": st.column_config.NumberColumn(
                "OHV (mg KOH/g)", min_value=0.0, format="%.1f"),
            "Viscosity": st.column_config.NumberColumn(
                "Viscosity", min_value=0.0, format="%.3f",
                help="Any consistent unit (poise, cP, bubble seconds) at a "
                     "fixed temperature."),
            "Temp (°C)": st.column_config.NumberColumn(
                format="%.0f", help="Stored per point; the fit uses ONE "
                                    "temperature (no cross-temp prediction)."),
        }, key=f"visc_cal_v{st.session_state.visc_cal_ver}")

    if st.button("📌 Pin points to formula", key="visc_pin",
                 help="Save these calibration points into the formula JSON. "
                      "The fit and prediction below already use the table "
                      "live; pinning is only needed before you Save."):
        st.session_state.visc_cal = ed.copy()
        st.session_state.visc_cal_ver += 1
        st.success("Calibration points pinned — they'll save with the formula.")
        st.rerun()

    valid = ed.dropna(subset=["OHV", "Viscosity"])
    valid = valid[valid["Viscosity"] > 0]
    temps = sorted({float(t) for t in valid["Temp (°C)"].dropna()})
    if len(temps) > 1:
        pick = st.selectbox(
            "Temperature to fit (no cross-temp prediction)", temps,
            format_func=lambda t: f"{t:.0f} °C", key="visc_temp")
        fitrows = valid[valid["Temp (°C)"] == pick]
        st.caption(f"Fitting {len(fitrows)} point(s) at {pick:.0f} °C; other "
                   f"temperatures excluded (viscosity is temperature-"
                   f"dependent).")
    else:
        fitrows = valid
        if temps:
            st.caption(f"All points at {temps[0]:.0f} °C.")

    pts = list(zip(fitrows["OHV"].astype(float),
                   fitrows["Viscosity"].astype(float)))
    fit = fit_log_viscosity_ohv(pts)
    if "a" not in fit:
        st.info(fit["warning"])
        return
    if fit["warning"]:
        st.warning(fit["warning"])

    c1, c2 = st.columns([1, 3])
    c1.metric("Fit R²", f"{fit['r2']:.4f}")
    c2.caption(f"ln(viscosity) = {fit['a']:.5g}·OHV + {fit['b']:.4g}  "
               f"(n = {fit['n']})")

    p1, p2 = st.columns(2)
    with p1:
        meas = st.number_input("Measured viscosity → predict OHV",
                               min_value=0.0, value=0.0, format="%.3f",
                               key="visc_meas",
                               help="Same unit & temperature as the "
                                    "calibration above.")
    with p2:
        conf = st.radio("Confidence", [0.95, 0.90],
                        format_func=lambda c: f"{int(c * 100)}%",
                        horizontal=True, key="visc_conf")
    if meas > 0:
        pr = predict_ohv_from_viscosity(fit, meas, conf)
        if "ohv" in pr:
            st.metric("Predicted OHV", f"{pr['ohv']:.1f} mg KOH/g",
                      delta=f"± {pr['half_width']:.1f} ({int(conf * 100)}% PI)",
                      delta_color="off")
            vmin = float(fitrows["Viscosity"].min())
            vmax = float(fitrows["Viscosity"].max())
            if meas < vmin or meas > vmax:
                st.warning(f"⚠️ {meas:g} is outside the calibrated viscosity "
                           f"range ({vmin:g}–{vmax:g}) — this is extrapolation; "
                           f"the predicted OHV and interval are unreliable.")
            st.caption("Empirical fit on your data — not a CoA, not "
                       "transferable to another resin or temperature. Confirm "
                       "against a titrated OHV.")
        else:
            st.warning(pr["warning"])
