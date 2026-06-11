"""
formulation_core.py -- pure chemistry/math logic for the formulation calculator.
No UI, no network. Fully unit-testable.

Component dict fields:
  name        str   display name
  cas         str   optional CAS number
  mw          float molecular weight (g/mol) of the ACTIVE substance
  assay       float assay / purity as fraction (0-1]; as-is grams = active/assay
  functionality float reactive groups per molecule (e.g. diacid = 2)
  group       str   'acid' | 'amine' | 'hydroxyl' | 'inert'
  ratio       float molar ratio (relative moles)

Reaction types supported for condensation accounting:
  'amidation'      acid + amine  -> amide + H2O
  'esterification' acid + hydroxyl -> ester + H2O
  'none'           no condensation math
"""

WATER_MW = 18.015
KOH_MW_MG = 56100.0  # mg KOH per equivalent, for AV/AmV/OHV


def component_calcs(components, scale=1.0):
    """Per-component moles/equivalents/grams at a given scale factor."""
    rows = []
    for c in components:
        moles = c["ratio"] * scale
        eq = moles * c["functionality"]
        g_active = moles * c["mw"]
        assay = c.get("assay", 1.0) or 1.0
        g_asis = g_active / assay
        rows.append({**c, "moles": moles, "eq": eq,
                     "g_active": g_active, "g_asis": g_asis})
    return rows


def totals(rows):
    return {
        "g_asis": sum(r["g_asis"] for r in rows),
        "g_active": sum(r["g_active"] for r in rows),
        "eq_acid": sum(r["eq"] for r in rows if r["group"] == "acid"),
        "eq_amine": sum(r["eq"] for r in rows if r["group"] == "amine"),
        "eq_oh": sum(r["eq"] for r in rows if r["group"] == "hydroxyl"),
    }


def solve_scale(components, anchor_type, anchor_value, anchor_component=None):
    """Find the scale factor so the batch hits the anchor.
    anchor_type: 'total_asis_mass' | 'component_asis_mass' | 'component_moles'
    """
    base = component_calcs(components, scale=1.0)
    if anchor_type == "total_asis_mass":
        per_unit = sum(r["g_asis"] for r in base)
        return anchor_value / per_unit if per_unit else 0.0
    idx = next(i for i, c in enumerate(components)
               if c["name"] == anchor_component)
    if anchor_type == "component_asis_mass":
        per_unit = base[idx]["g_asis"]
        return anchor_value / per_unit if per_unit else 0.0
    if anchor_type == "component_moles":
        per_unit = base[idx]["moles"]
        return anchor_value / per_unit if per_unit else 0.0
    raise ValueError(anchor_type)


def condensation(rows, reaction, extent=1.0):
    """Water of condensation + residual equivalents at conversion `extent`
    (extent = fraction of the LIMITING group that reacts)."""
    t = totals(rows)
    if reaction == "amidation":
        a, b = t["eq_acid"], t["eq_amine"]
    elif reaction == "esterification":
        a, b = t["eq_acid"], t["eq_oh"]
    else:
        return {"bonds": 0.0, "water_g": 0.0, "p_limiting": 0.0,
                "residual": {"acid": t["eq_acid"], "amine": t["eq_amine"],
                             "hydroxyl": t["eq_oh"]}}
    limiting = min(a, b)
    bonds = limiting * extent
    water_g = bonds * WATER_MW
    residual = {"acid": t["eq_acid"], "amine": t["eq_amine"],
                "hydroxyl": t["eq_oh"]}
    if reaction == "amidation":
        residual["acid"] -= bonds
        residual["amine"] -= bonds
    else:
        residual["acid"] -= bonds
        residual["hydroxyl"] -= bonds
    return {"bonds": bonds, "water_g": water_g, "p_limiting": extent,
            "residual": residual}


def batch_summary(components, scale, reaction="none", extent=1.0):
    """Everything the bench needs, at the solved scale."""
    rows = component_calcs(components, scale)
    t = totals(rows)
    cond = condensation(rows, reaction, extent)

    charge_mass = t["g_asis"]                       # what goes in the kettle
    inerts = t["g_asis"] - t["g_active"]            # solvent/water-of-dilution etc.
    final_mass = charge_mass - cond["water_g"]      # after condensate removed
    resin_mass = t["g_active"] - cond["water_g"]    # active solids produced
    solids_pct = 100.0 * resin_mass / final_mass if final_mass else 0.0

    # ---- end-group values (mg KOH / g of final resin solids) ----
    ev = {}
    if resin_mass > 0:
        r = cond["residual"]
        ev["acid_value"] = KOH_MW_MG * r["acid"] / resin_mass
        ev["amine_value"] = KOH_MW_MG * r["amine"] / resin_mass
        ev["hydroxyl_value"] = KOH_MW_MG * r["hydroxyl"] / resin_mass

    # ---- Carothers / step-growth stats (reacting monomers only) ----
    car = {}
    if reaction in ("amidation", "esterification"):
        partner = "amine" if reaction == "amidation" else "hydroxyl"
        ea = t["eq_acid"]
        eb = t["eq_amine"] if reaction == "amidation" else t["eq_oh"]
        if ea > 0 and eb > 0:
            r_ratio = min(ea, eb) / max(ea, eb)
            p = extent
            denom = 1 + r_ratio - 2 * r_ratio * p
            xn = (1 + r_ratio) / denom if denom > 1e-12 else float("inf")
            reacting = [x for x in rows if x["group"] in ("acid", partner)]
            n0 = sum(x["moles"] for x in reacting)
            mass_react = sum(x["g_active"] for x in reacting) - cond["water_g"]
            chains = n0 - cond["bonds"]
            mn = mass_react / chains if chains > 1e-12 else float("inf")
            excess_grp = "acid" if ea > eb else partner
            car = {"r": r_ratio, "Xn": xn, "Mn": mn, "excess_group": excess_grp,
                   "excess_pct": 100.0 * (max(ea, eb) / min(ea, eb) - 1.0)}

    return {"rows": rows, "totals": t, "cond": cond,
            "charge_mass": charge_mass, "inerts": inerts,
            "final_mass": final_mass, "resin_mass": resin_mass,
            "solids_pct": solids_pct, "end_values": ev, "carothers": car}


def charge_sheet_text(summary, meta):
    """Plain-text printable charge sheet with header block."""
    L = []
    ap = L.append
    ap("=" * 70)
    ap("BATCH CHARGE SHEET")
    ap("=" * 70)
    ap(f"Project: {meta.get('project','________')}    "
       f"Experiment ID: {meta.get('exp_id','________')}")
    ap(f"Chemist: {meta.get('chemist','________')}    "
       f"Date: {meta.get('date','________')}")
    ap(f"Description: {meta.get('description','')}")
    ap("-" * 70)
    ap(f"{'#':<3}{'Component':<28}{'CAS':<13}{'As-is g':>10}{'Active g':>10}"
       f"{'  [ ] added'}")
    for i, r in enumerate(summary["rows"], 1):
        ap(f"{i:<3}{r['name'][:27]:<28}{(r.get('cas') or ''):<13}"
           f"{r['g_asis']:>10.2f}{r['g_active']:>10.2f}      [ ]")
    ap("-" * 70)
    ap(f"Total charge (as-is): {summary['charge_mass']:.2f} g")
    if summary['cond']['water_g'] > 0:
        ap(f"Predicted condensate (H2O): {summary['cond']['water_g']:.2f} g  "
           f"-> collect & log actual: ________ g")
        ap(f"Theoretical final mass: {summary['final_mass']:.2f} g   "
           f"theoretical solids: {summary['solids_pct']:.1f}%")
    ev = summary["end_values"]
    if ev:
        ap(f"Theoretical AV: {ev['acid_value']:.1f}   "
           f"AmV: {ev['amine_value']:.1f}   OHV: {ev['hydroxyl_value']:.1f}"
           f"   (mg KOH/g resin)")
        ap("Titrated AV: ________   AmV: ________   OHV: ________")
    car = summary["carothers"]
    if car:
        ap(f"Stoich: r = {car['r']:.3f} ({car['excess_pct']:.1f}% excess "
           f"{car['excess_group']});  Xn = {car['Xn']:.1f};  "
           f"Mn ~ {car['Mn']:.0f} g/mol at stated conversion")
    ap("=" * 70)
    ap("Notes / observations:")
    ap("\n\n\n")
    return "\n".join(L)
