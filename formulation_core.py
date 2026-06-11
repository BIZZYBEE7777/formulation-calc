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
    t = {
        "g_asis": sum(r["g_asis"] for r in rows),
        "g_active": sum(r["g_active"] for r in rows),
        "eq_acid": sum(r["eq"] for r in rows if r["group"] in ("acid", "anhydride")),
        "eq_amine": sum(r["eq"] for r in rows if r["group"] == "amine"),
        "eq_oh": sum(r["eq"] for r in rows if r["group"] == "hydroxyl"),
        "eq_capper": sum(r["eq"] for r in rows if r["group"] == "capper"),
        # one water-free ring-opening per anhydride MOLECULE:
        "anh_moles": sum(r["moles"] for r in rows if r["group"] == "anhydride"),
    }
    return t


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


def condensation(rows, reaction, extent=1.0, cyc_extent=0.0):
    """Bonds + water of condensation + residual equivalents at conversion
    `extent` (fraction of the LIMITING side that reacts).

    Water-free bonds:
      - addition cappers (e.g. DCPD across a COOH): never release water
      - each anhydride molecule's FIRST ring-opening releases no water
    Imidazoline cyclization: `cyc_extent` (0-1) of eligible amide pairs
    release a SECOND water on ring closure. Eligible pairs are capped at
    the molar amount of amine components (one ring per polyamine molecule,
    DETA-style).
    """
    t = totals(rows)
    if reaction == "amidation":
        acid_side, nuc_side = t["eq_acid"], t["eq_amine"] + t["eq_capper"]
    elif reaction == "esterification":
        acid_side, nuc_side = t["eq_acid"], t["eq_oh"] + t["eq_capper"]
    else:
        return {"bonds": 0.0, "water_g": 0.0, "water_amide_g": 0.0,
                "water_cyc_g": 0.0, "cyc_rings": 0.0, "p_limiting": 0.0,
                "acid_side": t["eq_acid"], "nuc_side": 0.0,
                "residual": {"acid": t["eq_acid"], "amine": t["eq_amine"],
                             "hydroxyl": t["eq_oh"], "capper": t["eq_capper"]}}

    limiting = min(acid_side, nuc_side)
    bonds = limiting * extent

    # capper bonds happen first (no water), then anhydride first-openings
    capper_bonds = min(t["eq_capper"], bonds)
    anh_free = min(t["anh_moles"], bonds)
    water_free = min(bonds, capper_bonds + anh_free)
    # NOTE: assumes capper bonds and anhydride first-openings are distinct
    # bonds when possible; when a capper opens an anhydride this slightly
    # over-credits water-free bonds. Conservative for water COLLECTION
    # planning; flag in UI.
    water_g = (bonds - water_free) * WATER_MW

    # imidazoline cyclization: second water per ring, one ring max per
    # amine-component molecule
    amine_moles = sum(r["moles"] for r in rows if r["group"] == "amine")
    cyc_eligible = min(bonds, amine_moles)
    cyc_rings = cyc_eligible * max(0.0, min(1.0, cyc_extent))
    water_cyc_g = cyc_rings * WATER_MW
    water_total = water_g + water_cyc_g

    residual = {"acid": t["eq_acid"] - bonds,
                "amine": t["eq_amine"], "hydroxyl": t["eq_oh"],
                "capper": t["eq_capper"] - capper_bonds}
    nuc_bonds_remaining = bonds - capper_bonds
    if reaction == "amidation":
        residual["amine"] -= nuc_bonds_remaining
    else:
        residual["hydroxyl"] -= nuc_bonds_remaining

    return {"bonds": bonds, "water_g": water_total,
            "water_amide_g": water_g, "water_cyc_g": water_cyc_g,
            "cyc_rings": cyc_rings, "p_limiting": extent,
            "acid_side": acid_side, "nuc_side": nuc_side,
            "residual": residual}


def batch_summary(components, scale, reaction="none", extent=1.0, cyc_extent=0.0):
    """Everything the bench needs, at the solved scale."""
    rows = component_calcs(components, scale)
    t = totals(rows)
    cond = condensation(rows, reaction, extent, cyc_extent)

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
        ea, eb = cond["acid_side"], cond["nuc_side"]
        if ea > 0 and eb > 0:
            r_ratio = min(ea, eb) / max(ea, eb)
            p = extent
            denom = 1 + r_ratio - 2 * r_ratio * p
            xn = (1 + r_ratio) / denom if denom > 1e-12 else float("inf")
            react_groups = ("acid", "anhydride", partner, "capper")
            reacting = [x for x in rows if x["group"] in react_groups]
            n0 = sum(x["moles"] for x in reacting)
            mass_react = sum(x["g_active"] for x in reacting) - cond["water_g"]
            chains = n0 - cond["bonds"]
            mn = mass_react / chains if chains > 1e-12 else float("inf")
            excess_grp = "acid" if ea > eb else partner + "/capper"
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


# ---------------------------------------------------------------------------
# Inverse mode + reformulation helpers
# ---------------------------------------------------------------------------

def grams_to_ratios(components):
    """Inverse mode: each component dict carries 'g_asis_in' (grams as weighed).
    Sets 'ratio' = moles so the normal pipeline at scale=1.0 reproduces the
    batch exactly. Returns (components, normalized_ratio_list) where the
    normalized ratios are relative to the smallest nonzero molar amount."""
    for c in components:
        active = c["g_asis_in"] * (c.get("assay", 1.0) or 1.0)
        c["ratio"] = active / c["mw"] if c["mw"] else 0.0
    nz = [c["ratio"] for c in components if c["ratio"] > 0]
    base = min(nz) if nz else 1.0
    normalized = [c["ratio"] / base if base else 0.0 for c in components]
    return components, normalized


def weight_percents(rows):
    """As-is and active wt% for a computed rows list."""
    tot_asis = sum(r["g_asis"] for r in rows) or 1.0
    tot_act = sum(r["g_active"] for r in rows) or 1.0
    return [(100.0 * r["g_asis"] / tot_asis,
             100.0 * r["g_active"] / tot_act) for r in rows]


def adjust_to_yield(amounts, changed_idx, new_amount, absorber_idxs):
    """Reformulation: change one ingredient's as-is grams while holding total
    batch mass constant. The difference is absorbed by `absorber_idxs`
    proportionally to their current amounts.

    Returns (new_amounts, warning_or_None)."""
    amounts = list(amounts)
    delta = new_amount - amounts[changed_idx]
    absorbers = [i for i in absorber_idxs if i != changed_idx]
    pool = sum(amounts[i] for i in absorbers)
    if not absorbers or pool <= 0:
        return amounts, "No valid absorber ingredients selected."
    warning = None
    new_amounts = list(amounts)
    new_amounts[changed_idx] = new_amount
    for i in absorbers:
        share = amounts[i] / pool
        new_amounts[i] = amounts[i] - delta * share
        if new_amounts[i] < 0:
            warning = (f"Adjustment drives an absorber below zero grams "
                       f"(index {i}). Reduce the change or add absorbers.")
            new_amounts[i] = 0.0
    return new_amounts, warning


def cold_blend_solids(rows):
    """Cold blend (no reaction): solids = assay-corrected active grams of all
    non-inert components; inerts are the diluent/volatile side.
    Returns (solids_g, total_g, solids_pct)."""
    solids = sum(r["g_active"] for r in rows if r["group"] != "inert")
    total = sum(r["g_asis"] for r in rows)
    return solids, total, (100.0 * solids / total if total else 0.0)


def dilute_to_solids(amounts, rows, diluent_idx, target_pct):
    """Adjust ONE inert diluent's as-is grams so the blend hits target %
    solids. Solids mass is held constant; only the diluent changes.
    Returns (new_amounts, info_or_warning)."""
    if rows[diluent_idx]["group"] != "inert":
        return list(amounts), "Diluent must be a component marked 'inert'."
    if not (0 < target_pct <= 100):
        return list(amounts), "Target % solids must be between 0 and 100."
    solids = sum(r["g_active"] for r in rows if r["group"] != "inert")
    if solids <= 0:
        return list(amounts), "No non-inert solids in the formula."
    others = sum(a for i, a in enumerate(amounts) if i != diluent_idx)
    new_total = solids / (target_pct / 100.0)
    new_diluent = new_total - others
    if new_diluent < 0:
        max_pct = 100.0 * solids / others if others else 100.0
        return list(amounts), (f"Target unreachable by dilution alone: even "
                               f"at zero diluent the blend is {max_pct:.1f}% "
                               f"solids. Pick a lower target or reduce other "
                               f"components.")
    new_amounts = list(amounts)
    new_amounts[diluent_idx] = new_diluent
    return new_amounts, None


# ---------------------------------------------------------------------------
# Cold blend mode: weight-parts + per-ingredient % solids (no MW/chemistry)
# ---------------------------------------------------------------------------

US_FLOZ_ML = 29.5735
BLEND_PRESETS = [("4 oz", 4 * US_FLOZ_ML), ("8 oz", 8 * US_FLOZ_ML),
                 ("Pint (16 oz)", 16 * US_FLOZ_ML),
                 ("Quart (32 oz)", 32 * US_FLOZ_ML),
                 ("Gallon (128 oz)", 128 * US_FLOZ_ML),
                 ("100 g", None), ("1 kg", None), ("1 L", 1000.0)]


def blend_summary(parts, solids_pcts):
    """parts = weight parts; solids_pcts = % solids of each ingredient.
    Returns (wt_fractions, blend_solids_pct, solids_per_100g)."""
    total = sum(parts)
    if total <= 0:
        return [0.0] * len(parts), 0.0, 0.0
    wt = [p / total for p in parts]
    blend_pct = sum(w * s for w, s in zip(wt, solids_pcts))
    return wt, blend_pct, blend_pct  # per 100 g, grams of solids = pct value


def blend_solve_diluent(parts, solids_pcts, diluent_idx, target_pct):
    """Solve the diluent's parts so the blend hits target % solids.
    Works for any diluent solids content (water=0, a 30% dispersion, etc.).
    Returns (new_parts, warning_or_None)."""
    s_d = solids_pcts[diluent_idx]
    P_other = sum(p for i, p in enumerate(parts) if i != diluent_idx)
    S_other = sum(p * s for i, (p, s) in enumerate(zip(parts, solids_pcts))
                  if i != diluent_idx)
    t = target_pct
    if abs(s_d - t) < 1e-12:
        return list(parts), ("Diluent solids equals the target — its amount "
                             "can't move the blend toward that target.")
    x = (t * P_other - S_other) / (s_d - t)
    if x < 0:
        cur = S_other / P_other if P_other else 0.0
        direction = "below" if t < min(cur, s_d) else "above"
        return list(parts), (f"Target unreachable with this diluent: at zero "
                             f"{'' } diluent the blend is {cur:.1f}% solids "
                             f"and the diluent is {s_d:.0f}% — {t:.1f}% lies "
                             f"{direction} the reachable range.")
    new_parts = list(parts)
    new_parts[diluent_idx] = x
    return new_parts, None


def blend_batch_grams(parts, total_grams):
    """Scale weight parts to grams for a target batch mass."""
    total = sum(parts)
    return [total_grams * p / total for p in parts] if total else parts


def blend_sheet_text(names, solids_pcts, grams, batch_label, blend_pct, meta):
    L = []
    ap = L.append
    ap("=" * 64)
    ap(f"COLD BLEND CHARGE SHEET — {batch_label}")
    ap("=" * 64)
    ap(f"Project: {meta.get('project','____')}   Batch ID: "
       f"{meta.get('exp_id','____')}   By: {meta.get('chemist','____')}   "
       f"Date: {meta.get('date','')}")
    ap("-" * 64)
    ap(f"{'#':<3}{'Ingredient':<28}{'% solids':>9}{'grams':>10}{'  [ ] added'}")
    for i, (n, s, g) in enumerate(zip(names, solids_pcts, grams), 1):
        ap(f"{i:<3}{n[:27]:<28}{s:>8.1f}%{g:>10.2f}      [ ]")
    ap("-" * 64)
    ap(f"Total: {sum(grams):.2f} g    Blend solids: {blend_pct:.2f}%    "
       f"Solids mass: {sum(grams)*blend_pct/100:.2f} g")
    ap("Notes:")
    ap("\n\n")
    return "\n".join(L)


def solve_blend(specs, target_pct):
    """Constraint solver for cold blends, per 100 g of final product.

    specs: list of dicts {name, solids_pct, role, value} where role is:
      'fixed'   -- ingredient is a fixed wt% of the final blend (value = wt%)
      'ratio'   -- solids-carrier group held at relative weight ratios
                   (value = ratio number, e.g. 2 and 1 for a 2:1 pair)
      'balance' -- exactly ONE ingredient that fills to 100% (value ignored)
    target_pct: target % solids of the final blend.

    Returns (parts_per_100g_list, warning_or_None).
    """
    n = len(specs)
    parts = [0.0] * n
    fixed_idx = [i for i, s in enumerate(specs) if s["role"] == "fixed"]
    ratio_idx = [i for i, s in enumerate(specs) if s["role"] == "ratio"]
    bal_idx = [i for i, s in enumerate(specs) if s["role"] == "balance"]

    if len(bal_idx) != 1:
        return parts, ("Mark exactly ONE ingredient as 'balance' (usually "
                       "water) — it fills the blend to 100%.")
    b = bal_idx[0]
    s_w = specs[b]["solids_pct"] / 100.0

    F = sum(specs[i]["value"] for i in fixed_idx)
    Sf = sum(specs[i]["value"] * specs[i]["solids_pct"] / 100.0
             for i in fixed_idx)
    R_m = sum(specs[i]["value"] for i in ratio_idx)
    R_s = sum(specs[i]["value"] * specs[i]["solids_pct"] / 100.0
              for i in ratio_idx)

    if F > 100:
        return parts, "Fixed wt% ingredients alone exceed 100%."

    T = target_pct  # g solids per 100 g
    if ratio_idx:
        denom = R_s - R_m * s_w
        if abs(denom) < 1e-12:
            return parts, ("Solids-ratio group and balance ingredient have "
                           "the same effective solids — target can't be set "
                           "by trading between them.")
        k = (T - Sf - (100.0 - F) * s_w) / denom
        if k < 0:
            return parts, (f"Unreachable: hitting {T:.2f}% solids would need "
                           f"a negative amount of the ratio group. Check the "
                           f"target vs. fixed ingredients' solids.")
        for i in ratio_idx:
            parts[i] = k * specs[i]["value"]
    else:
        k = 0.0

    w = 100.0 - F - k * R_m
    if w < -1e-9:
        return parts, (f"Unreachable: fixed + solids-carrier masses exceed "
                       f"100 g per 100 g of blend (over by {-w:.2f} g). "
                       f"Lower the target % solids or a fixed wt%.")
    parts[b] = max(w, 0.0)
    for i in fixed_idx:
        parts[i] = specs[i]["value"]

    # verify
    tot = sum(parts)
    solids = sum(p * s["solids_pct"] / 100.0 for p, s in zip(parts, specs))
    if abs(tot - 100.0) > 1e-6 or abs(solids - T) > 1e-6:
        return parts, (f"Solver sanity check failed (mass {tot:.3f}, solids "
                       f"{solids:.3f}) — please report this combination.")
    return parts, None


# ---------------------------------------------------------------------------
# Cook tracker + spec-targeting solvers (1-D bisection on monotone responses)
# ---------------------------------------------------------------------------

def _value_at(components, reaction, key, extent, cyc_extent):
    out = batch_summary(components, 1.0, reaction=reaction, extent=extent,
                        cyc_extent=cyc_extent)
    return out["end_values"].get(key, 0.0), out


def p_from_measured(components, reaction, key, measured, cyc_extent=0.0):
    """Cook tracker: invert a measured end value (key in 'acid_value',
    'amine_value', 'hydroxyl_value') to the implied conversion p.
    Returns (p, summary_at_p, warning_or_None)."""
    lo, hi = 0.0, 1.0
    v_lo, _ = _value_at(components, reaction, key, lo, cyc_extent)
    v_hi, _ = _value_at(components, reaction, key, hi, cyc_extent)
    # end values fall as p rises
    if not (min(v_hi, v_lo) - 1e-9 <= measured <= max(v_hi, v_lo) + 1e-9):
        return None, None, (f"Measured {measured:.1f} is outside the "
                            f"theoretical range for this charge "
                            f"({min(v_lo, v_hi):.1f} at p=1 to "
                            f"{max(v_lo, v_hi):.1f} at p=0). Check the value, "
                            f"the charge entries, or sample dilution.")
    for _ in range(80):
        mid = (lo + hi) / 2.0
        v, _ = _value_at(components, reaction, key, mid, cyc_extent)
        if (v - measured) * (v_lo - measured) > 0:
            lo, v_lo = mid, v
        else:
            hi = mid
    p = (lo + hi) / 2.0
    _, summary = _value_at(components, reaction, key, p, cyc_extent)
    return p, summary, None


def solve_ratio_for_target(components, vary_name, reaction, key, target,
                           extent=1.0, cyc_extent=0.0, max_ratio=1000.0):
    """Spec targeting: solve ONE component's molar ratio so the end value
    `key` hits `target` at conversion `extent`. Other ratios held fixed.
    Returns (ratio, summary, warning_or_None)."""
    idx = next((i for i, c in enumerate(components)
                if c["name"] == vary_name), None)
    if idx is None:
        return None, None, "Component not found."

    def val(r):
        comps = [dict(c) for c in components]
        comps[idx]["ratio"] = r
        out = batch_summary(comps, 1.0, reaction=reaction, extent=extent,
                            cyc_extent=cyc_extent)
        return out["end_values"].get(key, 0.0), out

    lo, hi = 1e-6, 1.0
    v_lo, _ = val(lo)
    v_hi, _ = val(hi)
    # expand upper bracket until target enclosed or limit hit
    while (v_lo - target) * (v_hi - target) > 0 and hi < max_ratio:
        hi *= 2.0
        v_hi, _ = val(hi)
    if (v_lo - target) * (v_hi - target) > 0:
        return None, None, (f"No amount of '{vary_name}' between ~0 and "
                            f"{max_ratio:g} (molar ratio) reaches "
                            f"{key.replace('_', ' ')} = {target:.1f} at "
                            f"p = {extent:.3f}. Vary a different component "
                            f"or revisit the target.")
    for _ in range(90):
        mid = (lo + hi) / 2.0
        v, _ = val(mid)
        if (v - target) * (v_lo - target) > 0:
            lo, v_lo = mid, v
        else:
            hi = mid
    r = (lo + hi) / 2.0
    _, summary = val(r)
    return r, summary, None
