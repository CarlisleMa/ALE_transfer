"""Verifier for the high-content-screening dose-response family.

Truth is planted by the generator (reference/truth.json):
  * wells[well].n_cells    - number of cells the simulator rendered into the field
  * wells[well].qc_fail    - 1 if the field was rendered out of focus
  * compounds[c].marker_ic50 / count_gi50 - concentration at which the planted response curve of
    that readout equals 0.5 of vehicle, or null when the planted curve does not reach 0.5 within
    the tested range (no effect, plateau above 0.5, potency beyond the top concentration)
No segmentation result and no curve fit is used as reference, so tolerances only have to cover
segmentation error of correct methods and the precision of an 8-point curve.
"""
import math
import re
import sys

sys.path.insert(0, "/tests")
from _scoring import CandidateError, ramp, run, to_float  # noqa: E402

NA = {"", "na", "nan", "n/a", "none", "null"}


def well_id(s):
    m = re.fullmatch(r"\s*([A-Za-z])0*(\d{1,2})\s*", str(s))
    return f"{m.group(1).upper()}{int(m.group(2)):02d}" if m else None


def potency(s, what):
    s = str(s if s is not None else "").strip()
    if s.lower() == "inactive":
        return None
    v = to_float(s, what)
    if v <= 0:
        raise CandidateError(f"{what} must be a positive concentration or the literal 'inactive'")
    return v


def mcc(tp, fp, tn, fn):
    den = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return (tp * tn - fp * fn) / den if den else 0.0


def main(ctx):
    truth = ctx.ref_json("truth.json")
    wells = truth["wells"]
    comps = truth["compounds"]
    tol = truth["tolerance"]
    controls = {c.lower() for c in truth["control_names"]}

    # ---------------- well_measurements.tsv --------------------------------------------
    rows = ctx.table("well_measurements.tsv", ["well", "n_cells", "marker_per_cell", "qc_pass"],
                     max_rows=len(wells) + 10)
    got = {}
    for r in rows:
        w = well_id(r["well"])
        ctx.gate(w in wells, f"well_measurements.tsv: well {str(r['well'])[:12]!r} is on the plate layout")
        ctx.gate(w not in got, f"well_measurements.tsv: well {w} listed once")
        n = to_float(r["n_cells"], f"n_cells of {w}")
        if n < 0 or abs(n - round(n)) > 1e-6:
            raise CandidateError("gate failed: well_measurements.tsv: n_cells is a non-negative integer")
        q = to_float(r["qc_pass"], f"qc_pass of {w}")
        if q not in (0.0, 1.0):
            raise CandidateError("gate failed: well_measurements.tsv: qc_pass is 0 or 1")
        m = str(r["marker_per_cell"] or "").strip()
        if m.lower() not in NA:
            to_float(m, f"marker_per_cell of {w}")       # must be a finite number when given
        got[w] = (int(round(n)), int(q))
    ctx.gate(len(got) == len(wells), "well_measurements.tsv lists every well of the layout exactly once")

    usable = [w for w in wells if not wells[w]["qc_fail"]]
    good, rel = 0, []
    for w in usable:
        t, n = wells[w]["n_cells"], got[w][0]
        good += abs(n - t) <= max(tol["count_abs"], tol["count_rel"] * t)
        rel.append((n - t) / max(t, 1))
    frac = good / len(usable)
    rel.sort()
    ctx.info["count_fraction_within_tolerance"] = round(frac, 4)
    ctx.info["count_median_rel_error"] = round(rel[len(rel) // 2], 4)
    counts_value = ctx.component(
        "cell_counts", ramp(frac, 0.5, 0.9), 0.20,
        "fraction of in-focus wells with |n_cells - planted| <= max(abs, rel) tolerance; full credit >= 0.90")

    failed = {w for w in wells if wells[w]["qc_fail"]}
    flagged = {w for w in wells if got[w][1] == 0}
    tp, fp = len(flagged & failed), len(flagged - failed)
    fn, tn = len(failed - flagged), len(wells) - len(flagged | failed)
    ctx.info.update({"qc_tp": tp, "qc_fp": fp, "qc_fn": fn})
    if failed:
        m = mcc(tp, fp, tn, fn)
        ctx.info["qc_mcc"] = round(m, 4)
        ctx.component("qc_flags", ramp(m, 0.5, 0.9), 0.10,
                      "Matthews correlation of qc_pass = 0 with the wells rendered out of focus; full credit >= 0.90")
    else:
        # nothing to flag on this plate: credit for not flagging, but only as far as the wells were measured
        ctx.component("qc_flags", (1.0 - min(1.0, fp / 5.0)) * counts_value, 0.10,
                      "no well is out of focus on this plate: (1 - false flags / 5) x cell_counts credit")

    # ---------------- dose_response.tsv ------------------------------------------------
    drows = ctx.table("dose_response.tsv", ["compound", "ic50_marker_uM", "gi50_count_uM"], max_rows=200)
    names = {c.lower(): c for c in comps}
    dr = {}
    for r in drows:
        key = str(r["compound"]).strip().lower()
        if key in controls:
            continue                                  # rows for vehicle / reference are ignored
        ctx.gate(key in names, f"dose_response.tsv: compound {str(r['compound'])[:20]!r} is on the plate layout")
        c = names[key]
        ctx.gate(c not in dr, f"dose_response.tsv: compound {c} listed once")
        dr[c] = (potency(r["ic50_marker_uM"], f"ic50_marker_uM of {c}"),
                 potency(r["gi50_count_uM"], f"gi50_count_uM of {c}"))
    ctx.gate(len(dr) == len(comps), "dose_response.tsv lists every test compound exactly once")

    for idx, (key, comp_name, weight) in enumerate((("marker_ic50", "marker_potency", 0.45),
                                                    ("count_gi50", "count_potency", 0.25))):
        credits, worst, wrong_class = [], 0.0, []
        for c in sorted(comps):
            t, v = comps[c][key], dr[c][idx]
            if t is None or v is None:
                cr = 1.0 if (t is None and v is None) else 0.0
                if cr == 0.0:
                    wrong_class.append(c)
            else:
                err = abs(math.log10(v / t))
                worst = max(worst, err)
                cr = ramp(err, tol["log10_zero"], tol["log10_full"])
            credits.append(cr)
        mean = sum(credits) / len(credits)
        ctx.info[f"{comp_name}_mean_credit"] = round(mean, 4)
        ctx.info[f"{comp_name}_max_abs_log10_error"] = round(worst, 4)
        ctx.info[f"{comp_name}_wrong_activity_class"] = wrong_class
        ctx.component(comp_name, ramp(mean, 0.8, 1.0), weight,
                      "per compound: 'inactive' must match; active: |log10(est/planted)| <= full tolerance "
                      "gives 1, >= zero tolerance gives 0; component = ramp(mean credit, 0.8 -> 1.0)")


run(main, pass_threshold=0.9)
