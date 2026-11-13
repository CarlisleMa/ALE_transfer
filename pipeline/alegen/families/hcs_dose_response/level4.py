"""Level 4 of the hcs_dose_response family: the level-3 plate plus

  Q_validated_counts   the images are hard enough that the counting method matters (overlapping nuclei
                       without an intensity dip, chromatin texture, bright nuclear fragments the SOP does
                       not count, nuclei cut by the field border under a counting-frame rule, nuclei that
                       grow under treatment). The lab's manual point annotation of eight fields is among
                       the inputs and the SOP states the acceptance criterion on them.
  P_potency_shape      two curve shapes for which the textbook four-parameter fit with free plateaus
                       reports a wrong potency: a deep partial count response (the SOP's potency is the
                       absolute 0.5 crossing, far from the inflection point) and a marker signal per cell
                       that rises as cells are lost (not an inhibition at all).

Truth by construction as at levels 1-3: sim4.py plants every nucleus, fragment and curve before a pixel
exists; the reference count of a field is the number of rendered nuclei that the SOP counts.
Levels 1-3 never import this module, so their tasks stay byte-identical.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from ...spec import Attack, Probe, TaskBuild, Trap
from . import image, sim, sim4

ASSETS = Path(__file__).parent / "assets"
TOL = {"count_abs": 3, "count_rel": 0.10, "count_frac_full": 0.85, "count_frac_zero": 0.65,
       "log10_full": 0.20, "log10_zero": 0.40}
CONTROL_NAMES = ["DMSO", "REF-INH"]

INSTRUCTION = """\
You are the image-analysis scientist of a compound-profiling group. A 96-well plate of cells was
treated with eight test compounds in dose series, fixed, stained for nuclei and for a cytoplasmic
marker protein, and imaged on an automated microscope. The project team needs per-well
measurements and, for every test compound, its potency on the marker and on cell number. The
group's SOP requires that the cell-counting procedure is checked against the lab's manual counts
of this plate before it is used.

Task directory:
- `/task`

Inputs (read-only) under `/task/input`:
- `images/<well>_ch1.tif`, `images/<well>_ch2.tif` - one field per well ({n_wells} wells), 16-bit
  grayscale TIFF, {size} x {size} pixels. Channel 1 = nuclear stain, channel 2 = marker.
- `flatfield/flat_ch1.tif`, `flatfield/flat_ch2.tif`, `flatfield/dark.tif` - reference frames
  recorded with the same optical path (uniform fluorescent slide per channel; closed shutter).
- `plate_layout.tsv` - columns `well`, `compound`, `dose_uM`, `role` (`compound`, `vehicle` or
  `reference`).
- `annotations/<well>.tsv` - the lab's manual point annotation of {n_ann} fields of this plate, one
  file per annotated well: tab-separated, header `x_px`, `y_px`, one row per cell that was counted
  in that field under SOP section 1 (a click on its nucleus). Pixel coordinates in the image of
  that well: `x_px` = column, `y_px` = row, 0-based, origin at the top-left pixel. The manual count
  of a field is the number of rows.
- `assay_notes.md` - hand-over notes of the wet lab and of the imaging operator.
- `analysis_sop.md` - the group's SOP. It defines the two readouts and what is counted as a cell,
  the validation of the counting procedure, image QC, the response scale, the curve model, the
  potency and the meaning of `inactive`. These definitions are binding.

Available software (on `PATH`): Python 3.12 with numpy, scipy, pandas, scikit-image, tifffile,
matplotlib and Cellpose 3 (CPU only; the pretrained `nuclei` and `cyto3` models are installed,
nothing can be downloaded). There is no internet access.

Your task:
1. Quantify every well: the number of cells in the field and the marker signal per cell (SOP
   section 1), and decide whether the well passes image QC (SOP section 3).
2. Validate your counting procedure against the manual counts as SOP section 2 requires, and
   report the comparison for the procedure you finally used.
3. For every test compound, determine the potency on the marker (IC50) and on cell number (GI50)
   as defined in SOP sections 4 and 5, or report that the compound is inactive on that readout.

Required outputs under `/task/output`:
- `well_measurements.tsv` - tab-separated, header `well`, `n_cells`, `marker_per_cell`, `qc_pass`;
  one row for every well of `plate_layout.tsv`. Extra columns are allowed.
  - `n_cells`: integer, SOP readout R1 of the field as counted by your procedure (for the
    annotated wells too).
  - `marker_per_cell`: SOP readout R2 of the well in grey levels (`NA` if the field holds no cell).
  - `qc_pass`: `1` if the well passes the SOP's image QC, `0` if it fails. Wells that fail are
    still listed, with the values your analysis produced for them.
- `validation.tsv` - tab-separated, header `field`, `manual_count`, `pipeline_count`; one row for
  each annotated field and no other rows. `field` is the well id, `manual_count` the number of
  annotated cells, `pipeline_count` the count of your final procedure for that field (the same
  number as `n_cells` of that well in `well_measurements.tsv`).
- `dose_response.tsv` - tab-separated, header `compound`, `ic50_marker_uM`, `gi50_count_uM`; one
  row per test compound (no rows for vehicle or reference). Each potency is either a number in
  micromolar or the literal `inactive`, as defined in SOP section 5.

Constraints:
- Do not modify anything under `/task/input`; write deliverables only under `/task/output`.
- Do not ask for confirmation. Execute directly.
"""

SOP = """\
# SOP IA-07 (v5): image analysis of two-channel potency plates (nuclear stain / cytoplasmic marker)

This SOP applies to every plate of this assay format. It defines what is reported and which
criteria a result has to meet, not how a particular plate has to be processed.

## 1 Readouts (per well)
- **R1, cell number**: the number of cells in the imaged field. One nucleus is one cell; nuclei
  that touch or overlap each other are separate cells. Two rules decide what is counted:
  - *Fragments are not cells.* Nuclear fragments, apoptotic bodies and other stained debris are
    compact objects of a few micrometres (the largest are about 6.5 um across). Intact nuclei of
    this cell line measure at least 11 um in their shortest dimension, with or without treatment.
  - *Counting frame.* A nucleus that is cut by the top or the left edge of the image is counted.
    A nucleus that is cut by the bottom or the right edge of the image is not counted (also when
    it is cut by the top or the left edge as well).
- **R2, marker per cell**: for each counted cell, the mean marker-channel intensity of its
  cytoplasm that is due to the stain bound in the cell. Marker-channel signal that does not come
  from the cells is background and is not part of R2. The well value is the mean over the counted
  cells of the field.

## 2 Validation of the counting procedure
On every plate the lab counts a number of fields by hand under the rules of section 1 (point
annotations). Before a counting procedure is used on a plate, its counts are compared with these
manual counts. The procedure is fit for use when, for every annotated field,

    | procedure count - manual count |  <=  max(3 cells, 10 % of the manual count)

A procedure that does not meet this on every annotated field is not used: it is changed and
compared again. The comparison is reported for the procedure that was finally used. Annotated
fields are counted by the procedure like every other field; manual counts are not entered into
the well table.

## 3 Image QC
A field that was acquired out of focus cannot be quantified. Such a well fails QC (`qc_pass` = 0).
It stays in the well table but is not used for anything downstream (neither as a control nor as a
point of a curve). Every other well passes (`qc_pass` = 1). Few cells, a weak signal or any other
compound effect is a result, not a QC failure.

## 4 Response scale
Responses are relative to vehicle-treated cells:

    response = readout of the well / readout the same well would have given with vehicle only

so that vehicle-treated cells have response 1 for both readouts. The reference inhibitor (single
concentration) is a plate control; it does not define the scale.

## 5 Curve model, potency, `inactive`
For each compound and readout, the responses of the usable wells are described by a four-parameter
logistic function of the concentration c whose upper asymptote is fixed at the vehicle level:

    response(c) = bottom + (1 - bottom) / (1 + (c / EC50)^h),   0 <= bottom < 1,  h > 0

The reported potency is the concentration at which the fitted curve equals 0.5, i.e. the readout
is half of that of vehicle-treated cells (IC50 for R2, GI50 for R1). This equals EC50 only when
bottom = 0.

A compound is **inactive** on a readout when the fitted curve does not come down to 0.5 at or
below the highest concentration tested in a usable well (no effect, a plateau above 0.5, a
potency beyond the tested range, or a readout that does not fall at all). In that case the
literal `inactive` is reported, never a number.
"""

JUDGE = """\
**Where the truth comes from.** The simulator fixes, before rendering, each compound's response
curve on both readouts, each well's cell number, every nucleus (position, size, outline), every
nuclear fragment, each cell's marker level, the plate-position factors and which fields are out of
focus. `reference/truth.json` stores per well the number of rendered nuclei that the SOP counts
(fragments are not nuclei; nuclei cut by the bottom or right image edge are not counted), the focus
state and whether the field was annotated, and per compound and readout the concentration at which
the planted curve equals 0.5 (or null = `inactive`). The manual annotation of a field is the list of
the counted nuclei of that field (centres jittered by about a pixel), so the manual count equals
the planted count. Definitions are unambiguous by construction: fragments are 2.6-5.0 px
(3.4-6.5 um) across and lie at least 3.5 px from any nuclear outline, 7 px from the image edge and
8 px from each other; nuclei are at least 8.8 px (11.5 um) in their shortest dimension; a nucleus
either stays at least 3 px inside an image edge or extends at least 2.5 px beyond it, and a cut
nucleus shows at least 5 px of its extent. No segmentation and no fit enters the reference.

**Hard gates (score 0).** Missing/unparseable file; `well_measurements.tsv` not listing every well
exactly once; `n_cells` not a non-negative integer; `qc_pass` not 0/1; `validation.tsv` not listing
every annotated field exactly once, listing another field, or holding a count that is not a
non-negative integer; `dose_response.tsv` not listing every test compound exactly once; a potency
that is neither a positive number nor `inactive`.

**Components.**

| weight | component | metric | full credit | zero |
|---|---|---|---|---|
| 0.25 | cell_counts | fraction of in-focus wells *without annotation* with abs(n_cells - planted) <= max(3, 10%) | >= 0.85 | <= 0.65 |
| 0.05 | validation | fraction of annotated fields whose row is truthful (manual count = number of annotation points, pipeline count = `n_cells` of that well) and inside the SOP criterion | all | <= half |
| 0.05 | qc_flags | Matthews correlation of `qc_pass = 0` with the out-of-focus wells | >= 0.90 | <= 0.50 |
| 0.45 | marker_potency | mean per-compound credit (inactive must match; active: abs(log10(est/planted)) <= 0.20 -> 1, >= 0.40 -> 0) | mean = 1.0 | mean <= 0.8 |
| 0.20 | count_potency | same for GI50 | mean = 1.0 | mean <= 0.8 |

**Why these tolerances.** The count tolerance per field is the SOP's own criterion, max(3 cells,
10%). Measured in the gate runs of seeds 1 / 2, wells without annotation inside it: reference solution
(outline voting, no Cellpose) 100% / 100%; second solution (Cellpose `nuclei` model with the setting
chosen on the annotated fields) 96.4% / 100%; Cellpose with the level-3 settings 43% / 52% (median
error -11% / -10%); threshold + distance-transform watershed 49% / 53%; LoG-seeded watershed 55% / 52%;
the last two after a plate-wide factor fitted to the annotated fields 23% / 28% and 4% / 1%. The band
0.65-0.85 lies between the two groups. Potency tolerances are those of levels 1-3: the two solutions
are within 0.09 log10 units on every active entry; each ignored decision point gives a wrong activity
class or an error of 0.56-1.33 log10 units in at least one entry. The exception is Cellpose with the
level-3 settings: it loses nuclei almost uniformly in dense fields, the response ratios survive
(largest potency error 0.285 / 0.148) and it fails on the count and validation components alone
(0.647 / 0.700). That is why the count component weighs 0.25: a submission whose counts fail cannot
exceed 0.70. The steep ramp on the mean potency credit (one of eight compounds wrong -> component
0.375) is kept: a potency table with one wrong compound is a wrong deliverable. MCC is used for the QC
flags because 5 of 96 wells fail.

**What the validation component can and cannot do.** It pays 0.05 for a truthful table whose
counts meet the criterion. Copying the manual counts into the annotated wells earns these 0.05 and
nothing else: the annotated wells are excluded from `cell_counts`, and the potencies are fitted
from all wells.

**What is deliberately not scored.** The absolute value of `marker_per_cell`, the values reported
for out-of-focus wells, Hill slopes and plateaus, any plot, the method used for counting.
"""


def _fmt(v):
    return "inactive" if v is None else f"{v:.5g}"


def _wm_tsv(wells, n_fn, q_fn) -> str:
    out = ["well\tn_cells\tmarker_per_cell\tqc_pass"]
    for w in wells:
        m = f"{w.marker_mean:.2f}" if w.n_cells else "NA"
        out.append(f"{w.well}\t{int(n_fn(w))}\t{m}\t{int(q_fn(w))}")
    return "\n".join(out) + "\n"


def _val_tsv(wells, annotated, n_fn) -> str:
    return "field\tmanual_count\tpipeline_count\n" + "".join(
        f"{w.well}\t{w.n_cells}\t{int(n_fn(w))}\n" for w in wells if w.well in annotated)


def _dr_tsv(cs, m_fn, c_fn) -> str:
    return "compound\tic50_marker_uM\tgi50_count_uM\n" + "".join(
        f"{c.name}\t{_fmt(m_fn(c))}\t{_fmt(c_fn(c))}\n" for c in sorted(cs, key=lambda c: c.name))


def _notes(base: str, n_ann: int) -> str:
    extra = [
        "",
        "## Manual counts",
        f"- {n_ann} fields of this plate were counted by hand by two technicians with the point tool of the "
        "image viewer, following SOP IA-07 section 1 (annotations/). Fields were picked to cover untreated "
        "wells of the inner and outer plate area and wells of the highest and of intermediate "
        "concentrations.",
        "- Wells were not washed between treatment and fixation; floating material was fixed with the cells.",
    ]
    return base + "\n".join(extra) + "\n"


def generate_l4(seed: int, level: int, workdir: Path, notes_fn, plate_out: dict | None = None) -> TaskBuild:
    """`plate_out` is None at levels 4-6. Level 7 (l7.py) passes a dict that receives the planted plate
    (wells, compounds, nucleus positions of every field), so that it can file the images under other names."""
    pl = sim4.build_plate(seed, workdir)
    if plate_out is not None:
        plate_out.update(pl)
    cs, wells, plate, shapes = pl["cs"], pl["wells"], pl["plate"], pl["shapes"]
    blurred, annotated = pl["blurred"], set(pl["annotated"])

    input_files = dict(pl["files"])
    input_files["plate_layout.tsv"] = "well\tcompound\tdose_uM\trole\n" + "".join(
        f"{w.well}\t{w.compound}\t{w.dose:.4g}\t{w.role}\n" for w in wells)
    input_files["analysis_sop.md"] = SOP
    input_files["assay_notes.md"] = _notes(notes_fn(3, plate, np.random.default_rng(seed + 431)), len(annotated))

    truth = {
        "wells": {w.well: {"n_cells": w.n_cells, "qc_fail": int(w.defocus > 0), "annotated": int(w.well in annotated),
                           "role": w.role, "compound": w.compound, "dose_uM": w.dose,
                           "planted_marker_mean": round(w.marker_mean, 3),
                           "stain_factor": round(w.stain_factor, 4), "count_factor": round(w.count_factor, 4),
                           "background": round(w.background, 2), "defocus_sigma_px": round(w.defocus, 2),
                           "nuclear_size_factor": round(w.extra["size"], 3), "n_fragments": w.extra["n_debris"],
                           "n_nuclei_outside_counting_frame": w.extra["n_cut_uncounted"]}
                  for w in wells},
        "compounds": {c.name: {"marker_ic50": c.marker.ic50, "count_gi50": c.count.ic50, "class": c.cls,
                               "row": sim.ROWS[c.row],
                               "marker_curve": {"kind": c.marker.kind, "ec50": c.marker.ec50, "h": c.marker.h,
                                                "bottom": c.marker.bottom},
                               "count_curve": {"kind": c.count.kind, "ec50": c.count.ec50, "h": c.count.h,
                                               "bottom": c.count.bottom},
                               "autofluorescence_top_dose_x_M0": c.autofluor,
                               "marker_rise_per_fraction_lost": shapes["marker_rise"].get(c.name, 0.0),
                               "nuclear_size_factor_max": shapes["size_max"].get(c.name, 1.0)} for c in cs},
        "plate": {**plate, "dose_ascending_with_column": pl["ascending"], "out_of_focus": blurred,
                  "annotated_fields": sorted(annotated), "image_model": sim4.PARAMS},
        "tolerance": TOL, "control_names": CONTROL_NAMES,
    }

    # ---- battery: truth-derived probes / attacks -------------------------------------------------
    prng = np.random.default_rng(seed + 9409)
    bad_focus = {w.well for w in wells if w.defocus > 0}
    count_jit = {w.well: (1.0 + 0.07 * float(prng.choice([-1, 1]))) for w in wells}
    pot_jit = {c.name: (10 ** (0.15 * float(prng.choice([-1, 1]))), 10 ** (0.15 * float(prng.choice([-1, 1]))))
               for c in cs}
    cand = [w.well for w in wells if w.well not in bad_focus]
    extra_flag = cand[int(prng.integers(0, len(cand)))]

    def scaled(v, f):
        return None if v is None else v * f

    def n_true(w):
        return w.n_cells

    def q_true(w):
        return w.well not in bad_focus

    def n_jit(w):
        return round(w.n_cells * count_jit[w.well])

    def n_low(w):
        return round(w.n_cells * 0.8)

    true_wm = _wm_tsv(wells, n_true, q_true)
    true_val = _val_tsv(wells, annotated, n_true)
    true_dr = _dr_tsv(cs, lambda c: c.marker.ic50, lambda c: c.count.ic50)
    first_marker_active = sorted(c.name for c in cs if c.marker.ic50 is not None)[0]

    l3 = {"hcs_pipeline.py": ASSETS / "hcs_pipeline.py", "l3_solve_alt.py": ASSETS / "solve_alt.py",
          "l3_validation.py": ASSETS / "l3_validation.py"}
    l3_cp = "python3 /alegen_script/hcs_pipeline.py --segmenter cellpose"
    l3_ws = "python3 /alegen_script/hcs_pipeline.py --segmenter watershed"
    l3_log = "python3 /alegen_script/l3_solve_alt.py"
    val = " && python3 /alegen_script/l3_validation.py"
    ref = {"hcs_l4_voting.py": ASSETS / "hcs_l4_voting.py"}
    base = "python3 /alegen_script/hcs_l4_voting.py "

    probes = [
        Probe("planted_truth", "submission written directly from the planted cell numbers, focus states and "
              "potencies; validation table with the planted counts of the annotated fields", 0.999,
              files={"well_measurements.tsv": true_wm, "validation.tsv": true_val, "dose_response.tsv": true_dr}),
        Probe("within_tolerance", "every count off by 7%, every potency off by 0.15 log10 units, one in-focus "
              "well flagged", 0.99,
              files={"well_measurements.tsv": _wm_tsv(wells, n_jit, lambda w: q_true(w) and w.well != extra_flag),
                     "validation.tsv": _val_tsv(wells, annotated, n_jit),
                     "dose_response.tsv": _dr_tsv(cs, lambda c: scaled(c.marker.ic50, pot_jit[c.name][0]),
                                                  lambda c: scaled(c.count.ic50, pot_jit[c.name][1]))}),
        Probe("independent_solution", "second implementation, no code shared with the reference solution: "
              "Cellpose nuclei model with settings chosen on the annotated fields (the level-3 settings are tried "
              "first and rejected), fragments by mask area, median background far from nuclei, edge-gradient "
              "focus metric, least-squares plate model, scipy least_squares fit", 0.99,
              script="SOLUTION=2 bash /solution/solve.sh"),
    ]
    attacks = [
        # --- mandatory: the level-3 methods on the level-4 inputs (each with an honest validation table)
        Attack("l3_reference", "the level-3 reference solution unchanged (Cellpose nuclei model with the "
               "level-3 settings, every mask >= 20 px counted, level-3 plate model / QC / fit) plus a truthful "
               "validation table", 0.75, script=l3_cp + val, script_files=l3),
        Attack("l3_reference_classical", "the level-3 classical route unchanged (threshold + distance-transform "
               "watershed; the route the pilot agent took) plus a truthful validation table", 0.75,
               script=l3_ws + val, script_files=l3),
        Attack("l3_independent_solution", "the level-3 second solution unchanged (LoG-seeded watershed) plus a "
               "truthful validation table", 0.75, script=l3_log + val, script_files=l3),
        # --- Q_validated_counts
        Attack("counts_scaled_to_annotation", "level-3 classical route; the counts are then multiplied by one "
               "plate-wide factor that makes the annotated fields agree in total", 0.75,
               script=l3_ws + val + " --scale", script_files=l3),
        Attack("counts_scaled_to_annotation_log", "level-3 LoG route with the same plate-wide factor", 0.75,
               script=l3_log + val + " --scale", script_files=l3),
        Attack("manual_counts_copied", "level-3 classical route; the counts of the annotated wells are replaced "
               "by the manual counts so that the validation table looks perfect", 0.75,
               script=l3_ws + val + " --copy-manual", script_files=l3),
        Attack("naive_fragments_counted", "reference solution, everything else correct, but nuclear fragments "
               "are counted as cells", 0.75, script=base + "--counting fragments", script_files=ref),
        # --- P_potency_shape
        Attack("naive_free_plateau_fit", "reference solution, everything else correct, but the curves are fitted "
               "with free upper and lower plateaus and the fitted EC50 is reported", 0.75,
               script=base + "--fit free", script_files=ref),
        # --- level-3 decision points, each ignored on the correct level-4 base
        Attack("naive_no_focus_qc", "reference solution, but no focus QC: every well passes and is used", 0.75,
               script=base + "--qc none", script_files=ref),
        Attack("naive_whole_well_readout", "reference solution, but the marker readout is the mean "
               "background-subtracted intensity of the field instead of the intensity per cell", 0.75,
               script=base + "--readout wholewell", script_files=ref),
        Attack("naive_no_plate_model", "reference solution, but responses are normalised to the plate mean of the "
               "vehicle wells", 0.75, script=base + "--plate none", script_files=ref),
        Attack("naive_no_background", "reference solution, but the marker per cell is the raw ring intensity above "
               "the camera offset", 0.75, script=base + "--background none", script_files=ref),
        Attack("naive_global_background", "reference solution, but one plate-wide background constant is "
               "subtracted from every well", 0.75, script=base + "--background global", script_files=ref),
        # --- truth-derived
        Attack("prior_only", "no analysis: every well 100 cells and QC pass, every compound inactive, validation "
               "table claims 100 cells per field", 0.1,
               files={"well_measurements.tsv": _wm_tsv(wells, lambda w: 100, lambda w: 1),
                      "validation.tsv": _val_tsv(wells, annotated, lambda w: 100),
                      "dose_response.tsv": _dr_tsv(cs, lambda c: None, lambda c: None)}),
        Attack("all_potencies_beyond_tolerance", "planted counts and QC flags, every active potency 0.45 log10 "
               "units too high", 0.40,
               files={"well_measurements.tsv": true_wm, "validation.tsv": true_val,
                      "dose_response.tsv": _dr_tsv(cs, lambda c: scaled(c.marker.ic50, 10 ** 0.45),
                                                   lambda c: scaled(c.count.ic50, 10 ** 0.45))}),
        Attack("one_potency_beyond_tolerance", "planted truth, but one marker IC50 0.45 log10 units too high", 0.75,
               files={"well_measurements.tsv": true_wm, "validation.tsv": true_val,
                      "dose_response.tsv": _dr_tsv(
                          cs, lambda c: scaled(c.marker.ic50, 10 ** 0.45 if c.name == first_marker_active else 1.0),
                          lambda c: c.count.ic50)}),
        Attack("counts_beyond_tolerance", "planted potencies and QC flags, every cell count 20% too low, "
               "validation table truthful about it", 0.75,
               files={"well_measurements.tsv": _wm_tsv(wells, n_low, q_true),
                      "validation.tsv": _val_tsv(wells, annotated, n_low), "dose_response.tsv": true_dr}),
    ]

    by_cls = {c.cls: c for c in cs}
    traps = [
        Trap("Q_validated_counts",
             "Counting quality decides the task. Nuclei overlap without an intensity dip (centre distance down "
             f"to {sim4.PARAMS['dmin']} px for 12-15 x 10-12 px nuclei), chromatin is textured, every field holds "
             "bright nuclear fragments (about 12% of the seeded cells, many more where cells were lost) that the "
             "SOP does not count, nuclei are cut by the field border under a counting-frame rule, and two "
             "compounds enlarge the nuclei up to 1.3-1.5x. The level-3 methods fail the SOP criterion on the "
             "annotated fields and on the other wells; the errors have different signs in dense, sparse and "
             "treated fields, so no plate-wide factor repairs them.",
             "SOP section 2 states the acceptance criterion; annotations/ holds the manual point annotation of "
             "eight fields (dense and sparse untreated fields, fragment-rich and enlarged-nuclei fields); "
             "comparing any counting procedure with them, or overlaying its objects on the annotated points, "
             "shows which nuclei are lost and which objects are counted wrongly",
             "Verification / result not checked against the available control"),
        Trap("P_potency_shape",
             f"The count response of {by_cls['marker_potent'].name} is a deep partial response (plateau "
             f"{by_cls['marker_potent'].count.bottom:.2f}): the SOP potency (0.5 crossing) is "
             f"{math.log10(by_cls['marker_potent'].count.ic50 / by_cls['marker_potent'].count.ec50):.2f} log10 "
             "units above the inflection point that a free four-parameter fit reports. The marker per cell of "
             f"{by_cls['cytotoxic'].name} rises with the fraction of cells lost (factor "
             f"{1 + shapes['marker_rise'][by_cls['cytotoxic'].name]:.2f} when all are lost); a fit with free "
             "plateaus returns an EC50 for it although nothing is inhibited.",
             "SOP section 5 fixes the upper asymptote at the vehicle level and defines the potency as the 0.5 "
             "crossing and `inactive`; plotting the responses of each series shows the plateau at about a third "
             "and the rising series",
             "Understanding / textbook default instead of the written definition"),
        Trap("touching_nuclei",
             "Level 2 decision point, sharpened by Q_validated_counts: a global threshold merges the nuclei of "
             "a colony.",
             "overlay of the threshold mask on any vehicle field; comparison with the manual counts",
             "Approach / wrong strategy"),
        Trap("focus_qc",
             f"{len(blurred)} fields ({', '.join(blurred)}) are out of focus (extra blur sigma 4.5-6 px): two "
             "neighbouring series wells just below the marker potency of one compound, two just below the count "
             "potency of another, one vehicle well.",
             "SOP section 3 defines the QC failure; a signal-normalised sharpness metric of channel 1 puts these "
             "wells far below all others",
             "Verification / skipped QC step"),
        Trap("readout_confound",
             "Two compounds reduce cell number without reducing the marker per cell; the field intensity "
             "follows cell number.",
             "SOP section 1 defines R2 per cell; field intensity tracks cell number along those series",
             "Understanding / wrong readout"),
        Trap("plate_position_effects",
             f"Marker staining falls log-linearly to {plate['stain_end']:.2f}x along the columns (dose runs along "
             f"the columns) and the outer wells hold {plate['edge_count']:.2f}x the cells of inner wells.",
             "plate map of the 24 vehicle wells (3 per row, 2 per column) for each readout; assay notes give the "
             "pipetting order and the incubation conditions as facts",
             "Verification / skipped QC step"),
        Trap("marker_background",
             "Diffuse marker-channel background of about a third of the vehicle signal that differs between "
             "wells; one marker inhibitor fluoresces in the marker channel at high concentration.",
             "SOP section 1 excludes signal that does not come from the cells from R2; cell-free regions of "
             "channel 2 are far above the dark frame, differ between wells and rise with concentration in one "
             "series",
             "Approach / skipped QC step"),
    ]

    assert sum(c.marker.ic50 is not None for c in cs) == 4 and sum(c.count.ic50 is not None for c in cs) == 3
    assert len(wells) == 96 and len(annotated) == 8 and not (annotated & bad_focus)
    for c in cs:   # activity classes are separated by construction
        for curve in (c.marker, c.count):
            top = float(curve(sim.TOP_DOSE))
            assert (curve.ic50 is not None and sim.DOSES[0] * 3.2 < curve.ic50 <= sim.TOP_DOSE / 5.4 and top < 0.4) \
                or (curve.ic50 is None and top >= 0.59), (c.cls, curve)
    for w in wells:
        if w.well in annotated:
            assert len(pl["points"][w.well]) == w.n_cells

    return TaskBuild(
        task_id=f"hcs-dose-response-s{seed:04d}-l{level}", family="hcs_dose_response", seed=seed, level=level,
        title="High-content screening plate: validated cell counts and compound potencies",
        summary="Count the cells of a two-channel 96-well imaging plate with a procedure validated against the "
                "lab's manual counts, quantify the marker per cell, apply image QC and derive per-compound "
                "IC50/GI50 values under the SOP's definitions.",
        subdomain="3.3 Cell & Imaging Biology",
        software=["Cellpose", "scikit-image", "scipy", "tifffile", "Python"],
        instruction_md=INSTRUCTION.format(n_wells=len(wells), size=sim.SIZE, n_ann=len(annotated)),
        agent_dockerfile=image.agent_dockerfile(),
        verifier_dockerfile=image.verifier_dockerfile(),
        input_files=input_files,
        solution_files={"solve.sh": ASSETS / "solve_l4.sh", "hcs_l4_voting.py": ASSETS / "hcs_l4_voting.py",
                        "hcs_l4_cellpose.py": ASSETS / "hcs_l4_cellpose.py"},
        tests_files={"score.py": ASSETS / "score_l4.py", "reference/truth.json": json.dumps(truth, indent=1)},
        traps=traps, attacks=attacks, probes=probes,
        capabilities=["validating a segmentation against manual annotation before using it",
                      "tuning a pretrained segmentation model or building a classical alternative to a stated "
                      "quality bar", "applying written counting rules (fragments, counting frame)",
                      "plate-level QC before curve fitting", "per-well image QC",
                      "dose-response modelling under a written potency definition", "exact output contracts"],
        judge_explanation_md=JUDGE, pass_threshold=0.9,
        oracle_steps=1, oracle_runtime_min=0.2, cpus=4, memory_mb=8000, storage_mb=10000,
    )
