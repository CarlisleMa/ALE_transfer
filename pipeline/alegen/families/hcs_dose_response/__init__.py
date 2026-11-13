"""Family: high-content screening plate, 2-channel well images -> per-compound potencies.

Truth-by-construction: a forward simulator decides, before any pixel exists, the potency of every
compound on both readouts (marker per cell, cell number), the number of cells of every well, the
marker level of every cell, the plate-position factors, the background of every well and which
wells are acquired out of focus; the images are rendered from these values (sim.py). The reference
is the planted potency / planted cell count / planted focus state. No segmentation and no curve
fit is used as truth.

The agent image ships Cellpose 3 (CPU) with pretrained weights. Running the model is necessary for
nothing and sufficient for nothing: a classical watershed pipeline reaches full score as well, and
the Cellpose default pipeline without plate-level reasoning scores near zero.

Level 4 (level4.py, sim4.py) adds Q_validated_counts and P_potency_shape; there the counting method
does matter and has to be validated against manual annotation.
Level 7 (l7.py) is the level-4 plate of the seed with one silent decision point, T_plate_orientation:
the plate sat on the stage turned by 180 degrees, so every image is filed under the opposite well name.

Decision points (level 1: none; level 2: the first three; level 3: all five):
  touching_nuclei          cells grow in colonies, part of the nuclei abut their neighbours
  focus_qc                 five fields are out of focus, placed just below a potency
  readout_confound         cytotoxic compounds lower the field signal but not the marker per cell
  plate_position_effects   staining gradient along the columns (marker) and outer-well effect (cells)
  marker_background        diffuse marker-channel background that differs from well to well, and one
                           compound that fluoresces in the marker channel at high concentration
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from ...spec import Attack, Probe, TaskBuild, Trap
from . import image, sim

ASSETS = Path(__file__).parent / "assets"
# scoring tolerances (kept here so that probes, attacks and the verifier use the same numbers)
TOL = {"count_abs": 3, "count_rel": 0.10, "log10_full": 0.20, "log10_zero": 0.40}
CONTROL_NAMES = ["DMSO", "REF-INH"]

INSTRUCTION = """\
You are the image-analysis scientist of a compound-profiling group. A 96-well plate of cells was
treated with eight test compounds in dose series, fixed, stained for nuclei and for a cytoplasmic
marker protein, and imaged on an automated microscope. The project team needs per-well
measurements and, for every test compound, its potency on the marker and on cell number.

Task directory:
- `/task`

Inputs (read-only) under `/task/input`:
- `images/<well>_ch1.tif`, `images/<well>_ch2.tif` - one field per well ({n_wells} wells), 16-bit
  grayscale TIFF, {size} x {size} pixels. Channel 1 = nuclear stain, channel 2 = marker.
- `flatfield/flat_ch1.tif`, `flatfield/flat_ch2.tif`, `flatfield/dark.tif` - reference frames
  recorded with the same optical path (uniform fluorescent slide per channel; closed shutter).
- `plate_layout.tsv` - columns `well`, `compound`, `dose_uM`, `role` (`compound`, `vehicle` or
  `reference`).
- `assay_notes.md` - hand-over notes of the wet lab and of the imaging operator.
- `analysis_sop.md` - the group's SOP. It defines the two readouts, image QC, the response scale,
  the curve model, the potency and the meaning of `inactive`. These definitions are binding.

Available software (on `PATH`): Python 3.12 with numpy, scipy, pandas, scikit-image, tifffile,
matplotlib and Cellpose 3 (CPU only; the pretrained `nuclei` and `cyto3` models are installed,
nothing can be downloaded). There is no internet access.

Your task:
1. Quantify every well: the number of cells in the field and the marker signal per cell, and
   decide whether the well passes image QC (SOP sections 1 and 2).
2. For every test compound, determine the potency on the marker (IC50) and on cell number (GI50)
   as defined in SOP sections 3 and 4, or report that the compound is inactive on that readout.

Required outputs under `/task/output`:
- `well_measurements.tsv` - tab-separated, header `well`, `n_cells`, `marker_per_cell`, `qc_pass`;
  one row for every well of `plate_layout.tsv`. Extra columns are allowed.
  - `n_cells`: integer, SOP readout R1 of the field.
  - `marker_per_cell`: SOP readout R2 of the well in grey levels (`NA` if the field holds no cell).
  - `qc_pass`: `1` if the well passes the SOP's image QC, `0` if it fails. Wells that fail are
    still listed, with the values your analysis produced for them.
- `dose_response.tsv` - tab-separated, header `compound`, `ic50_marker_uM`, `gi50_count_uM`; one
  row per test compound (no rows for vehicle or reference). Each potency is either a number in
  micromolar or the literal `inactive`, as defined in SOP section 4.

Constraints:
- Do not modify anything under `/task/input`; write deliverables only under `/task/output`.
- Do not ask for confirmation. Execute directly.
"""

SOP = """\
# SOP IA-07 (v4): image analysis of two-channel potency plates (nuclear stain / cytoplasmic marker)

This SOP applies to every plate of this assay format. It defines what is reported, not how a
particular plate has to be processed.

## 1 Readouts (per well)
- **R1, cell number**: the number of cells in the imaged field. One nucleus is one cell.
- **R2, marker per cell**: for each cell, the mean marker-channel intensity of its cytoplasm that
  is due to the stain bound in the cell. Marker-channel signal that does not come from the cells
  is background and is not part of R2. The well value is the mean over the cells of the field.

## 2 Image QC
A field that was acquired out of focus cannot be quantified. Such a well fails QC (`qc_pass` = 0).
It stays in the well table but is not used for anything downstream (neither as a control nor as a
point of a curve). Every other well passes (`qc_pass` = 1). Few cells, a weak signal or any other
compound effect is a result, not a QC failure.

## 3 Response scale
Responses are relative to vehicle-treated cells:

    response = readout of the well / readout the same well would have given with vehicle only

so that vehicle-treated cells have response 1 for both readouts. The reference inhibitor (single
concentration) is a plate control; it does not define the scale.

## 4 Curve model, potency, `inactive`
For each compound and readout, the responses of the usable wells are described by a four-parameter
logistic function of the concentration c whose upper asymptote is fixed at the vehicle level:

    response(c) = bottom + (1 - bottom) / (1 + (c / EC50)^h),   0 <= bottom < 1,  h > 0

The reported potency is the concentration at which the fitted curve equals 0.5, i.e. the readout
is half of that of vehicle-treated cells (IC50 for R2, GI50 for R1). This equals EC50 only when
bottom = 0.

A compound is **inactive** on a readout when the fitted curve does not come down to 0.5 at or
below the highest concentration tested in a usable well (no effect, a plateau above 0.5, or a
potency beyond the tested range). In that case the literal `inactive` is reported, never a number.
"""

JUDGE = """\
**Where the truth comes from.** The simulator fixes, before rendering, each compound's response
curve on both readouts, each well's cell number, each cell's marker level, the plate-position
factors and which fields are out of focus. `reference/truth.json` stores the planted cell number
and focus state per well and, per compound and readout, the concentration at which the planted
curve equals 0.5 (or null = `inactive`: no effect, plateau >= 0.60, or potency >= 5x the top
concentration). Active potencies lie between 3.3x the lowest and 1/5.5 of the highest tested
concentration, so "active" and "inactive" are separated by construction. No segmentation and no
fit enters the reference.

**Hard gates (score 0).** Missing/unparseable file; `well_measurements.tsv` not listing every well
exactly once; `n_cells` not a non-negative integer; `qc_pass` not 0/1; `dose_response.tsv` not
listing every test compound exactly once; a potency that is neither a positive number nor `inactive`.

**Components.**

| weight | component | metric | full credit | zero |
|---|---|---|---|---|
| 0.20 | cell_counts | fraction of in-focus wells with abs(n_cells - planted) <= max(3, 10%) | >= 0.90 | <= 0.50 |
| 0.10 | qc_flags | Matthews correlation of `qc_pass = 0` with the out-of-focus wells | >= 0.90 | <= 0.50 |
| 0.45 | marker_potency | mean per-compound credit (inactive must match; active: abs(log10(est/planted)) <= 0.20 -> 1, >= 0.40 -> 0) | mean = 1.0 | mean <= 0.8 |
| 0.25 | count_potency | same for GI50 | mean = 1.0 | mean <= 0.8 |

**Why these tolerances.** Measured with two classical pipelines (distance-transform watershed,
LoG-seeded watershed) on six seeds and with Cellpose on two: the median per-well count error is
between -1.7% and 0%, the worst in-focus well is 11% low (16 of 18 cells; worst dense well 108 of 118), every in-focus
well is inside max(3, 10%); the QC flags are exact; every active potency is recovered within 0.11 log10
units (typically 0.01-0.07).
Ignoring a decision point moves the affected potencies by 0.25-1.3 log10 units or flips the activity
class. The steep ramp on the mean credit (one of eight compounds wrong -> component 0.375) is
deliberate: a potency table with one wrong compound is a wrong deliverable, and with a flat mean a
single-compound error (the fluorescent compound under a plate-wide background constant) would cost
only 0.06. MCC is used for the QC flags because 5 of 96 wells fail: "all pass" has accuracy 0.95 but
MCC 0. On a plate without out-of-focus wells (level 1) the QC component is (1 - false flags / 5)
times the cell-count credit, so that "flag nothing" earns nothing by itself.

**What is deliberately not scored.** The absolute value of `marker_per_cell` (it depends on the
legitimate choice of the cytoplasm region; only the potencies derived from it are scored), the
values reported for out-of-focus wells, Hill slopes and plateaus, any plot.
"""


def _fmt(v: float | None) -> str:
    return "inactive" if v is None else f"{v:.5g}"


def _wm_tsv(wells, n_fn, q_fn) -> str:
    out = ["well\tn_cells\tmarker_per_cell\tqc_pass"]
    for w in wells:
        m = f"{w.marker_mean:.2f}" if w.n_cells else "NA"
        out.append(f"{w.well}\t{int(n_fn(w))}\t{m}\t{int(q_fn(w))}")
    return "\n".join(out) + "\n"


def _dr_tsv(cs, m_fn, c_fn) -> str:
    return "compound\tic50_marker_uM\tgi50_count_uM\n" + "".join(
        f"{c.name}\t{_fmt(m_fn(c))}\t{_fmt(c_fn(c))}\n" for c in sorted(cs, key=lambda c: c.name))


def _notes(level: int, plate: dict, rng: np.random.Generator) -> str:
    first, last = (1, 12) if plate["stain_dir"] >= 0 else (12, 1)
    day = int(rng.integers(3, 27))
    lines = [
        f"# Assay notes - plate HCS-{int(rng.integers(2000, 9000))} (imaged 2025-{int(rng.integers(2, 11)):02d}-{day:02d})",
        "",
        "## Wet lab",
        "- Cells: adherent reporter line (identity blinded), 96-well imaging plate with optical bottom, "
        "all 96 wells in use.",
        "- Compounds: eight blinded test compounds, 8-point 3-fold series from 30 uM, one well per "
        "concentration; vehicle = DMSO (0.3% in every well); REF-INH = reference inhibitor of the marker "
        "pathway at 10 uM. Dispensed acoustically according to plate_layout.tsv, 48 h treatment.",
        ("- Incubation: 37 C / 5% CO2, standard polystyrene lid, plate at the front of the top shelf; the "
         "humidity tray of that incubator was found dry when the plate was taken out."
         if level >= 3 else
         "- Incubation: 37 C / 5% CO2 in a humidity cassette with a gas-permeable seal."),
        "- Fixation 4% PFA, permeabilisation 0.1% Triton X-100.",
        ("- Staining: Alexa Fluor 488 conjugate against the marker protein, added by hand with an 8-channel "
         f"pipette column by column from column {first} to column {last}; the stain was removed from the whole "
         "plate at once with the plate washer 45 min after the first column had been filled. Two PBS "
         "washes, wells left in 100 uL PBS."
         if level >= 3 else
         "- Staining: Alexa Fluor 488 conjugate against the marker protein, added to all wells at once with "
         "the 96-channel head, 45 min; three PBS washes on the plate washer, wells left in 100 uL PBS."),
        "- Nuclear counterstain: Hoechst 33342 in the last wash.",
        "",
        "## Imaging",
        "- Automated widefield microscope, 10x / NA 0.45 air objective, sCMOS camera with 2 x 2 binning "
        "(1.3 um per pixel), 16 bit.",
        "- Channel 1: Hoechst (ex 385 nm). Channel 2: Alexa Fluor 488 (ex 475 nm). One field per well at the "
        "well centre, cropped to 320 x 320 pixels for transfer. Same exposure time for all wells of a channel.",
        "- Focus: hardware autofocus on the plate bottom with one fixed offset for the whole plate; images "
        "were not reviewed at the microscope.",
        "- Reference frames (flatfield/): average of 50 exposures of a uniform dye slide per channel and of "
        "50 exposures with the shutter closed, recorded the same day.",
        "- Files are raw camera frames: no illumination, background or any other correction was applied.",
    ]
    return "\n".join(lines) + "\n"


def generate(seed: int, level: int, workdir: Path) -> TaskBuild:
    if level == 7:      # level 7 lives in l7.py: the level-4 plate, turned by 180 degrees on the stage
        from .l7 import generate_l7
        return generate_l7(seed, level, workdir, _notes)
    if level >= 4:      # level 4 has its own image model, SOP, verifier and battery (level4.py, sim4.py)
        from . import level4
        return level4.generate_l4(seed, level, workdir, _notes)

    import tifffile

    rng = np.random.default_rng(88000 + seed)
    cs = sim.draw_compounds(rng, level)
    sim.assign_rows(rng, cs, level)
    wells, ascending = sim.layout(rng, cs)
    plate = sim.plant_wells(rng, wells, cs, level)
    blurred = sim.plant_defocus(rng, wells, cs, level)

    img_dir = workdir / "images"
    ff_dir = workdir / "flatfield"
    for d in (img_dir, ff_dir):
        if d.exists():
            for f in d.glob("*.tif"):
                f.unlink()
        d.mkdir(parents=True, exist_ok=True)
    il1, il2 = sim.illumination(0.30), sim.illumination(0.26)
    input_files: dict = {}
    for w in wells:
        c1, c2 = sim.render_well(rng, w, level >= 2, il1, il2)
        for ch, im in ((1, c1), (2, c2)):
            p = img_dir / f"{w.well}_ch{ch}.tif"
            tifffile.imwrite(p, im, compression="zlib", compressionargs={"level": 9}, predictor=True)
            input_files[f"images/{p.name}"] = p
    f1, f2, dk = sim.reference_images(rng, il1, il2)
    for name, im in (("flat_ch1", f1), ("flat_ch2", f2), ("dark", dk)):
        p = ff_dir / f"{name}.tif"
        tifffile.imwrite(p, im, compression="zlib", predictor=True)
        input_files[f"flatfield/{p.name}"] = p

    layout_tsv = "well\tcompound\tdose_uM\trole\n" + "".join(
        f"{w.well}\t{w.compound}\t{w.dose:.4g}\t{w.role}\n" for w in wells)
    input_files["plate_layout.tsv"] = layout_tsv
    input_files["analysis_sop.md"] = SOP
    input_files["assay_notes.md"] = _notes(level, plate, np.random.default_rng(seed + 31))

    truth = {
        "wells": {w.well: {"n_cells": w.n_cells, "qc_fail": int(w.defocus > 0), "role": w.role,
                           "compound": w.compound, "dose_uM": w.dose,
                           "planted_marker_mean": round(w.marker_mean, 3),
                           "stain_factor": round(w.stain_factor, 4), "count_factor": round(w.count_factor, 4),
                           "background": round(w.background, 2), "defocus_sigma_px": round(w.defocus, 2)}
                  for w in wells},
        "compounds": {c.name: {"marker_ic50": c.marker.ic50, "count_gi50": c.count.ic50, "class": c.cls,
                               "row": sim.ROWS[c.row],
                               "marker_curve": {"kind": c.marker.kind, "ec50": c.marker.ec50, "h": c.marker.h,
                                                "bottom": c.marker.bottom},
                               "count_curve": {"kind": c.count.kind, "ec50": c.count.ec50, "h": c.count.h,
                                               "bottom": c.count.bottom},
                               "autofluorescence_top_dose_x_M0": c.autofluor} for c in cs},
        "plate": {**plate, "dose_ascending_with_column": ascending, "out_of_focus": blurred},
        "tolerance": TOL, "control_names": CONTROL_NAMES,
    }

    # ---- battery: truth-derived probes / attacks -------------------------------------------------
    prng = np.random.default_rng(seed + 9109)
    bad_focus = {w.well for w in wells if w.defocus > 0}
    count_jit = {w.well: (1.0 + 0.07 * float(prng.choice([-1, 1]))) for w in wells}
    pot_jit = {c.name: (10 ** (0.15 * float(prng.choice([-1, 1]))), 10 ** (0.15 * float(prng.choice([-1, 1]))))
               for c in cs}
    extra_flag = None
    if level >= 2:
        cand = [w.well for w in wells if w.well not in bad_focus]
        extra_flag = cand[int(prng.integers(0, len(cand)))]

    def scaled(v, f):
        return None if v is None else v * f

    true_wm = _wm_tsv(wells, lambda w: w.n_cells, lambda w: w.well not in bad_focus)
    true_dr = _dr_tsv(cs, lambda c: c.marker.ic50, lambda c: c.count.ic50)
    first_marker_active = sorted(c.name for c in cs if c.marker.ic50 is not None)[0]
    pipe = {"hcs_pipeline.py": ASSETS / "hcs_pipeline.py"}
    run = "python3 /alegen_script/hcs_pipeline.py "
    ws = run + "--segmenter watershed "

    probes = [
        Probe("planted_truth", "submission written directly from the planted cell numbers, focus states and "
              "potencies", 0.999, files={"well_measurements.tsv": true_wm, "dose_response.tsv": true_dr}),
        Probe("within_tolerance", "every count off by 7%, every potency off by 0.15 log10 units"
              + (", one in-focus well flagged" if extra_flag else ""), 0.99,
              files={"well_measurements.tsv": _wm_tsv(wells, lambda w: round(w.n_cells * count_jit[w.well]),
                                                      lambda w: w.well not in bad_focus and w.well != extra_flag),
                     "dose_response.tsv": _dr_tsv(cs, lambda c: scaled(c.marker.ic50, pot_jit[c.name][0]),
                                                  lambda c: scaled(c.count.ic50, pot_jit[c.name][1]))}),
        Probe("independent_solution", "second implementation without Cellpose: LoG-seeded watershed, histogram-"
              "mode background, power-spectrum focus metric, Theil-Sen plate model, grid-search fit", 0.99,
              script="SOLUTION=2 bash /solution/solve.sh"),
        Probe("classical_watershed_pipeline", "reference pipeline with distance-transform watershed instead of "
              "Cellpose (the correct base on which the naive variants below are built)", 0.99,
              script=ws.strip(), script_files=pipe),
    ]
    attacks = [
        Attack("prior_only", "no analysis: every well 100 cells and QC pass, every compound inactive", 0.1,
               files={"well_measurements.tsv": _wm_tsv(wells, lambda w: 100, lambda w: 1),
                      "dose_response.tsv": _dr_tsv(cs, lambda c: None, lambda c: None)}),
        Attack("all_potencies_beyond_tolerance", "planted counts and QC flags, every active potency 0.45 log10 "
               "units too high", 0.35,
               files={"well_measurements.tsv": true_wm,
                      "dose_response.tsv": _dr_tsv(cs, lambda c: scaled(c.marker.ic50, 10 ** 0.45),
                                                   lambda c: scaled(c.count.ic50, 10 ** 0.45))}),
        Attack("one_potency_beyond_tolerance", "planted truth, but one marker IC50 0.45 log10 units too high", 0.75,
               files={"well_measurements.tsv": true_wm,
                      "dose_response.tsv": _dr_tsv(
                          cs, lambda c: scaled(c.marker.ic50, 10 ** 0.45 if c.name == first_marker_active else 1.0),
                          lambda c: c.count.ic50)}),
        Attack("counts_beyond_tolerance", "planted potencies and QC flags, every cell count 20% too low", 0.82,
               files={"well_measurements.tsv": _wm_tsv(wells, lambda w: round(w.n_cells * 0.8),
                                                       lambda w: w.well not in bad_focus),
                      "dose_response.tsv": true_dr}),
    ]
    traps: list[Trap] = []
    if level >= 2:
        attacks += [
            Attack("naive_threshold_counting", "everything else correct, but nuclei counted as connected "
                   "components of a global Otsu threshold", 0.75,
                   script=run + "--segmenter threshold", script_files=pipe),
            Attack("naive_no_focus_qc", "everything else correct, but no focus QC: every well passes and is "
                   "used", 0.75, script=ws + "--qc none", script_files=pipe),
            Attack("naive_whole_well_readout", "everything else correct, but the marker readout is the mean "
                   "background-subtracted intensity of the field instead of the intensity per cell", 0.75,
                   script=ws + "--readout wholewell", script_files=pipe),
        ]
        traps += [
            Trap("touching_nuclei",
                 "Cells grow in colonies; a global threshold merges abutting nuclei, undercounts dense wells by "
                 "about half but sparse wells hardly at all, which compresses the count response: the GI50s move "
                 "up by 0.3-0.9 log10 units or the compound is reported inactive.",
                 "any overlay of the threshold mask on a vehicle field; object-area histogram with a long tail of "
                 "multiples of the single-nucleus area; counts that do not scale with the summed nuclear area",
                 "Approach / wrong strategy"),
            Trap("focus_qc",
                 f"{len(blurred)} fields ({', '.join(blurred)}) are out of focus (extra blur sigma 4.5-6 px): "
                 "two neighbouring series wells just below the marker potency of one compound, two just below "
                 "the count potency of another, one vehicle well. Used unflagged they pull both potencies down "
                 "by about one log10 unit.",
                 "SOP section 2 defines the QC failure; a per-well sharpness metric of channel 1 (edge gradient "
                 "relative to nuclear signal, or high/low spatial-frequency power) is bimodal with these wells "
                 "2-100x below the rest, while sparse wells are not low on a signal-normalised metric",
                 "Verification / skipped QC step"),
            Trap("readout_confound",
                 "Two compounds reduce cell number without reducing the marker per cell (one of them inhibits "
                 "the marker only at >= 12x higher concentration). The field intensity follows cell number, so "
                 "a whole-field readout reports a marker IC50 for a compound that is inactive on the marker and "
                 "one that is 1 log10 unit too low.",
                 "SOP section 1 defines R2 per cell; plotting field intensity against cell number across the "
                 "series of a compound; per-cell and whole-field curves disagree for these compounds",
                 "Understanding / wrong readout"),
        ]
    if level >= 3:
        attacks += [
            Attack("naive_no_plate_model", "everything else correct, but responses are normalised to the plate "
                   "mean of the vehicle wells", 0.75, script=ws + "--plate none", script_files=pipe),
            Attack("naive_no_background", "everything else correct, but the marker per cell is the raw ring "
                   "intensity above the camera offset (no field background)", 0.75,
                   script=ws + "--background none", script_files=pipe),
            Attack("naive_global_background", "everything else correct, but one plate-wide background constant "
                   "(median field background of the vehicle wells) is subtracted from every well", 0.75,
                   script=ws + "--background global", script_files=pipe),
            Attack("naive_cellpose_defaults", "Cellpose nuclei masks, mean marker intensity around them, "
                   "normalised to the vehicle mean, no QC: the professional model without plate-level reasoning",
                   0.30, script=run + "--segmenter cellpose --plate none --background none --qc none",
                   script_files=pipe),
            Attack("naive_default_pipeline", "threshold counting, whole-field marker readout above the camera "
                   "offset, vehicle-mean normalisation, no QC", 0.20,
                   script=run + "--segmenter threshold --plate none --readout wholewell --background none --qc none",
                   script_files=pipe),
        ]
        traps += [
            Trap("plate_position_effects",
                 f"Marker staining falls log-linearly to {plate['stain_end']:.2f}x from column "
                 f"{1 if plate['stain_dir'] >= 0 else 12} to column {12 if plate['stain_dir'] >= 0 else 1} "
                 f"(dose runs along the columns), and the outer wells hold {plate['edge_count']:.2f}x the cells "
                 "of inner wells (the pure cytotoxic and the cytostatic compound sit in rows A and H). "
                 "Normalising to the plate mean of the vehicle wells shifts potencies by 0.25-0.7 log10 units "
                 "and makes a plateau compound look active.",
                 "plate map of the 24 vehicle wells (3 per row, 2 per column) for each readout: marker per cell "
                 "changes monotonically along the columns by a factor ~2.5, cell number is ~45% lower in the "
                 "outer ring; assay notes give the pipetting order and the incubation conditions as facts",
                 "Verification / skipped QC step"),
            Trap("marker_background",
                 "The marker channel has a diffuse background of about a third of the vehicle signal that "
                 "differs between wells, and one marker inhibitor adds fluorescence of 1.3-2x the vehicle "
                 "signal at the top concentration. Without a per-field background the fluorescent compound "
                 "looks inactive and the other marker potencies move by 0.2-0.6 log10 units.",
                 "SOP section 1 excludes signal that does not come from the cells from R2; cell-free regions of channel 2 "
                 "are far above the dark frame and differ between wells; in the series of one compound the "
                 "cell-free background rises with concentration",
                 "Approach / skipped QC step"),
        ]

    n_active_m = sum(c.marker.ic50 is not None for c in cs)
    n_active_c = sum(c.count.ic50 is not None for c in cs)
    assert n_active_m == (4 if level >= 2 else 5) and n_active_c == 3 and len(wells) == 96
    for c in cs:   # activity classes are separated by construction
        for curve in (c.marker, c.count):
            top = float(curve(sim.TOP_DOSE))
            assert (curve.ic50 is not None and sim.DOSES[0] * 3.2 < curve.ic50 <= sim.TOP_DOSE / 5.4 and top < 0.4) \
                or (curve.ic50 is None and top >= 0.59), (c.cls, curve)
    assert not math.isnan(sum(w.n_cells for w in wells))

    return TaskBuild(
        task_id=f"hcs-dose-response-s{seed:04d}-l{level}", family="hcs_dose_response", seed=seed, level=level,
        title="High-content screening plate: well images to compound potencies",
        summary="Segment a two-channel 96-well imaging plate, quantify cell number and marker per cell, apply "
                "image QC and derive per-compound IC50/GI50 values in the presence of plate-level confounders.",
        subdomain="3.3 Cell & Imaging Biology",
        software=["Cellpose", "scikit-image", "scipy", "tifffile", "Python"],
        instruction_md=INSTRUCTION.format(n_wells=len(wells), size=sim.SIZE),
        agent_dockerfile=image.agent_dockerfile(),
        verifier_dockerfile=image.verifier_dockerfile(),
        input_files=input_files,
        solution_files={"solve.sh": ASSETS / "solve.sh", "hcs_pipeline.py": ASSETS / "hcs_pipeline.py",
                        "solve_alt.py": ASSETS / "solve_alt.py"},
        tests_files={"score.py": ASSETS / "score.py", "reference/truth.json": json.dumps(truth, indent=1)},
        traps=traps, attacks=attacks, probes=probes,
        capabilities=["plate-level QC before curve fitting", "choosing the readout the SOP defines",
                      "per-well image QC", "using a pretrained segmentation model as one step of a workflow",
                      "dose-response modelling with an explicit activity criterion", "exact output contracts"],
        judge_explanation_md=JUDGE, pass_threshold=0.9,
        oracle_steps=8, oracle_runtime_min=3.0, cpus=4, memory_mb=8000, storage_mb=10000,
    )
