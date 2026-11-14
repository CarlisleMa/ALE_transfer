# Family `hcs_dose_response` - high-content screening plate, well images to compound potencies

ALE subdomain 3.3 Cell & Imaging Biology. Uses a professional model (Cellpose 3, CPU, weights baked
into the agent image). Code: `pipeline/alegen/families/hcs_dose_response/` (`__init__.py`, `sim.py`,
`image.py`, `assets/`).

## Scenario
An image-analysis scientist receives one 96-well plate: 8 blinded compounds in 8-point 3-fold series
(0.0137-30 uM, one well per concentration), 24 vehicle wells, 8 reference-inhibitor wells, one
320 x 320 px field per well in two channels (nuclei, cytoplasmic marker), flat-field and dark frames,
the plate layout, hand-over notes and the group's SOP. Deliverables: `well_measurements.tsv`
(`well, n_cells, marker_per_cell, qc_pass`) and `dose_response.tsv` (`compound, ic50_marker_uM,
gi50_count_uM`, a number or `inactive`).

## Truth model
Everything is planted before a pixel exists (`sim.py`):
- per compound and readout a response curve relative to vehicle (top = 1): no effect, a logistic
  with a plateau, or a logistic beyond the tested range. For active curves the concentration at
  which the curve equals 0.5 is drawn first and EC50 is derived from it, so the reference potency
  is a planted number. Classes are separated by construction and asserted in `generate()`: actives
  lie between 3.3x the lowest and 1/5.5 of the top concentration and fall below 0.4 at the top
  concentration; inactives are flat, have a plateau >= 0.60, or a potency >= 5x the top concentration;
- per well the number of cells (vehicle level x outer-well factor x count curve x 3% lognormal),
  per cell its marker level (vehicle level x staining factor x marker curve x 3% well noise x 25%
  cell-to-cell), the diffuse marker-channel background, compound fluorescence, extra blur;
- rendering: nuclei as textured ellipses composited by maximum (abutting nuclei form one blob with
  a notch), cytoplasm as a larger ellipse with nearest-cell-wins territories and nuclear exclusion,
  Gaussian PSF, illumination roll-off (30% / 26% at the corners), Poisson and read noise, 16-bit.
`reference/truth.json` holds the planted cell number and focus state per well, the planted 0.5
crossing (or null) per compound and readout, and all plate parameters for audit. No segmentation
and no curve fit enters the reference. The plate is synthetic, so nothing is retrievable from the
internet; compound names are random `CPD-nnnn`.

Compound classes at levels 2-3 (rows and all parameters vary with the seed): potent marker inhibitor
that is cytotoxic 25-120x higher; weak marker inhibitor; pure cytotoxic compound; compound that is
cytotoxic at low dose and inhibits the marker only at >= 12x higher dose; marker inhibitor that
fluoresces in the marker channel (level 3); inactive; partial marker inhibitor (plateau 0.62-0.72);
cytostatic compound (plateau 0.60-0.65) whose marker potency lies beyond the range. Marker: 4 active
/ 4 inactive. Count: 3 active / 5 inactive.

## Decision points and why they are discoverable
None of them is named in the instruction, the SOP or the notes. The SOP defines the measurands
(R1, R2, QC failure = out of focus, response = well / what the same well would have given with
vehicle only, potency = 0.5 crossing, `inactive`); the notes state handling facts (pipetting order,
dry humidity tray, fixed autofocus offset) without consequences.

| trap (level) | what the naive path gets wrong (measured, seeds 1 / 2) | evidence available to the agent |
|---|---|---|
| `touching_nuclei` (2) | global threshold + connected components: median count -47.5% / -47.5%, 7% / 5% of the wells inside tolerance; largest GI50 error 0.70 / 0.58 log10 | overlay of the mask on any vehicle field; object-area histogram with multiples of the single-nucleus area |
| `focus_qc` (2) | 5 blurred fields (2 just below one marker IC50, 2 just below one GI50, 1 vehicle); used unflagged: that IC50 off by 1.19 / 1.09 log10, that GI50 by 1.15 / 1.37 | SOP section 2 defines the failure; a signal-normalised sharpness metric of channel 1 puts these wells 2-100x below all others, sparse wells stay normal |
| `readout_confound` (2) | whole-field intensity: the cytotoxic compound gets a marker IC50, the dual compound is off by 1.02 / 1.12 log10 | SOP section 1 defines R2 per cell; field intensity tracks cell number along the series of those compounds |
| `plate_position_effects` (3) | vehicle-mean normalisation: one wrong activity class per readout on both seeds, largest errors 0.62 / 0.59 (marker) and 0.66 / 0.76 (count), mean credit 0.64 / 0.71 and 0.67 / 0.63 | plate map of the 24 vehicle wells: planted factors: marker per cell changes 2.6x / 2.8x along the columns, cell number is 42% / 45% lower in the outer ring |
| `marker_background` (3) | camera offset only: fluorescent compound `inactive`, largest other IC50 error 0.53 / 0.50; plate-wide constant: fluorescent compound `inactive`, other IC50s up to 0.21 / 0.12 off | cell-free regions of channel 2 are far above the dark frame, differ between wells and rise with concentration in one series |

Things that look like traps but are not (measured): skipping the flat-field correction changes
nothing (seeds 1 and 2: all entries right, largest error 0.065); Cellpose with its default
per-image normalisation loses nuclei in fields with fewer than ~20 cells (5 / 3 of 91 wells outside
the count tolerance on seeds 1 / 2, one field with no cell found), everything else unchanged; scored
offline that is 0.945 / 0.967 of the wells inside tolerance, i.e. still full count credit.

## Output contract and scoring
Gates: files parse; every well exactly once; `n_cells` non-negative integer; `qc_pass` 0/1; every
test compound exactly once; potency = positive number or `inactive`. Well ids `A1`/`A01`, extra
columns, any row order and rows for DMSO / REF-INH are accepted.

| weight | component | metric | full | zero |
|---|---|---|---|---|
| 0.20 | cell_counts | fraction of in-focus wells with abs(n - planted) <= max(3, 10%) | >= 0.90 | <= 0.50 |
| 0.10 | qc_flags | MCC of `qc_pass = 0` vs out-of-focus wells | >= 0.90 | <= 0.50 |
| 0.45 | marker_potency | mean per-compound credit (class must match; active: abs log10 error <= 0.20 -> 1, >= 0.40 -> 0), then ramp 0.8 -> 1.0 | 1.0 | <= 0.8 |
| 0.25 | count_potency | same | 1.0 | <= 0.8 |

Pass threshold 0.9. `marker_per_cell` is parsed but its absolute value is not scored (it depends on
the legitimate choice of the cytoplasm region).

### Measured values
Seed 1, level 3 (`results/work/hcs-dose-response-s0001-l3/*/logs/details.json`). Count column: fraction of in-focus wells inside tolerance / median relative count error. Potency columns: compounds with the wrong activity class / largest abs log10 error among compounds called active on both sides / mean per-compound credit (8 compounds).

| submission | score | counts | QC MCC | marker IC50 | count GI50 |
|---|---|---|---|---|---|
| oracle (Cellpose nuclei model, plate-wide intensity scale) | 1.000 | 1.000 / +0.000 | 1.00 | 0 / 0.052 / 1.000 | 0 / 0.029 / 1.000 |
| independent solution (LoG-seeded watershed, no Cellpose) | 1.000 | 1.000 / +0.000 | 1.00 | 0 / 0.069 / 1.000 | 0 / 0.057 / 1.000 |
| reference pipeline with distance-transform watershed | 1.000 | 1.000 / -0.015 | 1.00 | 0 / 0.055 / 1.000 | 0 / 0.065 / 1.000 |
| truth, counts +-7%, potencies +-0.15, one false QC flag | 1.000 | 1.000 / +0.066 | 0.91 | 0 / 0.150 / 1.000 | 0 / 0.150 / 1.000 |
| naive: threshold + connected components | 0.550 | 0.066 / -0.475 | 1.00 | 0 / 0.056 / 1.000 | 0 / 0.696 / 0.714 |
| naive: no focus QC | 0.463 | 1.000 / -0.015 | 0.00 | 0 / 1.193 / 0.875 | 0 / 1.151 / 0.875 |
| naive: whole-field marker readout | 0.550 | 1.000 / -0.015 | 1.00 | 1 / 1.016 / 0.750 | 0 / 0.065 / 1.000 |
| naive: vehicle-mean normalisation | 0.300 | 1.000 / -0.015 | 1.00 | 1 / 0.622 / 0.637 | 1 / 0.658 / 0.666 |
| naive: camera offset only | 0.550 | 1.000 / -0.015 | 1.00 | 1 / 0.533 / 0.514 | 0 / 0.065 / 1.000 |
| naive: plate-wide background constant | 0.703 | 1.000 / -0.015 | 1.00 | 1 / 0.211 / 0.868 | 0 / 0.065 / 1.000 |
| Cellpose masks, no plate model, no background, no QC | 0.200 | 1.000 / +0.000 | 0.00 | 3 / 0.451 / 0.500 | 1 / 1.077 / 0.547 |
| fully naive default pipeline | 0.000 | 0.066 / -0.475 | 0.00 | 4 / 0.000 / 0.500 | 0 / 1.380 / 0.750 |
| prior only (all inactive, 100 cells) | 0.000 | 0.022 / -0.099 | 0.00 | 4 / 0.000 / 0.500 | 3 / 0.000 / 0.625 |
| truth, all potencies +0.45 log10 | 0.300 | 1.000 / +0.000 | 1.00 | 0 / 0.450 / 0.500 | 0 / 0.450 / 0.625 |
| truth, one marker IC50 +0.45 log10 | 0.719 | 1.000 / +0.000 | 1.00 | 0 / 0.450 / 0.875 | 0 / 0.000 / 1.000 |
| truth, all counts -20% | 0.800 | 0.044 / -0.200 | 1.00 | 0 / 0.000 / 1.000 | 0 / 0.000 / 1.000 |

Seed 2, level 3 (`results/work/hcs-dose-response-s0002-l3/*/logs/details.json`). Count column: fraction of in-focus wells inside tolerance / median relative count error. Potency columns: compounds with the wrong activity class / largest abs log10 error among compounds called active on both sides / mean per-compound credit (8 compounds).

| submission | score | counts | QC MCC | marker IC50 | count GI50 |
|---|---|---|---|---|---|
| oracle (Cellpose nuclei model, plate-wide intensity scale) | 1.000 | 1.000 / +0.000 | 1.00 | 0 / 0.027 / 1.000 | 0 / 0.063 / 1.000 |
| independent solution (LoG-seeded watershed, no Cellpose) | 1.000 | 1.000 / +0.000 | 1.00 | 0 / 0.068 / 1.000 | 0 / 0.066 / 1.000 |
| reference pipeline with distance-transform watershed | 1.000 | 1.000 / -0.016 | 1.00 | 0 / 0.033 / 1.000 | 0 / 0.070 / 1.000 |
| truth, counts +-7%, potencies +-0.15, one false QC flag | 1.000 | 1.000 / -0.064 | 0.91 | 0 / 0.150 / 1.000 | 0 / 0.150 / 1.000 |
| naive: threshold + connected components | 0.550 | 0.055 / -0.475 | 1.00 | 0 / 0.030 / 1.000 | 0 / 0.575 / 0.739 |
| naive: no focus QC | 0.463 | 1.000 / -0.016 | 0.00 | 0 / 1.093 / 0.875 | 0 / 1.374 / 0.875 |
| naive: whole-field marker readout | 0.550 | 1.000 / -0.016 | 1.00 | 1 / 1.124 / 0.750 | 0 / 0.070 / 1.000 |
| naive: vehicle-mean normalisation | 0.300 | 1.000 / -0.016 | 1.00 | 1 / 0.588 / 0.714 | 1 / 0.762 / 0.633 |
| naive: camera offset only | 0.550 | 1.000 / -0.016 | 1.00 | 1 / 0.495 / 0.647 | 0 / 0.070 / 1.000 |
| naive: plate-wide background constant | 0.719 | 1.000 / -0.016 | 1.00 | 1 / 0.116 / 0.875 | 0 / 0.070 / 1.000 |
| Cellpose masks, no plate model, no background, no QC | 0.200 | 1.000 / +0.000 | 0.00 | 3 / 0.358 / 0.526 | 1 / 1.341 / 0.503 |
| fully naive default pipeline | 0.000 | 0.055 / -0.475 | 0.00 | 4 / 0.000 / 0.500 | 2 / 0.037 / 0.750 |
| prior only (all inactive, 100 cells) | 0.000 | 0.033 / +0.205 | 0.00 | 4 / 0.000 / 0.500 | 3 / 0.000 / 0.625 |
| truth, all potencies +0.45 log10 | 0.300 | 1.000 / +0.000 | 1.00 | 0 / 0.450 / 0.500 | 0 / 0.450 / 0.625 |
| truth, one marker IC50 +0.45 log10 | 0.719 | 1.000 / +0.000 | 1.00 | 0 / 0.450 / 0.875 | 0 / 0.000 / 1.000 |
| truth, all counts -20% | 0.800 | 0.033 / -0.200 | 1.00 | 0 / 0.000 / 1.000 | 0 / 0.000 / 1.000 |

Reading of the bands:
- counts: correct pipelines have 100% of the in-focus wells inside max(3, 10%) (worst single well
  over six seeds and two classical pipelines: 16 of 18 cells; worst dense well 108 of 118);
  threshold counting has 7% / 5%. Band 0.5-0.9.
- potency: correct pipelines are within 0.07 log10 on the gated seeds (0.106 over six seeds and two
  classical pipelines). Every ignored trap produces a wrong activity class and/or an error of
  0.50-1.37 log10 in at least one entry. Full credit <= 0.20, zero >= 0.40. The 0.12-0.21 shifts of
  the other IC50s under a plate-wide background constant fall partly inside the tolerance; that
  variant fails because of the fluorescent compound alone.
- QC: exact for all correct pipelines (MCC 1.0); one false flag gives MCC 0.908 (still full),
  no flags gives 0.

## Levels
(Level 4 has its own section at the end of this file.)
1: uniform plate, separated nuclei, no blurred field, no compound fluorescence, cytotoxic compounds
inhibit the marker at lower dose than they kill (whole-field and per-cell readouts agree).
2: + touching nuclei, out-of-focus fields, readout confound. 3: + plate-position effects and
well-specific background / fluorescent compound. Gate runs: seed 1 level 3 with Harbor (oracle 1.0,
nop 0.0); seed 2 level 3, seed 1 level 1 and seed 1 level 2 without Harbor
(`results/validation/hcs-dose-response-*.json`). Level 2, seed 1: oracle, independent solution and
watershed pipeline 1.0; threshold counting 0.55, no focus QC 0.4625, whole-field readout 0.55.
Level 1, seed 1: the three correct pipelines 1.0 (largest potency error 0.057), prior-only 0.0.

## Image
`python:3.12-slim` base via `dockerlib.agent_dockerfile` (conda: python 3.12, numpy 2.0.2, scipy
1.14.1, pandas 2.2.3), pip: torch 2.5.1+cpu, cellpose 3.1.1.2, scikit-image 0.25.2, tifffile
2025.3.30, matplotlib 3.10.1. Weights fetched at build time from `https://www.cellpose.org/models/`
into `/opt/cellpose_models` (`CELLPOSE_LOCAL_MODELS_PATH`); the build verifies the sha256 sums pinned
in `image.py` (files as served on 2026-11-11):
- `nucleitorch_0` 89ca45e4a45048d5010d29621466b1274b91b3dcb9714dce9f9e90a9a8671303
- `size_nucleitorch_0.npy` a79107c9f284b569dd65f3a7363014d613a3403241c062822f07b6f026b4d3ef
- `cyto3` 2dc3087a8abd7da46d1ab0ddd5824639933cc3ff63b382af3fa1939a392db93c
- `size_cyto3.npy` 87f91feb48019b4dfb0400e7bbb89aca5873368686815aa7b4431c0faee97fc5

Host venv: `tifffile` (pure Python) was added for writing the TIFFs; the renderer needs numpy only.
Size: 2.82 GB as reported by `docker images` on the build host
(`docker images` on this host reports more than the file system holds: the untrimmed tools-only
image was 2.97 GB by `docker images` and 2.15 GB by `du -sx /`; removing PyTorch's C++ test binaries
and headers, 176 MB on disk, took the task image from 3.05 GB to 2.82 GB). Oracle runtime
2.6-4.5 min in the four gate runs (2 threads, shared host under load; one early version that used
the Cellpose size model and rescaling took 14 min); the classical pipelines take 6-10 s.

## Known limitations and things I am unsure about
- **No agent pilot was run.** Whether the obstacles are hard for a model is untested. The SOP has to
  define QC failure as "out of focus" and the response as relative to what the same well would give
  with vehicle; a model that reads these two sentences carefully is pointed at "compute a focus
  metric" and "position matters". What remains undisclosed is which wells, that the two readouts
  have different plate structure (gradient vs outer ring), that the background is well-specific, and
  the readout confound.
- **Steep potency scoring.** One wrong compound in a readout costs 0.28 (marker) or 0.16 (count), so
  a pass needs all 16 entries right. This is what makes single-compound traps bite, but it also
  means little partial credit inside the potency components.
- **`naive_global_background` has little margin**: 0.703 / 0.719 against a 0.75 ceiling. With
  exactly one compound wrong the score is 0.71875 by construction; it cannot exceed that unless the
  fluorescent compound is recovered, which the planted fluorescence (1.3-2x the vehicle signal at
  the top dose, potency 2.5-5.5 uM) rules out.
- **Plate model is a judgement.** The reference fits log(readout) ~ column + row + outer-ring
  indicator on the vehicle wells. Column-matched or nearest-vehicle normalisation handles the
  marker gradient but not the outer ring for cell number, and a smooth quadratic surface
  approximates the ring poorly. I did not measure those alternatives; an analyst who chooses one of
  them without looking at residuals may fail on the compounds in rows A / H.
- **Out-of-focus fields are strongly blurred** (extra sigma 4.5-6 px). Detection is easy once a
  normalised sharpness metric exists. I did not measure an un-normalised metric (e.g. Laplacian
  variance), which I expect to flag sparse wells as well.
- **Images are synthetic.** No debris, no mitotic or apoptotic morphology, nuclei lie fully inside
  the field, cell size does not change with treatment, the defocus is a Gaussian blur. Cellpose
  and classical methods both count these images almost exactly, which real plates would not allow.
- **Cellpose robustness** was measured on seeds 1 and 2 at level 3 and seed 1 at levels 1 and 2 only
  (classical pipelines on seeds 1-6 at level 3).
- Tolerances were derived with 3% planted well-to-well noise. With the first setting (5% / 4%) two
  correct pipelines were 0.13-0.14 log10 off on one compound; the noise was lowered rather than
  the tolerance widened, because the trap-induced shifts start at 0.21.
- The notes state the pipetting order and the dry humidity tray. They are facts a hand-over would
  contain and are not linked to any consequence, but they are hints.

# Level 4 - counts that have to be validated (Q + P)

Code: `level4.py` (instruction, SOP v5, battery), `sim4.py` (image model), `assets/hcs_l4_voting.py`
(reference solution), `assets/hcs_l4_cellpose.py` (second solution), `assets/score_l4.py`,
`assets/solve_l4.sh`, `assets/l3_validation.py` (battery helper). `generate()` only dispatches to
`level4.generate_l4` for `level >= 4`; levels 1-3 use the old code path and the old files
(`make --seed 1 --level 3 --out <scratch> --no-validate`, `diff -r` against
`tasks/hcs-dose-response-s0001-l3`: no difference, checked after the last code change).
Level 4 uses its own random stream (`884000 + seed`), so the level-4 plate of a seed is not the
level-3 plate of that seed.

## Why (pilot finding)
The pilot model found Cellpose slow on CPU, counted with a classical segmentation and scored 1.0 in
14 min 50 s. That was the predicted outcome: at level 3 every reasonable counting method is exact,
so nothing made the analyst look at the counts. Level 4 changes two things: the images are such
that the level-3 methods miscount, and the inputs contain the control (manual counts of eight
fields) with which the analyst can see that.

## What was added

| element | what is planted | tied to the pilot finding |
|---|---|---|
| `Q_validated_counts` | images on which the counting method matters, manual point annotation of 8 fields, an acceptance criterion in the SOP, a new deliverable `validation.tsv` | "switched to classical segmentation and still scored full": on these images the level-3 classical route, the level-3 second solution and Cellpose with the level-3 settings all fail the stated criterion and the count component (table below) |
| `P_potency_shape` | a deep partial count response and a marker signal per cell that rises as cells are lost | every level-3 decision point that a note or a routine QC step reveals was handled; these two shapes are only handled by fitting the SOP's model instead of the textbook default |

Discoverability (as written into `Trap.discoverable_by`): Q - SOP section 2 states the criterion and
`annotations/` holds the control; comparing any counting procedure with the manual counts, or
overlaying its objects on the annotated points, shows which nuclei are lost and which objects are
counted wrongly. P - SOP section 5 fixes the upper asymptote, defines the potency as the 0.5
crossing and lists the `inactive` cases; a plot of the responses of each series shows the plateau at
about a third and the rising series.

Image model (`sim4.py`, all planted before rendering; parameters in `truth.json` under
`plate.image_model`):
- **Overlapping nuclei without an intensity dip.** Colony growth with a minimal centre distance of
  9.5 px (x the size factor of the well) for nuclei of 12-15 x 10-12 px; 55% of the colony steps
  place the new nucleus 9.5-12 px from an existing one. Flat-topped profile, composited by maximum.
- **Chromatin.** Two texture scales (sigma 1.0 and 2.4 px, amplitudes 0.20 and 0.18) and 0-2 dark
  nucleoli per nucleus (depth 0.3-0.6, sigma 0.8-1.3 px). Synthetic. Real nuclear crops were not used
  (see "Dropped").
- **Fragments.** Bright compact objects, radius 1.3-2.5 px (3.4-6.5 um across), 1.3-2.8x the nuclear
  level, Poisson number with mean 0.12 x seeded cells + 0.40 x lost cells. Fragment centres lie at least
  3.5 px outside every nuclear outline (so the edge of a large fragment can come within 1 px of a nucleus),
  7 px from the image edge and 8 px from the next fragment centre. The SOP does not count them.
  Rendered: 0.117 / 0.127 fragments per cell in the median vehicle well (seeds 1 / 2), up to 3.2 per
  cell at the top dose of the cytotoxic compound.
- **Field border.** Nuclei are cut by the image edge. SOP rule (counting frame): cut by the top or
  left edge = counted, cut by the bottom or right edge = not counted. By construction a nucleus either
  stays >= 3 px inside an edge or extends >= 2.5 px beyond it, and a cut nucleus shows >= 5 px.
  Rendered: 1.8 nuclei per field lie outside the counting frame, 2.1% of the counted cells (both seeds).
- **Nuclear size changes with treatment.** The potent marker inhibitor and the cytostatic compound
  enlarge the nuclei with their count effect, up to 1.35-1.5x and 1.3-1.45x (linear). Enlargement
  only: no shrinking phenotype, so that the smallest nucleus (8.8 px = 11.5 um in its shortest
  dimension) stays far from the largest fragment and the SOP's definition has one reading.
- **Cytoplasm** wider and flatter than at level 3 (see worklog L4-5: with the level-3 profile the ring
  readout depended on the packing of the colony).
- **Manual annotation.** Eight in-focus fields: two inner and one outer vehicle well, the top dose
  and the dose nearest the GI50 of the cytotoxic compound, the top dose of the potent inhibitor and
  of the cytostatic compound (enlarged nuclei), the second-highest dose of the dual compound. One
  point per nucleus that the SOP counts (planted centre, 0.8 px jitter). The manual count equals the
  planted count.

Curve shapes (`sim4.reshape_compounds`):
- the count response of the potent marker inhibitor has a plateau of 0.28-0.34 and a GI50 of
  1.2-3.0 uM: the SOP potency (0.5 crossing) lies 0.32 (seed 1, plateau 0.30) / 0.36 (seed 2, plateau
  0.32) log10 units above the inflection point;
- the marker per cell of the pure cytotoxic compound is multiplied by 1 + k x (fraction of cells
  lost), k = 0.5-0.8 (0.76 on seed 1, 0.53 on seed 2). Planted as a real property of the surviving
  cells. Under SOP section 5 the compound is `inactive` on the marker.

SOP v5 (differences to v4): section 1 defines what is counted (fragments by size with a stated gap:
largest fragment about 6.5 um, smallest nucleus 11 um; counting frame; overlapping nuclei are
separate cells; R2 over the counted cells); new section 2 states the validation and its criterion,
`|procedure count - manual count| <= max(3 cells, 10% of the manual count)` on every annotated
field, and that manual counts are not entered into the well table; section 5 adds "a readout that
does not fall at all" to the list of `inactive` cases. No method, threshold of a method or
tolerance of the verifier is named. The notes gain two facts: how the fields were annotated, and
that wells were not washed before fixation.

## Output contract and scoring at level 4
As level 3 plus `validation.tsv` (`field`, `manual_count`, `pipeline_count`; one row per annotated
field, no other rows; gate).

| weight | component | metric | full | zero |
|---|---|---|---|---|
| 0.25 | cell_counts | fraction of in-focus wells **without annotation** with abs(n - planted) <= max(3, 10%) | >= 0.85 | <= 0.65 |
| 0.05 | validation | fraction of annotated fields whose row is truthful (manual count = number of points, pipeline count = `n_cells` of that well) and inside the SOP criterion | all | <= half |
| 0.05 | qc_flags | MCC of `qc_pass = 0` vs out-of-focus wells | >= 0.90 | <= 0.50 |
| 0.45 | marker_potency | as level 3 | 1.0 | <= 0.8 |
| 0.20 | count_potency | as level 3 | 1.0 | <= 0.8 |

Pass threshold 0.9.

## Measured values (level 4)
Gate runs `make hcs_dose_response --seed 1 --level 4 --keep-work` (with Harbor) and
`--seed 2 --level 4 --no-harbor --keep-work`; raw values in
`results/work/hcs-dose-response-s000{1,2}-l4/*/logs/details.json`, gate records in
`results/validation/hcs-dose-response-s000{1,2}-l4.json`.

Columns. *Hidden wells*: the 83 in-focus wells without annotation: fraction with
abs(n - planted) <= max(3, 10%) / median / minimum / maximum relative count error. *Annotated
fields*: number of the 8 fields inside the SOP criterion, then the relative count error of each
field in alphabetical order of the wells (seed 1: A06 A11 A12 C12 D10 E05 G07 H12; seed 2: A08 A12
C10 D01 D12 G03 G10 H12; a field with 8-10 cells is inside the criterion at +12%, the limit there is
3 cells). Potency columns as in the level-3 tables: compounds with the wrong activity class /
largest abs log10 error among compounds called active on both sides / mean per-compound credit.
Rows "by hand" were produced with the same verifier files outside the gate (they are not part of the
battery).

Seed 1 (`hcs-dose-response-s0001-l4`):

| submission | score | hidden wells: fraction inside / median / min / max rel. error | annotated fields inside criterion; rel. error per field | QC MCC | marker IC50: wrong class / max err / credit | count GI50: wrong class / max err / credit |
|---|---|---|---|---|---|---|
| oracle | 1.000 | 1.000 / -0.008 / -0.062 / +0.067 | 8/8; +0% -8% +0% +0% +0% -1% -2% +0% | 1.00 | 0 / 0.032 / 1.000 | 0 / 0.065 / 1.000 |
| probe: planted_truth | 1.000 | 1.000 / +0.000 / +0.000 / +0.000 | 8/8; +0% +0% +0% +0% +0% +0% +0% +0% | 1.00 | 0 / 0.000 / 1.000 | 0 / 0.000 / 1.000 |
| probe: within_tolerance | 1.000 | 1.000 / -0.066 / -0.078 / +0.091 | 8/8; +6% +8% +10% -6% -6% -7% +7% -8% | 0.91 | 0 / 0.150 / 1.000 | 0 / 0.150 / 1.000 |
| probe: independent_solution | 1.000 | 0.964 / -0.025 / -0.124 / +0.020 | 8/8; +0% -6% +0% +0% +0% -1% -7% -3% | 1.00 | 0 / 0.044 / 1.000 | 0 / 0.036 / 1.000 |
| attack: empty_output | 0.000 | gate: required output file missing: well_measurements.tsv | | | | |
| attack: inputs_copied | 0.000 | gate: required output file missing: well_measurements.tsv | | | | |
| attack: zero_byte_files | 0.000 | gate: well_measurements.tsv is empty | | | | |
| attack: header_only | 0.000 | gate: well_measurements.tsv lists every well exactly once | | | | |
| attack: l3_reference | 0.647 | 0.434 / -0.113 / -0.273 / +0.073 | 2/8; -15% -21% -10% +22% -10% -15% -19% -13% | 1.00 | 0 / 0.074 / 1.000 | 0 / 0.285 / 0.947 |
| attack: l3_reference_classical | 0.309 | 0.494 / -0.068 / -0.217 / +1.273 | 3/8; +9% -15% +80% +50% +65% -3% -8% +26% | 1.00 | 0 / 0.336 / 0.915 | 0 / 1.425 / 0.635 |
| attack: l3_independent_solution | 0.050 | 0.554 / +0.092 / -0.017 / +2.546 | 3/8; +38% +9% +300% +172% +184% +8% +2% +129% | 1.00 | 1 / 0.689 / 0.750 | 3 / 0.000 / 0.625 |
| attack: counts_scaled_to_annotation | 0.309 | 0.229 / -0.119 / -0.265 / +1.182 | 2/8; +3% -20% +70% +39% +55% -8% -13% +18% | 1.00 | 0 / 0.336 / 0.915 | 0 / 1.425 / 0.635 |
| attack: counts_scaled_to_annotation_log | 0.050 | 0.036 / -0.250 / -0.328 / +1.454 | 1/8; -6% -24% +180% +89% +97% -26% -30% +58% | 1.00 | 1 / 0.689 / 0.750 | 3 / 0.000 / 0.625 |
| attack: manual_counts_copied | 0.359 | 0.494 / -0.068 / -0.217 / +1.273 | 8/8; +0% +0% +0% +0% +0% +0% +0% +0% | 1.00 | 0 / 0.336 / 0.915 | 0 / 1.425 / 0.635 |
| attack: naive_fragments_counted | 0.312 | 0.879 / +0.024 / -0.074 / +1.546 | 5/8; +18% -8% +230% +0% +129% +6% +2% +0% | 1.00 | 1 / 0.556 / 0.750 | 2 / 0.073 / 0.750 |
| attack: naive_free_plateau_fit | 0.350 | 1.000 / -0.008 / -0.062 / +0.067 | 8/8; +0% -8% +0% +0% +0% -1% -2% +0% | 1.00 | 2 / 0.100 / 0.750 | 1 / 0.438 / 0.585 |
| attack: naive_no_focus_qc | 0.544 | 1.000 / -0.008 / -0.062 / +0.067 | 8/8; +0% -8% +0% +0% +0% -1% -2% +0% | 0.00 | 0 / 1.058 / 0.875 | 0 / 1.331 / 0.875 |
| attack: naive_whole_well_readout | 0.550 | 1.000 / -0.008 / -0.062 / +0.067 | 8/8; +0% -8% +0% +0% +0% -1% -2% +0% | 1.00 | 1 / 0.856 / 0.750 | 0 / 0.065 / 1.000 |
| attack: naive_no_plate_model | 0.350 | 1.000 / -0.008 / -0.062 / +0.067 | 8/8; +0% -8% +0% +0% +0% -1% -2% +0% | 1.00 | 1 / 0.477 / 0.750 | 1 / 0.823 / 0.555 |
| attack: naive_no_background | 0.550 | 1.000 / -0.008 / -0.062 / +0.067 | 8/8; +0% -8% +0% +0% +0% -1% -2% +0% | 1.00 | 1 / 0.343 / 0.766 | 0 / 0.065 / 1.000 |
| attack: naive_global_background | 0.719 | 1.000 / -0.008 / -0.062 / +0.067 | 8/8; +0% -8% +0% +0% +0% -1% -2% +0% | 1.00 | 1 / 0.165 / 0.875 | 0 / 0.065 / 1.000 |
| attack: prior_only | 0.000 | 0.024 / -0.115 / -0.213 / +8.091 | 0/8; +194% +52% +900% +456% +223% -12% -21% +163% | 0.00 | 4 / 0.000 / 0.500 | 3 / 0.000 / 0.625 |
| attack: all_potencies_beyond_tolerance | 0.350 | 1.000 / +0.000 / +0.000 / +0.000 | 8/8; +0% +0% +0% +0% +0% +0% +0% +0% | 1.00 | 0 / 0.450 / 0.500 | 0 / 0.450 / 0.625 |
| attack: one_potency_beyond_tolerance | 0.719 | 1.000 / +0.000 / +0.000 / +0.000 | 8/8; +0% +0% +0% +0% +0% +0% +0% +0% | 1.00 | 0 / 0.450 / 0.875 | 0 / 0.000 / 1.000 |
| attack: counts_beyond_tolerance | 0.700 | 0.024 / -0.200 / -0.212 / -0.182 | 1/8; -21% -20% -20% -22% -19% -20% -20% -21% | 1.00 | 0 / 0.000 / 1.000 | 0 / 0.000 / 1.000 |
| by hand: counting frame ignored (not in the battery) | 1.000 | 0.988 / +0.000 / -0.042 / +0.133 | 8/8; +0% -6% +0% +11% +0% +0% +1% +3% | 1.00 | 0 / 0.032 / 1.000 | 0 / 0.041 / 1.000 |
| by hand: level-4 Cellpose pipeline, level-3 settings, SOP rules (not in the battery) | 0.689 | 0.349 / -0.123 / -0.273 / +0.000 | 3/8; -15% -23% -10% +0% -10% -16% -21% -13% | 1.00 | 0 / 0.072 / 1.000 | 0 / 0.217 / 0.989 |

Seed 2 (`hcs-dose-response-s0002-l4`):

| submission | score | hidden wells: fraction inside / median / min / max rel. error | annotated fields inside criterion; rel. error per field | QC MCC | marker IC50: wrong class / max err / credit | count GI50: wrong class / max err / credit |
|---|---|---|---|---|---|---|
| oracle | 1.000 | 1.000 / -0.015 / -0.075 / +0.100 | 8/8; +0% +12% -1% -3% +0% +0% +0% +0% | 1.00 | 0 / 0.070 / 1.000 | 0 / 0.037 / 1.000 |
| probe: planted_truth | 1.000 | 1.000 / +0.000 / +0.000 / +0.000 | 8/8; +0% +0% +0% +0% +0% +0% +0% +0% | 1.00 | 0 / 0.000 / 1.000 | 0 / 0.000 / 1.000 |
| probe: within_tolerance | 1.000 | 1.000 / -0.062 / -0.100 / +0.077 | 8/8; +6% -12% -7% -7% +5% -7% +6% +8% | 0.91 | 0 / 0.150 / 1.000 | 0 / 0.150 / 1.000 |
| probe: independent_solution | 1.000 | 1.000 / -0.025 / -0.094 / +0.042 | 8/8; +0% +12% +1% -6% +10% -5% +0% +0% | 1.00 | 0 / 0.087 / 1.000 | 0 / 0.061 / 1.000 |
| attack: empty_output | 0.000 | gate: required output file missing: well_measurements.tsv | | | | |
| attack: inputs_copied | 0.000 | gate: required output file missing: well_measurements.tsv | | | | |
| attack: zero_byte_files | 0.000 | gate: well_measurements.tsv is empty | | | | |
| attack: header_only | 0.000 | gate: well_measurements.tsv lists every well exactly once | | | | |
| attack: l3_reference | 0.700 | 0.518 / -0.100 / -0.256 / +0.015 | 4/8; -9% -38% -13% -19% +0% -11% -18% +2% | 1.00 | 0 / 0.082 / 1.000 | 0 / 0.148 / 1.000 |
| attack: l3_reference_classical | 0.226 | 0.530 / -0.078 / -0.202 / +0.900 | 3/8; +21% +88% -4% -7% +24% -8% +55% +15% | 1.00 | 0 / 0.395 / 0.878 | 0 / 0.914 / 0.625 |
| attack: l3_independent_solution | 0.050 | 0.518 / +0.097 / -0.058 / +2.400 | 3/8; +79% +350% +9% -1% +210% +3% +164% +95% | 1.00 | 1 / 0.730 / 0.750 | 2 / 1.345 / 0.625 |
| attack: counts_scaled_to_annotation | 0.226 | 0.277 / -0.122 / -0.246 / +0.800 | 2/8; +15% +75% -9% -12% +19% -13% +45% +10% | 1.00 | 0 / 0.395 / 0.878 | 0 / 0.914 / 0.625 |
| attack: counts_scaled_to_annotation_log | 0.050 | 0.012 / -0.250 / -0.350 / +1.300 | 0/8; +21% +212% -26% -32% +114% -29% +82% +32% | 1.00 | 1 / 0.730 / 0.750 | 2 / 1.345 / 0.625 |
| attack: manual_counts_copied | 0.276 | 0.530 / -0.078 / -0.202 / +0.900 | 8/8; +0% +0% +0% +0% +0% +0% +0% +0% | 1.00 | 0 / 0.395 / 0.878 | 0 / 0.914 / 0.625 |
| attack: naive_fragments_counted | 0.481 | 0.879 / +0.024 / -0.071 / +2.000 | 5/8; +52% +150% +0% +1% +0% +3% +79% +0% | 1.00 | 0 / 0.558 / 0.875 | 0 / 0.893 / 0.750 |
| attack: naive_free_plateau_fit | 0.350 | 1.000 / -0.015 / -0.075 / +0.100 | 8/8; +0% +12% -1% -3% +0% +0% +0% +0% | 1.00 | 2 / 0.103 / 0.750 | 1 / 0.490 / 0.636 |
| attack: naive_no_focus_qc | 0.544 | 1.000 / -0.015 / -0.075 / +0.100 | 8/8; +0% +12% -1% -3% +0% +0% +0% +0% | 0.00 | 0 / 1.094 / 0.875 | 0 / 1.114 / 0.875 |
| attack: naive_whole_well_readout | 0.550 | 1.000 / -0.015 / -0.075 / +0.100 | 8/8; +0% +12% -1% -3% +0% +0% +0% +0% | 1.00 | 1 / 0.958 / 0.750 | 0 / 0.037 / 1.000 |
| attack: naive_no_plate_model | 0.513 | 1.000 / -0.015 / -0.075 / +0.100 | 8/8; +0% +12% -1% -3% +0% +0% +0% +0% | 1.00 | 1 / 0.204 / 0.872 | 1 / 0.726 / 0.616 |
| attack: naive_no_background | 0.550 | 1.000 / -0.015 / -0.075 / +0.100 | 8/8; +0% +12% -1% -3% +0% +0% +0% +0% | 1.00 | 1 / 0.439 / 0.731 | 0 / 0.037 / 1.000 |
| attack: naive_global_background | 0.550 | 1.000 / -0.015 / -0.075 / +0.100 | 8/8; +0% +12% -1% -3% +0% +0% +0% +0% | 1.00 | 2 / 0.192 / 0.750 | 0 / 0.037 / 1.000 |
| attack: prior_only | 0.000 | 0.012 / -0.115 / -0.194 / +9.000 | 0/8; +203% +1150% -14% +47% +376% -19% +203% +150% | 0.00 | 4 / 0.000 / 0.500 | 3 / 0.000 / 0.625 |
| attack: all_potencies_beyond_tolerance | 0.350 | 1.000 / +0.000 / +0.000 / +0.000 | 8/8; +0% +0% +0% +0% +0% +0% +0% +0% | 1.00 | 0 / 0.450 / 0.500 | 0 / 0.450 / 0.625 |
| attack: one_potency_beyond_tolerance | 0.719 | 1.000 / +0.000 / +0.000 / +0.000 | 8/8; +0% +0% +0% +0% +0% +0% +0% +0% | 1.00 | 0 / 0.450 / 0.875 | 0 / 0.000 / 1.000 |
| attack: counts_beyond_tolerance | 0.700 | 0.024 / -0.200 / -0.222 / -0.176 | 1/8; -21% -25% -20% -21% -19% -20% -21% -20% | 1.00 | 0 / 0.000 / 1.000 | 0 / 0.000 / 1.000 |
| by hand: counting frame ignored (not in the battery) | 1.000 | 1.000 / +0.008 / -0.060 / +0.100 | 8/8; +0% +12% +0% +3% +5% +1% +3% +2% | 1.00 | 0 / 0.070 / 1.000 | 0 / 0.069 / 1.000 |

Who is who: `oracle` = outline-voting pipeline without Cellpose (`hcs_l4_voting.py`);
`probe:independent_solution` = Cellpose `nuclei` model, settings tried on the annotated fields,
12 px / flow check off chosen on both seeds (`hcs_l4_cellpose.py`); `attack:l3_reference` = the
level-3 oracle unchanged (Cellpose, level-3 settings); `attack:l3_reference_classical` = the level-3
watershed route unchanged; `attack:l3_independent_solution` = the level-3 LoG solution unchanged;
all `naive_*` rows = the level-4 oracle with exactly one switch.

Reading of the numbers:
- **Counts, hidden wells.** Reference solution 1.000 / 1.000, second solution 0.964 / 1.000. The
  three level-3 methods 0.434-0.554 / 0.518-0.530, after a plate-wide factor 0.229 and 0.036 / 0.277
  and 0.012. Band of the component: zero <= 0.65, full >= 0.85.
- **Counts, annotated fields.** Both solutions 8/8 on both seeds (worst field relative to the limit:
  reference 0.76 / 0.33, second solution 0.71 / 0.67). Level-3 methods 2-4 of 8.
  The errors of the level-3 methods have both signs within one plate. Over the 16 annotated fields of the
  two seeds: Cellpose with the level-3 settings -9..-38% on 13 fields, 0..+2% on two, +22% on one field with
  enlarged nuclei; the watershed route -15..-3% on the six untreated fields and +9..+88% on the ten treated
  ones; the LoG route -1..+9% on the untreated and +38..+350% on the treated fields.
- **`l3_reference` (0.647 / 0.700)** loses the count and the validation component and nothing else on
  seed 2 (all 16 potency entries right, largest error 0.148); on seed 1 one GI50 is 0.285 off. 0.70 is
  the most a submission with failed counts can get; the ceiling is 0.75.
- **Fragments counted** (0.312 / 0.481): 3 of 8 annotated fields fail by +18..+230% and the potencies
  break (seed 1: one marker and two count entries with the wrong activity class; seed 2: errors up to
  0.56 and 0.89). On the hidden wells this variant still gets full count credit (0.879 / 0.879 inside,
  median error +2.4%) although the median vehicle well holds 12% fragments: the switch only removes the
  size test, and the voting step gives a seed to a minority of the fragments (I assume because the votes
  of a small object overshoot its centre; not examined). So this variant understates what counting every
  fragment costs. The level-3 LoG route does count them: median +9.2% / +9.7%, score 0.05.
- **Free-plateau fit** (0.350 / 0.350): the rising marker series and the partial marker inhibitor get a
  number instead of `inactive`, the cytostatic compound gets a GI50, and the largest GI50 error among the
  active compounds is 0.44 / 0.49.
- **Counting frame ignored**: 1.000 / 1.000. Not penalised (see limitations).
- **Level-3 decision points on the level-4 base** (no focus QC, whole-field readout, no plate model, no
  background, plate-wide background): 0.35-0.719 / 0.513-0.55, all below 0.75.
  `naive_global_background` is at 0.719 on seed 1 (one compound wrong, as at level 3).
- Potency errors of the two solutions: <= 0.065 (seed 1) and <= 0.087 (seed 2) log10 units.

Gate facts (final runs, family code unchanged since 02:48 UTC): seed 1 all 31 checks ok incl. Harbor
oracle 1.0 / nop 0.0, 11.07 min; seed 2 all 29 checks ok, 9.92 min. Oracle runtime 0.14 / 0.15 min.
Second solution (three Cellpose settings on the annotated fields, then 96 fields at 12 px, 2 threads)
5.6 / 5.8 min, `l3_reference` 1.9 / 1.9 min (differences of the file times under `results/work/`,
load average 1-4). The same gates had been run once before on the same inputs, solutions and verifier
(only the judge text of the task card and `oracle_runtime_min` changed in between, worklog L4-8):
every score was identical; the seed-1 second solution took about 16 min there under a load average of
7-11 (estimated from the container start time and the gate log). Agent image 2.82 GB (`docker images`),
same environment layers as levels 1-3; input 20 MB.

### Does the criterion on the annotated fields predict the hidden wells?
Outline-voting counter with 16 parameter sets per seed (the 4 the oracle chooses from and 12
deliberately detuned ones: larger or smaller vote smoothing, peak distance, radii, edge threshold,
no size adaptation), counted on the final inputs of seeds 1 and 2:

| annotated fields inside the criterion | settings (seed 1 + seed 2) | other in-focus wells inside the tolerance |
|---|---|---|
| 8 of 8 | 7 + 11 | 1.000 in 16 settings, 0.976 and 0.867 in the other two |
| 7 of 8 | 8 + 2 | 0.747-1.000 (below 0.85 in one setting) |
| 5 or 6 of 8 | 1 + 3 | 0.566, 0.590, 0.627, 0.952 |

Other methods (tables above): tuned Cellpose 8/8 -> 0.964 / 1.000; level-3 methods 2-4 of 8 ->
0.434-0.554; fragments counted 5/8 -> 0.879. So in every measured case a procedure that meets the
criterion on all eight fields gets full count credit (>= 0.85), the closest case being 0.867. The
reverse does not hold: a procedure that misses one field by a little can still be inside on most
hidden wells.

### Development measurements that shaped the design (seed 1, earlier image parameters)
These plates differ from the released ones (spacing, nucleoli, cytoplasm), so the numbers are not
comparable with the tables above; they are the reason for the released parameters.

| image parameters | method | result |
|---|---|---|
| spacing 9.0 px, deep nucleoli | Cellpose `nuclei`, annotated fields, mean / worst field: 17 px flow 0.4 (level 3) | -20.8% / 35% |
| | 17 px, flow check off | -8.7% / 22% |
| | 12 px, flow check off | -5.3% / 13% |
| | `cyto3` 12 px | -3.0% / 8.0% |
| spacing 9.5 px, light nucleoli (released) | 17 px flow 0.4 | -17.1% / 26% |
| | 12 px, flow check off | -1.7% / 6.1% |
| | `cyto3` 14 px | -1.7% / 5.6% |
| spacing 10.0 px, light nucleoli | 17 px flow 0.4 | -6.4% / 13% |
| | 12 px, flow check off | -0.8% / 2.7% |
| spacing 9.5 px | classical, hidden wells inside 10% / median error: LoG seeds on the smoothed image, per-well scale | 0.24-0.28 / -13..-15% |
| | distance-transform maxima, several thresholds | 0.12-0.19 / -16..-33% |
| | LoG on the raw image, sigma 2.4 px, fragment test, counting frame | 0.92 / -1.6%, one annotated field outside by one cell |
| | outline voting (the oracle's method), five parameter sets | 0.96-1.00 / -0.8..-3.5% |

Spacing 9.5 px is the window in which the level-3 Cellpose settings fail clearly and a tuned setting
passes with some margin; at 10.0 px the level-3 reference was only just outside 10%, at 9.0 px no
setting of the `nuclei` model passed.

Seeds 3 and 4 (generated into scratch, not gated, classical pipelines only): oracle 1.0 / 1.0 (hidden
wells 1.000 / 1.000, worst annotated field at 0.54 / 0.67 of the limit, potency errors <= 0.076),
level-3 watershed 0.219 / 0.340, level-3 LoG 0.05 / 0.219, fragments counted 0.481 / 0.478,
free-plateau fit 0.35 / 0.35. The Cellpose routes were not run on these seeds.

## Deviations from the level-4 spec, and what was dropped

| spec | what I did | why |
|---|---|---|
| "potencies keep their weight" (0.45 / 0.25), counts on their own component | counts 0.25, validation 0.05, QC 0.05, marker 0.45, count potency 0.20 | measured: Cellpose with the level-3 settings loses nuclei almost uniformly in dense fields, the response ratios survive, all 16 potency entries stay right or nearly right. With counts at 0.20 that submission scored 0.736 and would have scored 0.75 with one more potency inside the tolerance - no margin under the 0.75 ceiling. Now a submission whose counts fail cannot exceed 0.70. |
| reference solution = tuned Cellpose and/or careful classical | the classical outline-voting pipeline is the **oracle**, tuned Cellpose is the second solution | build-host load (4 cores, load average 12-14): the Cellpose route needs minutes of CPU per plate and would have run twice in every gate run (G2 and Harbor). Both routes are validated. Consequence: `oracle_steps = 1`, `oracle_runtime_min = 0.2`, a-priori difficulty 0.564 (level 3: 0.665). AUTHORING asks for >= 20 steps at level 4; the honest number is 1 (the Cellpose route has 8 invocations). |
| "nuclear size that changes with treatment" | enlargement only (two compounds) | a shrinking phenotype would put small nuclei next to large fragments and make the SOP's fragment definition ambiguous |
| "use real nuclear texture if you can do it with exact truth" | synthetic texture and nucleoli | not attempted: it needs an external image set with instance masks, a crop library and a second look at the truth of composited crops; no time left after the difficulty window had been found. Nothing was downloaded, `.cache/hcs_dose_response/` does not exist. |
| optional `I_time_budget` | dropped | not attempted. The tuned Cellpose route already takes several minutes for 96 fields of 320 x 320 px; larger fields would exceed the 25 MB input budget (the plate is 20 MB) |
| "one variant per new element" | Q: `naive_fragments_counted`, the three level-3 routes, the plate-wide factor (two bases), manual counts copied. P: `naive_free_plateau_fit` | the counting frame alone is **not** an attack: ignoring it changes the counts by about 2% and the submission scores 1.0 on both seeds (measured, table above) |
| `l3_reference` for both level-3 routes | `l3_reference` (Cellpose, the level-3 oracle), `l3_reference_classical` (watershed), `l3_independent_solution` (LoG, the level-3 second solution) | all three exist at level 3; each gets a truthful `validation.tsv` from `assets/l3_validation.py`, otherwise it would fail on a missing file and test nothing |
| "untuned Cellpose" as a careless variant on the level-4 pipeline | measured by hand on seed 1 only, not in the gated battery | one Cellpose plate run per gate saved; `l3_reference` uses the same segmentation |

## Known limitations and things I am unsure about (level 4)
- **No agent pilot.** Whether the task is hard for a model is untested. What the analyst is told: the
  criterion, the counting rules, that nuclei overlap. What is not: that the default settings fail, in
  which direction, and what repairs it.
- **Two ways through, both narrow.** (a) Cellpose: the level-3 settings fail; switching the flow-error
  check off and setting the diameter to the nuclei (12 px) passes; the fragments then appear as objects
  and have to be removed by size. During development I tried ten settings of the `nuclei` model and
  three of `cyto3`, all on earlier image parameters and on the annotated fields only (worklog L4-2,
  L4-3). On the released images I ran 12 px / flow off on the annotated fields of both seeds and
  17 px / 0.4 on seed 1; what the other settings do there I do not know. (b) Classical: threshold +
  distance transform and LoG seeds failed in every variant I tried; outline voting passes. An analyst
  who knows neither way has to find one, with the annotated fields as the only guide.
- **The second solution has little margin on the annotated fields**: worst field at 0.71 / 0.67 of the
  limit (seeds 1 / 2). On another seed it may fail a field; then the oracle still passes but the
  `independent_solution` probe would not run to the end (`choose` stops).
- **Meeting the criterion is sufficient in every case I measured, but not strictly necessary.** See
  the parameter sweep: 18 of 32 settings meet it on all eight fields and all of those have >= 0.867 of the
  other wells inside (full credit starts at 0.85); settings that miss one field by a little can still
  have 0.90-1.00 of the other wells inside and then lose only part of the 0.05 validation weight.
- **The counting frame is stated and annotated but not separately enforced** (2% effect, inside the
  tolerance).
- **`naive_fragments_counted` is a weak version of its own idea**: it disables the size test of the
  oracle, but the oracle's voting step finds only part of the fragments in the first place (hidden wells
  0.879 inside). It fails through the annotated fields and the potencies. A pipeline that counts every
  fragment is represented by `l3_independent_solution` (0.05 on both seeds).
- **`l3_reference` fails on the count and validation components only** (see its row): the potencies
  survive an almost uniform undercount. Its distance to the 0.75 ceiling comes from the weights.
- **`naive_global_background`** keeps the thin margin it had at level 3 (one compound wrong).
- **Synthetic images.** Elliptical nuclei of similar size with an unbroken outline are what makes
  outline voting work so well; I expect it to do much worse on real plates (not tested). Fragments are
  round, bright, never overlap a nucleus and never reach the border. Nucleoli and texture are much simpler
  than real chromatin.
- **Wording in the task card.** `README.md` of the two tasks says fragments "lie at least 3.5 px from any
  nuclear outline, 7 px from the image edge and 8 px from each other". These are distances of the fragment
  centre, as stated above. Not corrected, because the text is generated from `level4.py` and the gate
  records belong to the code as it is. The agent-facing SOP makes no statement about these distances.
- **Fragment definition.** The SOP gives sizes in micrometres of the object. My first classical attempts
  tested the area of the thresholded mask (>= 50 px) and over-counted fragment-rich wells (worst well +18%
  to +63%, development plates); I attribute that to the blurred image of a bright fragment being larger
  at a low threshold, without having checked it object by object. The half-maximum area and the Cellpose
  mask area both separate fragments from nuclei on both seeds; other size measures were not tried.
- **The rising marker series** is planted as a real effect and the SOP now lists "a readout that does not
  fall at all" under `inactive`. Whether that sentence makes the element too easy I cannot say.
- **Cellpose runtime.** 8 fields at 12 px took 21-23 s on the quiet host and 2.7-3.5 min under load;
  the second solution took 5.6 / 5.8 min in the final gate runs. On a loaded machine an agent that tunes
  Cellpose on the whole plate instead of the annotated fields will spend most of its time waiting.
- **Level 4 was validated on seeds 1 and 2 only.**

# Level 7 - the plate lay on the stage the other way round (added 2026-11-12, annotated fields revised 2026-11-14)

Spec: `docs/specs/l7/hcs_dose_response.md`. Code: `l7.py` (generator, battery, task card),
`assets/hcs_l7_voting.py` (reference solution and the by-name variants), `assets/hcs_l7_cellpose.py`
(independent second solution), `assets/solve_l7.sh`, `assets/score_l7.py` (wrapper around the level-4
verifier, which is shipped unchanged as `tests/score_l4.py`). `level4.generate_l4` got one optional
argument (`plate_out`: a dict that receives the planted plate; `None` at levels 4-6) and `__init__.py`
sends `level == 7` to `l7.py`. The task is built by `python -m alegen make hcs_dose_response --level 7`:
the family returns the level-7 build with the level-4 texts and `descaffold.apply(build, 7)` removes the
four remarks of the assay notes exactly as for level 6. `descaffold.py` was not changed and every edit
applied.

"Gate" = inside the task images during `python -m alegen make`. "Dev run" = the same script in a container
of the level-4 agent image (`--network none`, input directory of the generated task mounted) and scored
with the task's own `tests/score.py` on the host; where both exist they agree.

**Revision of 2026-11-14 (maintainer decision).** In the first version (2026-11-12) the annotation files
kept their level-4 names and held the cells of whatever image lay under those names. That changed which
phenotypes the manual counts cover, i.e. it made the counting part easier to pass on the annotated fields
and harder to get right elsewhere. The maintainer reported that the yardstick model, in two attempts on
each of seeds 1 and 2, never questioned the orientation (both potency components 0 every time) but met the
SOP criterion on those annotated fields with a quick Cellpose setting and then failed on the other wells
(`cell_counts` 0.0-0.73); I have not read the trial records myself. Now the lab annotated the same physical
fields as at level 4, so that level 7 differs from level 6 only by the re-filing. Everything in this
section that depends on the annotated fields was measured again on 2026-11-14 and both seeds were gated
again; numbers that do not depend on them are those of 2026-11-12 (the planted control statistics; the
potency tables of the by-name analyses, which I checked to be the same files on seeds 1, 2, 3 and 10).

## What was planted

One decision point, `T_plate_orientation`. The image pair named for well (row r, column c) shows plate
well (9 - r, 13 - c): `A01` shows the well the layout calls `H12`. Nothing else differs from the level-4
task of the seed:

- `plate_layout.tsv`, `analysis_sop.md`, `assay_notes.md`, the instruction, the flat-field and dark frames:
  the level-4 files (the notes then lose the four level-6 remarks).
- Images: the level-4 pixel data, filed under the opposite well name (`images/A01_ch1.tif` of level 7 is
  byte-identical to `images/H12_ch1.tif` of level 4; the TIFF tags hold no name). A field is turned with
  its plate. I did not turn the pixel data: the simulated cell layer has no direction (nuclei are placed
  uniformly on a range that is symmetric about the image centre, the rules at the two low and the two high
  edges mirror each other, fragments keep the same margin on all sides) and the illumination falls off
  radially, so the level-4 rendering of a well is a draw from the same distribution as its turned image.
  What is not symmetric is the SOP's counting frame, and that is applied to the image as it is filed.
- Annotations: the lab counted the same eight physical fields as at level 4, the plate wells that
  `sim4.choose_annotated` picks for coverage (two dense and one sparse untreated field, the top and a middle
  concentration of the cytotoxic compound, the top concentrations of the two compounds that enlarge the
  nuclei, a high concentration of the dual compound). The acquisition software names by stage position, so
  each annotation file carries the name of the image that shows its field: level-7 `annotations/H01.tsv`
  is level-4/6 `annotations/A12.tsv`, byte for byte (same clicks, same manual count), and its coordinates
  are those of level-7 `images/H01_ch*.tif`, which is level-4 `images/A12_ch*.tif`. Level 4 never
  annotates an out-of-focus field (`choose_annotated` draws from the in-focus wells and `generate_l4`
  asserts it); `l7.py` raises if it ever did. The replacement rule of the first version is gone.
- Truth: `truth.wells[<image name>]` holds what that image shows (cell number, focus state, annotated or
  not, plus `shows_plate_well` and the layout entry of the name); `truth.compounds` is untouched.

## Truth

Unchanged in meaning. Potencies are the planted 0.5 crossings of the compounds, as at level 4. The three
per-image quantities (cell number, focus state, manual count) are properties of an image and are stored
under the image name. The orientation itself is not a deliverable and is not scored.

## What goes wrong without the step

Measured with the level-4 reference solution, unchanged, on the level-7 inputs (`l6_reference`), seeds
1 / 2 (gate and dev run), 3 and 10 (dev run). Measured again on 2026-11-14 with the revised annotation
files: same scores, and `dose_response.tsv` of all four seeds is identical to that of 2026-11-12:

| seed | score | counts within tolerance | validation | focus flags | marker entries right | count entries right | reported |
|---|---|---|---|---|---|---|---|
| 1 | 0.35 | 100% | 8 of 8 | exact | 4 of 8 | 5 of 8 | 16 x `inactive` |
| 2 | 0.35 | 100% | 8 of 8 | exact | 3 of 8 | 4 of 8 | 12 x `inactive`; marker 0.0013, 0.0027 and 17.1 uM; count 0.0019 uM |
| 3 | 0.35 | 100% | 8 of 8 | exact | 4 of 8 | 5 of 8 | 16 x `inactive` |
| 10 | 0.35 | 100% | 8 of 8 | exact | 3 of 8 | 5 of 8 | 14 x `inactive`; marker 20.0 and 21.5 uM |

- Everything that is per image is right, which is what makes the result look sound: the validation
  table meets the SOP criterion on all eight fields.
- Both potency components are zero. A component needs a mean credit above 0.8; `inactive` everywhere is
  right for 4 of 8 marker and 5 of 8 count entries on every plate of this family (the generator asserts
  four marker-active and three count-active compounds), which is what seeds 1 and 3 deliver.
- How wrong the numbers are (seed 2): 0.0013 uM reported for a planted 5.34 uM (3.6 log10 units),
  0.0027 uM for 3.53 uM (3.1), 17.1 uM for a compound that is inactive on the marker, 0.0019 uM for a
  compound that is inactive on cell number. Three of the four are below the lowest tested concentration
  (0.0137 uM). The planted actives that come out `inactive`: seed 1 all seven active entries, seed 2 five
  of seven.
- Would the curves look acceptable? The table does; the curves do not, to anyone who plots them. A
  series holds the wells of the opposite row in reverse order with that row's vehicle and reference wells
  in between. Seed 2 (few vehicle wells are off, so the plate model is nearly right): smooth curves that
  rise with the nominal concentration (0.04, 0.08, 0.23, 0.54, 0.81, 0.95, 1.02, 0.90; reported
  `inactive`) and single wells at 0.05 inside otherwise flat series. Seed 1 (13 of 23 vehicle wells off,
  plate model distorted): jagged series with responses up to 2.0 (0.09, 1.47, 0.75, 0.08, 1.59, 1.53,
  1.44, 0.94; reported `inactive`; values of 2026-11-14, on 2026-11-12 three of them were 0.01 different
  and the largest response was 1.9). The SOP model has its upper asymptote fixed at 1 and falls with the
  concentration, so a rising series is fitted by a flat line (`inactive`) or by a step at the lower end
  of the search range. A pipeline that only writes the table never shows this.

## How it is discovered

By looking at the control wells before they are used. A plate fits into a stage holder in two
orientations, so there is exactly one alternative to try. Numbers for both orientations, seeds 1 / 2.

On the planted values (generator, `plants.plate_orientation.evidence_on_planted_values`; out-of-focus
images left out):

| statistic | by name, s1 | turned, s1 | by name, s2 | turned, s2 |
|---|---|---|---|---|
| marker per cell of the wells called reference / mean of the two vehicle wells of their column | 0.29-2.23 (median 1.74) | 0.05-0.05 (median 0.05) | 0.05-1.16 (median 1.04) | 0.05-0.06 (median 0.05) |
| wells called reference below half | 1 of 7 | 8 of 8 | 1 of 7 | 8 of 8 |
| the two wells called vehicle of a column, larger / smaller marker per cell: median / maximum | 1.60 / 22.01; 6 of 11 beyond 1.3 | 1.04 / 1.14; 0 of 11 beyond 1.3 | 1.04 / 1.44; 3 of 12 beyond 1.3 | 1.05 / 1.17; 0 of 11 beyond 1.3 |
| cell number of the wells called vehicle / median of their plate area (inner, outer ring) | 0.17-1.13; 6 of 23 beyond 1.3 | 0.92-1.06; 0 of 23 beyond 1.3 | 0.13-1.08; 3 of 24 beyond 1.3 | 0.95-1.05; 0 of 23 beyond 1.3 |
| wells called vehicle with a planted response beyond a factor 1.3 on either readout | 11 of 23 | 0 of 23 | 4 of 24 | 0 of 23 |
| dose series, cell number: falling / rising (rank correlation beyond 0.7) | 1 / 1 | 4 / 0 | 0 / 1 | 4 / 0 |
| dose series, marker per cell: falling / rising | 0 / 0 | 5 / 1 | 1 / 3 | 5 / 0 |
| balance of falling over rising series (16 series-readouts) | 0 | 8 | -3 | 9 |

Measured on the images by the reference solution (`control_check.tsv` in its output; plate model on
the wells called vehicle, responses relative to it):

| statistic | by name, s1 | turned, s1 | by name, s2 | turned, s2 |
|---|---|---|---|---|
| wells called reference, median marker response | 1.32 | 0.044 | 1.00 | 0.046 |
| wells called reference below half | 1 of 7 | 8 of 8 | 1 of 7 | 8 of 8 |
| wells called vehicle more than a factor 1.3 off the plate model | 13 of 23 | 0 of 23 | 4 of 24 | 0 of 23 |
| marker response of the wells called vehicle, range | 0.05-1.46 | 0.91-1.08 | 0.72-1.19 | 0.91-1.09 |
| count response of the wells called vehicle, range | 0.32-1.75 | 0.93-1.07 | 0.15-1.04 | 0.94-1.06 |

(2026-11-14, revised annotation files. On seed 1 the voting parameters are now those of level 4, which
moves three of the ranges in the second decimal against 2026-11-12; seed 2 is unchanged.)

Seeds 3 / 10 (dev run): reference median 2.06 / 0.89 by name and 0.046 / 0.046 turned; vehicle wells off
15 of 22 / 11 of 22 by name and 0 of 23 turned.

The second solution does not use the controls. It normalises the marker per cell by a robust line
through log(marker) against the image column over all in-focus fields and the cell number by the median
of the outer ring / of the inner fields (both sets are the same in either orientation), and counts dose
series whose readout falls or rises with the concentration (Spearman rho beyond 0.7, 16 series-readouts).
Dev runs of 2026-11-14 on the gated tasks, falling / rising: seed 1 by name 0 / 2 (balance -2), turned
9 / 1 (balance 8); seed 2 by name 1 / 5 (balance -4), turned 10 / 1 (balance 9). Seed 1 gave the same
counts on 2026-11-12.

A second sign in the inputs, not a pointer in any text: the manual counts. Since the revision the
annotation files carry names under which the layout lists something else than what was counted. The
assay notes say how the fields were picked (untreated wells of the inner and outer plate area, wells of
the highest and of intermediate concentrations); the files of level 7 (`truth.json`:
`plants.plate_orientation.annotated_fields_as_the_layout_names_them`; "untreated" = median planted cell
number of the in-focus vehicle wells of the same plate area):

Seed 1:

| file | manual count | the layout calls this name | the field counted (plate well) | count / untreated |
|---|---|---|---|---|
| `A01.tsv` | 38 | CPD-8550 at 0.0137 uM | CPD-7827 at 30 uM (H12) | 0.58 of 66 |
| `B06.tsv` | 127 | reference | vehicle (G07) | 1.06 of 120 |
| `D08.tsv` | 113 | vehicle | vehicle (E05) | 0.94 of 120 |
| `E03.tsv` | 31 | reference | CPD-4184 at 10 uM (D10) | 0.26 of 120 |
| `F01.tsv` | 18 | CPD-9526 at 0.0137 uM | CPD-9251 at 30 uM (C12) | 0.27 of 66 |
| `H01.tsv` | 10 | vehicle | CPD-8550 at 30 uM (A12) | 0.15 of 66 |
| `H02.tsv` | 66 | CPD-7827 at 0.0137 uM | vehicle (A11) | 1.00 of 66 |
| `H07.tsv` | 34 | CPD-7827 at 0.37 uM | CPD-8550 at 0.37 uM (A06) | 0.52 of 66 |

3 of 8 files below half of an untreated field, 2 of them under a name the layout calls vehicle or reference; 7 of 8 names stand for something else in the layout than the field that was counted.

Seed 2:

| file | manual count | the layout calls this name | the field counted (plate well) | count / untreated |
|---|---|---|---|---|
| `A01.tsv` | 40 | reference | CPD-5484 at 30 uM (H12) | 0.61 of 66 |
| `B03.tsv` | 33 | CPD-2467 at 0.0412 uM | CPD-1987 at 10 uM (G10) | 0.28 of 118 |
| `B10.tsv` | 123 | vehicle | vehicle (G03) | 1.04 of 118 |
| `E01.tsv` | 21 | CPD-1142 at 0.0137 uM | CPD-1550 at 30 uM (D12) | 0.32 of 66 |
| `E12.tsv` | 68 | vehicle | vehicle (D01) | 1.03 of 66 |
| `F03.tsv` | 116 | CPD-7090 at 0.123 uM | vehicle (C10) | 0.98 of 118 |
| `H01.tsv` | 8 | vehicle | CPD-3353 at 30 uM (A12) | 0.12 of 66 |
| `H05.tsv` | 33 | vehicle | CPD-3353 at 1.11 uM (A08) | 0.50 of 66 |

3 of 8 files below half of an untreated field, 1 of them under a name the layout calls vehicle or reference; 6 of 8 names stand for something else in the layout than the field that was counted.

An analyst who puts the eight manual counts next to the layout sees untreated-looking counts under
concentrations and very low counts under vehicle or reference names. I count this as evidence in the data.
It is a second route to the decision point besides the control wells. In the first version the names
were those of level 4 and the files held what lay under them (seed 1: `A12`, which the layout calls 30 uM
of the cytotoxic compound, held 68 cells).

The two mirror maps of the grid (rows only, columns only) are not orientations a plate can take. On the
planted values of seeds 1 / 2 they leave 6 of 6 / 6 of 7 (rows) and 5 of 7 / 5 of 7 (columns) wells called
reference uninhibited and 14 / 11 and 9 / 14 of 23 vehicle wells off.

## Seeds that cannot have a level 7

The layout puts the vehicle wells of a row into every fourth column with a per-row offset and the
reference well two columns further. If the offsets of opposite rows are complementary the vehicle wells
fall on vehicle wells under the turn, and the reference wells can do the same; such a plate shows
well-behaved controls in both orientations and the turn would be visible in the dose series only. The
generator therefore checks the planted values and raises for a seed unless, with the layout by name,
fewer than half of the in-focus wells called reference are below half of the vehicle wells of their column
and their median is at least 0.5, at least four wells called vehicle are off by more than a factor 1.3 on
either readout, and the balance of falling over rising dose series is at least 5 turned and at least 5
lower by name (and, turned, all reference wells are below 0.2 and no vehicle well is off). Survey of seeds
1-400 on the planned values (no rendering, `l7.survey`): 41 refused (31 on the reference criterion, 12 of
them also on the vehicle criterion; 10 on the vehicle criterion alone; none on the series criterion alone).
Among seeds 1-40: 4, 8, 12, 15, 26, 32, 33 and 36 are refused. Seed 2 passes the vehicle criterion with
exactly four wells.

## The two solutions

- Reference (`hcs_l7_voting.py`): the level-4 voting pipeline per image. The typical nuclear signal is
  taken over all fields instead of over the vehicle wells, so that nothing in the image analysis depends on
  the layout. Then `control_check` for both orientations (table above); the orientation in which the
  median reference response is below 0.5 and at most 15% of the vehicle wells are more than a factor 1.3
  off the plate model is used; if neither or both qualify the script stops. The plate model is fitted in
  plate coordinates. Per-well rows carry the image name; two extra columns (`image`, `plate_well`) say
  which plate well each image shows.
- Second (`hcs_l7_cellpose.py`): the level-4 Cellpose pipeline (intensity scale from all fields instead
  of the vehicle wells), a new stage `orient` (series direction, see above; stops unless the balance is
  at least 4 and at least 4 above the other orientation), then the level-4 statistics.
- No code is shared. Gate of 2026-11-14, seeds 1 / 2: reference 1.000 / 1.000, second solution 1.000 / 1.000.
- Tuning on the annotated fields (dev runs of 2026-11-14 on the gated tasks; seeds 1 and 2 also passed
  in the gate). Reference solution: outline threshold 0.12, vote sigma 1.2 on both seeds, worst annotated
  field at 0.76 / 0.33 of the SOP limit; every well without annotation inside the count tolerance. On
  seed 1 that is the level-4 result to the digit. Second solution: `nuclei`, diameter 12, flow check off
  on both seeds, worst annotated field at 0.63 / 0.67 of the limit (the two other settings it tries fail
  with 1.97-2.50); 96.4% / 100% of the wells without annotation inside the tolerance; largest potency
  error 0.053 / 0.085 log10. The level-4 notes give 0.71 / 0.67 for the worst field of the second solution;
  the seed-1 difference goes with the intensity scale, which is now taken over all fields.
- One variant shows the same code difference in the reference pipeline: `naive_fragments_counted` scores
  0.481 on seed 1 where the level-4 gate has 0.3125. Measured on the level-4 inputs of seed 1: level-4
  code 0.3125, level-7 code with the layout by name 0.481 (worst annotated field at 12.9 against 12.3 of
  the limit; both choose threshold 0.15, sigma 1.5). The full solution and the variant without focus QC
  give identical scores and potency errors with both codes on those inputs.

## Scorer

`score_l7.py` imports the level-4 verifier (`score_l4.py`, byte-identical to `tests/score.py` of levels 4
and 6) and runs its `main` under four namings of the per-well rows: each of `well_measurements.tsv` and
`validation.tsv` by image name or by the plate well the image shows. The best total is the score.
Tolerances, weights and gates are untouched.

Why two namings. The contract says "one row for every well of `plate_layout.tsv`" and "`field` is the
well id". Until the turn is known, the well of an image, of an annotation file and of the layout are the
same thing. Afterwards an analyst has to choose: keep each measurement under the name of the file it was
made on, or file it under the plate well it belongs to, which is what the project team needs when it joins
the table with the layout. I could not find a sentence that fixes it. The validation table ties `field`
to the annotation file and to `n_cells` of "that well", which speaks for the image name, but the same
sentence is satisfied when both tables are re-filed together. The spec asks to accept both in that case.
The mixed combination (well table by plate well, validation table by annotation file name) is accepted
as well, because the validation row is then still compared with the count reported for the same image.

What this cannot do: lift a by-name analysis. Potencies carry no well names and weigh 0.65; measured, the
by-name submissions score 0.35 under the naming they use and 0.0 under the three others (seeds 1 and 2).

## Battery

Pass threshold 0.9. Gate runs of 2026-11-14 with the revised annotation files, both seeds with Harbor
(`python -m alegen make hcs_dose_response --seed N --level 7`). The records of the first version are
kept in `<scratch>/hcs_l7_check/prev_l7/`.

| check | bound | seed 1 | seed 2 |
|---|---|---|---|
| G1 images_build | | ok | ok |
| G7 tools_in_login_shell | | ok | ok |
| oracle | >= 0.99 | 1.000 | 1.000 |
| G5 verifier_deterministic | | ok | ok |
| planted_truth | >= 0.999 | 1.000 | 1.000 |
| within_tolerance | >= 0.99 | 1.000 | 1.000 |
| independent_solution | >= 0.99 | 1.000 | 1.000 |
| per_well_tables_by_plate_well | >= 0.999 | 1.000 | 1.000 |
| well_table_by_plate_well_validation_by_image | >= 0.999 | 1.000 | 1.000 |
| G4 has_l3_reference_attack | | ok | ok |
| G4 has_l6_reference_attack | | ok | ok |
| empty_output | <= 0.0 | 0.000 | 0.000 |
| inputs_copied | <= 0.1 | 0.000 | 0.000 |
| zero_byte_files | <= 0.0 | 0.000 | 0.000 |
| header_only | <= 0.1 | 0.000 | 0.000 |
| l6_reference | <= 0.75 | 0.350 | 0.350 |
| layout_by_name_everything_else_right | <= 0.75 | 0.350 | 0.350 |
| controls_swapped_only | <= 0.75 | 0.350 | 0.350 |
| l3_reference | <= 0.75 | 0.050 | 0.050 |
| l3_reference_classical | <= 0.75 | 0.050 | 0.050 |
| l3_independent_solution | <= 0.75 | 0.050 | 0.050 |
| counts_scaled_to_annotation | <= 0.75 | 0.050 | 0.050 |
| counts_scaled_to_annotation_log | <= 0.75 | 0.050 | 0.050 |
| manual_counts_copied | <= 0.75 | 0.100 | 0.100 |
| naive_fragments_counted | <= 0.75 | 0.481 | 0.481 |
| naive_free_plateau_fit | <= 0.75 | 0.350 | 0.350 |
| naive_no_focus_qc | <= 0.75 | 0.544 | 0.544 |
| naive_whole_well_readout | <= 0.75 | 0.550 | 0.550 |
| naive_no_plate_model | <= 0.75 | 0.350 | 0.513 |
| naive_no_background | <= 0.75 | 0.550 | 0.550 |
| naive_global_background | <= 0.75 | 0.719 | 0.550 |
| prior_only | <= 0.1 | 0.000 | 0.000 |
| all_potencies_beyond_tolerance | <= 0.4 | 0.350 | 0.350 |
| one_potency_beyond_tolerance | <= 0.75 | 0.719 | 0.719 |
| counts_beyond_tolerance | <= 0.75 | 0.700 | 0.700 |
| G6 harbor_oracle | | ok [1.0] | ok [1.0] |
| G6 harbor_nop | | ok [0.0] | ok [0.0] |
| all gates | | ok (37 checks, 12.36 min) | ok (37 checks, 12.99 min) |

- `l6_reference`, `layout_by_name_everything_else_right` and `controls_swapped_only` are the three
  additions the spec asks for. `layout_by_name_everything_else_right` hands the planted per-image values
  to the reference plate model and fit (`--measurements`), layout by name; its potency table is identical
  to that of `l6_reference` on both seeds, i.e. the measurement error of the image analysis plays no role
  in the failure. `controls_swapped_only` exchanges the two control roles and nothing else (seed 1: 14
  numeric potencies, 1 of 8 marker entries right; seed 2: 3 of 8 marker and 5 of 8 count entries right).
- The seven "ignores exactly one element" attacks of level 4 run on `hcs_l7_voting.py` (controls checked,
  layout turned, that element ignored). The six attacks built on level-3 routes run as they were written,
  with the layout by name; they are wrong twice.
- Truth-derived probes and attacks are re-filed: every per-well row moves with its image.
- Two new probes file the planted truth by plate well (both tables; well table only).

## `naive_no_focus_qc`: ceiling 0.75 again, and what the first version showed

With the revised annotation files the reference solution tunes on the level-4 fields and selects the
level-4 parameters. Measured on 2026-11-14 (dev runs; seeds 1 and 2 also in the gate):

| inputs (revised annotation files) | outline threshold chosen by the tuning | worst annotated field / SOP limit | `--qc none` | full solution | `l6_reference` |
|---|---|---|---|---|---|
| level 7, seed 1 | 0.12 | 0.76 | 0.544 | 1.000 | 0.350 |
| level 7, seed 2 | 0.12 | 0.33 | 0.544 | 1.000 | 0.350 |
| level 7, seed 3 | 0.12 | 0.54 | 0.669 | 1.000 | 0.350 |
| level 7, seed 10 | 0.12 | 0.79 | 0.581 | 1.000 | 0.350 |

All four are under 0.75, so the attack has its level-4 ceiling again and the task card no longer carries
an exception. Seed 1 reproduces the level-4 numbers exactly (0.12, 0.76, 0.544).

History, because it says something about levels 4 and 6 too. In the first version the first gate run of
seed 1 failed on exactly this check (0.825 against 0.75; record in
`<scratch>/hcs_l7_check/first_gate_run_s0001_l7_ceiling_0.75.json`) and I had raised the ceiling to 0.85.
What I measured then (dev runs, 2026-11-12, first-version annotation files unless the row says level 4):

| inputs | outline threshold of the voting step | how it was set | worst annotated field / SOP limit | `--qc none` | full solution |
|---|---|---|---|---|---|
| level 7 first version, seed 1 | 0.15 | chosen by the tuning | 0.25 | 0.825 | 1.000 |
| level 7 first version, seed 1 | 0.12 | forced | 0.31 | 0.544 | 1.000 |
| level 7 first version, seed 2 | 0.12 | chosen by the tuning | 0.29 | 0.544 | 1.000 |
| level 7 first version, seed 10 | 0.15 | chosen by the tuning | 0.50 | 0.825 | 1.000 |
| level 4, seed 1 | 0.12 | chosen by the tuning (first of a tie) | 0.76 | 0.544 | 1.000 |
| level 4, seed 1 | 0.15 | forced | 0.76 | 0.825 | 1.000 |

With the threshold at 0.15 the voting step finds no nucleus in any of the five out-of-focus fields
(`n_cells` 0, marker `NA`); with 0.12 it finds 19-43 objects per field with a marker value that is too
low. In the first case the wells drop out of the marker fit by themselves and enter the count fit with
zero cells: the marker IC50s stay right (largest error 0.03 log10), one GI50 is 1.4 log10 units low, and
the submission loses the QC flags (0.05) and 0.125 of the count component. In the second case both
potencies of the affected compounds are pulled down, as the level-3 design intended.

What remains true after the revision: on the level-4 fields of seed 1 both thresholds are equally good
(0.76 of the limit), the reference tuning takes the first, and a counting procedure that finds nothing
in an out-of-focus field, combined with no focus QC, scores 0.825. Measured with the threshold forced to
0.15 on the level-4 inputs (2026-11-12) and on the revised level-7 task (2026-11-14) of seed 1; the full
solution with that threshold scores 1.000 on both. Level 6 has the images and annotation files of level
4. That is under the pass line by 0.075 and it is a property of the level-4 weights (QC flags 0.05, count
potency 0.20), not of level 7. The battery does not contain that submission any more; I list it under
the limitations.

## Levels 1-6 after the change

`__init__.py` checks `level == 7` before the existing `level >= 4` branch; `level4.generate_l4` has one
more optional argument. Levels 1-6 do not import `l7.py` (checked: importing the family loads neither
`l7` nor `level4`). The revision of 2026-11-14 touched `l7.py` only. With the final code, regenerated
into a scratch directory with `--no-validate` and compared with `diff -r`, nothing ignored (2026-11-14):

| comparison | files | result |
|---|---|---|
| `hcs-dose-response-s0001-l4` regenerated against `tasks/` | 219 | no difference |
| `hcs-dose-response-s0001-l6` regenerated against `tasks/` | 219 | no difference |
| `hcs-dose-response-s0001-l7` regenerated after its gate against `tasks/` (the gated task is what the final code writes) | 220 | no difference |
| `hcs-dose-response-s0002-l7` regenerated after its gate against `tasks/` | 220 | no difference |

The revised level-7 tasks in `tasks/` against the first version (sha256 manifests of the first version taken
before anything was changed, `<scratch>/hcs_l7_check/prev_l7/`) and against level 6 of the same seed
(`tasks/hcs-dose-response-s0001-l6`; `tasks/hcs-dose-response-s0002-l6` as the maintainer's gate of
2026-11-14 wrote it), seeds 1 / 2:

| what | seed 1 | seed 2 |
|---|---|---|
| files identical to the first version | 210 of 220 | 210 of 220 |
| images identical to the first version | 192 of 192 | 192 of 192 |
| files with other content | `README.md`, `tests/reference/truth.json` | `README.md`, `tests/reference/truth.json` |
| annotation files, first version -> revised | A06 A11 A12 C12 D10 E05 G07 H12 -> A01 B06 D08 E03 F01 H01 H02 H07 | A08 A12 C10 D01 D12 G03 G10 H12 -> A01 B03 B10 E01 E12 F03 H01 H05 |
| each level-7 annotation file byte-identical to the level-6 file of the turned name (e.g. level-7 `H01.tsv` = level-6 `A12.tsv`) | 8 of 8 | 8 of 8 |
| level-6 images found byte-identical under the turned name | 192 of 192 | 192 of 192 |
| layout, SOP, notes, instruction, Dockerfile, the three reference frames identical to level 6 | 8 of 8 | 8 of 8 |

2026-11-12, first version: `hcs-dose-response-s0001-l3` (211 files), `-s0001-l4`, `-s0001-l6` and
`-s0002-l6` (219 files each) regenerated without a difference. Levels 1 and 2 have no task directory in
`tasks/` to compare with; their code path is that of level 3. The agent image of both level-7 tasks is
derived from `alegen-hcs-dose-response-s0001-l4-agent:local` (`environment/Dockerfile` is byte-identical
to that of levels 3-6).

One slip of mine on 2026-11-12: while scoring dev runs on the host I ran `tests/score.py` of
`tasks/hcs-dose-response-s0001-l4` in place, which left a `tests/__pycache__/_scoring.cpython-313.pyc`
there from 22:11 to 22:33. I removed it; the directory is byte-identical to a fresh level-4 build again. A
verifier image built from that directory in those 22 minutes would contain the stray file (it is not read by
the Python 3.12 of the image).

## Known limitations and things I am unsure about (level 7)

- **Model runs.** I ran none. The maintainer's report on the first version is quoted at the top of this
  section; the revised tasks have not been run by a model as far as I know.
- **The manual counts are a second route to the decision point** (see "How it is discovered"): eight
  files whose names mean something else in the layout than what they hold, three of them per seed with
  fewer than half the cells of an untreated field. This is evidence in the data and follows from the
  decision that the lab counted the same physical fields. It makes the turn easier to notice than in the
  first version for an analyst who compares the manual counts with the layout; whether a model does that
  is untested.
- **Rejection rate.** About one seed in ten has no level 7 (see above). The thresholds of the check (half
  of the reference wells, four vehicle wells, factor 1.3, five series) are mine; seed 2 passes the vehicle
  criterion at exactly four. The check does not use the annotated fields.
- **A reference inhibitor that "did not work" is a competing explanation** for reference wells at vehicle
  level. The generator requires at least four off vehicle wells in addition, which that explanation does
  not cover, but an analyst may still stop at it. `controls_swapped_only` is the nearest thing in the
  battery; "exclude the odd controls and go on by name" is what `l6_reference` already does through its
  robust plate model.
- **Which naming an agent will choose is unknown,** and the verifier accepts both for each table. I see no
  wrong submission that gains from it, but I only tried the ones in the battery.
- **By-name output differs a lot between seeds**: all `inactive` (seeds 1, 3) or a few numbers, some below
  the tested range (seed 2). The second kind is easier to distrust from the table alone.
- **Pixel data are not turned** (see "What was planted"). The argument is distributional; the level-4
  counting frame is applied to the image as filed.
- **`task.toml` says "same inputs and planted truth as level 7" and `derived_from = ...-l7`.** That text
  comes from `descaffold.apply`, which treats the level-7 build as its own source; the same happens for
  the other two families and I did not change the core.
- **Focus QC is enforced through a tie-break** (levels 4, 6 and 7 alike; see the section on
  `naive_no_focus_qc`). With the outline threshold 0.15, which meets the SOP criterion on the annotated
  fields of seed 1 as well as 0.12 does, a submission without focus QC scores 0.825 instead of 0.544. It
  does not pass, but it is nearer to the pass line than the 0.15 margin of AUTHORING.md. Not in the
  battery; not fixable without touching the level-4 weights.
- **Gated on seeds 1 and 2 only**; seeds 3 and 10 were run by hand with the reference solution, its
  variant without focus QC and `l6_reference`. The second solution was not run on seeds 3 and 10.
