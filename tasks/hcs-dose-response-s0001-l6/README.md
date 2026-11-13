# High-content screening plate: validated cell counts and compound potencies (level 6: pointers removed)

Count the cells of a two-channel 96-well imaging plate with a procedure validated against the lab's manual counts, quantify the marker per cell, apply image QC and derive per-compound IC50/GI50 values under the SOP's definitions. Level 6: same inputs and planted truth as level 4; the text no longer says what to look out for (see the task card).

- family: `hcs_dose_response`  seed: `1`  level: `6`
- ALE subdomain: 3.3 Cell & Imaging Biology
- software: Cellpose, scikit-image, scipy, tifffile, Python
- a-priori difficulty: 0.564 (full-spectrum)  |  pass threshold: 0.9
- capabilities trained: validating a segmentation against manual annotation before using it, tuning a pretrained segmentation model or building a classical alternative to a stated quality bar, applying written counting rules (fragments, counting frame), plate-level QC before curve fitting, per-well image QC, dose-response modelling under a written potency definition, exact output contracts

## Undisclosed decision points (not shown to the agent)

- **Q_validated_counts** - Counting quality decides the task. Nuclei overlap without an intensity dip (centre distance down to 9.5 px for 12-15 x 10-12 px nuclei), chromatin is textured, every field holds bright nuclear fragments (about 12% of the seeded cells, many more where cells were lost) that the SOP does not count, nuclei are cut by the field border under a counting-frame rule, and two compounds enlarge the nuclei up to 1.3-1.5x. The level-3 methods fail the SOP criterion on the annotated fields and on the other wells; the errors have different signs in dense, sparse and treated fields, so no plate-wide factor repairs them.
  - discoverable by: SOP section 2 states the acceptance criterion; annotations/ holds the manual point annotation of eight fields (dense and sparse untreated fields, fragment-rich and enlarged-nuclei fields); comparing any counting procedure with them, or overlaying its objects on the annotated points, shows which nuclei are lost and which objects are counted wrongly
  - trains against ALE failure mode: Verification / result not checked against the available control
- **P_potency_shape** - The count response of CPD-9251 is a deep partial response (plateau 0.30): the SOP potency (0.5 crossing) is 0.32 log10 units above the inflection point that a free four-parameter fit reports. The marker per cell of CPD-8550 rises with the fraction of cells lost (factor 1.76 when all are lost); a fit with free plateaus returns an EC50 for it although nothing is inhibited.
  - discoverable by: SOP section 5 fixes the upper asymptote at the vehicle level and defines the potency as the 0.5 crossing and `inactive`; plotting the responses of each series shows the plateau at about a third and the rising series
  - trains against ALE failure mode: Understanding / textbook default instead of the written definition
- **touching_nuclei** - Level 2 decision point, sharpened by Q_validated_counts: a global threshold merges the nuclei of a colony.
  - discoverable by: overlay of the threshold mask on any vehicle field; comparison with the manual counts
  - trains against ALE failure mode: Approach / wrong strategy
- **focus_qc** - 5 fields (C04, C05, E09, G06, G08) are out of focus (extra blur sigma 4.5-6 px): two neighbouring series wells just below the marker potency of one compound, two just below the count potency of another, one vehicle well.
  - discoverable by: SOP section 3 defines the QC failure; a signal-normalised sharpness metric of channel 1 puts these wells far below all others
  - trains against ALE failure mode: Verification / skipped QC step
- **readout_confound** - Two compounds reduce cell number without reducing the marker per cell; the field intensity follows cell number.
  - discoverable by: SOP section 1 defines R2 per cell; field intensity tracks cell number along those series
  - trains against ALE failure mode: Understanding / wrong readout
- **plate_position_effects** - Marker staining falls log-linearly to 0.40x along the columns (dose runs along the columns) and the outer wells hold 0.54x the cells of inner wells.
  - discoverable by: plate map of the 24 vehicle wells (3 per row, 2 per column) for each readout; assay notes give the pipetting order and the incubation conditions as facts
  - trains against ALE failure mode: Verification / skipped QC step
- **marker_background** - Diffuse marker-channel background of about a third of the vehicle signal that differs between wells; one marker inhibitor fluoresces in the marker channel at high concentration.
  - discoverable by: SOP section 1 excludes signal that does not come from the cells from R2; cell-free regions of channel 2 are far above the dark frame, differ between wells and rise with concentration in one series
  - trains against ALE failure mode: Approach / skipped QC step

## Judge

**Level 6.** Inputs and planted truth are those of `hcs-dose-response-s0001-l4`; only text was taken away:

- S5: assay notes without the four remarks that explain where the plate artefacts come from (dry humidity tray and shelf position; stain added column by column and removed at once; one fixed focus offset and no review at the microscope; unwashed wells with floating material). The artefacts are visible in the vehicle wells of the plate map and in the images. The SOP, the annotated fields and the scorer are unchanged

Fairness no longer rests on a stated criterion but on three things that are checked for the level-4 task this one is derived from: the correct handling is standard practice of the field, it is discoverable from the data with a standard QC step (see the decision points above), and two independent solutions agree with the planted truth.

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

## Verifier attack suite (must score low)

- `l3_reference` (max 0.75): the level-3 reference solution unchanged (Cellpose nuclei model with the level-3 settings, every mask >= 20 px counted, level-3 plate model / QC / fit) plus a truthful validation table
- `l3_reference_classical` (max 0.75): the level-3 classical route unchanged (threshold + distance-transform watershed; the route the pilot agent took) plus a truthful validation table
- `l3_independent_solution` (max 0.75): the level-3 second solution unchanged (LoG-seeded watershed) plus a truthful validation table
- `counts_scaled_to_annotation` (max 0.75): level-3 classical route; the counts are then multiplied by one plate-wide factor that makes the annotated fields agree in total
- `counts_scaled_to_annotation_log` (max 0.75): level-3 LoG route with the same plate-wide factor
- `manual_counts_copied` (max 0.75): level-3 classical route; the counts of the annotated wells are replaced by the manual counts so that the validation table looks perfect
- `naive_fragments_counted` (max 0.75): reference solution, everything else correct, but nuclear fragments are counted as cells
- `naive_free_plateau_fit` (max 0.75): reference solution, everything else correct, but the curves are fitted with free upper and lower plateaus and the fitted EC50 is reported
- `naive_no_focus_qc` (max 0.75): reference solution, but no focus QC: every well passes and is used
- `naive_whole_well_readout` (max 0.75): reference solution, but the marker readout is the mean background-subtracted intensity of the field instead of the intensity per cell
- `naive_no_plate_model` (max 0.75): reference solution, but responses are normalised to the plate mean of the vehicle wells
- `naive_no_background` (max 0.75): reference solution, but the marker per cell is the raw ring intensity above the camera offset
- `naive_global_background` (max 0.75): reference solution, but one plate-wide background constant is subtracted from every well
- `prior_only` (max 0.1): no analysis: every well 100 cells and QC pass, every compound inactive, validation table claims 100 cells per field
- `all_potencies_beyond_tolerance` (max 0.4): planted counts and QC flags, every active potency 0.45 log10 units too high
- `one_potency_beyond_tolerance` (max 0.75): planted truth, but one marker IC50 0.45 log10 units too high
- `counts_beyond_tolerance` (max 0.75): planted potencies and QC flags, every cell count 20% too low, validation table truthful about it

## Must-pass probes (guard against false negatives)

- `planted_truth` (min 0.999): submission written directly from the planted cell numbers, focus states and potencies; validation table with the planted counts of the annotated fields
- `within_tolerance` (min 0.99): every count off by 7%, every potency off by 0.15 log10 units, one in-focus well flagged
- `independent_solution` (min 0.99): second implementation, no code shared with the reference solution: Cellpose nuclei model with settings chosen on the annotated fields (the level-3 settings are tried first and rejected), fragments by mask area, median background far from nuclei, edge-gradient focus metric, least-squares plate model, scipy least_squares fit
