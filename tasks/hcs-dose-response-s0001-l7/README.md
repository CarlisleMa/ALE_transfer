# High-content screening plate: validated cell counts and compound potencies (level 7: pointers removed, new silent decision points)

Count the cells of a two-channel 96-well imaging plate with a procedure validated against the lab's manual counts, quantify the marker per cell, apply image QC and derive per-compound IC50/GI50 values under the SOP's definitions. Level 7: same inputs and planted truth as level 7; the text no longer says what to look out for (see the task card).

- family: `hcs_dose_response`  seed: `1`  level: `7`
- ALE subdomain: 3.3 Cell & Imaging Biology
- software: Cellpose, scikit-image, scipy, tifffile, Python
- a-priori difficulty: 0.576 (full-spectrum)  |  pass threshold: 0.9
- capabilities trained: validating a segmentation against manual annotation before using it, tuning a pretrained segmentation model or building a classical alternative to a stated quality bar, applying written counting rules (fragments, counting frame), plate-level QC before curve fitting, per-well image QC, dose-response modelling under a written potency definition, exact output contracts, checking the plate controls against the layout before using them

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
- **T_plate_orientation** - The plate lay on the microscope stage turned by 180 degrees. The acquisition software names the images by stage position, so the image pair named for well (row r, column c) shows plate well (row 9 - r, column 13 - c): A01 shows the well the layout calls H12. plate_layout.tsv is what the wet lab dispensed and is unchanged; no text mentions orientation, the stage holder, checking the controls or a mismatch. An analysis that joins images and layout by well name uses compound wells as vehicle wells and mid-series wells as reference wells, reads every dose series backwards and credits it to the compound of the opposite row. Nothing fails: counting, its validation against the manual counts and the focus QC are per image and are right, curves are fitted and potencies or `inactive` are reported. The lab counted the same eight physical fields as at level 4; each annotation file is the level-4 file under the name of the image that shows its field: A01.tsv, 38 cells (the layout calls A01 CPD-8550 at 0.0137 uM; the field is plate well H12, CPD-7827 at 30 uM); B06.tsv, 127 cells (the layout calls B06 reference; the field is plate well G07, vehicle); D08.tsv, 113 cells (the layout calls D08 vehicle; the field is plate well E05, vehicle); E03.tsv, 31 cells (the layout calls E03 reference; the field is plate well D10, CPD-4184 at 10 uM); F01.tsv, 18 cells (the layout calls F01 CPD-9526 at 0.0137 uM; the field is plate well C12, CPD-9251 at 30 uM); H01.tsv, 10 cells (the layout calls H01 vehicle; the field is plate well A12, CPD-8550 at 30 uM); H02.tsv, 66 cells (the layout calls H02 CPD-7827 at 0.0137 uM; the field is plate well A11, vehicle); H07.tsv, 34 cells (the layout calls H07 CPD-7827 at 0.37 uM; the field is plate well A06, CPD-8550 at 0.37 uM).
  - discoverable by: looking at the plate controls before they are used: a plate map of both readouts with the roles of the layout, or the control statistics (reference against vehicle wells, spread of the vehicle wells). A plate fits into a stage holder in exactly two orientations, and the same statistics with the layout turned by 180 degrees settle it. Measured on the planted values (cell number and mean marker level of every image; out-of-focus images left out). Layout joined by well name: marker per cell of the 7 in-focus wells called reference relative to the vehicle wells of their column 0.29-2.23 (median 1.74; 1 of 7 below half); the two wells called vehicle in a column differ in marker per cell by a factor of 1.60 (median) and at most 22.01 (6 of 11 columns beyond 1.3); cell number of the wells called vehicle 0.17-1.13 of the median of their plate area (inner wells / outer ring; 6 of 23 beyond a factor 1.3); wells called vehicle whose planted response is beyond a factor 1.3 on either readout 11 of 23; dose series whose response falls / rises with the concentration (rank correlation beyond 0.7; eight series per readout): cell number 1 / 1, marker per cell 0 / 0. Layout turned by 180 degrees: marker per cell of the 8 in-focus wells called reference relative to the vehicle wells of their column 0.05-0.05 (median 0.05; 8 of 8 below half); the two wells called vehicle in a column differ in marker per cell by a factor of 1.04 (median) and at most 1.14 (0 of 11 columns beyond 1.3); cell number of the wells called vehicle 0.92-1.06 of the median of their plate area (inner wells / outer ring; 0 of 23 beyond a factor 1.3); wells called vehicle whose planted response is beyond a factor 1.3 on either readout 0 of 23; dose series whose response falls / rises with the concentration (rank correlation beyond 0.7; eight series per readout): cell number 4 / 0, marker per cell 5 / 1. The reference wells are the plainest sign: the reference inhibitor is described in the assay notes as an inhibitor of the marker pathway at 10 uM. A second sign, also in the data and in no text: the manual counts. The annotation files carry names that mean something else in the layout than the fields they hold; 3 of the 8 files hold fewer than half the cells of an untreated field of their plate area (120 cells inner wells, 66 outer ring), 2 of them under a name the layout calls vehicle or reference (E03: 31 cells, layout reference; F01: 18 cells, layout CPD-9526 at 0.0137 uM; H01: 10 cells, layout vehicle).
  - trains against ALE failure mode: Verification / skipped QC step

## Judge

**Level 7.** Inputs and planted truth are those of `hcs-dose-response-s0001-l7`; only text was taken away:

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

**Level-7 decision point (the plate lay on the stage the other way round).** The image pair named for a well
shows the well at the opposite position of the plate (row 9 - r, column 13 - c); `plate_layout.tsv` is in
plate coordinates and unchanged. The planted truth is that of the level-4 plate of the same seed with nothing
redefined: potencies are per compound and are the level-4 values; cell number, focus state and manual
annotation belong to an image and are stored under the name of the image. Gates, components, weights and
tolerances are those of level 4: `tests/score_l4.py` is the level-4 verifier, byte for byte, and
`tests/score.py` only runs it. No deliverable asks for the orientation.

**The two per-well tables are accepted under either naming of their rows.** An analyst who found the turn can
file a measurement under the name of the image (and annotation file) it was made on, which is how the truth
is stored, or under the name of the plate well that image shows, which is what a reader who joins the table
with the layout needs. "One row for every well of `plate_layout.tsv`" and "`field` is the well id" can be
read either way once the turn is known, no text says which is meant, and none may be added at this level. The
verifier therefore runs the level-4 logic under the four combinations (each of `well_measurements.tsv` and
`validation.tsv` by image or by plate well) and the best total is the score; under every combination a
validation row is compared with the planted count of the image it refers to and with the `n_cells` reported
for that same image. `info.well_naming` in details.json records the combination used and all four totals.
Potencies (0.65 of the score) carry no well names, so the naming cannot rescue an analysis made with the
layout the wrong way round. Measured: the reference result scores 1.0 by image name and 1.0 when both tables
are filed by plate well; the by-name submissions below score 0.35 under the naming they use and 0.0 under
the three others.

**What the layout by name costs** (measured, seed 1 / seed 2). `l6_reference`, the level-4 reference
solution unchanged: 0.35 / 0.35. It keeps the three per-image components (every well without annotation
within the count tolerance, validation table 8 of 8, focus flags exact) and earns nothing on either potency
component: marker 4 / 3 of 8 compounds right, cell number 5 / 4 of 8, where a component needs a mean credit
above 0.8. On seed 1 it reports `inactive` for all 16 entries, so the entries it has right are the
compounds that are inactive anyway. On seed 2 it reports 12 times `inactive` and four numbers: marker IC50s
of 0.0013 and 0.0027 uM for compounds whose planted IC50s are 5.3 and 3.5 uM, 17 uM for a compound that is
inactive, and a GI50 of 0.0019 uM for a compound that is inactive on cell number; the three small values lie
below the lowest tested concentration (0.0137 uM). The curves behind these tables are not smooth inhibition
curves: a series holds the wells of the opposite row in reverse order, with that row's vehicle and reference
wells in between, so responses rise with the nominal concentration (seed 2, one series: 0.04, 0.08, 0.23,
0.54, 0.81, 0.95, 1.02, 0.90) or jump (seed 1, one series: 0.09, 1.47, 0.75, 0.08, 1.59, 1.53, 1.44, 0.94);
the SOP model, whose upper asymptote is fixed at the vehicle level, turns a rising series into `inactive`
or into a potency at the lower end of its search range. The table alone does not show any of this.
`layout_by_name_everything_else_right` (planted cell numbers, marker levels and focus flags; reference plate
model and fit with the layout by name): 0.35 / 0.35, the same potency tables as `l6_reference` on both
seeds. `controls_swapped_only` (layout by name, the roles of vehicle and reference wells exchanged): 0.35 /
0.35 (seed 1: 14 numeric potencies, marker 1 of 8 right; seed 2: marker 3 of 8, cell number 5 of 8).
`prior_only` shows the floor from the other side: `inactive` everywhere is right for four of eight marker
and five of eight count entries on every plate of this family, which earns no potency credit.

**Routes that reach the truth** (probes): the reference solution finds the orientation from the control wells
(reference wells inhibited, vehicle wells alike), the second solution from the direction of the dose series
without using the controls; both then apply the layout turned and report the per-well tables by image name.
`per_well_tables_by_plate_well` and `well_table_by_plate_well_validation_by_image` are the planted truth
filed the other way.

**The annotated fields are those of level 4.** The lab counted the same eight physical fields (the plate
wells chosen at level 4 to cover dense and sparse untreated fields, fragment-rich fields and fields with
enlarged nuclei). The annotation files are the level-4 files, byte for byte, under the names of the images
that show those fields; `validation.tsv` lists exactly these eight fields (by image name or by plate well,
see above). The counting part of the task is therefore that of levels 4 and 6. Measured on seeds 1, 2, 3 and
10: the reference solution selects the outline threshold it selects on the level-4 inputs (0.12; seed 1 was
compared directly, same worst annotated field at 0.76 of the SOP limit), and `naive_no_focus_qc` scores
0.544 / 0.544 / 0.669 / 0.581, under its level-4 ceiling of 0.75. (In the first version of this level the
annotation files kept their level-4 names and held whatever image lay under them; the tuning then chose
another threshold on some seeds, that attack reached 0.825 and its ceiling had been raised to 0.85. Both
are gone.)

**A side effect that is evidence, not text.** The annotation files now carry names under which the layout
lists something else than what was counted: the entry of the decision point above gives, for the eight
files of this seed, the manual count, what the layout calls the name and which plate well the field is.
Seed 1 / seed 2: 3 / 3 of the 8 files hold fewer than half the cells of an untreated field of their plate
area (about 120 cells in inner wells, 66 in the outer ring), 2 / 1 of them under a name that the layout
calls vehicle or reference (seed 1: `H01` 10 cells, layout vehicle; `E03` 31 cells, layout reference; seed
2: `H01` 8 cells, layout vehicle). An analyst who compares the manual counts with the layout sees it; the
assay notes say how the fields were chosen (untreated wells, highest and intermediate concentrations) and
are the level-4 notes.

**Seeds without a level 7.** A seed whose control wells sit symmetrically under the turn would show
well-behaved controls in both orientations. The generator refuses such seeds (41 of 400 surveyed): it
requires, on the planted values, that with the layout by name fewer than half of the wells called reference
are inhibited and at least four wells called vehicle are off by more than a factor 1.3, and that the
balance of falling over rising dose series separates the two orientations by at least five series. The
check does not use the annotated fields.

## Verifier attack suite (must score low)

- `l6_reference` (max 0.75): the reference solution of levels 4-6, unchanged, on the level-7 inputs: every decision point of the lower levels is handled, images and layout are joined by well name
- `layout_by_name_everything_else_right` (max 0.75): T_plate_orientation ignored and nothing else wrong: planted cell numbers, planted marker levels and planted focus flags of every image, validation table from the planted counts, potencies as the reference plate model and curve fit give them when the layout is joined by well name
- `controls_swapped_only` (max 0.75): the odd controls are noticed but the layout is not turned: the reference solution with the layout joined by well name and the roles of the vehicle and the reference wells exchanged (the wells called reference look untreated and are used as vehicle wells)
- `l3_reference` (max 0.75): the level-3 reference solution unchanged (Cellpose nuclei model with the level-3 settings, every mask >= 20 px counted, level-3 plate model / QC / fit) plus a truthful validation table (level-7 inputs, layout joined by well name)
- `l3_reference_classical` (max 0.75): the level-3 classical route unchanged (threshold + distance-transform watershed; the route the pilot agent took) plus a truthful validation table (level-7 inputs, layout joined by well name)
- `l3_independent_solution` (max 0.75): the level-3 second solution unchanged (LoG-seeded watershed) plus a truthful validation table (level-7 inputs, layout joined by well name)
- `counts_scaled_to_annotation` (max 0.75): level-3 classical route; the counts are then multiplied by one plate-wide factor that makes the annotated fields agree in total (level-7 inputs, layout joined by well name)
- `counts_scaled_to_annotation_log` (max 0.75): level-3 LoG route with the same plate-wide factor (level-7 inputs, layout joined by well name)
- `manual_counts_copied` (max 0.75): level-3 classical route; the counts of the annotated wells are replaced by the manual counts so that the validation table looks perfect (level-7 inputs, layout joined by well name)
- `naive_fragments_counted` (max 0.75): reference solution, everything else correct, but nuclear fragments are counted as cells (plate controls checked, layout applied turned)
- `naive_free_plateau_fit` (max 0.75): reference solution, everything else correct, but the curves are fitted with free upper and lower plateaus and the fitted EC50 is reported (plate controls checked, layout applied turned)
- `naive_no_focus_qc` (max 0.75): reference solution, but no focus QC: every well passes and is used (plate controls checked, layout applied turned)
- `naive_whole_well_readout` (max 0.75): reference solution, but the marker readout is the mean background-subtracted intensity of the field instead of the intensity per cell (plate controls checked, layout applied turned)
- `naive_no_plate_model` (max 0.75): reference solution, but responses are normalised to the plate mean of the vehicle wells (plate controls checked, layout applied turned)
- `naive_no_background` (max 0.75): reference solution, but the marker per cell is the raw ring intensity above the camera offset (plate controls checked, layout applied turned)
- `naive_global_background` (max 0.75): reference solution, but one plate-wide background constant is subtracted from every well (plate controls checked, layout applied turned)
- `prior_only` (max 0.1): no analysis: every well 100 cells and QC pass, every compound inactive, validation table claims 100 cells per field
- `all_potencies_beyond_tolerance` (max 0.4): planted counts and QC flags, every active potency 0.45 log10 units too high
- `one_potency_beyond_tolerance` (max 0.75): planted truth, but one marker IC50 0.45 log10 units too high
- `counts_beyond_tolerance` (max 0.75): planted potencies and QC flags, every cell count 20% too low, validation table truthful about it

## Must-pass probes (guard against false negatives)

- `planted_truth` (min 0.999): submission written directly from the planted cell numbers, focus states and potencies; validation table with the planted counts of the annotated fields (per-well rows under the image names)
- `within_tolerance` (min 0.99): every count off by 7%, every potency off by 0.15 log10 units, one in-focus well flagged (per-well rows under the image names)
- `independent_solution` (min 0.99): second implementation, no code shared with the reference solution: Cellpose nuclei model with settings chosen on the annotated fields, as at level 4; the orientation of the plate is found from the direction of the dose series (rank correlation of both readouts with the concentration, readouts normalised without the control wells), not from the controls
- `per_well_tables_by_plate_well` (min 0.999): the turn is recognised and both per-well tables are filed by the plate well each image shows: planted truth with every row under the opposite well name. Accepted like the tables by image name: the verifier scores each table under both namings
- `well_table_by_plate_well_validation_by_image` (min 0.999): the well table is filed by plate well, the validation table keeps the names of the annotation files: planted truth. The tie between the two tables (pipeline_count = n_cells of the same image) holds under that naming
