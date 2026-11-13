# Level-7 spec: `hcs_dose_response` - the plate sat on the stage the other way round

Why: gemini-3.8-flash scored 1.0 on level 4 and, given enough turns (300), 1.0 on the pointer-free level 6 of
this family; its two level-6 attempts at 150 turns ended at the turn limit without a submission, which is a
workload effect and not a missed decision. The shape that did stop it elsewhere (crispr_screen level 6,
flow_compensation level 7) is: a table the lab wrote down does not describe the data, no text says so, and
the data show it to anyone who looks at the controls before using them. This family has that lever in its
plate layout.

## New decision point (silent)
**T_plate_orientation.** The plate was put on the microscope stage turned by 180 degrees. The acquisition
software names images by stage position, so the image pair named for well (row r, column c) shows plate well
(row 9 - r, column 13 - c) (rows A-H = 1-8): `A01` shows the well the layout calls `H12`, and so on.
`plate_layout.tsv` is what the wet lab dispensed (plate coordinates) and is the level-4 file unchanged; so
are the notes.
- Naive outcome: an analysis that joins images and layout by well name uses compound wells as vehicle wells
  and mid-series wells as reference wells, reads every dose series backwards and credits each series to the
  compound of the opposite row. Nothing fails: counting, its validation against the manual counts and image
  QC are per image and are unaffected, and curves are fitted and potencies (or `inactive`) are reported.
- Expert practice: plate controls are looked at before they are used (plate heat map of both readouts,
  control statistics). Wells labelled reference that look like vehicle, wells labelled vehicle with few cells
  and little marker, and readouts that rise with the nominal concentration are found at once; a plate has
  exactly two orientations in a stage holder, and under the other one every vehicle well is at vehicle
  level, every reference well is inhibited and every series falls with concentration. The layout is then
  applied in that orientation.
- Truth: unchanged in meaning. The plate is the level-4 plate of the same seed (same cells, same compounds,
  same planted potencies, same position effects and out-of-focus wells - all of these belong to the plate
  and turn with it). Per-image truths (`n_cells`, `marker_per_cell`, `qc_pass`, manual annotations) are
  keyed by image name and therefore move with the image; `dose_response.tsv` truth is per compound and is
  the level-4 truth. Flat-field and dark frames belong to the optical path and do not change.
- Annotated fields: the lab picks fields by the layout it believes in, and counts what the image shows.
  Keep the level-4 choice of well names; the annotation files hold the cells of the images that now carry
  those names.
- Nothing in the instruction, the SOP or the notes may mention orientation, the stage holder, checking the
  controls or any mismatch between layout and images. No deliverable asks for the orientation.
- `Trap.discoverable_by`: name the QC and give measured numbers for both orientations (for example: spread
  of the vehicle wells, separation of reference from vehicle wells on the marker readout, sign of the rank
  correlation between concentration and response per series).

## What stays
Everything of level 4 (image model, all level-4 decision points, deliverables, scorer, tolerances, software
in the image - the Dockerfile must stay byte-identical so that the task is built on top of the existing
image). The text of level 7 is produced by the pipeline: the family returns a level-7 build
(`generate(seed, 7, workdir)`, `level == 7`, `task_id` ending in `-l7`) whose instruction, SOP and notes are
the level-4 texts unchanged, and `alegen/descaffold.py` then applies the level-6 edits. Do not edit those
texts and do not change `descaffold.py` unless an edit fails to apply (then report it).

## Battery (in addition to the level-4 battery, on the level-7 data)
- `l6_reference` (required by the validator at level 7): the level-4 reference solution unchanged (layout
  joined by well name) on the level-7 inputs. Must score at most pass - 0.15.
- `layout_by_name_everything_else_right`: planted per-image truths, potencies as the reference curve fit
  gives them with the layout joined by name - attack, must fail.
- `controls_swapped_only`: a submission that notices the odd controls but only re-labels vehicle and
  reference wells without turning the layout - attack, must fail; report its score.
- Scorer: unchanged if possible. If per-well rows can be read two ways (by image name, or by plate well),
  decide from the output contract which reading the contract fixes; if the contract does not fix it, accept
  both (as amplicon_16s level 7 does) and say so in the task card. Do not change tolerances.
- The independent second solution must find the orientation by a different statistic than the oracle.

## Measurements to report (two seeds)
The numbers behind `discoverable_by`; the score of `l6_reference` and which components it loses; what the
by-name analysis outputs (how many compounds `inactive`, how far the reported potencies are from the
planted ones, whether its curves would look acceptable to a reader); whether a symmetric layout or a
symmetric position effect could make the two orientations indistinguishable for some seed (then the
generator must reject such draws); that levels 1-4 and 6 regenerate byte-identically.
