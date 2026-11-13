# Level-7 spec: `flow_compensation` - a disturbance during acquisition that only the Time axis shows

Why: gemini-3.8-flash scored 1.0 on levels 4, 5 and 6 of this family. With every pointer removed it still
read the detector voltages from the FCS headers and the reagent lots from the control table, and derived the
right matrices. What is left to exploit is the shape that defeated every model on somatic_tn and
trio_liftover: the routine workflow runs without an error, its result looks plausible, and it is wrong
because one non-default professional step is missing that nothing in the text asks for.

## New decision point (silent)
**T_acquisition_disturbance.** In two of the six samples (one in each voltage block), a contiguous part of
the run (12-25 % of the events, somewhere in the middle of the acquisition as given by the `Time` parameter)
was recorded during a fluidics disturbance. During that window the detectors of one laser line (choose the
line that carries a lineage marker and the HLA-DR reagent, e.g. the yellow-green detectors) record a
fraction of the true signal (about 0.3-0.45, with extra noise); scatter and the detectors of the other
lasers are unaffected, and the event rate changes little. Cells in the window therefore look negative or
dim for the markers on that laser.
- Naive outcome: gating all events gives lineage and activation percentages that are outside the tolerance
  in the affected samples; nothing fails, and the numbers are still plausible for PBMC.
- Expert practice: signal-versus-Time QC is routine (flowAI, PeacoQC, the Time gate of every analysis
  package). The step in the affected detectors over a time window is obvious in a plot or a rolling median;
  the window is excluded (or its signal corrected) before gating.
- Truth: population membership is the planted label of every event, disturbed or not. Percentages estimated
  from the undisturbed part of the file are unbiased and must be within tolerance; a submission that repairs
  the window is equally right. Populations that do not depend on the affected detectors (cells, singlets,
  live, lymphocytes, monocytes) can still be counted on all events; make sure the scored counts are of that
  kind, and say so in the judge explanation.
- Nothing in the instruction, the SOP or the notes may mention time, stability, the disturbance or the step.
  `Time` is already listed among the recorded parameters; that stays as it is.
- `Trap.discoverable_by`: name the QC (median of the affected detectors in consecutive time bins) and give
  the measured size of the step and the number of bins it spans.

## What stays
Everything of level 4 (data model, the other decision points, deliverables, scorer, tolerances). The text
of level 7 is produced by the pipeline: the family returns a level-7 build whose instruction, SOP and notes
are the level-4 texts unchanged, and `alegen/descaffold.py` removes the scaffolding and the incident reports
exactly as for level 6. Do not edit those texts and do not add sentences to them.

## Battery (in addition to the level-4 battery, on the level-7 data)
- `l6_reference` (required by the validator at level 7): the reference solution of level 4 without the new
  step, on the level-7 inputs. Must score at most pass - 0.15.
- `time_gate_everything_after`: discards everything from the start of the disturbance to the end of the file
  (valid: unbiased) - a probe, must pass.
- `repair_instead_of_exclude`: rescales the affected detectors inside the window - a probe, must pass if you
  can make it exact enough; otherwise document why not.
- `exclude_wrong_window`: a time gate of the right length in the wrong place - attack, must fail.
- The independent second solution must find the window by a different statistic than the oracle.

## Measurements to report (two seeds)
Step size and width of the window; score of `l6_reference`; how many of the scored items leave tolerance
without the time gate; that levels 1-6 of the family regenerate byte-identically.
