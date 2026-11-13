# 8-colour PBMC flow cytometry: per-sample spillover matrices, lineage and FMO-gated activation frequencies (level 7: pointers removed, new silent decision points)

Derive the spillover matrix that applies to every sample of an 8-colour PBMC session whose detector settings and control reagents changed, and report the lineage and FMO-gated HLA-DR population frequencies of every sample from raw FCS 3.1 files. Level 7: same inputs and planted truth as level 7; the text no longer says what to look out for (see the task card).

- family: `flow_compensation`  seed: `2`  level: `7`
- ALE subdomain: 3.3 Cell & Imaging Biology
- software: FlowKit, FlowIO, fcsparser, scikit-learn, Python, numpy/scipy/pandas
- a-priori difficulty: 0.7 (last-exam)  |  pass threshold: 0.9
- capabilities trained: instrument-file metadata inspection, compensation from single-stain controls, reasoning about detector gain and spillover, control selection by reagent lot, FMO-based gating, self-check against an in-data control, hierarchical gating with per-sample QC, exact output contracts, acquisition QC along the Time parameter

## Undisclosed decision points (not shown to the agent)

- **compensation_required** - Nothing was compensated at acquisition; uncompensated, CD3-BV421 spills into the viability detector, CD4-PE into the CD8 detector and CD56-APC into the CD19 detector at the level of true positives.
  - discoverable by: acquisition_notes.md: compensation setup not run; diagonal populations in any V450-A vs V525-A or YG586-A vs B710-A plot
  - trains against ALE failure mode: Understanding / domain knowledge gap
- **doublets** - 6-12% of cell events are aggregates of two cells; their summed fluorescence inflates DP_T and the lymphocyte count.
  - discoverable by: gating_sop.md lists a singlets population; FSC-A vs FSC-H shows a second population off the diagonal; notes say tubes were not strained
  - trains against ALE failure mode: Approach / skipped QC step
- **dead_cells** - 8-18% of the single cells are dead and take up every antibody non-specifically.
  - discoverable by: gating_sop.md defines live as viability dye negative; notes report 80-92% post-thaw viability
  - trains against ALE failure mode: Approach / skipped QC step
- **per_sample_shift** - PBMC_D507 was recorded with the FSC voltage lowered (gain x0.24).
  - discoverable by: acquisition_notes.md states the re-acquisition; $P1V/$P2V differ in that file
  - trains against ALE failure mode: Verification / per-sample QC
- **P_embedded_matrix** - Every file carries a well-formed, non-identity $SPILLOVER (APPLY COMPENSATION = TRUE) from the experiment template of 14-MAR-2024; the voltages were re-optimised later, so it is wrong for this session (PerCP-Cy5.5 -> YG780-A 0.093 instead of 0.205).
  - discoverable by: acquisition_notes.md gives the dates of the template and of the voltage re-optimisation and says that the compensation setup was not run; the single-stain controls of the session give other coefficients; with the embedded matrix the lineage gates break down (a population of the SOP's acceptance check is missing in the FMO tube)
  - trains against ALE failure mode: Prior decoupled from truth / trusting embedded metadata
- **P_voltage_change** - The YG780-A voltage was raised from 641 V to 716 V (gain x2.20) after PBMC_D507: PBMC_D124, PBMC_D507 need the session-start matrix, PBMC_D523, PBMC_D583, PBMC_D769, PBMC_D975 a second one. Only the PE-Cy7 bead tube and the FMO tube were recorded again; the spill of every other fluorochrome into YG780-A scaled by the same gain factor and has to be derived (bead brightness ratio of the two PE-Cy7 tubes, ratio of the two PE-Cy7 rows, or the FMO tube itself).
  - discoverable by: acquisition_notes.md states the voltage change and which tubes were recorded again; $P9V differs between the files; the SOP's acceptance check on the re-recorded FMO tube fails when only the PE-Cy7 row is replaced (CD8_T median far above DN_T)
  - trains against ALE failure mode: Prior decoupled from truth / one matrix per experiment
- **P_tandem_lot** - The PerCP-Cy5.5 bead control was stained from lot 3339183, the samples from lot 3651405; the bead lot spills 0.124 less into YG780-A. The cell-based single stain from the cocktail vial gives the row that applies to the samples.
  - discoverable by: reagent_lot columns of panel.tsv and controls.tsv differ for that bead control and agree for the cell-based one; acquisition_notes.md says when and from which vials the beads and the extra tube were stained; the two PerCP-Cy5.5 controls give different rows; the SOP's acceptance check on the FMO tube fails with the bead-derived row
  - trains against ALE failure mode: Prior decoupled from truth / beads as default compensation control
- **Q_fmo_gate** - HLA-DR on T cells has no valley; the SOP defines the gate on the FMO control (99th percentile, same parent, compensated) and states an acceptance criterion for the compensation on the same tube (medians of five lymphocyte populations within 50 units). A gate from the unstained control, from uncompensated data or from a valley search gives a plausible but out-of-band percentage; a residual spillover of several hundred units widens the FMO background and loses dim positives. The default matrices (embedded, bead lot, session-start column) miss the criterion by a factor 16 or more; the analyst has to find the cause or tune the matrix on the FMO tube.
  - discoverable by: gating_sop.md sections 1-3; the FMO tubes are listed in controls.tsv; the acceptance criterion lets the analyst measure the residual spillover on cells and says that a failing compensation has to be corrected
  - trains against ALE failure mode: Quality bar with an in-data control
- **T_acquisition_disturbance** - A contiguous part of two runs was recorded during a fluidics disturbance: PBMC_D124 (events 24065-30707 of 41059, 16% of the file, Time 7523-9791, signal x0.41, event rate x0.92) and PBMC_D523 (events 20870-27714 of 36071, 19% of the file, Time 6523-8733, signal x0.35, event rate x0.97). Inside the window the two detectors of the 561-nm laser, YG586-A (CD4 PE) and YG780-A (HLA-DR PE-Cy7), read a fraction of the signal (event-to-event jitter 1-2%); scatter and the other six detectors are unaffected. Gated with all events, the PE spill into B710-A is under-compensated inside the window, so its CD4 T cells are counted as CD4+CD8+ (CD4_T low and DP_T high by about the share of the window), and the full PerCP-Cy5.5 spill is subtracted from a YG780-A signal that was recorded too low, so no CD8 T cell of the window passes the HLA-DR gate (activated_CD8_T low by the share of the window). Nothing fails. No text mentions time, stability or the disturbance; `Time` is one of the recorded parameters.
  - discoverable by: signal-versus-Time QC as done by flowAI, PeacoQC or a manual Time gate: median of every fluorescence detector in consecutive bins of 500 events ordered by Time. PBMC_D124: bin medians YG780-A 1618 -> 680 (x0.42, 15 SD) and YG586-A 1790 -> 767 (x0.43, 7 SD) in 12 consecutive bins of 82 (SD = robust standard deviation of the bin medians outside the window); bins outside the window whose median is nearer to the level inside: 0 of 68 (YG780-A), 3 of 68 (YG586-A); 90th percentile per bin YG586-A 33115 -> 13502 (x0.41, 19 SD); the bin medians of the other six detectors differ between inside and outside by at most 0.3 SD. PBMC_D523: bin medians YG780-A 3898 -> 1411 (x0.36, 16 SD) and YG586-A 1794 -> 619 (x0.35, 9 SD) in 13 consecutive bins of 72 (SD = robust standard deviation of the bin medians outside the window); bins outside the window whose median is nearer to the level inside: 0 of 57 (YG780-A), 0 of 57 (YG586-A); 90th percentile per bin YG586-A 33406 -> 11626 (x0.35, 29 SD); the bin medians of the other six detectors differ between inside and outside by at most 0.4 SD. The same window shows in any CD4-versus-CD8 plot of the T cells of these samples as a third population (CD4 lower by the factor of the window, CD8 intermediate) that exists only in that Time range.
  - trains against ALE failure mode: Approach / skipped QC step

## Judge

**Level 7.** Inputs and planted truth are those of `flow-compensation-s0002-l7`; only text was taken away:

- S1: the compensation acceptance criterion of the gating SOP (section 3: medians of five lymphocyte populations in the FMO tube within 50 detector units)
- S5: acquisition notes without the three incident reports (the template compensation that was never re-run; the detector voltage raised in mid-session and which controls were recorded again; the sample re-acquired with a lower FSC voltage) and without the remarks on which vials the controls were stained from. Voltages are in the FCS headers, reagent lots in controls.tsv and panel.tsv

Fairness no longer rests on a stated criterion but on three things that are checked for the level-4 task this one is derived from: the correct handling is standard practice of the field, it is discoverable from the data with a standard QC step (see the decision points above), and two independent solutions agree with the planted truth.

**Where the truth comes from.** Every event is simulated from a hidden class label (debris, aggregate
of two cells, dead cell, CD4/CD8/DP/DN T cell, B cell, NK cell, other lymphocyte, monocyte). The
reference `truth.json` holds, per sample, the planted spillover matrix of that sample (sample reagent
lots, detector settings of the sample's acquisition block; the block-2 matrix is the block-1 matrix
with the HLA-DR detector column multiplied and the PE-Cy7 row divided by the planted gain factor) and
the number of events of each lineage population counted from the labels. The two FMO-gated
populations are defined operationally by the SOP; their reference is that rule applied to planted
quantities only: signals compensated with the planted matrix, parent population = hidden label,
threshold = 99th percentile of the label-defined parent in the FMO tube of the same block. No gate or
matrix estimate of the reference pipeline is used as reference.

**Hard gates (score 0).** Missing/unparseable file; `spillover.csv` without the columns `sample`,
`fluorochrome` and the eight detectors, not exactly one row per sample x fluorochrome, non-numeric
entries; `frequencies.tsv` not holding exactly one row per sample x SOP population, wrong `parent`,
`count` not a non-negative integer, `percent_of_parent` outside [0, 100].

**Components.**

| weight | component | metric | full credit | zero |
|---|---|---|---|---|
| 0.20 | spillover_matrices | maximum absolute error over the 8 x 8 entries of every sample's matrix | <= 0.02 | >= 0.05 |
| 0.25 | lineage_every_sample | per sample, fraction of the 12 lineage populations whose `percent_of_parent` is within max(0.5 points, 3% of the reference); **minimum** over samples | 1.0 | <= 10/12 |
| 0.10 | lineage_overall | fraction of all sample x lineage-population items within the same tolerance | 1.0 | <= 0.5 |
| 0.10 | lymphocyte_counts | fraction of samples whose `lymphocytes` count is within 3% of the planted count | 1.0 | <= 0.5 |
| 0.25 | activation_every_sample | per sample, fraction of the two FMO-gated populations within max(0.75 points, 6% of the reference); **minimum** over samples | 1.0 | <= 0.5 |
| 0.10 | activation_overall | fraction of all sample x FMO-gated-population items within the same tolerance | 1.0 | 0 |

**Why these tolerances.** Lineage populations and matrix: as at level 3 (two unrelated correct
strategies agree with the planted values to a small fraction of the band). FMO-gated populations: the
band is tied to the SOP's own acceptance criterion. In a sweep of deliberately mis-set coefficients on
eight seeds, a compensation whose FMO medians agreed within the stated 50 units stayed within 0.58 of
the band (25 block results), between 50 and 100 units one of 47 block results left the band, and each
careless matrix misses the criterion by a factor 16 or more (the residual signal widens the FMO
background of CD8 T cells and loses a large part of the dim positives). Measured values are in
docs/families/flow_compensation.md ("Level 4").

**What is deliberately not scored.** Counts of populations other than `lymphocytes`, the SOP's
acceptance check itself (it is the analyst's own control; only its consequences are scored), any
plot, any additional file.

**Level-7 decision point (acquisition disturbance).** In two samples a contiguous time window of the file
(12-25% of its events) was recorded with the two detectors of the 561-nm laser reading 0.30-0.45 of the
signal. The reference is that of the same simulated session without the disturbance: lineage membership is
the hidden label of every event, disturbed or not, and the FMO-gated populations are the SOP rule applied to
the signals as they would have been recorded. The scorer is the level-4 scorer, unchanged; the reference of
a level-7 task is identical to that of the level-4 task of the same seed.

Three treatments reach the reference and are probes of the battery: the window is set aside and the
populations that need the affected detectors are estimated from the other events; everything from the start
of the window to the end of the file is set aside; the affected detectors are rescaled inside the window.
The generator only plants a window if the reference percentages of the events that remain, in both of the
first two cases, are within half the tolerance band for all 14 populations, so the estimate from the
undisturbed part is within tolerance by construction and not by sampling luck. Measured on twelve seeds
(reference pipeline, worst item as a fraction of its band): window set aside 0.48 lineage / 0.41 FMO-gated,
everything after the start set aside 0.64 / 0.57, rescaled 0.02 / 0.37.

Without the step (the level-4 reference pipeline on these inputs) six of the 84 scored percentages leave the
band, the same six on twelve seeds. In each of the two samples CD4_T is 8-14 points low (4.4-8.0 bands) and
DP_T 8-14 points high (15.8-28.8 bands; reported 9-17% instead of 1-3%), because the CD4 T cells of the
window are counted as CD4+CD8+, and activated_CD8_T is 1.0-5.8 points low (1.3-4.1 bands). The lymphocyte
counts are unaffected. Score 0.472 on every seed; the lineage effect alone would give 0.739, the effect on
activated_CD8_T alone 0.733.

**Scored count.** The only scored count is that of `lymphocytes`. The population is defined by scatter, the
viability dye and CD14, none of which is recorded by the disturbed detectors, so it can be counted on all
events of the file (both solutions report it within 0.3% of the planted count). A submission that removes
the window and reports the count of the events it kept is 19-26% low in two samples, loses two thirds of
the `lymphocyte_counts` component and passes with 0.933 if everything else is right (probe
`time_gate_counts_of_kept_events`).

## Verifier attack suite (must score low)

- `l6_reference` (max 0.75): the reference pipeline of levels 4-6, unchanged, on the level-7 inputs: every decision point of the lower levels is handled, the files are not checked along Time and all events are gated
- `l6_second_solution` (max 0.75): the independent second solution of levels 4-6, unchanged, on the level-7 inputs (mixture-model gates instead of density valleys, no check along Time)
- `exclude_wrong_window` (max 0.75): T_acquisition_disturbance mis-handled: the level-7 reference pipeline, but the time gate (right length) is set on another part of the file: the start of the run if there is more room before the disturbed stretch than after it, otherwise the end
- `l3_reference` (max 0.75): the level-3 reference method on the level-7 inputs: one matrix for all samples from the first-listed single-stain control of each fluorochrome, level-3 gates, HLA-DR + by the density-valley rule used for every other marker (no FMO gate, no acceptance check)
- `l3_reference_last_control` (max 0.75): as l3_reference, but with the last-listed single-stain control of each fluorochrome (what the unmodified level-3 script does when a fluorochrome is listed twice)
- `careless_embedded_matrix` (max 0.75): P_embedded_matrix ignored: samples and FMO tubes compensated with the $SPILLOVER keyword of the files, which is also reported; FMO gate applied, acceptance check not acted on (signal-versus-time check done, disturbed stretch set aside)
- `careless_one_matrix` (max 0.75): P_voltage_change ignored: the session-start matrix (correct controls) and the session-start FMO tube are used for every sample; acceptance check not acted on (signal-versus-time check done, disturbed stretch set aside)
- `careless_row_only` (max 0.75): P_voltage_change half handled: for the later samples only the PE-Cy7 row is replaced from the re-recorded bead tube (and the re-recorded FMO is used); the spill of the other fluorochromes into the HLA-DR detector is left at its session-start value; acceptance check not acted on (signal-versus-time check done, disturbed stretch set aside)
- `careless_bead_tandem` (max 0.75): P_tandem_lot ignored: PerCP-Cy5.5 row from the bead control stained with the other lot (detector settings handled correctly); acceptance check not acted on (signal-versus-time check done, disturbed stretch set aside)
- `careless_no_fmo` (max 0.75): Q_fmo_gate ignored: correct matrices, HLA-DR threshold = 99th percentile of the unstained control instead of the FMO control (signal-versus-time check done, disturbed stretch set aside)
- `careless_fmo_uncompensated` (max 0.75): Q_fmo_gate mis-applied: correct matrices, FMO percentile and sample gate taken on the uncompensated HLA-DR-detector signal (lineage gates on compensated data) (signal-versus-time check done, disturbed stretch set aside)
- `textbook_prior` (max 0.1): identity matrix and textbook PBMC percentages, no data analysed
- `block1_matrix_true_frequencies` (max 0.85): true frequencies, planted block-1 matrix reported for every sample
- `all_percent_beyond_tolerance` (max 0.45): true matrices and counts, every percentage 1.4 x tolerance off
- `one_population_beyond_tolerance` (max 0.89): true matrices and counts, one lineage population per sample 1.5 x tolerance off
- `activation_beyond_tolerance` (max 0.7): true matrices and lineage, both FMO-gated populations of every sample 1.5 x tolerance off
- `one_activation_item_beyond_tolerance` (max 0.8): true matrices and lineage, activated_CD8_T of one sample 1.5 x tolerance off
- `matrix_entry_beyond_tolerance` (max 0.85): true frequencies, one spillover coefficient of one sample off by 0.06
- `naive_matrix_not_applied` (max 0.45): matrices computed and reported, samples and FMO gated on the uncompensated signals (signal-versus-time check done, disturbed stretch set aside)
- `naive_no_doublet_gate` (max 0.45): correct compensation, no aggregate exclusion (singlets = cells) (signal-versus-time check done, disturbed stretch set aside)
- `naive_no_viability_gate` (max 0.45): correct compensation and singlets, dead cells not excluded (signal-versus-time check done, disturbed stretch set aside)
- `naive_global_thresholds` (max 0.6): correct compensation and gates, but one gate template set on the first sample and applied unchanged to all samples (signal-versus-time check done, disturbed stretch set aside)
- `naive_default_pipeline` (max 0.2): embedded matrix, no doublet gate, no viability gate, one gate template, no FMO (signal-versus-time check done, disturbed stretch set aside)

## Must-pass probes (guard against false negatives)

- `planted_truth` (min 0.999): submission written directly from the planted matrices, the event labels and the SOP rule applied to them
- `within_tolerance` (min 0.99): matrix entries +-0.008, every percentage moved by 0.6 x its tolerance, counts +-1.5%
- `independent_solution` (min 0.99): second implementation: disturbed stretch located by an event-level rank scan (no bins, no medians), cut out of the tube, counts scaled back to the file; fcsparser reader, regression spillover, gain factor from mixture-model bead peaks, mixture-model gates
- `time_gate_everything_after` (min 0.99): reference pipeline, but everything from the start of the disturbed stretch to the end of the file is set aside (35-60% of the events remain for the populations that need the affected detectors)
- `repair_instead_of_exclude` (min 0.99): reference pipeline, but the stretch is not set aside: inside it the affected detectors are divided by their factor (position of the peaks of the log signal inside against outside), then all events are gated
- `time_gate_counts_of_kept_events` (min 0.9): reference pipeline, but the disturbed stretch is removed from the file and the sizes of what is left are reported. Floor = pass threshold instead of 0.99: the percentages are right, the lymphocyte count of the two samples is low by the share of the window, which costs two thirds of the lymphocyte_counts component (0.933)
