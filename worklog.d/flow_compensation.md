# Worklog - family `flow_compensation`

All entries 2026-11-11 (single session). Times are approximate order, not wall-clock stamps.

## 1. Reading and design choices
- Read AUTHORING.md, the spec, the crispr_screen reference family and the core modules.
- Panel fixed so that the three large spillovers land on gated markers: CD3-BV421 -> viability
  detector, CD4-PE -> CD8 detector, CD56-APC -> CD19 detector. First draft had spill x source
  brightness well below the true positives of the receiving channel; a per-sample density-valley
  gate would then still have separated "spill" from "true positive" without compensation, i.e. the
  trap would not have bitten. Brightness was therefore set so that the spill lands at the level of
  the true positives (CD3 40 000 x 0.18-0.32 vs dead 12 000; CD4 30 000 x 0.18-0.32 vs CD8 8 000;
  CD56 20 000 x 0.20-0.38 vs CD19 5 000).
- `per_sample_shift`: decided against a fluorescence PMT change. It rescales every spillover
  coefficient into and out of that detector, so the controls recorded at the panel voltages no
  longer describe that sample; a correct answer would need a second control set or an estimate of
  the gain, and I could not guarantee recoverability in the time available. Implemented as a
  lowered FSC voltage (affects FSC-A and FSC-H equally, compensation untouched).
- FCS: wrote a small FCS 3.1 writer inside the family instead of adding a host dependency, so the
  generator needs numpy only. No package was added to the pipeline venv.
- flowkit and flowutils are not on conda-forge; fcsparser 0.2.8 pins numpy<2 while FlowKit 1.3.2
  needs numpy>=2. Image: conda numpy/scipy/pandas/scikit-learn/flowio/six/pip (pinned), pip
  flowkit==1.3.2, then `pip install --no-deps fcsparser==0.2.8`.

## 2. First host test (seed 1, level 3) - independent solution failed
- Oracle 1.000. Independent solution 0.517: in one sample `live` = 100% and NK = 0%.
- Cause: the 2-component Gaussian mixture on arcsinh(x/150) converged to "narrow core + broad
  component" for the viability channel (means 1.16 / 1.74, sd 0.52 / 2.50) because the compensated
  negative population is wide and non-Gaussian on that scale. With the planted matrix the same fit
  was fine; with the regression-estimated matrix it was not, so the fit was fragile, not wrong by
  luck.
- Change: data-driven arcsinh cofactor (2 x |1st percentile| of the parent population, >= 150)
  and a fully specified start (means at the 5th / 99.9th percentile, equal start variances,
  weights 0.8/0.2). After that 1.000 on seeds 1-4.

## 3. Level 1 - oracle failed
- Oracle 0.25 at level 1 (independent solution 1.000): `live` = 17% in four of six samples.
- Cause: with only 0.4-1.2% dead cells the valley finder took the dip at zero that arcsinh with
  cofactor 150 creates inside a spread negative population (T cells spread by the CD3 spill) as the
  boundary. At level 3 the dead-cell mode was prominent enough to hide the problem.
- Change: cofactor = 3 x rms of the below-zero values of the parent population (>= 150); peak
  prominence floor lowered from 0.4% to 0.1% of the main peak so that a 0.4% dead population is
  still found. Oracle 1.000 on level 1 seeds 1-2 and level 3 seeds 1-6.
- Note: the two solutions now both use a data-driven cofactor, but estimate it differently.

## 4. Naive variants - two traps were too weak, retuned
- First measurement (tolerance 0.75 points / 5%):
  - no viability gate 0.50: almost all dead cells were CD14-"positive" and went to the monocyte
    gate; lymphocyte count only +1.1%, T-cell subsets unchanged. That contradicts the spec ("lineage
    frequencies shift").
  - one gate template 0.53 (seed 1) and 0.67 (seed 2): with FSC gain 0.38-0.48 the copied debris
    gate removed the dead cells of the shifted sample but kept its lymphocytes; only `cells` and
    `live` were off. 0.67 is under the 0.75 ceiling but with little margin.
  - no doublet gate 0.35-0.37, uncompensated 0.25 / 0.00: fine.
- Changes to the simulator (not to the gates or ceilings): non-specific binding of dead cells
  lowered from 0.12 to 0.06 of the positive level with a wider per-cell spread, so they smear
  across the thresholds instead of all being positive; FSC gain range changed to 0.20-0.27 so that
  the lymphocytes of the shifted sample fall below a debris gate set on any other sample.
- Tolerance tightened to max(0.5 points, 3%) and 3% for the lymphocyte count: correct solutions
  were within 0.03 bands either way, and the wrong side gained margin (the `singlets` item of a
  6%-aggregate sample was only 1.3 bands off at 5%).
- After the change (level-3 seeds 1-6, host): oracle and independent solution 1.000; no viability
  gate 0.250; one gate template 0.483; no doublet gate 0.272-0.328; matrix not applied 0.250;
  embedded matrix 0.000; all four ignored 0.000.
- Ceilings were then set 0.15-0.2 above these measurements (0.40, 0.20, 0.50, 0.45, 0.65, 0.20).

## 5. Reader validation
- All 15 files of seed 1 level 3: FlowKit `Sample` (version 3.1, event counts), flowio and
  fcsparser return identical arrays; FlowKit differs only in the Time channel, which it multiplies
  by `$TIMESTEP`; FlowKit `apply_compensation` with the planted matrix equals `X @ inv(S)`.
  Repeated inside the agent image (python 3.12.16, numpy 2.3.5, FlowKit 1.3.2, flowio 1.4.0,
  fcsparser 0.2.8).
- A first comparison asserted FlowKit raw events == flowio raw events and failed; that was my
  assertion (Time scaling), not the files.

## 6. Docker validation
- `make --seed 1 --level 3 --keep-work`: all gates ok including Harbor oracle 1.0 / nop 0.0 (3.64 min).
- `make --seed 2 --level 3 --no-harbor`: ok (1.15 min). `make --seed 1 --level 1 --no-harbor`: ok (1.19 min).
- Extra: `make --seed 1 --level 2 --no-harbor`: ok (0.85 min).
- Seed 2 L3 and seed 1 L1 were run a second time with `--keep-work` to read the raw metrics from
  `details.json`; same gate results.
- Docker scores equal the host scores for every variant.

## 7. Wording fix after validation
- Trap text said the doublets "double to triple" DP_T; measured 1.9x-3.9x. Changed the wording in
  `__init__.py` and regenerated the four task directories with `--no-validate`. Checked with a
  checksum over every file except README.md that the regenerated directories are byte-identical to
  the validated ones (also confirms that generation is deterministic).

## 8. Housekeeping (coordinator request: small disk footprint)
- No tools image was built. Removed the images of each instance once the next one was validated,
  and all `alegen-flow-compensation-*` images at the end.
- Deleted `results/gen/flow_compensation-*` (regenerated in < 1 s per instance) and the copied
  inputs / Harbor task copies under `results/work/flow-compensation-*`; the `details.json` and
  output files of every probe and attack are kept there.
- The development environment used for host tests lived in the session scratch directory and was
  deleted.

## Not done / open
- No fluorescence-detector voltage change (see 1).
- The case "gate template taken from the shifted sample" was not measured.
- "Not subtracting bead background costs <= 0.004" is a calculation, not a measurement.
- Harbor (G6) was only run for seed 1 level 3.
- No real agent was run on the task, so the instruction has not been tested for ambiguity by
  anything other than my own two solutions.

---

# Level 4 (2026-11-12, single session, 00:47-02:10 UTC)

Scratch: `/tmp/claude-0/-home-user/5a931900-2990-511a-8e45-6939bdda2f59/scratchpad/flow_compensation_l4/`
(host dev environment with the image's package versions, host replay of the battery, sweeps).

## L4-1. Reading, baseline, design (00:47-00:52)
- Read AUTHORING.md "Level 4", the level-4 spec, the level-3 notes/limitations, README section 6.
- Before touching the code: generated levels 1, 2, 3 (seed 1) and level 3 (seed 2) into the scratch
  directory as a byte baseline; `diff -r` of level 3 seed 1 against `tasks/flow-compensation-s0001-l3`
  was empty.
- Design decision: everything new goes through one detector, YG780-A (HLA-DR PE-Cy7), which level 3
  recorded but did not score. HLA-DR on T cells becomes a dim, valley-free signal gated on an FMO
  tube (Q); the voltage of that detector is raised mid-run (P); the PerCP-Cy5.5 bead control comes
  from another lot and spills less into that detector (P); the files carry an old template matrix
  (P). Reason: one coherent story, one in-data control (the FMO tube) that exposes all three
  matrix defaults.
- Level 4 lives in a new module `l4.py` and new assets (`flow_pipeline_l4.py`, `solve_alt_l4.py`,
  `score_l4.py`, `solve_l4.sh`); `__init__.py` only got a three-line dispatch and a docstring line,
  the level 1-3 assets were not touched. This is what keeps levels 1-3 byte-identical.
- Matrix deliverable: the spec suggested one file per acquisition block; I chose one matrix per
  sample in a long table, because a per-block file tells the agent how many blocks there are.

## L4-2. First implementation and host run (00:52-01:00)
- Both solutions scored 1.000 on the first host run of seed 1; all careless variants <= 0.65.
- `one_population_beyond_tolerance` scored exactly 0.900 = pass with my first level-4 weights (one
  lineage item out per sample cost only 0.10). Level 3 fails such a submission (0.833). Changed the
  ramp of `lineage_every_sample` to zero at 10/12: now 0.858.
- Six seeds on the host: bead-lot variant left two of six `activated_CD8_T` items in band on seed 1
  (lot difference 0.082 at the low end of my 0.08-0.12 range). Raised the planted difference to
  0.10-0.14; afterwards every CD8 item of the bead-lot variant was 13-74% off on eight seeds.
- The reference pipeline's largest matrix error (0.0065, seed 4) was the PerCP-Cy5.5 -> YG780-A entry
  of block 2: sampling noise of the cell-based single stain (16 000 events) times the gain factor.
  Raised that tube to 30 000 events (input +0.7 MB): <= 0.0041.

## L4-3. Docker validation of intermediate versions (01:06-01:28)
- Seed 1 passed all gates including Harbor three times (3.54, 3.86 and 5.07 min) with three
  intermediate versions (acceptance criterion 150 units on three, then five populations). None of
  these is the shipped version; the logs are kept in the scratch directory.
- Seed 2 was not built in the second and third round: my own guard refused (8149 MiB and 7658 MiB
  free, limit 8192 MiB). Disk was taken by other authors' builds; I removed my stale seed-1 images
  before every rebuild and deleted scratch copies I no longer needed.
- FlowKit 1.3.2 inside the agent image reads all 19 files, parses the embedded matrix and applies it.

## L4-4. Notes wording (01:10)
- First draft of the notes put both PerCP-Cy5.5 lots side by side in one bullet. That is the kind of
  hand-over sentence the pilot showed to be sufficient for a small model. Now the notes say when and
  from which vials the beads and the extra cell tube were stained; the lots are in `panel.tsv` and
  `controls.tsv`. Still facts only.
- Added `l3_reference_last_control`: the unmodified level-3 script overwrites a fluorochrome's row
  with the last-listed control. Scores 0.45 like `l3_reference`.

## L4-5. Acceptance criterion: what failed and what changed (01:17-01:50)
- First version: medians of CD4_T, CD8_T, DN_T in the FMO tube within 150 units. Two problems.
- (a) An analyst who tunes only the two coefficients that the three T-cell populations see would pass
  the check while the APC-Cy7 -> YG780-A entry of block 2 is still at its session-start value
  ((gain - 1) x coefficient; calculated, not measured). Added B_cells and NK_cells to the criterion,
  so that every coefficient of the column that a lymphocyte population carries is covered.
- (b) Sweep without fine-tuning (seeds 1, 2): an error of 0.0025 in PE -> YG780-A gave a spread of 78
  units in block 1 (criterion met) and moved `activated_CD4_T` of one sample by 0.76 points (1.5 bands
  of 0.5 points): score 0.742. The criterion was not sufficient for the band, which the Q element
  requires. Cause: PE is bright (30 000), the tubes differ by +-8% in staining, and in a sample with
  4% activated CD4 T cells a quarter of the reported events are the 1% tail of the negatives.
- Changes: criterion 50 units; absolute band of the FMO-gated populations back to 0.75 points (I had
  tightened it to 0.5 in between); sentence "has to be corrected before the HLA-DR populations are
  reported" instead of "are not reported" (the old sentence could be read as "report nothing").
- Re-measured on eight seeds, seven mis-settings each, no fine-tuning (112 block results): spread
  <= 50: worst item 0.58 band, none out (n = 25); 50-100: one of 47 out (1.13 band); 100-150: none of
  24 out (0.96); > 150: one of 16 out (1.56).
- With a 50-unit criterion the bead-based matrices of both solutions pass on eight seeds (spread
  2.7-35.8) but without much margin, so both solutions now fine-tune the YG780-A column on the FMO
  tube.

## L4-6. Fine-tuning: first method was wrong (01:35-01:45)
- First attempt in the reference pipeline: event-level least squares of the compensated YG780-A
  signal of the FMO lymphocytes on their other compensated signals. It made the result worse: spread
  13 -> 20 and 36 -> 43 (seed 1), largest matrix error 0.004 -> 0.007.
- Cause, checked with the planted matrix: the fit returns slopes of -0.0020 (PerCP-Cy5.5) and -0.0021
  (APC-Cy7) although the matrix is exact. The compensated signal contains minus coefficient x noise
  of the regressor, so the slope is biased by -coefficient x Var(noise) / Var(signal). Population
  medians of the same data agree within 4 units.
- Replaced by a weighted least-squares fit across seven populations on population medians
  (reference pipeline) and by per-population median matching against DN_T (second solution).
  Eight seeds: spread after fine-tuning 0.2-12.2 (reference), 0.0-0.2 (second solution); largest
  matrix error 0.0017-0.0022 and 0.0021-0.0025; FMO-gated items within 0.12 and 0.10 of the band.
- Consequence I had not planned: the check plus fine-tuning repairs the bead-lot row and the
  unscaled column without any reasoning about lots or gain (0.959-1.000 on eight seeds, see the
  design notes). The careless variants of the battery therefore do not act on the check.

## L4-7. Final validation and wrap-up (01:53-02:10)
- `make --seed 1 --level 4 --keep-work`: all 34 checks ok, Harbor oracle 1.0 / nop 0.0, 6.75 min
  (host load average about 13). `make --seed 2 --level 4 --no-harbor --keep-work`: all 32 checks ok,
  3.01 min. Logs `results/flow_s1_l4.log`, `results/flow_s2_l4.log`.
- Host replay of both batteries with the shipped code: 24 submissions per seed, every score equal
  to the Docker score.
- After the validation I found that the judge text quoted 0.61 for "spread <= 50"; 0.61 belongs to a
  spread of 51.7, the value for <= 50 is 0.58. I also found that the evidence line of
  `P_embedded_matrix` said the acceptance check "fails" with the embedded matrix, while the measured
  fact is that it cannot be evaluated (a population is missing). Corrected both sentences in `l4.py`,
  regenerated both task directories with `--no-validate` and compared them with copies of the
  validated directories (`diff -rq`): only `README.md` differs. I did not re-run the gates for this.
- Level 3 seed 1 regenerated into a temporary directory at the end: `diff -r` against
  `tasks/flow-compensation-s0001-l3` empty.
- Housekeeping: seed-2 images removed after seed 2 passed; the seed-1 images
  (`alegen-flow-compensation-s0001-l4-agent:local` 1.55 GB, `-verifier:local` 178 MB) are kept for
  the pilot; `results/gen/flow_compensation-s000{1,2}-l4`, the Harbor task copy
  `results/work/flow-compensation-s0001-l4/harbor_task` and the copied inputs of the `inputs_copied`
  attack removed (regenerated in seconds), `details.json` of every check and the outputs of every
  other probe and attack kept; scratch directory removed.
  `tasks/flow-compensation-s0002-l4` (20 MB) is left in place next to the seed-1 task.

## Level 4: not done / open
- The rare double-positive population of the spec was not built.
- Harbor (G6) only for seed 1. No agent was run.
- Not measured: an analyst who takes the gain from the voltage ratio with a textbook exponent; an
  analyst who fine-tunes on the FMO tubes but ignores the re-recorded PE-Cy7 bead tube (the PE-Cy7
  row of the later samples would then be off by (1 - 1/gain) x coefficient, calculated 0.03-0.09).
- "Tuning only the T-cell coefficients leaves APC-Cy7 -> YG780-A wrong" (L4-5a) is a calculation.
- The two solutions fine-tune with related methods (population medians).

---

# Level 7 (2026-11-12, 16:10-17:45 UTC)

Spec `docs/specs/l7/flow_compensation.md`: one new silent decision point (acquisition disturbance that
only the Time parameter shows), level-4 texts unchanged, text removal by `descaffold.apply(build, 7)`.

## L7-1. Reading and baseline (16:10-16:16)
- Before touching anything: regenerated seed 1 of levels 3, 4, 5, 6 with `--no-validate` into the
  scratch directory and compared with `tasks/`. Levels 3 and 4: `diff -r` empty. Levels 5 and 6:
  `README.md` and `task.toml` differ (title and summary wording, a `source_level` line), everything
  else identical. Cause: `alegen/descaffold.py` was changed (commits f09a54d, 671d503) after those two
  task directories were written (6fae8fb, b1d5918). Not a family matter and not mine to fix; noted so
  that the later comparison is read correctly.
- The host venv has no scipy / pandas / flowio, so the pipelines cannot be run on the host. Used a
  private container of the level-4 agent image (`flowl7-dev`, `--network none --cpus 2`, removed at the
  end) with the scripts called with `--input/--output`, and the task's `score.py` on the host.

## L7-2. Design (16:16-16:22)
- Disturbance on the two detectors of the 561-nm laser (YG586-A: CD4 PE, YG780-A: HLA-DR PE-Cy7), as
  the spec suggests. Worked out beforehand what an all-events analysis would report: with a factor of
  0.30-0.45 CD4 stays positive (10 000 against a negative limit of 540), so the lineage error cannot
  come from the CD4 gate; it comes from the compensation (PE spill left in B710-A). Predicted items
  out of band per sample: CD4_T, DP_T, activated_CD8_T; predicted score 0.472. Measured later: exactly
  those six items and 0.472 on all sixteen seeds tried.
- Truth: the level-4 reference of the same seed. `generate_l4` got an optional `disturb` function that
  returns what is written to a sample file; the reference is computed from the simulated recording as
  before. With `disturb=None` the function is unchanged.
- The disturbance has its own random stream, so the level-7 task of a seed is the level-4 session of
  that seed with two files changed.

## L7-3. Where the stretch can be (16:22-16:30): first version failed
- First version: one pass, the two samples drawn at random, then the first of 400 random positions
  for which the reference percentages of the remaining events (stretch removed; everything from its
  start removed) are within 0.5 band for all 14 populations. Seed 1 raised: for PBMC_D858 the best of
  400 was 0.79.
- Survey of 2000 positions per sample (seed 1): "stretch removed" is nearly always fine (73-100% of
  the positions within 0.5 band). "Everything after the start removed" is the problem: the kept set
  is a prefix of the file, so all candidate positions share most of their events, and a sample whose
  first 25-60% happen to be short of NK cells stays short whatever the start (within 0.5 band in
  0-40% of the positions, within the band in 40-91%).
- Change: two passes over the same simulated session. Pass 1 looks for a position in every sample;
  per block one sample is drawn among those that have one; pass 2 writes the files. If a block has no
  such sample the bound goes to 0.6, then 0.7; beyond that the generator raises.
- Scan of seeds 1-48: 90 stretches at 0.5, two at 0.6, two at 0.7; seed 32 raises (best 0.75). Side
  effect I had not planned: 38 of 94 stretches start after 55% of the run (a uniform draw would give
  14%).

## L7-4. The statistic that finds the stretch (16:24-16:32)
- First criterion in the generator's own check: bins of 500 events, a bin is flagged if its median
  is more than 6 robust SD from the median of all bins. It broke the stretch into pieces: in seed 1 /
  PBMC_D735 only 5 of 12 inside bins were flagged for YG780-A although the level drops from 1255 to
  447, and in seed 2 / PBMC_D124 none for YG586-A.
- Looked at the series. The median of YG586-A lies between populations and moves between 770 and
  2150 from bin to bin without any disturbance; the median of V450-A jumps between 1200 and 22 800 (it
  sits in the gap below the T cells). The median of YG780-A is stable (1500-1750), and so are the 75th
  and 90th percentiles of YG586-A.
- Reference pipeline now: median and 90th percentile per bin on an arcsinh scale, 4 robust SD, a run
  of at least 5 bins with at most 10% flagged bins elsewhere, then a confirmation (factor inside /
  outside over the 50th-95th percentiles differs from 1 by more than 15% in at least one parameter).
  No false report in 72 undisturbed tubes of twelve seeds.
- Generator check now: every inside bin and no outside bin nearer to the inside level, in at least
  one of the two detectors (YG780-A: 94 of 94 stretches; YG586-A: 35 of 94).
- FlowKit 1.3.2 has no `Sample.filter_anomalous_events` (AttributeError); the image offers
  `plot_channel` only.

## L7-5. Second solution: first version was wrong (16:27-16:30)
- Rank scan per parameter (cumulative sum of centred ranks, stretch from the highest to the lowest
  point), taken for both signs, union of all significant stretches. Scores 0.614 (seed 1) and 0.872
  (seed 2): 70-78% of the events of the disturbed files were cut. A stretch of low values makes the
  events before it look high, and that "high" stretch was significant too (z = 16-25 against 31-48
  for the real one).
- Fix: per parameter only the reading with the larger z. For a stretch of share w that starts after
  a share b of the run, z(prefix) / z(stretch) = sqrt(w b / ((1 - w)(1 - b))) < 1 as long as
  w + b < 1, which the generator guarantees (stretch ends before 85%). After the fix 1.000 on twelve
  seeds; z = 27.6-51.7 in disturbed tubes, at most 4.8 in undisturbed ones; threshold 10.

## L7-6. The repair treatment: three things had to change (16:33-16:45)
- First run on twelve seeds (jitter 2-4%, factor = median of percentile ratios, edges from a
  least-squares change point on the summed arcsinh signal): 11 pass, seed 4 scores 0.742
  (`activated_CD4_T` 10.62 for 9.26, 1.8 bands), seeds 10 and 12 at 0.79 and 0.87 band.
- Diagnosis on seed 4 / PBMC_D665 against the undisturbed recording. With the planted factor and the
  planted edges the repaired file still gives 12.8% instead of 9.4% HLA-DR + CD4 T cells inside the
  stretch: the jitter of the repaired PE signal (3.5%) leaves +-250 units in B710-A after
  compensation, and the PerCP-Cy5.5 -> YG780-A coefficient of block 2 (0.456) carries that into the
  HLA-DR detector. With the estimated factors (0.397 / 0.403 for 0.406) 15.8%.
- Changes: (1) jitter 1-2% instead of 2-4% (generator); (2) factor from the shift that brings the
  log-signal densities inside and outside into register, instead of percentile ratios (errors 5% ->
  at most 2.2%); (3) edges from an event-wise log ratio of the two histograms instead of the mean
  (seed 12 / PBMC_D801: the old change point was 178 events late because the 180 clean events after
  the stretch happened to hold few CD4 T cells; now within 7 events on 24 stretches, 21 of them
  within 2).
- After: repair 1.000 on twelve seeds, worst item 0.02 band (lineage) and 0.37 band (FMO-gated).
- The exclusion treatments never depended on any of this (guard of 250 events).

## L7-7. Battery and host replay (16:45-17:10)
- Battery: level-4 battery with the "ignores one element" scripts moved to the level-7 pipeline
  (time check done), the lower-level reference methods left as they are, plus `l6_reference`,
  `l6_second_solution`, `exclude_wrong_window` and three probes. Reasons in the design notes.
- Host replay: seeds 1-12, 30 checks each: 360 checks, none outside its bound. Seeds 13, 20 (bound
  0.7) and 33, 48 (bound 0.6), level-7 checks: all inside; `time_gate_everything_after` worst item
  0.79, 0.78, 0.58, 0.55 band.
- A background job started with `&` inside a background command was killed with its shell; noticed
  from an empty log, restarted as a plain background command.

## L7-8. Levels 1-6 unchanged (17:08 and again 17:36 with the final code)
```bash
S=<scratch>/flow_l7_check
git archive HEAD pipeline/alegen | tar -x -C $S/head_src          # family code before the change
(cd $S/head_src/pipeline && for L in 1 2 3 4 5 6; do /home/user/ale-lifesci-gym/.venv/bin/python -m alegen \
    make flow_compensation --seed 1 --level $L --no-validate --out $S/head --results $S/res_head; done)
cd pipeline && for L in 1 2 3 4 5 6; do ../.venv/bin/python -m alegen make flow_compensation \
    --seed 1 --level $L --no-validate --out $S/after --results $S/res_after; done
for L in 1 2 3 4 5 6; do diff -r $S/head/flow-compensation-s0001-l$L $S/after/flow-compensation-s0001-l$L; done
for L in 3 4 5 6;     do diff -r $S/after/flow-compensation-s0001-l$L ../tasks/flow-compensation-s0001-l$L; done
```
- Code before the change against final code, levels 1-6: all six `diff -r` empty (rc 0).
- Final code against `tasks/`: levels 3 and 4 empty; levels 5 and 6 `README.md` and `task.toml` only
  (12 and 16 diff lines), byte for byte the same diff as in L7-1 before any edit.

## L7-9. Docker validation (17:10-17:35)
- Disk checked before each run (4.5-6.2 GB free); `docker builder prune -af` after each.
- Round 1: seed 1 with Harbor, 41 checks ok, 4.08 min; seed 2 without Harbor, 39 checks ok, 2.38
  min. Every score equal to the host replay.
- Then found that the judge text gave the deviations of the naive analysis as "22-29 bands" for
  CD4_T and DP_T and "2.7-4.1 bands" for activated_CD8_T; those were the worst item per seed, not the
  range per item (CD4_T 4.4-8.0, DP_T 15.8-28.8, activated_CD8_T 1.3-4.1). Corrected in `l7.py`,
  validated both seeds again. `diff -rq` of the task directories of the two rounds: `README.md` only.
- Final round, seed 1 (`--keep-work`, with Harbor): G1 images_build ok (agent image derived from
  alegen-flow-compensation-s0001-l4-agent:local); G7 ok; G2 oracle 1.000 (0.04 min); G5 ok; G3
  planted_truth 1.000 (>=0.999), within_tolerance 1.000, independent_solution 1.000,
  time_gate_everything_after 1.000, repair_instead_of_exclude 1.000 (>=0.99),
  time_gate_counts_of_kept_events 0.933 (>=0.9); G4 has_l3_reference_attack ok,
  has_l6_reference_attack ok, empty_output 0, inputs_copied 0, zero_byte_files 0, header_only 0,
  l6_reference 0.472, l6_second_solution 0.472, exclude_wrong_window 0.472, l3_reference 0.189,
  l3_reference_last_control 0.189, careless_embedded_matrix 0.000, careless_one_matrix 0.517,
  careless_row_only 0.533, careless_bead_tandem 0.500, careless_no_fmo 0.650,
  careless_fmo_uncompensated 0.650 (all <=0.75), textbook_prior 0.000 (<=0.1),
  block1_matrix_true_frequencies 0.800 (<=0.85), all_percent_beyond_tolerance 0.300 (<=0.45),
  one_population_beyond_tolerance 0.858 (<=0.89), activation_beyond_tolerance 0.650 (<=0.7),
  one_activation_item_beyond_tolerance 0.742 (<=0.8), matrix_entry_beyond_tolerance 0.800 (<=0.85),
  naive_matrix_not_applied 0.200, naive_no_doublet_gate 0.281, naive_no_viability_gate 0.267 (<=0.45),
  naive_global_thresholds 0.417 (<=0.6), naive_default_pipeline 0.008 (<=0.2); G6 harbor_oracle [1.0],
  harbor_nop [0.0]. `{"task": "flow-compensation-s0001-l7", "ok": true, "minutes": 4.54}`
- Final round, seed 2 (`--keep-work --no-harbor`): same checks without G6, all ok; scores as seed 1
  except careless_embedded_matrix 0.175, careless_one_matrix 0.483, careless_row_only 0.483,
  textbook_prior 0.008, naive_no_doublet_gate 0.258, naive_default_pipeline 0.000.
  `{"task": "flow-compensation-s0002-l7", "ok": true, "minutes": 2.79}`

## L7-10. Housekeeping
- Left in place: `tasks/flow-compensation-s0001-l7`, `tasks/flow-compensation-s0002-l7`,
  `results/validation/flow-compensation-s000{1,2}-l7.json`, `results/work/flow-compensation-s000{1,2}-l7`
  (61 MB, details of every check), `results/gen/flow_compensation-s000{1,2}-l7`, and the four images
  `alegen-flow-compensation-s000{1,2}-l7-{agent,verifier}:local` (the agent images are a 20 MB layer
  on top of the level-4 image).
- Removed: my dev container and the development data in the session scratch directory. The
  regenerated levels 1-6 (`flow_l7_check/after`) and the four validation logs (`flow_l7_logs/`) are
  still there. No file outside the family, its design notes and this log was edited; nothing was
  committed.

## Level 7: not done / open
- No agent has been run on a level-7 task. Harbor (G6) only for seed 1.
- Seeds 13-48 were generated, not validated (level-7 checks replayed for 13, 20, 33, 48).
- Seed 32 has no level-7 instance (generator raises). Stretch positions cluster late in the run.
- The disturbance model is multiplicative with a 1-2% jitter; a physical noise model would make the
  repair treatment inexact for the FMO-gated populations (design notes).
- Not measured: an analyst who repairs with one common factor for both detectors; an analyst who
  counts B and NK cells on all events; a stretch at the very start or end of a file; two stretches in
  one file (the reference pipeline reports the longest run only).
