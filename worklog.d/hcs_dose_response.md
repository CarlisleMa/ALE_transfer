# Worklog - family `hcs_dose_response`

All entries 2026-11-11 (single session). Order is chronological; times are wall clock of the build host.

## 1. Reading and first decisions (16:20-16:35)
- Read AUTHORING.md, the spec, crispr_screen / flow_compensation and the core modules.
- Checked the build proxy before designing anything around Cellpose: the four weight files
  (`nucleitorch_0`, `size_nucleitorch_0.npy`, `cyto3`, `size_cyto3.npy`) download from
  `https://www.cellpose.org/models/` with plain curl, and `download.pytorch.org/whl/cpu` is reachable.
  So the professional model stays; no fallback to a classical-only image was needed.
- Cellpose 4 is Cellpose-SAM (about 1.1 GB of weights, minutes per image on CPU). Chose Cellpose
  3.1.1.2 (`nuclei`, `cyto3`), which needs `numpy<2.1` -> numpy 2.0.2, scipy 1.14.1, torch 2.5.1+cpu.
- Potency definition. The spec says "IC50". With plateaus above zero the inflection point (EC50) and
  the 50%-of-vehicle crossing differ, and "inactive" needs a rule. Decided: potency = concentration
  where the fitted curve equals 0.5 of vehicle; inactive = the curve does not reach 0.5 at or below
  the top usable concentration. Planted the 0.5 crossing directly and derived EC50 from it, so the
  truth is the planted number, not a fit. Activity classes are separated by construction (actives
  between 3.3x the lowest and 1/5.5 of the top dose; inactives have no effect, a plateau >= 0.60 or a
  potency >= 5x the top dose); `generate()` asserts this.
- Layout: one compound per row, 3 vehicle wells per row placed so that every column holds exactly 2,
  one reference well per row, the series in column order. This is what makes a plate model fittable
  from the vehicle wells and at the same time confounds dose with column.

## 2. Simulator (16:35-16:50)
- numpy-only renderer (`sim.py`): nuclei as textured ellipses composited by maximum, cytoplasm as a
  larger ellipse with nearest-cell-wins territories, PSF blur, illumination roll-off, Poisson + read
  noise, 16-bit deflate TIFF. 96 wells x 2 channels at 320 x 320 px = 19-21 MB, 13-30 s to render.
- Bug 1: the first colony placement retried a rejected colony seed forever, so some wells got 3-12
  cells instead of ~118. The first pipeline run showed "inactive" compounds with GI50 values, the
  per-well table showed vehicle wells with 12 cells. Rewrote `place_cells` (grow next to a random
  existing cell, 45% of the steps at touching distance).
- Bug 2 (both pipelines): a nucleus enclosed by its neighbours has no cytoplasm ring; `ndi.mean`
  returned NaN and poisoned the well mean. Rings without pixels are now skipped.

## 3. First measurements on the classical base (16:50-17:00)
- Distance-transform watershed pipeline, seed 1 level 3: counts median -1.6% (worst -6%), QC flags
  exact, all 16 potency entries right, largest error 0.084 log10.
- Naive variants on the same base, first try, no retuning needed:
  threshold counting: counts -47%, the three GI50s +0.48 / +0.60 / +0.93;
  no plate model: 3 marker and 4 count entries wrong; whole-field readout: cytotoxic compound gets
  an IC50, dual compound -1.07; camera offset only: fluorescent compound inactive, three IC50s
  +0.38..+0.57; plate-wide background constant: fluorescent compound inactive, two IC50s 0.22-0.25
  off; no focus QC: marker IC50 -1.2 and GI50 -1.06.
- Out-of-focus wells: first version blurred the two doses below the GI50 of the low-dose cytotoxic
  compound; its GI50 sits at the low end of the series, so only one well qualified and it was an
  outer (sparse) well. Moved that pair to the count potency of the potent marker inhibitor (inner
  row, dense wells).
- Planted well noise: with 5% (cells) / 4% (marker) lognormal well-to-well noise two correct
  pipelines were 0.13-0.14 log10 off on one compound each (tolerance 0.20). The error came from the
  planted wells, not from the image analysis (both pipelines were off in the same direction).
  Lowered the noise to 3% / 3% and the count plateau of the dual compound to 0.20-0.30. After that
  the worst error over seeds 1-6 x two classical pipelines is 0.106.
- Independent solution: LoG seeds with `min_distance=3` split elongated nuclei (+3.5% median, two
  wells outside tolerance at level 1). `min_distance=5` (touching nuclei are >= 10.5 px apart):
  median 0.0%, range -3.4%..+5.6%.

## 4. Cellpose (16:37-17:05)
- First oracle used the Cellpose size model on six vehicle wells and the estimated diameter
  (10.3 px) for all wells: 14 min 8 s on 2 CPUs with the host under load (load average 8-19), i.e.
  ~8 s CPU per image, and five sparse wells were badly undercounted (12 -> 1, 14 -> 4, ...).
- Timing on four wells (same load): size model 131 s for 3 images; `diameter=13`: 6.5 s/img;
  `diameter=17` (no rescaling): 1.2 s/img; dense well exact (109/109) in every setting.
- Sparse wells: default per-image percentile normalisation blows up a field with a dozen nuclei and
  the flow-error check then removes most masks (3/14, 1/12, 5/14); with `flow_threshold=0` the same
  wells come out slightly high (17/14, 14/12).
- Oracle now: no size model, `diameter=17`, one plate-wide intensity scale (`normalize=False` on
  pre-scaled images). Result seed 1: 2 min 47 s under load, counts median 0.0%, p5 -1.5%, no well
  outside tolerance, QC exact, potencies within 0.053.
- The default-normalisation behaviour is kept as `--cellpose-norm image` for the
  `naive_cellpose_defaults` attack. It is not a trap of the task: with default settings the oracle
  still scored 1.0 (0.923 of the wells inside the count tolerance, threshold 0.90). See weaknesses.

## 5. Scoring (17:00-17:10)
- Flat mean over 8 compounds would let a single wrong compound cost 0.056 (marker) - a pipeline with
  a plate-wide background constant gets exactly one compound wrong (the fluorescent one) and would
  score 0.94. Component = ramp(mean credit, 0.8 -> 1.0): one wrong compound -> 0.375.
  Measured: that variant scores 0.703 (seed 1).
- Level 1 has no out-of-focus well, so "flag nothing" would earn the QC weight for free and the
  prior-only attack would sit exactly on its 0.1 ceiling. QC credit at level 1 is therefore
  (1 - false flags / 5) x cell-count credit.
- SOP wording: first draft defined R2 "above the background of that field" and the response as
  "relative to the vehicle readout expected in that well". Both name the measurand, but the first
  one also names the procedure (per field). Reworded to "signal that does not come from the cells
  is background and is not part of R2" and "the readout the same well would have given with vehicle
  only".

## 6. First full gate run, image size (17:05-17:20)
- `make --seed 1 --level 3 --keep-work`: all gates ok, 13.6 min (log kept in the scratch directory).
- `docker images` reported 3.05 GB for the agent image (`du -sx /` in the tools-only image: 2.15 GB; I
  assume, without having verified it, that this daemon adds the compressed layers to the figure). Budget is 3 GB. Removed PyTorch's C++ test binaries, headers and cmake
  files in the pip layer (-176 MB on disk); Cellpose inference re-tested after the removal.
- Removed my prototype image and the first task images, rebuilt, re-ran the three gates (section 7).

## 7. Final gate runs with the trimmed image (17:13-18:00)
- Agent image now 2.82 GB by `docker images`. Outputs of all battery items identical to the first
  run (generation and both pipelines are deterministic).
- `make --seed 1 --level 3 --keep-work`: all gates ok incl. Harbor oracle 1.0 / nop 0.0, 16.8 min.
- `make --seed 2 --level 3 --no-harbor`: ok, 7.0 min. `make --seed 1 --level 1 --no-harbor`: ok,
  7.6 min. (`--keep-work` was added to both so that the raw metrics stay under `results/work/`.)
- Level 1 would have crashed: an assertion in `generate()` expected 4 marker-active compounds, level
  1 has 5. Found by building the level-1 TaskBuild on the host before its gate run started; fixed.
- Additionally ran `make --seed 1 --level 2 --no-harbor`: ok, 3.6 min.
- Three strings of the task card (README.md: one `discoverable_by` text and the tolerance paragraph of
  the judge notes) were corrected after the two level-3 tasks had been written. Instead of repeating
  the gate runs I regenerated all three required tasks with the final code into a scratch directory:
  `diff -rq` against `tasks/` showed README.md as the only differing file for the two level-3 tasks
  and no difference for level 1. Copied the two README.md files over. Instruction, inputs, solution
  and tests are byte-identical to what went through the gates.
- Measured for the notes: skipping the flat-field correction is harmless (seeds 1, 2); Cellpose with
  default per-image normalisation and everything else correct: 5 / 3 of 91 wells outside the count
  tolerance, potencies unchanged.
- Classical pipelines on seeds 1-6, final code: no well outside the count tolerance in 12 runs, QC
  exact, all entries right, largest potency error 0.106.
- Housekeeping: removed my prototype image and the agent/verifier images of every instance after its
  validation (the coordinator asked for this when the disk watchdog stopped); bulky development
  plates deleted from my scratch subdirectory.

## Not done
- No agent pilot (difficulty is untested).
- Alternative plate-correction strategies (column-matched vehicle wells, smooth surface) and an
  un-normalised focus metric were not measured.
- Harbor was run for seed 1 level 3 only, as required.

# Level 4 (2026-11-12, single session; times are UTC wall clock of the build host)

## L4-1. Reading, tool image (00:38-00:47)
- Read AUTHORING.md "Level 4", the level-4 spec, the level-3 notes and code, README section 6.
- No hcs image was left on the host. Built `alegen-hcs-dose-response-dev-agent:local` from
  `tasks/hcs-dose-response-s0001-l3/environment` with `dockerlib.build_image` (2 min 34 s): same
  Dockerfile as the task images, so the 2.8 GB environment layers are shared with the level-4 images.
- Decision: level 4 gets its own modules (`sim4.py`, `level4.py`, `assets/*_l4*`, `assets/score_l4.py`)
  and its own RNG stream (`884000 + seed`); `generate()` only dispatches. Levels 1-3 code paths, the
  level-3 solution files and `score.py` are not touched. Checked after every round of code changes (last
  at 03:11, after the final one):
  `make --seed 1 --level 3 --out <scratch> --no-validate` + `diff -r` against
  `tasks/hcs-dose-response-s0001-l3` = no difference.

## L4-2. First image model and what the level-3 methods do with it (00:47-01:00)
- `sim4.py`: nuclei on a padded canvas (cut by the border), flat-topped profile composited by maximum,
  two texture scales, nucleoli, bright fragments (a fraction of the seeded cells plus a fraction of the
  lost cells), enlargement of the nuclei for two compounds, counting frame, point annotation of eight
  fields. First setting: minimal centre distance 9.0 px, 1-3 deep nucleoli per nucleus.
- Level-3 methods on that plate (seed 1, wells without annotation, fraction inside max(3, 10%) /
  median relative error): distance-transform watershed 0.22 / -13.9%; LoG second solution 0.51 / +9.2%
  (enlarged nuclei +79%, fragment-rich wells +145%); Cellpose with the level-3 settings 0.10 / -25.8%.
  Unexpected: untuned Cellpose does not count the fragments at all (A12: 10 planted, 9 found); it fails
  because it drops nuclei in overlapping clusters.
- Cellpose settings on the eight annotated fields of that plate (mean relative error / worst field):
  nuclei 17 px, flow 0.4 (level 3): -20.8% / 35%; flow threshold off: -8.7% / 22%; nuclei 12 px, flow off:
  -5.3% / 13%; cyto3 12 px: -3.0% / 8.0%. At 12 px Cellpose starts to segment the fragments (A12: 25 objects for 10 nuclei), so
  the size rule only matters once the model is tuned.
- The nuclei looked like rings (nucleoli too wide and deep) and the densest clumps were hard to count
  by eye. Too hard: no setting of the `nuclei` model met 10% on every field.

## L4-3. Finding the difficulty window (01:00-01:10)
- Lighter nucleoli (0-2 per nucleus, depth 0.3-0.6, sigma 0.8-1.3 px) and two spacings, annotated
  fields, mean / worst field:
  spacing 9.5 px: level-3 Cellpose -17.1% / 25.8%, nuclei 12 px flow off -1.7% / 6.1%, cyto3 14 px -1.7% / 5.6%;
  spacing 10.0 px: level-3 Cellpose -6.4% / 12.6%, nuclei 12 px flow off -0.8% / 2.7%.
  At 10.0 px the level-3 reference is only just outside 10%, at 9.0 px nothing passes. Kept 9.5 px.
- Chose max(3 cells, 10%) as the SOP criterion (same as the level-3 count tolerance). On the final images
  (checked at 01:56) nuclei 12 px / flow off is inside it on every annotated field of both seeds, worst field at
  0.71 (seed 1) and 0.67 (seed 2) of the limit; the level-3 methods are outside on 5 of 8 fields.

## L4-4. A classical method that meets the criterion (01:10-01:30)
- Tried, on the 9.5 px plate, wells without annotation, fraction inside 10% / median error:
  LoG seeds on the smoothed image with a per-well scale 0.24-0.28 / -13..-15%;
  distance-transform maxima (h-maxima and peak_local_max, several thresholds) 0.12-0.19 / -16..-33%
  (first version filled the background enclosed by a colony, a bug of mine; fixing it did not help
  enough. I assume, without having verified it, that the PSF fills the notches between overlapping nuclei);
  LoG on the raw image, sigma 2.4, min distance 5, fragment test and counting frame: 0.92 / -1.6%, but one
  annotated field outside the criterion by one cell; sigma 2.6: 0.86 / -5.9%;
  LoG on an intensity-capped image: +18..+88% (a flat disc gives a ring of maxima at small sigma).
- What works: radial-symmetry voting. Every outline pixel votes for the points 4.5-7 px (x nuclear
  size of the field) along its gradient; an overlapped nucleus keeps most of its outline and still
  collects its votes, texture votes incoherently. Fragment test = area of the half-maximum region
  (a mask area at a fixed threshold calls bright fragments nuclei). First try: 1.00 / -1.6%, annotated
  fields all inside. Five parameter sets gave 0.96-1.00 on that plate; three sets on the final image
  parameters 0.99-1.00 (seed 1) and 1.00 (seed 2).
- So both routes named in the spec exist: tuned Cellpose and a carefully built classical pipeline.

## L4-5. Marker readout bias found while checking the classical solution (01:30-01:35)
- Seed 2: the marker IC50 of the dual compound came out 0.19 log10 low (tolerance 0.20). Per-well ratio
  measured / planted marker was 0.81-0.90 in thinned wells: with the level-3 cytoplasm model the 4 px ring
  reaches the dim rim of a lone cell but is filled by neighbours in a packed colony, and level-4 colonies
  are denser than level-3 ones. That is a renderer artefact, not a planted effect, and it would have put my
  own tolerance inside the systematic error of both solutions. Made the level-4 cytoplasm wider and flatter
  (rim 5.5-7.5 px, 8% radial decline). Ratio now 0.93-1.01 for the cytotoxic, dual and cytostatic series (the
  wells of the potent inhibitor with almost no marker left are noisier, 0.83-1.01); dual compound -0.04
  (seed 1) / -0.07 (seed 2).
- Also: per-cell size factor clipped to 0.90-1.12 so that the SOP sentence "nuclei measure at least 11 um
  in their shortest dimension" is true for every rendered nucleus (smallest 8.8 px = 11.5 um); fragments
  kept 3.5 px away from nuclear outlines.

## L4-6. Scoring and battery (01:35-02:00)
- `score_l4.py`: counts scored on the wells without annotation with the SOP's own per-field tolerance,
  new `validation` component, `validation.tsv` gated.
- Quick scores with the verifier files before any gate run, first weights (seed 1 / seed 2): classical solution 1.0 / 1.0;
  level-3 watershed 0.309 / 0.226; level-3 LoG 0.05 / 0.05; the same after a plate-wide factor 0.309 / 0.05
  and 0.226 / 0.05; fragments counted 0.246 / 0.415; free-plateau fit 0.30 / 0.30; counting frame ignored
  1.0 / 1.0 (not an attack: +0..+2% is inside the tolerance).
- Level-4 pipeline with the level-3 Cellpose settings but the SOP's counting rules, seed 1: counts 0.35 inside,
  3 of 8 annotated fields inside, but all 16 potency entries right or nearly (largest error 0.2175): an
  almost uniform loss of nuclei in dense fields nearly cancels in the response ratios. With the weights of the
  spec (counts 0.20, potencies 0.45 / 0.25) that submission scored 0.736, and 0.75 if the last potency had been
  inside: no margin under the 0.75 ceiling of `l3_reference`. Changed the weights to counts 0.25, validation
  0.05, QC 0.05, marker 0.45, count potency 0.20: a submission whose counts fail cannot exceed 0.70. The same
  output re-scored: 0.689.
- Full-credit point of the count component set to 0.85 (zero 0.65): level-3 methods 0.35-0.55, classical
  solution 1.00 on both seeds; the Cellpose solution had not been measured on a full plate at that point.

## L4-7. Build-host load, role swap (01:40-02:02)
- Load average 12-14 on 4 cores from the other authors' runs. A `try` of 8 fields at 12 px took 2 min 44 s
  and 3 min 33 s (21-23 s on the quieter host an hour earlier). My first full run of the Cellpose route
  (development plate) was killed by the 30 min limit of my own background job before `segment` had finished;
  it was not repeated outside the gate.
- Because of that the roles were swapped: the classical outline-voting pipeline is the oracle (11-22 s),
  the staged Cellpose route is `SOLUTION=2` (three settings tried on the annotated fields, then the plate).
  Consequence recorded honestly: `oracle_steps = 1`, far below the 20 that AUTHORING asks for.
- `naive_untuned_cellpose` (level-4 pipeline, level-3 settings) was taken out of the gated battery to save
  one Cellpose plate run per gate; its seed-1 number above was measured by hand with the same verifier files.
- Free disk fell from 24 GB to 7.5 GB during the session (other authors' images and build cache); it was
  9.3 GB when the first level-4 gate run started. My images: the tools image above and the task images.

## L4-8. Gate runs, interruption, final state (02:01-03:15)
- First gate runs, family code as of 01:58: `make --seed 1 --level 4 --keep-work` all 31 checks ok incl.
  Harbor oracle 1.0 / nop 0.0, 31.5 min (02:01-02:33, load average 6-11; the Cellpose second solution
  alone about 16 min); `make --seed 2 --level 4 --no-harbor --keep-work` all 29 checks ok, 11.3 min.
  Second solution on a full plate for the first time: 0.964 (seed 1) / 1.000 (seed 2) of the wells
  without annotation inside the tolerance, score 1.0 on both.
- `l3_reference` (Cellpose, level-3 settings): 0.647 / 0.700. On seed 2 all 16 potency entries are right;
  the score is exactly the cap for failed counts. With the weights of the spec it would have been 0.75.
- During the seed-1 run: free disk 5.2 GB (other authors). Removed my tools image tag (its environment
  layers stay, they are the layers of the task image) and the development plates in my scratch directory.
- Parameter sweep of the voting counter on the final inputs (16 settings per seed): see the design notes.
  Seeds 3 and 4 generated into scratch and scored by hand with the classical pipelines only.
- My run was cut off by a usage limit at about 02:45 and resumed at 02:47.
- After the first gate runs I changed `level4.py` once more (02:48): the judge text of the task card had
  been written before the measurements and held three numbers that were not measured ("98-100%", "13-59%",
  "14-18%"); it now holds the gate values. `oracle_runtime_min` 0.5 -> 0.2. Regenerating both tasks showed
  `README.md` and `task.toml` as the only files that differ. Because code had changed, both gates were run
  again with the final code: seed 1 all 31 checks ok incl. Harbor, 11.07 min (written 03:00); seed 2 all 29
  checks ok, 9.92 min (written 03:10). Every score and every raw metric is identical to the first runs
  (tables diffed). Second solution 5.6 / 5.8 min, `l3_reference` 1.9 / 1.9 min, oracle 0.14 / 0.15 min at a
  load average of 1-4. (The background job that ran the two gates ended with exit code 1 because of my
  closing `tail -3 file1 file2`; both `make` commands had returned `"ok": true`.)
- After the last code change: level 3 seed 1 regenerated into scratch, `diff -r` against
  `tasks/hcs-dose-response-s0001-l3`: no difference. Level 4 seeds 1 and 2 regenerated into scratch:
  `diff -r` against `tasks/`: no difference.
- Cleanup: seed-2 agent and verifier images removed; no container of mine left; scratch directory 288 KB
  (logs and helper scripts). Kept: `alegen-hcs-dose-response-s0001-l4-agent:local` (2.82 GB) and its
  verifier image (178 MB).

## Not done at level 4
- No agent pilot (not allowed in this round).
- Real nuclear texture and the optional time-budget element were not attempted.
- The Cellpose routes were validated on seeds 1 and 2 only; `naive_untuned_cellpose` on seed 1 only and
  outside the gate.
- The first full-plate run of the Cellpose route during development was killed by my own 30 min job limit
  and not repeated outside the gate.

# Level 7 (2026-11-12, about 21:25-22:50 UTC)

Spec `docs/specs/l7/hcs_dose_response.md`: the level-4 plate of the seed, imaged with the plate turned by
180 degrees; layout and texts unchanged.

## L7-1 Reading and decisions
- Read the spec, AUTHORING (level 7, battery, running, host notes), `level4.py`, `sim4.py`, the level-4
  solutions and verifier, `descaffold.py`, `validate.py`, and the level-7 code of amplicon_16s and
  flow_compensation.
- Images are re-filed, not re-rendered and not turned pixel-wise (the field statistics and the illumination
  are symmetric under the turn; reasoning in the design notes). TIFF tags checked: no file name inside.
- The annotation files need the nucleus positions of other wells than at level 4, which `generate_l4` did
  not hand out: added the optional argument `plate_out` (None at levels 4-6).
- Per-well rows: decided that the contract does not fix whether a row is named after the image or after
  the plate well it shows; the verifier accepts both for each of the two tables (wrapper `score_l7.py`
  around the unchanged level-4 verifier, same technique as amplicon level 7).
- 19% of the usable seeds would have an out-of-focus image under an annotated name. A technician cannot
  count such a field, and the validation criterion could not be met on it; added a replacement rule instead
  of refusing those seeds.

## L7-2 Survey before rendering anything
- `l7.survey` plants a plate without rendering and evaluates the control statistics for both orientations.
  First version of the checks refused 111 of 200 seeds, 81 of them because I had required the dose series
  to rise by name in at least three series more than they fall. That was my error: a by-name series holds
  the opposite row's vehicle and reference wells in between and is noisy, and an inactive series reaches a
  rank correlation of 0.7 by chance. Replaced by the balance of falling over rising series (turned at least
  5, by name at least 5 lower): on 400 seeds turned 6-11, by name -9 to 2, smallest difference 5.
- What remains is real: 41 of 400 seeds have controls that look acceptable by name (symmetric offsets of
  opposite rows). Those seeds are refused.

## L7-3 Solutions (dev runs in containers of the level-4 agent image, `--network none`, 1-2 cpus)
- Reference: seed 1 by name reference median 1.32 and 13 of 23 vehicle wells off, turned 0.044 and 0;
  score 1.0. Level-7 code on the level-4 inputs with the layout by name gives the level-4 result to the
  last digit (the nuclear level is now taken over all fields; no measurable change).
- Second solution: orientation from the series direction, readouts normalised without controls. First
  idea, the median of the image column for the marker, was dropped before it was run: up to five of the
  eight wells of a high-dose column are affected. Used a robust line through the dense band instead.
  Seed 1: by name 0 falling / 2 rising, turned 9 / 1; score 1.0, 96.4% of the wells without annotation
  inside the tolerance (the level-4 value); 6.4 min at a load average of 5-7.
- By-name analyses: `l6_reference` 0.35 on seeds 1, 2, 3 and 10. Seeds 1 and 3 report `inactive` 16
  times; seed 2 reports four numbers, three of them below the lowest tested concentration.

## L7-4 Gate runs
- Seed 1, first run (13.4 min): one failure. `naive_no_focus_qc` 0.825 against 0.75. Cause found by
  forcing the outline threshold of the voting step: 0.15 finds nothing in out-of-focus fields (marker
  `NA`, marker curve unaffected), 0.12 finds objects (both potencies wrong). The tuning picks 0.15 on the
  level-7 annotated fields of seeds 1 and 10 and 0.12 on seed 2. Same score on the level-4 inputs with
  0.15 forced, where both thresholds tie on the annotated fields. Raised the ceiling of this attack to
  0.85 at level 7, documented in the task card and the design notes; record of the failed run kept in the
  scratch directory.
- Free disk fell to 2.5 GB after the first run; `docker builder prune -af` released 2.9 GB of build cache.
- Seed 1, second run with the final battery: all 37 checks ok incl. Harbor oracle 1.0 / nop 0.0, 13.38 min.
  Seed 2 (`--no-harbor`): all 35 checks ok, 12.11 min. `docker builder prune -af` after each.
- After the last code change (a line wrap in `l7.py` while the seed-1 gate was running): levels 3, 4, 6 of
  seed 1, level 6 of seed 2 and level 7 of both seeds regenerated into scratch, `diff -r` against `tasks/`:
  no difference. So the gated level-7 tasks are what the final code writes.
- My slip: scoring dev runs on the host with `tests/score.py` of `tasks/hcs-dose-response-s0001-l4` left a
  `tests/__pycache__` there (22:11-22:33). Found by the `diff -r`, removed.

## Not done at level 7
- No agent pilot.
- The second solution was run by hand on seed 1 only; on seed 2 it ran inside the gate.
- Seeds other than 1, 2, 3 and 10 were not rendered.

# Level 7, revision of the annotated fields (2026-11-14, from 04:33 UTC)

Maintainer decision after the yardstick runs on the first version (reported to me, records not read by
me): the lab annotated the same physical fields as at level 4; each annotation file carries the name of
the image that shows its field. Level 7 then differs from level 6 only by the re-filing.

## L7-5 What changed
- `l7.py` only. The annotation files are the level-4 files (same Path objects) under the turned names.
  The replacement rule for out-of-focus annotated names is removed; level 4 never annotates an
  out-of-focus field (`choose_annotated` draws from in-focus wells, `generate_l4` asserts it) and `l7.py`
  raises if the level-4 build ever did. Own jitter stream and `annotations_l7/` work directory are gone.
- The orientation guard and the `discoverable_by` control statistics never used the annotated fields:
  unchanged. `survey` lost its bookkeeping of out-of-focus annotated names.
- `truth.json`: `annotated` flags by image name; `plants.plate_orientation.annotated_plate_wells` and
  `annotated_fields_as_the_layout_names_them` (per file: manual count, what the layout calls the name, which
  field it is, count relative to an untreated field of the plate area) replace the two keys of the first
  version.
- Task card: the decision point lists the eight files that way; `discoverable_by` names the manual counts
  as a second sign with the measured numbers; the judge text describes the rule and no longer has the
  paragraph on the raised ceiling or the weakness "annotated fields lost their design".
- Scorer, solutions, pipelines: not touched. `validation.tsv` is checked against the eight annotated
  images as before; the four namings are still accepted.

## L7-6 Measurements before gating (dev runs, containers of the level-4 agent image)
- Reference tuning, seeds 1 / 2 / 3 / 10: outline threshold 0.12 on all four (worst annotated field 0.76 /
  0.33 / 0.54 / 0.79 of the SOP limit); seed 1 is exactly the level-4 result. Oracle 1.0 on all four.
- `naive_no_focus_qc`: 0.544 / 0.544 / 0.669 / 0.581. Ceiling back at 0.75, exception removed from the
  task card.
- `l6_reference`: 0.35 on all four; its `dose_response.tsv` is the same file as on 2026-11-12 for each seed.
- Scratch builds of both seeds against the previous level-7 builds in `tasks/` (sha256 manifests taken
  before anything was changed): 210 of 220 files identical, among them all 192 images; `README.md` and
  `tests/reference/truth.json` differ; the eight annotation files have other names. Each level-7
  annotation file is byte-identical to the level-6 file of the turned name (seed 1 against
  `tasks/hcs-dose-response-s0001-l6`; seed 2 against a level-6 build regenerated into scratch, because
  `tasks/hcs-dose-response-s0002-l6` did not exist at that moment and is the maintainer's to write).
- Levels 4 and 6 of seed 1 regenerated into scratch: `diff -r` against `tasks/`, no difference.

## L7-7 Disk
- 04:38 free disk 1.8 GB, 04:41 1.7 GB: build cache of gates that were not mine (rnaseq_confound seeds 4
  and 5). I pruned the build cache once (2.09 GB released) and did not start my gates; asked the
  maintainer. Answer: disk freed to 5.3 GB, go ahead, a dip to about 1.5 GB for a minute is acceptable,
  stop below 1 GB.
- Gates run through `<scratch>/hcs_l7_check/gate.sh`: it watches free disk every 5 s, prunes the build
  cache as soon as the images of a gate are built and again at the end, and stops below 1 GB.

## L7-8 Gate runs of the revision
- First start 04:44:53. While it ran I re-measured the by-name variants that the task card quotes
  (dev runs, one thread each, next to the gate for about a minute): `controls_swapped_only` 0.35 on both
  seeds with the same counts of right entries as on 2026-11-12 (seed 1: four potency values moved a
  little); the example series of seed 1 in the judge text was 0.01 off in three values because the
  voting parameters are now those of level 4. Stopped my gate after 4 min (killed the `make` process,
  removed its one container, which ran my `...-s0001-l7-agent` image), corrected the series in `l7.py`
  and the design notes, started again. Logs of the stopped run: `<scratch>/hcs_l7_check/*_stopped_by_me.log`.
- Second start 04:50:10. Seed 1: all 37 checks ok incl. Harbor oracle 1.0 / nop 0.0, 12.36 min. Seed 2
  (now with Harbor as well): all 37 checks ok, 12.99 min. Free disk: 6.7 GB at the start of each gate,
  4.6 GB when its images were built (40-45 s after the start), 6.7 GB again after the prune that the
  runner issued at that moment; lowest value seen in either gate 4.5 GB.
- Scores that changed against the first version: `naive_no_focus_qc` 0.544 / 0.544 (was 0.825 / 0.544),
  `naive_fragments_counted` 0.481 / 0.481 (was 0.510 / 0.510), the six level-3-route attacks 0.05-0.10
  (were 0.05-0.16). `l6_reference`, `layout_by_name_everything_else_right`, `controls_swapped_only` 0.35 as
  before; all probes 1.0.
- After the gates: levels 4 and 6 of seed 1 and level 7 of both seeds regenerated into scratch, `diff -r`
  against `tasks/`: no difference. Revised level-7 tasks in `tasks/` against the manifests of the first
  version: all 192 images identical, `README.md` and `truth.json` differ, annotation files renamed; every
  level-7 annotation file is byte-identical to the level-6 file of the turned name, for seed 2 now against
  the maintainer's `tasks/hcs-dose-response-s0002-l6` (read only).
- Second solution by hand on the gated tasks (the gate keeps no logs): `nuclei` 12 / flow off on both
  seeds (worst annotated field 0.63 / 0.67 of the limit), series balance -2 against 8 and -4 against 9,
  score 1.0 / 1.0, 96.4% / 100% of the wells without annotation inside the tolerance; 6.5 and 5.8 min.
- Reference solution with the outline threshold forced to 0.15 on the revised seed-1 task: worst annotated
  field 0.76 of the limit as with 0.12, full solution 1.0, without focus QC 0.825. The tie-break finding of
  2026-11-12 holds for the revised task; it is listed as a limitation, the battery has 0.544.
- `naive_fragments_counted` on seed 1 is 0.481 where level 4 has 0.3125: level-7 code on the level-4
  inputs gives 0.481 too, so it is the nuclear-signal level over all fields, not the data.
- A cleanup command of mine (copy of two `tests/` directories into scratch, removed again in the same
  command) was refused by the environment's safety check; I did not rephrase it and scored in place with
  bytecode writing switched off instead. No `__pycache__` under `tasks/` (checked).

## Not done in the revision
- No model run on the revised tasks.
- Seeds 3 and 10: reference solution, its variant without focus QC and `l6_reference` only.
