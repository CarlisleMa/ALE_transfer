# Family `flow_compensation` - 8-colour PBMC flow cytometry: spillover, gating, frequencies

ALE subdomain 3.3 (Cell & Imaging Biology, cytometry). Code: `pipeline/alegen/families/flow_compensation/`.
Spec: `docs/specs/flow_compensation.md`. All numbers below were measured on 2026-11-11 with the
commands in "How the numbers were obtained".

## Scenario

An immunology core receives the raw export of one acquisition session: six PBMC samples stained
with an 8-colour lineage panel, eight single-stain controls (seven antibody conjugates on
capture beads, the viability dye on a live/heat-killed cell mix) and an unstained tube, plus the
operator's notes, the panel sheet and the lab's gating SOP. The analyst has to deliver the
spillover matrix of the panel and, for every sample, count and percent-of-parent of the twelve SOP
populations (cells, singlets, live, lymphocytes, monocytes, T cells, CD4/CD8/DP/DN T cells, B
cells, NK cells).

| detector | fluorochrome | marker | role |
|---|---|---|---|
| V450-A | BV421 | CD3 | T cells |
| V525-A | Zombie Aqua | Viability | dead-cell exclusion |
| B530-A | FITC | CD14 | monocytes |
| B710-A | PerCP-Cy5.5 | CD8 | T subsets |
| YG586-A | PE | CD4 | T subsets |
| YG780-A | PE-Cy7 | HLA-DR | recorded, not part of the reportable list |
| R670-A | APC | CD56 | NK cells |
| R780-A | APC-Cy7 | CD19 | B cells |

Agent image: Python 3.12, numpy 2.3.5, scipy 1.16.3, pandas 2.3.3, scikit-learn 1.7.2, FlowIO
1.4.0 (conda-forge), FlowKit 1.3.2 (pip), fcsparser 0.2.8 (pip, `--no-deps`, see limitations).
1.53 GB as reported by `docker images`; `/opt/tools` is 809 MB. Verifier image: `python:3.12-slim`,
no extra packages (score.py is stdlib only).

## Truth model (why the reference is correct)

Pure forward simulation in `__init__.py`; nothing is downloaded and no analysis tool is run at
generation time.

1. A spillover matrix `S` (rows = fluorochromes, columns = detectors, diagonal 1) is drawn from
   seed-dependent ranges: BV421->V525 0.18-0.32, PE->B710 0.18-0.32, APC->R780 0.20-0.38,
   PerCP-Cy5.5->YG780 0.15-0.30, 18 further entries 0.004-0.15, the rest 0-0.003.
2. Every sample gets seed-dependent proportions (debris 6-15% of events, aggregates 6-12% of cell
   events, dead 8-18% of single cells, monocytes 10-25% of live, T/B/NK/other and CD4/CD8/DP/DN by
   multinomial draws). Every event is created from its class: fluorochrome amounts (log-normal,
   positives 5 000-40 000, i.e. more than 1.5 decades above the unstained level of 40-120),
   autofluorescence, FSC-A/FSC-H/SSC-A. Aggregates are the sum of two independently drawn cells
   (fluorescence and SSC add, FSC-A = big + 0.8-1.0 x small, FSC-H ~ the bigger cell). Dead cells
   are viability-dye bright and carry a per-cell non-specific smear in every antibody channel.
3. Detector signal = amounts x `S` + autofluorescence + photon noise (variance = 6-12 x signal)
   + electronic noise (sd 20-35). Values are written unclipped below zero, clipped at 262 143.
4. Controls use the same `S`, noise and detector model (beads: 50% positive capture beads, own
   autofluorescence, 3% bead aggregates; viability control: 42% heat-killed cells).
5. Files are written by a small FCS 3.1 writer in the family (`write_fcs`: float32, little
   endian, `$PnN/$PnS/$PnR/$PnB/$PnE/$PnV`, `$CYT`, `$DATE`, `$BTIM/$ETIM`, `$TIMESTEP`,
   `$SPILLOVER` at level >= 2, FACSDiva-style form-feed delimiter).

`reference/truth.json` holds `S` (rounded to 4 decimals before it is used for the simulation) and,
per sample, the number of events of each SOP population **counted from the hidden class labels**.
Population membership by label: cells = everything but debris; singlets = cells minus aggregates;
live = singlets minus dead; lymphocytes = T + B + NK + lineage-negative lymphocytes; monocytes;
T = CD4 + CD8 + DP + DN. By construction every live single cell is either CD14+ with high scatter
(monocyte) or CD14- with low side scatter (lymphocyte), so the SOP's marker logic and the labels
describe the same sets.

Separation actually obtained after compensation with the planted matrix (one level-3 sample drawn
with the generator functions, live single cells; detector units):

| channel | negative median | negative 99.9th pct | positive 0.1th pct | positive median |
|---|---|---|---|---|
| CD3 | 100 | 399 | 16230 | 34932 |
| CD4 | 78 | 540 | 13070 | 30371 |
| CD8 | 49 | 937 | 3794 | 8549 |
| CD14 | 117 | 407 | 5516 | 14430 |
| CD19 | 24 | 874 | 2404 | 5659 |
| CD56 (NK vs CD56- non-T) | 63 | 308 | 8001 | 22001 |
| Viability (live vs dead) | 188 | 1089 | 3535 | 12793 |

The gap between the 99.9th percentile of the negatives and the 0.1th percentile of the positives
is at least a factor 2.7 in every gated channel, which is why two unrelated ways of placing the
boundary give the same counts to within 11 events.

Reader check (host dev environment and inside the agent image, all 15 files of seed 1 level 3):
FlowKit 1.3.2 `Sample` reports version 3.1 and the right event counts; flowio and fcsparser return
bit-identical arrays; FlowKit's raw events equal flowio's except the Time channel, which FlowKit
scales by `$TIMESTEP`; `Sample.apply_compensation` with the planted matrix equals `X @ inv(S)`.

## Decision points and how they can be discovered

| trap (level) | what goes wrong if ignored | evidence available to the agent |
|---|---|---|
| `compensation_required` (2, 3) | The exported files carry an identity `$SPILLOVER`. Uncompensated, CD3-BV421 puts live T cells at the level of dead cells in V525-A, CD4-PE puts CD4 T cells at the level of CD8+ cells in B710-A, CD56-APC puts NK cells at the level of B cells in R780-A. | `acquisition_notes.md` ("compensation setup ... was not run", "template compensation", "without any correction of the data"); keywords `$SPILLOVER` = identity and `APPLY COMPENSATION` = `FALSE`; single-stain controls show the slopes; any V450-A vs V525-A plot of a sample shows a diagonal population |
| `doublets` (2, 3) | 6-12% of the cell events are two cells. Their summed fluorescence roughly doubles to quadruples DP_T and inflates the lymphocyte count by 7-11%. | SOP lists `singlets` ("events that consist of exactly one cell"); FSC-A/FSC-H ratio is bimodal (single cells median 1.09, maximum 1.19; aggregates minimum 1.38, median 1.86 in a test draw); notes: tubes not strained |
| `dead_cells` (3) | 8-18% of single cells are dead and stain non-specifically; left in, they add 24-27% to the lymphocyte count and move T/CD4/CD8/DN/B/NK by several points. | SOP defines `live` as viability dye negative; notes: cryopreserved, 80-92% trypan-blue viability; compensated V525-A of singlets is bimodal |
| `per_sample_shift` (3) | One sample (never the first of the sheet, index drawn from the seed) was re-acquired with the FSC voltage lowered (gain x0.20-0.27). A debris/cell gate copied from another sample discards its lymphocytes. | notes state the re-acquisition and both voltages; `$P1V`/`$P2V` differ in that file only; its FSC-A histogram is shifted |

Level 1 has no aggregates, 0.4-1.2% dead cells, no voltage change and no `$SPILLOVER` keyword
(the notes say that nothing was compensated). Compensation itself is part of the workflow at every
level; the level >= 2 trap is the misleading identity keyword.

## Output contract

- `spillover.csv`: header `fluorochrome,<8 detectors>`, one row per fluorochrome, entry = signal of
  the fluorochrome in that detector relative to its primary detector. The verifier matches rows
  and columns by label (case/punctuation-insensitive) and also accepts the transposed layout.
- `frequencies.tsv`: `sample`, `population`, `parent`, `count`, `percent_of_parent`; one row per
  sample x SOP population; extra columns and rows for other populations are ignored.

## Scoring (pass threshold 0.9)

Hard gates: files present and parseable; matrix labelled with the panel's fluorochromes and
detectors, numeric; every sample x population exactly once; `parent` as in the SOP; `count` a
non-negative integer; `percent_of_parent` in [0, 100].

| weight | component | metric | full | zero |
|---|---|---|---|---|
| 0.25 | spillover_matrix | max abs error over the 64 entries | <= 0.02 | >= 0.05 |
| 0.40 | frequencies_every_sample | per sample, fraction of the 12 populations with `percent_of_parent` within max(0.5 points, 3% of truth); minimum over samples | 1.0 | <= 0.75 |
| 0.20 | frequencies_overall | fraction of all 72 items within the same band | 1.0 | <= 0.5 |
| 0.15 | lymphocyte_counts | fraction of samples with the `lymphocytes` count within 3% | 1.0 | <= 0.5 |

Consequence of the weights: one item out of band in one sample gives 0.833 (fail); a submission
needs every one of the 72 percentages inside the band to pass.

### Measured raw metrics (verifier `details.json`, runs inside the task images)

`worst/band` = largest |percent - truth| divided by its tolerance band over the 72 items;
`min frac` = smallest per-sample fraction of populations in band; `lymph err` = largest relative
error of the lymphocyte count.

| submission | seed 1 L3: spill max err | worst/band | min frac | lymph err | score | seed 2 L3: spill max err | worst/band | min frac | lymph err | score |
|---|---|---|---|---|---|---|---|---|---|---|
| oracle (valley gates, flowio) | 0.00042 | 0.030 | 1.0 | 0.0003 | 1.000 | 0.00086 | 0.017 | 1.0 | 0.0003 | 1.000 |
| independent (mixture gates, fcsparser) | 0.00058 | 0.017 | 1.0 | 0.0005 | 1.000 | 0.00053 | 0.012 | 1.0 | 0.0002 | 1.000 |
| probe within_tolerance | 0.00788 | 0.601 | 1.0 | 0.015 | 1.000 | 0.00788 | 0.601 | 1.0 | 0.015 | 1.000 |
| naive_matrix_not_applied | 0.00042 | 170.3 | 0.167 | 0.768 | 0.250 | 0.00086 | 164.9 | 0.167 | 0.704 | 0.250 |
| naive_embedded_matrix | 0.3019 | 170.3 | 0.167 | 0.768 | 0.000 | 0.3090 | 164.9 | 0.167 | 0.704 | 0.000 |
| naive_no_doublet_gate | 0.00042 | 5.46 | 0.500 | 0.094 | 0.306 | 0.00086 | 4.21 | 0.583 | 0.083 | 0.322 |
| naive_no_viability_gate | 0.00042 | 49.5 | 0.250 | 0.253 | 0.250 | 0.00086 | 47.9 | 0.167 | 0.244 | 0.250 |
| naive_global_thresholds | 0.00042 | 169.4 | 0.000 | 1.000 | 0.483 | 0.00086 | 122.4 | 0.000 | 1.000 | 0.483 |
| naive_default_pipeline (all four) | 0.3019 | 136.3 | 0.000 | 0.916 | 0.000 | 0.3090 | 140.3 | 0.000 | 0.928 | 0.000 |
| textbook_prior | 0.3019 | 19.2 | 0.167 | 0.292 | 0.000 | 0.3090 | 13.3 | 0.250 | 0.155 | 0.000 |
| all_percent_beyond_tolerance (1.4 x band) | 0 | 1.40 | 0.000 | 0 | 0.400 | 0 | 1.40 | 0.000 | 0 | 0.400 |
| one_population_beyond_tolerance (1.5 x band) | 0 | 1.50 | 0.917 | 0 | 0.833 | 0 | 1.50 | 0.917 | 0 | 0.833 |
| matrix_entry_beyond_tolerance (+0.06) | 0.06 | 0.001 | 1.0 | 0 | 0.750 | 0.06 | 0.001 | 1.0 | 0 | 0.750 |

Seed 1 level 1: oracle 0.00042 / 0.010 / 1.0 / 0.0002 -> 1.000; independent 0.00058 / 0.002 / 1.0 /
0.0000 -> 1.000; naive_matrix_not_applied 0.00042 / 176.5 / 0.167 / 0.723 -> 0.250.

Host-side repeat of the two solutions and the six naive variants on level-3 seeds 1-6 (same code,
dev environment with numpy 2.5.4 / scikit-learn 1.9.1): oracle and independent solution 1.000 on
all six seeds with worst/band <= 0.030; naive_matrix_not_applied 0.250 (6/6), naive_embedded_matrix
0.000 (6/6), naive_no_doublet_gate 0.272-0.328, naive_no_viability_gate 0.250 (6/6),
naive_global_thresholds 0.483 (6/6), naive_default_pipeline 0.000 (6/6).

Examples of what the trap-ignoring analyses report (seed 1 level 3, first affected sample; value
vs planted): uncompensated T_cells 0.00 vs 76.49, B_cells 94.38 vs 11.78, NK_cells 0.00 vs 9.10;
no doublet gate singlets 100.00 vs 88.21, DP_T 3.66 vs 0.93; no viability gate live 100.00 vs 85.30,
CD4_T 49.36 vs 61.45, DN_T 22.48 vs 3.85; one gate template (shifted sample) cells 7.67 vs 90.33,
lymphocytes 0.00 vs 84.72.

### Why the bands sit where they do

- Percentages: the two correct solutions deviate from the planted value by at most 3% of the band
  (largest absolute deviation 0.038 percentage points, largest count difference 11 events, over
  seed 1 L3, seed 2 L3 and seed 1 L1). The mildest trap-ignoring analysis is the missing doublet gate:
  in every sample 3-6 of the 12 items are out of band and the worst item of each sample is 2.2-5.5
  bands off (seeds 1 and 2; no viability gate: 8-10 items per sample, worst 7.9-49 bands; one gate
  template: all 12 items of the shifted sample, worst 122-169 bands; no compensation: 10 items per
  sample, worst 157-170 bands). The band max(0.5 points, 3%) therefore leaves a factor > 30 towards
  the correct side and > 2 towards the wrong side. It was not made tighter because a human-style
  gate that trims 1-2% of a population is still a correct analysis.
- Matrix: both estimators are within 0.0012 on six seeds. Not subtracting the bead background
  should cost at most ~0.004 (autofluorescence 35-70 against 20 000-70 000 bead signal; calculated,
  not measured). The identity matrix is off by 0.29-0.37. Full credit <= 0.02, none >= 0.05.
- Lymphocyte count: correct solutions within 0.05%; no doublet gate +7-10%, no viability gate
  +24-27%, shifted sample under one gate template -100% (six level-3 seeds). Band 3%.
- Ceilings of the trap-ignoring attacks are set 0.15-0.2 above their measured scores and are all
  <= 0.65, i.e. below pass_threshold - 0.15 = 0.75. The three truth-derived "just beyond tolerance"
  attacks have ceilings 0.45, 0.85 and 0.80 against measured 0.400, 0.833 and 0.750; the last two
  are above 0.75 on purpose: they show that a single out-of-band item, or a matrix that is wrong
  with perfect frequencies, is enough to fail (pass threshold 0.9), not that such a submission is
  worthless.

## The two solutions

| step | oracle `flow_pipeline.py` | independent `solve_alt.py` |
|---|---|---|
| reader | flowio | fcsparser |
| spillover | median(positive) - median(negative) per detector, populations split at a density valley | least-squares slope over all events of the tube |
| compensation | `X @ inv(S)` | `numpy.linalg.solve` |
| transform | arcsinh, cofactor = 3 x rms of the below-zero values (>= 150) | arcsinh, cofactor = 2 x abs(1st percentile) (>= 150) |
| debris | valley of log10 FSC-A | 2-component Gaussian mixture on log10 FSC-A |
| singlets | FSC-A/FSC-H below median + 6 robust SD | 2-component mixture on log(FSC-A/FSC-H) |
| marker +/- | valley between the two dominant modes of the parent population | posterior boundary of a 2-component mixture |
| lymphocytes / monocytes | CD14 valley and SSC-A valley | 2-D mixture on (CD14, log SSC-A) |

## Known limitations and open questions

- **The data are clean.** Populations are separated by design (spec: >= 1.5 decades), so once the
  four decisions are right any reasonable gating lands within 0.03 bands. The difficulty is in the
  decisions, not in gating finesse; real PBMC data are messier (CD56 dim NK cells, CD4 dim
  monocytes, CD8 dim NK cells are not simulated).
- **`per_sample_shift` is a scatter-detector change only.** A changed fluorescence PMT voltage
  would invalidate the single-stain controls for that sample (every spillover coefficient into and
  out of that detector scales with the gain), and the truth would then not be recoverable without a
  second control set. I did not have time to build and validate that variant, so the voltage
  change is on FSC, where compensation is unaffected. Fluorescence thresholds therefore do not
  need to move between samples (only a +-8% donor/staining factor).
- **The shifted sample is never the first sample of the sheet**, because the
  `naive_global_thresholds` attack takes its template from the first sample. I expect a template
  taken from the shifted sample to fail as well (its debris threshold lets the debris of the other
  samples through), but I did not measure that case.
- **Dead cells lose their lineage identity** in the simulator (viability dye plus a non-specific
  smear). Real dead cells keep part of their specific staining.
- **A correct analysis without compensation is conceivable**: in two dimensions the uncompensated
  populations still form separate diagonal clusters, so a multivariate clustering could assign
  cells correctly without subtracting spillover. Such a submission would be scored on its merits
  (it still needs the matrix for 0.25). Only one-dimensional thresholds on raw signals are covered
  by the attack suite.
- **Strictness.** One percentage out of band in one sample caps the score at 0.833. This is
  intended (every trap leaves several items out of band), but an otherwise careful submission
  with a single slip fails.
- **fcsparser 0.2.8 declares `numpy<2`**; it is installed with `--no-deps` next to numpy 2.3.5
  because FlowKit 1.3.2 needs numpy >= 2. It parsed all files in the image and the independent
  solution uses it, but the combination is not supported upstream.
- Spillover of beads and cells is identical by construction (no tandem-dye lot effects); the
  unstained control is provided but not needed by either solution.
- `lymphocyte_counts` is partly redundant with the percentages of `cells`...`lymphocytes`; it
  mainly penalises compensating errors along the hierarchy.
- The textbook-prior attack reaches up to 40% of the items in band on seed 2 (score still 0); with
  the 0.5 lower end of the `frequencies_overall` ramp a luckier prior could earn a few hundredths.
- Level 2 was validated once (seed 1, all gates except Harbor, all ok); levels 1 and 3 as required.
  Harbor (G6) was run for seed 1 level 3 only.

## How the numbers were obtained

```bash
cd /home/user/ale-lifesci-gym/pipeline && export PATH="$HOME/.local/bin:$PATH"
../.venv/bin/python -m alegen make flow_compensation --seed 1 --level 3 --keep-work
../.venv/bin/python -m alegen make flow_compensation --seed 2 --level 3 --no-harbor   # repeated with --keep-work for details.json
../.venv/bin/python -m alegen make flow_compensation --seed 1 --level 1 --no-harbor   # repeated with --keep-work for details.json
../.venv/bin/python -m alegen make flow_compensation --seed 1 --level 2 --no-harbor
```
Raw metrics are the `info` block of `results/work/<task_id>/{oracle,probe/*,attack/*}/logs/details.json`.
No external data are downloaded by this family (nothing under `.cache/`).

---

## Level 4 (added 2026-11-12)

Spec: `docs/specs/l4/flow_compensation.md`. Code: `l4.py` (generator, texts, battery),
`assets/flow_pipeline_l4.py` (reference pipeline and careless variants), `assets/solve_alt_l4.py`
(independent second solution), `assets/score_l4.py`, `assets/solve_l4.sh`. Levels 1-3 are generated
by the unchanged code path: level 3 seed 1 regenerated after the change is identical to
`tasks/flow-compensation-s0001-l3` under `diff -r` (nothing ignored); levels 1 and 2 (seed 1) and
level 3 (seed 2) are identical to copies generated before the change.

All numbers in this section were measured on 2026-11-12. "Docker" = inside the task images during
`python -m alegen make`; "host" = the same scripts in a scratch environment with the image's
package versions (numpy 2.3.5, scipy 1.16.3, pandas 2.3.3, scikit-learn 1.7.2, flowio 1.4.0,
fcsparser 0.2.8). Host and Docker scores agree for every submission of seeds 1 and 2.

### Why, and what was added

The pilot (README section 6) solved level 3 seed 1 with 1.0 in 2 min 40 s and 20 tool calls. The
level-3 limitations above say why: populations 1.5 decades apart (any gate lands in band), the
voltage change on a scatter detector only, bead and cell spillover identical. Level 4 keeps the
four level-3 decision points and adds four elements, all of them in the one detector that level 3
recorded but did not score (YG780-A, HLA-DR PE-Cy7).

| element | what is planted | pilot finding it answers | evidence available to the agent |
|---|---|---|---|
| `P_embedded_matrix` | Every file carries a well-formed, non-identity `$SPILLOVER` (`APPLY COMPENSATION` = `TRUE`) from the experiment template: the planted matrix re-scaled as if every detector had another gain (factors exp(U(-0.35, 0.35)), B710-A x1.6-2.4, YG780-A x0.7-1.0) plus 12% jitter. PerCP-Cy5.5 -> YG780-A is 0.093 instead of 0.254 (seed 1), 0.093 instead of 0.205 (seed 2). | "A plausible shortcut that the generator had closed": the level-3 keyword was an identity matrix, which nobody mistakes for a calibration. | Notes: date of the template, date of the later voltage re-optimisation, compensation setup not run. Single-stain controls give other coefficients. With the embedded matrix at least one population of the SOP's check is missing in the FMO tubes of all eight host seeds (CD8_T in seven), so the check cannot even be evaluated. |
| `P_voltage_change` | The YG780-A voltage is raised after 2-4 of the 6 samples (gain x1.8-2.6; the voltage exponent is drawn from 6.8-7.9, so the voltages alone do not give the factor). Recorded again: the PE-Cy7 bead tube and the FMO tube, nothing else. Planted matrix of the later samples: YG780-A column x gain, PE-Cy7 row / gain. | Level-3 limitation "voltage change on a scatter detector only". Shortcut that is closed: replace the one row that the new control gives. | Notes state the change and which tubes were recorded again; `$P9V` differs; the SOP's check on the re-recorded FMO gives a spread of 1413-3434 units (criterion 50) when only the row is replaced. |
| `P_tandem_lot` | The PerCP-Cy5.5 bead control was stained from another lot: PerCP-Cy5.5 -> YG780-A lower by 0.10-0.14, -> R670-A and -> R780-A x0.6-0.9. A cell-based single stain from the cocktail vial (30 000 events) is provided. | Level-3 limitation "spillover of beads and cells identical". Default that is wrong here: beads are the compensation control. | `reagent_lot` in `panel.tsv` and `controls.tsv`; notes say when and from which vials beads and the extra tube were stained (no sentence compares the lots); the two PerCP-Cy5.5 controls give different rows; SOP's check with the bead row: spread 828-1290 (session start), 1554-2620 (later setting). |
| `Q_fmo_gate` | HLA-DR PE-Cy7 on a per-sample fraction of the T cells (CD4 3-12%, CD8 10-35%), median 700 units, log-sd 0.75: a shoulder, no valley. Two new SOP populations (`activated_CD4_T`, `activated_CD8_T`) whose gate the SOP defines as the 99th percentile of the same parent population in the FMO tube of the same fluorescence-detector voltages, on compensated data. The SOP also states an acceptance criterion on the FMO tube: medians of the compensated YG780-A signal of CD4_T, CD8_T, DN_T, B_cells, NK_cells within 50 units, otherwise the compensation "has to be corrected". | "A final quantity that had to be good enough and needed self-checking and method tuning"; level-3 limitation "any reasonable gate lands in band". | `gating_sop.md` sections 1-3; the FMO tubes in `controls.tsv`. |

Inputs that changed against level 3: 13 control tubes instead of 9 (`ctrl_PerCPCy55_cells.fcs`,
`ctrl_PECy7_rerun.fcs`, `fmo_HLADR.fcs`, `fmo_HLADR_rerun.fcs`), a `reagent_lot` column in
`panel.tsv` and `controls.tsv`, control type `fmo`, SOP sections 1-3, notes. Input size 19.85 MB
(seed 1), 19.82 MB (seed 2). Agent image as at level 3 (same packages): 1.55 GB (`docker images`).
Oracle runtime in the agent image, wall time including container start: 0.13 min (seed 1) and
0.08 min (seed 2) while the build host had a load average of about 13, 0.03 min for seed 1 in an
earlier round on a quieter host; the reference pipeline itself takes 2 s on the host, the second
solution 9-22 s.

Deliverables that changed:
- `spillover.csv` is a long table with one 8 x 8 matrix per sample (`sample`, `fluorochrome`, eight
  detectors). The spec suggested one file per acquisition block; a per-sample table states the same
  thing without telling the agent how many blocks there are. The entry is defined "as it is present
  on the cells of that sample" and "at the detector settings with which that sample was recorded".
- `frequencies.tsv` has 14 populations per sample.

### Truth model of the new parts

- Matrix of a sample: the planted matrix of its block (sample lots; block 2 = block 1 with the
  YG780-A column multiplied and the PE-Cy7 row divided by the planted gain). The simulator applies
  the gain to optical signal and photon noise of YG780-A, not to electronic noise.
- Lineage populations: counted from the hidden labels, as at level 3.
- FMO-gated populations: there is no label that says "above the FMO gate"; the populations are
  operational by the SOP's definition. The reference applies that definition to planted quantities
  only: signals compensated with the planted matrix, parent population = hidden label, threshold =
  99th percentile (numpy, linear interpolation) of the label-defined parent in the FMO tube of the
  same block. It is not an output of either solution. For orientation: reference count divided by
  the number of cells that carry HLA-DR in the simulation is 0.75-1.05 (seed 1) and 0.87-1.33
  (seed 2); values above 1 occur in samples with few activated cells, where the 1% tail of the
  negative cells is a large share.

### Scoring at level 4 (pass threshold 0.9)

Hard gates as at level 3 (adapted to the long `spillover.csv`). Components:

| weight | component | metric | full | zero |
|---|---|---|---|---|
| 0.20 | spillover_matrices | max abs error over the 64 entries of every sample's matrix | <= 0.02 | >= 0.05 |
| 0.25 | lineage_every_sample | per sample, fraction of the 12 lineage populations within max(0.5 points, 3%); minimum over samples | 1.0 | <= 10/12 |
| 0.10 | lineage_overall | fraction of the 72 lineage items within the same band | 1.0 | <= 0.5 |
| 0.10 | lymphocyte_counts | fraction of samples with the `lymphocytes` count within 3% | 1.0 | <= 0.5 |
| 0.25 | activation_every_sample | per sample, fraction of the two FMO-gated populations within max(0.75 points, 6%); minimum over samples | 1.0 | <= 0.5 |
| 0.10 | activation_overall | fraction of the 12 FMO-gated items within the same band | 1.0 | 0 |

The spec asked to keep the level-3 structure. The components are the same kind, but the two new
populations got their own pair of components and the lineage ramp was steepened: with my first
weights a submission with one lineage item out of band in every sample scored exactly 0.900 and
passed (level 3: 0.833). Now 0.858; one FMO-gated item out in one sample: 0.742.

### Measured values (Docker, both seeds)

Columns per seed: largest matrix error | worst lineage item (error / band) | worst FMO-gated item
(error / band) | FMO-gated items in band | largest relative error of the lymphocyte count | score.

| submission | s1: matrix | lineage | FMO-gated | in band | lymph | score | s2: matrix | lineage | FMO-gated | in band | lymph | score |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| oracle (reference pipeline) | 0.0017 | 0.02 | 0.09 | 12/12 | 0.0004 | 1.000 | 0.0020 | 0.01 | 0.07 | 12/12 | 0.0002 | 1.000 |
| independent second solution | 0.0022 | 0.02 | 0.06 | 12/12 | 0.0003 | 1.000 | 0.0023 | 0.01 | 0.05 | 12/12 | 0.0004 | 1.000 |
| planted_truth | 0.0000 | 0.00 | 0.00 | 12/12 | 0.0000 | 1.000 | 0.0000 | 0.00 | 0.00 | 12/12 | 0.0000 | 1.000 |
| within_tolerance | 0.0079 | 0.60 | 0.60 | 12/12 | 0.0150 | 1.000 | 0.0080 | 0.60 | 0.60 | 12/12 | 0.0150 | 1.000 |
| l3_reference | 0.3298 | 0.06 | 16.67 | 0/12 | 0.0004 | 0.450 | 0.3696 | 0.13 | 16.67 | 0/12 | 0.0002 | 0.450 |
| l3_reference_last_control | 0.2292 | 0.02 | 16.67 | 0/12 | 0.0004 | 0.450 | 0.2460 | 0.01 | 16.67 | 0/12 | 0.0002 | 0.450 |
| careless_embedded_matrix | 0.3889 | 93.14 | 16.67 | 0/12 | 0.7635 | 0.000 | 0.3577 | 78.58 | 16.67 | 0/12 | 0.0002 | 0.175 |
| careless_one_matrix | 0.2292 | 0.06 | 123.68 | 8/12 | 0.0004 | 0.517 | 0.2460 | 0.14 | 127.55 | 4/12 | 0.0002 | 0.483 |
| careless_row_only | 0.2292 | 0.02 | 8.12 | 10/12 | 0.0004 | 0.533 | 0.2460 | 0.01 | 7.78 | 4/12 | 0.0002 | 0.483 |
| careless_bead_tandem | 0.1944 | 0.02 | 8.71 | 6/12 | 0.0004 | 0.500 | 0.2723 | 0.01 | 8.27 | 6/12 | 0.0002 | 0.500 |
| careless_no_fmo | 0.0017 | 0.02 | 42.96 | 0/12 | 0.0004 | 0.650 | 0.0020 | 0.01 | 47.50 | 0/12 | 0.0002 | 0.650 |
| careless_fmo_uncompensated | 0.0017 | 0.02 | 14.36 | 0/12 | 0.0004 | 0.650 | 0.0020 | 0.01 | 11.37 | 0/12 | 0.0002 | 0.650 |
| textbook_prior | 0.4818 | 12.40 | 5.96 | 0/12 | 0.1876 | 0.000 | 0.4509 | 19.78 | 7.34 | 1/12 | 0.2902 | 0.008 |
| block1_matrix_true_frequencies | 0.2281 | 0.00 | 0.00 | 12/12 | 0.0000 | 0.800 | 0.2456 | 0.00 | 0.00 | 12/12 | 0.0000 | 0.800 |
| all_percent_beyond_tolerance | 0.0000 | 1.40 | 1.40 | 0/12 | 0.0000 | 0.300 | 0.0000 | 1.40 | 1.40 | 0/12 | 0.0000 | 0.300 |
| one_population_beyond_tolerance | 0.0000 | 1.50 | 0.00 | 12/12 | 0.0000 | 0.858 | 0.0000 | 1.50 | 0.00 | 12/12 | 0.0000 | 0.858 |
| activation_beyond_tolerance | 0.0000 | 0.00 | 1.50 | 0/12 | 0.0000 | 0.650 | 0.0000 | 0.00 | 1.50 | 0/12 | 0.0000 | 0.650 |
| one_activation_item_beyond_tolerance | 0.0000 | 0.00 | 1.50 | 11/12 | 0.0000 | 0.742 | 0.0000 | 0.00 | 1.50 | 11/12 | 0.0000 | 0.742 |
| matrix_entry_beyond_tolerance | 0.0600 | 0.00 | 0.00 | 12/12 | 0.0000 | 0.800 | 0.0600 | 0.00 | 0.00 | 12/12 | 0.0000 | 0.800 |
| naive_matrix_not_applied | 0.0041 | 166.67 | 16.67 | 0/12 | 0.7865 | 0.200 | 0.0020 | 166.73 | 16.67 | 0/12 | 0.7684 | 0.200 |
| naive_no_doublet_gate | 0.0017 | 4.20 | 2.69 | 6/12 | 0.0833 | 0.278 | 0.0025 | 4.78 | 3.24 | 3/12 | 0.1026 | 0.253 |
| naive_no_viability_gate | 0.0017 | 46.88 | 29.04 | 8/12 | 0.2554 | 0.267 | 0.0020 | 50.34 | 23.59 | 8/12 | 0.2576 | 0.267 |
| naive_global_thresholds | 0.0017 | 120.01 | 16.67 | 10/12 | 1.0000 | 0.417 | 0.0020 | 33.33 | 15.78 | 10/12 | 1.0000 | 0.417 |
| naive_default_pipeline | 0.3889 | 61.05 | 16.67 | 1/12 | 0.9149 | 0.008 | 0.3577 | 51.87 | 33.32 | 0/12 | 0.9536 | 0.000 |

The four generic attacks of the validator (`empty_output`, `inputs_copied`, `zero_byte_files`,
`header_only`) score 0 at a hard gate on both seeds.

Ceilings: 0.75 (= pass threshold - 0.15) for `l3_reference`, `l3_reference_last_control` and the
six `careless_*` variants; 0.1 `textbook_prior`; 0.85 / 0.45 / 0.89 / 0.70 / 0.80 / 0.85 for the six
truth-derived perturbations in table order; 0.45 / 0.45 / 0.45 / 0.60 / 0.20 for the level-3
variants. The largest careless score is 0.667 (host, eight seeds), 0.083 below the ceiling.

What the submissions report for the FMO-gated populations (host, eight seeds = 48 samples;
deviation from the reference in percent of the reference value; the band is 6% or 0.75 points):
- reference pipeline -1.4% to +1.0%, second solution -0.8% to +0.6% (96 items each, all in band).
- bead lot for PerCP-Cy5.5: `activated_CD8_T` -74% to -17% in the earlier and -66% to +6.9% in the
  later samples, 0 of 48 in band; `activated_CD4_T` -3.5% to +9.6%, 48 of 48 in band.
- only the PE-Cy7 row replaced: in the later samples `activated_CD8_T` -62% to +8.6%, 0 of 24 in
  band, `activated_CD4_T` -46% to +46%, 9 of 24 in band; the earlier samples are unaffected.
- one matrix and one FMO for all samples: later samples +270% to +2490%.
- gate from the unstained control: +61% to +887%. FMO gate on uncompensated signals:
  `activated_CD8_T` -87% to -12%, 0 of 48 in band; `activated_CD4_T` -70% to +197%, 8 of 48 in band.
- level-3 method (valley rule, first-listed controls): no HLA-DR + event at all in 79 of 96 items,
  none in band.
- embedded matrix: fewer than 100 CD8_T events in 40 of 48 samples; lineage percentages 66-143
  bands off.

### The acceptance criterion, measured

The Q element requires that meeting the stated criterion is necessary and close to sufficient.

- Spread of the five FMO medians (units), host, eight seeds, session-start tube / re-recorded tube:
  reference pipeline before fine-tuning 5.6-13.1 / 12.5-35.8, after 0.2-6.9 / 0.3-12.2; second
  solution before 2.7-13.0 / 8.4-33.3, after 0.0-0.0 / 0.0-0.2; bead lot 828-1290 / 1554-2620;
  row-only 5.6-13.1 / 1413-3434; embedded matrix: not evaluable in 16 of 16 tubes (a population of
  the check is empty). In Docker (`qc_fmo.tsv` of the oracle run): seed 1 13.09 -> 0.59 and
  35.81 -> 0.35, seed 2 5.70 -> 1.32 and 12.46 -> 4.62.
- Sufficiency: sweep without fine-tuning, eight seeds x seven deliberate mis-settings of one
  coefficient of the session-start matrix (PE 0.0008 / 0.0016 / 0.0032, PerCP-Cy5.5 0.003 / 0.006 /
  0.012, APC-Cy7 0.01), evaluated per acquisition block (112 block results):

  | spread of the block's FMO tube | block results | worst FMO-gated item (error / band) | results with an item out of band |
  |---|---|---|---|
  | <= 50 | 25 | 0.58 | 0 |
  | 50-100 | 47 | 1.13 | 1 |
  | 100-150 | 24 | 0.96 | 0 |
  | > 150 (up to 248) | 16 | 1.56 | 1 |

  A first version of the criterion (150 units) was not sufficient: an error of 0.0025 in
  PE -> YG780-A gave a spread of 78 units and moved `activated_CD4_T` of one seed-2 sample by 0.76
  points (worklog L4-5). That is why the criterion is 50 units and the absolute band 0.75 points.
- Necessity: every default matrix misses the criterion by a factor 16 or more and fails. Between 50
  and about 150 units the result is usually still in band, so the criterion is conservative, not
  sharp.
- The check repairs what it finds (host, eight seeds, not part of the battery): bead-lot row plus
  fine-tuning on the FMO tubes scores 0.959-1.000 (the two PerCP-Cy5.5 entries that the FMO tube
  cannot see stay 0.008-0.026 off); row-only plus fine-tuning 1.000; both together 0.959-1.000.
  The correct matrices without fine-tuning also score 1.000 on eight seeds (worst FMO-gated item
  0.35 band). The `careless_*` variants of the battery therefore compute the SOP's check and do not
  act on it.

### Is the truth recoverable from the inputs? Three routes to the later matrix

Gain factor of YG780-A against the planted value, host, eight seeds: ratio of the
background-subtracted positive-bead medians of the two PE-Cy7 tubes (reference pipeline) -0.41% to
+0.28%; ratio of the two PE-Cy7 spillover rows -0.35% to +0.26%; ratio of the positive-bead peaks
from a two-component mixture on the log signal (second solution) -0.40% to +0.48%. Third route
without any bead tube: an event-level regression of the raw YG780-A signal of the re-recorded FMO
tube on the unmixed amounts of the other seven fluorochromes gives the seven column entries within
0.0077 of the planted block-2 values (seeds 1 and 2; worst entry FITC, +0.006 and +0.008, because
monocytes are more autofluorescent; PerCP-Cy5.5 -0.003 and -0.002).

### The two solutions at level 4

| step | reference `flow_pipeline_l4.py` | independent `solve_alt_l4.py` |
|---|---|---|
| reader, voltages | flowio keywords | fcsparser metadata |
| control choice | lot of `panel.tsv`, then session-start settings | lot filter, then the one settings group that has a control for every fluorochrome |
| spillover rows | median(positive) - median(negative) | least-squares slope over all events |
| gain factor | ratio of background-subtracted bead medians | ratio of mixture-model bead peaks |
| later matrix | column x gain, PE-Cy7 row from the re-recorded tube | same construction, own estimates |
| gates | density valleys | Gaussian mixtures |
| fine-tuning on the FMO tube | weighted least squares across seven populations (medians), five coefficients | per-population median matching against DN_T, four coefficients |
| FMO percentile | numpy, interpolated | nearest rank |

### What was dropped from the spec, and why

- **The rare double-positive population of `Q_dim_and_rare` was not built.** The three P elements
  act on YG780-A; none of them moves a lineage population (worst lineage item under the bead-lot and
  row-only variants: 0.01-0.03 band on eight seeds). A rare CD4/CD8 double-positive population that
  under-compensation inflates would have needed a fourth mis-compensation path into the CD4 or CD8
  detector, with its own control and its own careless variant. I preferred to spend the time on
  measuring the acceptance criterion. DP_T (0.8-2.5% of T cells) is still scored as a lineage
  population.
- **One file per acquisition block** became one matrix per sample in one file (see above).
- **Scoring weights** differ from level 3 (see above).
- `P_voltage_change` follows the spec with one addition: the FMO tube is re-recorded as well,
  otherwise the FMO gate of the later samples would have no control.

### Known limitations and doubts (level 4)

- **One detector, one population.** All four elements meet in YG780-A, and the population that fails
  reliably under the bead-lot and row-only variants is `activated_CD8_T`; `activated_CD4_T` stays in
  band under the bead-lot variant (48 of 48). Those two variants fail through the matrix component
  (0.20) and through `activation_every_sample` (0.25, minimum over samples); the lineage components
  give them full credit.
- **The SOP's check is a way around the reasoning.** An agent that runs the check and tunes the
  YG780-A column on the FMO tubes passes without understanding reagent lots or detector gain
  (measured, see above). I think this is legitimate, because it is what the SOP asks for and what a
  cytometrist does, but it means the difficulty of level 4 is: gate the FMO tube like a sample,
  compute the check, and be able to correct a failing compensation. Whether a small model does that
  is not known; no agent was run.
- **The embedded matrix is not subtle.** It wrecks the lineage gates (hardly any CD8_T). It closes the
  "use the keyword" shortcut but will not mislead anyone who looks at one plot.
- **The criterion is in absolute detector units.** The later block has 1.8-2.6 times the gain, so 50
  units is about twice as strict there. A criterion relative to the FMO threshold would be
  scale-free; I kept units because an SOP would state it that way.
- **Sensitivity of the FMO-gated populations is a property of the simulator.** The worst case of the
  sweep (0.97 points per 100 units of spread) comes from a sample with 4.0% activated CD4 T cells and
  a staining factor different from the FMO tube (tubes differ by 8% log-sd). Real data have other
  sources of mismatch between an FMO tube and a sample.
- **The two solutions fine-tune in related ways** (both on population medians). An event-level
  regression would have been a more different second method, but it is biased by the measurement
  noise of the regressor (-0.0020 for PerCP-Cy5.5 with the exact matrix, worklog L4-6). The other
  steps are independent.
- **Reference solution: 17 steps, 2 s.** Below the 20 steps AUTHORING asks for, recorded as 17. There
  is no element C or I; the task is harder by reasoning and self-checking, not by length. If the
  pilot shows that this is still short, a convergence element would have to be added.
- **Strictness.** One lineage item out of band in every sample: 0.858; one FMO-gated item out in one
  sample: 0.742. Both fail.
- **Gain from the voltages.** An agent that derives the gain from the voltage ratio with a textbook
  exponent is off by an amount I did not measure; the data (bead tubes, FMO) give it within 0.5%.
- **Realism.** A lot difference of 0.10-0.14 in one spillover coefficient is large for a real tandem
  lot; PerCP-Cy5.5 in a 561-nm 780 detector is inherited from the level-3 panel; dead cells and
  aggregates as at level 3.
- **Not validated:** Harbor (G6) for seed 2; levels other than 4 were only checked for byte
  identity, not re-validated; the instruction has been read by no agent.

### How the level-4 numbers were obtained

```bash
cd /home/user/ale-lifesci-gym/pipeline && export PATH="$HOME/.local/bin:$PATH"
../.venv/bin/python -m alegen make flow_compensation --seed 1 --level 4 --keep-work      # with Harbor
../.venv/bin/python -m alegen make flow_compensation --seed 2 --level 4 --no-harbor --keep-work
```
Result of the shipped version: seed 1 all 34 checks ok including Harbor oracle 1.0 / nop 0.0
(6.75 min); seed 2 all 32 checks ok, Harbor not run (3.01 min). Gate lines: `results/flow_s1_l4.log`,
`results/flow_s2_l4.log`; per-check details: `results/validation/flow-compensation-s000{1,2}-l4.json`.
After this validation two sentences of the task card were corrected in `l4.py`: one number in the
judge text (0.61 -> 0.58, the sweep value for a spread <= 50) and the evidence line of
`P_embedded_matrix` (the check is not evaluable with the embedded matrix, rather than "fails"). Both
task directories were regenerated with `--no-validate`; `diff -rq` against copies of the validated
directories shows `README.md` as the only file that differs, so the images, the instruction, the
inputs, the solutions and the verifier are the validated bytes. The gates were not re-run for this.
Three earlier validation rounds of seed 1 (all gates ok) were run on intermediate versions and are
not the shipped task (worklog L4-3).

Raw metrics: `info` block of `results/work/flow-compensation-s000{1,2}-l4/{oracle,probe/*,attack/*}/logs/details.json`;
the acceptance check of every pipeline run is in `output/qc_fmo.tsv` next to it. Host repeats,
sweep and rescue runs used the same scripts with `--input/--output` (flags `--refine off|on`, and
an environment variable that mis-sets one coefficient in a scratch copy of the pipeline); their
logs were in the session scratch directory and are not kept.

---

## Level 7 (added 2026-11-12)

Spec: `docs/specs/l7/flow_compensation.md`. Code: `l7.py` (generator, battery, task card),
`assets/flow_pipeline_l7.py` (reference pipeline and alternative treatments of the disturbed stretch),
`assets/solve_alt_l7.py` (independent second solution), `assets/solve_l7.sh`. Scorer:
`assets/score_l4.py`, unchanged. `l4.py` got one optional argument (`disturb`: a function that
returns the events written to a sample file; `None` at levels 4-6) and `__init__.py` sends
`level == 7` to `l7.py`. The task is built by `python -m alegen make flow_compensation --level 7`:
the family returns the level-7 build with the level-4 texts, and `descaffold.apply(build, 7)` takes
the acceptance criterion and the incident reports out exactly as for level 6.

All numbers in this section were measured on 2026-11-12. "Docker" = inside the task images during
`python -m alegen make`. "Host replay" = the same battery scripts run with `--input/--output` in a
container of the level-4 agent image (same package versions) and scored with the task's own
`score.py`; Docker and host replay give the same score, to six decimals, for all 30 non-generic checks
of seed 1 and of seed 2.

### What was planted

One decision point, `T_acquisition_disturbance`. Two of the six samples, one from each voltage block,
contain a contiguous stretch of events (contiguous in `Time`) that was recorded during a fluidics
disturbance:

- in the two detectors of the 561-nm laser, YG586-A (CD4 PE) and YG780-A (HLA-DR PE-Cy7), the file
  holds the recorded value times one factor per event: attenuation x exp(jitter x z). Attenuation
  U(0.30, 0.45) and jitter U(0.01, 0.02) are drawn per sample, z is standard normal and the same for
  both detectors of an event;
- FSC-A, FSC-H, SSC-A and the six detectors of the other three lasers are untouched;
- the time stamps inside the stretch are spread so that the event rate is U(0.90, 0.97) of the rate
  before it;
- the stretch holds U(12%, 25%) of the events, starts at U(25%, 60%) of the run and ends before 85%.

Everything else is the level-4 session of the same seed. The generator simulates the session twice
with the level-4 code and the level-4 random stream; the disturbance has its own stream
(`default_rng(71000 + seed)`). Measured for seed 1: against `tasks/flow-compensation-s0001-l6`, the
level-7 task differs in two input files (`PBMC_D735.fcs`, `PBMC_D960.fcs`: identical keywords, the
YG586-A / YG780-A columns differ in the planted rows only, `Time` from the first planted row on), in
`solution/`, in the task card and in the `plants` block of `truth.json`. `instruction.md`,
`gating_sop.md`, `acquisition_notes.md`, `score.py`, the Dockerfiles, the other 17 FCS files and the
scored part of `truth.json` are identical. Input size 19.84 MB.

| seed | sample (block) | planted events | share | start (fraction of run) | Time | signal | jitter | event rate |
|---|---|---|---|---|---|---|---|---|
| 1 | PBMC_D735 (1) | 22784-29651 of 38220 | 18.0% | 0.596 | 7104-9391 | x0.356 | 1.9% | x0.927 |
| 1 | PBMC_D960 (2) | 19748-25472 of 36413 | 15.7% | 0.542 | 6177-8154 | x0.304 | 1.7% | x0.919 |
| 2 | PBMC_D124 (1) | 24065-30707 of 41059 | 16.2% | 0.586 | 7523-9791 | x0.412 | 1.5% | x0.922 |
| 2 | PBMC_D523 (2) | 20870-27714 of 36071 | 19.0% | 0.579 | 6523-8733 | x0.354 | 1.7% | x0.967 |

In seed 1 the sample with the lowered FSC voltage (PBMC_D858) is not one of the two; in seed 2
neither (PBMC_D507). Nothing forbids the coincidence.

### Truth

The reference of a level-7 task is the reference of the level-4 task of the same seed, number for
number (checked for seed 1: `truth.json` without `plants.acquisition_disturbance` equals that of
level 6). `generate_l4` computes it from the recording as simulated; what `disturb` returns is only
written to the file. Membership in a lineage population is the hidden label of the event, inside the
stretch or not; the FMO-gated populations are the SOP rule applied to the signals as they would have
been recorded. The existing definitions already say this: `count` is the "number of events of the
sample's FCS file that belong to the population", and a CD4 T cell does not stop being one because a
detector read low while it passed the laser. No sentence was added to any text.

An analyst cannot see the label of a disturbed event; what is available is an unbiased estimate from
the undisturbed events (or a correction of the signal). With a tolerance of max(0.5 points, 3%) that
estimate is not automatically in band: over 2000 random positions in the six samples of seed 1, the
reference percentages of the events *before* the stretch were within the band for all 14 populations
in 40-91% of the positions, and within half the band in 0-40% (one sample: none). The generator
therefore chooses the position: per sample, the first of up to 600 seeded candidates for which the
reference percentages of the remaining events are within 0.5 band for all 14 populations, both when
the stretch alone is removed and when everything from its start to the end of the file is removed;
per block, one sample at random among those that have such a position. If no sample of a block has
one, the bound is relaxed to 0.6, then 0.7 (recorded in `truth.json`); beyond that the generator
raises. Scan of seeds 1-48: 90 of the 94 planted stretches met 0.5, two 0.6, two 0.7, and seed 32
cannot be generated (best position of block 1: 0.75).

### What goes wrong without the step

The level-4 reference pipeline runs through and reports numbers of the usual shape. Two things are
wrong in each of the two samples (12 seeds, 24 samples, host replay; Docker for seeds 1 and 2):

- The compensation subtracts PE spill from B710-A in proportion to the YG586-A signal. Inside the
  stretch YG586-A is a third of what it should be, so two thirds of the PE spill stay in the CD8
  detector (planted PE -> B710-A 0.18-0.32, CD4 at 30 000: 3000-6700 units against 8000 for real CD8)
  and the CD4 T cells of the stretch are counted as CD4+CD8+. `CD4_T` is 7.9-14.4 points low (4.4-8.0
  bands), `DP_T` 7.9-14.4 points high (15.8-28.8 bands; reported 9.3-16.5% of T cells, reference
  0.8-2.7%).
- The compensation subtracts the full PerCP-Cy5.5 spill from a YG780-A signal that was recorded too
  low, and the surplus in B710-A is carried on into YG780-A. The compensated HLA-DR signal of the T
  cells of the stretch has a median of -1100 to -1300 units in a block-1 sample and -2100 to -2500 in
  three block-2 samples (seeds 1, 4, 12), so 0.0-0.7% of them pass the FMO gate instead of 4-22%.
  `activated_CD8_T` is 1.0-5.8 points low (1.3-4.1 bands).
  `activated_CD4_T` stays in band because its parent population has lost the cells of the stretch to
  `DP_T`.

Six of the 84 scored percentages leave the band, the same six on all twelve seeds; the six scored
lymphocyte counts are unaffected (error at most 0.04%). Score 0.472 on every seed
(`lineage_every_sample` 0, `activation_every_sample` 0). Each effect alone would also fail, but
narrowly: lineage alone 0.739, `activated_CD8_T` alone 0.733 (scorer arithmetic, not a measured
submission).

### How it is discovered

No text says anything; `Time` is listed among the recorded parameters, as at every level. The
evidence is the signal along `Time`, the first check of flowAI or PeacoQC and what a Time gate is
for. As computed by the generator on the files as written (median of each fluorescence detector in
consecutive bins of 500 events ordered by `Time`; SD = robust standard deviation of the bin medians
outside the stretch):

| seed | sample | bins inside / in file | YG780-A bin median | YG586-A bin median | YG586-A 90th percentile per bin | other six detectors |
|---|---|---|---|---|---|---|
| 1 | PBMC_D735 | 12 / 76 | 1255 -> 447 (x0.36, 10 SD) | 1511 -> 536 (x0.35, 9 SD) | 33436 -> 11980 (23 SD) | at most 0.5 SD |
| 1 | PBMC_D960 | 10 / 72 | 2487 -> 768 (x0.31, 12 SD) | 1507 -> 447 (x0.30, 12 SD) | 36176 -> 11039 (26 SD) | at most 0.4 SD |
| 2 | PBMC_D124 | 12 / 82 | 1618 -> 680 (x0.42, 15 SD) | 1790 -> 767 (x0.43, 7 SD) | 33115 -> 13502 (19 SD) | at most 0.3 SD |
| 2 | PBMC_D523 | 13 / 72 | 3898 -> 1411 (x0.36, 16 SD) | 1794 -> 619 (x0.35, 9 SD) | 33406 -> 11626 (29 SD) | at most 0.4 SD |

The generator refuses an instance unless, in at least one of the two detectors, every bin inside the
stretch and no bin outside it has a median nearer to the inside level than to the outside level.
Scan of seeds 1-48 (94 stretches): YG780-A meets this in 94 of 94 and its ratio of bin medians is
within 5% of the planted attenuation; YG586-A meets it in 35 of 94. The median of YG586-A sits between
populations (CD4-negative lymphocytes below, monocytes and dead cells around it), so its trace is
noisy and up to 25 of about 60 outside bins come as low as the stretch; its 90th percentile, which
lies in the CD4 peak, is clean. The stretch is 8-20 bins long (72-87 bins per file). It also shows
without any time plot as a third population in a CD4-versus-CD8 plot of T cells that exists in one
Time range only.

The image has no ready-made tool for this: flowAI and PeacoQC are R packages, and FlowKit 1.3.2
offers `Sample.plot_channel` (a channel against the event index) but no automatic check
(`filter_anomalous_events` does not exist in this version; tried). The check is fifteen lines of
numpy.

### The two solutions and the three treatments

| step | reference `flow_pipeline_l7.py` | independent `solve_alt_l7.py` |
|---|---|---|
| statistic that finds the stretch | median and 90th percentile of every detector and of FSC-A / SSC-A in bins of 500 events; a run of at least 5 bins beyond 4 robust SD (arcsinh scale) with at most 10% such bins elsewhere | no bins: ranks of every parameter in Time order, cumulative sum of the centred ranks, stretch = highest to lowest point of the curve, z = rank sum / sqrt(length x (1 - length / n)), accepted from z = 10 |
| which detectors | factor inside / outside per parameter (median over the 50th-95th percentiles), affected if it differs from 1 by more than 15% | every parameter with z >= 10 |
| edges | change points of an event-wise score (log ratio of the histograms of the affected detectors outside and inside) | the two extremes of the cumulative sum; union over the parameters; 150 events margin |
| treatment | events of the stretch (and 250 on either side) are not used for thresholds; populations whose gates need neither an affected detector nor one that is compensated for it by 2% or more are counted on all events; the others are estimated as share of the parent among the undisturbed events x size of the parent in the file | the stretch is cut out of the tube before anything else; all percentages from the remaining events; every count scaled by events in the file / events used |
| rest | as at level 4 | as at level 4 |

Measured (12 seeds, 24 stretches, host replay). Reference pipeline: first and last event within 7
events of the planted ones (21 of 24 stretches within 2); factors within 2.2% of the planted attenuation
(Docker: 0.3555 / 0.3551 for 0.3560, 0.3048 / 0.3044 for 0.3043, 0.4122 / 0.4112 for 0.4122, 0.3512 /
0.3577 for 0.3544); no stretch reported in any of the 72 undisturbed sample and FMO tubes or in the
cell-based single stain. Second solution: z = 27.6-51.7 for the two affected detectors of the 24
disturbed tubes, at most 4.8 for any parameter of the 72 undisturbed tubes; its stretch is up to 523
events too long on one side (seed 9, PBMC_D415: 17919-27230 for 17935-26707), which costs clean
events and nothing else, and at most 5 events too short (seed 2, PBMC_D523), which its margin of 150
events covers. With the reference rule, `B_cells` and `NK_cells` are estimated from the
undisturbed events as well (PE-Cy7 -> R780-A is 0.06-0.15); `cells` to `T_cells` are counted on all
events.

Three treatments of a stretch reach the reference and are probes: set aside (the two solutions), set
aside everything from its start to the end of the file, and rescale the affected detectors inside it.

### Battery (pass threshold 0.9; ceiling of the "ignores one element" attacks 0.75)

Columns per seed: score | worst lineage item (error / band) | worst FMO-gated item (error / band).
Docker, seeds 1 and 2; last column: range of the score in the host replay of seeds 1-12.

| check | bound | s1 score | lineage | FMO-gated | s2 score | lineage | FMO-gated | seeds 1-12 |
|---|---|---|---|---|---|---|---|---|
| oracle (reference pipeline) | >= 0.99 | 1.000 | 0.24 | 0.41 | 1.000 | 0.35 | 0.26 | 1.000 |
| planted_truth | >= 0.999 | 1.000 | 0.00 | 0.00 | 1.000 | 0.00 | 0.00 | 1.000 |
| within_tolerance | >= 0.99 | 1.000 | 0.60 | 0.60 | 1.000 | 0.60 | 0.60 | 1.000 |
| independent_solution | >= 0.99 | 1.000 | 0.25 | 0.33 | 1.000 | 0.30 | 0.22 | 1.000 |
| time_gate_everything_after | >= 0.99 | 1.000 | 0.41 | 0.57 | 1.000 | 0.38 | 0.35 | 1.000 |
| repair_instead_of_exclude | >= 0.99 | 1.000 | 0.02 | 0.12 | 1.000 | 0.01 | 0.10 | 1.000 |
| time_gate_counts_of_kept_events | >= 0.90 | 0.933 | 0.24 | 0.40 | 0.933 | 0.35 | 0.25 | 0.933 |
| **l6_reference** | <= 0.75 | 0.472 | 23.46 | 3.12 | 0.472 | 23.85 | 3.38 | 0.472 |
| l6_second_solution | <= 0.75 | 0.472 | 23.47 | 3.18 | 0.472 | 23.85 | 3.38 | 0.472 |
| exclude_wrong_window | <= 0.75 | 0.472 | 29.15 | 3.77 | 0.472 | 29.87 | 3.86 | 0.472 |
| l3_reference | <= 0.75 | 0.189 | 23.46 | 59.49 | 0.189 | 23.85 | 88.04 | 0.186-0.189 |
| l3_reference_last_control | <= 0.75 | 0.189 | 23.46 | 16.67 | 0.189 | 23.85 | 34.33 | 0.189 |
| careless_embedded_matrix | <= 0.75 | 0.000 | 93.14 | 16.67 | 0.175 | 78.58 | 16.67 | 0.000-0.317 |
| careless_one_matrix | <= 0.75 | 0.517 | 0.22 | 123.68 | 0.483 | 0.35 | 127.55 | 0.372-0.517 |
| careless_row_only | <= 0.75 | 0.533 | 0.24 | 8.12 | 0.483 | 0.35 | 7.78 | 0.483-0.533 |
| careless_bead_tandem | <= 0.75 | 0.500 | 0.24 | 8.71 | 0.500 | 0.35 | 8.27 | 0.500 |
| careless_no_fmo | <= 0.75 | 0.650 | 0.24 | 42.88 | 0.650 | 0.35 | 47.50 | 0.650 |
| careless_fmo_uncompensated | <= 0.75 | 0.650 | 0.24 | 14.36 | 0.650 | 0.35 | 11.37 | 0.650-0.667 |
| textbook_prior | <= 0.10 | 0.000 | 12.40 | 5.96 | 0.008 | 19.78 | 7.34 | 0.000-0.033 |
| block1_matrix_true_frequencies | <= 0.85 | 0.800 | 0.00 | 0.00 | 0.800 | 0.00 | 0.00 | 0.800 |
| all_percent_beyond_tolerance | <= 0.45 | 0.300 | 1.40 | 1.40 | 0.300 | 1.40 | 1.40 | 0.300 |
| one_population_beyond_tolerance | <= 0.89 | 0.858 | 1.50 | 0.00 | 0.858 | 1.50 | 0.00 | 0.858 |
| activation_beyond_tolerance | <= 0.70 | 0.650 | 0.00 | 1.50 | 0.650 | 0.00 | 1.50 | 0.650 |
| one_activation_item_beyond_tolerance | <= 0.80 | 0.742 | 0.00 | 1.50 | 0.742 | 0.00 | 1.50 | 0.742 |
| matrix_entry_beyond_tolerance | <= 0.85 | 0.800 | 0.00 | 0.00 | 0.800 | 0.00 | 0.00 | 0.800 |
| naive_matrix_not_applied | <= 0.45 | 0.200 | 166.67 | 16.67 | 0.200 | 166.70 | 16.67 | 0.200 |
| naive_no_doublet_gate | <= 0.45 | 0.281 | 4.20 | 2.99 | 0.258 | 4.78 | 3.11 | 0.236-0.281 |
| naive_no_viability_gate | <= 0.45 | 0.267 | 46.88 | 29.04 | 0.267 | 50.34 | 23.21 | 0.250-0.283 |
| naive_global_thresholds | <= 0.60 | 0.417 | 120.01 | 16.67 | 0.417 | 33.33 | 15.78 | 0.417 |
| naive_default_pipeline | <= 0.20 | 0.008 | 61.05 | 16.67 | 0.000 | 52.17 | 32.87 | 0.000-0.008 |

The four generic attacks of the validator score 0 at a hard gate on both seeds. Further numbers of
the level-7 checks over seeds 1-12 (host replay), worst item as a fraction of its band, lineage /
FMO-gated: reference pipeline 0.48 / 0.41, second solution 0.52 / 0.48, everything after the start
set aside 0.64 / 0.57, rescaled 0.02 / 0.37. Largest relative error of the lymphocyte count: 0.04%
(reference pipeline), 0.26% (second solution, scaled counts), 19-26% (counts of the kept events).

How the battery was assembled:
- `l6_reference`, `l6_second_solution`, `l3_reference`, `l3_reference_last_control` run the level-4
  scripts as they are (`flow_pipeline_l4.py`, `solve_alt_l4.py`), i.e. without the new step.
  `l3_reference` scored 0.450 at level 4; here it also loses the two lineage items.
- The eleven "ignores exactly one element" attacks of level 4 (`careless_*`, `naive_*`) run
  `flow_pipeline_l7.py` with their level-4 flags, i.e. with the signal-versus-time check done. Ten of
  them score exactly what they scored at level 4 on both seeds; `naive_no_doublet_gate` is 0.003 /
  0.006 higher (level 4: 0.278 / 0.253). Handling the new step does not rescue an old mistake.
  Ceilings unchanged.
- The truth-derived probes and attacks are the level-4 files (same reference).
- `time_gate_counts_of_kept_events` has the pass threshold as floor instead of 0.99, on purpose: the
  percentages are right and two of six lymphocyte counts are 19-26% low, which costs two thirds of
  the `lymphocyte_counts` component. It documents that the right QC with this reading of `count`
  passes (0.933), not that it is flawless.

### Scorer and judge text

No change to the scorer. The spec asks that the scored counts be countable on all events: the only
scored count is `lymphocytes`, defined by scatter, the viability dye (405 nm) and CD14 (488 nm).
Planted spill of PE and PE-Cy7 into B530-A and V525-A is at most 0.015, so the residual in those
detectors inside the stretch stays far below the gates; both solutions report the count within 0.3%.
The judge explanation says so, describes the three accepted treatments and gives the numbers of the
naive analysis.

### Levels 1-6 after the change

`diff -r` of seed-1 tasks generated with the family code at git HEAD (before the change) and with the
changed code, levels 1, 2, 3, 4, 5, 6: no difference. Against the task directories in `tasks/`:
levels 3 and 4 identical; levels 5 and 6 differ in `README.md` and `task.toml` (title and summary
wording, `source_level` line). That difference is there with the unchanged family code as well: it
comes from changes to `alegen/descaffold.py` made after those two directories were written. Every
other file of levels 5 and 6 is identical.

### Known limitations and doubts (level 7)

- **The disturbance is a multiplication of the recorded value, with a jitter of 1-2%.** A physical
  loss of light would leave fewer photoelectrons (relatively more photon noise) and an unscaled
  electronic noise, and rescaling would then not restore the background width of YG780-A. I chose the
  multiplicative model so that the repair treatment can be exact. The spec's "with extra noise" is
  therefore implemented weakly. With a jitter of 2-4% (first version) the repair probe failed on one
  of twelve seeds (`activated_CD4_T` 1.8 bands, seed 4) and was at 0.79 and 0.87 band on two more: a
  3.5% event-to-event error of the repaired PE signal leaves about +-250 units in B710-A, which the
  compensation carries into YG780-A with the PerCP-Cy5.5 coefficient (0.44-0.46 in block 2 of those
  seeds) and widens the HLA-DR background of the CD4 T cells of the stretch.
- **Repair needs the factor to about 1%.** The first estimator (median of percentile ratios) was off
  by up to 5% and failed the probe; the present one (shift between the peaks of the log signal) is
  within 2.2%. An agent that rescales with a rough factor can miss an FMO-gated population.
  Exclusion is the robust route; repair is accepted, not easy.
- **The position of the stretch is not free.** The condition on the remaining events favours late
  starts: 38 of 94 stretches (seeds 1-48) start after 55% of the run, where a uniform draw would put
  14%; 15 start before 40%. One seed in 48 has no instance. With the relaxed bound (4 of 94
  stretches) the margin of `time_gate_everything_after` is thin: host replay of seeds 13, 20, 33, 48
  passes with the worst item at 0.79, 0.78, 0.58, 0.55 band.
- **The naive tables are not subtle to an immunologist.** 9-17% CD4+CD8+ T cells in blood is a number
  a specialist would question. Nothing fails and a program sees tables of the usual shape, which is
  what the spec asks for, but "plausible for PBMC" holds for CD4_T and activated_CD8_T and only
  loosely for DP_T. A smaller lineage effect would cost the margin of `l6_reference`.
- **`activated_CD8_T` is out by as little as 1.3 bands** (a sample with 8.8% reference). The 0.472
  rests on both effects; see the two "alone" scores above.
- **The median of the CD4 detector is a poor indicator** in 59 of 94 stretches (see above). An agent
  that plots only YG586-A medians sees a noisy trace with a dip; YG780-A, an upper percentile or a
  rank statistic show a clean step.
- **Easy once looked at:** one stretch per file, sharp edges, the same factor in both detectors, the
  other detectors exactly untouched, undisturbed files exactly stationary (events are exchangeable in
  time by construction, so any drift test has an exact null). Real files drift.
- **Counts of kept events pass (0.933).** With four of six lymphocyte counts right the component
  still gives a third of its weight. I did not change the scorer (spec: the scorer stays).
- **The 2% rule** of the reference pipeline for "counted on all events" is a rule of thumb; I did not
  measure what counting `B_cells` and `NK_cells` on all events would give.
- **The two solutions share the level-4 gating code paths** (valleys against mixtures, as at level 4)
  and differ in the time check and in the treatment.
- **Reference solution: 21 steps** as listed in its docstring, 2-3 s; Docker wall time including container
  start 0.04 min (seed 1, both rounds) and 0.05 / 0.14 min (seed 2, two rounds, host load 10-17).
- **Not validated:** Harbor (G6) for seed 2; generation only for seeds 13-48 (battery replay for
  seeds 1-12, level-7 checks for 13, 20, 33, 48); no agent has been run on a level-7 task.

### How the level-7 numbers were obtained

```bash
cd /home/user/ale-lifesci-gym/pipeline && export PATH="$HOME/.local/bin:$PATH"
../.venv/bin/python -m alegen make flow_compensation --seed 1 --level 7 --keep-work               # with Harbor
../.venv/bin/python -m alegen make flow_compensation --seed 2 --level 7 --no-harbor --keep-work
```
Result of the shipped version: seed 1 all 41 checks ok including Harbor oracle 1.0 / nop 0.0 (4.54 min),
seed 2 all 39 checks ok, Harbor not run (2.79 min). Per-check records:
`results/validation/flow-compensation-s000{1,2}-l7.json`; raw metrics: `info` block of
`results/work/flow-compensation-s000{1,2}-l7/{oracle,probe/*,attack/*}/logs/details.json`; what the
reference pipeline found in every run: `output/qc_time.tsv` next to it. A first round of both
validations (all checks ok, same scores) was run before two sentences of the judge text were
corrected (ranges of the naive deviations per item); `diff -rq` of the task directories of the two
rounds shows `README.md` as the only file that differs, and both seeds were validated again after
the correction.

Host replay (seeds 1-12 all 30 checks, seeds 13 / 20 / 33 / 48 the level-7 checks), the generation
scan of seeds 1-48, the 2000-position survey and the first versions that failed (worklog L7-3 to
L7-6) used scripts in the session scratch directory; they are not kept. The byte-identity check of
levels 1-6 is described in the worklog (L7-8).
