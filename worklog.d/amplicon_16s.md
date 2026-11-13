# Worklog: family `amplicon_16s`

All entries 2026-11-11 (UTC), single authoring session.

## 15:30 Reading
AUTHORING.md, docs/specs/amplicon_16s.md, the crispr_screen reference family, core modules.

## 15:35 Reference material
- Downloaded `https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz`
  to `.cache/amplicon_16s/` (sha256 `c5dc51b877cfed87b526b5cb9425f05d59ce96805ce4d4eb2015d65d75436e7e`,
  26,888 records).
- Exploration (`515F GTGYCAGCMGCCGCGGTAA`, `806R GGACTACNVGGGTWTCTAAT`, exact IUPAC match):
  24,344 records give an amplicon, 21,472 of them 253 nt; 13,679 unique amplicons, 406 shared by
  more than one genus (dropped); 2,316 same-genus pairs at Hamming distance 1 and 3,483 at distance 2.
  Enough real material for strain pairs, no need to invent them.
- Built the tools image `alegen-tools-1d9afa22bd69:local` (vsearch=2.30.1, cutadapt=5.2; 1.0 GB).
  Note: the session scratchpad is shared between authors (another author's `toolimg.log` overwrote
  mine), so everything of mine lives in `scratchpad/amp16s/`.

## 15:40 Design decisions (before code)
- Truth = planted sequences + per-library read-origin counts; recombinant pairs count for nobody
  unless the recombinant is identical to a planted sequence.
- Reference set shipped V4-trimmed (like the 515F/806R-trimmed SILVA classifiers), ids blinded
  (`REFnnnnn`), planted sources hidden among thousands of distractors; 6 planted strains carry 1-2
  substitutions relative to their reference (not in the set), 2 have every close reference removed
  (expected genus `Unclassified`).
- Genus truth without an aligner: 8-mer pigeonhole bound (an edit destroys at most 8 k-mers), used
  to remove confusable references. First idea was vsearch in a tool image at generation time;
  dropped because the bound is provable, needs no Docker and is independent of the agent's tool.
- Reagent contaminants are picked at random among the planted strains, not from the usual
  suspects (Ralstonia, ...), so prior knowledge of genus names does not help.
- Counts are scored as counts (read pairs that derive from the ASV), not as relative abundance:
  otherwise a workflow that silently loses one read orientation is indistinguishable from a
  correct one. The instruction defines the count accordingly.
- Contaminants.txt carries a header line so that "no contaminant" is not an empty file (the core's
  generic `zero_byte_files` attack must stay at 0).

## 15:47 First generation + oracle (seed 1, level 3)
- Generation 6-9 s, inputs 13 MB.
- Oracle (UNOISE3 + uchime3_denovo): 0.984. Three non-planted sequences survived, all exact
  two-parent recombinants with a short (37-38 nt) tail from the second parent, 4-5 nt from the main
  parent. One genus mismatch.
- Bug found through the genus mismatch: `protected = set(db)` was taken after the distractors had
  been added, so no reference was ever removed (`db_removed: 0`) and the `Unclassified` plant had a
  close reference. Fixed (protected = planted sources only). After the fix 2,383 of 4,000
  distractors were removed (the bound for the `Unclassified` plants is conservative), so the
  number of distractors drawn was raised to 9,000 (3.8k references remain for seed 1).
- Chimera filters compared on the 163 UNOISE3 centroids of seed 1 (45 planted + 118 recombinants):

  | filter | kept | planted kept | non-planted kept |
  |---|---|---|---|
  | none | 163 | 45 | 118 |
  | `--uchime_denovo` | 45 | 45 | 0 |
  | `--uchime_denovo --mindiffs 1 --mindiv 0.3 --minh 0.1` | 45 | 45 | 0 |
  | `--uchime2_denovo` | 48 | 45 | 3 |
  | `--uchime3_denovo` | 48 | 45 | 3 |
  | `--uchime3_denovo --abskew 2` | 48 | 45 | 3 |
  | `--uchime3_denovo --mindiffs 1 --mindiv 0.3` | 48 | 45 | 3 |
  | `--uchime3_denovo --abskew 2 --mindiffs 1 --mindiv 0.3 --minh 0.1` | 48 | 45 | 3 |

  Oracle switched to `--uchime_denovo`. The precision band was moved from [0.70, 0.95] to
  [0.60, 0.90]: the stock UNOISE3 + uchime3 route is a legitimate workflow whose leftover
  recombinants should cost little, while "no chimera removal" sits at precision 0.28.

## 15:53 Second solution and full dev battery (tool image, 2 threads, verifier run on the host)
- `solve_alt.py` written: merge first, primer search in python on the merged read (both strands),
  per-library Poisson abundance test against a data-estimated substitution rate, exact two-parent
  bimera test, numpy Hamming read assignment, python k-mer + banded-edit-distance taxonomy,
  contaminants by constancy of relative abundance x qPCR copies. Seed 1: 1.000 on the first run
  (45 sequences, all planted; the oracle reported 44, both miss only gray-zone strains).
- The read-mapping step of the reference pipeline (`usearch_global --maxaccepts 16 --maxrejects 64`)
  took 146 s for the "no primer removal" variant (every read is a non-exact query). Reduced to
  `--maxaccepts 4 --maxrejects 8`; oracle counts unchanged to the third decimal.
- Dev results, level 3 (score | recall | precision | min study count agreement | mock | MCC | genus):

  | run | seed 1 | seed 2 |
  |---|---|---|
  | oracle | 1.000 / 1.0 / 1.0 / 0.977 / 0.991 / 1.0 / 1.0 | 1.000 / 1.0 / 1.0 / 0.979 / 0.988 / 1.0 / 1.0 |
  | solve_alt | 1.000 / 1.0 / 1.0 / 0.978 / 0.991 / 1.0 / 1.0 | 1.000 / 1.0 / 1.0 / 0.979 / 0.987 / 1.0 / 1.0 |
  | UNOISE3 + uchime3 | 0.977 / 1.0 / 0.88 / 0.975 / 0.991 / 1.0 / 1.0 | 1.000 / 1.0 / 1.0 / 0.979 / 0.988 / 1.0 / 1.0 |
  | no primer removal | 0.000 / 0 / 0 / 0 / 0 / 0 / 0 | 0.000 / 0 / 0 / 0 / 0 / 0 / 0 |
  | fixed-length trim | 0.000 / 0.07 / 0.04 / 0 / 0.10 / 0.43 / 0.07 | 0.000 / 0.05 / 0.02 / 0 / 0.07 / 0 / 0.05 |
  | single orientation, untrimmed kept | 0.250 / 1.0 / 0.47 / 0.449 / 0.439 / 1.0 / 1.0 | 0.600 / 1.0 / 0.98 / 0.512 / 0.525 / 1.0 / 1.0 |
  | single orientation, untrimmed discarded | 0.600 / 1.0 / 1.0 / 0.449 / 0.438 / 1.0 / 1.0 | 0.600 / 1.0 / 1.0 / 0.512 / 0.525 / 1.0 / 1.0 |
  | no chimera removal | 0.650 / 1.0 / 0.30 / 0.945 / 0.973 / 1.0 / 1.0 | 0.650 / 1.0 / 0.31 / 0.932 / 0.972 / 1.0 / 1.0 |
  | 97 % OTUs | 0.327 / 0.86 / 1.0 / 0.439 / 0.857 / 1.0 / 0.86 | 0.327 / 0.86 / 1.0 / 0.472 / 0.876 / 1.0 / 0.86 |
  | no contaminant handling | 0.550 / 1.0 / 1.0 / 0.407 / 0.978 / 0 / 1.0 | 0.550 / 1.0 / 1.0 / 0.482 / 0.975 / 0 / 1.0 |
  | everything in a blank is a contaminant | 0.450 / 1.0 / 1.0 / 0.020 / 0.268 / 0.16 / 1.0 | 0.450 / 1.0 / 1.0 / 0.025 / 0.000 / 0.12 / 1.0 |
  | fully naive | 0.040 / 0.86 / 0.17 / 0.178 / 0.425 / 0 / 0.86 | 0.040 / 0.86 / 0.13 / 0.135 / 0.504 / 0 / 0.86 |

  Observation: "single orientation, untrimmed kept" differs between the seeds (precision 0.47 vs
  0.98). The reverse-strand reads carry 16 spacer combinations x 96 primer instances, so their
  copies are diluted about 1,500-fold; only the most abundant strains of seed 1 reach the UNOISE
  minimum size as primer-carrying junk. The variant fails on the count component in both seeds.
- Level 1, seed 1: oracle 1.000, solve_alt 1.000 (count agreement 0.999).
- Generation probed for seeds 2, 3, 4 at level 3: no assertion failure; smallest scored strain
  135-246 read pairs, rare (gray) strains 8-27.

## 16:08 Battery additions before the gate runs
- `closed_reference_best_case` attack (truth-derived upper bound of reference-only picking).
- `stock_unoise3_uchime3` probe with floor = pass threshold (written reason in the code): the
  stock recipe must pass but is allowed to lose precision credit for leftover recombinants.

## 16:10-16:25 First gate runs
- `make --seed 1 --level 3 --keep-work`: all gates ok (6.6 min), Harbor oracle 1.0 / nop 0.0.
  Oracle count agreement after the `--maxaccepts 4 --maxrejects 8` change: 0.9768 (min study) and
  0.9907 (mock), identical to the values before the change.
- `make --seed 2 --level 3 --no-harbor`: all gates ok (7.1 min).
- `make --seed 1 --level 1 --no-harbor`: **G4 FAIL** `reference_as_asvs score=0.15 ceiling=0.1`.
  Cause: a real hole in the verifier, not an attack that was too strict. Without reagent background
  (levels 1-2) the contaminant component gave 1.0 whenever no planted strain was flagged, so a
  header-only `contaminants.txt` earned 0.15 for a submission that recovered nothing.
  Fix (ceiling untouched): in batches without reagent background the component is the fraction of
  scored planted sequences that are *reported and not flagged* (ramp 0.80 -> 0.95, same convention
  as the genus component: an unreported strain cannot earn credit). The level-3 branch (MCC) was
  already safe: MCC is 0 without a true positive.
  Rerun level 1: all gates ok, `reference_as_asvs` 0.0, `within_tolerance` 0.993.
- Because `score.py` changed, both level-3 runs were repeated with the final verifier (the gate
  lines in the report are from the repeated runs). `naive_everything_in_blanks` now runs at every
  level (blanks hold cross-talk at every level). A level-2 run was added.
- Checked that planted sources are not clustered in `reference_db.fasta` (seed 1: 38 planted among
  3,782 records at positions 18 ... 3723, uniformly spread); the 2 hits of `reference_as_asvs`
  among the first 40 records are chance.
- Measured separation of the two contaminant statistics on the planted counts: blank/other
  relative-abundance ratio 4.9-6.5 (reagent) vs <= 0.27 (others), seed 1; relative abundance x
  qPCR copies, blank/other, 0.70-1.04 (reagent) vs <= 0.004 (others), seeds 1-2.
- Footprint: validated task images are removed after each run (`docker rmi` of my own
  `alegen-amplicon-16s-*` tags), `results/work` kept only for seed 1 level 3; one tools image.

## 16:28-16:50 Final gate runs (final `score.py`, final battery)
- `make amplicon_16s --seed 1 --level 3 --keep-work`: all gates ok incl. Harbor (oracle 1.0, nop 0.0), 8.1 min.
- `make amplicon_16s --seed 2 --level 3 --no-harbor` (run with `--keep-work` added to read the raw
  metrics, work dir deleted afterwards): all gates ok, 7.0 min; metric table identical to the first run.
- `make amplicon_16s --seed 1 --level 1 --no-harbor` (also with `--keep-work`): all gates ok, 1.4 min.
- Extra: `make amplicon_16s --seed 1 --level 2 --no-harbor`: all gates ok, 5.5 min.
  Level 2 note: single orientation with discarded pairs reaches a count agreement of 0.5646
  (forward share 0.563), 0.035 below the zero point of the count band (0.60). With the allowed
  range of the orientation split (42-58 %) the worst case is about 0.58: still zero credit, but
  the margin is thin.
- Measured in the dev harness: exact-match-only counting (`--counts exact`) gives a min study count
  agreement of 0.917 / 0.914 (seeds 1 / 2), full credit.
- Measured on the dev outputs: without chimera removal the mock library shows 30 / 24 sequences
  with >= 8 read pairs (15 / 11 not planted) where the standard has 8 strains.
- Logs: `results/amplicon_16s_s1_l3.log`, `_s2_l3.log`, `_s1_l1.log`, `_s1_l2.log`;
  per-check JSON in `results/validation/amplicon-16s-*.json`; kept work dir:
  `results/work/amplicon-16s-s0001-l3/`.
- Cleanup: removed all my task images and my tools image (`alegen-tools-1d9afa22bd69:local`),
  `results/gen/amplicon_16s-*`, the other work dirs and the scratch data.

## Not done / open
- Only seeds 1 and 2 went through the full battery (generation alone was probed for seeds 1-4).
- Why `uchime3_denovo` misses short-tail recombinants was not investigated.
- `score.py` has an unused local (`seq_of`); left untouched so that the verifier in the reported
  runs is byte-identical to the file in the repository.

## 16:52 Byte-level determinism, last code change, all gate runs repeated
- `git status` showed the regenerated `*.fastq.gz` of seed 2 as modified against the committed
  first run. Decompressed content was identical (md5 of S01_R1 and S11_R2 checked); the
  difference was the timestamp `gzip.open` writes into the header.
  Fix: `gzip.GzipFile(..., mtime=0)`. Verified by generating seed 2 / level 3 twice into different
  directories: 42 files, identical sha256; decompressed FASTQ identical to the validated task;
  the only other difference to the validated task dir was `tests/score.py`, from which the unused
  local `seq_of` was removed in the same step.
- Because two files changed, all gate runs were repeated once more with the final code
  (seed 1 L3 with Harbor, seed 2 L3, seed 1 L1, seed 1 L2). The metric tables of the repeated runs
  were compared with `diff` against the previous runs. The gate lines in the report are from these
  last runs.
- Considered and rejected: narrowing the orientation split from 42-58 % to 45-55 % to widen the
  margin below the count band. It would redraw every instance and invalidate the development
  measurements (exact-count yield, mock library evidence); recorded as a limitation instead.

## 17:08 Final state
- Final runs, all with the final code: seed 1 L3 (Harbor) ok 7.3 min; seed 2 L3 ok 3.7 min;
  seed 1 L1 ok 0.8 min; seed 1 L2 ok 2.5 min. Metric tables identical to the runs before the last
  code change for all four.
- Removed all `alegen-amplicon-16s-*` images, the tools image, `results/gen/amplicon_16s-*`, all
  work dirs except `results/work/amplicon-16s-s0001-l3/` (requested `--keep-work`), and my scratch
  directory.


# 2026-11-12 Level 7: M_control_swap (a blank and a specimen that changed places)

Spec `docs/specs/l7/amplicon_16s.md`. Files: `families/amplicon_16s/l7.py` (new), `__init__.py` (dispatch for
`level == 7`, docstring), `assets/amplicon_pipeline_l7.py`, `assets/solve_alt_l7.py`, `assets/solve_l7.sh` (new).
Level-3 assets, scorer and texts untouched.

## 19:35 Design
- Level 7 = `generate(seed, 3, workdir)` unchanged, then post-processing in `l7.py`: two FASTQ pairs filed under each
  other's name, truth permuted, battery rebuilt. No change to the level-3 code path, so levels 1-6 cannot move.
- The FASTQ pairs are decompressed and written again under the other name. Mapping the level-3 files to the other
  name would have left the old file name in the gzip header (`GzipFile` stores it), i.e. a label that says what
  the file "really" is. Checked on the built tasks: header name = file name, mtime 0, reads identical to the
  level-6 file of the other sample; read names carry no sample information.
- Truth: counts follow the files (spec). `samples[].type` in truth.json = kind of template the file holds, so the
  unchanged scorer takes `sample_counts` over the libraries that hold a specimen. Alternatives considered: type as
  in the sheet (then the scored "specimen" column would hold ~70 cross-talk read pairs of a blank: noise); neither
  of the two files scored (accepts a repaired table, but no longer sees a zeroed specimen column, and is not what
  the spec's sentence about the count columns implies). Kept the first; consequence recorded as a weakness.
- Host venv has no scipy: Poisson tail written in `l7.py`, compared with `scipy.stats.poisson.sf` inside the tools
  image on 90 (k, mu) points: same side of 1e-6 everywhere, relative error < 1e-6 where the value is not tiny.

## 19:45 Solutions
- Oracle: Bray-Curtis between libraries + ASVs above 0.5 %, no qPCR. First rule ("labelled blank at least as rich
  as the poorest other labelled specimen") was wrong before it ran: on seed 1 the planted numbers showed a real
  specimen within 0.245 of the exchanged specimen (below the 0.25 limit) and the blank-like library among the
  labelled specimens (richness 6), so the true blank qualified as well and the result would have been "ambiguous".
  Changed to: labelled blank >= half the median richness of the other labelled specimens > labelled specimen.
- Second solution: mass balance against the qPCR column. ASVs whose log10(relative abundance x copies) has a robust
  spread <= 0.3 decades over the libraries (present in >= 80 % of them); a library is off when the median over
  these ASVs is >= 1 decade from their common level; exchange accepted when exactly two libraries are off in
  opposite directions and both come back within 0.5 decades with the rows exchanged. No sample types, no profile
  comparison.
- Caught on re-reading before the first run: the refile loop read the second source from the mapping it had just
  overwritten; sources are now collected first.
- Both level-7 solutions on the level-3 inputs of seeds 1 and 2 (no exchange): nothing reported, 1.000.

## 19:50-20:05 Development harness (existing agent image, inputs mounted, score.py on the host)
- Whole battery, seeds 1 and 2: oracle 1.000, second solution 1.000, `l6_reference` 0.55, `l6_second_solution` 0.55,
  `sheet_roles_everything_else_right` 0.55, `drop` 1.000, `drop_table` 0.70; level-3 entries reproduce the level-3
  numbers to the last digit.
- Core entries, seeds 3-6: same picture (oracle and second solution 1.000 four times; `l6_reference` 0.55 four
  times; uchime3 probe 1.000 / 0.982 / 0.922 / 0.979).
- `--roles drop_blank` (only the odd blank left out), seeds 1 and 2: 1.000. Not put into the battery; recorded as
  a weakness (enrichment test needs one clean blank only).
- Planted-count survey of all 2 x 5 exchanges per seed, seeds 1-16: level-3 generation fails on seeds 7, 12, 14, 15
  ("pair split did not converge", code of 2026-11-11, not touched); on the other 12 seeds all 120 exchanges give
  5 of 5 reagent ASVs missed and score 0.55 with the sheet's roles, exact contaminants with the real roles and
  with both libraries left out.
- Generator condition on richness was first ">= the poorest other specimen": rejected 2 of 120 exchanges (seed 5,
  the exchanged specimen is itself the poorest). Relaxed to >= 0.75 x / <= 0.5 x of the poorest other specimen.
- Added `counts_filed_under_the_other_sample_id` after asking what a solver who found the exchange would write
  into the table: a repaired table scores 0.70. That is what the spec's definition gives; reported, not changed.

## 20:09-20:45 Gate runs
- First `make amplicon_16s --seed 1 --level 7`: all gates ok incl. Harbor, 7.12 min. Afterwards the level-7 judge
  paragraph was corrected ("six / five" profiles below 0.80 instead of "five / six") and the trap sentence about
  the second method reworded, so both gate runs were repeated with the final code:
  - `make amplicon_16s --seed 1 --level 7`: `{"task": "amplicon-16s-s0001-l7", "ok": true, "minutes": 7.08}`
    (34 checks, Harbor oracle [1.0], nop [0.0], agent image derived from the level-3 image).
  - `make amplicon_16s --seed 2 --level 7 --no-harbor`: `{"task": "amplicon-16s-s0002-l7", "ok": true, "minutes": 6.01}`.
  - Scores of every entry equal those of the development harness.
- Disk: 3.8 / 4.5 GB free before the two builds; `docker builder prune -af` after each make.
- Byte identity: `git archive HEAD pipeline/alegen` into scratch, seed 1 at levels 3, 6, 1, 2 and seed 2 at levels
  3, 6 generated with HEAD code and with the final code (`--no-validate --out ... --results ...`), `diff -r`:
  no difference (42 files per task). `tasks/amplicon-16s-s0001-l3`, `-s0001-l6`, `-s0002-l6` are identical to what
  the final code writes. Level-7 tasks regenerated into scratch: identical to `tasks/`.
- Agent-visible content of level 7 against the level-6 task of the seed: four FASTQ files differ, nothing else.

## Left in place
- `tasks/amplicon-16s-s0001-l7`, `tasks/amplicon-16s-s0002-l7`, `results/validation/amplicon-16s-s000{1,2}-l7.json`,
  `results/gen/amplicon_16s-s000{1,2}-l7` (13-14 MB each, written by the cli), images
  `alegen-amplicon-16s-s000{1,2}-l7-{agent,verifier}:local`.

## Not done / open
- Whether a repaired count table should be accepted (see docs, level 7, weaknesses).
- Gate runs for seeds other than 1 and 2; an exchange with a low-biomass specimen or the mock; more than one
  exchange.
- The level-3 generator's failure on seeds 7, 12, 14, 15.


## 20:44-21:04 Maintainer decision: the count table is accepted under either filing
- Decision received after the first hand-over: accept both readings of the count table; scorer = the unchanged
  level-3 logic run on the table as submitted and with the two exchanged columns swapped (both together), better
  total counts; `counts_filed_under_the_other_sample_id` becomes a probe; add an attack with only one column moved.
- `assets/score_l7.py` (new): imports the level-3 verifier as a module (`tests/score_l3.py` = `assets/score.py`,
  byte-identical, checked with `cmp` in both task directories) and runs its `main()` once per reading. The
  level-3 file is a script whose last line is `run(main, ...)`; it is imported with `_scoring.run` replaced by a
  function that records the call, so the file itself needs no change. The second reading swaps the values of the
  two columns in the rows the level-3 code gets from `ctx.table("asv_table.tsv", ...)`; the pair of ids comes from
  `truth.json` (`plants.control_swap`). `details.json` carries gates and components of the reading used and
  `info.count_table_reading` (which one, both totals). Equal totals: the table as submitted.
- `l7.py`: `tests_files` gets `score_l3.py` and the new `score.py`; probe `counts_filed_under_the_other_sample_id`
  (>= 0.99); attack `one_column_moved_specimen_reads_lost` (<= 0.75); judge paragraph on the two filings.
- Measured with the new verifier (planted truth, right contaminants, the two columns filed in seven ways; seeds 1
  and 2 identical): as filed 1.00; both moved 1.00; blank library in both columns 0.70; blank moved and the other
  column zero 0.70; both zero 0.70; **specimen library in both columns 1.00; specimen moved and the other column
  zero 1.00**. The last two are inconsistent tables that the two-reading rule cannot reject, because the column
  that would have to hold the blank library is a blank column and blank columns are not scored at any level.
  The request was that a table with only one column moved stays below the pass line; that holds for one of the
  two possible moves only. Not tuned (a check on the blank column would be new scoring); the attack in the battery
  is the move that loses the specimen's read pairs, the other one is in the docs table and in the judge paragraph.
- `l6_reference`, `l6_second_solution`, `sheet_roles_everything_else_right`: 0.55 under both readings, both seeds.
  No other battery entry changed its score.
- Gate runs with the final code:
  - `make amplicon_16s --seed 1 --level 7`: `{"task": "amplicon-16s-s0001-l7", "ok": true, "minutes": 6.48}`, 35 checks,
    Harbor oracle [1.0], nop [0.0].
  - `make amplicon_16s --seed 2 --level 7 --no-harbor`: `{"task": "amplicon-16s-s0002-l7", "ok": true, "minutes": 5.24}`,
    33 checks.
  - 5.0 GB free before the builds, `docker builder prune -af` after each; the previous level-7 images were replaced
    under the same tags (no untagged leftovers).
- Byte identity re-checked for seed 1, levels 3 and 6: `git archive HEAD` (now `ba7d1b3`, which changed no file the
  generator or the cli uses) against the final code, `diff -r`: no difference, 42 files each; equal to
  `tasks/amplicon-16s-s0001-l3` and `-l6`. Level-7 tasks regenerated into scratch: identical to `tasks/` (43 files).
- Superseded by this entry: the "Not done / open" point about the repaired table above, and the level-7 scores of
  `counts_filed_under_the_other_sample_id` (0.70 as an attack) in the entries above.
