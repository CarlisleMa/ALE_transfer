# Worklog: family `rnaseq_confound`

All entries 2026-11-11, build-host clock (UTC); times are approximate, reconstructed from file and log timestamps.

## 15:33 - reading, plan
- Read AUTHORING.md, the spec, the `crispr_screen` reference family and the core modules.
- Plan: forward simulator on real Ensembl transcripts (sequence realism only), negative-binomial
  counts from planted effects, single-end 50 nt reads; oracle = salmon + pyDESeq2; independent
  solution = own k-mer counter + voom-style linear model; R/Bioconductor in the image if it fits.

## 15:35 - reference download
- `curl -fsSL` of Ensembl release 113 `Homo_sapiens.GRCh38.cdna.all.fa.gz` (79.7 MB) and
  `Homo_sapiens.GRCh38.ncrna.fa.gz` (41.0 MB; XIST is a lncRNA and not in the cDNA file) into
  `.cache/rnaseq_confound/`. sha256 recorded in `sim.py` and in the design notes.

## 15:37-15:52 - tool image: four failed attempts, one finding about the core
- `salmon=1.10.3` + `r-base=4.5.3` + Bioconductor 3.22: unsolvable (salmon 1.x builds pin old
  icu/libdeflate/zlib). `salmon=1.12.1`: same. `salmon=2.11.0` (Rust rewrite, same CLI: `index -t -i
  -k`, `quant -i -l A -r -o --validateMappings`) solves.
- `dockerlib.tool_image` then died with `UnicodeDecodeError: 'utf-8' codec can't decode byte 0xc0`
  inside `dockerlib.build_image` (`subprocess.run(..., text=True)` on the captured build log). The
  core is not mine to edit, so the family got its own Dockerfile string (`image.py`): same base and
  tail layers as `dockerlib`, but `micromamba create -q ... > /tmp/mm.log`, and the log is printed
  ASCII-filtered only on failure.
- Next attempt: "No space left on device" (shared quota; the coordinator pruned the BuildKit cache
  and asked for a small footprint). Retried, built in 202 s.
- Result: with the R stack the conda layer is 2.63 GB and `docker images` reports 3.82 GB
  (`du -sh /opt/tools` = 2.5 GB). Above the 3 GB budget -> **the shipped image is Python only**
  (salmon, pyDESeq2, statsmodels). The R image was used as a development image for calibration
  with real DESeq2 / limma-voom / edgeR and removed once that was done (about 16:10).

## 15:45-15:54 - simulator
- `sim.py`: eligible-transcript cache (17,135 protein-coding genes, <= 3 transcripts each, plus XIST
  and five chrY genes), seeded gene subset with a 31-mer uniqueness filter, design, planted effects,
  NB counts, vectorised read sampling. First version hashed k-mers with Python `hash()` (salted per
  process); replaced by exact 2-bit uint64 codes.
- First instance: 149-198k reads per sample = 35.9 MB of input. Too big. Measured 19.7 gzip bytes
  per 50 nt read (15.1 for the bases, 3 for the numeric name), not the 14 I had assumed. Fixes:
  expected library size is normalised after planting (up-regulation adds more reads than
  down-regulation removes) and the target went to 98k reads per library.

## 15:50 - design decision: 5+1 / 1+5 instead of the spec's 4+2 / 2+4
- Worked out on paper first: with 4/2 vs 2/4 a condition-only model sees a batch-only gene with an
  apparent effect of b/3 against a within-group SD of 0.47 b, t ~ 1.2 whatever b is. With 5/1 vs 1/5
  the apparent effect is 2b/3, t ~ 3. Confirmed by measurement at 16:45 (below).
- Consequence for the mix-up: both exchanged libraries have to come from a 5-sample cell, otherwise
  excluding them leaves a session with a single condition. They are therefore in different prep
  sessions and the scenario is "two FASTQ files carry each other's library".

## 15:56 - first end-to-end run on real reads (dev image, 0.5% error rate, 148k reads)
- oracle (salmon + pyDESeq2): count agreement 0.936, 78/79 planted genes, 0 false. 25 s.
- alt (k-mer counter + voom-style, two samples excluded): agreement 0.985, 77/79, 4 false. 9 s.
- salmon leaves 6.4% of the reads unassigned. A default run must earn full credit, so the simulated
  error rate went from 0.5% to 0.3% (agreement 0.962 afterwards; full credit from 0.90).

## 15:58-16:07 - calibration on true counts, 2 seeds, 22 analyses each (98k reads)
- First run of the harness died in `solve_alt.py --counts` (it still opened the transcript FASTA,
  and `set -e` then skipped the R part); fixed the development switch and repeated.
- All correct analyses (pyDESeq2, R DESeq2, R limma-voom, R edgeR, own voom; records exchanged or
  samples excluded) score 1.000: FDP 0-0.062, sensitivity 0.924-0.975, median |LFC error| 0.166-0.241.
- No batch term: FDP 0.38 (edgeR QL) to 0.60 (DESeq2), scores 0.45-0.63.
- Labels as given: sensitivity 0.10-0.23, scores 0.27-0.38. Both ignored: 0.20-0.22.
- Exchanging only the condition of the two samples (sheet batch kept): FDP 0.24-0.42, score
  0.60-0.80. Below the pass mark; reported, not gated.

## 16:05 - scoring changes after thinking about trivial answers
- `sample_sex` / `sample_labels` are credited only if `gene_counts` earns credit; otherwise a QC
  table copied from the sample sheet earned 0.15 at levels 1-2 with no analysis at all.
- FDP credit is scaled by ramp(sensitivity, 0, 0.5) and sensitivity credit by ramp(FDP, 0.60, 0.30),
  so "nothing significant" and "everything significant" do not collect a rate for free.
- A planted gene called significant with the wrong sign counts as a false discovery.

## 16:07-16:21 - first complete `make` runs (98k reads per library) - all gates ok, later superseded
- seed 1 level 3 with Harbor (5.4 min), seed 2 level 3 (2.4 min), seed 1 level 1 (2.7 min).
- Shipped (Python-only) agent image: 1.98 GB by `docker images` (`du -sh /opt/tools` = 1.1 GB).
- Removed the R development image and the finished task images (coordinator's request).
- Discoverability measured on the oracle output of seed 1: XIST 1337-2044 CPM in female libraries
  vs <= 11 in male ones, chrY genes 5936-7892 CPM vs 0; PC1 of log-CPM (57% of variance) splits the
  two 5-sample cells, PC2 (17%) isolates the two single-sample cells, and both mixed-up samples sit
  in the other majority cluster.

## 16:21-16:24 - problems found in my own first version; all runs repeated afterwards
1. Level-1 input was 25.41 MB, over the 25 MB budget (the level-3 instances were 24.5 and 23.4 MB).
   Library size varies by a few percent with the seed, so the target went from 98k to 93k reads per
   library (same genes, effects and design; about 5% fewer reads). Sizes now 22.2-22.9 MB.
2. Hostile-input pass over the verifier on the host (short rows, binary garbage, CSV instead of
   TSV, inf, negative counts, duplicated genes, gzip'ed tables, counts x1000): a short row made
   `r["inferred_sex"].strip()` hit `None` -> AttributeError -> verifier error instead of score 0.
   All cells now go through a `cell()` helper. The `significant`/`padj` consistency gate got 0.005
   slack because a padj written with three decimals can land on the other side of 0.05.
3. The same pass showed that perfect counts and DE results with a wrong `condition_used` column
   scored exactly 0.90 = pass. Weights changed from 0.05/0.10 to 0.04 (sample_sex) / 0.11
   (sample_labels) so that a wrong sample table ends at 0.89.
4. Median |LFC error| of correct analyses was 0.17-0.24; the full-credit edge moved from 0.30 to
   0.33 (zero stays at 0.55). The one-third attenuation attack sits at 0.71-0.75.
- The real-R calibration was NOT repeated at 93k (the R image is gone; rebuilding costs 3.8 GB of
  the shared disk). The design notes label those rows as 98k / first weights.

## 16:24-16:35 - final `make` runs (93k reads, final verifier) - the ones reported
- `--seed 1 --level 3 --keep-work`: ok, Harbor oracle [1.0], nop [0.0], 3.36 min.
- `--seed 2 --level 3 --no-harbor` (with `--keep-work` added to read the raw metrics): ok, 4.81 min.
- `--seed 1 --level 1 --no-harbor` (same): ok, 2.41 min.
- Generator determinism checked: two independent processes give the same md5 of the truth for
  seed 1 level 3, identical to the shipped `truth.json`.

## 16:35-16:48 - supplementary measurements with the final code
- Five seeds on true counts, Python analyses: correct ones 1.000 on every seed (FDP <= 0.051,
  sensitivity >= 0.924, LFC error <= 0.227); no batch term 0.45-0.53; labels as given 0.25-0.41;
  both ignored 0.19-0.22; condition-only exchange 0.56-0.82.
- The spec's 4/2 vs 2/4 design, three seeds: the condition-only analysis scores 1.000 every time
  (FDP 0.013-0.052). This is the measurement behind the 5/1 deviation.
- `make --seed 1 --level 2 --no-harbor`: all gates ok (4.03 min). Own voom on all 12 samples has 8
  false calls there (FDP 0.096, still full credit; thinnest margin seen on the correct side).

## 16:50 - clean-up
- Removed my task images and the FASTQ copies under `results/gen` and in the work directories
  (`harbor_task`, `attack/inputs_copied`); kept every `details.json`.

# Level 4 (2026-11-12, UTC, build-host clock)

## 00:37 - reading, plan
- Read AUTHORING.md ("Level 4"), `docs/specs/l4/rnaseq_confound.md`, the level-3 notes and code, the core modules.
- Constraint that shaped the code: `tasks/rnaseq-confound-s0001-l3` must regenerate byte-identically, including
  `solution/` and `tests/`. Nothing of level 1-3 is edited: level 4 lives in a new module `l4.py` and new assets
  under `assets/l4/` (own verifier, own oracle, own second solution). `__init__.py` got three small branches
  (`level >= 4`): another base seed, the mixed library, dispatch to `l4.generate`.
- Level 4 uses its own base seed (64000 + seed instead of 52000 + seed), so the level-4 task of a seed is a
  different instance from the level-3 task of that seed (other genes, effects, design, mix-up).

## 00:50 - tool image for generation
- The spec wants the predecessor's files "produced by running mutated copies of the reference pipeline". That
  needs salmon and pyDESeq2 at generation time. `dockerlib.tool_image` would create a second 1.9 GB conda layer
  (its RUN line differs from the family's own Dockerfile string). Instead `l4.ensure_gentools()` builds
  `alegen-rnaseq-confound-gentools:local` from the *same Dockerfile text as the task's agent image* with an empty
  input directory: the conda layer is shared with every task image of the family, and the predecessor's scripts
  run on exactly the software the agent gets. Built in 99 s, 1.88 GB by `docker images`.
- First attempt failed (`ModuleNotFoundError: alegen`, script started without PYTHONPATH). No other build problem.
- Probed salmon 2.11.0 in that image on one level-3 FASTQ: `-l A`/`-l U` assign 89,247 of 92,685 reads, `-l SF`
  44,671 and `-l SR` 44,582 (48.2% / 48.1%), `quant.sf` is byte-identical between two runs and between 2 and 4
  threads, the first 30,000 reads give 28,857 assigned reads.

## 00:55-01:06 - count-level measurements that decided the design (planted counts, seeds 1-3)
- Mixed library kept in the model (male/female typing only), pyDESeq2 `~batch + condition`, records exchanged:
  sensitivity 0.899 / 0.911 / 0.962, FDP 0.014 / 0 / 0, median LFC error 0.222 / 0.191 / 0.200. All inside the
  level-3 bands. **One mixed library of twelve does not damage the DE table**; the Q element can only bite through
  the sample table. Consequence: weight on `sample_labels` plus a stated release rule (see design notes).
- Correct analyses: oracle chain 0.987 / 0.937 / 0.987 sensitivity (FDP <= 0.037), with all three doubtful
  libraries excluded (9 samples left) 0.949 / 0.899 / 0.987, voom-style second solution 0.962 / 0.949 / 0.962.
- No session term: FDP 0.359 / 0.263 / 0.427 (lower than at level 3; seed 2 is the thin one).
- Raw p-values written as `padj` (spec candidate): FDP only 0.194 / 0.084 / 0.225 -> not a decision point. The
  checklist item stays, the defect is never planted.
- Predecessor's "first pass" (condition-only, sample sheet as delivered), top 40 by adjusted p-value: 33 / 35 / 36
  planted responders and only 7 / 5 / 4 others. An erratum scored by MCC on 4-7 negatives is a lottery. The mailed
  list is therefore composed by the generator from that ranking (16-24 responders, the rest non-responders, both
  taken from the top of the first-pass ranking).

## 01:07-01:15 - generator, verifier, oracle, second solution
- `l4.py`: mixture (male library + 25-45% of a female library of the other condition and session; fraction raised
  when needed so that the minor marker share is >= 2.5 x the SOP limit), defect plan, script templates, the
  container run, SOP, instruction, battery.
- First generation run failed after 74 s: `build_previous` still read `mail_top_responders.tsv` from the container
  although the list had moved to the host. Fixed.
- Oracle `takeover_pipeline.py` (fresh salmon, SOP identity control, pyDESeq2, audit by re-running the shipped
  scripts, erratum) and `solve_alt_l4.py` (k-mer counter, voom-style model, exclusion, audit by profile distance to
  the k-mer counts plus static inspection) both scored 1.000 on seed 1 at the first run.

## 01:15-01:32 - careless variants on seed 1 (dev harness, tool image), seeds 2-3 for the two solutions
- seed 1 (planted: TPM matrix, one column from another FASTQ, stale draft): level-3 reference 0.470, inherited
  matrix trusted 0.497, identity criterion not applied 0.670, sample sheet as delivered 0.324, TPM kept 0.660,
  wrong-file column kept 0.646.
- Variants of defects that are *not* planted on seed 1 (so the forced `ok` is the right verdict and only the
  downstream damage shows): stranded library type 0.830, partial file 0.830, no session term 0.765, reversed
  contrast 0.520. With the audit miss they would lose a further 0.13-0.17.
- Raw p-values as padj: 0.913 (FDP 0.202); with an audit miss about 0.75-0.78. Confirms dropping it.
- Seeds 2 and 3: oracle, second solution and oracle-with-exclusion all 1.000 with a complete audit (seed 2: stranded
  `SR`, partial file, no session term, stale draft; seed 3: stranded `SF`, wrong file, partial file, stale draft).
- Hostile-input pass over the level-4 verifier on the host (23 cases: missing/short/duplicated/unknown rows in
  `audit.tsv` and `erratum.tsv`, CSV instead of TSV, random bytes, upper-case verdicts, gzip'ed table, `inf`):
  no verifier error; every malformed file scores 0, upper-case verdicts are accepted.
- Level 3 seed 1 regenerated into a scratch directory: `diff -r` against `tasks/rnaseq-confound-s0001-l3` empty.
- Level 4 seed 1 generated twice: `diff -r` empty apart from a metadata line I had edited in between
  (`oracle_steps`), i.e. the salmon / pyDESeq2 products are byte-reproducible.

## 01:32-01:53 - blocked by the disk rule, used the time for more instances
- Free disk fell to 7.4-7.7 GB (other authors' builds). I had no stale image to remove (only the tool image, which
  shares its layers with the task images), so no `make` was started for about 20 minutes. Removed my own scratch
  copies; that freed 0.1 GB.
- Meanwhile, without any build: seeds 4-6 generated into scratch and run through the dev harness. Oracle and
  second solution 1.000 on all three with a complete audit; identity criterion not applied 0.670 on all three;
  inherited matrix trusted 0.490 / 0.630 / 0.344.
- Spec element (a), paired donors with donor partly aliased with session, simulated at count level (numpy on the
  host, pyDESeq2 in the tool image, three instances, donor SD 0.5 and 1.0 log2). Unpaired model: sensitivity
  0.81-0.91 (SD 0.5) and 0.39-0.42 (SD 1.0). Correct paired model: 0.90-0.94 and 0.785-0.899. No setting has margin
  on both sides with the level-3 bands -> not implemented; numbers are in the design notes.
- Wording change after thinking about another CPU type: the predecessor's DE script now prints six significant
  digits and the SOP item `draft_reproducible` says "the same numbers, apart from rounding in the last printed
  digit".

## 01:53-02:08 - first complete `make` of seed 1 (earlier defect draw): all gates ok, then superseded
- One `make` was started and killed by me after two minutes (still generating): I had found a wrong statement in
  the judge text (MCC of one wrong verdict on eight items is 0.745-0.775, not 0.76-0.80; I had computed it for nine
  items). Fixed in `l4.py`, `score.py` comment and the notes.
- `make --seed 1 --level 4 --keep-work`: ok, Harbor oracle [1.0], nop [0.0], 10.46 min. Planted then: TPM matrix,
  wrong-file column, stale draft. Raw metrics archived in my scratch directory and quoted in the notes as "d1".
- While it ran I tabulated the defect draw over 30 seeds: every recompute item was defective in 80-87% of the
  tasks and `merge_read_counts` in 1 of 30 (the rejection rule "no TPM together with a partial file" removed it).
  A policy could rate the three recompute items `defective` unseen. Changed the draw: exactly two of the three
  recompute items, the TPM defect with probability 0.75 when the partial file is not planted. Over 24 seeds:
  13-18 tasks per recompute item, 15-17 per code item, 6 with TPM. Also raised the target for the mixed library's
  expected minor share from 2.5 to 3.5 x the limit (seeds 4, 25, 27 had 0.046-0.049 realised; none of 40 seeds is
  refused, minimum now 0.055).
- The new draw changes the plans of seeds 1 and 2, so both were validated again.

## 02:08-02:30 - second complete runs (final draw)
- seed 1 with Harbor: all gates ok, 11.41 min. seed 2 without Harbor: all gates ok, 6.26 min. Seed-2 images removed.
- Seed 2 is the first instance with a *sound* draft: both solutions reproduce it exactly (maximum difference
  0.0000) and rate it `ok`.
- Thin spot found in these numbers: `defect_kept_matrix_sample_files` scored 0.714 on seed 2 (ceiling 0.75). The
  copied column agrees 0.719 with the truth of the column it sits in, which the level-3 band (zero at 0.50) still
  rewards with half of the count credit.

## 02:30-02:35 - verifier: a count column has to be the column of its own file
- Added to `assets/l4/score.py`: `gene_counts` earns nothing when a column agrees better with the planted counts
  of another FASTQ file than with those of its own. Checked on the host against every output kept from the two
  validation runs before shipping it: the wrong-file variant of seed 2 goes from 0.714 to 0.580, no other score
  changes (oracle, second solution, probes, all other attacks).
- Level 3 seed 1 regenerated with the final code: `diff -r` against `tasks/rnaseq-confound-s0001-l3` empty.
  Level 4 seeds 1 and 2 regenerated: identical to the validated task directories except `tests/score.py` and the
  judge table in `README.md`, i.e. the edit above; all inputs, including the files made with salmon and pyDESeq2,
  are byte-identical between the two generations.
- Third complete runs started (both seeds), because the verifier changed.

## 02:35-02:44 - third complete runs = the ones reported (final code)
- `make rnaseq_confound --seed 1 --level 4 --keep-work`: ok, 32 checks, Harbor oracle [1.0], nop [0.0], 5.09 min,
  oracle 0.39 min. Log: `results/rnaseq_confound_s1_l4.log`.
- `make rnaseq_confound --seed 2 --level 4 --no-harbor --keep-work`: ok, 31 checks, 3.38 min, oracle 0.29 min. Log:
  `results/rnaseq_confound_s2_l4.log`. `defect_kept_matrix_sample_files` now 0.580 (was 0.714 before the column
  rule); every other score of both seeds is the same as in the second runs.
- Seed-2 images removed right after the run.

## 02:45-02:55 - closing checks and clean-up (after an interruption of a few minutes by a usage limit)
- No code was changed after 02:31:51 (`l4.py`, `assets/l4/score.py`); the validation records are from 02:39:58 and
  02:44:13. Confirmed by regenerating with the code as it is now: level 4 seed 1 and seed 2 into a scratch
  directory, `diff -r` against `tasks/rnaseq-confound-s0001-l4` and `tasks/rnaseq-confound-s0002-l4` empty; level 3
  seed 1, `diff -r` against `tasks/rnaseq-confound-s0001-l3` empty.
- Design notes: "Level 4" section filled from the `details.json` files of the final runs (`dev/l4_table.py`).
- `l4.ensure_gentools()` had never run with the image absent (I had built the image by hand first). Removed the
  image, called the function: rebuilt from the build cache in under a second. Then removed the image again;
  generation of a level-4 task rebuilds it on demand (needs the conda download when the build cache is gone,
  99 s at the first build).
- Removed the FASTQ copies in my own result directories (`results/gen/rnaseq_confound-s000{1,2}-l4/fastq`,
  `l4_stage`, `transcripts.fa.gz`; `results/work/rnaseq-confound-s000{1,2}-l4/harbor_task` and
  `attack/inputs_copied/output/fastq`). Every `details.json` and every other output is kept.
- Left in place: `alegen-rnaseq-confound-s0001-l4-agent:local` (1.97 GB) and `-verifier:local` (178 MB);
  `tasks/rnaseq-confound-s0002-l4` (written by the required `make` command for seed 2).

## Not done / open at hand-over
- `oracle_runtime_min` in the task metadata is 1.5; measured 0.29-0.61 min. Not corrected (would need a fourth
  validation of both tasks).
- Element (a) of the spec (paired donors) is measured at count level only and not implemented.
- No agent pilot (not allowed in this round), so whether level 4 is hard for a small model is unknown.
- Real R/Bioconductor methods were not run on the level-4 instances (the R development image was not rebuilt).
