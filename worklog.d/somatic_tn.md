# Worklog: family `somatic_tn`

All entries 2026-11-11, build-host clock (UTC); times are approximate, reconstructed from file and
log timestamps.

## 17:03 - reading, plan
- Read AUTHORING.md, the spec, the `crispr_screen`, `rnaseq_confound` and `variant_review` families
  and the core modules (`spec`, `dockerlib`, `harbor`, `validate`, `scoring`).
- Plan: host-side forward simulator (numpy only) on an hg38 chr22 slice; paired 2 x 100 FASTQ for
  three tumour/normal pairs; oracle = bwa + pysam fragment counting; independent solution =
  minimap2 + bcftools mpileup; no GATK (image size), no turnkey somatic caller.
- Design decisions taken before writing code:
  - Orientation artefacts are planted as read-pair-orientation artefacts (F1R2 / F2R1, the real
    OxoG / FFPE-deamination signature), not as forward/reverse strand bias. Alt reads are then
    balanced over the two alignment strands, so the classic strand-bias filter does not find them.
  - The duplication is exact, and the duplicated-segment variants sit more than the largest insert
    size away from the copy edges, so that *no* read pair can place them. Their true position is
    therefore gray (unscored), the mirror position in the other copy is a false positive.
  - VCFs are named by patient and the patient is anchored to the tumour library in the output
    definition, because "normals exchanged" and "tumours exchanged" cannot be told apart from
    genotypes.
  - Scoring uses precision and recall separately (min over patients) instead of one F1: with ~30
    true variants, 9 mirror calls only lower F1 to 0.87, which a ramp cannot separate cleanly from
    a good analysis with two slips.

## 17:08 - tools image
- Bioconda versions looked up through the anaconda.org API. `dockerlib.tool_image(conda=[bwa=0.7.19,
  minimap2=2.30, samtools=1.25, bcftools=1.25, pysam=0.24.1])` built in 90 s on the first try:
  `alegen-tools-1c4673c508ab:local`, 1.91 GB disk usage / 447 MB content as listed by `docker images`.
- The chr22 FASTA was copied from `.cache/variant_review/` into `.cache/somatic_tn/` (sha256 checked
  by `sim.load_chr22`); nothing was downloaded.

## 17:15 - simulator
- `sim.py`: window + target choice using the UCSC soft-mask (targets <= 35% repeat), exact
  duplication inserted into an inter-target gap, uniqueness mask from 35-mers (both strands),
  planting, fragment simulation with per-fragment haplotype / orientation / duplicate count,
  vectorised sequencing errors, deterministic gzip (mtime 0).
- First smoke test: one tumour library (46.9 k pairs at 180x) is 2 x 2.32 MB of FASTQ, i.e. 49
  bytes per read; six libraries at 180x / 90x come to 21 MB. Generation takes 11 s.

## 17:21 - first end-to-end run (seed 1, level 3, in the tools image)
- Harness mistake: my `run_local.py` took a stray `-v` as its scratch directory and wrote into
  `pipeline/-v/` and `results/dev/`; both removed, scratch moved to the session scratch directory.
- Oracle: all 95 scored variants found, pairing concordance 1.00 vs 0.41-0.48, **one false positive**:
  `27111 G>C`, 6 reads at AF 0.05, right next to the planted deletion `27110 CG>C`. bwa draws a
  1-bp deletion that lies within a few bases of the read end as a mismatch. This is the kind of
  alignment artefact the instruction excludes, not a planted trap.
  - Change: the oracle no longer uses the first/last 5 bases of an alignment as SNV evidence. The
    independent solution keeps BAQ on and drops an SNV within 10 bp of a stronger indel call.
  - Cost: the oracle's allele-fraction error rises from 0.003 to 0.009 (10% fewer reads per site);
    the full-credit bound is 0.03.
- Independent solution (minimap2 + bcftools mpileup on orientation-split BAMs, no dedup): 95/95, no
  false positive, on the first run. bcftools' indel allele depths were good enough at this depth
  (all planted indels found, AF error 0.009), so no CIGAR-walker fallback was needed.
- Naive variants (same instance): sheet pairing 0.34, no MAPQ filter 0.61, no orientation filter
  0.51, tumour-only 0.40, everything ignored 0.20.

## 17:27 - rehearsal `make` (seed 1, level 3, Harbor) and more seeds
- Rehearsal: all gates ok in 6.4 min, `harbor_oracle rewards=[1.0]`, `harbor_nop rewards=[0.0]`.
- Agent image: 2.0 GB disk usage / 491 MB content (`docker images`).
- Added `--strand-instead-of-orientation` to the pipeline and a `naive_strand_bias_only` attack: a
  forward/reverse filter leaves every artefact in (score identical to "no orientation filter").
- Verifier hardening after re-reading it against hostile inputs: `patient` cell may be `None` for
  a short row (would have raised `TypeError` = verifier error), `int(inf)` raises `OverflowError`,
  VCF size cap lowered from 20 MB to 3 MB compressed.
- **Seed 3 failed to generate** (`RuntimeError: could not place a site`): 13,962 eligible positions
  and ~270 sites with 40-bp exclusion each. Changes: spacing 25 bp (germline) / 30 bp (somatic SNV,
  artefact) / 60 bp (indel), and `build_reference` now rejects a window with fewer than 15,000
  eligible positions. Seeds 1-8 generate; this changed the random stream, so every number measured
  before 17:38 was re-measured.
- Added generation-time assertions (>= 20 scored variants and >= 3 scored indels per patient, >= 9
  mirror sites per patient, enough strong artefact sites) so that a weak instance fails loudly
  instead of shipping.

## 17:38 - first complete set of runs (single damage chemistry per seed; superseded below)
- `make --seed 1 --level 3 --keep-work`: all gates ok incl. Harbor (6.7 min); `--seed 2 --level 3`,
  `--seed 1 --level 1`, `--seed 1 --level 2` (all `--no-harbor`): ok.
- Determinism: `truth.json` and all six R1 files of seed 1 level 3 had identical md5 in two
  independent generator processes.
- Discoverability evidence (per-target MAPQ, orientation and strand tables, concordance) computed
  from the oracle's BAMs with `dev/evidence.py`.
- The coordinator asked for finished instances' images to be removed: `docker rmi` of the
  `somatic-tn-s0001-l3` agent/verifier tags; from here on every instance's tags are removed right
  after its run.
- Seed 2 showed only 7 sheet-pairing false positives in one of the two affected patients (the other
  had 16). Added a generation-time assertion that at least one affected patient has >= 12 germline
  sites it does not share with the owner of the exchanged normal. `oracle_runtime_min` set to the
  measured 0.4 (was a guess of 2.0), `oracle_steps` to 33.

## 17:52 - attacking my own verifier: a blanket spectrum filter passed
- Question I had left open in the first draft of the design notes: what does "drop every low-fraction
  call of the damage substitution class in the FFPE tumours" score? Measured with
  `dev/spectrum_filter.py` on the `naive_no_orientation_filter` outputs: **0.91-0.97 on both seeds**
  (cut-off 0.20-0.25: precision 1.0, recall 0.91). The orientation trap could be bypassed without
  ever looking at a read. Reported here because it is the most important finding of the session.
- Change 1 (generator): a damaged tumour now carries artefacts of *both* chemistries (deamination
  C>T/G>A and oxidation G>T/C>A, 5-11 sites each), so the artefacts are not one substitution class,
  and the artefact rate goes up to 56% of the affected orientation (apparent allele fraction 8-30%).
  The chemistry no longer varies with the seed; which tumours are damaged, how many sites, the mix
  and the positions do. Notes: damaged tumours are FFPE sections; the shearing-batch story is gone.
  All five `make` runs passed again (seed 1 level 3 with Harbor).
- Re-measured: best blanket filter 0.94 on seed 1 (FFPE tumour only, cut-off 0.30: the damaged
  tumour happened to have only 6 of 23 true SNVs in the two classes), 0.87 on seed 2. Still a pass
  on one seed.
- Change 2 (generator): five SNVs of the first subclone of a damaged tumour are planted in the damage
  classes (allele fraction 0.13-0.29, inside the artefact range). Change 3 (verifier): recall band
  moved from 0.80-0.92 to 0.85-0.95; every correct route measured so far has recall 1.000, one miss
  in 28 is 0.964.
- Re-measured with the pipeline's new `--no-orientation --spectrum-filter X` switch (tools image):
  0.60 / 0.70 / 0.70 for X = 0.25 / 0.30 / 0.35 on seed 1 and 0.64 / 0.70 / 0.70 on seed 2 (8-13
  true variants lost in the damaged tumour; before change 2 and 3 these were passes). Added to the battery as `shortcut_spectrum_filter`
  with the naive ceiling 0.75.

## 18:10 - final runs (finished 18:25)
- `make --seed 1 --level 3 --keep-work` (all gates incl. Harbor), `--seed 2 --level 3 --no-harbor`,
  `--seed 1 --level 1 --no-harbor`, plus `--seed 3 --level 3` and `--seed 1 --level 2` as extra
  checks. Gate lines are in `results/somatic_tn_s*_l*.log`.
- All five passed: seed 1 level 3 in 4.4 min incl. `harbor_oracle rewards=[1.0]` and
  `harbor_nop rewards=[0.0]`; `shortcut_spectrum_filter` 0.700 / 0.700 / 0.589 on seeds 1-3.
- Analyst-variation runs on the final seed 1 and seed 2 instances (dev harness, tools image): oracle
  without read-end margin, with stricter thresholds (8 reads / 6%), with looser thresholds (4 reads /
  3%); independent solution without the indel-proximity rule, and additionally without BAQ: all 1.0
  with 0 FP / 0 FN.
- Determinism re-checked on the final code: `truth.json` and all twelve FASTQ files of seed 1 level 3
  have the same md5 in the `make` run and in a separate `dev/score_files.py` process.
- Images: every instance's agent/verifier tags were removed after its run; the tools image is
  removed at the end.
- No pip packages were added to the venv (pysam 0.24.1 was already there; it is used only by
  `dev/evidence.py`).
