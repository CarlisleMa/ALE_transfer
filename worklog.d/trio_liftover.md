# Worklog: trio_liftover

All entries 2026-11-11 (single session). Times are approximate, in order.

## 1. Reading and plan
- Read AUTHORING.md, docs/specs/trio_liftover.md, the crispr_screen reference family and the core modules.
- Decided up front: plant everything on GRCh38, derive the hg19 input from it; tools only as guards.

## 2. Reference material
- Downloaded UCSC hg19/hg38 chr22 FASTA, both liftOver chains, Ensembl 113 chr22 GFF3 into
  `.cache/trio_liftover/dl/` (URL + sha256 in `refdata.py` and in docs/families/trio_liftover.md).
- Chain analysis: main chr22->chr22 chain (id 23) is plus strand, 1129 blocks, 839 hg19-side gaps; offsets
  hg19->hg38 are +12477 (22q11.21), about -354289, -342187/-396448, -396004 and -438429 further down.
- Picked five 600-kb GRCh38 windows by protein-coding gene density: 88 genes with an Ensembl-canonical
  mRNA, 87 with a clean CDS (ATG..stop, no internal stop). Window 3 turned out to contain the 54-kb hg19
  segment chr22:24,343,276-24,397,527 (GSTT region) that only aligns to an ALT contig in the full chain;
  after restricting the chain to chr22->chr22 it has no GRCh38 image, which is a natural "hg19-only" region.
- Masked FASTA bgzips to ~1 MB per assembly; total input is 2.5 MB.

## 3. Tool image: three attempts
- `picard=3.3.0` pulls r-base: image 3.52 GB (over the 3 GB budget). Switched to `picard-slim=3.3.0`: 1.98 GB.
- **CrossMap 0.7.3 corrupts deletions.** Test VCF with `AT>A`: CrossMap fetches a single reference base for
  REF and writes `REF=A ALT=<empty>`; a padded SNV `CAT>AAT` became `C>CAT`. Read `cmmodule/mapvcf.py`: only
  the first base of REF is lifted and `refFasta.fetch(target_start, target_end)` is one base long. A tool
  that silently damages every deletion cannot be one of the two consensus tools and should not be the only
  lift-over tool besides Picard in the agent image. Replaced it by `bcftools-liftover-plugin=1.21`
  (bcftools +liftover, Genovese et al. 2024). Deviation from the spec, see design notes.
- Removed my two superseded tool images after the coordinator's disk-space notice
  (`alegen-tools-f6d11e242257`, `alegen-tools-cbf5a389fb92`); kept `alegen-tools-5886c792c5f1`.
- Other measured tool facts that shaped the design:
  - `bcftools norm -m-any -f` trims padded SNVs (`CAT>AAT` -> `C>A`), left-aligns right-shifted indels and
    splits `TGCA>T,TGCC` into `TGCA>T` + `A>C` at pos+3 with AD split correctly.
  - Picard LiftoverVcf needs an uncompressed FASTA + .dict; on the contig name `22` it rejects 100% as
    `NoTarget`; **it writes the samples in alphabetical order** (my first consensus check compared GT
    columns by position and flagged the causal de novo as "not reproduced" - fixed by `bcftools view -S`).
  - `bcftools csq` stops with "the fasta reference does not match the VCF REF allele" when hg19
    coordinates are annotated against GRCh38 (also with `--force`).
  - bcftools +liftover lifts some alleles that sit in 1-2 bp hg19-side chain gaps by adding a reference
    allele / swapping alleles, Picard rejects them: these are the gray zone (5-7 alleles per task).

## 4. Generator
- First end-to-end generation: 5-9 s per attempt including the docker guard step.
- Guards that fired during development: seed 2 / attempt 0 "could not place compound-het LoF" (gene without
  a right-shiftable frameshift site that passes all placement rules) -> automatic retry with the next
  sub-seed, attempt 1 passed. No csq disagreement with a designed consequence so far.
- Added `level` to the RNG seed so that level 1 and level 3 of one seed do not share the answer.

## 5. Solutions and naive variants (measured before the first `make`, tools image, seed 1 level 3)
- Oracle (Picard + bcftools csq) 1.000; independent pure-Python solution 1.000 on the first run.
- `--no-lift` crashed at bcftools csq (REF mismatch) -> a crash tests nothing. `bcftools norm -c x` as the
  "drop the offending records" reaction produced an unreadable file, so the careless path now drops
  REF-inconsistent records in Python and carries on; it then reports exactly the planted wrong-build decoy
  as a de novo missense candidate.
- Order-of-operations check (not part of the battery): bcftools +liftover on the raw records, +liftover
  after pre-normalisation, and Picard on un-split records, each followed by `bcftools norm -m-any` on
  GRCh38, all give lifted F1 1.000 / GT concordance 1.000 / unmapped F1 1.000 against the planted truth.

## 6. Scoring calibration
- Weights moved from the spec sketch (0.20/0.10/0.10/0.40/0.20) to 0.15/0.10/0.10/0.45/0.20 so that the
  "earlier artefacts" carry 0.35, and the candidate ramp was set to [0.70, 1.0] so that one false or one
  missing candidate fails (see the table in the design notes).
- `liftover_rejects` is not declared as a trap: both tools hand the rejects over, the instruction asks for
  them explicitly, and ignoring them would cost only the 0.10 of that component (0.90 = pass), so a
  per-trap naive variant could not stay under the required ceiling without distorting the weights.

## 7. First validation round (code state before the fixes below)
- `make --seed 1 --level 3 --keep-work`: all gates ok incl. G6 (oracle 1.0, nop 0.0), 5.95 min (first
  agent-image build), 4.49 min on the repeat.
- `make --seed 2 --level 3 --no-harbor`: all gates ok, 1.12 min.
- `make --seed 1 --level 1 --no-harbor`: **generator crash** (`TypeError: 'NoneType' object does not support
  item assignment` in `fuse`). Cause: the second allele of a shared record was planted with `force=True`
  under a stricter block-edge margin (150 bp) than the first allele (>= 100 bp + span), so `add` could
  return None and nothing checked for it. It had not fired on the level-3 seeds by luck. Fix: a
  `forceable()` pre-check before the first allele of a pair is planted, and a clean `_Retry` if the
  partner still cannot be placed.

## 8. Retry rate
- Measured the pre-docker retry rate over seeds 11-40 at level 3: 20 of 30 seeds needed at least one
  redraw and 2 of 30 would have exhausted the 6 attempts. All reasons were "could not place <role>":
  a fixed gene had been assigned to every role before checking that it can host the variant.
- Changed role assignment to `host()`: walk through the unused panel genes and take the first one that
  can host the variant (all-or-nothing planting functions). Re-measured over seeds 1-40 at levels 3 and 1:
  79 of 80 succeed on attempt 0, one on attempt 1.
- Full generation including the docker guards for seeds 1-12 (level 3), seeds 1-2 (level 1) and seeds 3-4
  (level 2): 16 of 16 on attempt 0, i.e. no csq disagreement, no gray causal/decoy, SOP-on-csq == planted
  set every time. Gray zone 4-7 alleles per instance, always hg19-only alleles in small chain gaps.
- This changed the RNG stream, so every instance differs from the first validation round; all three
  `make` runs were repeated.

## 9. Verifier hardening
- `_scoring.Ctx.text()` decompresses a .gz completely before checking its size. Replaced it for the VCF
  by a streaming, capped decompression of the multi-member BGZF (300 MB cap) inside score.py. Scores of
  the kept submissions are unchanged (checked before the final runs).

## 10. Second validation round (adaptive role assignment, verifier before hardening)
- s1-l3 with Harbor ok (4.49 min), s2-l3 ok (2.62 min), s1-l1 ok (0.41 min), s1-l2 `--no-harbor` ok (0.75 min).
  Log of the level-2 run kept as `results/trio_s1_l2.pre_hardening.log`.
- Regression check of the hardened score.py on the 18 kept submissions of s1-l3: identical rewards.
  Hostile inputs: truncated BGZF -> 0 ("truncated"), 400 MB gzip bomb -> 0 after 0.4 s / 308 MB RSS with a
  300 MB cap (cap then lowered to 50 MB), fake gzip magic -> 0.

## 11. Final validation round (final code)
- `make --seed 1 --level 3 --keep-work`: all gates ok, G6 oracle 1.0 / nop 0.0, 4.71 min.
- `make --seed 2 --level 3 --no-harbor`: all gates ok, 4.0 min (host under load, image rebuilt).
- `make --seed 1 --level 1 --no-harbor`: all gates ok, 2.89 min.
- Free disk dropped to 4.6 GB during this round. I stopped my batch before the optional level-2 repeat,
  removed the images of finished tasks, and deleted the level-2 task directory / validation JSON so that
  no task validated with the older verifier is left behind.
- Both required level-3 seeds drew the same configuration (3 causal variants, compound-het missense as
  2nd ALT), so I added `make --seed 3 --level 3 --no-harbor` (4 causal variants, de novo as 2nd ALT):
  all gates ok, 3.29 min; naive_no_normalise 0.20, quality/phase naive 0.55.
- Cleanup: all `alegen-trio-liftover-*` images removed (rebuildable from the task directories); kept the
  tools image `alegen-tools-5886c792c5f1:local` (1.95 GB), which generation needs; scratch files removed.
- Raw metric values: table in docs/families/trio_liftover.md.

---

# Level 4

All entries 2026-11-12 UTC (single session, 00:38 - 02:00). Times are wall clock from `date -u`.

## L4-1. Reading and state of the build host (00:38)
- Read AUTHORING.md "Level 4", docs/specs/l4/trio_liftover.md, the level-3 notes and code, README section 6.
- No `alegen-*` image of this family was left (cleaned up at the end of level 3). Rebuilt the generation tools
  image `alegen-tools-5886c792c5f1:local` through `dockerlib.tool_image` (87 s, 1.95 GB); the agent images share
  its layers.
- Free disk 26 GB at the start, 16 GB ten minutes later, 7.4-8.0 GB from 01:19 to 01:47, 13 GB at 01:53 (other
  authors building and cleaning up).

## L4-2. What the real sequences offer (00:45)
- Counted the positions inside aligned main-chain blocks of the five windows where hg19 and GRCh38 carry different
  bases: **173** (9 / 13 / 98 / 52 / 1 per window), one pair 2 bp apart, one pair 9 bp apart, 13 within 100 bp of a
  block edge, none within 5 bp, 7 inside a CDS. Enough for the call-set side of the reference-difference element
  without planting; planting is only needed where a difference has to sit in a chosen codon.
- The two minus-strand chains of chr22 (ids 135 and 227) lie outside the windows: an inverted block has to be made.

## L4-3. Design decisions before coding
- Truth stays "family alleles on GRCh38"; at a differing base the provider record is derived from the allele
  *pair* of each member, so sites without a record, swapped records, vanishing records and three-allele sites all
  come out of one function (`rd_site`).
- Genotype in a one-ALT-per-record file had to be defined (copies of that ALT; other alleles count as 0), and the
  inheritance models and the allele-balance rule had to be restated in copies against GRCh38. Without that a
  three-allele site has no single answer.
- Haplotype rule restricted to SNVs in one codon on one copy, with an explicit copy-assignment rule. Brute force
  over the genetic code: no two individually synonymous changes give a stop together (the spec's second example
  does not exist); 18 "silent + missense -> stop" combinations; only four "both non-synonymous -> silent together"
  combinations (Ser TCT/TCC <-> AGT/AGC).
- Homozygous partners: a position without a record is homozygous in the child, hence in cis with a de novo
  variant whatever its copy. That links the two elements the coordinator cares most about (and is the C element).

## L4-4. First level-4 generation and tool behaviour (00:56 - 01:01)
- First end-to-end generation of seed 1 worked on attempt 0 (6.7 s before the docker guard, 14 s in total).
- bcftools +liftover on the level-4 records: `Lines total/swapped/reference added/rejected: 3008/72/26/24` on input
  that had been split with `bcftools norm -m-any` on the provider's reference first. **Splitting first is wrong at
  three-allele sites**: a sample that is GRCh38-ref/third-allele becomes `1/0` in the record "provider base > GRCh38
  base", and after the swap it carries a copy of the provider's base that it does not have. Lifting the records as
  delivered and splitting on GRCh38 gives the planted genotypes. The level-4 oracle therefore does not split
  before the lift-over, and `split_before_liftover` became a measured careless variant.
- For a new reference allele +liftover writes `.` as its AD; the oracle treats it as 0 reads.
- Picard LiftoverVcf (defaults) rejects every record whose REF is not the GRCh38 base: on the final seed-1 instance
  the level-3 method reports 97 such alleles as unmapped next to the 26 scored alleles in unaligned sequence
  (unmapped F1 0.349).
- Oracle (bcftools +liftover + added sites + per-copy `bcftools csq`) scored 1.000 on the first run, the
  pure-Python second solution 1.000 on the first run, i.e. three implementations of rule 5 (generator, csq-based,
  own codon logic) agreed on the first attempt.
- First careless variants (10 candidates, weights 0.10/0.15/0.10/0.45/0.20): l3_reference 0.22, references
  ignored 0.22, records only 0.53, per-variant 0.53, codon merge ignoring phase 0.63, provider genotypes 0.49.

## L4-5. First `make` (01:02 - 01:06)
- `make --seed 1 --level 4 --keep-work` on that version: all gates ok including G6 (oracle 1.0, nop 0.0), 3.66 min.
  Superseded by the versions below; kept only as evidence that the two core elements validated before anything
  else was added.
- Level 3 seed 1 regenerated byte-identical (`diff -r`, exit 0).

## L4-6. Three-allele candidate (01:07)
- Estimated that a Picard `RECOVER_SWAPPED_REF_ALT` route plus added sites would lose only the three-allele sites
  and end near the pass mark. Planted a de novo third allele on a planted difference (provider: parents 1/1, child
  1/2, REF = a base nobody carries, which would itself be a rare qualifying allele), so that mishandling
  three-allele sites changes the candidate list. Measured `split_before_liftover`: 0.777-0.812 on seeds 1-10 with
  this 11-candidate version, 0.798-0.839 on seeds 1-8 with the final version (after the split the record of the
  third allele has no read count for the GRCh38 base, the de novo then fails the allele-balance rule and is lost;
  18 duplicate records appear at three-allele sites).
- Real differences within 30 bp of an exon are now always "no record" sites, so that no SOP rule can depend on how
  reads of a third allele are counted.

## L4-7. Robustness run without `make` (01:09 - 01:14)
- Seeds 3-10 generated into scratch space and scored inside the tools image: 8 of 8 on generation attempt 0,
  oracle 1.000 and second solution 1.000 on all eight, l3_reference 0.209 on all eight (11-candidate version).

## L4-8. Inverted block (01:09 - 01:17)
- Prototype outside the package: reverse-complemented a 4.3-kb segment of the GRCh37 slice, split the main chain
  block in both chain files and appended a minus-strand chain; 73 test alleles (39 indels) written on the other
  strand and left-aligned there. Picard and bcftools +liftover both returned all 73 planted GRCh38 alleles with
  genotypes (after `bcftools norm`), 0 rejects. So the element can be planted with exact truth and both guard
  tools agree; it is a trap for a home-made lift-over, not for a tool user.
- Integrated: segment chosen around an exon (>= 150 bp) of an unused panel gene, 80-bp exclusion band at both
  ends, compound-het pair (frameshift + missense) inside the exon, 32-43 background alleles inside the block, no
  multi-allelic record and no planted difference inside. Third call-set category in the verifier
  (0.08/0.12/0.07/0.08), candidate ramp 0.87-1.0 for 13 candidates.
- First run of the integrated version (seed 1): oracle 1.000, second solution 1.000 (after adding strand handling
  and the minus-strand comparison to both), l3_reference 0.273.
- Sanity test of the edited *reverse* chain (not used by any solution): the 46 inverted truth alleles of seed 1
  lifted GRCh38 -> provider reference with +liftover and the shipped `hg38ToHg19` file all land on provider
  records (46 of 46).

## L4-9. A gate failure and its cause (01:18 - 01:22)
- `make --seed 1 --level 4` on the inverted-block version: G1-G5, G7 ok, **G6 failed** (`harbor run`:
  "ValueError: Either datasets or tasks must be provided", no reward for oracle and nop). Cause: I had put a nested
  dict (`inverted_block = {...}`) into `extra_metadata`; `harbor.write_task` serialises values with `json.dumps`,
  which is not a TOML inline table, so Harbor could not load `task.toml`. Flattened the metadata to scalars and
  lists. The core was not touched; a dict value in `extra_metadata` is a trap for other authors too.
- The generator now deletes the 51-MB plain FASTA it writes for the tools image after the guard step (ten scratch
  generations had used 600 MB).

## L4-10. Final validation
- `make --seed 1 --level 4 --keep-work` (final code, 01:24 - 01:29): all 35 checks ok, G6 oracle 1.0 / nop 0.0,
  5.01 min, oracle 0.09 min.
- Levels 1-3 after all edits: level 3 seed 1 `diff -r` identical to `tasks/trio-liftover-s0001-l3`; level 3 seeds
  2 and 3 and level 1 seed 1 regenerate the `truth.json`, `trio.vcf` and `stats.json` kept in `results/gen/`.
- `make --seed 2 --level 4 --no-harbor --keep-work`: the build waited from 01:32 to 01:47 for the host to have
  8 GB free again (a background loop checked `df` every 20 s and started `make` at 8.9 GB). All 33 checks ok,
  6.11 min, oracle 0.11 min (`results/trio_s2_l4.log`). Seeds 1 and 2 drew different level-3 configurations this
  time (causal SNV as 2nd ALT: compound-het missense / de novo).
- While waiting I ran the same checks for seed 2 without building images (`results/trio_l4_emu_validate.py`: oracle,
  determinism, probes and attacks in the generation tools image with the input mounted, scoring in the host venv;
  log `results/trio_s2_l4.no_image_build.log`). 31 checks ok, every score identical to the later `make` run.
- Robustness with the final code, seeds 3-8 in scratch space: 5 of 6 on generation attempt 0, seed 8 on attempt 1
  ("too few alleles fit into the inverted block"); oracle 1.000 and second solution 1.000 on all six; ranges of
  the careless variants in the design notes.
- Level-4 generation is byte-deterministic: seed 1 regenerated into a scratch directory is `diff -r`-identical to
  `tasks/trio-liftover-s0001-l4`.
- Outside the battery: Picard with `RECOVER_SWAPPED_REF_ALT=true` + added sites + rule 5 scores 0.702 / 0.708 on
  seeds 1 / 2 (three-allele records rejected, three-allele de novo lost). Verifier checked by hand on equivalent
  encodings of the oracle output (phased, other allele order, GT-only, non-carrier records, gray positions as one
  record or none: all 1.000) and on a truncated BGZF file (0, gate).
- Found while reading the split-before-liftover details: duplicate records are counted as *standard* false
  positives whatever their position (inherited from the level-3 verifier). Not changed, because it does not move
  any credit in the measured cases and a change would have meant a third full validation of seed 1; listed as a
  limitation.

## L4-11. Cleanup
- Removed the seed-2 images (`alegen-trio-liftover-s0002-l4-agent/-verifier`) after seed 2 had passed, and the
  seed-2 task directory (only seed 1 is kept under `tasks/`, as at level 3; `results/validation/trio-liftover-s0002-l4.json`
  and `results/work/trio-liftover-s0002-l4/` are kept).
- Kept: `alegen-trio-liftover-s0001-l4-agent:local` (1.96 GB) and `-verifier:local` (178 MB) and the generation tools
  image `alegen-tools-5886c792c5f1:local` (1.95 GB); the agent image shares the 1.95-GB conda layer with the tools
  image, so the three together occupy about 2 GB.
- Scratch directory emptied of the large files; nothing of other authors was touched; no `docker ... prune`.
- My hand-run verifier tests (01:33) had executed `tasks/trio-liftover-s0001-l4/tests/score.py` with the host
  Python, which left `tests/__pycache__/_scoring.cpython-313.pyc` inside the shipped task directory (after its
  images had been built and validated). Removed it at 01:57 and re-checked: level 4 seed 1 and level 3 seed 1
  regenerated into scratch are both `diff -r`-identical to the directories under `tasks/`.

## L4-12. What I did not do
- No agent pilot (by instruction).
- No third seed through `make` (seeds 3-8 only through the scratch route above).
- No change to the core. Two things I worked around or noticed there: `harbor.write_task` writes `extra_metadata`
  values with `json.dumps` (a dict value yields an invalid `task.toml`, see L4-9), and `validate.py` has no way to
  run the gates without building images (hence the small driver in `results/`).
