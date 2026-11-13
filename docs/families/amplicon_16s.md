# Family `amplicon_16s`: 16S V4 amplicon reads to exact sequence variants, abundances, taxonomy and reagent contaminants

ALE subdomain 3.2 (Genomics & Sequence Analysis, microbiome). Code:
`pipeline/alegen/families/amplicon_16s/`. Spec: `docs/specs/amplicon_16s.md`.

## Scenario

A microbiome core receives demultiplexed paired-end 2 x 250 reads of 11 libraries: 8 study
specimens, 2 extraction blanks and one mock community standard. Deliverables: `asv.fasta`
(primer-free exact variants in the 515F -> 806R direction), `asv_table.tsv` (read pairs per ASV and
library), `taxonomy.tsv` (genus by the SOP rule), `contaminants.txt` (reagent contaminants).
Tools in the image: `vsearch=2.30.1`, `cutadapt=5.2`, Python 3.12 with numpy/scipy/pandas.

Inputs (13 MB): `fastq/<id>_R{1,2}.fastq.gz`, `samples.tsv` (type, specimen, qPCR 16S copies),
`primers.tsv`, `reference_db.fasta` + `reference_taxonomy.tsv` (V4-trimmed reference set, ids
blinded), `mock_community.tsv` (genus and nominal fraction of the 8 mock strains, no sequences),
`library_prep_notes.md`, `sop.md`.

## Reference material

| file | URL | sha256 |
|---|---|---|
| `.cache/amplicon_16s/bacteria.16SrRNA.fna.gz` | https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz (downloaded 2026-11-11, 26,888 records) | `c5dc51b877cfed87b526b5cb9425f05d59ce96805ce4d4eb2015d65d75436e7e` |

The generator refuses to run if the cached file has another checksum (NCBI updates this file in
place; a different release would silently change every instance). Derived caches in the same
directory: `v4_pool.tsv` (13,183 unique amplicons of 251-255 nt that occur in exactly one genus,
cut out by exact IUPAC match of `GTGYCAGCMGCCGCGGTAA` / `GGACTACNVGGGTWTCTAAT`) and
`v4_close_pairs.json` (same-genus pairs at Hamming distance 1 or 2).

## Truth model (what is planted and why it is correct)

Everything scored is decided before a read exists and stored in `tests/reference/truth.json`.

1. **Planted sequences.** Level 3: 30 specimen strains (5 real same-genus pairs at Hamming distance
   1 or 2, 20 singles), 3 rare specimen strains, 5 reagent contaminants, 8 mock strains (one pair
   at distance 1-2). Every unit comes from a different genus and all non-pair sequences are
   provably more than 16 edits apart when they are selected (8-mer pigeonhole bound), so nothing
   but the designated pairs is close. 6 specimen singles then receive 1-2 substitutions (a strain
   that is not in the reference set); the provable separation after that is more than 12 edits. No planted sequence is an exact prefix+suffix combination of two others (checked), so a
   bimera test cannot legitimately remove one. Which strains are contaminants is random: genus
   names carry no information about it.
2. **Read-origin counts.** Per library the generator draws the composition (specimen strains
   x (1 - c) + reagent strains x c, with c = 0.4-2 % for normal specimens and the mock, 25-60 % for
   three low-biomass specimens, 100 % for the blanks), the number of non-recombinant pairs per
   strain (multinomial), ~70 cross-talk pairs drawn from the other libraries of the run, and 5-8 %
   recombinant pairs. `counts[strain][library]` = own + cross-talk pairs + recombinant pairs whose
   sequence happens to be identical to a planted sequence (as almost every recombinant of the two
   members of a distance-1 pair is). Other recombinant pairs belong to nobody.
3. **Genus.** The genus of the RefSeq record, or `Unclassified` for two specimen strains for which
   every reference within a provable 22 edits was removed from the shipped set. The reference set
   is the planted sources (except those two) plus up to 9,000 random pool sequences, minus every
   other-genus sequence that might be within 15 edits of a planted sequence. Consequence, proved
   by the k-mer bound and independent of any aligner: a classified planted sequence has a
   same-genus reference at >= 99.2 % and no other-genus reference above 95.3 %; an `Unclassified`
   one has no reference above 95.3 % (no distractor above 92 %).
4. **Scored subset (the "what denoising can resolve" rule).** A planted sequence is scored for
   recall, genus and contaminant status iff the batch contains >= 100 read pairs of it (the counts
   of point 2). The three rare strains (5-35 pairs in one or two libraries) are a gray zone:
   reporting them is not a false positive, missing them is not a miss. Measured gap around the
   limit: the smallest scored strain has 135-246 read pairs, the largest gray one 27 (level 3,
   seeds 1-4). Both solutions recover every scored strain on seeds 1 and 2; of the gray ones the
   oracle finds 1 of 3 (seed 1) and 3 of 3 (seed 2), the second solution 2 of 3 and 3 of 3. Pair members stay within a factor 3
   of each other in the pooled data (asserted), so pooled UNOISE3 (skew limit 1/8 at distance 1,
   1/32 at distance 2) can separate them.

Read simulator: construct = spacer(0-3 nt, fixed per primer version) + 515F instance + amplicon +
rc(806R instance) + rc(spacer); degenerate primer positions are instantiated at random; 42-58 % of
the pairs are read from the 806R end (ligation library); three quality classes (Q37/Q25/Q12) whose
frequency rises towards the 3' end (faster in R2) and substitution errors at exactly the nominal
rate of the class; recombinants switch template inside a shared 12-mer.

## Traps and how each can be discovered

| trap (level) | what goes wrong | discoverable by |
|---|---|---|
| `primer_spacers` (2+) | fixed-length trimming or no trimming leaves primer/spacer bases: no reported sequence is a real amplicon | notes: "four phased versions per primer, carrying 0, 1, 2 or 3 extra nucleotides"; first 25 bases of the reads; `primers.tsv` has IUPAC codes |
| `mixed_orientation` (2+) | half of the pairs start with 806R on R1 | notes: adapters added by ligation; cutadapt reports ~50 % untrimmed pairs; R1 heads |
| `chimeras` (2+) | 101 of 145 UNOISE3 centroids (seed 1; seed 2: 103 of 149) are recombinants | the mock library yields more sequences than the standard has strains: without chimera removal 30 (seed 1) / 24 (seed 2) sequences have >= 8 read pairs in the mock library, 15 / 11 of them not planted; they are exact two-parent combinations of more abundant sequences |
| `strain_variants` (3) | 97 % OTUs (or coarse error merging) report one member per pair; the minor member dominates in two specimens | `sop.md` item 1; notes: "the two <genus> strains of the standard differ at 1 position(s) of the V4 amplicon", so the mock shows whether a workflow resolves them |
| `reagent_contaminants` (3) | five reagent taxa are 25-60 % of the low-biomass specimens; abundant specimen taxa also occur in the blanks (cross-talk), so "remove everything seen in a blank" is wrong too | blanks in `samples.tsv`; reagent taxa have their highest relative abundance in the blanks and rise as `qpcr_16s_copies_per_ul` falls. Measured on the planted counts of seed 1: relative abundance in the blanks divided by pooled relative abundance in the other libraries is 4.9-6.5 for the reagent strains and at most 0.27 for specimen and mock strains (31 of which have at least one read pair in a blank); relative abundance x qPCR copies, blanks over other libraries, is 0.70-1.04 for reagent strains and at most 0.004 for the others (seeds 1 and 2) |

Level 1 has none of them (primers without spacers, one orientation, no recombinants, no pairs, no
reagent background; the blanks contain only ~70 cross-talk pairs). Level 2 adds the first three.

## Output contract

See `instruction.md` of a generated task. Points that matter for scoring: sequences are compared
exactly (upper-cased); ids are free; `asv_table.tsv` and `taxonomy.tsv` must list exactly the ASVs
of `asv.fasta`; `contaminants.txt` has a header line `asv_id` (accepted with or without it) so that
"no contaminant" is not an empty file; contaminant ASVs stay in the table and are removed by the
verifier before the specimen profiles are compared (stated in the instruction).

## Scoring

Hard gates: files present/parseable, unique ids and sequences, ACGT only, table/taxonomy list
exactly the ASVs, integer counts, known ids in `contaminants.txt`.

| weight | component | raw metric | full | zero |
|---|---|---|---|---|
| 0.35 | `asv_sequences` | recall of scored planted sequences x precision of reported sequences (exact) | recall >= 0.95, precision >= 0.90 | recall <= 0.85 or precision <= 0.60 |
| 0.30 | `sample_counts` | min over study libraries of `sum min(obs,true) / max(sum obs, sum true)` after removing listed contaminants | >= 0.80 | <= 0.60 |
| 0.15 | `contaminants` | MCC of the flags over scored planted sequences | >= 0.85 | <= 0.50 |
| 0.10 | `genus` | fraction of scored planted sequences with the planted genus | >= 0.95 | <= 0.80 |
| 0.10 | `mock_counts` | the count metric on the mock library | >= 0.80 | <= 0.60 |

Pass threshold 0.9.

## Measured values

All numbers below are from `python -m alegen make` (agent image, verifier image), final code,
2026-11-11. Columns: score | recall | precision | min study count agreement | mock count
agreement | contaminant MCC | genus accuracy | number of reported sequences.

### Level 3

| submission | seed 1 | seed 2 |
|---|---|---|
| oracle (`amplicon_pipeline.py`) | 1.000 / 1.000 / 1.000 / 0.9768 / 0.9907 / 1.00 / 1.000 / 44 | 1.000 / 1.000 / 1.000 / 0.9788 / 0.9875 / 1.00 / 1.000 / 46 |
| independent solution (`solve_alt.py`) | 1.000 / 1.000 / 1.000 / 0.9776 / 0.9906 / 1.00 / 1.000 / 45 | 1.000 / 1.000 / 1.000 / 0.9786 / 0.9874 / 1.00 / 1.000 / 46 |
| probe `stock_unoise3_uchime3` | 0.977 / 1.000 / 0.880 / 0.9750 / 0.9907 / 1.00 / 1.000 / 50 | 1.000 / 1.000 / 1.000 / 0.9788 / 0.9875 / 1.00 / 1.000 / 46 |
| probe `planted_truth` | 1.000 / 1.000 / 1.000 / 1.0000 / 1.0000 / 1.00 / 1.000 / 46 | 1.000 / 1.000 / 1.000 / 1.0000 / 1.0000 / 1.00 / 1.000 / 46 |
| probe `within_tolerance` | 1.000 / 0.977 / 0.978 / 0.8836 / 0.9028 / 1.00 / 0.953 / 46 | 1.000 / 0.977 / 0.978 / 0.8845 / 0.9045 / 1.00 / 0.953 / 46 |
| `naive_no_primer_removal` (trap primer_spacers) | 0.000 / 0.000 / 0.000 / 0.0000 / 0.0000 / 0.00 / 0.000 / 66 | 0.000 / 0.000 / 0.000 / 0.0000 / 0.0000 / 0.00 / 0.000 / 5 |
| `naive_fixed_length_trim` (trap primer_spacers) | 0.000 / 0.070 / 0.035 / 0.0000 / 0.1000 / 0.43 / 0.070 / 85 | 0.000 / 0.047 / 0.023 / 0.0000 / 0.0729 / 0.00 / 0.047 / 86 |
| `naive_single_orientation` (trap mixed_orientation, untrimmed pairs kept) | 0.250 / 1.000 / 0.468 / 0.4485 / 0.4386 / 1.00 / 1.000 / 94 | 0.600 / 1.000 / 0.979 / 0.5123 / 0.5253 / 1.00 / 1.000 / 47 |
| `naive_single_orientation_discard` (trap mixed_orientation, untrimmed pairs dropped) | 0.600 / 1.000 / 1.000 / 0.4485 / 0.4384 / 1.00 / 1.000 / 44 | 0.600 / 1.000 / 1.000 / 0.5123 / 0.5253 / 1.00 / 1.000 / 46 |
| `naive_no_chimera_removal` (trap chimeras) | 0.650 / 1.000 / 0.303 / 0.9446 / 0.9734 / 1.00 / 1.000 / 145 | 0.650 / 1.000 / 0.309 / 0.9322 / 0.9716 / 1.00 / 1.000 / 149 |
| `naive_otu97` (trap strain_variants) | 0.327 / 0.860 / 1.000 / 0.4388 / 0.8571 / 1.00 / 0.860 / 40 | 0.327 / 0.860 / 1.000 / 0.4720 / 0.8759 / 1.00 / 0.860 / 40 |
| `minor_variants_merged` (truth-derived, trap strain_variants) | 0.327 / 0.860 / 1.000 / 0.4608 / 0.8653 / 1.00 / 0.860 / 40 | 0.327 / 0.860 / 1.000 / 0.4819 / 0.8870 / 1.00 / 0.860 / 40 |
| `naive_no_contaminant_handling` (trap reagent_contaminants) | 0.550 / 1.000 / 1.000 / 0.4068 / 0.9776 / 0.00 / 1.000 / 44 | 0.550 / 1.000 / 1.000 / 0.4820 / 0.9747 / 0.00 / 1.000 / 46 |
| `naive_everything_in_blanks` (trap reagent_contaminants) | 0.450 / 1.000 / 1.000 / 0.0197 / 0.2680 / 0.16 / 1.000 / 44 | 0.450 / 1.000 / 1.000 / 0.0250 / 0.0003 / 0.12 / 1.000 / 46 |
| `naive_default_pipeline` (all traps ignored) | 0.040 / 0.860 / 0.175 / 0.1779 / 0.4253 / 0.00 / 0.860 / 223 | 0.040 / 0.860 / 0.135 / 0.1353 / 0.5035 / 0.00 / 0.860 / 297 |
| `half_of_the_read_pairs` (beyond tolerance) | 0.600 / 1.000 / 1.000 / 0.4982 / 0.5043 / 1.00 / 1.000 / 46 | 0.600 / 1.000 / 1.000 / 0.4905 / 0.5091 / 1.00 / 1.000 / 46 |
| `closed_reference_best_case` | 0.559 / 0.814 / 1.000 / 0.8015 / 0.9993 / 1.00 / 0.814 / 38 | 0.483 / 0.814 / 1.000 / 0.7492 / 0.9986 / 1.00 / 0.814 / 38 |
| `reference_as_asvs` (prior only) | 0.000 / 0.047 / 0.050 / 0.0000 / 0.0000 / 0.00 / 0.047 / 40 | 0.000 / 0.000 / 0.000 / 0.0000 / 0.0000 / 0.00 / 0.000 / 40 |
| `empty_output`, `inputs_copied`, `zero_byte_files`, `header_only` | 0.000 (gate) | 0.000 (gate) |

Oracle runtime in the agent image: 0.2-0.6 min (2 threads, shared host). Generation: 6-9 s.
Inputs: 13 MB. Whole validation: 3.7-8.1 min per level-3 instance (shared host).
The gate runs were repeated after the last code change; the metric tables of the repeated runs
are identical to the ones before it (compared with `diff`).

### Level 1 (seed 1) and level 2 (seed 1)

| submission | level 1, seed 1 | level 2, seed 1 |
|---|---|---|
| oracle | 1.000 / 1.000 / 1.000 / 0.9991 / 0.9994 / n.a. / 1.000 / 35 | 1.000 / 1.000 / 1.000 / 0.9786 / 0.9873 / n.a. / 1.000 / 36 |
| independent solution | 1.000 / 1.000 / 1.000 / 0.9993 / 0.9994 / n.a. / 1.000 / 36 | 1.000 / 1.000 / 1.000 / 0.9799 / 0.9901 / n.a. / 1.000 / 36 |
| probe `stock_unoise3_uchime3` | not run (no recombinants) | 1.000 / 1.000 / 0.947 / 0.9779 / 0.9873 / n.a. / 1.000 / 38 |
| probe `within_tolerance` | 0.993 / 0.970 / 0.972 / 0.8751 / 0.9060 / n.a. / 0.939 / 36 | 0.993 / 0.970 / 0.972 / 0.8828 / 0.9048 / n.a. / 0.939 / 36 |
| `naive_no_primer_removal` | not a trap at level 1 | 0.000 / 0.000 / 0.000 / 0.0000 / 0.0000 / n.a. / 0.000 / 57 |
| `naive_fixed_length_trim` | not a trap at level 1 | 0.000 / 0.061 / 0.030 / 0.0000 / 0.0816 / n.a. / 0.061 / 66 |
| `naive_single_orientation` | not a trap at level 1 | 0.438 / 1.000 / 0.761 / 0.5646 / 0.5699 / n.a. / 1.000 / 46 |
| `naive_single_orientation_discard` | not a trap at level 1 | 0.600 / 1.000 / 1.000 / 0.5646 / 0.5698 / n.a. / 1.000 / 35 |
| `naive_no_chimera_removal` | not a trap at level 1 | 0.650 / 1.000 / 0.221 / 0.9523 / 0.9654 / n.a. / 1.000 / 163 |
| `naive_everything_in_blanks` | 0.450 / 1.000 / 1.000 / 0.0006 / 0.0000 / n.a. / 1.000 / 35 | 0.450 / 1.000 / 1.000 / 0.0000 / 0.0000 / n.a. / 1.000 / 36 |
| `half_of_the_read_pairs` | 0.600 / 1.000 / 1.000 / 0.4859 / 0.5024 / n.a. / 1.000 / 36 | 0.600 / 1.000 / 1.000 / 0.4828 / 0.5061 / n.a. / 1.000 / 36 |
| `closed_reference_best_case` | 0.100 / 0.758 / 1.000 / 0.4503 / 0.9974 / n.a. / 0.758 / 28 | 0.100 / 0.758 / 1.000 / 0.4316 / 0.9970 / n.a. / 0.758 / 28 |
| `reference_as_asvs` | 0.000 / 0.000 / 0.000 / 0.0000 / 0.0000 / n.a. / 0.000 / 40 | 0.000 / 0.000 / 0.000 / 0.0000 / 0.0000 / n.a. / 0.000 / 40 |

At levels 1 and 2 there is no reagent background, so there is no MCC. The contaminant component is
then the fraction of scored planted sequences that are reported and not flagged: 1.000 for oracle
and second solution, 0.970 for `within_tolerance`, 0.758 for `closed_reference_best_case`, 0.030 for
`naive_everything_in_blanks`, 0.000 for `reference_as_asvs`. (The first level-1 run failed G4
because this component used to give 1.0 to any submission that flagged nothing; see the worklog.)

### Why the bands sit where they do

| raw metric | correct workflows (oracle, second solution; seeds 1, 2) | stock UNOISE3 + uchime3 | closest wrong workflow | band (zero -> full) |
|---|---|---|---|---|
| recall of scored planted sequences | 1.000 | 1.000 | 0.8605 (97 % OTUs: 37 of 43, the six minor pair members are gone) | 0.85 -> 0.95 |
| precision of reported sequences | 1.000 | 0.88 (seed 1), 1.00 (seed 2) | 0.30 / 0.31 (no chimera removal), 0.47 (single orientation, seed 1) | 0.60 -> 0.90 |
| min study count agreement | 0.977 - 0.979 (level 1: 0.999) | 0.975 / 0.979 | 0.449 / 0.512, level 2: 0.565 (one orientation; bounded by the share of forward pairs <= 0.58, so the worst case is only about 0.02 below the zero point), 0.407 / 0.482 (contaminants kept; bounded by 1 - 0.45), 0.439 / 0.472 (OTUs) | 0.60 -> 0.80 |
| mock count agreement | 0.987 - 0.991 | 0.988 / 0.991 | 0.438 / 0.525 (one orientation) | 0.60 -> 0.80 |
| contaminant MCC | 1.000 | 1.000 | 0.00 (none listed), 0.16 / 0.12 (everything seen in a blank) | 0.50 -> 0.85 |
| genus accuracy | 1.000 | 1.000 | 0.8605 (OTUs), 0.047 (reference sequences as ASVs) | 0.80 -> 0.95 |

- One missed scored sequence (42/43 = 0.977), one extra non-planted sequence and one wrong genus
  (41/43 = 0.953) stay at full credit (`within_tolerance` probe: 1.000 on both seeds, 0.993 at
  level 1 where only 33 sequences are scored).
- One wrong contaminant flag among 43 scored sequences (5 true) gives MCC 0.88 (one missed) or 0.90
  (one extra): full credit. (Computed from the formula, not run.)
- Counting only error-free read pairs (no assignment of reads with errors) was measured in the
  development harness (tools image, verifier on the host) with `amplicon_pipeline.py --counts exact`: min study count agreement 0.917 / 0.914, mock 0.920 / 0.923 (seeds 1 / 2), i.e.
  still full credit and score 1.000. The convention "reads with errors count" is therefore not what
  separates pass from fail; losing an orientation or keeping the reagent taxa is.
- No naive variant exceeds 0.65 (ceiling 0.75 = pass threshold - 0.15); both correct solutions are
  at 1.000 on seeds 1 and 2 of level 3 and on seed 1 of levels 1 and 2.

## Deviations from the spec (and why)

1. **Reference database is V4-trimmed** (primer-free amplicon region, ids blinded, genus in a
   separate TSV) instead of full-length records. This is how 515F/806R-trimmed classifiers are
   shipped in practice, keeps the input small and lets the genus truth be proved with a k-mer
   bound instead of an aligner run.
2. **46 planted sequences** (43 scored) instead of "~40".
3. **Sequence recovery**: recall is unweighted over the scored planted sequences (not
   abundance-weighted) and precision counts sequences, not reads. Abundance weighting would make
   the minor strain variants and the recombinants (many sequences, few reads) almost free, which
   are exactly the two things the component is meant to see.
4. **"Per-sample composition distance" became a count agreement** (`sum min / max`, minimum over
   study libraries, contaminants removed first). A distance between relative profiles cannot see a
   workflow that drops one read orientation (measured: recall and precision stay 1.0, only the
   counts halve), and it would not connect the contaminant list with the specimen profiles.
5. **Mock component** uses the same count metric on the mock library against the emitted read
   pairs, not against the nominal fractions (the real fractions deviate from nominal by design,
   sd 15 %).
6. **Extra inputs**: `mock_community.tsv`, `reference_taxonomy.tsv`, qPCR copies in the sample
   sheet (needed for the discoverability of the reagent background).
7. **Extra battery entries** beyond the five naive variants of the spec: fixed-length primer
   trimming, single orientation with discarded pairs, "everything seen in a blank", the fully naive
   pipeline, the truth-derived `minor_variants_merged`, `half_of_the_read_pairs` and
   `closed_reference_best_case`, and the `stock_unoise3_uchime3` probe (floor 0.9, see below).
8. **Level 2** was not among the requested runs; it was run once in addition
   (`make amplicon_16s --seed 1 --level 2 --no-harbor`: all gates ok, numbers above).

## Known limitations and things I am unsure about

- **The community is cleaner than a real one.** Apart from the designated pairs every planted
  sequence is more than 12 edits (< 95.3 % identity) from every other. Real gut communities contain
  many taxa at 95-99 %; denoising, chimera detection and genus assignment are easier here.
- **Error model.** Substitutions only (no indels), three quality classes drawn independently per
  base, errors exactly at the nominal rate. Real quality scores are autocorrelated and
  miscalibrated; expected-error filtering works better here than on real runs.
- **Chimera model.** Two parents, switch point uniform over shared 12-mers, no errors-on-chimeras
  beyond the read error model. `--uchime_denovo` removed every recombinant and no planted sequence
  on both seeds; `--uchime3_denovo` left 6 (seed 1), 0 (seed 2) and 2 (level 2, seed 1). The three
  leftovers I inspected (an earlier draw of seed 1) were exact two-parent recombinants with a
  37-38 nt tail of the second parent, 4-5 substitutions from the main parent. The oracle uses `--uchime_denovo` because of this measurement on seed 1; I did not
  investigate why uchime3 misses them (parameter changes `--abskew`, `--mindiffs`, `--mindiv`,
  `--minh` made no difference).
- **Precision band history.** The band started at [0.70, 0.95]. After seeing the stock recipe at
  precision 0.94 (3 leftovers, an earlier draw of seed 1) I changed two things at once: the oracle moved to `--uchime_denovo`
  (precision 1.0, so it passes under either band) and the band became [0.60, 0.90]. With the final
  data the stock recipe is at 0.88 on seed 1, i.e. it still loses credit (score 0.977). The
  `stock_unoise3_uchime3` probe has floor 0.9 for that reason; on a seed with more leftovers it
  could drop below 0.9 and G3 would fail. Only two seeds were measured.
- **The second solution is not fully independent**: it uses `vsearch --fastq_mergepairs` like the
  oracle. Primer handling, denoising, chimera test, read assignment, taxonomy and contaminant
  logic are independent code and different methods.
- **Contaminants are easy to separate** once the right comparison is made: each reagent strain is
  at least 8 % of the reagent background and the blanks have 1.5-2.5k read pairs. There is no taxon
  that is both a specimen member and a reagent contaminant.
- **Count component and yield.** With this (too clean) error model exact matches alone account for
  about 92 % of the planted read pairs, so even exact-match-only counting gets full credit
  (measured above). On real
  runs with 60-80 % error-free reads the 0.80 limit would start to matter; here it does not.
  Recombinant pairs are not mentioned in the count definition ("read pairs that derive from no ASV
  are not counted"): the oracle assigns those within 97 % of a parent to that parent, the truth
  counts them for nobody; the oracle's count agreement of 0.977-0.979 includes that difference.
- **What the reference set gives away.** It contains 38 of the 46 planted sequences verbatim among
  3.8k records (planted sources are spread uniformly over the file; positions were checked, the 2
  planted records among the first 40 of seed 1 are chance). `closed_reference_best_case` shows
  that even perfect reference-only picking stays at 0.56 / 0.48 (level 3) and 0.10 (level 1),
  because the 6 strains with substitutions and the 2 `Unclassified` strains are not in the set.
  An agent that combines reference matching with abundance reasoning is doing a legitimate
  closed-reference analysis, not exploiting a leak, but it does shorten the denoising step.
- **Thin lower margin of the count band for one variant.** A workflow that drops one orientation
  reaches a count agreement equal to the share of forward pairs (level 2, seed 1: 0.565); the
  generator draws that share from 42-58 %, so the worst case (about 0.58) is only 0.02 below the
  zero point of the band (0.60). The attack ceilings are not at risk (0.60 total against 0.75), but
  narrowing the draw to 45-55 % would give a cleaner margin. Not changed, because it would have
  redrawn every instance after the measurements were taken.
- **Determinism.** The same seed gives byte-identical task directories (checked by generating
  seed 2 / level 3 twice: 42 files, identical sha256). This needed `mtime=0` in the gzip writer;
  before that fix the FASTQ content was identical but the gzip headers were not.
- **Generator depends on one RefSeq release** (checksum-pinned); if the cache is lost and NCBI has
  updated the file, generation stops with an explicit error.
- **Robustness beyond the tested seeds**: generation was probed on seeds 1-4 (no assertion
  failure); the full battery only on seeds 1 and 2.

## Level 7

Spec: `docs/specs/l7/amplicon_16s.md`. Code: `families/amplicon_16s/l7.py` (generator, dispatched from
`generate()` for `level == 7`), `assets/amplicon_pipeline_l7.py` (oracle), `assets/solve_alt_l7.py` (second
solution), `assets/solve_l7.sh`, `assets/score_l7.py` (verifier: the level-3 verifier run for both filings of
the two exchanged columns). All numbers below are from 2026-11-12, final code, i.e. after the maintainer's
decision of the same day to accept the count table under either filing. Gate runs:
`make amplicon_16s --seed 1 --level 7` (all gates incl. Harbor, 35 checks, 6.5 min) and `--seed 2 --level 7 --no-harbor`
(33 checks, 5.2 min). Raw metrics behind the scores come from a development harness (same scripts in the existing
agent image with the level-7 inputs mounted, the task's `score.py` on the host); its scores agree with the
gate runs for every battery entry of both seeds.

### What is planted

**M_control_swap.** The level-3 batch of the seed is generated unchanged (same strains, read pairs, sample
sheet, texts). Then the FASTQ pair of one extraction blank and of one study specimen are filed under each
other's `sample_id` (same reads, same order; written again so that the name in the gzip header is the name of
the file; read names carry no sample information). `samples.tsv` is the level-3 sheet. Instruction, SOP and
notes are the level-3 texts; `make --level 7` removes the three pointer sentences as for level 6, so the
agent-visible files differ from the level-6 task of the seed in four FASTQ files and nothing else (checked
with `diff -r`).

| seed | labelled blank | holds | labelled specimen | holds |
|---|---|---|---|---|
| 1 | S04 (182 copies/uL) | library of the stool specimen S05: 9,301 read pairs, 2.2 % reagent reads | S05 (9.33e3 copies/uL) | library of blank S04: 1,731 read pairs, 96.4 % reagent reads |
| 2 | S03 (146 copies/uL) | library of the stool specimen S05: 8,547 read pairs, 1.1 % reagent reads | S05 (1.94e4 copies/uL) | library of blank S03: 2,518 read pairs, 97.4 % reagent reads |

Truth (`tests/reference/truth.json`), unchanged in meaning: sequences, genus, contaminant status and the
scored subset are those of level 3; `counts` follow the files (the two columns change places);
`samples[].type` is the kind of template the file holds, with the sheet's entry kept as `sheet_type`.
`sample_counts` is the minimum over the libraries that hold a specimen (seven unaffected specimens and the
file labelled as a blank); the columns of the two libraries that are blanks are not scored, as at level 3.

Verifier: `tests/score_l3.py` is the level-3 `score.py`, byte-identical (checked with `cmp`); every gate,
component, weight and tolerance comes from it. `tests/score.py` (`assets/score_l7.py`) imports it and runs its
`main()` twice, on the table as submitted and with the values of the two exchanged columns swapped (both
together; the pair of ids is read from `truth.json`), and takes the better total. `details.json` holds the
gates and components of the reading that was used and, under `info.count_table_reading`, which one that was
and both totals. See "Either filing of the count table" below.

The pair is drawn from the seed (own random stream) among 2 blanks x the specimens with less than 5 %
reagent reads (10 candidates). A draw is used only if, on the planted counts as filed: (1) the QC below
separates the two libraries (Bray-Curtis labelled blank to other blank >= 0.75, labelled specimen to other
blank <= 0.15, every other specimen to that blank >= 0.33; richness of the labelled blank >= 0.75 x and of
the labelled specimen <= 0.5 x that of the poorest other specimen; reagent share x qPCR copies off by at
least 10 x for both and within 2.5 x of the median for the nine others); (2) the reference contaminant test
with the sheet's roles gives a submission that scores <= 0.75; (3) the same test with the real roles returns
exactly the planted contaminants; (4) so does the test with both libraries left out. On the 12 of the
seeds 1-16 that the level-3 generator can build, all 120 candidate exchanges satisfy all four (first draw
used on seeds 1-6).

### Discoverability (`Trap.discoverable_by`), measured on the planted counts as filed

| QC | seed 1 | seed 2 |
|---|---|---|
| read pairs: labelled blank / other blank | 9,301 / 1,769 | 8,547 / 2,317 |
| read pairs: labelled specimen / other specimens | 1,731 / 7,974-9,195 | 2,518 / 7,652-9,506 |
| ASVs above 0.5 % of the library: labelled blank / other blank / other specimens / labelled specimen | 23 / 6 / 17-26 / 6 | 25 / 5 / 19-29 / 5 |
| share of reads in the five reagent ASVs: labelled blank / other blank | 0.022 / 0.965 | 0.011 / 0.972 |
| the same: labelled specimen / other specimens | 0.964 / 0.007-0.586 | 0.974 / 0.009-0.509 |
| Bray-Curtis, labelled blank to: other blank / nearest other specimen (median) | 0.952 / 0.245 (0.649) | 0.965 / 0.314 (0.562) |
| Bray-Curtis, labelled specimen to other blank; nearest that any other specimen comes to it | 0.080; 0.388 | 0.056; 0.465 |
| reagent share x qPCR copies of the row: nine other libraries (median) | 173-256 (204) | 139-204 (167) |
| the same: labelled specimen / labelled blank | 8,989 (44 x above) / 4.1 (50 x below) | 18,899 (113 x above) / 1.5 (108 x below) |

Over the 120 candidate exchanges: Bray-Curtis labelled blank to other blank 0.95-0.97, labelled specimen to
other blank 0.04-0.08, other specimens to that blank >= 0.38; copies relation off by 37-262 x (labelled
blank) and 43-261 x (labelled specimen), the other libraries within 1.7 x of their median.

The same seen by the two solutions on the reads (seeds 1 / 2; seeds 3-6 in brackets, development harness):

- Oracle (profiles only, no qPCR): Bray-Curtis between the two labelled blanks 0.952 / 0.965 (0.95-0.97;
  limit > 0.50); labelled specimen to the other blank 0.080 / 0.054 (0.04-0.07; limit <= 0.25); nearest any
  other specimen comes to that blank 0.39 / 0.47 (0.40-0.44); ASVs above 0.5 %: labelled blank 23 / 25,
  labelled specimen 6 / 5, half the median of the other specimens 12 / 12.
- Second solution (qPCR mass balance, no sample types, no profile comparison): 5 / 6 ASVs with a constant
  absolute abundance; the two libraries are -1.68 and +1.63 / -2.04 and +2.06 decades off the common amount
  (seeds 3-6: 1.68-2.27; limit 1.0), all other libraries within 0.09 / 0.12 (seeds 3-6: 0.26); after
  exchanging the two rows within 0.08 / 0.09 (0.25).

### What the sheet's roles cost

| | seed 1 | seed 2 |
|---|---|---|
| `l6_reference` (reference pipeline of levels 3-6, unchanged) | **0.55** | **0.55** |
| components lost | `sample_counts` 0 (min agreement 0.317), `contaminants` 0 (MCC -0.099) | `sample_counts` 0 (0.353), `contaminants` 0 (MCC -0.132) |
| components kept | `asv_sequences` 1, `genus` 1, `mock_counts` 1 (0.978) | the same (mock 0.975) |
| contaminant calls that change (of 43 scored sequences) | 8: all 5 reagent ASVs missed, 3 ASVs of the specimen listed | 10: 5 missed, 5 listed |
| specimen profiles below 0.80 | 6 of 8 | 5 of 8 |
| `l6_second_solution` (second solution of levels 3-6, unchanged) | 0.55: lists nothing (5 calls change), min agreement 0.407 | 0.55: lists nothing, 0.482 |
| planted counts: enrichment of the reagent ASVs with the sheet's roles (test asks for 2) / real roles | 0.69-1.03 / 4.87-6.55 | 0.97-1.26 / 5.47-6.67 |

Why: the two labelled blanks together are about four fifths stool community, so the reagent ASVs are not
enriched in "the blanks", while ASVs that are more abundant in the exchanged specimen than elsewhere are.
Over the 120 candidate exchanges (planted counts, sheet's roles): all five reagent ASVs missed in every
case (enrichment 0.58-1.52), 1-9 ASVs of the specimen listed, score 0.55 in every case. Development
harness, seeds 3-6: `l6_reference` 0.55 (3-6 specimen ASVs listed), `l6_second_solution` 0.55.

### Battery (gate runs; score, bound)

| entry | kind | seed 1 | seed 2 | bound |
|---|---|---|---|---|
| oracle | G2 | 1.000 | 1.000 | >= 0.99 |
| `planted_truth` | probe | 1.000 | 1.000 | >= 0.999 |
| `within_tolerance` | probe | 1.000 | 1.000 | >= 0.99 |
| `independent_solution` | probe | 1.000 | 1.000 | >= 0.99 |
| `stock_unoise3_uchime3` | probe | 0.977 | 1.000 | >= 0.9 |
| `drop_the_odd_blank` | probe | 1.000 | 1.000 | >= 0.99 |
| `counts_filed_under_the_other_sample_id` | probe | 1.000 | 1.000 | >= 0.99 |
| `l6_reference` | attack | 0.550 | 0.550 | <= 0.75 |
| `l6_second_solution` | attack | 0.550 | 0.550 | <= 0.75 |
| `sheet_roles_everything_else_right` | attack | 0.550 | 0.550 | <= 0.75 |
| `drop_the_odd_blank_and_its_counts` | attack | 0.700 | 0.700 | <= 0.75 |
| `one_column_moved_specimen_reads_lost` | attack | 0.700 | 0.700 | <= 0.75 |
| `reference_as_asvs` | attack | 0.000 | 0.000 | <= 0.1 |
| `half_of_the_read_pairs` | attack | 0.600 | 0.600 | <= 0.75 |
| `closed_reference_best_case` | attack | 0.559 | 0.483 | <= 0.75 |
| `naive_everything_in_blanks` | attack | 0.450 | 0.450 | <= 0.75 |
| `naive_no_primer_removal` | attack | 0.000 | 0.000 | <= 0.3 |
| `naive_fixed_length_trim` | attack | 0.000 | 0.000 | <= 0.5 |
| `naive_single_orientation` | attack | 0.250 | 0.600 | <= 0.75 |
| `naive_single_orientation_discard` | attack | 0.600 | 0.600 | <= 0.75 |
| `naive_no_chimera_removal` | attack | 0.650 | 0.650 | <= 0.75 |
| `minor_variants_merged` | attack | 0.327 | 0.327 | <= 0.75 |
| `naive_otu97` | attack | 0.327 | 0.327 | <= 0.75 |
| `naive_no_contaminant_handling` | attack | 0.550 | 0.550 | <= 0.75 |
| `naive_default_pipeline` | attack | 0.040 | 0.040 | <= 0.3 |
| `empty_output`, `inputs_copied`, `zero_byte_files`, `header_only` | attack | 0.000 | 0.000 | 0 / 0.1 |
| Harbor oracle / nop | G6 | 1.0 / 0.0 | not run | |

The level-3 entries score what they score at level 3 (table above): with the roles corrected the level-7
pipeline produces the level-3 result with two columns exchanged. The naive variants run the level-7
pipeline (controls checked) with their one flag, so each still ignores exactly its own decision point;
`naive_default_pipeline` also takes the roles from the sheet. Oracle and second solution: recall 1.0,
precision 1.0, min specimen agreement 0.977 / 0.979 (oracle) and 0.978 / 0.979 (second), mock 0.991 / 0.987,
MCC 1.0, genus 1.0, as at level 3. Development harness, seeds 3-6: oracle and second solution 1.000 on all
four, `drop_the_odd_blank` 1.000, `stock_unoise3_uchime3` 1.000 / 0.982 / 0.922 / 0.979.

**`drop_the_odd_blank` is a probe (valid route).** The submission recognises the odd blank, leaves both
libraries out of the contaminant identification and reports their counts like those of every other file.
Measured: 1.000 on seeds 1-6 (reference pipeline, `--roles drop`), contaminant list exact. Reasons: the
count table is defined per file and no deliverable asks for sample types, so nothing in the deliverables
depends on the role given to the two libraries; and one blank is enough, the reagent ASVs are enriched
4.0-6.5 / 5.0-6.6 fold in the remaining blank (planted counts, seeds 1 / 2; 3.8-7.2 over the 120 candidate
exchanges; the test asks for 2). The generator only plants an exchange for which this route returns exactly
the planted contaminants, so the route is valid by construction, not by luck. What is not valid is to
remove the two libraries from the table as well (`drop_the_odd_blank_and_its_counts`, 0.70: the profile of
one specimen is gone).

### Either filing of the count table (maintainer decision, 2026-11-12)

The first version of this level scored the count table as the files are and nothing else; a table that had
been "repaired" (the two columns moved to the `sample_id` of the extract they come from) scored 0.70. That
was reported as the one place where the level could fail somebody who did the QC, and the maintainer decided
to accept both readings. Implementation: see "Verifier" above; no tolerance and no weight was touched.

Planted truth with the right contaminants and the two columns filed in every possible way (level-7
verifier; identical on seeds 1 and 2; scored in the development harness, rows 1-3 are also battery entries
of the gate runs with the same scores):

| the two columns | as submitted | two columns exchanged | score | battery |
|---|---|---|---|---|
| as the files are (planted truth) | 1.00 | 0.70 | 1.00 | probe `planted_truth` |
| by the sample each library really is (both columns moved) | 0.70 | 1.00 | 1.00 | probe `counts_filed_under_the_other_sample_id` |
| only the blank's read pairs moved: blank library in both columns, specimen's read pairs nowhere | 0.70 | 0.70 | 0.70 | attack `one_column_moved_specimen_reads_lost` |
| blank's read pairs moved, the other column zero | 0.70 | 0.70 | 0.70 | - |
| both columns zero | 0.70 | 0.70 | 0.70 | attack `drop_the_odd_blank_and_its_counts` (pipeline run) |
| only the specimen's read pairs moved: specimen library in both columns, blank's read pairs nowhere | 1.00 | 1.00 | **1.00** | - |
| specimen's read pairs moved, the other column zero | 0.70 | 1.00 | **1.00** | - |

So the rule the verifier implements is: the specimen's library must be under one of the two `sample_id`s
with the right counts. The last two rows are inconsistent tables that pass; see the weaknesses. The
sheet-role submissions gain nothing from the second reading (`l6_reference`, `l6_second_solution`,
`sheet_roles_everything_else_right`: 0.55 as submitted, 0.55 with the columns exchanged, both seeds), and no
other battery entry changes its score.

### Byte identity of the lower levels

```
C=<scratchpad>/amplicon_l7_check
git archive HEAD pipeline/alegen | tar -x -C $C/head ; ln -s /home/user/ale-lifesci-gym/.cache $C/head/.cache
(cd $C/head/pipeline && python -m alegen make amplicon_16s --seed 1 --level L --no-validate --out $C/head_out --results $C/head_results)
(cd pipeline         && python -m alegen make amplicon_16s --seed 1 --level L --no-validate --out $C/final_out --results $C/final_results)
diff -r $C/head_out/amplicon-16s-s0001-lL $C/final_out/amplicon-16s-s0001-lL
```
for L = 3, 6: no difference, 42 files each (final code, i.e. with the level-7 verifier in place, against
HEAD `ba7d1b3`; `assets/score.py` is byte-identical to HEAD). Before the verifier change the same check was
also done for seed 1 at levels 1 and 2 and for seed 2 at levels 3 and 6 (against HEAD `4dec6ec`): no
difference; the level-3 code path has not been touched since. The level-3 and level-6 task directories of
seed 1 in `tasks/` are identical to what the final code writes, README.md and task.toml included. The
level-7 tasks are deterministic (regenerated into scratch, `diff -r` against `tasks/`: identical, 43 files).

### Weaknesses and open points (level 7)

- **Not every table with one column moved is rejected.** Asked for: an inconsistent table (only one of the
  two columns moved) stays below the pass line. That holds when the move loses the specimen's read pairs
  (blank library in both columns: 0.70, in the battery). It does not hold for the opposite move: with the
  specimen's library in both columns, or under the specimen's `sample_id` with the other column zero, the
  score is 1.00 on both seeds. The column that would have to hold the blank library is a blank column, and
  blank columns are not scored at any level ("the counts in the negative controls" are deliberately left
  out since level 3), so the level-3 logic has nothing to object to under either reading. Rejecting these
  tables needs a new check on the blank column (for instance that the two columns are not copies of each
  other, or a count agreement for blank libraries), which is a change of the scoring beyond "both readings,
  better total". Not done and not tuned; measured and reported. In such a table every scored quantity is
  right (all specimen profiles are available, the contaminant list is right); what is wrong is the column
  of a blank.
- **Easy to see once looked at.** The blanks have a quarter of the read pairs of the other libraries, so
  `ls -l` on the FASTQ directory already shows a "blank" of specimen size and a "specimen" of blank size.
  The level tests whether the controls are looked at before they are used, not whether subtle evidence is
  read correctly.
- **A half-noticed route passes.** Leaving out only the odd blank and keeping the blank-like library as a
  specimen scores 1.000 with the reference pipeline on seeds 1 and 2 (`--roles drop_blank`, development
  harness; exact on the planted counts of all 120 candidate exchanges): the enrichment test needs one clean
  blank and is not disturbed by a blank among the specimens. The absolute-abundance test of the second
  solution fails in that situation (0 of 120 exact on the planted counts). Not in the battery.
- **A method that never uses the sample types can pass without noticing anything** (for example a robust
  regression of relative abundance on qPCR copies). Not tested. The deliverables of such a submission are
  correct, so the scorer has nothing to object to.
- **The checks of the two solutions are written for this decision point.** One blank against one
  high-biomass stool specimen, exactly one exchange. The oracle looks for exactly one (labelled blank,
  labelled specimen) pair; the second solution for exactly two libraries that are at least a decade off in
  opposite directions. Two exchanges, an exchange of two specimens (oracle), or an exchange with a
  low-biomass specimen or the mock are not planted and not handled; in those cases both fall back to the
  sheet and say so on stdout. Thresholds (Bray-Curtis 0.25 / 0.50, half the median richness, one decade)
  have margins measured on seeds 1-6 only (0.08 against 0.39; 1.6 decades against 0.26).
- **Second solution, constant-amount ASVs.** On seed 2 six ASVs pass the spread filter (five reagent ASVs
  and one other); the median over them is not affected there, but that is measured, not proved.
- **Generator-side replica of the scorer.** Condition (2) uses a re-implementation of the level-3 component
  formulas (`l7._score`) and of the two contaminant tests. If `score.py` or the tests change, the replica
  must follow; the gate is the validator's run of the real scorer (`sheet_roles_everything_else_right`,
  whose score equals the replica's on both seeds).
- **Seeds.** The level-3 generator stops with "pair split did not converge" on seeds 7, 12, 14 and 15 (of
  1-16), so level 7 does not exist for them either. Gate runs: seeds 1 and 2 only; seeds 3-6 went through
  the development harness (no image build, no Harbor).
- **`stock_unoise3_uchime3`** keeps its level-3 fragility (floor 0.9): 0.922 on seed 5.
- **Judge text.** The level-7 paragraph of the task card quotes the measured numbers of seeds 1 and 2 on
  every instance; the instance's own numbers are in the trap entry. The sentence the pipeline puts in front
  ("same inputs and planted truth as level 7") comes from `descaffold.py` and is not accurate for this
  family; not edited (core).
- `task.toml` lists `M_control_swap` under `traps` and a capability about checking controls, like the other
  families do for their decision points; neither file is in the agent's container.
