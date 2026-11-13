# Tumour/normal pairs: FASTQ to filtered somatic variants

Align six panel libraries of three cancer patients, establish which normal belongs to which tumour, and deliver filtered somatic SNV/indel calls with allele fractions per patient.

- family: `somatic_tn`  seed: `1`  level: `3`
- ALE subdomain: 3.2 Genomics & Sequence Analysis
- software: bwa, minimap2, samtools, bcftools, pysam, Python
- a-priori difficulty: 0.702 (last-exam)  |  pass threshold: 0.9
- capabilities trained: sample identity QC before analysis, read-level artefact diagnosis (orientation, mapping quality), long tool chain from raw reads to filtered calls, tumour/normal comparison without a turnkey caller, exact output contracts (VCF normalisation)

## Undisclosed decision points (not shown to the agent)

- **sample_swap** - The normal libraries of PT-0452 and PT-0740 are exchanged in manifest.tsv. With the sheet pairs, every germline variant private to the patient is called somatic in two of three patients.
  - discoverable by: pairwise germline genotype concordance between the six libraries: each tumour agrees with exactly one normal at ~100% of variable sites and with the others at well under 60%; for two patients the matching normal is not the one in the sheet. Also visible downstream: dozens of 'somatic' calls at allele fraction ~0.5 and ~1.0 in two patients
  - trains against ALE failure mode: Understanding / unverified assumption about inputs
- **duplicated_segment** - A 3988 bp segment occurs twice in the panel reference (two targets). Read pairs from its interior align equally well to both copies; a tumour variant of one copy appears at both, with half the allele fraction, supported by mapping-quality-0 reads only.
  - discoverable by: mapping-quality distribution per target (two targets consist almost entirely of MAPQ 0 pairs); the candidate calls there come in pairs at identical offsets in the two targets; aligning the reference against itself shows the duplication
  - trains against ALE failure mode: Approach / missing QC of alignment evidence
- **orientation_artifacts** - The tumour libraries of PT-0760 carry 14-18 C>T / G>A (deamination) and G>T / C>A (oxidation) artefact sites each at 8-30% allele fraction. Every supporting fragment has the same read-pair orientation (F1R2 for one substitution, F2R1 for its complement); forward/reverse strand counts are balanced, so a plain strand-bias filter does not see them.
  - discoverable by: the unfiltered calls of the affected tumours contain an excess of two substitution types that the other tumours do not show; tabulating alt-supporting reads by read number and strand shows 100% F1R2 or 100% F2R1 at those sites while true variants are balanced; core_notes.md lists the specimen type of each library without comment
  - trains against ALE failure mode: Understanding / domain knowledge gap

## Judge

**Where the truth comes from.** Before any read exists the generator fixes (a) the germline
genotype of each patient at 110 population SNP sites, (b) the somatic SNVs and indels of each
tumour with their allele fractions, (c) artefact sites and the read-pair orientation that shows
them, (d) which library is written to which sample id, and then (e) which sequenced fragment
carries which allele. `reference/truth.json` holds the planted variants and, per variant, the
realised fraction of tumour fragments that carry it. No aligner or caller output is used.

**Scored set and gray zone.** A planted somatic variant is *scored* if the simulated tumour library
holds at least 8 alt-carrying fragments at the site, a realised allele fraction of at least
0.07, and at least max(3, 0.25 x alt fragments) alt fragments in the rarer
read-pair orientation. Planted variants below that floor, and planted variants in the interior of
the duplicated segment (every supporting read pair maps equally well to both copies), are *gray*:
calling them is neither rewarded nor penalised. Every other PASS call inside the targets is a false
positive; that includes the patient's germline variants, orientation artefacts and the mirror
position of a duplicated-segment variant in the other copy. Calls outside the targets are ignored.

**Hard gates (score 0).** Missing/unparseable file; a patient missing or duplicated in
`pairing.tsv`; a VCF record on another contig, outside the reference, or whose REF does not match
`ref.fa`; `summary.json` without a `patients` object.

**Components.** `recall` and `precision` are computed per patient and the **minimum** over the
three patients is used, so a patient whose pair was left wrong cannot be averaged away.

| weight | component | metric | full credit | zero |
|---|---|---|---|---|
| 0.20 | pairing | every patient's (tumour, normal) pair equals the planted one | all | any wrong |
| 0.30 | somatic_precision | min over patients of TP / (TP + FP) | >= 0.95 | <= 0.85 |
| 0.30 | somatic_recall | min over patients of TP / scored variants; multiplied by ramp(precision, 0.50, 0.85) | >= 0.95 | <= 0.85 |
| 0.10 | allele_fraction | worst patient's median abs(AF - realised fragment fraction) over true positives | <= 0.03 | >= 0.08 |
| 0.10 | summary | fraction of patients whose summary entry matches their own VCF counts and pairing row | all | none |

pairing, somatic_precision, allele_fraction and summary are multiplied by ramp(min recall, 0.2,
0.5): a deliverable without real calls earns nothing.

**Why these tolerances.** See `docs/families/somatic_tn.md` for the measured values of the oracle
(bwa + pysam fragment counting), the independent solution (minimap2 + bcftools mpileup on
orientation-split BAMs) and every trap-ignoring analysis on two seeds.

**What is deliberately not scored.** QUAL, genotype columns, variants below the scoring floor,
anything outside the targets, the numeric concordance values behind the pairing.

## Verifier attack suite (must score low)

- `sheet_pairing_no_calls` (max 0.1): prior-only answer: pairing copied from the sample sheet, VCFs without records
- `false_calls_beyond_tolerance` (max 0.75): planted truth plus 30% false calls per patient
- `missed_calls_beyond_tolerance` (max 0.75): planted truth with a quarter of the variants missing
- `naive_tumour_only` (max 0.5): correct alignment and filters, but no matched normal is used: germline variants are reported as somatic
- `naive_sheet_pairing` (max 0.75): all filters applied, but the pairs of the sample sheet are used as given
- `naive_no_mapq_filter` (max 0.75): verified pairs and orientation filter, but reads of any mapping quality are counted
- `naive_default_pipeline` (max 0.4): sample sheet pairs, no mapping-quality filter, no orientation filter
- `naive_no_orientation_filter` (max 0.75): verified pairs and mapping-quality filter, but no check of the read-pair orientation of the supporting reads
- `naive_strand_bias_only` (max 0.75): verified pairs and mapping-quality filter; supporting reads are required on both alignment strands (classic strand-bias filter) instead of in both read-pair orientations
- `shortcut_spectrum_filter` (max 0.75): verified pairs and mapping-quality filter; instead of looking at the supporting reads, every C>T, G>A, G>T and C>A call below 30% allele fraction is dropped in the tumours the notes list as FFPE (removes the artefacts and the true variants of those classes alike)

## Must-pass probes (guard against false negatives)

- `planted_truth` (min 0.999): submission written directly from the planted truth
- `within_tolerance` (min 0.99): per patient one scored variant missed, one false call, AF noise +-0.02
- `independent_solution` (min 0.99): second tool chain: minimap2 (no bwa), bcftools mpileup allele depths on orientation-split BAMs (no pysam), genotype discordance from bcftools call for the pairing
