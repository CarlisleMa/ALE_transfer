# Family spec: `somatic_tn` - tumour/normal pairs from FASTQ to somatic calls

ALE subdomain: 3.2 Genomics & Sequence Analysis. (ALE's WGS task is single-sample germline calling
with a prescribed tool chain; this is paired somatic calling where sample identity and artefacts
decide the result.)

## Scenario
A cancer genomics core receives tumour and normal FASTQ for three patients on a small capture
region. Deliver verified tumour/normal pairing and filtered somatic SNV/indel calls per patient.

## Truth model
- Real reference slice (UCSC hg38, ~120 kb; cache), with a planted intra-region duplication so that
  some reads are ambiguously placed.
- Per patient: private germline SNPs (a fingerprint), somatic SNVs/indels at planted allele
  fractions (clonal and subclonal), tumour purity.
- Read simulator (paired-end, written by the family): sequencing errors, duplicates, plus planted
  artefact classes with known sites: orientation-specific damage (one read orientation only),
  paralog-driven mismatches in the duplicated segment.
- Manifest mix-up: the normal files of two patients are exchanged in the manifest.

## Inputs
`fastq/*.fastq.gz`, `manifest.tsv` (patient, sample type, files), `reference/ref.fa` (+ indexes),
`targets.bed`, `core_notes.md`. Tools: bwa, samtools, bcftools, GATK4 (Mutect2), python (pysam).

## Outputs
`pairing.tsv` (`patient, tumour_sample, normal_sample` as verified), `somatic/<patient>.vcf.gz`
(PASS records only), `summary.json` (per patient counts).

## Traps
sample_swap (germline concordance between the manifest pairs is near random for two patients);
orientation_artifacts; duplicated_segment (low mapping quality evidence).

## Scoring sketch
gates: files parse. Components: pairing (0.20); per-patient somatic F1 against planted variants,
minimum over patients (0.60); allele-fraction error on true positives (0.10); summary consistency
(0.10). Naive variants: manifest pairing used as given; no artefact filtering; tumour-only calling.
