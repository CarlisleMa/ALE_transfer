# Family spec: `rnaseq_confound` - bulk RNA-seq from FASTQ to differential expression

ALE subdomain: 3.2 Genomics & Sequence Analysis. (ALE's TCGA task starts from a finished matrix and
prescribes the test; this one starts from reads and the difficulty is the experimental design.)

## Scenario
A pharmacology group treated donor-derived cells with a compound (treated vs vehicle, 6 vs 6) and
hands over FASTQ files, a sample sheet and wet-lab notes. Deliver gene-level counts, sample QC and
the differential-expression table.

## Truth model
- Real transcript sequences (Ensembl human cDNA; a subset of ~800 genes including sex-specific
  genes XIST, RPS4Y1, DDX3Y, KDM5D, UTY, EIF1AY and some multi-isoform genes). Cache the download.
- Planted effects per gene: treatment log2FC (about 10% of genes, |LFC| 1-3, both directions);
  library-prep batch log2FC on about 30% of genes (|LFC| 0.5-2); donor sex. Negative-binomial counts
  with gene-wise dispersion; fragments sampled from the transcripts into single-end reads with
  sequencing errors. True per-gene fragment counts are kept as hidden truth.
- Design: two library-prep batches, unbalanced across conditions (e.g. 4/2 vs 2/4) so the batch
  effect is estimable but strongly aliased with treatment when ignored.
- Planted sample mix-up: the labels of two samples (one treated, one vehicle, different donor sex)
  are exchanged in the sample sheet. Evidence: expression of sex-specific genes contradicts the
  recorded donor sex for exactly those two samples.

## Inputs
`fastq/*.fastq.gz` (12 samples), `samples.tsv` (sample_id, condition, donor_id, donor_sex,
library_prep_batch, fastq_file), `reference/transcripts.fa.gz`, `reference/tx2gene.tsv`,
`lab_notes.md` (facts: prep dates, a relabelling incident; no instructions).
Tools: salmon, python (pydeseq2, statsmodels), and R with DESeq2/edgeR/limma/tximport if the image
stays within budget.

## Outputs
- `gene_counts.tsv`: gene_id x sample estimated counts.
- `sample_qc.tsv`: `sample_id, inferred_sex, condition_used` where condition_used is the condition
  the analyst used in the model, or `excluded`.
- `de_results.tsv`: `gene_id, log2fc, padj, significant` (treated vs vehicle; significant = 1 at
  FDR 5%, no fold-change cutoff), one row per gene.

## Traps
batch_confound; sample_swap (identity check via sex genes). Accept both resolutions of the mix-up
(labels corrected, or the two samples excluded) because both give valid inference.

## Scoring sketch
gates: files parse, every gene and sample present. Components: quantification agreement with true
fragment counts (0.15); sample_qc (sex inferred correctly; the two swapped samples handled) (0.15);
false-discovery proportion among `significant` genes against planted truth (0.30); sensitivity on
planted DE genes (0.25); log2FC error on planted DE genes (0.15). Calibrate bands from two correct
analyses (e.g. DESeq2 with batch term vs limma-voom or pyDESeq2 with batch term) and from the naive
ones (no batch term; labels as given; both).
