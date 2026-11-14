You are the bioinformatician of a pharmacology group and take over an RNA-seq project from a
colleague who has left. The wet lab treated donor-derived primary cells with a candidate compound or
with vehicle and prepared poly-A RNA-seq libraries. Your predecessor quantified the libraries, ran a
differential-expression analysis and mailed a list of "top responder genes" to the wet lab. You
release the final result: gene-level counts, a per-sample QC table, the genes that respond to the
compound and an erratum for the mailed list.

Task directory:
- `/task`

Inputs (read-only) under `/task/input`:
- `fastq/<sample_id>.fastq.gz` - one single-end, unstranded FASTQ per sample (12 samples)
- `samples.tsv` - sample sheet as filled in by the wet lab (`sample_id`, `condition` = `treated` or
  `vehicle`, `donor_id`, `donor_sex` = `F` or `M`, `library_prep_batch`, `fastq_file`)
- `reference/transcripts.fa.gz` - transcript sequences of the gene panel that was sequenced
- `reference/tx2gene.tsv` - `transcript_id`, `gene_id`, `gene_name` for every transcript
- `lab_notes.md` - wet-lab notes for this experiment
- `previous_analysis/` - what your predecessor left: `HANDOVER.md`, `scripts/01_quantify.sh`,
  `scripts/02_merge_counts.py`, `scripts/03_de.py`, `gene_counts_merged.tsv`, `de_draft.tsv`,
  `mail_top_responders.tsv`

Available software (on `PATH`): `salmon`, Python 3.12 with numpy, scipy, pandas, statsmodels,
scikit-learn and pydeseq2. There is no internet access and no R installation.

Your task:
1. Quantify every sample against the provided transcripts and summarise to genes.
2. Decide how every sample enters the statistical analysis.
3. Test every gene for differential expression between compound-treated and vehicle-treated cells.
4. Tell the wet lab which of the mailed genes hold up.

Required outputs under `/task/output`:
- `gene_counts.tsv` - tab-separated, header `gene_id` followed by one column per `sample_id`
  (all 12 samples, named as in `samples.tsv`). One row per `gene_id` of `tx2gene.tsv`.
  Values are the estimated number of reads of that FASTQ file that originate from the gene
  (non-negative, not length- or depth-normalised; fractional values are allowed).
- `sample_qc.tsv` - tab-separated with at least the columns `sample_id`, `inferred_sex`,
  `condition_used`; one row per sample of `samples.tsv`.
  - `inferred_sex`: `F` or `M`, the sex of the donor of the library contained in that sample's
    FASTQ file, as determined from the sequencing data; `mixed` for a library that contains
    material of more than one donor.
  - `condition_used`: `treated` or `vehicle` - the condition under which the reads of that FASTQ
    file entered your differential-expression model - or `excluded` if you left the sample out.
    A library that contains material of more than one donor cannot be assigned to a record of the
    sample sheet and is excluded; every other library enters the model under the condition its
    cells actually received.
- `de_results.tsv` - tab-separated with at least the columns `gene_id`, `log2fc`, `padj`,
  `significant`; one row per `gene_id` of `tx2gene.tsv`.
  - `log2fc`: estimated log2 fold change of the gene's expression caused by the compound (treated
    over vehicle), on the scale on which a gene that does not respond has `log2fc` = 0. Must be a
    finite number for every gene (use 0 for a gene that cannot be estimated).
  - `padj`: multiple-testing-adjusted p-value of the test "the compound does not change the
    expression of this gene"; `NA` for a gene that was not tested.
  - `significant`: `1` if the gene is differentially expressed at a false-discovery rate of 5%
    (`padj` < 0.05, no fold-change cutoff), otherwise `0`.
- `erratum.tsv` - tab-separated with at least the columns `gene_id`, `verdict`; one row per gene of
  `previous_analysis/mail_top_responders.tsv`; `verdict` = `confirmed` if the gene is
  differentially expressed in your final analysis (`significant` = 1 in your `de_results.tsv`),
  otherwise `retracted`.

Constraints:
- Do not modify anything under `/task/input`; write deliverables only under `/task/output`.
- Extra columns are allowed in every output file.
- Do not ask for confirmation. Execute directly.
