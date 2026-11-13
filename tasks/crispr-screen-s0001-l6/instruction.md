You are the computational biologist of a functional-genomics core facility. You are taking over a
pooled CRISPR-Cas9 knockout drug-modifier screen (secondary validation library VAL-2) from a
colleague who has left. He had started the analysis and mailed a draft hit list to the PI. You
deliver the final results and an erratum for the mailed list.

Task directory:
- `/task`

Inputs (read-only) under `/task/input`:
- `fastq/<sample_id>.fastq.gz` - one single-end FASTQ per library (8 libraries)
- `samples.tsv` - sample sheet as it was submitted with the libraries (`sample_id`, `arm`,
  `replicate`, `fastq_file`, `sequencing_run`); `arm` is `T0`, `vehicle` or `drug`
- `library.tsv` - current release (r2) of the sgRNA library: `sgRNA`, `sequence`, `gene`.
  Non-targeting control sgRNAs have `gene = NON_TARGETING`. Gene identities were blinded by the client.
- `reference_genes.tsv` - library genes whose behaviour under this compound is established
  (`gene`, `expected_class`)
- `screen_design.md` - wet-lab notes on the library, the treatment arms, PCR and sequencing
- `screen_sop.md` - the facility SOP for drug-modifier screens: definitions of the reported
  quantities and of the erratum
- `handover/` - the colleague's project folder: `HANDOVER.md`, `count_guides.py`,
  `analyse_screen.py`, `counts_draft.tsv`, `gene_draft.tsv`, `hits_mailed.tsv`

Available software (on `PATH`): `cutadapt`, `mageck` (MAGeCK), Python 3.12 with numpy, scipy and
pandas. There is no internet access.

Your task:
1. Count how many reads carry each library sgRNA in every library.
2. For every targeted gene, estimate the fitness effect of the knockout and its effect on the
   response to the drug, and classify the gene as sensitiser, resistance gene or neither.
3. Write the erratum for the hit list that was mailed.

Required outputs under `/task/output`:
- `counts.tsv` - tab-separated, header `sgRNA`, `gene`, then one column per `sample_id` in the
  order of `samples.tsv`. One row per sgRNA of `library.tsv` (non-targeting sgRNAs included), `gene`
  as in `library.tsv`. Values are raw integer read counts of the FASTQ file listed for that sample.
- `gene_summary.tsv` - tab-separated with at least the columns `gene`, `lfc_fitness`, `lfc_drug`,
  `drug_class`; one row per targeted gene (no row for `NON_TARGETING`). Extra columns are allowed.
  - `lfc_fitness`: log2 fold change of the gene, vehicle-treated libraries of day 14 versus T0,
    defined as the median over the gene's sgRNAs of the sgRNA-level log2 ratio of depth-normalised
    abundance (mean of the replicates of one group over mean of the replicates of the other). The
    scale must be such that a knockout without any fitness effect has `lfc_fitness` = 0.
  - `lfc_drug`: the same quantity for drug-treated versus vehicle-treated libraries of day 14;
    negative means less abundant under drug. A knockout that does not change the drug response
    has `lfc_drug` = 0.
  - `drug_class`: `sensitiser` if the gene's sgRNAs are significantly less abundant under drug than
    under vehicle, `resistance` if they are significantly more abundant, otherwise `none`; false
    discovery rate 5%.
- `erratum.tsv` - tab-separated, columns `gene`, `verdict`; one row per gene of
  `handover/hits_mailed.tsv`. `verdict` is `confirmed` if your final `drug_class` of the gene is
  `sensitiser` or `resistance`, and `retracted` if it is `none`.
- `qc.json` - `{"samples": {"<sample_id>": {"total_reads": <int>, "mapped_reads": <int>}, ...}}`
  where `total_reads` is the number of reads in the FASTQ, `mapped_reads` the number of reads
  assigned to a library sgRNA.

Constraints:
- Do not modify anything under `/task/input`; write deliverables only under `/task/output`.
- Every library in the sample sheet belongs to the experiment and must be reported and used.
- Do not ask for confirmation. Execute directly.
