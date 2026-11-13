You are the bioinformatician of a cancer-genomics core facility. Tumour and normal DNA of three
patients was sequenced on a small hybrid-capture panel. The clinical team needs, for every patient,
the tumour/normal pair and the somatic variants of the tumour.

Task directory:
- `/task`

Inputs (read-only) under `/task/input`:
- `fastq/<sample_id>_R1.fastq.gz`, `fastq/<sample_id>_R2.fastq.gz` - paired-end reads (2 x 100 nt),
  one pair of files per library (six libraries)
- `manifest.tsv` - sample sheet exported from the LIMS: `patient`, `sample_type` (`tumour` or
  `normal`), `sample_id`, `fastq_1`, `fastq_2`
- `reference/ref.fa` (+ `.fai`) - the panel reference, one contig (`chr22_panel`); no aligner index is
  provided
- `targets.bed` - the capture targets on the panel reference (0-based, half-open)
- `core_notes.md` - hand-over notes of the sequencing core

Available software (on `PATH`): `bwa`, `minimap2`, `samtools`, `bcftools` (with `bgzip`, `tabix`),
Python 3.12 with `pysam`, numpy, scipy and pandas. There is no internet access.

Your task:
1. Align the reads of every library to the panel reference.
2. Establish, for every patient, which normal library was prepared from the same individual as the
   patient's tumour library.
3. Call the somatic single-nucleotide variants and short insertions/deletions of every tumour
   inside the capture targets.

Required outputs under `/task/output`:
- `pairing.tsv` - tab-separated with the columns `patient`, `tumour_sample`, `normal_sample`; one
  row per patient of the manifest. `tumour_sample` is the `sample_id` of the patient's tumour
  library as listed in the manifest; `normal_sample` is the `sample_id` of the normal library that
  comes from the same individual as that tumour library.
- `somatic/<patient>.vcf.gz` - one gzip- or bgzip-compressed VCF per patient holding the somatic
  variants of that patient's tumour library.
  - A somatic variant is a single-nucleotide variant or short insertion/deletion that is present in
    tumour cells of the patient and absent from the patient's germline. Germline variants, and
    differences from the reference that arise during library preparation, sequencing or read
    alignment, are not somatic variants.
  - Reported variants have `FILTER` = `PASS`; records with any other FILTER value are not part of
    the deliverable.
  - `CHROM`, `POS`, `REF`, `ALT` follow the VCF specification relative to `reference/ref.fa`; one
    alternate allele per record; insertions and deletions left-aligned and parsimonious (the
    representation `bcftools norm` produces).
  - Every record carries `AF=<number>` in the `INFO` column: the fraction (0-1) of the tumour
    library's DNA fragments covering the site that carry the alternate allele.
  - The reporting range of the assay starts at an allele fraction of 5%; variants below it may be
    left out.
- `summary.json` - `{"patients": {"<patient>": {"tumour_sample": "<sample_id>", "normal_sample":
  "<sample_id>", "n_snv": <int>, "n_indel": <int>}, ...}}` where `n_snv` and `n_indel` are the
  numbers of `PASS` single-nucleotide and insertion/deletion records in that patient's VCF and the
  two sample ids repeat `pairing.tsv`.

Constraints:
- Do not modify anything under `/task/input`; write deliverables only under `/task/output`.
- Extra columns, INFO fields and sample columns are allowed.
- Do not ask for confirmation. Execute directly.
