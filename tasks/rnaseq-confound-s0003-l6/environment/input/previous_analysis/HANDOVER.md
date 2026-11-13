# Hand-over note: compound vs vehicle RNA-seq

Written on my last day. Sorry for the state of things.

## What is in this folder
- `scripts/01_quantify.sh` - salmon index and quantification of the 12 FASTQ files
- `scripts/02_merge_counts.py` - gene-level matrix from the salmon output
- `scripts/03_de.py` - DESeq2 model (pydeseq2), compound vs vehicle
- `gene_counts_merged.tsv` - gene-level matrix (output of 01 + 02)
- `de_draft.tsv` - draft differential-expression table (output of 03)
- `mail_top_responders.tsv` - the 40 genes I mailed to the wet lab on 2 July as "top responder genes"

## Notes
- The salmon output directories were on the scratch disk, which has been purged since. Only the
  merged matrix is left.
- The mailed list came from my first differential-expression run (smallest adjusted p-values first,
  sex-chromosome control genes taken out). The table of that run was overwritten; `de_draft.tsv` is
  from a later run.
- Samples were used as listed in the sample sheet.
- The wet lab has started to order qPCR primers for the mailed genes, so they need to know soon
  which of them hold up.
- All scripts take their paths as arguments, see the usage line at the top of each.
