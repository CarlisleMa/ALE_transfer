# Level-4 spec: `rnaseq_confound` - taking over an RNA-seq analysis (H + Q)

Pilot at level 3 (seed 1): Claude Haiku scored 1.0 in 3 min 4 s with 16 tool calls. It quantified,
checked sex markers, found the mix-up, added the batch term.

## Scenario
Same experiment. A predecessor quantified the libraries, ran a differential-expression script,
mailed "top responder genes" to the wet lab and left. The analyst takes the project over.

## Inherited project: `/task/input/previous_analysis/`
Quantification script, merged gene count matrix, DE script, draft DE table, the mailed gene list,
a hand-over note. Produce them by running **mutated copies of the reference pipeline**. The SOP
holds the audit checklist: 8-10 items, 3-5 defective per seed, the rest genuinely fine (include
fine items that look unusual, e.g. an uncommon but correct parameter). Candidate operators:
count-matrix columns in lexicographic order (S1, S10, S11, S12, S2 ...) joined to a numerically
ordered sample sheet; transcript-to-gene join on versioned vs unversioned ids that silently drops
genes; library type forced to stranded on unstranded reads; one sample quantified from another
sample's FASTQ; design formula without the batch term; reversed contrast; unadjusted p-values
written as `padj`; a draft table older than the count matrix next to it. At least two defects must
be detectable only by recomputing (comparing the inherited matrix with a fresh quantification),
not by reading the code.

## New element in the analysis itself (Q) - choose with measurements
The level-3 analysis must not be sufficient for the final table. Preferred: (a) a strong donor
effect in a paired design in which donor and library batch are partly aliased, so that the model
has to be specified with care (an unpaired model loses sensitivity below the band; a naively
over-specified one is rank deficient); and/or (b) one library that is a cross-contaminated mixture
of two donors (identity markers intermediate), which the SOP's stated identity criterion excludes
and which cannot be repaired by relabelling. State every acceptance criterion in the SOP. Keep
whichever you can show to separate correct from careless analyses on two seeds with margin; say
which you dropped and why.

## Deliverables
The level-3 files plus `audit.tsv` (`item`, `verdict` in `ok`/`defective`) and `erratum.tsv` (one
row per gene of the mailed list: `gene_id`, `verdict` in `confirmed`/`retracted`).

## Scoring
Level-3 components rescaled to 0.65; audit MCC 0.20; erratum MCC 0.15. Pass 0.9.

## Battery
`l3_reference`; `trust_inherited` (finish from the inherited matrix); one variant per active defect
left unfixed; the careless variant of the Q element; audit "everything defective" / "everything
ok"; beyond tolerance. Probes as in levels 1-3.
