# SOP FG-017: pooled CRISPR drug-modifier screens, analysis and release (rev. 4)

## 1. Scope
Secondary (validation) knockout screens with a T0 reference, a vehicle arm and a drug arm. This
SOP defines the quantities of a released analysis. It does not prescribe software.

## 2. Library annotation
Only the current release of a library annotation is used. For VAL-2 this is release 2
(`library.tsv`). Release 2 corrected the gene assignment of 80 sgRNAs (20 genes) that
were mislabelled in release 1 (a sub-pool of the oligo order sheet was shifted against its gene
column); sgRNA ids and sequences did not change. Release 1 is withdrawn.

## 3. Counting
- A read is counted for an sgRNA if it carries that sgRNA's spacer in the amplicon context
  described in the wet-lab notes. Every library of the sample sheet is counted and reported.

## 4. Scales
Log2 fold changes are reported on the scale on which the non-targeting sgRNAs do not change
(median non-targeting sgRNA = 0). The fitness effect is vehicle day 14 against T0. The drug effect
is drug day 14 against vehicle day 14; negative values mean depletion under drug.

## 5. Erratum
If results have left the facility, an erratum lists every gene of the mailed list with one verdict:
`confirmed` if the gene is a modifier (either class) in the released analysis, `retracted` if not.
