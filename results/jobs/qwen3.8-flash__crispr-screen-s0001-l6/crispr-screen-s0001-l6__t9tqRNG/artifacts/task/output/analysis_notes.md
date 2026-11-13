# VAL-2 released analysis (SOP FG-017 rev. 4) - method notes and erratum background

## 1. Counting (`counts.tsv`, `qc.json`)
Reads were counted for an sgRNA when the read carried that sgRNA's 20-nt spacer directly in front of
the scaffold (`GTTTTAGAGCTA`), i.e. in the amplicon context given in the wet-lab notes
(`...ACACCG [spacer] GTTTTAGAGCTAGAAATAG...`). Both orientations were accepted: run A (NextSeq, 1x75)
reads are sense, run B (NovaSeq LIB15) reads are the reverse complement of the amplicon, which is why
the hand-over script reported **0** reads for LIB15 (its `count_guides.py` only looks for the sense
scaffold). With the orientation fix LIB15 yields 116805 mapped reads (90.3%), in line with the other
libraries. All eight libraries of the sample sheet are counted and reported. Gene labels are from the
current release (`library.tsv`, r2); the hand-over `counts_draft.tsv` still used the withdrawn release 1,
in which 80 sgRNAs (20 genes) were mislabelled.

## 2. Arm assignment - LIB10 and LIB15 are swapped in the submitted sample sheet
The submitted sheet labels LIB10 as vehicle and LIB15 as drug. The libraries themselves say the
opposite, and three independent lines of evidence agree:

1. **Reference genes.** Only the partition drug = {LIB10, LIB13, LIB14} vs vehicle = {LIB11, LIB12, LIB15}
   reproduces the compound classes: the five known resistance genes go up and the five known
   sensitisers go down by 1-2.7 log2 (signature gap 3.10 log2). Scanning all ten possible 3-vs-3
   splits of the six day-14 libraries, every other split gives a gap of <= 0.99.
2. **Template dropouts.** The wet-lab notes state ~150 cells/sgRNA for the vehicle arm and only
   ~25 cells/sgRNA for the drug arm, so the drug libraries must show Poisson zeros. LIB10, LIB13,
   LIB14 have 13, 11, 12 sgRNAs with zero reads and count CV ~1.12-1.15, while LIB11, LIB12, LIB15
   have 1, 0, 0 zeros and CV ~0.97. LIB10 therefore belongs to the drug arm, LIB15 to the vehicle arm.
3. **Replicate correlation.** In log space LIB15 correlates 0.954/0.955 with LIB11/LIB12 but only
   0.83 with LIB10, so LIB15 is a genuine vehicle replicate (not a re-sequencing of a drug library).

All reported quantities use the corrected arms (T0 = LIB8, LIB9; vehicle = LIB11, LIB12, LIB15;
drug = LIB10, LIB13, LIB14); no library is dropped. This also fixes the hand-over analysis, which
silently discarded LIB15 ("a library without reads is left out") and analysed 2 vehicle + 2 drug
libraries.

## 3. Gene table (`gene_summary.tsv`)
Counts were depth-normalised by the mapped-read total of each library and log2 ratios were centred on
the non-targeting guides so that the median non-targeting sgRNA is exactly 0 (SOP section 4).
`lfc_fitness` = median over the gene's sgRNAs of log2(mean vehicle / mean T0);
`lfc_drug` = median over the gene's sgRNAs of log2(mean drug / mean vehicle).
Using median-of-ratios size factors instead changes nothing (max 0.017 log2). The scale behaves as
the design predicts: 190 genes (the carried-forward fitness candidates) are depleted below -1 and the
40 putatively neutral genes sit within +-0.5.

Significance was taken from a count-aware negative-binomial model of the sgRNA abundances
(MAGeCK `test`, treatment LIB10/13/14 vs control LIB11/12/15), because the drug arm is strongly
under-templated and its rare guides have Poisson-dominated variance that a variance-homogeneous null
over-calls. Gene p-values were taken in the observed direction and adjusted with Benjamini-Hochberg
over the 230 genes; `drug_class` = sensitiser (FDR < 5%, depleted under drug), resistance
(FDR < 5%, enriched under drug), else none. Result: 18 sensitisers, 18 resistance genes, 194 none.
Validation on the reference genes: 10/10 established modifiers recovered with the expected class,
0/30 established neutral genes called. The FDR decision is not marginal - the largest FDR among the
calls is 0.028, the smallest among the non-calls 0.066 - and all 36 calls move in the same direction
in each of the three drug replicates. (A two-sided BH test gives the same list except TGT4571,
FDR 0.057; TGT4571 is not on the mailed list, so the erratum is unaffected.)

## 4. Erratum (`erratum.tsv`)
Every gene of `handover/hits_mailed.tsv` is listed with a verdict (SOP section 5):
`confirmed` where the released analysis classes the gene as sensitiser or resistance (26 genes),
`retracted` where it is classified as none (14 genes). No reagent had been ordered on the mailed
list, so the correction is editorial only. The mailed list has three systematic faults:

* it was computed on the withdrawn release-1 gene labels. Seven of the fourteen retracted rows name a
  gene that does not carry the sgRNAs which produced the mailed statistic at all; those guide pools
  belong, under release 2, to a different gene which *is* a modifier. The effect the colleague measured
  was real but misattributed, and the name that was mailed has to be withdrawn:

  | mailed (r1) name | verdict | r2 gene of the same 4 sgRNAs | released class |
  |---|---|---|---|
  | TGT4344 | retracted | TGT7853 | sensitiser |
  | TGT0462 | retracted | TGT7971 | sensitiser |
  | TGT3362 | retracted | TGT1631 | sensitiser |
  | TGT2478 | retracted | TGT4168 | sensitiser |
  | TGT4232 | retracted | TGT7793 | resistance |
  | TGT4662 | retracted | TGT5795 | sensitiser |
  | TGT3979 | retracted | TGT4331 | sensitiser |

  The remaining seven retracted rows (TGT7728, TGT5816, TGT1012, TGT2987, TGT0409, TGT1378, TGT4120)
  keep the correct name but are not modifiers at 5% FDR (their FDR is 0.49-0.84); their apparent
  depletion/enrichment comes from the under-templated drug libraries (guides with 3-25 reads, see
  e.g. TGT0409, whose only well-powered sgRNA does not change) and, for TGT1378/TGT3362-like rows,
  from the mailed ranking, which was cut at "top 40 by |z|" with a 3-sigma rule instead of a 5% FDR.
* the analysis behind the mailed list also dropped the drug replicate LIB15 (0 reads, counting bug),
  so it compared 2 vehicle against 2 drug libraries and never used the re-sequenced library;
* no confirmed row changes identity between releases (all four of its sgRNAs carry the same gene
  label in r1 and r2), and every mailed gene is either confirmed with FDR <= 0.023 or retracted with
  FDR >= 0.49 - there is no gene in the grey zone, so the erratum does not depend on the exact test.
