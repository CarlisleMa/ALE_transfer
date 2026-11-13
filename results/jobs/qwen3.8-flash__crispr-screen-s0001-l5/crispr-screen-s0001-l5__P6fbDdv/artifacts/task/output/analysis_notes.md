# VAL-2 drug-modifier screen - final analysis (replaces the mailed draft)

## Pipeline
1. **Counting** (`counts.tsv`, `qc.json`).
   Reads were counted for a library sgRNA when the read carries that 20-nt spacer in the
   amplicon context of the wet-lab notes (`...TCTTGTGGAAAGGACGAAACACCG` + spacer +
   `GTTTTAGAGCTAGAAATAGCAAGTT...`), searched in **both orientations**. Because the P5
   stagger (0-7 nt) and the constant regions carry sequencing errors, a read whose spacer
   matches a library sgRNA exactly but whose anchor region is error-free nowhere is still
   assigned (sliding 20-mer); no read ever matched two different sgRNAs, so there are no
   ambiguous assignments. Mapping rate is 90.5-92.6 % in all eight libraries.
   *LIB15* is the Read 2 file of the NovaSeq re-sequencing (run B) and is the reverse
   complement of the amplicon; the draft script only searched the forward scaffold, which is
   why it reported **0** reads for that library. It is now counted and used.
2. **Arms.** The day-14 tubes were re-racked and the sheet was filled from the new rack
   positions; two day-14 labels are wrong. Groups used (see `qc.json`, `arm`):
   - T0: LIB8, LIB9
   - vehicle: LIB11, LIB12, **LIB15**
   - drug: **LIB10**, LIB13, LIB14
   Evidence: (i) the drug arm was harvested at ~25 cells per sgRNA, so 11-13 library sgRNAs
   are absent from LIB10/LIB13/LIB14 while LIB11/LIB12/LIB15 lose none - LIB10 behaves like
   the drug arm and LIB15 like the vehicle arm; (ii) the clustering of log-ratio profiles puts
   LIB10 with LIB13/LIB14 and LIB15 with LIB11/LIB12; (iii) the established reference genes
   only separate in this direction - over the ten possible 3-vs-3 splits of the day-14
   libraries the split above gives the largest separation (resistance median +1.69 vs
   sensitiser median -1.42, score 3.11 against 1.02 for the next best split), and it
   reproduces **40/40** genes of `reference_genes.tsv` (the sheet's own labels reverse the
   sign of all ten reference modifiers).
3. **Effect sizes** (`gene_summary.tsv`), per SOP FG-017 section 4:
   raw counts were converted to abundance per library by dividing by the mapped-read depth
   (CPM); the sgRNA-level log2 ratio is log2((mean CPM group A + 1)/(mean CPM group B + 1));
   the gene value is the **median over its 4 sgRNAs**; both quantities are shifted so that the
   median non-targeting sgRNA is 0, i.e. a knockout with no effect reads 0.
   `lfc_fitness` = vehicle day 14 vs T0, `lfc_drug` = drug day 14 vs vehicle day 14
   (negative = depleted under drug).
4. **Significance / class.** sgRNA-level counts were tested with MAGeCK `test` (negative
   binomial, median normalisation) for drug vs vehicle; the gene-level two-sided p-value is
   2 x min(p_depleted, p_enriched) from the RRA output, and Benjamini-Hochberg over the 230
   targeted genes gives the FDR. `drug_class` = `sensitiser` if FDR < 0.05 and `lfc_drug` < 0,
   `resistance` if FDR < 0.05 and `lfc_drug` > 0, else `none`. 35 modifiers (18 sensitiser,
   17 resistance); none of the 30 reference "none" genes is called and all 10 reference
   modifiers are recovered with the right sign. A non-targeting-resampling null was rejected
   because it is too narrow for the low-depth drug arm (it would have called 8 reference
   "none" genes).

## Erratum (`erratum.tsv`)
The mailed list was produced with (a) withdrawn release-1 gene labels (80 sgRNAs / 20 genes
mislabelled), (b) LIB15 silently dropped because it counted as 0, and (c) the swapped day-14
arms, so its z-scores are not interpretable. 26 of the 40 mailed genes are confirmed as
modifiers, 14 are retracted - including the top two entries (TGT4344, TGT0462) and
TGT4120, a reference gene with established neutral behaviour.
