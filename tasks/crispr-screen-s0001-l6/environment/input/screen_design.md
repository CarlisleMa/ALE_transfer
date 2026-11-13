# Screen design notes (wet lab)

- Cells: Cas9-stable clone (identity blinded by the client), lentiviral infection at MOI ~0.3, 500x library coverage until the arms were split.
- Library: VAL-2 secondary validation library in lentiGuide-Puro. 190 genes carried forward as candidate fitness genes from the primary genome-wide dropout screen and 40 genes picked as putative neutral; 4 sgRNAs per gene, plus 80 non-targeting sgRNAs. The reference genes for the compound class are part of these genes.
- Arms: T0 = day 3 after puromycin selection (2 replicates). The culture was then split into a vehicle arm (DMSO, 3 replicates) and a drug arm (compound at IC70, 3 replicates), both harvested 14 days later.
- Harvest: vehicle arm about 150000 cells per replicate went into the gDNA prep (about 150 cells per sgRNA). The drug arm was limited by the treatment: all surviving cells were taken, about 25000 per replicate (about 25 per sgRNA).
- gDNA PCR: one-step PCR with the P5 stagger primer mix (8 forward primers, 0-7 nt stagger) and one P7 index primer per library.
- Amplicon, sense strand 5'->3': ...TCTTGTGGAAAGGACGAAACACCG [20-nt spacer] GTTTTAGAGCTAGAAATAGCAAGTTAAAAT... (sgRNA scaffold).
- Sequencing run A: NextSeq 550, 75-cycle single read.
- Sequencing run B: NovaSeq, paired-end lane shared with other users.
