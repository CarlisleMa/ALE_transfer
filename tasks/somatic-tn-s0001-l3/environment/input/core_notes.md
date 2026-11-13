# Sequencing core hand-over: run 214, panel CGP-22

## Assay

- Hybrid-capture panel CGP-22: 13 target regions, 46.3 kb in total (`targets.bed`). The panel
  reference `reference/ref.fa` (contig `chr22_panel`) is the hg38 chr22 slice the baits were designed on.
- Libraries: enzymatic end repair / A-tailing, ligation of dual-indexed adapters, PCR amplification,
  one pooled capture reaction for all six libraries.
- Sequencing: NovaSeq 6000, flow cell HKF3WDSX7 lane 2, 2 x 100 cycles. Demultiplexed with bcl2fastq;
  adapters were trimmed by the instrument software. FASTQ files are otherwise as produced.
- Requested mean depth on target: tumour 180x, normal 90x.

## Libraries

| sample_id | specimen | DNA input (ng) | shearing | PCR cycles | i7+i5 index |
|---|---|---|---|---|---|
| CG24-7410 | fresh-frozen biopsy | 80 | batch A | 8 | CTCTATAC+AAATAGAA |
| CG24-1366 | peripheral blood (EDTA) | 200 | batch A | 8 | CTAATGAC+TGTGTCGA |
| CG24-6425 | peripheral blood (EDTA) | 200 | batch A | 8 | GGCATGTG+TCGCCCTT |
| CG24-4187 | fresh-frozen biopsy | 100 | batch A | 8 | CATTAACT+CACGCTGT |
| CG24-9024 | FFPE section (block from 2018) | 80 | batch A | 10 | TAAGTGTC+TGCTGTGG |
| CG24-3973 | peripheral blood (EDTA) | 200 | batch A | 8 | GCCACTTC+TCACTCCC |

## Remarks

- Shearing: Covaris LE220, low-EDTA TE, standard setting, one batch.
- FFPE sections were deparaffinised and extracted with the FFPE DNA kit of the pathology lab.
- Tumour cell content estimated by the pathologist on the adjacent H&E section: PT-0452 70%, PT-0740 60%, PT-0760 70%.
- All libraries passed the core's yield and fragment-size checks (mean insert 260-280 bp).
