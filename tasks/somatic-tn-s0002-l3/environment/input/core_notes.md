# Sequencing core hand-over: run 214, panel CGP-22

## Assay

- Hybrid-capture panel CGP-22: 13 target regions, 47.5 kb in total (`targets.bed`). The panel
  reference `reference/ref.fa` (contig `chr22_panel`) is the hg38 chr22 slice the baits were designed on.
- Libraries: enzymatic end repair / A-tailing, ligation of dual-indexed adapters, PCR amplification,
  one pooled capture reaction for all six libraries.
- Sequencing: NovaSeq 6000, flow cell HKF3WDSX7 lane 2, 2 x 100 cycles. Demultiplexed with bcl2fastq;
  adapters were trimmed by the instrument software. FASTQ files are otherwise as produced.
- Requested mean depth on target: tumour 180x, normal 90x.

## Libraries

| sample_id | specimen | DNA input (ng) | shearing | PCR cycles | i7+i5 index |
|---|---|---|---|---|---|
| CG24-6922 | peripheral blood (EDTA) | 200 | batch A | 8 | AGACACTC+CTGTAACA |
| CG24-3081 | fresh-frozen biopsy | 150 | batch A | 8 | GGACACCC+GACATTCA |
| CG24-4358 | peripheral blood (EDTA) | 200 | batch A | 8 | TCTGTTTA+CGGACCCG |
| CG24-2309 | FFPE section (block from 2019) | 200 | batch A | 10 | GAAAGCAC+GACAACTT |
| CG24-1557 | peripheral blood (EDTA) | 200 | batch A | 8 | AAGCAACA+GGATAGCC |
| CG24-9370 | fresh-frozen biopsy | 50 | batch A | 8 | AATACACG+GCCGGTGG |

## Remarks

- Shearing: Covaris LE220, low-EDTA TE, standard setting, one batch.
- FFPE sections were deparaffinised and extracted with the FFPE DNA kit of the pathology lab.
- Tumour cell content estimated by the pathologist on the adjacent H&E section: PT-0204 70%, PT-0476 80%, PT-0789 60%.
- All libraries passed the core's yield and fragment-size checks (mean insert 260-280 bp).
