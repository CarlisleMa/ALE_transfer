# IMM-FC-012 PBMC lineage panel - analysis notes

## 0. State of the raw files
* FCS 3.1, list-mode, 12 parameters (FSC-A, FSC-H, SSC-A, 8 fluorescence areas, Time), `$MODE=L`.
* Detector voltages are stored per parameter in `$PnV`.
* The event data are **uncompensated** (single-stain tubes show their full spillover), even though every file
  carries `APPLY COMPENSATION = TRUE` and a `$SPILLOVER` keyword.
* The `$SPILLOVER` keyword is **identical in all 19 files** and is therefore *not* the matrix that describes
  each individual file: it lists PerCP-Cy5.5 -> YG780-A = 0.093, which matches neither the bead tube
  (0.152) nor the lot-corrected cell tube (0.253), and it does not change between the 641 V and the 703 V
  acquisitions. It was ignored; every matrix below was recomputed from the control tubes.
  Cross-check: regressing the raw YG780-A signal of the HLA-DR FMO tubes (which contain every reagent except
  PE-Cy7) on the compensated fluorochrome signals gives PerCP-Cy5.5 -> YG780-A = 0.252 (641 V tube) and
  0.476 (703 V tube), i.e. the control-derived values, not the keyword values.

## 1. Spillover matrices (`spillover.csv`)
Row *f*, column *d* = (mean signal of the f-positive control population in d) - (mean of the f-negative
population in d), divided by the same difference in f's primary detector, background-free, diagonal = 1.

Control selection (`controls.tsv`):
* one candidate per fluorochrome was dropped on purpose: `ctrl_PerCPCy55.fcs` is antibody-capture beads of
  reagent lot 2807268, while `panel.tsv` uses lot **3066970** for CD8 PerCP-Cy5.5. That lot was acquired on
  cells in `ctrl_PerCPCy55_cells.fcs`, which is what the samples contain, so the **cell** tube was used
  (spill into YG780-A 0.253 instead of the bead tube's 0.152).
* `ctrl_PECy7.fcs` (YG780-A = 641 V) is recorded at the wrong settings for the last two samples; the
  rerun tube `ctrl_PECy7_rerun.fcs` (703 V) is lot-identical and was used for those files.
* Zombie Aqua was controlled on cells (1:1 live/heat-killed PBMC), all other conjugates on beads.

Voltage groups of the fluorescence detectors (`$PnV`):
* group A, YG780-A = 641 V: PBMC_D617, PBMC_D685, PBMC_D735, PBMC_D829;
* group B, YG780-A = 703 V: PBMC_D858, PBMC_D960 (D858 additionally has different FSC/SSC voltages, which
  does not change the fluorescence spillover).
All other fluorescence detectors have identical voltages in both groups.

Group B matrix from the group-A matrix: the signal gain of YG780-A between the two settings was measured
from the PE-Cy7 tube pair (positive-minus-negative excess in YG780-A: 30977 at 641 V -> 58681 at 703 V,
k = 1.8944; the 703/641 = 1.097 voltage ratio is *not* the signal gain). Every signal recorded in YG780-A
is k times larger at 703 V, hence
* A_B[f, d!=YG780-A] = A_A[f, d] for f != PE-Cy7,
* A_B[f, YG780-A] = k * A_A[f, YG780-A] for f != PE-Cy7,
* A_B[PE-Cy7, d!=YG780-A] = A_A[PE-Cy7, d]/k, A_B[PE-Cy7, YG780-A] = 1.
The PE-Cy7 row of group B agrees with the directly measured rerun row to <1%.

## 2. Compensation and gating (`frequencies.tsv`)
Samples and both FMO tubes were compensated with the matrix of their own voltage group, then gated.
All gate boundaries were placed automatically in the empty band (or smoothed minimum) between the two
populations of the relevant channel, so the counts do not depend on the exact position of the cut.
* `cells`   : FSC-A above the debris gap (FSC voltages differ between files, so the cut is found per file).
* `singlets`: FSC-H/FSC-A above the gap (empty band 0.72-0.84 in every file).
* `live`    : compensated V525-A (Zombie Aqua) below the gap of the singlet population.
* `lymphocytes` = live AND SSC-A below the lymphocyte/monocyte gap AND CD14 -; `monocytes` = CD14 +.
* `T_cells` = lymphocytes AND CD3 +; `CD4_T` = CD3+ CD4+ CD8-; `CD8_T` = CD3+ CD4- CD8+;
  `DP_T` = CD4+ CD8+; `DN_T` = CD4- CD8-; `B_cells` = lymphocytes CD3- CD19+;
  `NK_cells` = lymphocytes CD3- CD19- CD56+. CD4/CD8/CD19/CD56 cuts are taken inside T_cells / CD3-
  lymphocytes (hierarchical), all in empty bands.
  CD4+CD8+DP+DN = T_cells exactly for every sample; lymphocytes+monocytes = live for every sample.
* `activated_CD4_T`, `activated_CD8_T`: FMO gate of section 2 of the SOP. The FMO tube was compensated with
  the matrix of its own voltage group (fmo_HLADR.fcs, 641 V, for the group-A samples; fmo_HLADR_rerun.fcs,
  703 V, for PBMC_D858 and PBMC_D960) and gated identically down to CD4_T and CD8_T; the threshold is the
  99th percentile of the compensated YG780-A signal of that population (linear units):
  group A: CD4_T 302.0, CD8_T 431.7; group B: CD4_T 599.9, CD8_T 792.9 (ratio 1.98 ~ the k gain factor).
  A sample event of P is HLA-DR + if its compensated YG780-A signal is strictly greater than that threshold.

* `count` = number of events of the sample's FCS file in the population; `percent_of_parent` =
  100 x count / parent count (`all_events` = `$TOT` of the file).
