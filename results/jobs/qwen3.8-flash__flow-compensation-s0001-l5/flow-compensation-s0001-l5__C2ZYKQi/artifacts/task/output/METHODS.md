# Method notes (IMM-FC-012_PBMC_lineage)

## Spillover
The compensation stored in the FCS headers is the 2024 experiment-template matrix and is stale
(voltages were re-optimised 04-JAN-2025, no new compensation setup was acquired), so every matrix was
recomputed from the single-stain controls of this session. Data are raw linear detector values
(PnE = 0,0; gain 1), so no untransform step is needed.

For each control tube the positive and negative carrier events were separated by Otsu on the primary
detector; row f = (median_pos - median_neg) / (median_pos - median_neg)_primary, i.e. background-free
and normalised to 1 in the primary detector (mean- and regression-based estimators agree to <=0.002).

Control choices:
* PerCP-Cy5.5: ctrl_PerCPCy55_cells.fcs (PBMC, reagent lot 3066970 = the lot of the sample cocktail).
  ctrl_PerCPCy55.fcs uses lot 2807268 and its spill into YG780-A is 0.152 instead of 0.253, so it is
  the wrong-lot tube and is not used.
* PE-Cy7: ctrl_PECy7_rerun.fcs for the samples recorded at 703 V, ctrl_PECy7.fcs for the 641 V samples.
* Zombie Aqua: cell-based single stain (as in the panel).

Voltage groups (from PnV in each file): all detectors are identical in every file except YG780-A,
641 V (PBMC_D617/685/735/829) and 703 V (PBMC_D858/960). From the same PE-Cy7 bead tube re-recorded
at 703 V, the YG780-A gain factor rho = 1.89141 ((703/641)^6.9); all other fluorochromes keep their
measured spill into the other detectors, only their YG780-A column is multiplied by rho, and the
PE-Cy7 row is taken from the rerun. Both matrices have the required diagonal of 1.

## Gating
Compensation: a = M . inv(S) with S the matrix above. FSC-A was brought to a common scale by
(420/P1V)^7 because PBMC_D858 was re-acquired at 348 V (fluorescence voltages unchanged).

Thresholds sit in wide low-density gaps that are present in all six samples and both FMO tubes:
* cells       FSC-A    > 40000   (valley empty 34.4k-47.6k; debris excluded)
* singlets    FSC-H/FSC-A > 0.80 (empty 0.743-0.838)
* live        V525-A   < 2500    (gap 1.5k-3.4k)
* lymphocytes SSC-A    < 50000 and CD14-
* monocytes   CD14 (B530-A) > 3000 (gap 0.42k-7.9k)
* T_cells     CD3 (V450-A) > 5000  (gap 0.48k-17.8k)
* CD4_T CD4+CD8-  CD4 (YG586-A) > 6000 (gap 0.6k-13.5k), CD8 (B710-A) < 2000 (gap 1.0k-4.0k)
* CD8_T, DP_T, DN_T from the same two thresholds; the four quadrants partition T_cells exactly.
* B_cells     CD3- CD19 (R780-A) > 1500 (gap 0.72k-2.6k)
* NK_cells    CD3- CD19- CD56 (R670-A) > 3000 (gap 1.3k-6.0k)
Moving any of these thresholds anywhere inside its gap changes the counts by <0.1 %.

## HLA-DR (activated_CD4_T, activated_CD8_T)
Per SOP section 2, the threshold is the 99th percentile of the compensated YG780-A signal over CD4_T
and CD8_T of the FMO tube recorded with the same fluorescence-detector voltages as the sample:
* fmo_HLADR.fcs      (641 V): CD4_T 299.52, CD8_T 424.43  -> PBMC_D617/685/735/829
* fmo_HLADR_rerun.fcs(703 V): CD4_T 598.21, CD8_T 790.27  -> PBMC_D858/960
The two thresholds scale like the YG780-A gain, which is why the SOP requires matching voltages.
