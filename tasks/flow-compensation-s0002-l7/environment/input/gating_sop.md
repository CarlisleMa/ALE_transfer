# Gating SOP IMM-FC-012 (PBMC lineage panel with T-cell activation read-out, 8 colours)

Marker names are those of `panel.tsv`. "+" means that the cell carries the marker, i.e. its signal
from that reagent is clearly above the signal of cells that do not carry it; "-" means that it does not.
Populations are nested: an event can only belong to a population if it belongs to the parent.

| population | parent | definition |
|---|---|---|
| cells | all_events | intact cells of any kind; sub-cellular debris is excluded |
| singlets | cells | events that consist of exactly one cell |
| live | singlets | Viability dye - (membrane-intact cells) |
| lymphocytes | live | mononuclear cells with low side scatter, CD14 - |
| monocytes | live | CD14 + |
| T_cells | lymphocytes | CD3 + |
| CD4_T | T_cells | CD4 + CD8 - |
| CD8_T | T_cells | CD4 - CD8 + |
| DP_T | T_cells | CD4 + CD8 + |
| DN_T | T_cells | CD4 - CD8 - |
| B_cells | lymphocytes | CD3 - CD19 + |
| NK_cells | lymphocytes | CD3 - CD19 - CD56 + |
| activated_CD4_T | CD4_T | HLA-DR + by the FMO gate of section 2 |
| activated_CD8_T | CD8_T | HLA-DR + by the FMO gate of section 2 |

`all_events` denotes every event recorded in the file.

## 1. Controls of the HLA-DR read-out

HLA-DR on T cells is a dim, continuous signal without a separate positive peak. Its gate is
therefore not set on the stained sample but on a fluorescence-minus-one control ("FMO HLA-DR"):
PBMC of the laboratory's reference pool, stained with the complete panel except the HLA-DR reagent
and processed like a sample.

## 2. FMO gate for HLA-DR + (populations activated_CD4_T and activated_CD8_T)

1. The FMO tube is analysed exactly like a sample: same compensation as the samples that were
   recorded with the same fluorescence-detector voltages, same gates down to CD4_T and CD8_T.
2. For each of the two parent populations P (CD4_T, CD8_T), the threshold is the 99th percentile of
   the compensated signal of the HLA-DR detector (linear scale, detector units) over the events of
   P in the FMO tube.
3. An event of P in a sample is HLA-DR + if its compensated HLA-DR-detector signal is greater than
   the threshold of P from the FMO tube that was recorded with the same fluorescence-detector
   voltages as the sample.
