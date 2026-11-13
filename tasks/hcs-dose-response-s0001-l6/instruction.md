You are the image-analysis scientist of a compound-profiling group. A 96-well plate of cells was
treated with eight test compounds in dose series, fixed, stained for nuclei and for a cytoplasmic
marker protein, and imaged on an automated microscope. The project team needs per-well
measurements and, for every test compound, its potency on the marker and on cell number. The
group's SOP requires that the cell-counting procedure is checked against the lab's manual counts
of this plate before it is used.

Task directory:
- `/task`

Inputs (read-only) under `/task/input`:
- `images/<well>_ch1.tif`, `images/<well>_ch2.tif` - one field per well (96 wells), 16-bit
  grayscale TIFF, 320 x 320 pixels. Channel 1 = nuclear stain, channel 2 = marker.
- `flatfield/flat_ch1.tif`, `flatfield/flat_ch2.tif`, `flatfield/dark.tif` - reference frames
  recorded with the same optical path (uniform fluorescent slide per channel; closed shutter).
- `plate_layout.tsv` - columns `well`, `compound`, `dose_uM`, `role` (`compound`, `vehicle` or
  `reference`).
- `annotations/<well>.tsv` - the lab's manual point annotation of 8 fields of this plate, one
  file per annotated well: tab-separated, header `x_px`, `y_px`, one row per cell that was counted
  in that field under SOP section 1 (a click on its nucleus). Pixel coordinates in the image of
  that well: `x_px` = column, `y_px` = row, 0-based, origin at the top-left pixel. The manual count
  of a field is the number of rows.
- `assay_notes.md` - hand-over notes of the wet lab and of the imaging operator.
- `analysis_sop.md` - the group's SOP. It defines the two readouts and what is counted as a cell,
  the validation of the counting procedure, image QC, the response scale, the curve model, the
  potency and the meaning of `inactive`. These definitions are binding.

Available software (on `PATH`): Python 3.12 with numpy, scipy, pandas, scikit-image, tifffile,
matplotlib and Cellpose 3 (CPU only; the pretrained `nuclei` and `cyto3` models are installed,
nothing can be downloaded). There is no internet access.

Your task:
1. Quantify every well: the number of cells in the field and the marker signal per cell (SOP
   section 1), and decide whether the well passes image QC (SOP section 3).
2. Validate your counting procedure against the manual counts as SOP section 2 requires, and
   report the comparison for the procedure you finally used.
3. For every test compound, determine the potency on the marker (IC50) and on cell number (GI50)
   as defined in SOP sections 4 and 5, or report that the compound is inactive on that readout.

Required outputs under `/task/output`:
- `well_measurements.tsv` - tab-separated, header `well`, `n_cells`, `marker_per_cell`, `qc_pass`;
  one row for every well of `plate_layout.tsv`. Extra columns are allowed.
  - `n_cells`: integer, SOP readout R1 of the field as counted by your procedure (for the
    annotated wells too).
  - `marker_per_cell`: SOP readout R2 of the well in grey levels (`NA` if the field holds no cell).
  - `qc_pass`: `1` if the well passes the SOP's image QC, `0` if it fails. Wells that fail are
    still listed, with the values your analysis produced for them.
- `validation.tsv` - tab-separated, header `field`, `manual_count`, `pipeline_count`; one row for
  each annotated field and no other rows. `field` is the well id, `manual_count` the number of
  annotated cells, `pipeline_count` the count of your final procedure for that field (the same
  number as `n_cells` of that well in `well_measurements.tsv`).
- `dose_response.tsv` - tab-separated, header `compound`, `ic50_marker_uM`, `gi50_count_uM`; one
  row per test compound (no rows for vehicle or reference). Each potency is either a number in
  micromolar or the literal `inactive`, as defined in SOP section 5.

Constraints:
- Do not modify anything under `/task/input`; write deliverables only under `/task/output`.
- Do not ask for confirmation. Execute directly.
