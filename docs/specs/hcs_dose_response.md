# Family spec: `hcs_dose_response` - high-content screening plate: images to potency values

ALE subdomain: 3.3 Cell & Imaging Biology. (ALE's translocation task prescribes the ratio and the
threshold; here the analyst must get from raw plate images to per-compound potencies and the
difficulty is plate-level confounding and readout choice.) Uses a professional model (Cellpose).

## Scenario
A screening group imaged a 2-channel (nuclei, cytoplasmic marker) plate with 8 compounds in
8-point dose series plus vehicle and reference-inhibitor wells. Deliver per-well measurements and,
per compound, the potency on the marker (IC50) and on cell number (GI50), or "inactive".

## Truth model (forward-rendered images; every cell and every well parameter is planted)
- Per compound: planted IC50 on per-cell marker intensity and GI50 on cell number (some compounds
  active on one readout only, some inactive, potencies drawn from the seed). 4-parameter logistic.
- Per well: planted cell count (dense wells contain touching nuclei), per-cell marker intensity.
- Rendering: nuclei as textured ellipses, cytoplasm around them, PSF blur, shot and read noise,
  field illumination roll-off, plate-level staining gradient and edge effect, a few out-of-focus
  wells, one compound that is autofluorescent in the marker channel at high dose.
- Plate layout places some top doses in edge wells; vehicle controls are spread over the plate.

## Inputs
`images/<well>_ch1.tif`, `_ch2.tif`, `plate_layout.tsv` (well, compound, dose_uM, role),
`flatfield/` reference images if used, `assay_notes.md` (stain, objective/pixel size, known plate
handling facts), `analysis_sop.md` (definitions of the two readouts, of "inactive", of the QC
exclusion rule, curve model). Tools: python with cellpose (CPU), scikit-image, scipy, tifffile.

## Outputs
`well_measurements.tsv` (`well, n_cells, marker_per_cell, qc_pass`),
`dose_response.tsv` (`compound, ic50_marker_uM, gi50_count_uM`; the literal `inactive` where the SOP's
activity criterion is not met).

## Traps
plate_position_effects (edge + gradient, must be corrected with the distributed vehicle wells);
readout_confound (a cytotoxic compound lowers total well signal but not per-cell marker);
autofluorescent_compound (needs local background handling); focus_qc (blurred wells excluded);
touching_nuclei (thresholding alone undercounts dense wells).

## Scoring sketch
gates: files parse, every well/compound present. Components: per-well cell counts (0.20); QC flags
(0.10); marker potency: activity classification and |log10(IC50_est/IC50_true)| (0.45); count
potency likewise (0.25). Bands from two correct pipelines (Cellpose-based and classical
watershed-based segmentation) and the naive ones (no plate correction; whole-well intensity
readout; global threshold counting; no QC exclusion).
