# Family spec: `flow_compensation` - multicolour flow cytometry: compensation, gating, frequencies

ALE subdomain: 3.3 Cell & Imaging Biology (cytometry). No counterpart among ALE public tasks.

## Scenario
An immunology core must turn raw 8-colour PBMC acquisitions into a spillover matrix and a table of
population frequencies following the lab's gating SOP.

## Truth model (pure forward simulation; every event carries a hidden label)
- Populations: debris, doublets, dead cells, live singlet lymphocytes -> T (CD4, CD8, DN, DP), B,
  NK; monocytes. Per-sample proportions drawn from the seed.
- True fluorochrome signals per event from population-specific log-normal distributions with
  positive/negative populations separated by >= 1.5 decades after compensation; scatter channels
  (FSC-A, FSC-H, SSC-A) with doublets off the FSC-A/FSC-H diagonal.
- Planted spillover matrix with realistic structure (e.g. PE -> PE-Cy7, APC -> APC-Cy7, BV421 ->
  BV510 up to 0.2-0.4). Detector data = true signal x spillover + photon and electronic noise.
- Files: FCS 3.1 for each sample (>= 5 samples), one single-stain control per fluorochrome
  (compensation beads or cells) and an unstained control.

## Inputs
`fcs/*.fcs`, `panel.tsv` (detector, fluorochrome, marker), `controls.tsv`, `acquisition_notes.md`
(facts only: compensation not applied at acquisition, which controls were recorded, any
re-acquisition with changed detector voltage), `gating_sop.md` (population definitions as marker
logic and the parent of each population). Tools: python with flowkit/flowio/fcsparser,
scikit-learn, numpy/scipy/pandas.

## Outputs
`spillover.csv` (fluorochrome x detector, diagonal 1), `frequencies.tsv` (`sample, population,
parent, count, percent_of_parent` for the fixed population list in the SOP).

## Traps (level 3 = all)
- compensation_required: the `$SPILLOVER` keyword in the sample files is absent or the acquisition
  default (identity); it must be computed from the single-stain controls.
- doublets: without FSC-A/FSC-H discrimination, doublets inflate double-positive populations.
- dead_cells: dead cells bind antibody non-specifically; without the viability gate lineage
  frequencies shift.
- per_sample_shift: one sample acquired with a different detector gain -> one global threshold set
  misclassifies it.

## Scoring sketch
gates: files parse, all samples x populations present, percentages in [0,100]. Components:
spillover max abs error (0.25); per (sample, population) agreement of percent_of_parent with the
planted truth within max(abs, rel) tolerance, aggregated so a wholly wrong sample is not averaged
away (0.60); live singlet lymphocyte counts within tolerance (0.15). Tolerances from the measured
spread of two correct, different gating strategies (threshold/manual-style vs model-based).

## Naive variants
no compensation; identity/embedded matrix; no doublet gate; no viability gate; one threshold set
for all samples.
