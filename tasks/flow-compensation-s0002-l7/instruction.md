You are the cytometry data analyst of an immunology core facility. A study team acquired PBMC
samples with an 8-colour immunophenotyping panel and needs, for every sample, the spillover matrix
that applies to it and the size of every population in the lab's gating SOP.

Task directory:
- `/task`

Inputs (read-only) under `/task/input`:
- `fcs/*.fcs` - FCS 3.1 list-mode files as exported by the acquisition software: 6 PBMC
  samples and the 13 control tubes of the session
- `samples.tsv` - sample sheet (`sample_id`, `fcs_file`, `donor`, `acquisition_order`)
- `controls.tsv` - control tubes (`fcs_file`, `control_type`, `fluorochrome`, `detector`, `carrier`,
  `reagent_lot`). For a single-stain control the columns name the one reagent in the tube; for an
  FMO control they name the one reagent of the panel that was left out.
- `panel.tsv` - staining panel of the samples (`detector`, `fluorochrome`, `marker`, `reagent_lot`)
- `gating_sop.md` - the lab's gating SOP: population names, parent of each population, definitions
  and the control-based HLA-DR gate
- `acquisition_notes.md` - the operator's notes from the acquisition session

Available software: Python 3.12 (on `PATH`) with numpy, scipy, pandas, scikit-learn, FlowKit,
FlowIO and fcsparser. There is no internet access.

Your task:
1. For every sample in `samples.tsv`, determine the spillover matrix that describes the signals in
   its file.
2. For every sample in `samples.tsv`, identify every population of `gating_sop.md` as defined there
   and report its size.

Required outputs under `/task/output`:
- `spillover.csv` - comma-separated. Header: `sample`, `fluorochrome`, followed by the eight
  fluorescence detectors named as in `panel.tsv`. Exactly one row for every `sample_id` of
  `samples.tsv` and every fluorochrome of `panel.tsv` (names as written there). The eight rows of a
  sample are the spillover matrix of that sample: the entry in row *f*, column *d* is the signal
  that fluorochrome *f*, as it is present on the cells of that sample, produces in detector *d* at
  the detector settings with which that sample was recorded, as a fraction of the signal it
  produces in its own primary detector at the same settings (background-free; the entry of each
  fluorochrome in its own detector is 1).
- `frequencies.tsv` - tab-separated with the columns `sample`, `population`, `parent`, `count`,
  `percent_of_parent`; exactly one row for every `sample_id` of `samples.tsv` and every population
  of `gating_sop.md` (names exactly as written there). Extra columns are allowed.
  - `parent`: the parent population named in the SOP.
  - `count`: number of events of the sample's FCS file that belong to the population (integer).
  - `percent_of_parent`: 100 x `count` / (count of the parent population). For the parent
    `all_events` the denominator is the total number of events in the FCS file.

Constraints:
- Do not modify anything under `/task/input`; write deliverables only under `/task/output`
  (additional files there are ignored).
- Every sample in the sample sheet belongs to the study and must be reported.
- Do not ask for confirmation. Execute directly.
