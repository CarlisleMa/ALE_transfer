# Level-4 spec: `crispr_screen` - taking over a drug-modifier screen (H + P + Q)

Pilot at level 3 (seed 1): Claude Haiku scored 1.0 in 2 min 1 s with 17 tool calls. The lab notes
named all three decision points and a per-sample mapping rate revealed two of them.

## Scenario
Same validation library and amplicon design. The experiment is now a chemogenomic screen: T0 x2,
vehicle day 14 x3, drug day 14 x3 (8 FASTQ files). A predecessor started the analysis, mailed a
draft hit list to the PI and left. The analyst takes the project over and delivers the genes that
modify the drug response.

## Truth (planted before any read exists)
Per gene two independent latent effects: fitness (vehicle vs T0) and drug effect (drug vs vehicle;
sensitiser < 0, resistance > 0, none). Drug effects include a moderate tier (|log2FC| about 0.8-1.5)
next to the strong one, so that power and false-positive control both matter (Q). Fitness genes and
drug genes are drawn independently, so "depleted in the drug arm compared with T0" does not
identify sensitisers (P). Keep blinded gene ids.

## Inherited project: `/task/input/handover/`
Count script, count table, analysis script, draft gene table, the draft hit list that was mailed,
a hand-over note. Produce them by running **mutated copies of the reference pipeline** on the task
data, so every intermediate file is exactly what the defective code yields. The SOP
(`screen_sop.md`) holds the audit checklist: 8-10 items, of which 3-5 are defective for a given
seed and the others genuinely fine (include fine items that *look* unusual). Candidate operators:
library annotation version (a v2 file with a changelog supersedes v1, ~8% of guides re-assigned);
sample-to-file assignment (lexicographic vs numeric column order, or two files exchanged between
arms); fixed-offset read parsing; the Read-2-only sample; a bottlenecked replicate that the SOP's
stated replicate criterion excludes; normalisation reference; contrast direction; a stale draft
table that was not produced from the count table next to it. At least two defects must be
detectable only by recomputing, not by reading the code. The three level-3 decision points stay.

## Deliverables
`counts.tsv`, `gene_summary.tsv` (`gene`, `lfc_fitness`, `lfc_drug`, `drug_class` in
`sensitiser`/`resistance`/`none`), `audit.tsv` (`item`, `verdict`), `erratum.tsv` (one row per gene
of the mailed draft list: `gene`, `verdict` in `confirmed`/`retracted`), `qc.json`.

## Scoring (set the bands from measurements on two seeds)
counts (minimum over samples) 0.15; `lfc_drug` error 0.20; `drug_class` (multi-class MCC) 0.30;
audit MCC 0.20; erratum MCC 0.15. Pass 0.9.

## Battery
`l3_reference` (level-3 method: counts + fitness test, no drug contrast logic, no audit);
`trust_inherited` (finish from the inherited count table); one variant per active defect left
unfixed; `drug_vs_t0`; audit "everything defective" and "everything ok"; a perturbation just
beyond tolerance. Probes: planted truth, within tolerance, independent second solution.
