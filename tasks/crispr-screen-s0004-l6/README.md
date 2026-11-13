# Pooled CRISPR drug-modifier screen: taking over an unfinished analysis (level 6: pointers removed)

Take over a pooled CRISPR drug-modifier screen from a predecessor: recount staggered amplicon reads, establish which libraries received drug, call sensitisers and resistance genes to a stated release criterion, audit the inherited analysis and write the erratum for a mailed hit list. Level 6: same inputs and planted truth as level 4; the text no longer says what to look out for (see the task card).

- family: `crispr_screen`  seed: `4`  level: `6`
- ALE subdomain: 3.2 Genomics & Sequence Analysis
- software: cutadapt, MAGeCK, Python, numpy/scipy/pandas
- a-priori difficulty: 0.675 (full-spectrum)  |  pass threshold: 0.9
- capabilities trained: read-structure reasoning from lab notes, per-sample QC before analysis, auditing inherited work by recomputation, trusting data over metadata, self-checking calls against in-data controls, choosing statistics that match the noise, exact output contracts

## Undisclosed decision points (not shown to the agent)

- **staggered_amplicon** - The spacer starts 24-31 nt into the read; a fixed-offset parser assigns ~1/8 of reads.
  - discoverable by: screen_design.md names the P5 stagger mix; per-library mapping rate collapses with a fixed offset
  - trains against ALE failure mode: Approach / wrong strategy
- **read2_only_sample** - LIB8 is the reverse complement of the amplicon; forward search maps ~0% of it.
  - discoverable by: screen_design.md says only Read 2 was delivered; read names end in ' 2:N:0'; mapping rate ~0 for one library
  - trains against ALE failure mode: Understanding / domain knowledge gap
- **hit_enriched_library** - 83% of genes are true fitness genes. Median/total normalisation shifts every lfc_fitness by > 1 log2 unit.
  - discoverable by: screen_design.md describes a validation library of carried-forward hits; the SOP states the scale (non-targeting sgRNAs at 0)
  - trains against ALE failure mode: Understanding / domain knowledge gap
- **H_inherited_project** - The predecessor's count table, gene table and mailed list carry 5 of 8 checklist defects (library_release, arm_labels, control_normalisation, drug_baseline, contrast_sign); two leave no trace in the scripts.
  - discoverable by: SOP section 6 lists the checklist; every item is decided by recounting / re-running the inherited scripts and comparing with the inherited tables
  - trains against ALE failure mode: Verification / unverified reuse of prior work
- **P_sample_sheet_arms** - LIB8 (sheet: vehicle) received drug and LIB11 (sheet: drug) received vehicle. With the sheet arms the drug contrast is diluted to one third.
  - discoverable by: SOP section 3 states the criterion (a library must agree better with its own arm); library-library correlation of centred log abundance, or the reference modifiers, separate the arms 3:3; notes mention the re-racked tubes
  - trains against ALE failure mode: Verification / trusting metadata over data
- **P_drug_vs_t0** - Fitness and drug effects are planted independently in a hit-enriched library: 'depleted in the drug arm compared with T0' returns the fitness genes, not the sensitisers.
  - discoverable by: instruction and SOP section 4 define the drug effect against vehicle; reference `none` genes with strong fitness effects are called by the T0 contrast
  - trains against ALE failure mode: Approach / wrong strategy
- **Q_modifier_calls** - Drug-arm libraries come from ~25 cells per sgRNA, so the noise of depleted sgRNAs is several times that of the non-targeting sgRNAs: a null taken from the non-targeting sgRNAs calls dozens of fitness genes, the level-3 rank test both over- and under-calls.
  - discoverable by: SOP section 5 states the release criterion on reference_genes.tsv (30 `none` genes over the whole fitness range, 10 modifiers); notes give the cell numbers at harvest; replicate scatter grows with depletion
  - trains against ALE failure mode: Verification / no self-check of result quality

## Judge

**Level 6.** Inputs and planted truth are those of `crispr-screen-s0004-l4`; only text was taken away:

- S2: the instruction to confirm the arm of every day-14 library from the data, and the sentence that libraries are grouped by the treatment they actually received (SOP section 3, instruction)
- S1: the release criterion on the reference genes (SOP section 5); reference_genes.tsv stays in the inputs
- S3: the eight-item hand-over audit checklist (SOP section 6) and the deliverable audit.tsv; its weight (0.20) is redistributed proportionally over the other four components
- S5: wet-lab notes without the two incident reports (the day-14 tubes re-racked before the sample sheet was filled in; the library that failed QC and came back as a Read-2 file). The sequencing run of every library is in the sample sheet
- the `arm` field of qc.json is no longer asked for (it was never scored; its description told the reader that the analysed group can differ from the sample sheet)

Fairness no longer rests on a stated criterion but on three things that are checked for the level-4 task this one is derived from: the correct handling is standard practice of the field, it is discoverable from the data with a standard QC step (see the decision points above), and two independent solutions agree with the planted truth.

**Where the truth comes from.** The generator fixes, before any read exists: each gene's fitness
effect and, independently, its drug effect (sensitiser, resistance or none); each sgRNA's
efficiency; the treatment every library really received (one vehicle and one drug library are
exchanged on the sample sheet); how many reads each sgRNA contributes to each library; which
checklist items are defective in the inherited project. `reference/truth.json` stores these planted
values. The inherited files are the output of the predecessor's scripts (seeded mutations of a
correct pipeline) run on the task data; the audit truth is the list of mutations, the erratum
truth is the planted class of every mailed gene. No tool output is used as reference. A numpy port
of the reference statistics is run at generation time only to confirm that the stated acceptance
criteria can be met (arms recoverable, SOP section 5 reachable); it never defines a truth value.
A drug effect is planted only where a sound analysis has power (expected gene z >= 6), so no scored
gene is a borderline case.

**Hard gates (score 0).** Missing/unparseable file; `counts.tsv` not listing every library sgRNA
exactly once with non-negative integers; `gene_summary.tsv` not listing every targeted gene exactly
once; `drug_class` outside {sensitiser, resistance, none}; `audit.tsv` not listing every checklist
item exactly once with `ok`/`defective`; `erratum.tsv` not listing every mailed gene exactly once
with `confirmed`/`retracted`; non-finite numbers.

**Components.**

| weight | component | metric | full credit | zero |
|---|---|---|---|---|
| 0.15 | guide_counts | per library `sum(min(obs,true)) / max(sum obs, sum true)`, then the **minimum** over libraries | >= 0.90 | <= 0.50 |
| 0.20 | gene_lfc | the largest of three median absolute errors against the planted latent values: `lfc_fitness` (all genes), `lfc_drug` (planted non-modifiers), `lfc_drug` (planted modifiers) | <= 0.25 | >= 0.60 |
| 0.30 | drug_class | three-class Matthews correlation with the planted class, over the genes that are **not** reference genes | >= 0.90 | <= 0.75 |
| 0.20 | handover_audit | Matthews correlation of `defective` verdicts with the planted defects (8 items) | 1.0 | <= 0.50 |
| 0.15 | erratum | balanced accuracy (mean of sensitivity and specificity) of the verdicts against "mailed gene is a planted modifier" | >= 0.85 | <= 0.60 |

`qc.json` is required and format-gated (an entry per library); its content is written to the
verifier details (`qc_arms_right`, `qc_read_numbers_right`) but carries no weight: the arms are
scored through `lfc_drug` and `drug_class`, the read numbers through `guide_counts`.

**Why these tolerances.** Counts and the LFC band are those of level 3 (exact-match counting keeps
about 96% of the reads; correct LFC estimates have median errors of 0.08-0.2, the wrong
normalisation, the wrong baseline, the sample-sheet arms or a reversed contrast give 0.7-2.4). The
modifier subset is scored separately because the median over all genes is blind to a reversed or
diluted drug contrast (85% of the genes have a drug effect of 0). For `drug_class` the stated
requirement is a false discovery rate of 5% (one or two false calls among ~26 scored modifiers,
MCC >= 0.96). Full credit starts at MCC 0.90, i.e. about five false calls or four misses, three
times the nominal error; two different abundance-aware tests reach 0.96-1.0 and MAGeCK RRA on the
right contrast 0.89-1.0 (it is a must-pass probe with floor 0.95). Calls against the
spread of the non-targeting sgRNAs give 20-55 false calls (MCC 0.55-0.8) and the level-3 rank test
both over- and under-calls (0.65-0.9); the band is put between the two groups. Reference genes are
excluded from this component because their class is given to the agent. Audit verdicts are facts,
so only a complete audit earns full credit; one wrong item costs about 0.09 of the total. The erratum
uses balanced accuracy because a mailed list can contain as few as four true modifiers (the generator
guarantees at least four of each verdict), where MCC moves by 0.1 per verdict; "all confirmed" and
"all retracted" both score exactly 0.5, i.e. zero credit. Every component is therefore blind to its
trivial answers ("all ok", "all defective", "all retracted", "nothing called").

**What is deliberately not scored.** P-values, the method used, the reference genes' own classes
(they are the agent's control), free-text explanations.

## Verifier attack suite (must score low)

- `l3_reference` (max 0.75): level-3 reference method on the level-4 inputs: correct counts, control normalisation, every contrast against T0 with the sample-sheet arms, sensitiser = level-3 hit in the drug arm that is no hit in the vehicle arm, no audit (all ok), nothing retracted
- `trust_inherited` (max 0.75): ignores H: finishes from the inherited count table with the reference statistics, audit all ok
- `trust_sample_sheet` (max 0.75): ignores P(a): everything as the reference, but the arms are taken from the sample sheet
- `drug_vs_t0` (max 0.75): ignores P(b): everything as the reference, but the drug effect is the drug arm against T0
- `calls_vs_nontargeting_spread` (max 0.75): ignores Q: everything as the reference, but genes are called against the spread of the non-targeting sgRNAs (the level-3 second solution's null)
- `calls_l3_rank_test` (max 0.89): ignores Q, second variant: the level-3 reference test (Mann-Whitney against the non-targeting sgRNAs), two-sided on the right contrast
- `unfixed_library_release` (max 0.75): one defect operator left in: gene annotation of release 1
- `unfixed_sample_columns` (max 0.75): one defect operator left in: FASTQ files taken in lexicographic order
- `unfixed_stagger_reads` (max 0.86): one defect operator left in (level-3 trap): fixed-offset parsing
- `unfixed_read2_library` (max 0.86): one defect operator left in (level-3 trap): the Read-2 library is not rescued, reported as zeros and left out of the statistics
- `unfixed_control_normalisation` (max 0.81): one defect operator left in (level-3 trap): median-of-ratios normalisation over the whole library
- `unfixed_contrast_sign` (max 0.75): one defect operator left in: vehicle against drug
- `naive_default_pipeline` (max 0.5): fully careless: forward strand only (the Read-2 library is lost and left out), sample-sheet arms, median-of-ratios normalisation, calls against the non-targeting spread, no audit (all ok)
- `erratum_nothing_retracted` (max 0.86): correct results, every mailed gene confirmed
- `prior_only` (max 0.05): no analysis: counts all 1, every gene a sensitiser with lfc -3 / -1, audit all defective, everything retracted
- `lfc_drug_shifted_beyond_tolerance` (max 0.81): planted truth, every lfc_drug shifted by +0.7
- `calls_beyond_tolerance` (max 0.75): planted truth with ten false and eight missed modifiers among the scored genes

## Must-pass probes (guard against false negatives)

- `planted_truth` (min 0.999): submission written directly from the planted truth
- `within_tolerance` (min 0.99): 4.5% of reads lost, LFC noise sd 0.08 on both scales, two false and one missed modifier among the scored genes, one wrong erratum verdict
- `independent_solution` (min 0.99): second implementation without cutadapt, MAGeCK or scipy: pure-python counting, arms from the reference-gene signature, calls from the spread of sgRNAs of similar abundance (no replicate variance model), the inherited tables reproduced numerically instead of re-running the inherited scripts
- `mageck_rra_calls` (min 0.95): third route to the calls: the reference pipeline with the calling step replaced by MAGeCK RRA (mageck test, control normalisation) on the drug-versus-vehicle contrast with the arms from the data. Floor 0.95 instead of 0.99: RRA estimates its variance from the vehicle arm and the drug arm is noisier, so its FDR is slightly anti-conservative (0-6 false calls measured, MCC 0.89-1.0)
