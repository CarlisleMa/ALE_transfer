# Bulk RNA-seq: taking over a half-finished differential-expression project (level 6: pointers removed)

Take over an RNA-seq project from a predecessor: re-quantify twelve donor samples, apply the SOP's identity control, deliver gene counts, sample QC and treated-vs-vehicle differential expression, audit the inherited scripts and tables against a checklist and issue an erratum for a gene list that was already mailed. Level 6: same inputs and planted truth as level 4; the text no longer says what to look out for (see the task card).

- family: `rnaseq_confound`  seed: `3`  level: `6`
- ALE subdomain: 3.2 Genomics & Sequence Analysis
- software: salmon, pyDESeq2, statsmodels, Python, numpy/scipy/pandas
- a-priori difficulty: 0.707 (last-exam)  |  pass threshold: 0.9
- capabilities trained: auditing inherited code and intermediate files by recomputation, applying a stated acceptance criterion to an in-data control, experimental-design reasoning (aliased batch), sample identity QC before statistics, long tool chain from raw reads to statistics, exact output contracts

## Undisclosed decision points (not shown to the agent)

- **batch_confound** - Two library-prep sessions, 5 treated + 1 vehicle in one and 1 treated + 5 vehicle in the other; 30% of genes differ between sessions by 0.5-2 log2 units. A condition-only model reports batch genes as treatment genes.
  - discoverable by: samples.tsv has a library_prep_batch column that is almost aliased with condition; lab_notes.md gives two prep sessions; the SOP checklist asks whether the predecessor's model contains the session
  - trains against ALE failure mode: Approach / wrong strategy
- **sample_swap** - The FASTQ files of RS43 and RS49 hold each other's library (one treated, one vehicle, opposite donor sex, different prep sessions).
  - discoverable by: the required inferred_sex column and SOP section 2 make the analyst type every library from XIST and chrY genes; exactly two libraries contradict donor_sex and they are a reciprocal pair; lab_notes.md records that all tubes were re-labelled from a box map
  - trains against ALE failure mode: Understanding / unverified assumption about inputs
- **H_inherited_project** - previous_analysis/ holds the predecessor's scripts, merged matrix, draft table and mailed gene list. Planted defects (5 of 8 checklist items): quant_library_type, matrix_sample_files, matrix_all_reads, de_session_term, de_contrast_direction. The other items are sound; `quant_library_type` uses `-l SR`, the contrast is written in the `plain` style and padj in the `deseq` style. 23 of the 40 mailed genes are not planted responders. Reusing the inherited matrix or draft, or rating the files without recomputing, fails.
  - discoverable by: SOP section 1 forbids unverified reuse and section 3 lists the eight requirements; code-level defects are visible in the three short scripts; file-level defects appear when the scripts are re-run on the FASTQ files and compared with the shipped matrix / draft (identical columns, a column with a third of the reads, a draft that the script does not reproduce)
  - trains against ALE failure mode: Verification / trusting inherited intermediate results
- **Q_identity_control_mixed_library** - The library in RS45.fastq.gz is a mixture: 69% of the recorded male donor's library and 31% of a female library of the other condition and session (RS56). Its Y-gene reads still outnumber XIST, so a male/female call agrees with the sample sheet; its minor marker share is 0.089 against the SOP limit of 0.02 (largest value among the other libraries: 0.0025). It cannot be relabelled and has to be excluded.
  - discoverable by: SOP section 2 states the acceptance criterion on the sex-chromosome control genes of the panel; computing X and Y per library, as the SOP defines them, shows one male library with XIST reads at a third to a half of the female level
  - trains against ALE failure mode: Verification / acceptance criterion on an in-data control not evaluated

## Judge

**Level 6.** Inputs and planted truth are those of `rnaseq-confound-s0003-l4`; only text was taken away:

- S4 + S2 + S1: the take-over SOP is not shipped: its principle about inherited files, the identity control with the marker-share formula and the 0.02 limit, and the sentence in the hand-over note that refers to it. The policy for a mixed library (reported as `mixed`, excluded) stays, as a definition of the deliverable; the marker genes stay in the gene panel
- S3: the eight-item audit checklist and the deliverable audit.tsv; its weight (0.20) is redistributed proportionally over the other seven components
- S5: wet-lab notes without the incident report (tube labels illegible, tubes re-labelled from the box map)

Fairness no longer rests on a stated criterion but on three things that are checked for the level-4 task this one is derived from: the correct handling is standard practice of the field, it is discoverable from the data with a standard QC step (see the decision points above), and two independent solutions agree with the planted truth.

**Where the truth comes from.** Everything scored is planted. Before any read exists the generator
fixes which genes respond and by how much, the library-prep batch effects, sex and treatment of
every library, which two FASTQ files hold each other's library, which library is a mixture of two
donors, and which checklist items of the SOP are defective in the predecessor's files.
`reference/truth.json` stores exactly these planted values. The predecessor's files are produced by
running mutated copies of the reference pipeline on the task's own reads; they are inputs, never
the reference. The erratum truth is the planted responder status of the mailed genes.

**Hard gates (score 0).** As at level 3 (missing/unparseable file, gene or sample missing or
duplicated, negative counts, invalid `condition_used`, `significant` inconsistent with `padj`,
non-finite `log2fc`), plus: `audit.tsv` must give `ok`/`defective` for every checklist item exactly
once; `erratum.tsv` must give `confirmed`/`retracted` for every mailed gene exactly once.

**Components.**

| weight | component | metric | full credit | zero |
|---|---|---|---|---|
| 0.08 | gene_counts | per sample `sum(min(obs,true)) / max(sum obs, sum true)`, minimum over samples; 0 when a column agrees better with the planted counts of another FASTQ file than with those of its own | >= 0.90 | <= 0.50 |
| 0.02 | sample_sex | `inferred_sex` is the planted sex for every library, `mixed` for the mixed one; only if gene_counts > 0 | all | any wrong |
| 0.07 | sample_labels | mixed library `excluded`; every other library under its planted condition (the two mixed-up ones may be `excluded`); nobody else dropped; only if gene_counts > 0 | all | any wrong |
| 0.19 | de_false_discoveries | as level 3 | <= 0.12 | >= 0.30 |
| 0.16 | de_sensitivity | as level 3 | >= 0.85 | <= 0.50 |
| 0.13 | de_log2fc | as level 3 | <= 0.33 | >= 0.55 |
| 0.20 | audit | Matthews correlation of the verdicts against the planted defects | 1.0 (complete audit) | <= 0.70 |
| 0.15 | erratum | Matthews correlation of `confirmed`/`retracted` against the planted responder status of the mailed genes | >= 0.80 | <= 0.40 |

The three DE components are multiplied by 0.5 when the sample table is not valid (SOP section 2: a
table from a model that contains the mixed library, or a library under the wrong condition, is not
released). One wrong audit verdict gives MCC 0.745-0.775 on eight items with three to five defects,
i.e. 0.15-0.25 of the audit credit and a task score of at most 0.85 with everything else perfect: a
single missed defect or false alarm ends below the pass mark. Two wrong verdicts give no audit credit.

**Why these tolerances.** See the "Level 4" section of `docs/families/rnaseq_confound.md` for the
measured values of the oracle, the independent solution, the level-3 reference method and every
careless variant on two seeds.

**What is deliberately not scored.** The numeric value of `padj`, `log2fc` of genes without a planted
effect, the DE status of the sex-chromosome control genes, the evidence column of the audit, which
of the two valid resolutions of the file mix-up was chosen.

## Verifier attack suite (must score low)

- `l3_reference` (max 0.75): the level-3 reference method on the level-4 inputs: fresh salmon quantification, male/female call, mix-up resolved, session term; no SOP identity criterion; the predecessor's work is taken as fine (audit all ok); erratum from its own DE table
- `trust_inherited` (max 0.8): H ignored: the analysis is finished from the inherited count matrix (identity control, mix-up and session term handled correctly), nothing is re-checked, audit all ok
- `q_identity_not_applied` (max 0.75): Q ignored: everything right (fresh quantification, complete audit, erratum) but libraries are only typed male/female, so the mixed library stays in the model
- `sheet_copied_no_analysis` (max 0.1): prior-only answer: QC table copied from the sample sheet, counts all 1, no gene significant, audit all ok, every mailed gene confirmed
- `everything_significant` (max 0.45): true counts, QC and audit, but every gene called significant with log2fc 1 and every mailed gene confirmed
- `erratum_all_confirmed` (max 0.89): everything perfect, but the mailed list is confirmed wholesale
- `fdp_beyond_tolerance` (max 0.85): everything perfect, but 40% of the significant genes are false
- `lfc_attenuated_beyond_tolerance` (max 0.89): everything perfect, but every log2fc shrunk by a third
- `naive_labels_as_given` (max 0.75): level-3 trap ignored: SOP identity criterion applied, fresh quantification, complete audit, but the sample sheet is used as delivered
- `defect_kept_quant_library_type` (max 0.8): planted defect `quant_library_type` neither reported nor repaired: the final analysis repeats it, everything else is right
- `defect_kept_matrix_sample_files` (max 0.8): planted defect `matrix_sample_files` neither reported nor repaired: the final analysis repeats it, everything else is right
- `defect_kept_matrix_all_reads` (max 0.8): planted defect `matrix_all_reads` neither reported nor repaired: the final analysis repeats it, everything else is right
- `defect_kept_de_session_term` (max 0.8): planted defect `de_session_term` neither reported nor repaired: the final analysis repeats it, everything else is right
- `defect_kept_de_contrast_direction` (max 0.8): planted defect `de_contrast_direction` neither reported nor repaired: the final analysis repeats it, everything else is right

## Must-pass probes (guard against false negatives)

- `planted_truth` (min 0.999): submission written directly from the planted truth
- `within_tolerance` (min 0.99): 6% of reads lost, 4 false discoveries, 8% of planted genes missed, log2fc noise sd 0.2, two erratum verdicts wrong (one each way); audit exact
- `independent_solution` (min 0.99): second tool chain: k-mer read counter (no salmon), voom-style weighted linear model (no DESeq2), mismatching and mixed libraries excluded, audit by profile comparison and static inspection instead of a salmon re-run
- `oracle_exclude_resolution` (min 0.99): oracle tool chain, but the two mixed-up samples are excluded instead of relabelled
