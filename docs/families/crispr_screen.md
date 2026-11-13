# Family `crispr_screen` - pooled CRISPR dropout screen, FASTQ to fitness genes

Reference implementation of the pipeline (written by the maintainer).

## Scenario
A core facility analyst gets single-end FASTQ files of a pooled CRISPR-Cas9 knockout screen run with
a secondary validation library (T0 x2, T14 x3), a library file, a sample sheet and wet-lab notes,
and must deliver guide counts, per-gene log2 fold changes and fitness-gene calls.

## Truth model
Everything is planted before a read exists: which genes have a fitness effect and its size
(log2FC -1.5 to -4), each sgRNA's efficiency (at most one inactive guide per gene), each guide's
abundance per sample, and therefore the exact number of reads per guide per sample. Reads are then
rendered from the amplicon structure (U6 3' end, 20-nt spacer, scaffold) with 0.2% substitution
errors and 3-6% unrelated reads. Gene identities are blinded (`TGT0421`), so knowledge of real
essential genes cannot replace the analysis. `truth.json` = planted counts, latent median-guide
LFC per gene, fitness flag.

## Decision points and why they are discoverable
| trap | effect if ignored | evidence available to the agent |
|---|---|---|
| `staggered_amplicon` | fixed-offset parsing assigns ~1/8 of reads | notes name the "P5 stagger primer mix"; mapping rate |
| `read2_only_sample` | one sample maps ~0% | notes say only Read 2 was delivered; read names end in ` 2:N:0`; mapping rate |
| `hit_enriched_library` | median normalisation shifts every LFC by ~+1.9; rank-based RRA finds almost no hits | notes describe a validation library of carried-forward hits; non-targeting guides get a large positive LFC under default normalisation |

## Output contract
`counts.tsv`, `gene_summary.tsv` (`gene, lfc, is_hit`), `qc.json`; definitions in the instruction.
The instruction defines the LFC scale ("a knockout without fitness effect has lfc 0") and the
meaning of a hit ("depleted compared with sgRNAs that have no fitness effect"); it does not say how
to normalise or which test to use.

## Scoring and measured values
Raw metrics are from seed 1, level 3 (`results/work/crispr-screen-s0001-l3/*/logs/details.json`);
the last column is the score on seed 2, level 3 (`results/validation/crispr-screen-s0002-l3.json`).

| submission | guide_counts (min over samples) | LFC median abs error | hit MCC | score (seed 1) | score (seed 2) |
|---|---|---|---|---|---|
| oracle (cutadapt + mageck count + control-based test) | 0.9553 | 0.0829 | 0.9848 | 1.0 | 1.0 |
| independent solution (pure python, bootstrap null) | 0.9747 | 0.0755 | 0.9848 | 1.0 | 1.0 |
| truth with 4.5% read loss, LFC noise 0.08, two false hits | 0.9544 | 0.0567 | 0.9696 | 1.0 | 1.0 |
| fixed-offset parser | 0.0 | 0.2144 | 0.9697 | 0.589 | 0.60 |
| failed sample dropped instead of rescued | 0.0 | 0.0928 | 0.9543 | 0.64 | 0.64 |
| correct counts, median normalisation + RRA | 0.9553 | 1.9269 | 0.0 | 0.40 | 0.40 |
| MAGeCK end to end with defaults | 0.0 | 1.9528 | 0.0527 | 0.04 | 0.04 |
| every gene a hit, lfc -3, counts all 1 | 0.0053 | 0.7318 | 0.0 | 0.0 | 0.0 |
| perfect counts and calls, every lfc shifted by +0.7 | 1.0 | 0.70 | 1.0 | 0.70 | 0.70 |

Bands: counts full >= 0.90 and zero <= 0.50 (exact-match counting keeps 95.5%, a 1-mismatch counter
97.5%); LFC full <= 0.20 and zero >= 0.60 (correct analyses 0.06-0.08, wrong normalisation 1.93);
MCC full >= 0.90 and zero <= 0.50 (correct 0.98, defaults 0.0). Pass threshold 0.9.
One number sits close to a band edge and is worth knowing: the fixed-offset parser, which uses
one eighth of the reads, still estimates the LFC with error 0.21, i.e. the LFC and hit components
are robust to an 8x loss of depth. That is why the counts component uses the minimum over samples.

## Levels
1: no stagger, all forward, library mostly neutral (default tools work). 2: stagger + Read-2
sample. 3: + hit-enriched library. All three validated (`results/validation/crispr-screen-*`).
4: the level-3 task inside a hand-over of a drug-modifier screen (H + P + Q); see "Level 4" at the
end of this file.

## Small-model pilot (not a Harbor run)
Two single attempts by Claude Haiku through `docker exec` (`alegen/pilot.py`), seed 1, level 3:

| run | isolation | reward | wall time | count agreement | LFC error | hit MCC |
|---|---|---|---|---|---|---|
| 16:05 | by instruction only (agent had host access); unaudited | 1.0 | 2 min 49 s | 0.999 | 0.0845 | 0.9848 |
| 19:17 | answers sealed, access-time audit clean (`results/pilot/AUDIT_batchA.txt`) | 1.0 | 2 min 1 s | 0.9987 | 0.0807 | 0.9848 |

Reading: this task is **easier than its prior (0.643) says**. The lab notes state all three decision
points plainly, the instruction's definitions of `lfc` and `is_hit` point at the right normalisation
and test, and a small model does the per-sample QC that reveals the first two anyway.

## Known limitations
- Difficulty is overstated by the prior; see the pilot above. The decision points are disclosed
  clearly enough that a small model handled all three.
- Reads are rendered from a clean amplicon model: no index hopping, no PCR chimeras, constant
  qualities. Real data are messier; the counting step here is easier than in production.
- The fitness effect is bimodal (clear hits vs neutral). Borderline genes, where methods
  legitimately disagree, are deliberately absent, so the task does not train threshold judgement.
- "Hit" is defined relative to non-targeting controls; an analyst who defines hits relative to the
  library median is graded wrong even though some labs do that for genome-wide screens. The
  instruction states the definition, so this is a specification, not a hidden requirement.

# Level 4 (2026-11-12): taking over a drug-modifier screen (H + P + Q)

Code: `alegen/families/crispr_screen/l4.py` (generator, own random stream; levels 1-3 are untouched
and regenerate byte-identically), `assets/score_l4.py`, `assets/crispr_l4_pipeline.py` (reference
solution; `--mut` / `--l3` replay the careless variants), `assets/solve_l4_alt.py` (independent
second solution), `assets/solve_l4.sh`. Spec: `docs/specs/l4/crispr_screen.md`.

## What was added, and which pilot finding each element answers

The pilot (README section 6) solved level 3 in 2 min 1 s with 17 tool calls: the lab notes named all
three decision points and a per-sample mapping rate revealed two of them. Level 4 keeps the three
level-3 decision points and adds three elements that are not of that kind.

| element | what the task contains | pilot finding it answers | where the acceptance criterion is stated |
|---|---|---|---|
| **H** inherited project | `/task/input/handover/`: the predecessor's count script, count table, analysis script, gene table, mailed hit list and a hand-over note. An 8-item checklist (SOP section 6); 3-5 items are defective per seed, the others are fine. Deliverables `audit.tsv` and `erratum.tsv`. | "a plausible shortcut that the generator had closed" (somatic-tn): finishing from the inherited count table is the shortcut, and every seed carries at least one count-table defect that leaves no trace in the scripts. Also lengthens the chain: the reference solution has 19 tool invocations (level 3: 9). | SOP section 6 (one statement per item), SOP section 7 and the instruction (erratum verdicts) |
| **P(a)** sample sheet != data | One vehicle and one drug library are exchanged on the sample sheet (which pair varies with the seed). With the sheet arms the drug contrast is diluted to one third. | "every decision point that lab notes or a routine QC step reveal was handled": here the authoritative file is wrong and has to be overruled from the data. The notes only say that the tubes were re-racked. | SOP section 3 (a library must agree better with its own arm; it is analysed in the arm its profile supports; the arm used goes into `qc.json`), instruction ("grouped by the treatment they actually received") |
| **P(b)** dropout contrast | Fitness and drug effects are drawn independently in a library where 83% of the genes are fitness genes, so "depleted in the drug arm compared with T0" (the contrast of every level-3 run and of the standard dropout workflow) returns the fitness genes. In some seeds the predecessor's script uses exactly this contrast. | same finding; the level-3 habit is the prior that fails | instruction and SOP section 4 define the drug effect against vehicle |
| **Q** quality bar with an in-data control | Drug-arm libraries come from about 25 cells per sgRNA (vehicle 150), so the noise of an sgRNA grows steeply with its depletion and the non-targeting sgRNAs (never depleted) no longer represent the null. 40 reference genes of known behaviour (10 modifiers, 30 `none` genes over the whole fitness range) are given. Modifiers include a moderate tier (planted \|log2FC\| 0.7-1.4) next to a strong one (1.5-3.0). | "a final quantity that had to be good enough and needed self-checking and method tuning" (amplicon-16s): every step can be done and the call set still be wrong | SOP section 5: at least 9 of 10 reference modifiers with the expected class, at most 1 of 30 reference `none` genes called; instruction: FDR 5% |

## Truth model (additions to level 3)
Planted before any read exists (`l4.plan`): fitness effect per gene (190 of 230, log2FC -1.5 to -4);
independently a drug effect for 18 sensitisers and 18 resistance genes; sgRNA efficiencies (at most
one inactive sgRNA per gene, 45% of genes); the treatment each library really received; cells per
sgRNA at harvest (T0 500, vehicle 150, drug 25: a Poisson bottleneck before the read multinomial);
read counts per sgRNA and library; the 20 genes whose labels are shifted in release 1 of the
annotation (10 scored modifiers alternating with 10 scored non-modifiers); the defective checklist
items. sgRNA ids are neutral (`VAL2_0001`), so an id does not name its gene.

A drug effect is planted only where a sound analysis has power: expected gene-level z >= 6 from the
planted read and cell numbers, otherwise the effect is redrawn from the strong tier or the gene is
skipped. There is therefore no gray zone and no scored borderline gene.

`l4.solvable` runs numpy ports of the reference statistics on the planted counts and re-draws the
instance (deterministically, `sub` = 0, 1, ...) unless the arms are recoverable from the first
principal component, the release criterion of SOP section 5 is met and the scored MCC is >= 0.93.
This is a solvability filter; no truth value comes from it. Of seeds 1-40, 39 pass at the first draw
(seed 29 is redrawn: the port calls two reference `none` genes).

## The inherited project
The predecessor's scripts are text templates of a correct stdlib/numpy pipeline; each seeded defect
is one mutation. The inherited tables are produced at generation time by **running those scripts on
the task data** (`sys.executable`, host venv), so every inherited file is exactly what the script in
the folder yields when it is invoked the way the defect implies.

| item | mutation operator | visible in the script? |
|---|---|---|
| `library_release` | count script called with release 1 of the annotation | no (the path is an argument) |
| `sample_columns` | FASTQ files passed in lexicographic (shell glob) order while the script takes them in sample-sheet order; library ids cross from one to two digits (`LIB8` ... `LIB15`) | no (the docstring states the expected order) |
| `stagger_reads` | spacer read at fixed offsets | yes |
| `read2_library` | forward strand only | yes |
| `arm_labels` | arms taken from the sample sheet | no (the sheet is wrong, not the script); defective in every seed |
| `control_normalisation` | median-of-ratios over the whole library, no centring on the non-targeting sgRNAs | yes |
| `drug_baseline` | drug arm against T0 | yes |
| `contrast_sign` | vehicle against drug | yes |

Seed 1: `library_release`, `read2_library`, `arm_labels`, `control_normalisation`. Seed 2:
`library_release`, `read2_library`, `arm_labels`, `drug_baseline`, `contrast_sign`. Over seeds 1-40
every operator is active in 11-25 seeds and a seed has 3, 4 or 5 defects (12 / 13 / 15 seeds).

Per seed: `arm_labels` plus 2-4 of the other seven, at least one of the two that leave no trace in
the script. Together with `arm_labels` every seed has at least two defects that can only be found
by recomputing. Fine items that look unusual: exact-match counting (about 4% fewer reads than a
mismatch-tolerant counter), size factors from 80 sgRNAs only, an analysis that silently leaves out a
library without reads.

The mailed list is the top 40 genes by |z| of the inherited gene table (gene names only); if fewer
than four of them are planted modifiers, or fewer than four are not, the predecessor's `TOP` is
raised to 60, 80, 120.

## Output contract
`counts.tsv`, `gene_summary.tsv` (`gene`, `lfc_fitness`, `lfc_drug`, `drug_class`), `audit.tsv`
(`item`, `verdict`), `erratum.tsv` (`gene`, `verdict`), `qc.json` (reads, mapped reads, arm used).
The instruction defines both LFC scales, the class definition (FDR 5%, release criterion of the SOP),
which FASTQ a count column refers to, and the erratum verdict as a function of the final classes.
It does not say how to find the arms, how to normalise or how to test.

## Scoring and measured values
Weights are those of the spec: `guide_counts` 0.15, `gene_lfc` 0.20, `drug_class` 0.30,
`handover_audit` 0.20, `erratum` 0.15; pass 0.9. Raw metrics below are from
`results/work/crispr-screen-s000{1,2}-l4/*/logs/details.json`, scores from
`results/validation/crispr-screen-s000{1,2}-l4.json` (2026-11-12, all gates green for both seeds;
seed 1 including Harbor oracle = 1.0 / nop = 0.0; gate logs in `results/crispr_s1_l4.log` and
`results/crispr_s2_l4.log`). Agent image 1.18 GB; inputs 21 MB; generation 13 s; oracle runtime in
the agent image 0.39-1.21 min (seed 1, idle and loaded host), 0.55-0.58 min (seed 2); validation
13.8 min (seed 1 with Harbor) and 7.9 min (seed 2).

Columns: counts = minimum over libraries of the count agreement; LFC = median absolute errors of
`lfc_fitness` / `lfc_drug` of planted non-modifiers / `lfc_drug` of planted modifiers; class =
three-class MCC on the 190 non-reference genes (false calls, missed or wrong direction; 26 scored
modifiers); ref = reference modifiers right of 10 / reference `none` genes called of 30 (the SOP
criterion is >= 9 and <= 1); audit = MCC over the 8 items; erratum = balanced accuracy.

| submission | counts | LFC errors | class MCC | ref | audit | erratum | score s1 | score s2 |
|---|---|---|---|---|---|---|---|---|
| **oracle** (cutadapt + mageck count, arms from PC1, abundance-aware z test) s1 | 0.957 | 0.092 / 0.121 / 0.150 | 1.0 (0, 0) | 10 / 0 | 1.0 | 1.0 | **1.0** | |
| oracle s2 | 0.956 | 0.095 / 0.111 / 0.117 | 1.0 (0, 0) | 10 / 0 | 1.0 | 1.0 | | **1.0** |
| **independent solution** (pure python, arms from reference signature, spread of similar-abundance sgRNAs, tables reproduced numerically) s1 | 0.975 | 0.088 / 0.121 / 0.147 | 1.0 (0, 0) | 10 / 0 | 1.0 | 1.0 | **1.0** | |
| independent solution s2 | 0.975 | 0.106 / 0.107 / 0.133 | 0.979 (1, 0) | 10 / 0 | 1.0 | 0.986 | | **1.0** |
| **MAGeCK RRA calls** (probe `mageck_rra_calls`, floor 0.95: the reference pipeline with `mageck test`, control normalisation, on drug vs vehicle) | 0.957 / 0.956 | as oracle | 0.923 (4, 0) / 1.0 (0, 0) | 10 / 0, 10 / 1 | 1.0 | 1.0 | **1.0** | **1.0** |
| planted truth | 1.0 | 0 / 0 / 0 | 1.0 | 10 / 0 | 1.0 | 1.0 | 1.0 | 1.0 |
| within tolerance (4.5% reads lost, LFC noise 0.08, 2 false + 1 missed, one wrong erratum verdict) | 0.954 | 0.05 / 0.055 / 0.06 | 0.937 (2, 1) | 10 / 0 | 1.0 | 0.964 / 0.986 | 1.0 | 1.0 |
| `l3_reference` (level-3 method: right counts, control normalisation, both arms against T0 with the sheet arms, no audit) | 0.957 | 0.13 / 2.15 / 1.41 (s1) | 0.04 (1, 26) / 0.0 (0, 26) | 0 / 0, 1 / 0 | 0.0 | 0.5 | 0.15 | 0.15 |
| `trust_inherited` (ignores H: reference statistics on the inherited count table, audit all ok) | 0.0 | 0.13 / 0.15 / 0.21 | 0.571 (10, 10) | 10 / 0 | 0.0 | 0.75 / 0.70 | 0.29 | 0.26 |
| `trust_sample_sheet` (ignores P(a)) | 0.957 | 0.13 / 0.13 / 1.30 (s1), 1.40 (s2) | 0.499 (0, 19) / 0.264 (0, 24) | 3 / 0 | 0.775 | 0.654 / 0.5 | 0.292 | 0.260 |
| `drug_vs_t0` (ignores P(b)) | 0.957 | 0.09 / 2.17 / 1.93 (s1) | 0.324 (0, 23) / 0.082 (18, 22) | 2 / 0, 2 / 3 | 1.0 / 0.775 | 0.596 / 0.6 | 0.35 | 0.260 |
| `calls_vs_nontargeting_spread` (ignores Q: level-3 second solution's null) | 0.957 | as oracle | 0.674 (25, 0) / 0.527 (48, 0) | 10 / 8, 10 / 10 | 1.0 | 0.786 / 0.757 | 0.661 | 0.644 |
| `calls_l3_rank_test` (ignores Q: level-3 reference test, two-sided) | 0.957 | as oracle | 0.811 (5, 4) / 0.827 (7, 2) | 10 / 2, 10 / 4 | 1.0 | 0.926 / 0.971 | 0.821 | 0.853 |
| `unfixed_library_release` | 0.957 | 0.10 / 0.13 / 0.19 | 0.571 (10, 10) | 10 / 0 | 0.775 | 0.75 / 0.70 | 0.550 | 0.520 |
| `unfixed_sample_columns` | 0.557 / 0.583 | 2.10 / 2.16 / 1.46 (s1) | 0.08 / 0.09 | 0 / 0, 0 / 3 | 1.0 | 0.558 / 0.6 | 0.222 | 0.231 |
| `unfixed_stagger_reads` (level-3 trap) | 0.118 / 0.119 | 0.21 / 0.16 / 0.23 (s1), 0.26 / 0.18 / 0.25 (s2) | 0.979 (0, 1) / 0.979 (1, 0) | 10 / 0 | 1.0 | 1.0 | 0.85 | 0.842 |
| `unfixed_read2_library` (level-3 trap) | 0.0 | 0.11 / 0.13 / 0.12 | 1.0 | 10 / 0 | 0.775 | 1.0 | 0.760 | 0.760 |
| `unfixed_control_normalisation` (level-3 trap) | 0.957 | 1.76 (s1), 1.87 (s2) / 0.12 / 0.16 | 1.0 | 10 / 0 | 0.775 (s1), 1.0 (s2) | 1.0 | 0.710 | 0.80 |
| `unfixed_contrast_sign` | 0.957 | 0.09 / 0.12 / 3.56 (s1), 3.78 (s2) | 0.443 (0, 26) | 0 / 0 | 1.0 / 0.775 | 1.0 | 0.50 | 0.410 |
| `naive_default_pipeline` (forward strand only, sheet arms, median normalisation, calls against the non-targeting spread, no audit) | 0.0 | 1.74 / 0.16 / 0.58 (s1) | 0.604 (31, 1) / 0.114 (19, 21) | 10 / 6 (s1) | 0.0 | 0.714 (s1) | 0.069 | 0.0 |
| `audit_all_defective`, `audit_all_ok` (planted truth otherwise) | 1.0 | 0 | 1.0 | | 0.0 | 1.0 | 0.80 | 0.80 |
| `erratum_nothing_retracted` (planted truth otherwise) | 1.0 | 0 | 1.0 | | 1.0 | 0.5 | 0.85 | 0.85 |
| `prior_only` (counts 1, every gene a sensitiser, all defective, all retracted) | 0.008 | 0.87 / 1.0 / 1.98 | 0.0 | | 0.0 | 0.5 | 0.0 | 0.0 |
| `lfc_drug_shifted_beyond_tolerance` (+0.7) | 1.0 | 0 / 0.7 / 0.7 | 1.0 | | 1.0 | 1.0 | 0.80 | 0.80 |
| `calls_beyond_tolerance` (10 false + 8 missed) | 1.0 | 0 | 0.630 | | 1.0 | 1.0 | 0.70 | 0.70 |

In an `unfixed_*` row the analyst repeats exactly one defect operator in the final analysis and
reports the corresponding checklist item `ok`; the audit column is 0.775 (one miss) when the
inherited project has that defect in this seed and 1.0 when it does not.

Measured outside the gates (scratch runs in the validated seed-1 image, same scorer; not kept in
`results/`):

| submission | class MCC | ref | score |
|---|---|---|---|
| reference pipeline with MAGeCK RRA and MAGeCK's default median normalisation, seeds 1 / 2 / 3 | 0.979 (1 missed) / 1.0 / 1.0 | 10 / 0, 10 / 0, 9 / 0 | 1.0 / 1.0 / 1.0 |
| MAGeCK RRA (control normalisation) on the second solution's counts (1-mismatch counting), seed 1 | 0.890 (6 false) | 10 / 1 | 0.988 |
| MAGeCK RRA (control normalisation) on scratch seed 3 | 1.0 | 10 / 0 | 1.0 |
| oracle and second solution on scratch seed 3 (defects `library_release`, `sample_columns`, `stagger_reads`, `arm_labels`, `contrast_sign`; mailed list of 80 with 4 true modifiers) | 1.0 / 0.960 (2 false) | 10 / 0 | 1.0 / 1.0 |
| numpy ports of both callers on the planted counts of seeds 1-20 (true arms) | 0.96-1.0 | >= 9 / <= 1 in all 20 | not scored |

Bands and why:
- `guide_counts`: full >= 0.90, zero <= 0.50, as level 3 (exact-match counting 0.957, 1-mismatch
  0.975; one eighth of the reads 0.12; a lost library 0.0).
- `gene_lfc`: the largest of the three errors, full <= 0.25, zero >= 0.60. Correct analyses
  0.12-0.15; the right method on one eighth of the reads 0.23-0.26 (the upper value just leaves the
  full-credit band: 0.96 credit); single wrong decisions (sheet arms, T0 baseline, reversed sign,
  library-wide normalisation, permuted files) 1.30-3.78; the lowest wrong value measured is 0.58
  (fully careless pipeline, seed 1), still 0.06 credit. The modifier subset is separate because the
  median over all genes cannot see a reversed or diluted contrast (85% of the genes are 0).
  `lfc_fitness` is included so that the level-3 normalisation decision is still scored.
- `drug_class`: full >= 0.90, zero <= 0.75. The stated requirement is FDR 5%; 0.90 is about five
  false calls or four misses. Sound methods 0.92-1.0 on the validated seeds (reference, second
  solution, MAGeCK RRA), careless ones 0.53-0.83. The margin on the careless side is thin for the
  level-3 rank test (0.81 / 0.83, 0.07 below full credit).
- `handover_audit`: full only at MCC 1.0, zero <= 0.5. Verdicts are facts and both solutions get all
  eight; one wrong item costs 0.09 (score 0.91 if nothing else is wrong).
- `erratum`: balanced accuracy, full >= 0.85, zero <= 0.60 (trivial answers are exactly 0.5).

## Deviations from the spec, and what was dropped
1. **Sample exchange moved from the predecessor to the sample sheet (P(a)).** The spec lists "two
   files exchanged between arms" as an H operator. With the instruction defining `lfc_drug` against
   vehicle, the spec's own P ("drug vs T0") is a definition the agent just follows, so it has no
   bite of its own; the exchange on the sheet is the P element that does. The SOP therefore has to
   state the sample-identity criterion (section 3). This is the element I am least sure is judged
   fair by everyone; see the doubts below.
2. **Dropped: bottlenecked replicate with an SOP replicate criterion.** Not built. With one library
   exchanged per arm and one library only as Read 2, a third special library would leave arms of two
   and interact with the arm detection; I had no time to validate that interaction. The heteroscedastic
   noise that Q needs comes from the drug arm's cell number instead (stated in the notes).
3. **Dropped: stale draft table.** A stale gene table changes nothing in the final deliverables, so
   its "left unfixed" variant can only lose audit credit (0.91, a pass). It would also force the
   analysis-stage items to be defined on the table instead of the script, which makes them
   forensic and ambiguous when several defects overlap. The two "only by recomputing" defects are
   `arm_labels` plus `library_release` or `sample_columns` instead.
4. **Weights as in the spec; two definitions changed.** `gene_lfc` (0.20) covers `lfc_fitness` as
   well as `lfc_drug`. The erratum uses balanced accuracy instead of MCC: on scratch seed 3 the
   mailed list of 80 held 4 planted modifiers, the second solution's two false calls (inside the FDR
   tolerance) were both on it, MCC was 0.81 and a correct solution scored 0.981. The erratum verdict
   is defined on "modifier of either class" and the mailed list carries gene names only, because
   with a reversed contrast every class-matched verdict would be `retracted` and the component
   degenerate. `qc.json` is required and format-gated but has no weight (the spec gives it none).
5. **Inherited files come from a stdlib/numpy re-implementation**, not from mutated copies of
   `crispr_l4_pipeline.py`: the generator runs on the host venv, which has neither cutadapt, MAGeCK
   nor scipy. The scripts in `handover/` are the ones that were run, with the arguments the defect
   implies, so the property the spec asks for (each inherited file is exactly what the defective
   code yields) holds.
6. **Ceilings above `pass - 0.15` (0.75)**, all below pass. Measured (seed 1 / seed 2): the three
   level-3 traps left unfixed one at a time, `unfixed_stagger_reads` 0.85 / 0.842 (ceiling 0.86),
   `unfixed_read2_library` 0.760 / 0.760 (ceiling 0.86; both seeds also lose audit credit because
   the inherited project has the same defect, the ceiling leaves room for a seed where it does not),
   `unfixed_control_normalisation` 0.710 / 0.80 (ceiling 0.81); constant audits 0.80 (0.81);
   constant erratum 0.85 (0.86); the LFC shift 0.80 (0.81); the level-3 rank test as second Q
   variant 0.821 / 0.853 (0.89). This follows from the spec weights: counting is 0.15 and
   `gene_lfc` 0.20, so a single level-3 mistake cannot cost 0.25. `l3_reference` and the
   element-ignoring attacks for H, P(a), P(b) and Q score 0.15-0.66.
7. **`oracle_steps` = 19** (8 cutadapt, 8 mageck count, the two inherited scripts, the pipeline),
   one short of the 20 the contract asks for. Not padded.
8. **Battery additions beyond the spec.** A must-pass probe `mageck_rra_calls` (floor 0.95 with the
   reason written in the probe: RRA takes its variance from the vehicle arm and is slightly
   anti-conservative here) and the fully careless pipeline `naive_default_pipeline` (ceiling 0.5).
   Per-operator `unfixed_*` attacks are run for all operators on every seed, not only for the
   defects that are active in the inherited project of that seed.

## Known limitations and doubts
- **Q does not bite an agent that runs `mageck test` on the right contrast.** MAGeCK RRA with the
  right arms scores full on seeds 1-3 (MCC 0.92-1.0; 0.89 on the 1-mismatch counts). Q stops the
  level-3 statistics (null from the non-targeting sgRNAs: 25-48 false calls; level-3 rank test: 5-7
  false and 2-4 missed) and ad hoc cut-offs, i.e. an agent that carries its level-3 method over.
  How often a small model does that is unknown until the pilot. MAGeCK was the wrong tool at
  level 3 (RRA collapses in the hit-enriched fitness contrast) and is a right one here.
- **The reference criterion is necessary, not sufficient.** "At most 1 of 30 `none` genes called"
  admits a false-call rate of about 3%, which on 164 scored non-modifiers is five false calls, the
  edge of the band (seed 1, MAGeCK on 1-mismatch counts: criterion met, MCC 0.89, score 0.988).
  A stricter criterion (0 of 30) would be failed by sound methods by chance.
- **P(a) is stated in the SOP.** Fairness requires the criterion to be written down, so a model
  that works through the SOP will run the check. The arms separate (seed 1, planted counts,
  Pearson of log2 normalised counts: vehicle-vehicle 0.956-0.958, drug-drug 0.847-0.885, between
  arms 0.813-0.827; after centring each sgRNA over the six day-14 libraries: vehicle-vehicle +0.58
  to +0.60, drug-drug -0.03 to +0.17, between arms -0.48 to -0.56). The drug libraries are the
  noisy ones, so on raw correlations their margin is small (0.85 against 0.83). Whether the model
  then dares to overrule the sheet is the open question. `arm_labels` is defective in every
  seed (only the pair varies), so it is not a seeded item.
- **Five of the eight audit items can be decided without recounting** (three by reading the
  analysis script, two from column totals). The hand-over note lists the reads assigned per library,
  so an empty column is visible at a glance.
- **Single level-3 traps cost less than at level 3** (see deviation 6).
- **A mailed list can hold few true modifiers** (seed 2: 5 of 40; scratch seed 3: 4 of 80). The
  generator guarantees four of each verdict; with four, one missed confirmation is still full credit
  (0.875), two are not.
- **The two validated seeds carry six of the eight operators as inherited defects.**
  `sample_columns` and `stagger_reads` are inherited defects only on scratch seed 3, where the
  oracle and the second solution were run and scored (both 1.0, audit complete) but the gates were
  not. Every operator is exercised as a careless variant on both validated seeds.
- Noise model: Poisson cell bottleneck plus read sampling; no jackpot clones, no index hopping,
  constant qualities (as level 3). No borderline modifiers by construction.
- No agent has run this task. The instruction and SOP were only tested for ambiguity by my own two
  solutions.
