# Family `xlms_model_selection` - choosing among predicted complex models with XL-MS restraints

ALE subdomain 3.1 (Biomolecular Structure & Design). Code: `pipeline/alegen/families/xlms_model_selection/`.
Task ids: `xlms-model-selection-s<seed>-l<level>`.

## Scenario

A structural-biology group has 22-28 predicted models of a two-protein complex and a table of
lysine-lysine cross-links measured on the complex in solution. The agent has to map every
cross-linked position onto the models, count per model how many inter- and intra-protein
cross-links are satisfied under the group's distance criterion, rank the models by how well the
data support the arrangement of the two proteins, and name the best model.

Inputs (`/task/input`): `models/model_NN.pdb`, `crosslinks.csv`, `sequences.fasta`,
`xlms_notes.md`, and at level 3 `ranking.json`. Tools: Python 3.12, numpy, scipy, pandas,
Biopython 1.85 (one conda env, nothing else). Input size 7.5-13.3 MB; generation 5-10 s on the host
(numpy only, no tool image); oracle runtime 3 s.

## Truth model (what is planted and why it is correct)

1. **Hidden native.** One of five RCSB heterodimers, chosen by `(seed - 1) % 5`. The file is read
   once on the build host and never shipped; `reference/truth.json` (verifier image only) is the
   only place that names it.

   | entry | chains used | URL | sha256 of the cached file |
   |---|---|---|---|
   | 1FIN | A, B | https://files.rcsb.org/download/1FIN.pdb | `12f7920d265a08c0145152e6f960d942a17dd947585fab0470869111e8136fca` |
   | 1EFU | A, B | https://files.rcsb.org/download/1EFU.pdb | `f621d20b05abfc7550a5250a971023cbfa11a49fefaf53946153674d11b687c0` |
   | 1KXP | A, D | https://files.rcsb.org/download/1KXP.pdb | `5d68ab0d9bec2230c7240a77988cd00a4cadff3dd97d814c60a85e6182d76418` |
   | 1IBR | A, B | https://files.rcsb.org/download/1IBR.pdb | `181e03e9c5b5122ed91b91600927f017286a0cabb43575cb0a36550aa6a65c14` |
   | 1I2M | A, B | https://files.rcsb.org/download/1I2M.pdb | `d11d60e277912fc62f2379ae4d2a64fb5d7970dbc9742197a11ea6f4cf4a096c` |

   Cache: `.cache/xlms_model_selection/`. The generator verifies the sha256 and refuses a file
   that differs. Selection criteria (measured in a survey of 30 candidates): both chains large
   enough that a wrong pose moves lysines out of range, >= 31 inter-chain Lys-Lys pairs with
   CA-CA <= 26 A, no insertion codes, no modified residues, residue numbers a constant shift of
   SEQRES.
2. **Full-length sequences.** `sequences.fasta` = a synthetic N-terminal region (14-94 residues,
   disordered-like composition, seeded) + the SEQRES sequence of the chain. Accessions are
   synthetic (`XL` + 4 characters). "Full-length position" = position in this sequence.
3. **Modelled constructs.** A model contains the residues that are resolved in the native with a
   complete backbone, minus resolved islands shorter than 12 residues, minus one planted internal
   deletion of 6-14 residues per protein (surface, > 14 A from the partner, containing a lysine
   when possible, placed after 15-60% of the chain, never where the deletion could be slid along a
   repeat). Lysines whose modelled/missing status would depend on where an aligner opens a gap
   (any gap, planted or native) are never used in a cross-link, so the position mapping is unique.
4. **Models.** Protein 1 is kept, protein 2 is moved as a rigid body.
   - 3 near-native models: rotation 1-9 degrees about the interface centre plus a shift <= 1.6 A,
     accepted at interface RMSD 0.5-1.9 A (targets 0.75 / 1.25 / 1.75 A).
   - 19-25 wrong models: random orientation, pushed against another part of the surface until the
     pose has 0.6-1.6 x the native number of CA-CA contacts < 8 A; accepted only with interface
     RMSD >= 6 A, and (checked afterwards) DockQ < 0.23 and fnat < 0.3. Nothing lies between the
     two classes.
   - Each model is written in its own random frame, heavy atoms only, `ATOM`/`TER`/`END` records
     only, B-factor column = a pLDDT-like value derived from the rank of the native B-factor
     (identical for all models up to noise, so it carries no quality information).
   - Level >= 2: residues are numbered 1..n per chain and the chain letters / chain order follow
     one of four layouts per model; all near-native models use a layout different from model_01.
     Level 1: full-length numbering, chains A/B everywhere.
   - Hidden quality of each model: interface RMSD (backbone atoms of residues within 10 A of the
     partner), fnat (5 A heavy-atom contacts), ligand RMSD after superposing the larger chain, and
     DockQ computed from these, all from the exact transform.
5. **Cross-links.** Sampled from the native geometry (CA-CA):
   - true inter-protein: 13-17 pairs with native distance <= cutoff - 4 A (<= 3 links per lysine);
   - true intra-protein: 24-32 pairs with native distance <= cutoff - 1.5 A;
   - up to 4 true pairs to a lysine of a planted deletion (in range in the native, but not evaluable);
   - 1 pair reaching into a region without any coordinates (not evaluable);
   - level 3 only: false positives = lysine pairs far apart in the native (inter >= cutoff + 8 A,
     intra >= cutoff + 3 A), 11.5-13% of the rows.
   The cross-linker and with it the cutoff vary with the seed (DSSO 30 A, BS3 30 A, DSBU 32 A).
6. **Reference counts.** After the model files are written the generator reads the CA atoms back
   from the text and counts with the planted mapping. It raises instead of emitting a task unless
   - every evaluated distance (all models x all evaluable links) is >= 0.5 A from the cutoff
     (poses are pre-filtered at 0.6 A);
   - every near-native model satisfies >= 3 more inter-protein links than any wrong model;
   - the inter-protein satisfied count takes >= 6 distinct values and no value is shared by more
     than 40% of the models;
   - the intra-protein counts are identical in all models (rigid chains);
   - at least two cross-linked positions are missing from the models.

Independent check of the labels (Biopython parser + SVDSuperimposer on CA atoms, starting from the
shipped model files and the cached native, seed 1 level 3): CA ligand RMSD agrees with the stored
backbone ligand RMSD within 0.58 A for all 24 models; near-native models 1.3-3.5 A, wrong models
>= 19.7 A; receptor fit 0.000-0.001 A.

## Decision points and why they are discoverable

| trap (level) | what goes wrong on the naive path | discoverable by |
|---|---|---|
| `numbering_offset_and_gap` (2) | positions are full-length positions; model residues are 1..n over trimmed constructs, so the shift changes after the trimmed loop | the residue at a cross-linked position is not a lysine (or does not exist) in the model; chain sequences are shorter than the FASTA sequences; notes: prediction input were trimmed constructs |
| `chain_identity_varies` (2) | chain letter and chain order of a protein differ between model files | chain lengths/sequences per letter differ between files; notes: models were pooled from several runs |
| `misleading_confidence` (3) | model_01 (top confidence) is wrong; the confidence order has AUC 0.30-0.60 for near-native vs wrong | satisfied inter-protein counts disagree with `ranking.json` |
| `false_positive_crosslinks` (3) | 11.5-13% of the pairs are false; the best model violates the 2-3 false inter-protein links and every model violates the false intra-protein links | notes state the FDR and that nothing was validated manually; no model satisfies everything |

The instruction names none of these. It defines the scored quantities (inter/intra by accession,
"evaluable", totals, satisfied "by the distance criterion that `xlms_notes.md` gives") and the
output files. The atoms (CA-CA) and the cutoff are stated in the notes as the group's convention;
the notes also restate that same-protein pairs are intra-protein (the sample is a 1:1 heterodimer)
and that a cross-link to a residue absent from a structure is left out of both the satisfied count
and the total. Those four statements are everything the exact counts depend on.

## Output contract

- `restraints.tsv`: `model`, `inter_satisfied`, `inter_total`, `intra_satisfied`, `intra_total`
  (extra columns and any row order allowed; `model_07.pdb` accepted for `model_07`).
- `selection.json`: `best_model`, `ranking` (permutation of all models, first element = best).
- `mapping.tsv`: `accession`, `position`, `model_chain`, `model_resnum` for every distinct
  cross-linked position, in the labelling of `model_01.pdb`; `NA` (any case, also empty, `-`,
  `None`, `nan`) when the residue is not in the model.

## Scoring (pass threshold 0.9)

Hard gates (score 0): missing/unparseable file; a model missing, duplicated or unknown in
`restraints.tsv`; non-integer counts; `selection.json` without a valid `best_model`, without a
permutation of all models as `ranking`, or with `ranking[0] != best_model`; duplicated positions
or non-integer positions in `mapping.tsv`.

| weight | component | metric | credit |
|---|---|---|---|
| 0.30 | restraint_counts | fraction of models with all four counts exact | 1.0 if all models exact, else 0.6 x fraction |
| 0.35 | best_model | hidden class of `best_model` | 1.0 if near-native (interface RMSD < 2 A), else 0 |
| 0.20 | ranking | ROC AUC of the submitted order, near-native vs wrong | linear from 0.70 (0) to 0.95 (1) |
| 0.15 | residue_mapping | fraction of cross-linked positions with exact chain + residue number (or NA) | 1.0 if all exact, else 0.6 x fraction |

The spec sketch had 0.35 / 0.35 / 0.20 / 0.10. With those weights a submission with a useless
`mapping.tsv` scored 0.35 + 0.35 + 0.20 = 0.8999999999999999, i.e. it failed only through floating
point rounding (measured on the host). Moving 0.05 from the counts to the mapping puts it at 0.85.

## Measured values (Docker runs of the final code, level 3, `details.json` of every run)

s1 = seed 1 (native 1FIN, DSBU, cutoff 32 A, 24 models, 48 cross-links); s2 = seed 2 (native 1EFU,
BS3, cutoff 30 A, 24 models, 52 cross-links). "models exact" = models whose four counts are all
exact; "iRMSD of best_model" = hidden interface RMSD of the submitted `best_model` (near-native
< 2 A); "positions exact" = mapping rows right.

| submission | s1: models exact | s1: iRMSD of best_model (A) | s1: AUC | s1: positions exact | s1: score | s2: models exact | s2: iRMSD of best_model (A) | s2: AUC | s2: positions exact | s2: score |
|---|---|---|---|---|---|---|---|---|---|---|
| `oracle` | 24/24 | 1.246 | 1.0 | 36/36 | **1.000** | 24/24 | 1.708 | 1.0 | 44/44 | **1.000** |
| `independent_solution` | 24/24 | 1.246 | 1.0 | 36/36 | **1.000** | 24/24 | 1.25 | 1.0 | 44/44 | **1.000** |
| `planted_truth` | 24/24 | 0.751 | 1.0 | 36/36 | **1.000** | 24/24 | 0.75 | 1.0 | 44/44 | **1.000** |
| `within_tolerance` | 24/24 | 0.751 | 0.9841 | 36/36 | **1.000** | 24/24 | 0.75 | 0.9841 | 44/44 | **1.000** |
| `naive_positions_as_resnums` | 0/24 | 18.2 | 0.9524 | 0/36 | **0.200** | 0/24 | 21.509 | 0.7143 | 0/44 | **0.011** |
| `naive_single_offset` | 0/24 | 10.783 | 0.7619 | 15/36 | **0.087** | 0/24 | 1.708 | 0.9841 | 14/44 | **0.579** |
| `naive_constant_chains` | 11/24 | 10.783 | 0.127 | 36/36 | **0.233** | 12/24 | 8.919 | 0.1429 | 44/44 | **0.240** |
| `naive_top_confidence` | 24/24 | 9.214 | 0.3175 | 36/36 | **0.450** | 24/24 | 8.919 | 0.4127 | 44/44 | **0.450** |
| `naive_prune_unsatisfied` | 0/24 | 1.246 | 1.0 | 36/36 | **0.700** | 0/24 | 1.708 | 1.0 | 44/44 | **0.700** |
| `naive_pooled_counts` | 0/24 | 1.246 | 1.0 | 36/36 | **0.700** | 0/24 | 1.708 | 1.0 | 44/44 | **0.700** |
| `naive_default_pipeline` | 0/24 | 9.214 | 0.3175 | 0/36 | **0.000** | 0/24 | 8.919 | 0.4127 | 0/44 | **0.000** |
| `prior_only_constant` | 0/24 | 9.214 | 0.3175 | 0/36 | **0.000** | 0/24 | 8.919 | 0.4127 | 0/44 | **0.000** |
| `one_count_off` | 23/24 | 0.751 | 1.0 | 36/36 | **0.873** | 23/24 | 0.75 | 1.0 | 44/44 | **0.873** |
| `best_is_best_wrong_model` | 24/24 | 9.214 | 0.9524 | 36/36 | **0.650** | 24/24 | 8.919 | 0.9524 | 44/44 | **0.650** |

Bands and where the measurements sit:

- restraint_counts: oracle / independent solution 24/24 on both seeds; every variant that gets the
  counts wrong at all gets 0/24 to 12/24. Credit is 1.0 only at 24/24, so one model off
  (`one_count_off`, 23/24) gives 0.575 and the task score 0.8725 < 0.9.
- best_model: oracle and independent solution pick different near-native models at seed 2 (1.71 A
  and 1.25 A) and both get full credit. In the 60-instance host sweep the worst near-native model
  is at 1.81 A (DockQ >= 0.54) and the best wrong model at 6.04 A (DockQ <= 0.18, fnat <= 0.25).
- ranking: correct solutions 1.0, the in-tolerance probe 0.984; the confidence order 0.32 / 0.41
  (planted in 0.30-0.60); band 0.70-0.95. Variants with a broken mapping land anywhere between 0.13
  and 0.98, which is why the ranking has only 0.20 weight.
- residue_mapping: correct solutions 36/36 and 44/44; single-offset 15/36 and 14/44; positions as
  residue numbers 0.

Release conditions measured by the generator on the written files:

| instance | min distance-to-cutoff (A) | near-native inter counts | best wrong inter count | distinct inter counts | false rows | confidence/name-order AUC | cross-linked positions (missing from models) |
|---|---|---|---|---|---|---|---|
| seed 1, level 3 | 0.630 | 15, 15, 15 of 17 | 12 | 14 | 12.5% | 0.317 | 36 (3) |
| seed 2, level 3 | 0.656 | 16, 16, 16 of 18 | 13 | 14 | 11.5% | 0.413 | 44 (3) |

Host sweep of the final generator over seeds 1-20 at levels 1, 2 and 3: 60 of 60 instances
generated. Level 3: min distance-to-cutoff 0.600-0.742 A, lead of the near-native models 3-6 links,
11-16 distinct inter counts, false rows 11.5-13.0%, confidence-order AUC 0.317-0.600, 1-9 wrong
models per instance satisfy at least one false inter-protein link, model_01 never near-native,
input 7.5-13.3 MB. Levels 1-2: min distance-to-cutoff 0.603-0.854 A, lead 3-5 links, name-order
AUC 0.302-0.591. (An earlier sweep of 40 instances had one failure, see the worklog.)

Contact statistics that motivated the pose generator (seed 1 / seed 2, level 3): CA-CA pairs < 8 A
across the interface, near-native 58, 83, 124 / 79, 84, 149; wrong 51-108 / 50-126. CA-CA pairs
< 4 A, near-native 0, 1, 10 / 0, 4, 8; wrong 0-20 / 2-21.

### Gate lines

`python -m alegen make xlms_model_selection --seed 1 --level 3 --keep-work`

```
wrote /home/user/ale-lifesci-gym/tasks/xlms-model-selection-s0001-l3  (prior difficulty 0.565)
  [ok] G1 images_build 
  [ok] G2 oracle score=1.0 runtime_min=0.04 failure= stderr=
  [ok] G5 verifier_deterministic first=1.0 second=1.0
  [ok] G3 probe:planted_truth score=1.0 floor=0.999 failure=
  [ok] G3 probe:within_tolerance score=1.0 floor=0.99 failure=
  [ok] G3 probe:independent_solution score=1.0 floor=0.99 failure=
  [ok] G4 attack:empty_output score=0.0 ceiling=0.0 failure=required output file missing: restraints.tsv
  [ok] G4 attack:inputs_copied score=0.0 ceiling=0.1 failure=required output file missing: restraints.tsv
  [ok] G4 attack:zero_byte_files score=0.0 ceiling=0.0 failure=restraints.tsv is empty
  [ok] G4 attack:header_only score=0.0 ceiling=0.1 failure=gate failed: restraints.tsv lists every model exactly once
  [ok] G4 attack:prior_only_constant score=0.0 ceiling=0.1 failure=
  [ok] G4 attack:one_count_off score=0.8725 ceiling=0.89 failure=
  [ok] G4 attack:best_is_best_wrong_model score=0.65 ceiling=0.7 failure=
  [ok] G4 attack:naive_pooled_counts score=0.7 ceiling=0.75 failure=
  [ok] G4 attack:naive_positions_as_resnums score=0.2 ceiling=0.75 failure=
  [ok] G4 attack:naive_single_offset score=0.087024 ceiling=0.75 failure=
  [ok] G4 attack:naive_constant_chains score=0.2325 ceiling=0.75 failure=
  [ok] G4 attack:naive_top_confidence score=0.45 ceiling=0.75 failure=
  [ok] G4 attack:naive_prune_unsatisfied score=0.7 ceiling=0.75 failure=
  [ok] G4 attack:naive_default_pipeline score=0.0 ceiling=0.6 failure=
  [ok] G6 harbor_oracle rewards=[1.0] tail=
  [ok] G6 harbor_nop rewards=[0.0] tail=
{"task": "xlms-model-selection-s0001-l3", "ok": true, "minutes": 2.93}
exit=0
```

`python -m alegen make xlms_model_selection --seed 2 --level 3 --no-harbor`

```
wrote /home/user/ale-lifesci-gym/tasks/xlms-model-selection-s0002-l3  (prior difficulty 0.565)
  [ok] G1 images_build 
  [ok] G2 oracle score=1.0 runtime_min=0.05 failure= stderr=
  [ok] G5 verifier_deterministic first=1.0 second=1.0
  [ok] G3 probe:planted_truth score=1.0 floor=0.999 failure=
  [ok] G3 probe:within_tolerance score=1.0 floor=0.99 failure=
  [ok] G3 probe:independent_solution score=1.0 floor=0.99 failure=
  [ok] G4 attack:empty_output score=0.0 ceiling=0.0 failure=required output file missing: restraints.tsv
  [ok] G4 attack:inputs_copied score=0.0 ceiling=0.1 failure=required output file missing: restraints.tsv
  [ok] G4 attack:zero_byte_files score=0.0 ceiling=0.0 failure=restraints.tsv is empty
  [ok] G4 attack:header_only score=0.0 ceiling=0.1 failure=gate failed: restraints.tsv lists every model exactly once
  [ok] G4 attack:prior_only_constant score=0.0 ceiling=0.1 failure=
  [ok] G4 attack:one_count_off score=0.8725 ceiling=0.89 failure=
  [ok] G4 attack:best_is_best_wrong_model score=0.65 ceiling=0.7 failure=
  [ok] G4 attack:naive_pooled_counts score=0.7 ceiling=0.75 failure=
  [ok] G4 attack:naive_positions_as_resnums score=0.011429 ceiling=0.75 failure=
  [ok] G4 attack:naive_single_offset score=0.578636 ceiling=0.75 failure=
  [ok] G4 attack:naive_constant_chains score=0.24 ceiling=0.75 failure=
  [ok] G4 attack:naive_top_confidence score=0.45 ceiling=0.75 failure=
  [ok] G4 attack:naive_prune_unsatisfied score=0.7 ceiling=0.75 failure=
  [ok] G4 attack:naive_default_pipeline score=0.0 ceiling=0.6 failure=
{"task": "xlms-model-selection-s0002-l3", "ok": true, "minutes": 0.71}
exit=0
```

`python -m alegen make xlms_model_selection --seed 1 --level 1 --no-harbor`

```
wrote /home/user/ale-lifesci-gym/tasks/xlms-model-selection-s0001-l1  (prior difficulty 0.165)
  [ok] G1 images_build 
  [ok] G2 oracle score=1.0 runtime_min=0.03 failure= stderr=
  [ok] G5 verifier_deterministic first=1.0 second=1.0
  [ok] G3 probe:planted_truth score=1.0 floor=0.999 failure=
  [ok] G3 probe:within_tolerance score=1.0 floor=0.99 failure=
  [ok] G3 probe:independent_solution score=1.0 floor=0.99 failure=
  [ok] G4 attack:empty_output score=0.0 ceiling=0.0 failure=required output file missing: restraints.tsv
  [ok] G4 attack:inputs_copied score=0.0 ceiling=0.1 failure=required output file missing: restraints.tsv
  [ok] G4 attack:zero_byte_files score=0.0 ceiling=0.0 failure=restraints.tsv is empty
  [ok] G4 attack:header_only score=0.0 ceiling=0.1 failure=gate failed: restraints.tsv lists every model exactly once
  [ok] G4 attack:prior_only_constant score=0.082286 ceiling=0.1 failure=
  [ok] G4 attack:one_count_off score=0.873077 ceiling=0.89 failure=
  [ok] G4 attack:best_is_best_wrong_model score=0.65 ceiling=0.7 failure=
  [ok] G4 attack:naive_pooled_counts score=0.7 ceiling=0.75 failure=
{"task": "xlms-model-selection-s0001-l1", "ok": true, "minutes": 0.4}
exit=0
```

Additional run, not required: `--seed 1 --level 2 --no-harbor`

```
wrote /home/user/ale-lifesci-gym/tasks/xlms-model-selection-s0001-l2  (prior difficulty 0.432)
  [ok] G1 images_build 
  [ok] G2 oracle score=1.0 runtime_min=0.03 failure= stderr=
  [ok] G5 verifier_deterministic first=1.0 second=1.0
  [ok] G3 probe:planted_truth score=1.0 floor=0.999 failure=
  [ok] G3 probe:within_tolerance score=1.0 floor=0.99 failure=
  [ok] G3 probe:independent_solution score=1.0 floor=0.99 failure=
  [ok] G4 attack:empty_output score=0.0 ceiling=0.0 failure=required output file missing: restraints.tsv
  [ok] G4 attack:inputs_copied score=0.0 ceiling=0.1 failure=required output file missing: restraints.tsv
  [ok] G4 attack:zero_byte_files score=0.0 ceiling=0.0 failure=restraints.tsv is empty
  [ok] G4 attack:header_only score=0.0 ceiling=0.1 failure=gate failed: restraints.tsv lists every model exactly once
  [ok] G4 attack:prior_only_constant score=0.0 ceiling=0.1 failure=
  [ok] G4 attack:one_count_off score=0.873077 ceiling=0.89 failure=
  [ok] G4 attack:best_is_best_wrong_model score=0.65 ceiling=0.7 failure=
  [ok] G4 attack:naive_pooled_counts score=0.7 ceiling=0.75 failure=
  [ok] G4 attack:naive_positions_as_resnums score=0.0 ceiling=0.75 failure=
  [ok] G4 attack:naive_single_offset score=0.238571 ceiling=0.75 failure=
  [ok] G4 attack:naive_constant_chains score=0.219231 ceiling=0.75 failure=
{"task": "xlms-model-selection-s0001-l2", "ok": true, "minutes": 1.41}
exit=0
```

## Known limitations and open points

- **Pooling inter + intra cannot change the order.** The models are strict rigid-body variants, so
  the intra-protein counts are identical in every model (asserted). The spec's naive variant
  "intra+inter pooled" is therefore only a reporting error: `naive_pooled_counts` loses the count
  component and nothing else (0.70, ceiling 0.75). Making pooling matter needs models that differ
  inside the chains (flexible tail or domain); not built.
- **The false-positive trap does not change the selection for a count-based analysis.** It adds
  noise (wrong models satisfy false inter-protein links: 2 of 21 wrong models at seed 1, 6 of 21
  at seed 2) and it
  makes "satisfies everything" unattainable, but the only naive variant I could make concrete
  (`naive_prune_unsatisfied`: discard the links no model satisfies) fails the count component only
  (0.70). Its margin to the 0.75 ceiling is 0.05.
- **Physical plausibility is a second route to `best_model`.** Wrong poses are matched to the
  native in CA-CA contact number and they clash about as much as the perturbed near-native poses,
  but their side chains were not repacked. A real interface-scoring function would probably pick
  the near-native models without any cross-link. That can give 0.55 (best_model + ranking); the
  counts and the mapping (0.45) still need the cross-link analysis, so it cannot pass.
- **Protein identity.** Accessions are synthetic and there is no header, but the sequences are the
  real SEQRES sequences behind a synthetic N-terminal region. A model that recognises the complex
  and remembers its interface could reason about `best_model` without the data (not tested).
- **Garbage mappings keep some signal.** With a small offset, "position = residue number" still
  points to residues near the right ones; at seed 10 level 2 (host run) that variant picked a
  near-native model and scored 0.55. Such variants are bounded by 0.55 (0.64 when part of the
  mapping is right) because the counts fail, but their best_model/ranking credit can be luck.
- **Ranking is scored only as near-native vs wrong.** The order among wrong models and among
  near-native models is not scored. The three near-native models always satisfy the same number of
  inter-protein links in the instances I looked at.
- **Exactness is strict.** One miscounted model fails the task (0.8725). This is intended and safe
  only because of the 0.5 A margin; a solver that uses another atom or cutoff than the notes state
  fails by design.
- **Level 1 is easy by construction**: identity mapping, constant chains, no false positives, no
  confidence file. A constant submission gets 0-0.084 there from the identity mapping (host runs,
  seeds 1-5; bounded by 0.09).
- **Realism.** The "full-length" sequence is not a real database sequence; chains have a physical
  break at the trimmed loop while numbering runs through; all models share identical chain
  geometry; `n_csm`/`score` are drawn from two overlapping distributions and are not used by the
  truth.
- **File-name order is planted at every level.** model_01 is never near-native and the name order
  separates near-native from wrong models with AUC 0.30-0.60 (this keeps a constant answer at 0).
  It is a weak learnable prior ("the best model is not the first file"); it does not identify the
  best model.
- **Level 2 went through Docker once** (seed 1, all gates ok, no Harbor run); it is otherwise
  covered by the host sweep only.
- **Five natives.** Seeds 6+ reuse a native with new trims, links, poses, labels and cutoff.
- The oracle reads the cutoff from the notes with a regular expression tied to the generator's
  wording.
- Hidden quality metrics come from my own numpy implementation (iRMSD / fnat / LRMSD / DockQ);
  only the ligand RMSD was cross-checked with a second tool (above). The class labels do not
  depend on the fine values (0.75-1.75 A vs >= 6 A).

---

# Level 4: cross-links that a straight line misjudges (P_surface_distance + P_homo_oligomer)

Code: `pipeline/alegen/families/xlms_model_selection/l4.py` (generator, imported by `generate` for
`level >= 4`), `assets/l4/` (`score.py`, `solve.sh`, `xl4_pipeline.py` = reference, `solve_alt.py` =
second solution and careless variants, `sasd_np.py` = lattice search used by the generator and the
second solution). Levels 1-3 are untouched apart from a three-line dispatch; `make --seed 1 --level 3
--out <scratch> --no-validate` followed by `diff -r` against `tasks/xlms-model-selection-s0001-l3`
shows no difference.

## What the pilot showed and what level 4 changes

Level 3, seed 1: Claude Haiku scored 1.0 in 1 min 35 s with 15 tool calls. Every decision point was
one that the notes or a sanity check reveal, and once the numbering was mapped the rest was a
straight-line distance count. README section 6 lists what the small model did fail on elsewhere; the
one that applies here is "a plausible shortcut that the generator had closed".

| element | pilot finding it answers | what it looks like here |
|---|---|---|
| `P_surface_distance` | The level-3 count needs one line of numpy per link. The obvious default (straight-line CA-CA distance with the stated cutoff) was also the right answer. | The group's criterion is a surface distance: shortest path through solvent on a 1 A lattice, written out step by step in `xlms_notes.md`. No program for it is installed. The straight line is now the closed shortcut: with the same cutoff it gets no model's counts right and ranks a wrong model (model_01, the top-confidence one) first. |
| `P_homo_oligomer` | "Which chain is which protein" was solved once per file and then never mattered again. | The complex is 2:2. A cross-link is satisfied if any assignment of its two positions to chains is within range; a row with the same position on both sides can only run between the two copies. The rule is in the instruction and in the notes. Looking two positions of one protein up in one chain ranks a wrong model first. |

Both are P elements (an authoritative-looking default that is wrong in this instance, with the right
rule stated in the SOP). The level-3 decision points (numbering, chain identity per file, misleading
confidence, false positives) are all still present.

## Scenario and truth model at level 4

- **Hidden native.** Two heterodimers of one RCSB entry (same cached files and sha256 as level 3):
  1FIN chains A,B + C,D for odd seeds, 1EFU chains A,B + C,D for even seeds. 1IBR (heterodimers do
  not touch), 1KXP (one heterodimer) and 1I2M (14 candidate between-copy pairs) are not used.
- **Constructs.** Residues resolved with a complete backbone in both copies, minus islands shorter
  than 12, minus one planted internal deletion per protein (same residues in both copies). Both
  copies are numbered 1..n identically.
- **Models.** Heterodimer 1 is kept, heterodimer 2 is moved as a rigid body. A pool of 6 near-native
  poses (two each at interface RMSD 0.75 / 1.25 / 1.75 A; one of each is kept) and 36 wrong poses
  (interface RMSD >= 6 A, DockQ < 0.23, fnat < 0.3, 0.6-1.6 x the native number of CA-CA contacts
  < 8 A between the heterodimers) is tabulated completely; 15-18 models are then chosen. Each model
  has its own random frame, chain order and chain letters (18 layouts).
- **Surface distance (the definition as shipped in `xlms_notes.md`).** Lattice points = points with
  whole-number x, y, z in the frame of the file; a point is blocked if it is less than 3.0 A from the
  centre of any non-hydrogen atom, otherwise open; a path is a sequence of open points that differ by
  at most 1 in each coordinate, its length is the sum of the straight-line step lengths; each CA atom
  is joined to the first / last point by a straight segment of at most 7.0 A that is not tested
  against the protein; the surface distance is the minimum of segment + path + segment over all
  choices, infinite if there is none. Cutoff by cross-linker: DSSO 32 A, BS3 34 A, DSBU 36 A.
- **Copies.** Satisfied = at least one assignment of the two positions to chains is within range;
  the two ends are always two different residues.
- **Cross-links** (52 at seed 1, 59 at seed 2):
  - *true, between the heterodimers* (8 / 13): within range in the native only through an assignment
    across the two heterodimers, and satisfied in all three near-native models. Includes rows with
    identical ends (5 / 2) and, at seed 2, pairs of two different positions of one protein across the
    copies; at seed 1 (1FIN) no such pair exists in the native.
  - *true, inside a heterodimer* (31 / 33): within one chain or between the two proteins of one
    heterodimer.
  - *false* (10 / 10, i.e. 19% / 17% of the rows; the notes state this rate): surface distance in the
    native > cutoff + 8 A and violated in all near-native models. Four kinds: 3 pairs that every model
    has within straight-line range but out of range along the surface (the path would pass through
    the protein); 1-2 hetero pairs that one wrong model satisfies along the surface; 3-4 pairs that
    are within straight-line range in model_01 and in no near-native model; 1 pair that is far by
    both measures; plus 1 row into a region without coordinates.
  - 3 true rows to a lysine of a trimmed loop (not evaluable).
- **Status truth.** Computed with the numpy implementation on the coordinates as written (rounded to
  0.001 A). This is a computed quantity, not a planted one, so three things are enforced:
  1. a pair is used only if its surface distance is <= cutoff - 3 A or >= cutoff + 3 A in every
     shipped model (no gray-zone call exists, so no tolerance is needed in the scorer);
  2. a model frame is redrawn if any lattice point is exactly 3.000 A from an atom or exactly 7.000 A
     from a lysine CA (integer arithmetic on the written coordinates), so "less than" / "at most"
     cannot depend on rounding;
  3. the reference solution is a separate implementation (Biopython, k-d tree, sparse-graph Dijkstra
     in scipy) and gate G2 requires it to reproduce every count and status. The generator does not
     run it; an instance on which the two disagree fails G2 and is not released.
- **Release conditions checked by the generator** (otherwise it redraws the whole instance, up to 6
  attempts): every near-native model satisfies >= 3 more cross-links than any wrong model; with
  straight-line distances (copies handled correctly) a wrong model has more satisfied cross-links
  than every near-native model; the same under "two positions of one protein are looked up in one
  chain"; the same for the inter-protein count of the level-3 method (one chain per protein, straight
  line), which is arranged through the chain layouts; >= 5 distinct totals and no total shared by
  more than 60% of the models; model_01 is the straight-line look-alike and the file order separates
  near-native from wrong models with AUC 0.30-0.60; >= 2 cross-linked positions missing from the
  models.

## Output contract and scoring (pass 0.9)

As level 3 plus `link_status.tsv` (`model`, `xl_id`, `status` in `satisfied` / `violated` /
`not_evaluable`, one row per model and cross-link). `ranking` is by the total number of satisfied
cross-links (stated in the instruction; ties in any order). `model_chain` in `mapping.tsv` lists both
chains of the protein in model_01, comma-separated (the scorer accepts `,` `;` `/` `|` or blanks and
any order).

| weight | component | metric | credit |
|---|---|---|---|
| 0.30 | restraint_counts | fraction of models with all four counts exact | 1.0 if all models exact, else 0.6 x fraction |
| 0.15 | link_status | three-class Matthews correlation over all model x cross-link statuses | linear from 0.95 (0) to 0.995 (1) |
| 0.30 | best_model | hidden class of `best_model` | 1.0 if near-native (interface RMSD < 2 A), else 0 |
| 0.10 | ranking | ROC AUC of the submitted order, near-native vs wrong | linear from 0.70 (0) to 0.95 (1) |
| 0.15 | residue_mapping | fraction of positions with the exact set of chains + residue number (or NA) | 1.0 if all exact, else 0.6 x fraction |

Hard gates as at level 3, plus: `link_status.tsv` must hold exactly one row per model and cross-link
with a valid status.

Band of the status component: correct implementations give exactly 1.0 (0 wrong statuses, both
seeds, both implementations); the two "ignore one element" variants give 0.82-0.86; the band is
0.95-0.995. A submission with zero status credit scores at most 0.55 + 0.18 x f (f = fraction of
models with exact counts) even if it picks a near-native model and has the mapping right; the
highest such value measured is 0.556 (`naive_constant_chains`, seed 2).


## Measured values (Docker runs of the final code, `results/work/<task>/**/details.json`)

s1 = seed 1 (native 1FIN, cutoff 34 A, 17 models, 52 cross-links, 36 cross-linked positions of which 3 are missing from the models); s2 = seed 2 (native 1EFU, cutoff 34 A, 17 models, 59 cross-links, 44 cross-linked positions of which 3 are missing from the models). "models exact" = models whose four counts are all exact; "statuses wrong" = rows of `link_status.tsv` that differ from the truth; "iRMSD of best_model" = hidden interface RMSD of the submitted `best_model` (near-native < 2 A).

| submission | s1: models exact | s1: statuses wrong | s1: status MCC | s1: iRMSD of best_model (A) | s1: AUC | s1: positions exact | s1: score | s2: models exact | s2: statuses wrong | s2: status MCC | s2: iRMSD of best_model (A) | s2: AUC | s2: positions exact | s2: score |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `oracle` | 17/17 | 0 | 1.0 | 1.25 | 1.0 | 36/36 | **1.000** | 17/17 | 0 | 1.0 | 1.248 | 1.0 | 44/44 | **1.000** |
| `independent_solution` | 17/17 | 0 | 1.0 | 1.25 | 1.0 | 36/36 | **1.000** | 17/17 | 0 | 1.0 | 1.248 | 1.0 | 44/44 | **1.000** |
| `planted_truth` | 17/17 | 0 | 1.0 | 0.749 | 1.0 | 36/36 | **1.000** | 17/17 | 0 | 1.0 | 0.75 | 1.0 | 44/44 | **1.000** |
| `within_tolerance` | 17/17 | 0 | 1.0 | 1.25 | 0.9762 | 36/36 | **1.000** | 17/17 | 0 | 1.0 | 1.248 | 0.9762 | 44/44 | **1.000** |
| `l3_reference` | 0/17 | 193 | 0.5529 | 12.222 | 0.2857 | 3/36 | **0.007** | 0/17 | 277 | 0.4957 | 41.749 | 0.5714 | 3/44 | **0.006** |
| `naive_straight_line` | 0/17 | 74 | 0.8397 | 12.222 | 0.9286 | 36/36 | **0.241** | 0/17 | 83 | 0.8455 | 9.533 | 0.9286 | 44/44 | **0.241** |
| `naive_same_protein_same_chain` | 0/17 | 82 | 0.8221 | 18.775 | 0.9286 | 36/36 | **0.241** | 2/17 | 72 | 0.8604 | 9.533 | 0.9286 | 44/44 | **0.263** |
| `naive_one_chain_per_protein` | 0/17 | 149 | 0.6647 | 12.222 | 0.2857 | 3/36 | **0.007** | 0/17 | 244 | 0.594 | 41.749 | 0.6429 | 3/44 | **0.006** |
| `naive_pooled_counts` | 0/17 | 0 | 1.0 | 1.25 | 1.0 | 36/36 | **0.700** | 0/17 | 0 | 1.0 | 1.248 | 1.0 | 44/44 | **0.700** |
| `naive_positions_as_resnums` | 0/17 | 573 | 0.1962 | 1.25 | 0.8571 | 0/36 | **0.363** | 0/17 | 712 | -0.0383 | 42.346 | 0.0476 | 0/44 | **0.000** |
| `naive_single_offset` | 0/17 | 288 | 0.4754 | 1.698 | 1.0 | 20/36 | **0.450** | 0/17 | 413 | 0.3717 | 0.75 | 1.0 | 14/44 | **0.429** |
| `naive_constant_chains` | 8/17 | 198 | 0.6425 | 15.744 | 0.619 | 36/36 | **0.235** | 10/17 | 126 | 0.7821 | 0.75 | 0.4762 | 44/44 | **0.556** |
| `naive_top_confidence` | 17/17 | 0 | 1.0 | 12.222 | 0.4286 | 36/36 | **0.600** | 17/17 | 0 | 1.0 | 9.533 | 0.4762 | 44/44 | **0.600** |
| `naive_prune_unsatisfied` | 0/17 | 102 | 0.8033 | 1.25 | 1.0 | 36/36 | **0.550** | 0/17 | 85 | 0.8543 | 1.248 | 1.0 | 44/44 | **0.550** |
| `naive_default_pipeline` | 0/17 | 576 | 0.0515 | 12.222 | 0.4286 | 0/36 | **0.000** | 0/17 | 639 | -0.0202 | 9.533 | 0.4762 | 0/44 | **0.000** |
| `prior_only_constant` | 0/17 | 630 | 0.0 | 12.222 | 0.4286 | 0/36 | **0.000** | 0/17 | 684 | 0.0 | 9.533 | 0.4762 | 0/44 | **0.000** |
| `one_count_off` | 16/17 | 0 | 1.0 | 0.749 | 1.0 | 36/36 | **0.869** | 16/17 | 0 | 1.0 | 0.75 | 1.0 | 44/44 | **0.869** |
| `best_is_best_wrong_model` | 17/17 | 0 | 1.0 | 12.222 | 0.9286 | 36/36 | **0.691** | 17/17 | 0 | 1.0 | 9.533 | 0.9286 | 44/44 | **0.691** |
| `empty_output` | - | - | - | - | - | - | **0.000** (gate) | - | - | - | - | - | - | **0.000** (gate) |
| `inputs_copied` | - | - | - | - | - | - | **0.000** (gate) | - | - | - | - | - | - | **0.000** (gate) |
| `zero_byte_files` | - | - | - | - | - | - | **0.000** (gate) | - | - | - | - | - | - | **0.000** (gate) |
| `header_only` | - | - | - | - | - | - | **0.000** (gate) | - | - | - | - | - | - | **0.000** (gate) |


All four generic attacks (`empty_output`, `inputs_copied`, `zero_byte_files`, `header_only`) fail a gate (score 0).

Reading the table:

- Reference and second solution: 0 wrong statuses and all counts exact on both seeds.
- `l3_reference` (straight line, one chain per protein, ranking by the inter-protein count): no model's counts right, a wrong model first, mapping wrong for every modelled position (it reports one chain).
- `naive_straight_line` (ignores only P_surface_distance): no model's counts right; it selects model_01, which is wrong.
- `naive_same_protein_same_chain` (ignores only P_homo_oligomer): selects a wrong model on both seeds.
- The level-2/3 variants (`naive_positions_as_resnums`, `naive_single_offset`, `naive_constant_chains`, `naive_top_confidence`, `naive_prune_unsatisfied`) can still pick a near-native model by luck; they are held below 0.75 by the counts and statuses (highest: 0.70 for `naive_pooled_counts`, which has the statuses right and only reports the counts in the wrong columns).


Release conditions measured by the generator (`reference/truth.json`, `measured`):

| quantity | seed 1 | seed 2 |
|---|---|---|
| near-native lead in satisfied cross-links (total) | 4 | 4 |
| lead of a wrong model with straight-line distances | 1 | 2 |
| lead of a wrong model with "two positions of one protein in one chain" | 1 | 2 |
| lead of a wrong model in the level-3 method's inter-protein count | 5 | 9 |
| distinct totals over the models | 6 | 5 |
| true links between the heterodimers (of which hetero / same position) | 8 (3 / 5) | 13 (4 / 2) |
| true links inside a heterodimer | 31 | 33 |
| false: through the protein / satisfied by a wrong model / within straight-line range in model_01 but in no near-native model / far | 3 / 1 / 4 / 1 | 3 / 2 / 3 / 1 |
| false rows | 19.2% | 16.9% |
| confidence-order AUC | 0.429 | 0.476 |
| smallest margin of a satisfied call, numpy (A) | 3.433 | 3.344 |
| satisfied / violated calls | 562 / 254 | 616 / 319 |
| generator attempt used | 0 | 0 |

Totals of satisfied cross-links: seed 1 near-native [39], wrong [35, 35, 34, 32, 31, 31, 31, 31, 31, 31, 31, 31, 31, 30]; seed 2 near-native [46], wrong [42, 36, 36, 34, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33].

### Margin of the scored surface distances to the cutoff, under both implementations

"Reference" = `xl4_pipeline.py` (Biopython reader, i.e. float32 coordinates; k-d tree; scipy
Dijkstra on a sparse graph). "Second" = `solve_alt.py` + `sasd_np.py` (own reader, float64; bucketed
search in numpy; the generator's truth uses the same module). Both write a diagnostic
`distances.tsv` (smallest surface distance over the allowed assignments per model and cross-link).
In the Docker runs the search stops at cutoff + 3.5 A, so the violated side is resolved only up to
there; for the full distribution both programs were rerun on the same inputs on the host with the
search limit at cutoff + 12 A (`XL_LIMIT_EXTRA=12`; reference in a scratch venv with scipy 1.18.1 and
Biopython 1.85). Margin = cutoff - distance for satisfied calls, distance - cutoff for violated calls.

| | seed 1 reference | seed 1 second | seed 2 reference | seed 2 second |
|---|---|---|---|---|
| scored calls (model x evaluable cross-link) | 816 | 816 | 935 | 935 |
| calls on the wrong side of the cutoff | 0 | 0 | 0 | 0 |
| calls with a margin below 3 A | 0 | 0 | 0 | 0 |
| satisfied calls | 562 | 562 | 616 | 616 |
| smallest margin, satisfied (A) | 3.433 | 3.433 | 3.344 | 3.344 |
| satisfied, margin [3,4) / [4,6) / [6,10) / >= 10 A | 9 / 59 / 133 / 361 | 9 / 59 / 133 / 361 | 6 / 73 / 150 / 387 | 6 / 73 / 150 / 387 |
| violated calls | 254 | 254 | 319 | 319 |
| smallest margin, violated (A) | 3.003 | 3.003 | 3.215 | 3.215 |
| violated, margin [3,4) / [4,6) / [6,10) / [10,12] / > 12 A or no path | 3 / 19 / 22 / 23 / 187 | 3 / 19 / 22 / 23 / 187 | 8 / 25 / 45 / 34 / 207 | 8 / 25 / 45 / 34 / 207 |
| calls finite in both programs | 629 | | 728 | |
| largest difference between the two programs on those (A) | 0.001 | | 0.000 | |
| calls finite in one program only | 0 | | 0 | |

The 0.001 A at seed 1 is one call whose two printed values (three decimals) round differently.
The `distances.tsv` files of the Docker runs (`results/work/<task>/oracle/output`,
`.../probe/independent_solution/output`) give the same picture within their shorter search limit:
satisfied side identical to the table; violated side 2 (seed 1) / 5 (seed 2) calls between 3.0 and
3.5 A beyond the cutoff, all others beyond the limit; largest difference between the programs
0.001 / 0.000 A on 564 / 621 calls finite in both.

The smallest violated margin at seed 1 is 3.003 A, i.e. the generator's 3 A filter is binding there.

Statuses that change when one parameter of the definition is perturbed (second solution, notes
edited, host): blocking radius 2.9 A: 0 of 884 (seed 1) / 0 of 1003 (seed 2); 3.1 A: 3 / 0; tether
6.5 A: 13 / 0; tether 7.5 A: 1 / 0.

Host generation of two further seeds (not validated in Docker): seed 3 (1FIN, cutoff 36 A) needed 4
attempts (431 s in total, three rejected for lack of compatible wrong poses, a total shared by 11
models, and no usable decoy pose) and then gave 17 models, 54 cross-links, lead 4, straight-line lead
of a wrong model 1, one-chain lead 2; seed 4 (1EFU, cutoff 32 A) passed on the first attempt (131 s):
16 models, 57 cross-links, lead 3, straight-line lead 2, one-chain lead 2.


### Runtime, sizes

| | seed 1 | seed 2 |
|---|---|---|
| reference solution in the agent image (G2, minutes; load average on the 4-core host was 10-15 during both runs) | 2.81 | 6.59 |
| reference solution on the host outside Docker, scratch venv (s) | 115 | 178 |
| second solution on the host (s) | 30 | 38 |
| whole validation (minutes) | 18.1 (with Harbor) | 13.73 (no Harbor) |
| generation on the host (s) | 143 | 162 |
| `environment/input` (MB) | 11.7 | 12.7 |
| agent image | 1.05 GB | 1.06 GB |

### Gate lines

`python -m alegen make xlms_model_selection --seed 1 --level 4 --keep-work`

```
wrote /home/user/ale-lifesci-gym/tasks/xlms-model-selection-s0001-l4  (prior difficulty 0.65)
  [ok] G1 images_build 
  [ok] G7 tools_in_login_shell detail=nonlogin=/opt/tools/bin/python3 login=/opt/tools/bin/python3
  [ok] G2 oracle score=1.0 runtime_min=2.81 failure= stderr=
  [ok] G5 verifier_deterministic first=1.0 second=1.0
  [ok] G3 probe:planted_truth score=1.0 floor=0.999 failure=
  [ok] G3 probe:within_tolerance score=1.0 floor=0.99 failure=
  [ok] G3 probe:independent_solution score=1.0 floor=0.99 failure=
  [ok] G4 has_l3_reference_attack 
  [ok] G4 attack:empty_output score=0.0 ceiling=0.0 failure=required output file missing: restraints.tsv
  [ok] G4 attack:inputs_copied score=0.0 ceiling=0.1 failure=required output file missing: restraints.tsv
  [ok] G4 attack:zero_byte_files score=0.0 ceiling=0.0 failure=restraints.tsv is empty
  [ok] G4 attack:header_only score=0.0 ceiling=0.1 failure=gate failed: restraints.tsv lists every model exactly once
  [ok] G4 attack:prior_only_constant score=0.0 ceiling=0.1 failure=
  [ok] G4 attack:one_count_off score=0.869412 ceiling=0.89 failure=
  [ok] G4 attack:best_is_best_wrong_model score=0.691429 ceiling=0.72 failure=
  [ok] G4 attack:l3_reference score=0.0075 ceiling=0.75 failure=
  [ok] G4 attack:naive_straight_line score=0.241429 ceiling=0.75 failure=
  [ok] G4 attack:naive_same_protein_same_chain score=0.241429 ceiling=0.75 failure=
  [ok] G4 attack:naive_one_chain_per_protein score=0.0075 ceiling=0.75 failure=
  [ok] G4 attack:naive_pooled_counts score=0.7 ceiling=0.75 failure=
  [ok] G4 attack:naive_positions_as_resnums score=0.362857 ceiling=0.75 failure=
  [ok] G4 attack:naive_single_offset score=0.45 ceiling=0.75 failure=
  [ok] G4 attack:naive_constant_chains score=0.234706 ceiling=0.75 failure=
  [ok] G4 attack:naive_top_confidence score=0.6 ceiling=0.75 failure=
  [ok] G4 attack:naive_prune_unsatisfied score=0.55 ceiling=0.75 failure=
  [ok] G4 attack:naive_default_pipeline score=0.0 ceiling=0.6 failure=
  [ok] G6 harbor_oracle rewards=[1.0] tail=
  [ok] G6 harbor_nop rewards=[0.0] tail=
{"task": "xlms-model-selection-s0001-l4", "ok": true, "minutes": 18.1}
exit=0
```

`python -m alegen make xlms_model_selection --seed 2 --level 4 --no-harbor --keep-work`

```
wrote /home/user/ale-lifesci-gym/tasks/xlms-model-selection-s0002-l4  (prior difficulty 0.655)
  [ok] G1 images_build 
  [ok] G7 tools_in_login_shell detail=nonlogin=/opt/tools/bin/python3 login=/opt/tools/bin/python3
  [ok] G2 oracle score=1.0 runtime_min=6.59 failure= stderr=
  [ok] G5 verifier_deterministic first=1.0 second=1.0
  [ok] G3 probe:planted_truth score=1.0 floor=0.999 failure=
  [ok] G3 probe:within_tolerance score=1.0 floor=0.99 failure=
  [ok] G3 probe:independent_solution score=1.0 floor=0.99 failure=
  [ok] G4 has_l3_reference_attack 
  [ok] G4 attack:empty_output score=0.0 ceiling=0.0 failure=required output file missing: restraints.tsv
  [ok] G4 attack:inputs_copied score=0.0 ceiling=0.1 failure=required output file missing: restraints.tsv
  [ok] G4 attack:zero_byte_files score=0.0 ceiling=0.0 failure=restraints.tsv is empty
  [ok] G4 attack:header_only score=0.0 ceiling=0.1 failure=gate failed: restraints.tsv lists every model exactly once
  [ok] G4 attack:prior_only_constant score=0.0 ceiling=0.1 failure=
  [ok] G4 attack:one_count_off score=0.869412 ceiling=0.89 failure=
  [ok] G4 attack:best_is_best_wrong_model score=0.691429 ceiling=0.72 failure=
  [ok] G4 attack:l3_reference score=0.006136 ceiling=0.75 failure=
  [ok] G4 attack:naive_straight_line score=0.241429 ceiling=0.75 failure=
  [ok] G4 attack:naive_same_protein_same_chain score=0.262605 ceiling=0.75 failure=
  [ok] G4 attack:naive_one_chain_per_protein score=0.006136 ceiling=0.75 failure=
  [ok] G4 attack:naive_pooled_counts score=0.7 ceiling=0.75 failure=
  [ok] G4 attack:naive_positions_as_resnums score=0.0 ceiling=0.75 failure=
  [ok] G4 attack:naive_single_offset score=0.428636 ceiling=0.75 failure=
  [ok] G4 attack:naive_constant_chains score=0.555882 ceiling=0.75 failure=
  [ok] G4 attack:naive_top_confidence score=0.6 ceiling=0.75 failure=
  [ok] G4 attack:naive_prune_unsatisfied score=0.55 ceiling=0.75 failure=
  [ok] G4 attack:naive_default_pipeline score=0.0 ceiling=0.6 failure=
{"task": "xlms-model-selection-s0002-l4", "ok": true, "minutes": 13.73}
exit=0
```

## What was dropped from the spec, and why

- **Q_two_states: not built.** The spec makes it conditional on keeping one answer. With two
  arrangements, which pair of models "jointly explains the most links" would have needed its own tie
  and margin analysis on top of the surface-distance margins; I did not have a design for which I
  could show a unique answer.
- **An established tool (Jwalk) as the professional software: not attempted.** I did not try to
  install it. Reason for not trying: the scored definition would then be the tool's behaviour
  (start voxels, radii, grid placement), and the second implementation would have to agree with a
  program rather than with a written definition. The definition is written out in the notes and the
  analyst implements it.
- **Models that differ inside the chains: not built.** Chains are still rigid. Under the surface
  rule the counts inside a heterodimer are nevertheless not constant any more (the other heterodimer
  can block a path or offer an assignment): the intra-protein satisfied count takes 6 / 4 distinct
  values at seed 1 / 2, the inter-protein count 4 / 4. Ranking is by the total, so pooling is now the
  stated criterion and `naive_pooled_counts` is only a reporting error (0.70).
- **Repacking the wrong poses: not built.** The shipped wrong models are far from the native anyway
  (see limitations), with similar CA contact numbers and clash counts as the near-native ones
  (CA-CA < 4 A: near-native 4, 3, 0 / 2, 9, 0; wrong 0-8 / 4-15), but no side chain was repacked.
- **Gray-zone tolerance in the count component: not needed.** The generator ships no scored call
  within 3 A of the cutoff, so counts stay exact everywhere.
- **Per-link classification** is scored on all models (three-class MCC), not only on the best model.

## Known limitations and doubts (level 4)

- **Few informative cross-links.** Only 8 of 52 (seed 1) and 13 of 59 (seed 2) rows are true links
  between the heterodimers. 9 of the 14 wrong models at seed 1 and 10 of 14 at seed 2 share one
  total (31 / 33). The near-native models lead the best wrong model by 4 at both seeds.
- **Wrong models are far wrong.** Poses with interface RMSD 6-9 A satisfy the same cross-links as
  the near-native ones under any distance definition (measured on the 1FIN pool: 18 of 19 hetero
  candidates, 4-5 of 5 same-position candidates), so the generator never ships them. Shipped wrong
  models: interface RMSD 12.2-44.2 A (seed 1), 9.5-52.4 A (seed 2).
- **The decoys are helped by planted false links.** For the straight line to rank model_01 first, 4
  (seed 1) / 3 (seed 2) false pairs are chosen that are within straight-line range in model_01 and in
  none of the near-native models, and 1 / 2 false hetero pairs that one wrong model satisfies along
  the surface. Random false
  positives would not line up behind one model like this. The three "through the protein" false
  pairs per instance are of the kind random false positives do produce. Residue-pair FDR is 19% /
  17%, higher than at level 3 (12%).
- **Small leads of the decoys.** Straight line: the wrong model leads by 1 (seed 1) / 2 (seed 2)
  cross-links. One-chain rule: 1 / 2. They are strict leads (the careless variants do put a wrong
  model first on both seeds), but one link wide at seed 1.
- **P_homo_oligomer is thinner on 1FIN.** In that crystal the two heterodimers touch through
  hetero contacts; there is no pair of two different positions of one protein that is in range only
  across the copies. The element rests on 5 same-position rows plus the assignment over copies for
  hetero links. On 1EFU (EF-Ts dimer in the middle) 9 of the 13 between-heterodimer links are
  same-protein links.
- **The 2:2 arrangement of 1FIN is the crystal's asymmetric unit**, not an assembly I can cite as
  biological. It is a real, clash-free arrangement and serves as hidden geometry only.
- **The surface distance is a lattice quantity.** 26-neighbour path lengths depend on the direction
  relative to the lattice, so the same residue pair has slightly different distances in the frames
  of different models. Tethers are not tested against the protein. One radius for all atoms. An
  XL-MS specialist may know the measure in a different form (CB atoms, element radii); here the
  notes are the definition.
- **Two implementations, one author.** The reference (scipy) and the numpy search agree to the
  printed precision, but the second solution shares the lattice module with the generator; the
  independent pair is scipy reference versus numpy generator. Agreement was measured on two seeds
  (docker) and is enforced per instance only by G2.
- **Float32 readers.** The reference takes its coordinates from Biopython (float32) and still
  reproduces every status on both seeds; lattice points that lie within about 1e-5 A of the blocking
  radius can differ for such a reader, and the 3 A margin is what absorbs it. Not proven in general.
- **Sensitivity to the parameters** (second solution with the notes edited): radius 2.9 / 3.1 A
  changes 0 / 3 of 884 statuses at seed 1; tether 6.5 / 7.5 A changes 13 / 1. A solver that uses a
  slightly different radius or tether fails the count component.
- **Compute.** Reference 2.8 / 6.6 min in the agent image while the build host was overloaded
  (1.9 / 3.0 min on the host), second solution 30 / 38 s on the host; `oracle_runtime_min` in
  `task.toml` says 2.0 / 3.0. A lattice search written with Python loops over a million points per
  source will be far slower. The agent limit is the default 2 h. I did not measure a slow
  implementation, so I do not know how much of the time limit a first naive attempt burns.
- **Reference solution = one script with 13 analysis steps** (`oracle_steps = 13`), below the
  level-4 target of 20.
- **Generator.** 2-3 min per attempt, up to 6 attempts; on 1FIN attempts are rejected often (seed 3
  needed four, 7 min on the loaded host), so some seeds may exhaust the attempts and the 10 min
  budget. Only seeds 1 and 2 went through the gates; only two natives; no host sweep beyond seeds
  1-4.
- **Both solutions read cutoff, radius and tether from the notes** with patterns tied to the
  generator's wording.
- **Not piloted.** Whether a small model still solves this in one attempt is unknown.

## Reproducibility checks (level 4)

- Generation is deterministic: after the last code edit both instances were regenerated into a
  scratch directory (`make --seed N --level 4 --out <scratch> --results <scratch> --no-validate`);
  `diff -r` against `tasks/xlms-model-selection-s0001-l4` and `-s0002-l4` shows no difference.
- Level 3 seed 1 regenerated with the final code: `diff -r` against
  `tasks/xlms-model-selection-s0001-l3` shows no difference.
