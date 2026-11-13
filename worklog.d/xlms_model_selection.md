# Worklog: xlms_model_selection

## 2026-11-11

**Reading.** AUTHORING.md, the family spec, the crispr_screen reference family and the core modules
(spec, dockerlib, harbor, validate, scoring). Core modules not touched.

**Native shortlist.** Downloaded 30 candidate heterodimers from RCSB into
`.cache/xlms_model_selection/` and surveyed chain sizes, lysine counts and the number of
inter-chain Lys-Lys pairs within 20/26/30 A (CA-CA). Small complexes (barnase/barstar class,
protease/inhibitor) were rejected: with a 30 A CA-CA cutoff almost any pose of a 100-residue ligand
satisfies many links. Kept five entries with >= 31 inter-chain K-K pairs <= 26 A, no insertion
codes, no modified residues, and a constant shift between residue numbers and SEQRES:
1FIN A/B, 1EFU A/B, 1KXP A/D, 1IBR A/B, 1I2M A/B.

**First generator.** Rigid-body poses of chain 2, links sampled from native geometry, own PDB
writer, hidden quality (iRMSD, fnat, LRMSD, DockQ) from the transform. Ran seeds 1, 2 (level 3) and
seed 1 (level 1) on the host: release conditions held (margin >= 0.6 A, near-native lead >= 3-4
links, 11-13 distinct inter counts).

**Verifier + two solutions + naive flags.** Host run of oracle (Biopython), independent solution
(stdlib + numpy, difflib) and all naive flags against score.py: oracle and independent solution 1.0
on all three instances; naive variants 0.0-0.65.

**Found: wrong decoys were recognisable without any cross-link.** The first decoys were slid in
until the first CA-CA contact (min CA-CA 4.6-6 A). Measured on seed 1: CA-CA pairs < 8 A across the
interface were 58/83/124 for the three near-native models and 0-34 for every wrong model. Even the
garbage numbering variant then ranked near-native models high (AUC 0.90/0.92) purely because they
are the compact ones. Changed the slide: protein 2 is pushed in until the pose has 0.6-1.6 x the
native number of CA-CA contacts < 8 A (closed form: k-th largest "reach" along the slide axis).
After the change, seed 1: wrong models 51-108 contacts (near-native 58-124), CA-CA pairs < 4 A 0-20
(near-native 0-10). Also fixed a slip found at the same time (the slide used residues of the
trimmed loops, so two poses had no contact at all).

**Notes wording.** "one flexible internal loop per protein" was not true for natives that have
additional unresolved segments (1EFU, 1KXP, 1I2M, ...): reworded to "disordered terminal regions and
flexible internal loops had been trimmed". FDR line reworded to a CSM-level cut with the estimated
residue-pair FDR (the planted rate).

**Robustness sweep.** Seeds 1-20 at levels 2 and 3 on the host (40 instances): 39 generated, one
failed (`seed 10, level 2`: fewer than two cross-linked positions missing from the models). Cause:
below level 3 the only links to unmodelled lysines came from the trimmed loops. Fix: at every level
one reported pair reaches into a region without coordinates (N-terminal region / segment disordered
in the native). Seed 10 level 2 then generates (2 unresolved positions) and passes the host battery.

**Verifier fuzzing on the host (19 malformed submissions).** No verifier crash; all contract
violations give reward 0. Finding: with the spec's weights (0.35/0.35/0.20/0.10) a submission whose
`mapping.tsv` is useless (header only, all NA, garbage residue numbers) scored
0.35 + 0.35 + 0.20 = 0.8999999999999999 and was reported as failed only because of floating point
rounding. Changed the weights to 0.30 (counts) / 0.35 / 0.20 / 0.15 (mapping) and made the mapping
credit "1.0 only if every position is right, else 0.6 x fraction" (same rule as the counts). The
missing-mapping case is now 0.85; the level-1 constant answer (identity mapping) is bounded by 0.09.

**Notes: conventions restated.** Added two lines to `xlms_notes.md` (same-protein pairs are
intra-protein because the sample is a 1:1 heterodimer; a cross-link to an absent residue is left out
of satisfied count and total) so that everything the exact counts depend on is in the notes as well
as in the instruction.

**Independent label check.** Biopython parser + SVDSuperimposer on CA atoms from the shipped model
files against the cached native: ligand RMSD within 0.58 A (seed 1) / 0.53 A (seed 2) of the stored
value for all models; near-native 1.3-4.0 A, wrong >= 17.5 A.

**Final sweep.** Seeds 1-20 x levels 1-3 on the final generator: 60/60 generated (numbers in the
design notes).

**First Docker validation (seed 1, level 3, before the two wording changes and the fix above).**
All gates passed including Harbor oracle = 1.0 / nop = 0.0, 3.47 min. Harbor ran without
`ALEGEN_NO_EGRESS_SIDECAR`.

**Final runs.** The three `make` commands were rerun on the final code (all gates ok; seed 1 level 3
including Harbor oracle = 1.0 / nop = 0.0); logs in `results/xlms_s1_l3.log`,
`results/xlms_s2_l3.log`, `results/xlms_s1_l1.log` (pasted in the design notes). For seed 2 (and
level 1) the command was run once more with `--keep-work` *before* the exact command, only to read
the raw metric values out of `details.json`. Extra run `--seed 1 --level 2 --no-harbor`: all gates
ok (`results/xlms_s1_l2.log`). A run of the three commands that was in progress when the verifier
weights changed was stopped and discarded.

**Disk hygiene (coordinator request).** No tools image was built (generation needs numpy only).
Own task images are removed with `docker rmi` as soon as the next instance is validated; the Biopython
scratch venv used for host tests lived in the session scratchpad and was deleted at the end;
`results/gen/xlms_model_selection-*` deleted after the task directories were written.

**Not done / cut.**
- No flexible regions: models are strict rigid-body variants, so intra-protein counts are identical
  in all models and pooling inter+intra cannot change the order. The spec's "intra+inter pooled"
  variant is therefore only a reporting error (fails the count component).
- No naive variant makes the false-positive trap change the *selection*; the variant that discards
  unsatisfied links fails the count component only.
- No coordinate noise inside the chains; no side-chain repacking of decoys.
- Cache trimmed to the five entries the generator uses (25 surveyed candidates deleted).

## 2026-11-12 (level 4)

**Reading.** AUTHORING "Level 4", `docs/specs/l4/xlms_model_selection.md`, README section 6, the level-3
code and notes. Core modules not touched. Decision: implement P_surface_distance and P_homo_oligomer,
drop Q_two_states (the spec makes it conditional on keeping one answer; no time to show that).

**Is there a tool to ship?** Not attempted. A surface distance from an existing program would make that
program's undocumented details (start voxels, radii, how the ends are attached) the scored definition,
and I would then need a second independent implementation that agrees with it on every scored link.
I wrote the definition myself instead (integer lattice in the frame of the file, one blocking radius,
26 neighbours, straight tethers of bounded length at both ends, minimum over everything) and
implemented it twice.

**Engine 1 (numpy, `assets/l4/sasd_np.py`).** Bucketed search (Dial): because every step is >= 1 A,
all points with a tentative length in [k, k+1) are final and can be relaxed together with array
operations. Checked against a heap Dijkstra written for the test on a random 260-atom system: 132
anchor pairs, 0 mismatches, asymmetry 7e-15. First speed measurement on the 1FIN tetramer (8994
atoms, 1.3 M lattice points): 0.14 s per source with a 60 A limit.

**Found: distances beyond the search limit are not symmetric.** With a limit, D[i,j] and D[j,i]
differed by up to 0.34 A when the total was above the limit (the best last point of one direction was
pruned). Fixed by reporting everything above the limit as infinite.

**Numbers that decided the parameters (1FIN, chains A-D, all 80 lysine CA atoms).** Blocking radius
3.0 A: every lysine has an open lattice point within 5.85 A (median 3.5 A); with 3.2 A two lysines
need more than 7 A. Chose radius 3.0 A and tether 7.0 A. In the same exploratory run (radius 3.2 A,
tether 8 A, i.e. before the parameters were fixed) 1051 lysine pairs were within straight-line range at
a 34 A cutoff and 498 within surface range; 287 pairs were >= 3 A inside by straight line and >= 3 A
outside along the surface; median surface/straight ratio 1.28. Not re-measured with the final values.

**Natives for a complex with copies.** Surveyed the five cached entries for a second heterodimer.
1FIN: A-D and B-C touch (17 CA-CA contacts < 8 A each). 1EFU: B-D 42 contacts, A-C 26 (the EF-Ts
dimer). 1I2M: thin (14 between-copy candidate pairs in total). 1IBR: the heterodimers do not touch.
1KXP: one heterodimer only. Level 4 uses 1FIN and 1EFU, heterodimer 1 fixed, heterodimer 2 moved.

**First attempt failed: links first, poses second (the level-3 order).** With a 3 A margin on both
sides of the cutoff, 22 native between-copy candidate pairs shrank to 3 once 16 poses were fixed, and
no same-position pair survived. Changed the order: a pool of 6 near-native and 36 wrong poses gets a
full table of surface distances (every lysine CA to every other, both copies, final frame, rounded
coordinates); models and links are then chosen together from the tables.

**Found: wrong poses with interface RMSD 6-9 A satisfy the same links as the near-native ones.**
(1FIN pool: 18 of 19 hetero between-copy candidates and 4-5 of 5 same-position candidates.) They
cannot be told apart by this data at any distance definition, so they are never shipped: every wrong
model is at least 3 satisfied links behind. Consequence: the shipped wrong models are far away
(seed 1: interface RMSD >= 12.2 A).

**Selection.** Version 1 (pick the links for one decoy pose, then look for further wrong poses in
which all of them are clean) found nothing for any of the three pools: too few compatible poses
(typically 5-11 of the 12-13 needed). Version 2 drops the between-copy links that sit in the gray zone
of many pool poses: 1EFU and 1I2M then gave instances, 1FIN still not (11-12 of 13). Version 3: pool
of 36 instead of 30 wrong poses, 15-18 models instead of 16-19. The 1FIN and 1EFU pools then gave an
instance for each of three random draws (1I2M was not rebuilt).

**Found: straight-line and surface statuses were nearly identical.** First complete instance (seed 1):
the straight-line variant had 13 wrong statuses of 748 and 11 of 17 models exact (score 0.44, status
MCC 0.968, inside the credit band). Cause: my false links were all chosen far apart by both measures.
Added three false pairs per instance that are within straight-line range in every model and >= 3 A
out of range along the surface (the path would go through the protein). After the change: 74 / 83
wrong statuses (seed 1 / 2), 0 of 17 models exact, MCC 0.84 / 0.85.

**Weights.** With 0.25 (counts) the pooled-counts variant scored exactly 0.75 = its ceiling. Moved to
0.30 counts / 0.15 statuses / 0.30 best / 0.10 ranking / 0.15 mapping: pooled 0.70, one count off
0.869, wrong best model 0.691.

**Host replay of the battery (both seeds, final code).** Reference solution in a scratch venv
(scipy 1.18.1, Biopython 1.85), second solution and all careless variants with the task's scorer:
reference 1.000 / 1.000, second solution 1.000 / 1.000, every careless variant <= 0.70 (table in the
design notes). Reference runtime on the host 115 s / 178 s, second solution 30 s / 38 s.

**Margins under both implementations (search limit raised to cutoff + 12 A).** Seed 1: 816 scored
calls, satisfied side min 3.433 A, violated side min 3.003 A, no call on the wrong side in either
implementation; 629 calls finite in both, largest difference 0.001 A (rounding of the printed value).
Seed 2 in the design notes.

**Sensitivity (seed 1, second solution, notes edited).** Blocking radius 2.9 A: 0 of 884 statuses
change; 3.1 A: 3; tether 6.5 A: 13; tether 7.5 A: 1.

**Levels 1-3 unchanged.** `make ... --seed 1 --level 3 --out <scratch> --results <scratch> --no-validate`
and `diff -r` against `tasks/xlms-model-selection-s0001-l3`: identical. Levels 1 and 2 generate.
The only edit to the level-3 module is a three-line dispatch at the top of `generate`.

**Disk.** Free space fell from 24 GB (start) to 7.4-8.0 GB while the host battery was running; no
image was built below 8 GB. Seed 1 `make` was started at 8.7 GB free (01:47 UTC); when its build step
began the disk showed 6.3 GB (another author's build was running), by then my agent image was already
built and only the verifier image (no new layers) followed, so I let it finish. Seed 2 `make` was
started at 8.6 GB and its build ran at 11 GB free.

**Docker validation, seed 1** (`results/xlms_s1_l4.log`): all gates ok including Harbor oracle = 1.0
and nop = 0.0; 18.1 min; reference solution 2.81 min in the agent image (load average 10-14 on the
host). Scores identical to the host replay. Task directory and truth are byte-identical to the host
generation.

**Docker validation, seed 2** (`results/xlms_s2_l4.log`, `--no-harbor`): all gates ok; 13.7 min;
reference solution 6.59 min (my own low-priority seed sweep and other authors' jobs were running).
Seed-2 images removed afterwards.

**Generator robustness (host, not Docker).** Seed 3 (1FIN): attempts 0-2 rejected (too few compatible
wrong poses; a total shared by 11 models; no pose that violates >= 5 same-protein between-copy links
while satisfying a false hetero pair), attempt 3 accepted, 431 s in total. Seed 4 (1EFU): first
attempt, 131 s. 1FIN is the weak native for this design: its heterodimers touch through hetero
contacts only.

**After both instances had been generated for their validation runs** (seed 2's battery was still
running) I removed unused imports from `l4.py` and corrected its docstring (it
still said the scipy implementation runs at generation time, which I had dropped: the generator uses
the numpy search only and G2 is the cross-check). Regenerated both seeds with the final file into a
scratch directory and compared with `diff -r` against `tasks/` (result in the design notes / report).

**Not done / cut at level 4.**
- Q_two_states; Jwalk or any other shipped tool for the distance; models that differ inside the
  chains; repacked wrong poses.
- No run of the scipy implementation inside the generator (it would need a tools image at generation
  time); G2 is the per-instance check that the two implementations agree.
- The residue-pair FDR ended at 17-19% because the decoy constructions need 4-6 planted false pairs;
  I did not rebalance it.
- No measurement of a slow (pure Python) lattice search, no pilot.
- Clean-up at the end: scratch venv (scipy/Biopython, 233 MB) and all generated scratch data deleted
  (the measurement scripts and small logs stay in the session scratch directory; the logs quoted in
  the design notes are copied to `results/xlms_l4_host_measurements.txt`); `results/gen` of the two
  level-4 instances deleted; `results/work/xlms-model-selection-s000{1,2}-l4` kept (`--keep-work`);
  seed-1 images kept, seed-2 images removed.

**Regeneration check with the final files.** Both level-4 instances regenerated into scratch are
identical to `tasks/xlms-model-selection-s0001-l4` and `-s0002-l4` (`diff -r`); level 3 seed 1 is
identical to `tasks/xlms-model-selection-s0001-l3`.
