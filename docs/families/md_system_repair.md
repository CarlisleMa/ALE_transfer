# Family `md_system_repair` - from a defective crystal structure to a simulation-ready Amber system

ALE subdomain 3.1 (Biomolecular Structure & Design). Code: `pipeline/alegen/families/md_system_repair/`.

## Scenario

A computational chemist receives a refined crystallographic model (`structure.pdb`) and the group's
SOP (`simulation_sop.md`) and has to hand over an Amber topology, the as-built coordinates, a
briefly equilibrated restart and a small JSON report. The SOP defines the target system completely
(what is simulated, ff14SB/TIP3P, rectangular cell with a stated padding, neutralisation plus NaCl
by an explicit formula, a three-stage equilibration) but is written as a general group document:
it does not say which of its clauses matter for this file. The file is what a crystallographer
exports: REMARK 465/470, SEQRES, MODRES, HETNAM, alternate locations, insertion codes, hetero groups.

Tools in the agent image (conda, pinned): `ambertools=25.3` (tleap, pdb4amber, cpptraj, sander,
ParmEd 4.3.1), `openmm=8.6.1`, `pdbfixer=1.12`, python 3.12 with numpy/scipy/pandas.
Verifier image: `openmm=8.6.1`, `parmed=4.3.1`, `scipy=1.18.1` (ParmEd's NetCDF reader), numpy.

## Truth model

The hidden truth is a complete, real protein structure; the defects are planted by the generator,
so every scored fact is known before any tool runs.

| source | URL | sha256 | used for |
|---|---|---|---|
| 1AKI (hen egg-white lysozyme, 129 aa, 4 disulfides, 2 Met) | https://files.rcsb.org/download/1AKI.pdb | `c75d7a689617248cdd92dc6633531d2506fb9bef1e6e21e26c8f579ae6955abb` | even seeds |
| 4PTI (BPTI, 58 aa, 3 disulfides, 1 Met) | https://files.rcsb.org/download/4PTI.pdb | `cd9b0faf2e66c0e031f30d1c84cb359691b8f5819f1b0f784be15f6005d58e08` | odd seeds |

Cached under `.cache/md_system_repair/`; the generator re-checks the sha256 on every run. Both
entries are complete (asserted against a heavy-atom table, OXT included), free of their own
alternate locations and of hydrogens. The deposited SSBOND pairs are cross-checked against SG-SG
distances (< 2.5 Å) at generation time, so the disulfide truth rests on two independent routes.
The entry id is blinded in the file (`HEADER ... XXXX`); the sequence itself is visible (SEQRES).

Planted per `(seed, level)` (`plant.py`):

| defect | how it is planted | what is true by construction |
|---|---|---|
| chain break (L3) | 4-6 consecutive internal residues deleted (no Cys/Met/Pro in the segment or its flanks; C-N distance of the flanks > 6 Å); listed in REMARK 465 | expected residue list, two fragments, caps at the break, net charge of the remaining residues |
| alternate conformers (L3) | 3 residues get a decoy side-chain conformer (side chain rotated +-120 degrees about CA-CB, accepted only if >= 3.0 Å from all other residues); occupancies 0.58-0.72 for the true one; for at least one residue the decoy is listed first as altloc A | true chi1 (the crystal conformer) and its coordinates |
| numbering (L3) | offset 17/100/201/300, two insertion-coded residues placed so that a disulfide cysteine carries a code, chain renamed | labels of every residue as written in the file |
| modified residue (L2+) | one Met written as `HETATM ... MSE` with `SE`, MODRES and SEQRES consistent | the residue is MET with atom SD |
| missing side chains (L2+) | 3-4 residues truncated after CB or CG, REMARK 470 written | full heavy-atom set of each residue |
| additives and waters (L2+) | glycerol, sulfate, sulfate/ethanediol, one Cl-, one Na+ placed at crystal-water sites (>= 3.1 Å from the protein); remaining crystal waters of the entry kept | positions of the crystal water and ion sites (55 for seed 1, 70 for seed 2) |
| hidden disulfides (L2+) | no SSBOND/CONECT records (level 1 keeps SSBOND) | true pairs |

SOP parameters vary with the seed: padding 8.0 or 9.0 Å, NaCl 0.10 or 0.15 M.

Level 1 is the long workflow on a clean file (complete protein, SSBOND present, original numbering,
no hetero groups). Level 2 adds hidden disulfides, additives/waters, the modified residue and
missing side chains. Level 3 adds the chain break, alternate conformers and the numbering scheme.

## Decision points and why each is discoverable

| trap | naive outcome | discoverable by |
|---|---|---|
| hidden_disulfides | tleap builds reduced CYS with clashing thiol hydrogens | SG-SG distances ~2.0 Å; `pdb4amber` renames to CYX and writes an sslink file; SOP 1.7 |
| additives_and_waters | waters and Na+/Cl- load into tleap without any message and become part of the system | HETATM/HETNAM records; tleap "Unknown residue: GOL/SO4/EDO" and `saveamberparm` "FATAL: Atom ... does not have a type"; SOP 1.8 |
| modified_residue | `grep ^ATOM` (or deleting what tleap rejects) removes a residue from the chain | MODRES/SEQRES, HETATM lines inside the chain, tleap "Unknown residue: MSE"; SOP 1.6 |
| missing_sidechains | stub kept as deposited (ALA/GLY) | REMARK 470; `pdb4amber` "misses n heavy atom(s)"; tleap "Added missing heavy atom"; SOP 1.4 |
| chain_break | no TER: tleap bonds across 13.7 Å (it only warns); `pdb4amber` TER: charged termini inside the protein | REMARK 465, SEQRES vs ATOM, C-N distance, `pdb4amber` "gap of ... A", tleap "There is a bond of 13.705 angstroms between C and N atoms"; SOP 1.1-1.3 |
| altlocs | `pdb4amber`, tleap and OpenMM keep the first listed conformer | altLoc column, occupancies < 1.00, `pdb4amber` "alternate locations" list, tleap "Atom names in each residue should be unique"; SOP 1.5 |
| numbering_insertion_codes | disulfides reported in the builder's 1..N numbering | column 27, `pdb4amber` renum file; output contract of `report.json` |

Measured tool behaviour that the traps rely on (AmberTools 25.3, on the seed-1 level-3 file):
`pdb4amber` keeps the first conformer, leaves MSE as HETATM, keeps waters and additives, puts a TER
at the gap and adds spurious TER records around the insertion-coded residues. tleap on the raw
file keeps the first of duplicate atoms, creates untyped units for MSE/GOL/SO4/EDO and then refuses
to save ("FATAL: Atom .R<MSE 65>.A<N 1> does not have a type"); with the HETATM records deleted it
saves a topology with "Errors = 0; Warnings = 4", two of the warnings being bonds of 13.705 Å (the
unresolved segment) and 3.291 Å (the deleted MSE); with only waters, Cl- and Na+ kept it saves
without any message about them. tleap rebuilds missing side-chain atoms ("Added missing heavy
atom"). `missing_sidechains` is therefore not a trap for a tleap or PDBFixer user
(the tools repair it by default); it is a defect whose repair the verifier confirms, and its naive
variant is the stub-preserving workflow.

## Output contract

`system.prmtop`; `system.inpcrd` (as built, with box; ASCII or NetCDF); `equilibrated.rst7` (last
equilibration frame with velocities and box; ASCII or NetCDF); `report.json` with `disulfides`
(pairs of residue labels as written in `structure.pdb`, insertion code included),
`protein_net_charge`, `n_na`, `n_cl`, `n_water`, `box_angstrom` (of the restart).

## Scoring

Pass threshold 0.90. Gates (score 0): files parse; atom counts agree; periodic box present; >= 70%
of the expected residues found in order; as-built CA RMSD to the hidden structure <= 1.0 Å.
Each component is the mean of its checks. `production_ready` (0.30) is 1 only when every check of
every other component passed, so one ignored requirement caps the score at 0.70 while partial work
keeps its per-requirement credit (0.70 in total).

| weight | component | checks |
|---|---|---|
| 0.07 | sequence | residue identity vs the resolved sequence (ramp 0.90-1.00); nothing missing/extra; no foreign residue |
| 0.07 | atoms_and_charge | fraction of residues with exactly the expected heavy atoms (ramp 0.90-1.00); solute charge = SOP value +-0.05 e |
| 0.08 | disulfides | Jaccard of S-S bonds vs true pairs; bonded SG carrying H counts as false |
| 0.08 | termini_and_breaks | N-/C-terminus charged; peptide bonds inside fragments (ramp 0.9-1.0); no other inter-residue bond; N-methylamide and acetyl at the break |
| 0.06 | conformers | fraction of resolved heavy atoms within 1.0 Å of the hidden structure (ramp 0.90-0.99); chi1 within 40 degrees of the true conformer for each planted residue |
| 0.10 | solvent_and_ions | neutral +-0.05 e; only Na+/Cl-; salt by the SOP formula +-0.02 M and |nCl-nNa| = |solute charge|; TIP3P charges; padding in [P-0.5, P+5.5] Å on each axis, orthorhombic; <= 20% of crystal water/ion sites occupied within 0.5 Å |
| 0.07 | physical_sanity | longest solute heavy-atom bond as built <= 2.6 Å; OpenMM builds the PME system; 60 L-BFGS iterations give a finite negative energy |
| 0.08 | equilibrated | restart with orthorhombic box; density 0.97-1.08 g/cm3; kinetic temperature of stored velocities 250-350 K (degrees of freedom counted with constraints when the stored velocities show constrained X-H bonds); CA RMSD <= 3.0 Å; 0.3 ps verifier MD from the restart finite, E < 0, 280-320 K, no atom moved >= 8 Å |
| 0.09 | report | six items, each equal to the recomputed value |
| 0.30 | production_ready | all checks passed |

## Why these tolerances (measured)

All numbers below were read from the verifier's `details.json` of the runs listed in the tables
(seed 1 = BPTI, padding 8 Å; seed 2 = lysozyme, padding 9 Å; level 1 and 3) plus the smoke tests listed
further down. "Naive" is the corresponding trap-ignoring build.

| check | tolerance | correct submissions (oracle, 2nd solution, truth build) | naive / wrong | comment |
|---|---|---|---|---|
| as-built CA RMSD (gate) | <= 1.0 Å | 0.000 | unrelated peptide: fails the 70%-of-residues gate before RMSD | both builders keep deposited coordinates |
| resolved heavy atoms within 1.0 Å | >= 99% | 100% | seed 1: truncated stubs 96.2%, dropped residue 98.0%, first conformers 97.5%, raw 95.4%; seed 2: 97.9%, 99.2% (passes), 99.7% (passes), 98.8% | generic guard only; the specific checks (sequence, heavy atoms, chi1) carry the decision in every one of these cases |
| chi1 of planted residues | within 40 degrees | 0.0 degrees | 120.0 degrees for every decoy | decoys are built 120 degrees away |
| S-S bonds | Jaccard = 1 | 3/3 (seed 1), 4/4 (seed 2) | 0/3, 0/4 | |
| solute net charge | +-0.05 e | exact integer | - | |
| salt concentration (SOP formula) | +-0.02 M | 0.1380 - 0.1512 (n=9) M for SOP 0.15 M | no ions 0.000; pairs removed 0.1109 - 0.1167 (n=3) | one NaCl pair is 0.016 M at 3500 waters, 0.011 M at 5200 |
| padding per axis | [P-0.5, P+5.5] Å | tleap `solvatebox` with buffer P: P+2.7..P+3.2; Modeller with explicit box P+0.25 | none measured below the band; with crystal waters kept the solute extent grows and the padding reads up to P+4.3 (still inside) | the SOP itself states "at least P, at most 5 Å excess" |
| crystal water/ion sites occupied (0.5 Å) | <= 20% | 0.0% | 98.2 - 100.0 (n=4) % | chance occupancy of a 0.5 Å sphere by bulk water is 1.7% per site |
| longest solute heavy-atom bond as built | <= 2.6 Å | 2.02-2.09 Å (S-S) | 13.7 Å (bonded across the break) | |
| minimised energy | finite, < 0 | -221288 - -131687 (n=9) kJ/mol | (not a discriminator; kept as a sanity floor) | 60 L-BFGS iterations |
| density of the restart | 0.97-1.08 g/cm3 | 0.9984 - 1.0295 (n=30) (every run whose restart passed, naive variants included) | tleap box without NPT 0.759-0.791; tleap box after 1 ps NPT 0.845-0.850; Modeller box after 30 minimiser iterations and 20 MD steps 0.939-0.947 | lower edge halfway between 0.947 and 0.998; nothing wrong was observed above the band, the upper edge leaves room for protein-rich cells |
| kinetic temperature of stored velocities | 250-350 K | 296.7 - 307.5 (n=30) K | 416-430 K after 20 MD steps (lysozyme), but 297.8 K in the same situation for BPTI; no velocities in a copied inpcrd; fabricated Maxwell-Boltzmann velocities pass (298.9 K) | X-H velocity along the bond 0.8-1.2 Å/ps for constrained runs, threshold 5 (thermal value ~16) |
| CA RMSD of the restart | <= 3.0 Å | 0.59 - 1.84 (n=30) Å | 1.7-2.5 Å after minimising a bond drawn across the break (passes this check, fails others) | generous: guards against unfolding/explosion only |
| temperature after 0.3 ps verifier MD from the restart | 280-320 K | 298.1 - 308.1 (n=30) K | unrelaxed tleap coordinates 342.7-389.1 K; Modeller build after 20 MD steps 330.9-357.6 K; minimised only with fabricated velocities 234.9 K | relaxed-but-cold and unrelaxed-hot configurations both leave the band; sigma of T for 12k atoms is ~2 K |
| reported box | +-0.02 Å | exact / rounded to 0.01 Å | - | restart files carry 7 decimals |

## Measured results

Raw values recomputed by the final verifier for every submission of the battery (`results/work/<task_id>/*/logs/details.json`). Gate lines are verbatim. For seed 1 level 3 the first block is the complete `make` with the verifier before the last band change and the second block is the re-validation with the final verifier (script outputs of the first block replayed, see limitation 0); the table shows the second block. Seed 2 level 3 and seed 1 level 1 are from-scratch `make` runs with the final verifier.

### Seed 1, level 3 (4PTI, padding 8.0 Å, 0.15 M)

```
# (1) ../.venv/bin/python -m alegen make md_system_repair --seed 1 --level 3 --keep-work
#     2026-11-11 16:26-17:13 UTC, complete run, verifier BEFORE the final band placement (density 0.95-1.06, verifier-MD 260-340 K, atoms in place 97%)
wrote /home/user/ale-lifesci-gym/tasks/md-system-repair-s0001-l3  (prior difficulty 0.675)
  [ok] G1 images_build 
  [ok] G2 oracle score=1.0 runtime_min=3.17 failure= stderr=
  [ok] G5 verifier_deterministic first=1.0 second=1.0
  [ok] G3 probe:within_tolerance score=1.0 floor=0.99 failure=
  [ok] G3 probe:independent_solution score=1.0 floor=0.99 failure=
  [ok] G4 attack:empty_output score=0.0 ceiling=0.0 failure=required output file missing: system.prmtop
  [ok] G4 attack:inputs_copied score=0.0 ceiling=0.1 failure=required output file missing: system.prmtop
  [ok] G4 attack:zero_byte_files score=0.0 ceiling=0.0 failure=system.prmtop is empty
  [ok] G4 attack:header_only score=0.0 ceiling=0.1 failure=report.json is not valid JSON: Expecting property name enclosed in double quotes: line 2 column 1 (char 2)
  [ok] G4 attack:prior_only_generic_system score=0.0 ceiling=0.05 failure=gate failed: topology contains a polypeptide
  [ok] G4 attack:salt_beyond_tolerance score=0.683333 ceiling=0.75 failure=
  [ok] G4 attack:naive_reduced_cys score=0.62 ceiling=0.75 failure=
  [ok] G4 attack:naive_no_ions score=0.666667 ceiling=0.75 failure=
  [ok] G4 attack:naive_crystal_solvent_kept score=0.683333 ceiling=0.75 failure=
  [ok] G4 attack:naive_mse_dropped score=0.660736 ceiling=0.75 failure=
  [ok] G4 attack:naive_sidechains_truncated score=0.62973 ceiling=0.75 failure=
  [ok] G4 attack:naive_gap_bonded score=0.620667 ceiling=0.75 failure=
  [ok] G4 attack:naive_gap_charged_termini score=0.666852 ceiling=0.75 failure=
  [ok] G4 attack:naive_altloc_first score=0.67 ceiling=0.75 failure=
  [ok] G4 attack:naive_internal_numbering score=0.685 ceiling=0.75 failure=
  [ok] G4 attack:naive_raw_tleap score=0.419733 ceiling=0.6 failure=
  [ok] G6 harbor_oracle rewards=[1.0] tail=
  [ok] G6 harbor_nop rewards=[0.0] tail=
{"task": "md-system-repair-s0001-l3", "ok": true, "minutes": 47.15}

# (2) ../.venv/bin/python ../results/md_system_repair_s1_l3_replay.py
#     2026-11-11 18:17-18:26 UTC, FINAL verifier: validate() with all gates; the MD outputs of run (1) are replayed (no solution/probe/attack
#     script changed in between); images, every verifier run, the generic attacks and the Harbor oracle/nop runs are fresh
  [ok] G1 images_build 
  [ok] G2 oracle score=1.0 runtime_min=0.0 failure= stderr=
  [ok] G5 verifier_deterministic first=1.0 second=1.0
  [ok] G3 probe:within_tolerance score=1.0 floor=0.99 failure=
  [ok] G3 probe:independent_solution score=1.0 floor=0.99 failure=
  [ok] G4 attack:empty_output score=0.0 ceiling=0.0 failure=required output file missing: system.prmtop
  [ok] G4 attack:inputs_copied score=0.0 ceiling=0.1 failure=required output file missing: system.prmtop
  [ok] G4 attack:zero_byte_files score=0.0 ceiling=0.0 failure=system.prmtop is empty
  [ok] G4 attack:header_only score=0.0 ceiling=0.1 failure=report.json is not valid JSON: Expecting property name enclosed in double quotes: line 2 column 1 (char 2)
  [ok] G4 attack:prior_only_generic_system score=0.0 ceiling=0.05 failure=gate failed: topology contains a polypeptide
  [ok] G4 attack:salt_beyond_tolerance score=0.683333 ceiling=0.75 failure=
  [ok] G4 attack:naive_reduced_cys score=0.62 ceiling=0.75 failure=
  [ok] G4 attack:naive_no_ions score=0.666667 ceiling=0.75 failure=
  [ok] G4 attack:naive_crystal_solvent_kept score=0.683333 ceiling=0.75 failure=
  [ok] G4 attack:naive_mse_dropped score=0.659027 ceiling=0.75 failure=
  [ok] G4 attack:naive_sidechains_truncated score=0.626776 ceiling=0.75 failure=
  [ok] G4 attack:naive_gap_bonded score=0.620667 ceiling=0.75 failure=
  [ok] G4 attack:naive_gap_charged_termini score=0.666852 ceiling=0.75 failure=
  [ok] G4 attack:naive_altloc_first score=0.667447 ceiling=0.75 failure=
  [ok] G4 attack:naive_internal_numbering score=0.685 ceiling=0.75 failure=
  [ok] G4 attack:naive_raw_tleap score=0.417141 ceiling=0.6 failure=
  [ok] G6 harbor_oracle rewards=[1.0] tail=
  [ok] G6 harbor_nop rewards=[0.0] tail=
{'task': 'md-system-repair-s0001-l3', 'ok': True, 'minutes': 8.85, 'mode': 'replayed MD outputs, fresh verifier + Harbor'}
```

| submission | score | failed checks | S-S | chi1 delta (deg) | atoms in place | salt (M) | padding x,y,z (Å) | xtal sites kept | max bond (Å) | density | T stored vel. (K) | CA RMSD eq. (Å) | T after verifier MD (K) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| oracle | 1.0 | - | 3/3 true, 0 false | 0.0, 0.0, 0.0 | 1.0 | 0.1426 | 11.02, 10.84, 10.78 | 0.0 | 2.088 | 1.0097 | 297.3 | 1.195 | 298.1 |
| independent_solution | 1.0 | - | 3/3 true, 0 false | 0.0, 0.0, 0.0 | 1.0 | 0.1422 | 8.25, 8.25, 8.25 | 0.0 | 2.088 | 1.0025 | 298.5 | 1.243 | 307.4 |
| within_tolerance | 1.0 | - | 3/3 true, 0 false | 0.0, 0.0, 0.0 | 1.0 | 0.1426 | 11.02, 10.84, 10.78 | 0.0 | 2.088 | 1.015 | 302.5 | 1.166 | 301.9 |
| naive_altloc_first | 0.6674 | conformers 0.457: resolved_atoms_in_place, major_conformer_32, major_conformer_34 | 3/3 true, 0 false | 120.0, 120.0, 0.0 | 0.9747 | 0.1518 | 10.83, 10.97, 10.78 | 0.0 | 2.088 | 1.0168 | 301.3 | 1.118 | 298.2 |
| naive_crystal_solvent_kept | 0.6833 | solvent_and_ions 0.833: no_crystal_solvent_or_ions | 3/3 true, 0 false | 0.0, 0.0, 0.0 | 1.0 | 0.1623 | 11.29, 11.2, 11.96 | 0.9818 | 2.088 | 1.0143 | 302.9 | 1.1 | 299.9 |
| naive_gap_bonded | 0.6207 | termini_and_breaks 0.5: no_bond_between_unlinked_residues, break1_c_side_n_methylamide, break1_n_side_acetyl; physical_sanity 0.667: as_built_bond_lengths; equilibrated 0.8: density | 3/3 true, 0 false | 0.0, 0.0, 0.0 | 1.0 | 0.1426 | 11.02, 10.84, 10.78 | 0.0 | 13.705 | 0.8502 | 289.4 | 2.505 | 293.8 |
| naive_gap_charged_termini | 0.6669 | atoms_and_charge 0.907: heavy_atoms_complete; termini_and_breaks 0.667: break1_c_side_n_methylamide, break1_n_side_acetyl | 3/3 true, 0 false | 0.0, 0.0, 0.0 | 1.0 | 0.1426 | 11.02, 10.84, 10.78 | 0.0 | 2.088 | 1.0084 | 306.5 | 1.842 | 301.7 |
| naive_internal_numbering | 0.685 | report 0.833: disulfides | 3/3 true, 0 false | 0.0, 0.0, 0.0 | 1.0 | 0.1426 | 11.02, 10.84, 10.78 | 0.0 | 2.088 | 1.0083 | 300.7 | 1.118 | 300.4 |
| naive_mse_dropped | 0.659 | sequence 0.605: residue_identity, no_extra_or_missing_residues; atoms_and_charge 0.907: heavy_atoms_complete; termini_and_breaks 0.936: peptide_bonds; conformers 0.972: resolved_atoms_in_place | 3/3 true, 0 false | 0.0, 0.0, 0.0 | 0.9797 | 0.1425 | 11.02, 10.84, 10.78 | 0.0 | 2.088 | 1.0139 | 298.6 | 0.966 | 301.3 |
| naive_no_ions | 0.6667 | solvent_and_ions 0.667: neutral, salt_concentration | 3/3 true, 0 false | 0.0, 0.0, 0.0 | 1.0 | 0.0 | 11.02, 10.84, 10.78 | 0.0 | 2.088 | 0.9984 | 298.9 | 0.996 | 304.0 |
| naive_raw_tleap | 0.4171 | sequence 0.605: residue_identity, no_extra_or_missing_residues; atoms_and_charge 0.907: heavy_atoms_complete; disulfides 0.0: ss_bonds; termini_and_breaks 0.436: peptide_bonds, no_bond_between_unlinked_residues, break1_c_side_n_methylamide, break1_n_side_acetyl; conformers 0.401: resolved_atoms_in_place, major_conformer_32, major_conformer_34; solvent_and_ions 0.667: salt_concentration, no_crystal_solvent_or_ions; physical_sanity 0.667: as_built_bond_lengths; equilibrated 0.8: density; report 0.833: disulfides | 0/3 true, 0 false | 120.0, 120.0, 0.0 | 0.9544 | 0.0146 | 11.26, 12.33, 11.96 | 1.0 | 13.705 | 0.8452 | 293.1 | 2.284 | 297.6 |
| naive_reduced_cys | 0.62 | disulfides 0.0: ss_bonds | 0/3 true, 0 false | 0.0, 0.0, 0.0 | 1.0 | 0.1426 | 11.02, 10.84, 10.78 | 0.0 | 1.815 | 1.0076 | 307.2 | 1.261 | 300.2 |
| naive_sidechains_truncated | 0.6268 | sequence 0.481: residue_identity, no_extra_or_missing_residues; atoms_and_charge 0.722: heavy_atoms_complete; termini_and_breaks 0.84: peptide_bonds; conformers 0.922: resolved_atoms_in_place | 3/3 true, 0 false | 0.0, 0.0, 0.0 | 0.962 | 0.1425 | 11.02, 10.84, 10.78 | 0.0 | 2.088 | 1.0084 | 299.0 | 1.255 | 304.3 |
| prior_only_generic_system | 0.0 | gate failed: topology contains a polypeptide | - | - | - | - | - | - | - | - | - | - | - |
| salt_beyond_tolerance | 0.6833 | solvent_and_ions 0.833: salt_concentration | 3/3 true, 0 false | 0.0, 0.0, 0.0 | 1.0 | 0.1109 | 11.02, 10.84, 10.78 | 0.0 | 2.088 | 1.014 | 296.7 | 1.388 | 304.2 |

Generic attacks: `empty_output` 0.0, `header_only` 0.0, `inputs_copied` 0.0, `zero_byte_files` 0.0.

### Seed 2, level 3 (1AKI, padding 9.0 Å, 0.15 M)

```
wrote /home/user/ale-lifesci-gym/tasks/md-system-repair-s0002-l3  (prior difficulty 0.675)
  [ok] G1 images_build 
  [ok] G2 oracle score=1.0 runtime_min=4.5 failure= stderr=
  [ok] G5 verifier_deterministic first=1.0 second=1.0
  [ok] G3 probe:within_tolerance score=1.0 floor=0.99 failure=
  [ok] G3 probe:independent_solution score=1.0 floor=0.99 failure=
  [ok] G4 attack:empty_output score=0.0 ceiling=0.0 failure=required output file missing: system.prmtop
  [ok] G4 attack:inputs_copied score=0.0 ceiling=0.1 failure=required output file missing: system.prmtop
  [ok] G4 attack:zero_byte_files score=0.0 ceiling=0.0 failure=system.prmtop is empty
  [ok] G4 attack:header_only score=0.0 ceiling=0.1 failure=report.json is not valid JSON: Expecting property name enclosed in double quotes: line 2 column 1 (char 2)
  [ok] G4 attack:prior_only_generic_system score=0.0 ceiling=0.05 failure=gate failed: topology contains a polypeptide
  [ok] G4 attack:salt_beyond_tolerance score=0.683333 ceiling=0.75 failure=
  [ok] G4 attack:naive_reduced_cys score=0.62 ceiling=0.75 failure=
  [ok] G4 attack:naive_no_ions score=0.666667 ceiling=0.75 failure=
  [ok] G4 attack:naive_crystal_solvent_kept score=0.683333 ceiling=0.75 failure=
  [ok] G4 attack:naive_mse_dropped score=0.66972 ceiling=0.75 failure=
  [ok] G4 attack:naive_sidechains_truncated score=0.647005 ceiling=0.75 failure=
  [ok] G4 attack:naive_gap_bonded score=0.620667 ceiling=0.75 failure=
  [ok] G4 attack:naive_gap_charged_termini score=0.670488 ceiling=0.75 failure=
  [ok] G4 attack:naive_altloc_first score=0.685 ceiling=0.75 failure=
  [ok] G4 attack:naive_internal_numbering score=0.685 ceiling=0.75 failure=
  [ok] G4 attack:naive_raw_tleap score=0.446772 ceiling=0.6 failure=
{"task": "md-system-repair-s0002-l3", "ok": true, "minutes": 49.47}
```

| submission | score | failed checks | S-S | chi1 delta (deg) | atoms in place | salt (M) | padding x,y,z (Å) | xtal sites kept | max bond (Å) | density | T stored vel. (K) | CA RMSD eq. (Å) | T after verifier MD (K) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| oracle | 1.0 | - | 4/4 true, 0 false | 0.0, 0.0, 0.0 | 1.0 | 0.1485 | 11.74, 12.01, 11.74 | 0.0 | 2.018 | 1.025 | 302.5 | 0.649 | 302.8 |
| independent_solution | 1.0 | - | 4/4 true, 0 false | 0.0, 0.0, 0.0 | 1.0 | 0.1512 | 9.25, 9.25, 9.25 | 0.0 | 2.018 | 1.0171 | 303.4 | 0.748 | 305.0 |
| within_tolerance | 1.0 | - | 4/4 true, 0 false | 0.0, 0.0, 0.0 | 1.0 | 0.138 | 11.74, 12.01, 11.74 | 0.0 | 2.018 | 1.0228 | 300.4 | 0.836 | 306.1 |
| naive_altloc_first | 0.685 | conformers 0.75: major_conformer_188 | 4/4 true, 0 false | 0.0, 0.0, 120.0 | 0.9968 | 0.1485 | 11.74, 12.01, 11.74 | 0.0 | 2.018 | 1.0191 | 300.4 | 0.907 | 302.7 |
| naive_crystal_solvent_kept | 0.6833 | solvent_and_ions 0.833: no_crystal_solvent_or_ions | 4/4 true, 0 false | 0.0, 0.0, 0.0 | 1.0 | 0.1598 | 13.31, 12.14, 11.74 | 1.0 | 2.018 | 1.0154 | 303.0 | 0.633 | 302.2 |
| naive_gap_bonded | 0.6207 | termini_and_breaks 0.5: no_bond_between_unlinked_residues, break1_c_side_n_methylamide, break1_n_side_acetyl; physical_sanity 0.667: as_built_bond_lengths; equilibrated 0.8: density | 4/4 true, 0 false | 0.0, 0.0, 0.0 | 1.0 | 0.1522 | 11.51, 12.01, 11.74 | 0.0 | 18.948 | 0.8825 | 306.4 | 1.979 | 305.4 |
| naive_gap_charged_termini | 0.6705 | atoms_and_charge 0.959: heavy_atoms_complete; termini_and_breaks 0.667: break1_c_side_n_methylamide, break1_n_side_acetyl | 4/4 true, 0 false | 0.0, 0.0, 0.0 | 1.0 | 0.1523 | 11.55, 12.01, 11.74 | 0.0 | 2.018 | 1.0295 | 300.6 | 0.636 | 301.8 |
| naive_internal_numbering | 0.685 | report 0.833: disulfides | 4/4 true, 0 false | 0.0, 0.0, 0.0 | 1.0 | 0.1485 | 11.74, 12.01, 11.74 | 0.0 | 2.018 | 1.026 | 303.0 | 0.846 | 301.5 |
| naive_mse_dropped | 0.6697 | sequence 0.64: residue_identity, no_extra_or_missing_residues; atoms_and_charge 0.959: heavy_atoms_complete; termini_and_breaks 0.972: peptide_bonds | 4/4 true, 0 false | 0.0, 0.0, 0.0 | 0.9915 | 0.1485 | 11.74, 12.01, 11.74 | 0.0 | 2.018 | 1.0244 | 301.2 | 0.976 | 300.3 |
| naive_no_ions | 0.6667 | solvent_and_ions 0.667: neutral, salt_concentration | 4/4 true, 0 false | 0.0, 0.0, 0.0 | 1.0 | 0.0 | 11.74, 12.01, 11.74 | 0.0 | 2.018 | 1.0101 | 304.8 | 0.787 | 303.7 |
| naive_raw_tleap | 0.4468 | sequence 0.64: residue_identity, no_extra_or_missing_residues; atoms_and_charge 0.959: heavy_atoms_complete; disulfides 0.0: ss_bonds; termini_and_breaks 0.472: peptide_bonds, no_bond_between_unlinked_residues, break1_c_side_n_methylamide, break1_n_side_acetyl; conformers 0.745: resolved_atoms_in_place, major_conformer_188; solvent_and_ions 0.667: salt_concentration, no_crystal_solvent_or_ions; physical_sanity 0.667: as_built_bond_lengths; equilibrated 0.8: density; report 0.833: disulfides | 0/4 true, 0 false | 0.0, 0.0, 120.0 | 0.9883 | 0.0099 | 13.76, 12.14, 11.74 | 1.0 | 18.948 | 0.873 | 285.6 | 1.744 | 291.2 |
| naive_reduced_cys | 0.62 | disulfides 0.0: ss_bonds | 0/4 true, 0 false | 0.0, 0.0, 0.0 | 1.0 | 0.1485 | 11.74, 12.01, 11.74 | 0.0 | 1.841 | 1.014 | 305.0 | 0.742 | 306.1 |
| naive_sidechains_truncated | 0.647 | sequence 0.558: residue_identity, no_extra_or_missing_residues; atoms_and_charge 0.837: heavy_atoms_complete; termini_and_breaks 0.89: peptide_bonds; conformers 0.969: resolved_atoms_in_place | 4/4 true, 0 false | 0.0, 0.0, 0.0 | 0.9787 | 0.1484 | 11.74, 12.01, 11.74 | 0.0 | 2.018 | 1.0204 | 300.8 | 0.589 | 299.3 |
| prior_only_generic_system | 0.0 | gate failed: topology contains a polypeptide | - | - | - | - | - | - | - | - | - | - | - |
| salt_beyond_tolerance | 0.6833 | solvent_and_ions 0.833: salt_concentration | 4/4 true, 0 false | 0.0, 0.0, 0.0 | 1.0 | 0.1167 | 11.74, 12.01, 11.74 | 0.0 | 2.018 | 1.0216 | 302.4 | 0.704 | 304.0 |

Generic attacks: `empty_output` 0.0, `header_only` 0.0, `inputs_copied` 0.0, `zero_byte_files` 0.0.

### Seed 1, level 1 (4PTI, clean file)

```
wrote /home/user/ale-lifesci-gym/tasks/md-system-repair-s0001-l1  (prior difficulty 0.282)
  [ok] G1 images_build 
  [ok] G2 oracle score=1.0 runtime_min=1.63 failure= stderr=
  [ok] G5 verifier_deterministic first=1.0 second=1.0
  [ok] G3 probe:within_tolerance score=1.0 floor=0.99 failure=
  [ok] G3 probe:independent_solution score=1.0 floor=0.99 failure=
  [ok] G4 attack:empty_output score=0.0 ceiling=0.0 failure=required output file missing: system.prmtop
  [ok] G4 attack:inputs_copied score=0.0 ceiling=0.1 failure=required output file missing: system.prmtop
  [ok] G4 attack:zero_byte_files score=0.0 ceiling=0.0 failure=system.prmtop is empty
  [ok] G4 attack:header_only score=0.0 ceiling=0.1 failure=report.json is not valid JSON: Expecting property name enclosed in double quotes: line 2 column 1 (char 2)
  [ok] G4 attack:prior_only_generic_system score=0.0 ceiling=0.05 failure=gate failed: topology contains a polypeptide
  [ok] G4 attack:salt_beyond_tolerance score=0.683333 ceiling=0.75 failure=
  [ok] G4 attack:naive_reduced_cys score=0.62 ceiling=0.75 failure=
  [ok] G4 attack:naive_no_ions score=0.666667 ceiling=0.75 failure=
{"task": "md-system-repair-s0001-l1", "ok": true, "minutes": 10.94}
```

| submission | score | failed checks | S-S | chi1 delta (deg) | atoms in place | salt (M) | padding x,y,z (Å) | xtal sites kept | max bond (Å) | density | T stored vel. (K) | CA RMSD eq. (Å) | T after verifier MD (K) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| oracle | 1.0 | - | 3/3 true, 0 false | - | 1.0 | 0.1428 | 11.02, 10.84, 10.78 | 0.0 | 2.088 | 1.0071 | 307.5 | 0.794 | 302.6 |
| independent_solution | 1.0 | - | 3/3 true, 0 false | - | 1.0 | 0.1421 | 8.25, 8.25, 8.25 | 0.0 | 2.088 | 1.0144 | 304.2 | 0.802 | 302.1 |
| within_tolerance | 1.0 | - | 3/3 true, 0 false | - | 1.0 | 0.1428 | 11.02, 10.84, 10.78 | 0.0 | 2.088 | 1.004 | 303.1 | 0.728 | 300.9 |
| naive_no_ions | 0.6667 | solvent_and_ions 0.667: neutral, salt_concentration | 3/3 true, 0 false | - | 1.0 | 0.0 | 11.02, 10.84, 10.78 | 0.0 | 2.088 | 1.0035 | 298.4 | 0.849 | 308.1 |
| naive_reduced_cys | 0.62 | disulfides 0.0: ss_bonds | 0/3 true, 0 false | - | 1.0 | 0.1428 | 11.02, 10.84, 10.78 | 0.0 | 1.815 | 1.0161 | 303.0 | 1.016 | 300.7 |
| prior_only_generic_system | 0.0 | gate failed: topology contains a polypeptide | - | - | - | - | - | - | - | - | - | - | - |
| salt_beyond_tolerance | 0.6833 | solvent_and_ions 0.833: salt_concentration | 3/3 true, 0 false | - | 1.0 | 0.111 | 11.02, 10.84, 10.78 | 0.0 | 2.088 | 1.0028 | 304.6 | 0.758 | 306.4 |

Generic attacks: `empty_output` 0.0, `header_only` 0.0, `inputs_copied` 0.0, `zero_byte_files` 0.0.

## Additional measurements outside the three validation runs

These were single manual runs in the task images (1-2 CPU threads); the files were not kept.

| what | result |
|---|---|
| Oracle output with `system.inpcrd` and `equilibrated.rst7` converted to NetCDF (ParmEd), scored in the verifier image | 1.0, identical raw values (the verifier reads both encodings) |
| Level 2, seeds 1 and 2: both builders without real equilibration (tleap: none; Modeller route: 30 minimiser iterations + 20 MD steps) | every as-built check passes for both builders and both seeds; only `equilibrated` fails (density 0.765 / 0.786 for tleap, 0.940 / 0.947 for Modeller; stored-velocity temperature 297.8 / 429.9 K; temperature after verifier MD 330.9 / 357.6 K) |
| Level 3 seed 2 and level 1 seed 1, same smoke test before the full runs | as above (density 0.791 / 0.764 tleap, 0.947 / 0.939 Modeller) |
| Seed 1 level 3: tleap build, minimised only, Maxwell-Boltzmann velocities written by hand | `velocities_thermal` passes (298.9 K, X-H velocity along the bond 16.3 Å/ps, so unconstrained counting is used), `density` fails (0.759), `md_continues_stably` fails (234.9 K) |
| Ten naive builds of the seed-1 file without equilibration (development, padding 9 Å) | same trap-specific failed checks as in the table for seed 1; as-built energies up to +1.3e6 kJ/mol minimise to -1.4e5..-1.8e5 kJ/mol within 60 iterations |
| Density of the tleap box during NPT (seed-1 file, MC barostat): frequency 25 after 2 ps NVT | 0.82, 0.85, 0.89, 0.91, 0.93, 0.935, 0.963, 0.973 g/cm3 after 1..8 ps |
| same, frequency 10 after 1 ps NVT | 0.835, 0.877, 0.920, 0.942, 0.975, 0.990 after 1..6 ps |
| same, frequency 3 after 0.5 ps NVT | 0.808, 0.832, 0.880, 0.910, 0.928, 0.949 after 0.5..3 ps |
| `pdb4amber` (default options) on the seed-1 level-3 file | keeps first conformers, leaves MSE as HETATM, keeps waters/additives, TER at the gap, spurious TER records at insertion-coded residues |

History of the bands: a first complete validation of seed 1 level 3 (all gates green, 47 min) used
density 0.95-1.06 g/cm3, 260-340 K for the verifier MD and 97% for atoms in place. The smoke tests
above then showed wrong-side values at 0.947 g/cm3 and 330.9 K, i.e. at or inside those edges, so
the bands were moved to the midpoints between measured correct and wrong values (0.97-1.08 g/cm3,
280-320 K, 99%) and all three validations were repeated with the final verifier. No band was
widened to obtain a pass.

## Solutions

- **Oracle** (`solve.sh` -> `build_tleap.py`): text-level repair in Python (highest-occupancy conformer,
  MSE -> MET, HETATM groups dropped, chain break found from C-N geometry, one anchor atom per cap,
  disulfides from SG-SG distances), `pdb4amber` as a consistency check of the edited file, `tleap`
  (ff14SB, TIP3P, `solvatebox`, `addIonsRand` with the SOP count), OpenMM for the SOP equilibration
  (200 L-BFGS iterations, 0.5 ps NVT, 6 ps NPT with a Monte-Carlo barostat), report recomputed from
  the written topology with ParmEd.
- **Second solution** (`SOLUTION=2`, `solve_alt.py`): no pdb4amber, no tleap. Own PDB parser, caps
  placed atom by atom (NeRF), PDBFixer for the missing heavy atoms, OpenMM Modeller for hydrogens and
  the water box (explicit `boxSize`), own ion placement by replacing waters, ff14SB from the OpenMM
  XML files, ParmEd export to prmtop/inpcrd (plus the H-H water bond of the Amber TIP3P library),
  equilibration with a different integrator, barostat frequency and length (3 ps NPT).
  The two solutions share OpenMM as MD engine and ParmEd as restart writer; everything that decides
  the scored structure (repair, capping, hydrogens, solvation, ions, parameters) is separate code and
  separate tools. `sander` was not used for MD: about 10 steps/s single-threaded on these systems.

## Battery

Probes: `within_tolerance` (reference builder on the hidden clean model, then NaCl pairs removed
while the deficit stays <= 0.015 M and the reported box rounded to 0.01 Å), `independent_solution`.
There is no separate `planted_truth` probe: the deliverables cannot be written from the truth
without a builder, and the unperturbed truth build is strictly closer to every target than the
`within_tolerance` build, so it was merged to save one MD run per validation.

Attacks: `prior_only_generic_system` (solvated ACE-ALA3-NME, constant report; ceiling 0.05),
`salt_beyond_tolerance` (truth build with NaCl pairs removed until the deficit is >= 0.025 M),
one naive variant per decision point, `naive_no_ions`, and `naive_raw_tleap` (ceiling 0.60).
Naive variants are flags of the reference builder, so each differs from the oracle in exactly one
decision and runs the full reference equilibration; the only failing checks are those of the
ignored requirement (tables below). Two exceptions: `naive_gap_bonded` and `naive_raw_tleap` fail
several components by construction and run only 1 ps of NPT to save time, so their density check
fails as well.

## Known limitations and open points

0. **The literal `make ... --seed 1 --level 3` was not repeated after the last verifier change.**
   It passed every gate including Harbor with the verifier as it was at 16:26 UTC. Three constants
   were then tightened (see "History of the bands"). With the final verifier, seed 2 level 3 and
   seed 1 level 1 were run through `make` from scratch; seed 1 level 3 was re-validated by
   `results/md_system_repair_s1_l3_replay.py`, i.e. the real `validate()` with all gates, fresh
   images, fresh verifier runs and a fresh Harbor oracle/nop run, but with the 17 script outputs of
   the earlier run replayed instead of recomputed (no solution, probe or attack script changed in
   between). A from-scratch `make` takes 35-50 min. The replayed outputs were deleted afterwards to
   free disk, so the replay script documents the procedure but cannot be re-run as is.
1. **Force-field family is not verified.** ParmEd's export of an OpenMM system renames atom types
   (`N1`, `C1`, ...), so type names are not comparable across tool chains, and ff14SB, ff99SB and
   ff19SB share the same charges. A system built with another Amber protein force field plus TIP3P
   passes. Ion parameters are not checked either. Only TIP3P charges, total charge and composition are.
2. **"Was it equilibrated" rests on three outcome checks, none of which sees the protocol.**
   Density separates a tleap box that never saw NPT (0.76-0.79 g/cm3) or saw 1 ps of it (0.85) and
   a Modeller box (0.94-0.95 before NPT) from relaxed cells (1.00-1.03). The stored velocities only
   show that a thermal velocity set is present: hand-written Maxwell-Boltzmann velocities pass.
   The temperature after 0.3 ps of verifier MD rejects unrelaxed (331-389 K) and minimised-only
   (235 K) configurations. A submission could still satisfy all three without following the three
   SOP stages literally (for example NPT from the first step).
3. **The equilibration protocol itself (thermostat, time step, cut-off, length) is not verifiable**
   from the hand-over files; only its outcome is scored.
4. **SOP 3.3 "until the density has levelled off"** is implemented in the oracle as a fixed 6 ps.
   Measured final densities are 1.00-1.03 g/cm3, but the oracle does not test for a plateau.
5. **Salt tolerance is about one NaCl pair for the small system** (0.016 M per pair at 3500 waters
   against a tolerance of 0.02 M). The SOP gives the count formula explicitly, and both solutions
   land within 0.008 M, but an agent that derives the count from the volume of the as-built tleap
   cell (density 0.77) would add ~45% more pairs and fail. I consider that a correct fail (the SOP
   formula is explicit) but it is the tightest tolerance in the task.
6. **The SOP states an acceptance window for the padding** ("at least P ... not more than 5 Å
   excess"). This is a threshold in an input document, on purpose: tleap's `solvatebox` buffer gives
   P+2.7..3.0 Å measured from atom centres, so "padding P" alone would be ambiguous. The verifier
   adds 0.5 Å slack on both sides.
7. **`missing_sidechains` is a defect, not a trap**, for anyone using tleap or PDBFixer.
8. **Two source proteins only**, and SEQRES shows the sequence, so a model can recognise BPTI or
   lysozyme. Knowing the entry gives no scored advantage beyond the disulfide pairing, which is also
   in the geometry; rebuilt atoms and the unresolved segment are not scored against the entry.
   Coordinates are in the frame of the RCSB entry (no rigid-body blinding).
9. **Histidine**: the SOP fixes neutral His; only the net charge is checked (lysozyme has one His,
   BPTI none), not the tautomer.
10. **`within_tolerance` is weak on the salt side for the small system**: removing a single pair
    already leaves the 0.015 M window, so for seed 1 the probe only exercises the rounding of the
    reported box. For seed 2 it removes pairs (see table).
11. **Level 2 was generated and smoke-tested (table above) but not run through `make`.**
12. **Image size and runtime.** Agent image: 2.24 GB conda layer (`docker images` reports 3.26 GB with
    the containerd store) after deleting the CUDA runtime compiler that the openmm package pulls in
    (`rm -rf /opt/tools/targets` in the same RUN; this is a small departure from the plain
    `dockerlib` Dockerfile strings). Verifier image 1.2 GB as reported. Oracle with 2 OpenMM threads:
    1.6 min (BPTI, idle host), 3.2 min (BPTI, loaded host), 3.8-4.5 min (lysozyme, loaded host).
    Verifier 20-45 s. A level-3 validation took 47-49 min on the loaded host because every naive
    variant needs its own equilibration; level 1 took 11 min on the idle host.
13. **The verifier depends on ParmEd/OpenMM parsers for untrusted files.** Nothing is executed or
    unpickled, sizes and atom counts are capped, and parser exceptions become score 0 (or a failed
    check for the restart), but the parsers themselves are not hardened code.
14. **Not exercised:** the 0.10 M salt option (both validated seeds drew 0.15 M; padding 8.0 and 9.0 Å
    were both exercised), seeds other than 1 and 2 (different gap/conformer/numbering draws on the
    same two proteins), and any agent run. The planting code raises if it cannot place a gap, decoy
    conformers or additives, so an unlucky seed fails at generation rather than producing a wrong task.
