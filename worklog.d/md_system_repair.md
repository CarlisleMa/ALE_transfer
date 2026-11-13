# Worklog: `md_system_repair`

All entries 2026-11-11 (UTC, build host shared with two other authors, 4 cores).

## 15:25 - reading, tool choice
- Read AUTHORING.md, the spec, the crispr_screen reference family and the core modules.
- conda-forge versions available for python 3.12: ambertools 25.3 / 26.0, openmm 8.6.1, pdbfixer 1.12,
  parmed 4.3.1. Chose `ambertools=25.3 openmm=8.6.1 pdbfixer=1.12` for the agent image (ParmEd 4.3.1 and
  pdb4amber come with AmberTools) and `openmm=8.6.1 parmed=4.3.1` for the verifier.
- Built a tools image with that list to experiment: 2.61 GB conda layer (`docker images` shows 3.87 GB
  because the containerd store counts compressed + unpacked). `/opt/tools/targets` (368 MB) is the CUDA
  runtime compiler pulled in by the openmm package; OpenMM only offers `Reference` and `CPU` here.

## 15:30 - source structures
- Downloaded 1AKI, 4PTI, 1BPI, 1CRN, 7RSA, 2LZT, 1HEL from `https://files.rcsb.org/download/<ID>.pdb`
  into `.cache/md_system_repair/`. Kept 1AKI (lysozyme, 129 aa, 4 SS, 2 Met) and 4PTI (BPTI, 58 aa, 3 SS,
  1 Met): both are complete (every heavy atom incl. OXT), have no alternate locations of their own and no
  hydrogens. Rejected 7RSA (hydrogens, 12 real altloc residues, no SSBOND records, 4 His), 1BPI (real
  altlocs, PO4), 1CRN (no Met, so the modified-residue defect could not be planted), 2LZT/1HEL (same
  protein as 1AKI).

## 15:35 - generator (pdbkit.py, plant.py)
- Truth = the complete entry. `read_entry` asserts completeness against a heavy-atom table and
  cross-checks the deposited SSBOND pairs against SG-SG geometry (< 2.5 A).
- Planting: internal gap (4-6 residues, no Cys/Met/Pro in the segment or its flanks, C-N distance of the
  flanks > 6 A), MSE, side chains truncated at CB/CG, decoy conformers (side chain rotated by +-120 degrees
  about CA-CB, accepted only if >= 3.0 A from every other residue), numbering offset + insertion codes
  placed so that a disulfide cysteine carries one, chain renamed, additives placed at crystal-water
  sites, crystal waters of the entry kept.
- First test of the planted files: OK for both entries, level 1 and 3.

## 15:40 - what the tools do with the defective file (measured, AmberTools 25.3)
- `pdb4amber` default: keeps the FIRST conformer ("Only the first occurrence for each atom was kept"),
  does not convert MSE, keeps waters and additives, inserts TER at the gap (-> charged termini in tleap),
  and is confused by insertion codes (spurious TER records around `67A`/`67B`). So a "pdb4amber then
  tleap" default run fails several SOP clauses at once; each became a naive variant.
- `tleap`: rebuilds missing side-chain atoms and the rest of ACE/NME from a single anchor atom
  ("Added missing heavy atom: .R<NME 20>.A<C 3>"), so the reference only places the acetyl C and the
  N-methylamide N explicitly.
- tleap `solvatebox m TIP3PBOX 10.0`: box 59.9 x 54.0 x 58.9 A for a solute of 33 x 28 x 34 A, i.e. the
  padding measured from atom centres is ~P + 3 A, and the as-built density is only 0.76-0.78 g/cm3.

## 15:45 - equilibration cost (the main budget problem)
- OpenMM CPU, 2 threads, 12.9-14.6 k atoms: ~37 MD steps/s. `sander` (serial) was not tried for MD: an
  estimated ~10 steps/s is too slow for a battery of ~15 runs per task.
- Density relaxation of the tleap box (MC barostat):
  - frequency 25, after 2 ps NVT: 0.82, 0.85, 0.89, 0.91, 0.93, 0.935, 0.963, 0.973 g/cm3 after 1..8 ps
  - frequency 10, after 1 ps NVT: 0.835, 0.877, 0.920, 0.942, 0.975, 0.990 after 1..6 ps
  - frequency 3, after 0.5 ps NVT: 0.808, 0.832, 0.880, 0.910, 0.928, 0.949 after 0.5..3 ps
  More barostat attempts help only a little; the cell needs ~6 ps of real time. Settled on 200 L-BFGS
  iterations, 0.5 ps NVT, 6 ps NPT (frequency 5) for the tleap route and 3 ps NPT for the Modeller route
  (its box starts at ~1.0 g/cm3).
- Considered and rejected: shrinking the tleap box by hand before MD (non-standard, risk of periodic
  clashes); a lower density bound of 0.90 with a 3 ps protocol (would let "not levelled off" pass and
  could not be justified for other barostats).

## 15:50 - verifier (score.py)
- ParmEd for topology/coordinates (ASCII and NetCDF restarts read back identically; garbage files raise
  and become CandidateError), OpenMM CPU for energies, minimisation and 0.3 ps MD.
- First run on the oracle output: 1.0 in 17 s (2 threads). Second solution: 1.0 in 28 s.
- All ten naive builds (without equilibration, to iterate quickly): every trap-specific check fired as
  intended; no verifier error. Numbers in docs/families/md_system_repair.md.
- NetCDF restarts need scipy (ParmEd fallback) -> added `scipy=1.18.1` to the verifier image.

## 15:57 - coordinator notice: disk quota
- Removed my verifier test image, later the tools image (the generator needs no tool image; development
  continued on the task images). Agent/verifier Dockerfiles now delete `/opt/tools/targets` in the same
  RUN as the conda install (-368 MB each). Only the images of the task being validated are kept.

## 16:00 - second solution
- `solve_alt.py`: hand-written repair (occupancy-based conformer choice, MSE->MET, explicit NeRF
  placement of all cap heavy atoms), PDBFixer for missing atoms, Modeller for hydrogens and the TIP3P
  box (explicit `boxSize` = extent + 2P + 0.5 A; the `padding=` argument of `addSolvent` builds a cube
  from the bounding sphere and misses the SOP), own ion placement, ParmEd export.
- The exported topology had only two O-H bonds per water, so OpenMM treated the H-O-H angle as flexible
  (7.7 k constraints instead of 11.4 k). Added the H-H bond the Amber water library carries.

## 16:12 - first `make` attempt for seed 1 level 3, stopped on purpose
- Oracle 1.0 (3.1 min on the loaded host), `within_tolerance` 1.0. While it ran, a smoke test of the
  second solution with 20 MD steps showed a 416 K restart passing `velocities_thermal`: the check
  accepted the temperature under either degree-of-freedom count. Stopped the run and changed the
  check: the verifier now decides from the stored velocities whether X-H bonds were constrained
  (relative velocity along the bond 0.8-1.2 A/ps when constrained, ~16 when not; threshold 5) and
  uses one count. Narrowed the verifier-MD temperature band from 200-400 to 260-340 K.

## 16:26-17:13 - complete `make` seed 1 level 3 (verifier before the last band change)
- All gates green including Harbor oracle = 1.0 / nop = 0.0, 47 min. Log kept as
  the first block of `results/md_system_repair_s1_l3.log`.
- Every single-decision naive variant failed only the checks of its own requirement while its
  equilibration passed; `naive_gap_bonded` and `naive_raw_tleap` (1 ps NPT) also failed density.

## 17:14-17:28 - seed 2 started; smoke tests between runs; bands re-placed
- Seed 2 (lysozyme) oracle 1.0, density 1.027.
- Level-2 smoke test (both builders, both seeds, no real equilibration) and a NetCDF copy of the
  oracle output: as-built checks all pass; NetCDF scores 1.0.
- The same smoke tests gave wrong-side values at or inside two band edges: an unequilibrated Modeller
  box at 0.947 g/cm3 (band started at 0.95) and 330.9 K after verifier MD for a 20-step "equilibration"
  (band ended at 340 K). Moved the bands to the midpoints between measured correct and wrong values:
  density 0.97-1.08 g/cm3, verifier-MD temperature 280-320 K, atoms in place >= 99% (was 97%, wrong
  values at 95.4-98.0%). Stopped the seed-2 run (oracle and `within_tolerance` had passed) and
  restarted it with the final verifier.
- Measured afterwards: minimised-only tleap build with hand-written Maxwell-Boltzmann velocities ->
  234.9 K after verifier MD (fails), stored-velocity check passes (298.9 K).

## 17:45 - coordinator notice (time and disk)
- Only the images of the task being validated exist; a background driver removes them after each run.
- To hand back sooner, seed 1 level 3 is not run a second time from scratch. Because no solution,
  probe or attack script changed after its complete run (only three verifier constants), the driver
  re-validates it with `results/md_system_repair_s1_l3_replay.py`: the normal `validate()` with all
  gates, fresh images, fresh verifier runs and a fresh Harbor oracle/nop run, but each script call is
  answered with the output directory the same script produced in the 16:26 run.
- Order after seed 2: seed 1 level 3 replay, then seed 1 level 1.

## 17:28-18:17 - `make` seed 2 level 3, final verifier
- All gates green (no Harbor, as specified), 49.5 min. Oracle 1.0 (4.5 min), second solution 1.0.
- Also recorded what tleap really prints for the raw seed-1 file (unknown residues, "does not have a
  type", bonds of 13.705 and 3.291 A as warnings only, silent acceptance of HOH/CL/NA) and replaced
  the assumed wording in the trap evidence with the measured messages.

## 18:17-18:26 - seed 1 level 3 re-validated with the final verifier (replay)
- All gates green, Harbor oracle = 1.0 / nop = 0.0 (fresh oracle run inside Harbor), 8.9 min.
- Scores moved only where the tightened "atoms in place" guard now also fires:
  `naive_altloc_first` 0.670 -> 0.667, `naive_mse_dropped` 0.661 -> 0.659,
  `naive_sidechains_truncated` 0.630 -> 0.627, `naive_raw_tleap` 0.420 -> 0.417.

## 18:27-18:38 - `make` seed 1 level 1, final verifier
- All gates green, 10.9 min on the now idle host (oracle 1.63 min).

## 18:40 - wrap-up
- Removed all my Docker images (the generator needs none; `make` rebuilds the task images), the
  unused cached PDB entries, the kept output directories under `results/work/md-system-repair-*`
  (the small `logs/details.json` files stay) and my scratch files.
- `tasks/md-system-repair-s0002-l3` was rewritten once more with `--no-validate` so that its README
  and task.toml carry the final trap wording and the measured oracle runtime; only those two files
  changed (checked with `diff -rq`).
- Not done: a from-scratch `make` of seed 1 level 3 with the final verifier (see limitation 0 in the
  design notes), level 2 through `make`, any seed beyond 1 and 2.

## Results
Gate lines and measured tables: docs/families/md_system_repair.md.
