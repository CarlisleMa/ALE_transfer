# Family spec: `md_system_repair` - from a defective crystal structure to a simulation-ready system

ALE subdomain: 3.1 Biomolecular Structure & Design. (ALE has an Amber MMGBSA scripting task; this is
the upstream, different skill: structure repair and system building, graded on the built system.)

## Scenario
A computational chemist receives a raw PDB entry-like file and the group's simulation SOP and must
deliver an Amber-format, solvated, neutralised system plus short equilibration, ready for
production MD.

## Truth model
Start from a small, complete, real protein with disulfides (e.g. lysozyme 1AKI, RNase A, BPTI,
crambin; choose by seed from a short list fetched from RCSB). The complete structure is the hidden
truth. The generator *plants defects* into the file the agent sees (level 3 = all; choice and
positions vary with seed):
- side-chain heavy atoms deleted for several residues;
- an internal segment of 4-6 residues removed (chain break, REMARK 465 listing them);
- alternate-location records (A/B) with occupancies for 2-3 residues;
- one MET written as selenomethionine (MSE, HETATM) when the protein has Met;
- crystallisation additives (GOL, SO4, ...) and crystal waters added;
- SSBOND/CONECT records removed (disulfides only detectable from SG-SG geometry);
- residues renumbered with an offset and an insertion code; chain renamed.

## Inputs
`structure.pdb` (the defective file), `simulation_sop.md` (force field family, water model, box
padding, salt concentration, how unresolved segments are to be treated: not modelled, neutral caps;
which heteroatoms to keep/convert; equilibration protocol), no reference structure.
Tools: AmberTools (tleap, pdb4amber, cpptraj, sander), OpenMM, PDBFixer, ParmEd, python.

## Outputs
`system.prmtop`, `system.inpcrd` (as built), `equilibrated.rst7` (after the SOP's short protocol),
`report.json` (disulfide pairs in input numbering, protein net charge before ions, counts of
Na+/Cl-/water, box dimensions).

## Verifier (runs OpenMM + ParmEd in the verifier image; never trusts numbers the agent reports
without recomputing where possible)
gate: files load; protein backbone RMSD to the hidden complete structure small (same protein, same
fold). Components (weights to calibrate): residue sequence as expected incl. caps at the break and
MSE->MET; all expected heavy atoms present; each true disulfide present as an SG-SG bond; no
covalent bond across the unresolved segment; total charge ~ 0 and ion count consistent with the SOP
concentration; water model and padding per SOP; physical sanity from the verifier's own short
minimisation + MD on the submitted system (finite, stable, no blow-up); report.json consistent with
recomputed values.

## Traps
hidden_disulfides; chain_break (must not be bonded across); altlocs; modified_residue;
additives_and_waters; missing_sidechains; numbering/insertion code.

## Naive variants
CYS left reduced; break bonded or left with charged termini contrary to SOP; additives kept or the
MSE residue dropped; no ions; default tleap on the raw file (make it produce what a careless user
ships after "fixing errors until it runs").
