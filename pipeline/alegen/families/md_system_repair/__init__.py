"""Family: from a defective crystal structure to a simulation-ready Amber system.

Truth by construction: the generator starts from a small, complete, real protein with disulfides
(RCSB entry, cached), which stays hidden, and *plants* defects into the file the agent receives.
Everything the verifier needs (expected residues and heavy atoms, true conformers, disulfide
pairs, where the chain break is, where the crystal solvent was) is known from the planting, never
from the output of a preparation tool.

Decision points (level 2 adds the first four, level 3 adds the last three):
  hidden_disulfides       SSBOND/CONECT removed, cysteines named CYS
  additives_and_waters    crystal waters, glycerol/sulfate/ethanediol and Na+/Cl- sites in the file
  modified_residue        one Met written as selenomethionine (HETATM MSE)
  missing_sidechains      side chains truncated at CB/CG for several residues
  chain_break             4-6 internal residues unresolved (REMARK 465): caps, no bond across
  altlocs                 A/B conformers; for at least one residue the minor conformer is listed first
  numbering               offset, insertion codes (one on a disulfide cysteine), chain renamed
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np

from ... import dockerlib
from ...spec import Attack, Probe, TaskBuild, Trap
from . import pdbkit, plant

ASSETS = Path(__file__).parent / "assets"
CACHE = Path(__file__).resolve().parents[4] / ".cache" / "md_system_repair"
ENTRIES = {   # RCSB entry -> sha256 of https://files.rcsb.org/download/<ID>.pdb
    "1AKI": "c75d7a689617248cdd92dc6633531d2506fb9bef1e6e21e26c8f579ae6955abb",   # hen egg-white lysozyme, 129 aa, 4 SS
    "4PTI": "cd9b0faf2e66c0e031f30d1c84cb359691b8f5819f1b0f784be15f6005d58e08",   # BPTI, 58 aa, 3 SS
}
ORDER = ["1AKI", "4PTI"]
AGENT_CONDA = ["ambertools=25.3", "openmm=8.6.1", "pdbfixer=1.12"]
VERIFIER_CONDA = ["openmm=8.6.1", "parmed=4.3.1", "scipy=1.18.1"]   # scipy: NetCDF restarts in ParmEd
# the CUDA runtime compiler pulled in by the openmm package is useless without a GPU (368 MB)
_TRIM = ("&& micromamba clean -a -y\n", "&& micromamba clean -a -y \\\n    && rm -rf /opt/tools/targets\n")


def _entry(pid: str) -> Path:
    CACHE.mkdir(parents=True, exist_ok=True)
    p = CACHE / f"{pid}.pdb"
    if not p.exists():
        subprocess.run(["curl", "-fsSL", f"https://files.rcsb.org/download/{pid}.pdb", "-o", str(p)], check=True)
    if hashlib.sha256(p.read_bytes()).hexdigest() != ENTRIES[pid]:
        raise RuntimeError(f"{p} does not match the recorded sha256")
    return p


SOP = """\
# SOP MD-07 - preparing crystallographic protein models for explicit-solvent MD (Amber)

Revision 4. Applies to every single-chain protein job routed to the production queue. The SOP is
general; not every clause applies to every entry.

## 1. What is simulated

1.1 The polypeptide is simulated exactly as far as it was resolved in the experimental model.
    Residues without coordinates are not modelled (no loop building, no template grafting).
1.2 Ends of the polypeptide that are true ends of the expressed construct (per the sequence
    record of the entry) are simulated in their charged forms.
1.3 Ends that exist only because neighbouring residues are unresolved are blocked with neutral
    groups: an acetyl group (ACE) on an N-terminal end, an N-methylamide group (NME) on a
    C-terminal end.
1.4 Every resolved residue is simulated with the complete side chain of the amino acid given by
    the sequence record; atoms without coordinates are rebuilt.
1.5 Where more than one conformer was modelled, the conformer with the highest occupancy is
    used. Apart from a rigid-body placement in the cell, resolved atoms start from their deposited
    coordinates.
1.6 Residues that were chemically modified only to solve the structure or as an expression
    artefact (selenomethionine and the like) are simulated as their parent standard amino acid.
1.7 Covalent links between residues that exist in the crystal (disulfide bridges) exist in the
    topology.
1.8 Nothing else from the crystal is simulated: ordered water, buffer and precipitant molecules,
    cryoprotectants and crystallographic ions are discarded. (Cofactors and ligands need a
    separate sign-off and do not occur in jobs routed through this SOP.)
1.9 Protonation corresponds to pH 7 with standard pKa values: Asp and Glu deprotonated, Lys and
    Arg protonated, His neutral (either tautomer), Cys in a disulfide without thiol hydrogen.

## 2. Force field, solvent, ions

2.1 Protein: Amber ff14SB. Water: TIP3P. Ions: Na+ / Cl- with parameters matched to TIP3P.
2.2 Periodic cell: rectangular (all angles 90 degrees). Solvent padding: {pad:.1f} Å. On every cell axis
    the edge is at least the extent of the solute along that axis (all solute atoms, hydrogens
    included) plus twice the padding. More than 5 Å of excess padding on a side is not accepted
    (queue cost).
2.3 The system is neutralised with Na+ or Cl- (whichever is needed, never both for this
    purpose). In addition {salt:.2f} M NaCl is added; the number of added NaCl pairs is
    N = {salt:.2f} x N_water / 55.5, rounded, where N_water is the number of water molecules in
    the final system.

## 3. Equilibration (short pre-production)

3.1 Energy minimisation until the worst contacts are removed (a few hundred steps).
3.2 Thermalisation at constant volume: 300 K, at least 0.5 ps.
3.3 Constant pressure: 300 K, 1 bar, at least 3 ps and until the density of the cell has levelled off.
3.4 Bonds to hydrogen constrained, 2 fs time step, particle-mesh Ewald, direct-space cut-off 8-9 Å.

## 4. Hand-over

The system as built (before any minimisation) and the last frame of step 3.3 (coordinates,
velocities, cell) are handed over together with a short machine-readable report. File names and
formats are fixed by the job ticket.
"""

INSTRUCTION = """\
You are the computational chemist of a structural-biology group. A crystallographer handed over a
refined protein model. Turn it into a simulation-ready, explicitly solvated Amber system that
follows the group's SOP, run the SOP's short equilibration, and report the key properties.

Task directory:
- `/task`

Inputs (read-only) under `/task/input`:
- `structure.pdb` - the refined crystallographic model as exported by the crystallographer (PDB format)
- `simulation_sop.md` - the group's SOP for preparing MD systems (what is simulated, force field,
  solvent, ions, equilibration)

Available software (on `PATH`): AmberTools 25 (`tleap`, `pdb4amber`, `cpptraj`, `sander`, `parmed`, ...),
OpenMM 8.6, PDBFixer, ParmEd, Python 3.12 with numpy, scipy and pandas. 4 CPU cores, no GPU, no
internet access.

Your task:
1. Build the system that the SOP defines for `structure.pdb`.
2. Run the SOP's equilibration on it.
3. Write the report.

Required outputs under `/task/output`:
- `system.prmtop` - Amber parm7 topology of the complete system (solute, water, ions).
- `system.inpcrd` - Amber coordinate/restart file (ASCII or NetCDF) of the system exactly as built:
  after solvation and ion placement, before any minimisation or dynamics. Same atom order as
  `system.prmtop`; it must contain the periodic box.
- `equilibrated.rst7` - Amber restart file (ASCII or NetCDF) holding the last frame of the SOP
  equilibration: coordinates, velocities and box, same atom order as `system.prmtop`.
- `report.json` - a JSON object with these keys:
  - `disulfides`: list of the disulfide bridges in the system. Each bridge is a list of two strings:
    the residue numbers of its two cysteines exactly as written in `structure.pdb`, including the
    insertion code where a residue has one (for example `["27", "112A"]`). Order does not matter.
  - `protein_net_charge`: net charge in units of e (integer) of the solute in `system.prmtop`, i.e.
    of everything except water and the added ions.
  - `n_na`, `n_cl`, `n_water`: number of Na+ ions, Cl- ions and water molecules in `system.prmtop`.
  - `box_angstrom`: `[a, b, c]`, the box edge lengths in Å of `equilibrated.rst7`.

Constraints:
- Do not modify anything under `/task/input`; write deliverables only under `/task/output`.
- The SOP is binding wherever it applies to this entry.
- Do not ask for confirmation. Execute directly.
"""

JUDGE = """\
**Where the truth comes from.** The generator starts from a complete RCSB structure (hidden) and
plants every defect itself, so `reference/truth.json` holds planted facts: the residues that remain
after the unresolved segment was cut out (names, heavy atoms, coordinates), which atoms were left
in the file, the true conformer of every residue that received a decoy conformer, the disulfide
pairs (deposited SSBOND records, cross-checked against SG-SG geometry at generation time), the
net charge implied by the SOP protonation rules, the positions of the crystal waters/ions written
into the file, and the SOP padding and salt concentration. No preparation tool is run to produce
the reference.

**Nothing in `report.json` is trusted.** The verifier parses `system.prmtop` with ParmEd, the
coordinates with ParmEd, superposes the as-built CA atoms on the hidden structure, and recomputes
every quantity; the physics checks run OpenMM (CPU) on the submitted files.

**Hard gates (score 0).** A file missing/unparseable; atom counts of topology and as-built
coordinates differ; no periodic box; fewer than 70% of the expected residues found in sequence
order; as-built CA RMSD to the hidden structure > 1.0 Å (measured for both solutions: 0.000 Å).

**Components.** Each component is the mean of its checks; `production_ready` is 1 only if every
check of every other component passed. A system that ignores one requirement therefore scores at
most 0.70 and cannot reach the pass threshold of 0.90, while partial work still earns the
per-requirement credit.

| weight | component | checks |
|---|---|---|
| 0.07 | sequence | residue identities equal the resolved sequence (modified residue as parent amino acid); no missing/extra residue; no foreign residue |
| 0.07 | atoms_and_charge | every expected residue has exactly its heavy atoms; solute net charge equals the SOP value |
| 0.08 | disulfides | Jaccard index of S-S bonds in the topology against the true pairs (a bonded SG that still carries H counts as false) |
| 0.08 | termini_and_breaks | true termini charged; all peptide bonds inside resolved fragments; no other inter-residue bond (none across the unresolved segment); break capped with N-methylamide / acetyl |
| 0.06 | conformers | >= 99% of resolved heavy atoms within 1.0 Å of the hidden structure in the as-built file; chi1 of every residue with a planted decoy conformer within 40 degrees of the true conformer |
| 0.10 | solvent_and_ions | total charge 0 (+-0.05 e); only Na+/Cl- ions; salt concentration by the SOP formula within +-0.02 M and excess ions equal to the solute charge; TIP3P charges on 3-site water; padding per axis within [P-0.5, P+5.5] Å in an orthorhombic cell; <= 20% of the planted crystal water/ion sites occupied (0.5 Å) in the as-built file |
| 0.07 | physical_sanity | no solute heavy-atom bond longer than 2.6 Å as built; OpenMM builds a PME system; 60 L-BFGS iterations reach a finite negative energy |
| 0.08 | equilibrated | restart readable with orthorhombic box; density in [0.97, 1.08] g/cm3; kinetic temperature of the stored velocities in [250, 350] K (degrees of freedom counted with constraints when the stored velocities show constrained X-H bonds); CA RMSD to the hidden structure <= 3.0 Å; 0.3 ps of verifier MD from the restart stays finite, negative in energy, between 280 and 320 K, no atom moving 8 Å |
| 0.09 | report | each of the six reported items equals the recomputed value (disulfide labels in input numbering; charge; counts; box of the restart within 0.02 Å) |
| 0.30 | production_ready | all of the above |

**What is deliberately not scored.** The histidine tautomer, rotamers of rebuilt side chains,
the positions of rebuilt atoms and caps, ion positions, the exact equilibration recipe, and the
force-field family beyond what the checks above imply (see design notes).
"""


def _clean_pdb(clean, chain: str) -> str:
    """The hidden model as a perfect repair would leave it (input numbering, true conformers, MET)."""
    out, serial = [], 0
    for r in clean:
        num, ic = pdbkit.split_label(r.label)
        for at in r.atoms:
            serial += 1
            out.append(pdbkit.atom_line("ATOM", serial, at.name, " ", r.name, chain, num, ic.strip(), at.xyz,
                                        1.0, at.b, at.elem))
    out.append("END")
    return "\n".join(x.rstrip() for x in out) + "\n"


def generate(seed: int, level: int, workdir: Path) -> TaskBuild:
    pid = ORDER[seed % len(ORDER)]
    entry = pdbkit.read_entry(_entry(pid))
    planted = plant.plant(entry, np.random.default_rng(5000 + seed), level)
    prng = np.random.default_rng(9000 + seed)
    pad = float(prng.choice([8.0, 9.0]))
    salt = float(prng.choice([0.10, 0.15]))
    truth = planted["truth"]
    truth["sop"] = {"padding": pad, "salt_molar": salt}
    truth["source_entry"] = pid
    sop_md = SOP.format(pad=pad, salt=salt)
    clean_pdb = _clean_pdb(planted["clean"], truth["plan"]["chain"])

    builder = {"build_tleap.py": ASSETS / "build_tleap.py"}
    truth_in = {**builder, "perturb.py": ASSETS / "perturb.py", "truth_in/structure.pdb": clean_pdb,
                "truth_in/simulation_sop.md": sop_md}
    pre = "export OPENMM_CPU_THREADS=2\n"
    run_truth = pre + "python3 /alegen_script/build_tleap.py --input /alegen_script/truth_in\n"

    def naive(flag: str, short: bool = False) -> str:
        # `short`: variants that fail several checks by construction do not need the long NPT stage
        return pre + ("export EQ_NPT_STEPS=500\n" if short else "") + f"python3 /alegen_script/build_tleap.py {flag}\n"

    probes = [
        Probe("within_tolerance", "planted truth inside the tolerances: reference builder applied to the hidden "
              "clean model (true conformers, complete residues, no hetero groups), then NaCl pairs removed while "
              "the concentration stays within 0.015 M of the SOP value and the reported box rounded to 0.01 Å", 0.99,
              script=run_truth + f"python3 /alegen_script/perturb.py --salt {salt} --max-deficit 0.015 --round-box 2\n",
              script_files=truth_in),
        Probe("independent_solution", "hand-written repair + PDBFixer + OpenMM Modeller + ParmEd export "
              "(no pdb4amber, no tleap)", 0.99, script="SOLUTION=2 bash /solution/solve.sh"),
    ]
    prior = ("cd /tmp && printf 'source leaprc.protein.ff14SB\\nsource leaprc.water.tip3p\\n"
             "m = sequence { ACE ALA ALA ALA NME }\\nsolvatebox m TIP3PBOX 10.0\\n"
             "saveamberparm m system.prmtop system.inpcrd\\nquit\\n' > leap.in && tleap -f leap.in > /dev/null\n"
             "cp system.prmtop system.inpcrd /task/output/ && cp system.inpcrd /task/output/equilibrated.rst7\n"
             "echo '{\"disulfides\": [], \"protein_net_charge\": 0, \"n_na\": 0, \"n_cl\": 0, \"n_water\": 5000, "
             "\"box_angstrom\": [60.0, 60.0, 60.0]}' > /task/output/report.json\n")
    attacks = [
        Attack("prior_only_generic_system", "a valid solvated Amber system of an unrelated peptide with a constant "
               "report", 0.05, script=prior),
        Attack("salt_beyond_tolerance", "system built from the hidden clean model, then NaCl pairs removed until the "
               "concentration is >= 0.025 M below the SOP value", 0.75,
               script=run_truth + f"python3 /alegen_script/perturb.py --salt {salt} --min-deficit 0.025\n",
               script_files=truth_in),
        Attack("naive_reduced_cys", "cysteines left as CYS with thiol hydrogens, no S-S bonds", 0.75,
               script=naive("--ss none"), script_files=builder),
        Attack("naive_no_ions", "solvated but neither neutralised nor salted", 0.75,
               script=naive("--ions none"), script_files=builder),
    ]
    traps = []
    if level >= 2:
        attacks += [
            Attack("naive_crystal_solvent_kept", "crystal waters and Na+/Cl- sites kept (only the residues tleap "
                   "cannot load were deleted)", 0.75, script=naive("--xtal keep"), script_files=builder),
            Attack("naive_mse_dropped", "`grep ^ATOM`: every HETATM record removed, so the modified residue is "
                   "missing", 0.75, script=naive("--het atom-records-only"), script_files=builder),
            Attack("naive_sidechains_truncated", "residues with unresolved side-chain atoms simulated as the "
                   "deposited stub (ALA/GLY)", 0.75, script=naive("--incomplete truncate"), script_files=builder),
        ]
        traps += [
            Trap("hidden_disulfides", "No SSBOND/CONECT records and cysteines are named CYS; tleap then builds "
                 "reduced cysteines with clashing thiol hydrogens.",
                 "SG-SG distances of ~2.0 A; pdb4amber renames the residues to CYX and writes an sslink file; SOP 1.7",
                 "Approach / missing verification of intermediate results"),
            Trap("additives_and_waters", "Crystal waters, glycerol, sulfate/ethanediol and Na+/Cl- sites are in the "
                 "file; water and ions load into tleap without any error and silently become part of the system.",
                 "HETATM/HETNAM records; tleap 'Unknown residue: GOL/SO4/EDO' and saveamberparm 'FATAL: Atom ... does not "
                 "have a type'; SOP 1.8",
                 "Understanding / instructions in provided documents not applied"),
            Trap("modified_residue", "One methionine is deposited as selenomethionine (HETATM MSE, atom SE); deleting "
                 "HETATM records removes a residue from the chain.",
                 "MODRES and SEQRES records; HETATM lines inside the chain; tleap 'Unknown residue: MSE'; SOP 1.6",
                 "Understanding / domain knowledge gap"),
            Trap("missing_sidechains", "Several residues are truncated at CB or CG. tleap/PDBFixer rebuild them, but "
                 "a stub-preserving workflow simulates the wrong amino acid.",
                 "REMARK 470; pdb4amber 'misses n heavy atom(s)'; tleap 'Added missing heavy atom'; SOP 1.4",
                 "Approach / wrong strategy"),
        ]
    if level >= 3:
        attacks += [
            Attack("naive_gap_bonded", "no TER at the unresolved segment: tleap bonds the flanking residues", 0.75,
                   script=naive("--gap bonded", short=True), script_files=builder),
            Attack("naive_gap_charged_termini", "TER at the unresolved segment (pdb4amber default) but no caps: "
                   "charged termini inside the protein", 0.75, script=naive("--gap charged"), script_files=builder),
            Attack("naive_altloc_first", "first listed conformer kept (pdb4amber / OpenMM default) instead of the "
                   "highest-occupancy one", 0.75, script=naive("--altloc first"), script_files=builder),
            Attack("naive_internal_numbering", "disulfides reported in the builder's own sequential numbering", 0.75,
                   script=naive("--numbering internal"), script_files=builder),
        ]
        traps += [
            Trap("chain_break", "4-6 internal residues are unresolved. Without a TER tleap bonds across the gap; with "
                 "pdb4amber's TER the break gets charged termini. The SOP wants neutral caps and no modelling.",
                 "REMARK 465 and SEQRES versus ATOM records; C-N distance > 6 A; pdb4amber 'gap of ... A'; tleap 'There is "
                 "a bond of 13.7 angstroms between C and N atoms'; SOP 1.1-1.3",
                 "Approach / wrong strategy"),
            Trap("altlocs", "Two or three residues carry A/B conformers; for at least one of them the minor "
                 "conformer is listed first, so 'keep altloc A' picks the wrong one.",
                 "altLoc column and occupancies < 1.00; pdb4amber lists residues with alternate locations; tleap 'Atom "
                 "names in each residue should be unique' (first occurrence kept); SOP 1.5",
                 "Understanding / tool default taken for the requirement"),
            Trap("numbering_insertion_codes", "Residue numbers carry an offset and insertion codes (one on a "
                 "disulfide cysteine) and the chain is renamed; builders renumber from 1.",
                 "column 27 of the ATOM records; pdb4amber renum file; output contract of report.json",
                 "Output / contract not followed"),
        ]
    if level >= 2:
        attacks.append(Attack(
            "naive_raw_tleap", "raw file into tleap, deleting whatever makes it fail: first conformers, reduced "
            "cysteines, HETATM residues dropped, break bonded, crystal water/ions kept, neutralised only, internal "
            "numbering", 0.6, script=naive("--raw", short=True), script_files=builder))

    agent_df = dockerlib.agent_dockerfile(conda=AGENT_CONDA).replace(*_TRIM)
    verifier_df = dockerlib.verifier_dockerfile(pip=[], conda=VERIFIER_CONDA).replace(*_TRIM)
    return TaskBuild(
        task_id=f"md-system-repair-s{seed:04d}-l{level}", family="md_system_repair", seed=seed, level=level,
        title="Defective crystal structure to simulation-ready Amber system",
        summary="Repair a crystallographic protein model according to a group SOP, build the solvated, neutralised "
                "Amber system, equilibrate it briefly and report its key properties.",
        subdomain="3.1 Biomolecular Structure & Design",
        software=["AmberTools (tleap, pdb4amber)", "OpenMM", "PDBFixer", "ParmEd", "Python"],
        instruction_md=INSTRUCTION,
        agent_dockerfile=agent_df, verifier_dockerfile=verifier_df,
        input_files={"structure.pdb": planted["pdb"], "simulation_sop.md": sop_md},
        solution_files={"solve.sh": ASSETS / "solve.sh", "build_tleap.py": ASSETS / "build_tleap.py",
                        "solve_alt.py": ASSETS / "solve_alt.py"},
        tests_files={"score.py": ASSETS / "score.py", "reference/truth.json": json.dumps(truth)},
        traps=traps, attacks=attacks, probes=probes,
        capabilities=["structure QC before modelling", "reading an SOP against a concrete data file",
                      "tool defaults versus stated requirements", "verifying a built system before hand-over",
                      "exact output contracts"],
        judge_explanation_md=JUDGE, pass_threshold=0.9,
        oracle_steps=9, oracle_runtime_min=4.0, cpus=4, memory_mb=8000, storage_mb=10000,
        extra_metadata={"source_entry_blinded": True},
    )
