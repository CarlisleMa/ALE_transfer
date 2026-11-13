#!/usr/bin/env python3
"""Independent second solution: no pdb4amber, no tleap.

hand-written PDB repair (conformer selection, modified residue, chain-break caps with explicit
geometry) -> PDBFixer (missing heavy atoms) -> OpenMM Modeller (hydrogens, TIP3P box, ions) ->
ff14SB system from the OpenMM force-field XML -> ParmEd export to Amber prmtop/inpcrd ->
equilibration with a different integrator/barostat schedule than the reference.
"""
import io
import json
import math
import os
import re
from pathlib import Path

import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as u
import parmed as pmd
from pdbfixer import PDBFixer

inp, out = Path(os.environ.get("ALT_INPUT", "/task/input")), Path(os.environ.get("ALT_OUTPUT", "/task/output"))
out.mkdir(parents=True, exist_ok=True)
sop = (inp / "simulation_sop.md").read_text()
PAD = float(re.search(r"padding[^\n]*?(\d+(?:\.\d+)?)\s*(?:Å|A\b|angstrom)", sop, re.I).group(1))
SALT = float(re.search(r"(\d+\.\d+)\s*(?:M|mol/L)\s+NaCl", sop).group(1))
STD = "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split()

# ---------------------------------------------------------------- parse, keep the polypeptide only
recs = []
for ln in (inp / "structure.pdb").read_text().splitlines():
    if ln.startswith(("ATOM", "HETATM")):
        recs.append({"name": ln[12:16].strip(), "alt": ln[16], "resn": ln[17:20].strip(),
                     "rid": (ln[21], ln[22:27]), "label": ln[22:27].strip(),
                     "xyz": np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]),
                     "occ": float(ln[54:60]), "elem": (ln[76:78].strip() or ln[12:14].strip()).upper()})
chain = []
for r in recs:
    if r["resn"] == "MSE":                                   # selenomethionine is simulated as Met
        r["resn"] = "MET"
        if r["name"] == "SE":
            r["name"], r["elem"] = "SD", "S"
    if r["resn"] not in STD or r["elem"] in ("H", "D"):
        continue                                             # waters, additives, ions: not simulated
    if not chain or chain[-1]["rid"] != r["rid"]:
        chain.append({"rid": r["rid"], "resn": r["resn"], "label": r["label"], "atoms": []})
    chain[-1]["atoms"].append(r)
for res in chain:                                            # one conformer per residue: highest occupancy
    tot = {}
    for at in res["atoms"]:
        if at["alt"] != " ":
            tot.setdefault(at["alt"], []).append(at["occ"])
    if tot:
        keep = max(sorted(tot), key=lambda k: sum(tot[k]) / len(tot[k]))
        res["atoms"] = [at for at in res["atoms"] if at["alt"] in (" ", keep)]
    res["xyz"] = {at["name"]: at["xyz"] for at in res["atoms"]}


def nerf(a, b, c, bond, angle, torsion):
    """Natural-extension-reference-frame placement of a new atom bonded to c."""
    bc = c - b
    bc /= np.linalg.norm(bc)
    n = np.cross(b - a, bc)
    n /= np.linalg.norm(n)
    frame = np.array([bc, np.cross(n, bc), n])
    th, ph = math.radians(180.0 - angle), math.radians(torsion)
    local = bond * np.array([math.cos(th), math.sin(th) * math.cos(ph), math.sin(th) * math.sin(ph)])
    return c + local @ frame


# ---------------------------------------------------------------- fragments and blocking groups
frags = [[chain[0]]]
for prev, cur in zip(chain, chain[1:]):
    if np.linalg.norm(prev["xyz"]["C"] - cur["xyz"]["N"]) > 1.8:
        frags.append([])
    frags[-1].append(cur)
everything = np.array([at["xyz"] for res in chain for at in res["atoms"]])


def clearance(points, anchor):
    far = np.linalg.norm(everything - anchor, axis=1) > 2.7
    return min(np.min(np.linalg.norm(everything[far] - p, axis=1)) for p in points)


pieces = []          # (resname, label or None, {atom: xyz})
for fi, frag in enumerate(frags):
    if fi > 0:       # N-terminal side of a break: acetyl, all heavy atoms placed explicitly
        r = frag[0]
        n_, ca, c = r["xyz"]["N"], r["xyz"]["CA"], r["xyz"]["C"]
        best = None
        for phi in range(-180, 180, 30):
            cc = nerf(c, ca, n_, 1.335, 121.9, phi)
            ch3 = nerf(ca, n_, cc, 1.522, 116.6, 180.0)
            oo = nerf(ca, n_, cc, 1.229, 122.9, 0.0)
            sc = clearance([cc, ch3, oo], n_)
            if best is None or sc > best[0]:
                best = (sc, {"CH3": ch3, "C": cc, "O": oo})
        pieces.append(("ACE", None, best[1]))
    for r in frag:
        pieces.append((r["resn"], r["label"], r["xyz"]))
    if fi < len(frags) - 1:   # C-terminal side of a break: N-methylamide
        r = frag[-1]
        n_, ca, c, o = r["xyz"]["N"], r["xyz"]["CA"], r["xyz"]["C"], r["xyz"]["O"]
        psi_o = math.degrees(math.atan2(
            np.dot(np.cross(np.cross(n_ - ca, c - ca), np.cross(c - ca, o - c)), (c - ca) / np.linalg.norm(c - ca)),
            np.dot(np.cross(n_ - ca, c - ca), np.cross(c - ca, o - c))))
        nn = nerf(n_, ca, c, 1.335, 116.6, psi_o + 180.0)
        cm = nerf(ca, c, nn, 1.449, 121.9, 180.0)
        pieces.append(("NME", None, {"N": nn, "C": cm}))
    pieces.append(None)       # chain terminator

lines, serial, resnum, cid = [], 0, 0, 0
labels = {}                   # (chain index, residue index in chain) is not needed: labels by running residue number
for p in pieces:
    if p is None:
        lines.append("TER")
        cid += 1
        continue
    resn, label, xyz = p
    resnum += 1
    labels[resnum] = label
    for name, pos in xyz.items():
        serial += 1
        el = "S" if name in ("SG", "SD") else name[0]
        lines.append(f"ATOM  {serial:>5} {name:<4}".replace(f" {name:<4}", f"  {name:<3}" if len(name) < 4 else f" {name}")
                     + f" {resn:>3} {chr(65 + cid)}{resnum:>4}    {pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}  1.00 20.00"
                     + f"          {el:>2}")
lines.append("END")
repaired = "\n".join(lines) + "\n"
Path("/tmp/alt_repaired.pdb").write_text(repaired)

# ---------------------------------------------------------------- PDBFixer: missing heavy atoms only
fixer = PDBFixer(pdbfile=io.StringIO(repaired))
fixer.missingResidues = {}                 # unresolved residues are not modelled
fixer.findMissingAtoms()
fixer.addMissingAtoms(seed=3)
ff = app.ForceField("amber14/protein.ff14SB.xml", "amber14/tip3p.xml")
mod = app.Modeller(fixer.topology, fixer.positions)
mod.addHydrogens(ff, pH=7.0)

# disulfides: OpenMM bonds SG atoms closer than 0.3 nm when it reads the PDB; record them from the topology
res_of = {}
for ch in mod.topology.chains():
    for r in ch.residues():
        res_of[r.index] = r
order = {r.index: k + 1 for k, r in enumerate(mod.topology.residues())}
ss = []
for b in mod.topology.bonds():
    if b[0].name == "SG" and b[1].name == "SG":
        ss.append(sorted([labels[order[b[0].residue.index]], labels[order[b[1].residue.index]]],
                         key=lambda s: (int(re.match(r"-?\d+", s).group()), s)))

# ---------------------------------------------------------------- solvent box and ions
pos = np.array(mod.positions.value_in_unit(u.nanometer))
ext = pos.max(axis=0) - pos.min(axis=0)
box = ext + 2.0 * (PAD / 10.0) + 0.05                         # nm; rectangular, padding on every side
mod.addSolvent(ff, model="tip3p", boxSize=mm.Vec3(*box) * u.nanometer, neutralize=False)
system = ff.createSystem(mod.topology, nonbondedMethod=app.PME, nonbondedCutoff=0.9 * u.nanometer,
                         constraints=None, rigidWater=False)
nb = [f for f in system.getForces() if isinstance(f, mm.NonbondedForce)][0]
n_solute = sum(1 for at in mod.topology.atoms() if at.residue.name != "HOH")
q = sum(nb.getParticleParameters(i)[0].value_in_unit(u.elementary_charge) for i in range(n_solute))
qi = int(round(q))
n_wat = sum(1 for r in mod.topology.residues() if r.name == "HOH")
n_pairs = int(round(SALT * (n_wat - abs(qi)) / (55.5 + 2 * SALT)))
n_na, n_cl = n_pairs + max(0, -qi), n_pairs + max(0, qi)
# replace randomly chosen waters (at least 5 A from the solute and from each other) by ions
rng = np.random.default_rng(12)
allpos = np.array(mod.positions.value_in_unit(u.angstrom))
solute_xyz = allpos[:n_solute]
wat = [r for r in mod.topology.residues() if r.name == "HOH"]
chosen, chosen_xyz = [], []
for wi in rng.permutation(len(wat)):
    o = [at for at in wat[wi].atoms() if at.name == "O"][0]
    p = allpos[o.index]
    if np.min(np.linalg.norm(solute_xyz - p, axis=1)) < 5.0:
        continue
    if chosen_xyz and np.min(np.linalg.norm(np.array(chosen_xyz) - p, axis=1)) < 5.0:
        continue
    chosen.append(wat[wi])
    chosen_xyz.append(p)
    if len(chosen) == n_na + n_cl:
        break
boxvec = mod.topology.getPeriodicBoxVectors()
mod.delete(chosen)
ion_top = app.Topology()
ich = ion_top.addChain()
for k, p in enumerate(chosen_xyz):
    name, el = ("NA", app.element.sodium) if k < n_na else ("CL", app.element.chlorine)
    ion_top.addAtom(name, el, ion_top.addResidue(name, ich))
mod.add(ion_top, [mm.Vec3(*p) for p in chosen_xyz] * u.angstrom)
mod.topology.setPeriodicBoxVectors(boxvec)

# ---------------------------------------------------------------- ff14SB system -> Amber files via ParmEd
system = ff.createSystem(mod.topology, nonbondedMethod=app.PME, nonbondedCutoff=0.9 * u.nanometer,
                         constraints=None, rigidWater=False)
struct = pmd.openmm.load_topology(mod.topology, system, xyz=mod.positions)
struct.box = [box[0] * 10, box[1] * 10, box[2] * 10, 90.0, 90.0, 90.0]
# Amber treats TIP3P as a rigid triangle: add the H-H bond the Amber water library carries
hh = pmd.BondType(553.0, 1.5136, list=struct.bond_types)
struct.bond_types.append(hh)
for res in struct.residues:
    if res.name == "HOH":
        hs = [at for at in res.atoms if at.atomic_number == 1]
        struct.bonds.append(pmd.Bond(hs[0], hs[1], type=hh))
struct.save(str(out / "system.prmtop"), overwrite=True)
struct.save(str(out / "system.inpcrd"), overwrite=True)

# ---------------------------------------------------------------- equilibration from the written Amber files
prm = app.AmberPrmtopFile(str(out / "system.prmtop"))
crd = app.AmberInpcrdFile(str(out / "system.inpcrd"))
sysm = prm.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=0.8 * u.nanometer, constraints=app.HBonds,
                        rigidWater=True)
integ = mm.LangevinIntegrator(300 * u.kelvin, 2.0 / u.picosecond, 0.002 * u.picoseconds)
integ.setRandomNumberSeed(99)
sim = app.Simulation(prm.topology, sysm, integ, mm.Platform.getPlatformByName("CPU"))
sim.context.setPositions(crd.positions)
sim.context.setPeriodicBoxVectors(*crd.boxVectors)
sim.minimizeEnergy(maxIterations=int(os.environ.get("EQ_MIN_ITER", "400")))
sim.context.setVelocitiesToTemperature(300 * u.kelvin, 5)
sim.step(int(os.environ.get("EQ_NVT_STEPS", "250")))
baro = mm.MonteCarloBarostat(1.0 * u.atmosphere, 300 * u.kelvin, 20)
baro.setRandomNumberSeed(3)
sysm.addForce(baro)
sim.context.reinitialize(preserveState=True)
sim.step(int(os.environ.get("EQ_NPT_STEPS", "1500")))
st = sim.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=False)
bv = st.getPeriodicBoxVectors(asNumpy=True).value_in_unit(u.angstrom)
final_box = [float(bv[0][0]), float(bv[1][1]), float(bv[2][2])]
rst = pmd.amber.Rst7(natom=sysm.getNumParticles())
rst.coordinates = st.getPositions(asNumpy=True).value_in_unit(u.angstrom)
rst.vels = st.getVelocities(asNumpy=True).value_in_unit(u.angstrom / u.picosecond)
rst.box = final_box + [90.0, 90.0, 90.0]
rst.time = 0.0
rst.write(str(out / "equilibrated.rst7"))

# ---------------------------------------------------------------- report
names = [r.name for r in mod.topology.residues()]
report = {"disulfides": sorted(ss), "protein_net_charge": qi,
          "n_na": names.count("NA"), "n_cl": names.count("CL"), "n_water": names.count("HOH"),
          "box_angstrom": [round(v, 3) for v in final_box]}
(out / "report.json").write_text(json.dumps(report, indent=1) + "\n")
print(json.dumps(report))
