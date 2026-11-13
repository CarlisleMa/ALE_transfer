#!/usr/bin/env python3
"""Reference solution: text-level repair -> pdb4amber -> tleap (ff14SB / TIP3P) -> OpenMM equilibration.

The non-default flags replay plausible-but-wrong preparations for the validator's attack suite;
the oracle (solve.sh) uses the defaults. Parameters of the target system are read from the SOP
values hard-wired below by the generator-independent defaults (padding, salt) or given as flags.
"""
import argparse
import json
import math
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--input", default="/task/input")
ap.add_argument("--output", default="/task/output")
ap.add_argument("--work", default="/tmp/mdprep_work")
ap.add_argument("--padding", type=float, default=None, help="Angstrom; default: read from the SOP")
ap.add_argument("--salt", type=float, default=None, help="mol/L; default: read from the SOP")
ap.add_argument("--altloc", choices=["occupancy", "first"], default="occupancy")
ap.add_argument("--ss", choices=["bond", "none"], default="bond")
ap.add_argument("--gap", choices=["cap", "charged", "bonded"], default="cap")
ap.add_argument("--het", choices=["convert", "atom-records-only"], default="convert")
ap.add_argument("--xtal", choices=["drop", "keep"], default="drop")
ap.add_argument("--ions", choices=["sop", "neutral", "none"], default="sop")
ap.add_argument("--numbering", choices=["input", "internal"], default="input")
ap.add_argument("--incomplete", choices=["rebuild", "truncate"], default="rebuild")
ap.add_argument("--raw", action="store_true", help="careless end-to-end: raw file straight into tleap")
ap.add_argument("--no-equil", action="store_true")
a = ap.parse_args()
if a.raw:
    a.altloc, a.ss, a.gap, a.het, a.xtal, a.ions, a.numbering = (
        "first", "none", "bonded", "atom-records-only", "keep", "neutral", "internal")

inp, out, work = Path(a.input), Path(a.output), Path(a.work)
out.mkdir(parents=True, exist_ok=True)
work.mkdir(parents=True, exist_ok=True)
os.chdir(work)

sop = (inp / "simulation_sop.md").read_text()
m = re.search(r"padding[^\n]*?(\d+(?:\.\d+)?)\s*(?:Å|A\b|angstrom)", sop, re.I)
padding = a.padding if a.padding is not None else float(m.group(1))
m = re.search(r"(\d+\.\d+)\s*(?:M|mol/L)\s+NaCl", sop)
salt = a.salt if a.salt is not None else float(m.group(1))


def sh(cmd, **kw):
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        raise SystemExit(f"command failed: {' '.join(map(str, cmd))}\n{r.stdout[-1500:]}\n{r.stderr[-1500:]}")
    return r


# ------------------------------------------------------------------ 1. read the entry
AA = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL MSE".split())
atoms = []           # dict(rec, name, alt, resn, chain, num, icode, xyz, occ, b, elem)
for line in (inp / "structure.pdb").read_text().splitlines():
    if line[:6] in ("ATOM  ", "HETATM"):
        atoms.append(dict(rec=line[:6], name=line[12:16].strip(), alt=line[16], resn=line[17:20].strip(),
                          chain=line[21], num=int(line[22:26]), icode=line[26].strip(),
                          xyz=np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
                          occ=float(line[54:60]), b=float(line[60:66]),
                          elem=(line[76:78].strip() or line[12:14].strip()).upper()))
missing_465 = []
for line in (inp / "structure.pdb").read_text().splitlines():
    mm = re.match(r"REMARK 465\s+([A-Z]{3})\s+(\S)\s+(-?\d+)(\w?)\s*$", line)
    if mm:
        missing_465.append((mm.group(1), f"{mm.group(3)}{mm.group(4)}"))

# residues in file order
residues, key = [], None
for at in atoms:
    k = (at["chain"], at["num"], at["icode"], at["resn"])
    if k != key:
        residues.append({"resn": at["resn"], "label": f"{at['num']}{at['icode']}", "atoms": [], "rec": at["rec"]})
        key = k
    residues[-1]["atoms"].append(at)

if a.het == "atom-records-only":
    # careless: `grep '^ATOM'` -- drops every HETATM record, the modified residue included
    protein = [r for r in residues if r["rec"] == "ATOM  " and r["resn"] in AA]
else:
    protein = [r for r in residues if r["resn"] in AA]
    for r in protein:
        if r["resn"] == "MSE":                       # selenomethionine -> methionine
            r["resn"] = "MET"
            for at in r["atoms"]:
                at["rec"] = "ATOM  "
                if at["name"] == "SE":
                    at["name"], at["elem"] = "SD", "S"
xtal = [r for r in residues if r["resn"] in ("HOH", "WAT", "NA", "CL")] if a.xtal == "keep" else []

NHEAVY = {"ALA": 5, "ARG": 11, "ASN": 8, "ASP": 8, "CYS": 6, "GLN": 9, "GLU": 9, "GLY": 4, "HIS": 10, "ILE": 8,
          "LEU": 8, "LYS": 9, "MET": 8, "PHE": 11, "PRO": 7, "SER": 6, "THR": 7, "TRP": 14, "TYR": 12, "VAL": 7}
if a.incomplete == "truncate":
    # careless: residues with unresolved side-chain atoms are simulated as the stub that was deposited
    for r in protein:
        names = {at["name"] for at in r["atoms"] if at["alt"] in (" ", "A")} - {"OXT"}
        if r["resn"] in NHEAVY and len(names) < NHEAVY[r["resn"]]:
            r["resn"] = "ALA" if "CB" in names else "GLY"
            r["atoms"] = [at for at in r["atoms"] if at["name"] in ("N", "CA", "C", "O", "CB", "OXT")]

# alternate conformers
for r in protein:
    alts = sorted({at["alt"] for at in r["atoms"] if at["alt"] != " "})
    if not alts:
        continue
    if a.altloc == "first":
        best = alts[0]
    else:
        occ = {al: np.mean([at["occ"] for at in r["atoms"] if at["alt"] == al]) for al in alts}
        best = max(alts, key=lambda al: occ[al])
    r["atoms"] = [at for at in r["atoms"] if at["alt"] in (" ", best)]
    for at in r["atoms"]:
        at["alt"], at["occ"] = " ", 1.0


# ------------------------------------------------------------------ 2. chain breaks
def get(r, name):
    for at in r["atoms"]:
        if at["name"] == name:
            return at["xyz"]
    return None


breaks = []          # index i: break between protein[i] and protein[i+1]
for i in range(len(protein) - 1):
    c, n = get(protein[i], "C"), get(protein[i + 1], "N")
    if np.linalg.norm(c - n) > 1.8:
        breaks.append(i)


def place(p1, p2, p3, bond, angle, torsion):
    """Position of an atom bonded to p3 with angle (p2,p3,new) and dihedral (p1,p2,p3,new)."""
    ang, tor = math.radians(angle), math.radians(torsion)
    bc = (p3 - p2) / np.linalg.norm(p3 - p2)
    nrm = np.cross(p2 - p1, bc)
    nrm /= np.linalg.norm(nrm)
    mvec = np.cross(nrm, bc)
    d = np.array([-bond * math.cos(ang), bond * math.sin(ang) * math.cos(tor), bond * math.sin(ang) * math.sin(tor)])
    return p3 + d[0] * bc + d[1] * mvec + d[2] * nrm


def cap_atom(name, resn, xyz, elem):
    return dict(rec="ATOM  ", name=name, alt=" ", resn=resn, chain="A", num=0, icode="", xyz=xyz, occ=1.0,
                b=30.0, elem=elem)


model = []           # ordered: residues and "TER" markers
allxyz = np.array([at["xyz"] for r in protein for at in r["atoms"]])
for i, r in enumerate(protein):
    if i - 1 in breaks and a.gap == "cap":
        # acetyl carbonyl carbon; tleap builds CH3/O/H from the template. Pick the phi that keeps the
        # cap furthest from the rest of the model.
        n_, ca, c = get(r, "N"), get(r, "CA"), get(r, "C")
        best = None
        for phi in (-60, -90, -120, -150, 180, 60):
            pos = place(c, ca, n_, 1.335, 121.9, phi)
            ch3 = place(ca, n_, pos, 1.522, 116.6, 180.0)
            d = min(np.min(np.linalg.norm(allxyz - pos, axis=1)[np.linalg.norm(allxyz - n_, axis=1) > 2.6]),
                    np.min(np.linalg.norm(allxyz - ch3, axis=1)[np.linalg.norm(allxyz - n_, axis=1) > 2.6]))
            if best is None or d > best[0]:
                best = (d, pos)
        model.append({"resn": "ACE", "label": None, "atoms": [cap_atom("C", "ACE", best[1], "C")]})
    model.append(r)
    if i in breaks:
        if a.gap == "cap":
            ca, c, o = get(r, "CA"), get(r, "C"), get(r, "O")
            u = (o - c) / np.linalg.norm(o - c) + (ca - c) / np.linalg.norm(ca - c)
            npos = c - 1.335 * u / np.linalg.norm(u)
            model.append({"resn": "NME", "label": None, "atoms": [cap_atom("N", "NME", npos, "N")]})
        if a.gap in ("cap", "charged"):
            model.append("TER")

# ------------------------------------------------------------------ 3. disulfides from geometry
seq_index, k = {}, 0             # id(residue dict) -> 1-based residue number in tleap
for item in model:
    if item != "TER":
        k += 1
        seq_index[id(item)] = k
cys = [r for r in protein if r["resn"] == "CYS" and get(r, "SG") is not None]
ss_pairs = []
for x in range(len(cys)):
    for y in range(x + 1, len(cys)):
        if np.linalg.norm(get(cys[x], "SG") - get(cys[y], "SG")) < 2.5:
            ss_pairs.append((cys[x], cys[y]))
if a.ss == "bond":
    for x, y in ss_pairs:
        x["resn"] = y["resn"] = "CYX"


def fmt(serial, at, resn, num):
    name = at["name"]
    nm = f"{name:<4}" if (len(name) >= 4 or len(at["elem"]) == 2) else f" {name:<3}"
    x, y, z = at["xyz"]
    return (f"ATOM  {serial:>5} {nm} {resn:>3} A{num:>4}    {x:8.3f}{y:8.3f}{z:8.3f}{at['occ']:6.2f}{at['b']:6.2f}"
            f"          {at['elem']:>2}")


lines, serial, num = [], 0, 0
for item in model:
    if item == "TER":
        lines.append("TER")
        continue
    num += 1
    for at in item["atoms"]:
        serial += 1
        lines.append(fmt(serial, at, item["resn"], num))
lines.append("TER")
for r in xtal:
    num += 1
    resn = {"HOH": "WAT", "NA": "Na+", "CL": "Cl-"}.get(r["resn"], r["resn"])
    for at in r["atoms"]:
        serial += 1
        at = dict(at, name={"NA": "Na+", "CL": "Cl-"}.get(at["name"], at["name"]))
        lines.append(fmt(serial, at, resn, num))
    lines.append("TER")
lines.append("END")
Path("prepared.pdb").write_text("\n".join(lines) + "\n")

# pdb4amber as an independent consistency check of the hand-edited file (gaps, disulfides, residues)
r4 = subprocess.run(["pdb4amber", "-i", "prepared.pdb", "-o", "p4a.pdb"], capture_output=True, text=True)
Path("pdb4amber.log").write_text(r4.stdout + r4.stderr)

# ------------------------------------------------------------------ 4. tleap
n_prot_res = sum(1 for item in model if item != "TER")
bonds = "".join(f"bond m.{seq_index[id(x)]}.SG m.{seq_index[id(y)]}.SG\n" for x, y in ss_pairs) if a.ss == "bond" else ""
head = ("source leaprc.protein.ff14SB\nsource leaprc.water.tip3p\n"
        "m = loadpdb prepared.pdb\n" + bonds)
# closeness 1.0 is the tleap default; the buffer is counted from the solute's van der Waals surface box
solv = f"solvatebox m TIP3PBOX {padding:.2f}\n"
Path("leap1.in").write_text(head + solv + "charge m\nsavepdb m probe.pdb\nquit\n")
r1 = sh(["tleap", "-f", "leap1.in"])
Path("leap1.log").write_text(r1.stdout)
q = float(re.search(r"Total (?:unperturbed|perturbed) charge:\s+(-?\d+\.\d+)", r1.stdout).group(1))
n_wat0 = sum(1 for ln in Path("probe.pdb").read_text().splitlines()
             if ln.startswith(("ATOM", "HETATM")) and ln[17:20] == "WAT" and ln[12:16].strip() == "O")
qi = int(round(q))
n_neut = abs(qi)
if a.ions == "sop":
    n_pairs = int(round(salt / 55.5 * (n_wat0 - n_neut) / (1.0 + 2.0 * salt / 55.5)))
else:
    n_pairs = 0
n_na = n_pairs + (n_neut if qi < 0 else 0)
n_cl = n_pairs + (n_neut if qi > 0 else 0)
ions = ""
if a.ions != "none":
    if n_na:
        ions += f"addIonsRand m Na+ {n_na}\n"
    if n_cl:
        ions += f"addIonsRand m Cl- {n_cl}\n"
Path("leap2.in").write_text(head + solv + ions + "charge m\n"
                            "saveamberparm m system.prmtop system.inpcrd\nquit\n")
r2 = sh(["tleap", "-f", "leap2.in"])
Path("leap2.log").write_text(r2.stdout)
if not Path("system.prmtop").exists() or Path("system.prmtop").stat().st_size == 0:
    raise SystemExit("tleap did not write a topology:\n" + r2.stdout[-3000:])

# ------------------------------------------------------------------ 5. equilibration (OpenMM)
import openmm as mm                                  # noqa: E402
import openmm.app as app                             # noqa: E402
import openmm.unit as u                              # noqa: E402
import parmed as pmd                                 # noqa: E402

prm = app.AmberPrmtopFile("system.prmtop")
crd = app.AmberInpcrdFile("system.inpcrd")
box_final = None
if not a.no_equil:
    system = prm.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=0.8 * u.nanometer,
                              constraints=app.HBonds, rigidWater=True)
    integ = mm.LangevinMiddleIntegrator(300 * u.kelvin, 1.0 / u.picosecond, 0.002 * u.picoseconds)
    integ.setRandomNumberSeed(20240611)
    plat = mm.Platform.getPlatformByName("CPU")
    sim = app.Simulation(prm.topology, system, integ, plat)
    sim.context.setPositions(crd.positions)
    sim.context.setPeriodicBoxVectors(*crd.boxVectors)
    total_mass = sum(system.getParticleMass(i).value_in_unit(u.dalton) for i in range(system.getNumParticles()))
    import time as _t
    _t0 = _t.time()
    sim.minimizeEnergy(maxIterations=int(os.environ.get("EQ_MIN_ITER", "200")))
    good = sim.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=False)
    print(f"minimisation {_t.time() - _t0:.0f} s", file=sys.stderr)
    try:
        sim.context.setVelocitiesToTemperature(300 * u.kelvin, 7)
        sim.step(int(os.environ.get("EQ_NVT_STEPS", "250")))                  # NVT
        baro = mm.MonteCarloBarostat(1.0 * u.bar, 300 * u.kelvin, int(os.environ.get("EQ_BARO_FREQ", "5")))
        baro.setRandomNumberSeed(11)
        system.addForce(baro)
        sim.context.reinitialize(preserveState=True)
        left = int(os.environ.get("EQ_NPT_STEPS", "3000"))                     # NPT
        while left > 0:
            sim.step(min(500, left))
            left -= 500
            st = sim.context.getState(getPositions=True, getVelocities=True, getEnergy=True,
                                      enforcePeriodicBox=False)
            if not math.isfinite(st.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)):
                raise ValueError("energy is not finite")
            good = st
            vol = st.getPeriodicBoxVolume().value_in_unit(u.angstrom ** 3)
            print(f"t={_t.time() - _t0:.0f}s npt steps_left={left} volume={vol:.0f} A^3 density={total_mass / vol * 1.66054:.4f} g/cc",
                  file=sys.stderr)
    except Exception as e:          # a careless run that blows up still ships its last finite frame
        print("equilibration stopped early:", e, file=sys.stderr)
    pos = np.array(good.getPositions(asNumpy=True).value_in_unit(u.angstrom))
    vel = np.array(good.getVelocities(asNumpy=True).value_in_unit(u.angstrom / u.picosecond))
    bv = good.getPeriodicBoxVectors(asNumpy=True).value_in_unit(u.angstrom)
    box_final = [float(bv[0][0]), float(bv[1][1]), float(bv[2][2])]
    rst = pmd.amber.Rst7(natom=len(pos))
    rst.coordinates = pos
    rst.vels = vel
    rst.box = box_final + [90.0, 90.0, 90.0]
    rst.time = 0.0
    rst.write("equilibrated.rst7")

# ------------------------------------------------------------------ 6. report (recomputed from the topology)
parm = pmd.load_file("system.prmtop")
names = [r.name for r in parm.residues]
n_w = sum(1 for x in names if x in ("WAT", "HOH"))
prot_q = sum(at.charge for at in parm.atoms if at.residue.name not in ("WAT", "HOH", "Na+", "Cl-"))
if a.numbering == "input":
    pairs = [[x["label"], y["label"]] for x, y in ss_pairs]
else:
    pairs = [[str(seq_index[id(x)]), str(seq_index[id(y)])] for x, y in ss_pairs]
report = {
    "disulfides": pairs,
    "protein_net_charge": int(round(prot_q)),
    "n_na": sum(1 for x in names if x == "Na+"),
    "n_cl": sum(1 for x in names if x == "Cl-"),
    "n_water": n_w,
    "box_angstrom": [round(v, 3) for v in (box_final or [float(x) for x in parm.box[:3]])],
}
for f in ("system.prmtop", "system.inpcrd"):
    (out / f).write_bytes(Path(f).read_bytes())
if box_final is not None:
    (out / "equilibrated.rst7").write_bytes(Path("equilibrated.rst7").read_bytes())
else:   # skipped equilibration: the as-built coordinates are shipped under the restart name
    (out / "equilibrated.rst7").write_bytes(Path("system.inpcrd").read_bytes())
(out / "report.json").write_text(json.dumps(report, indent=1) + "\n")
print(json.dumps({"padding": padding, "salt": salt, "charge_before_ions": q, "n_wat0": n_wat0,
                  "breaks": [protein[i]["label"] for i in breaks], "report": report}))
