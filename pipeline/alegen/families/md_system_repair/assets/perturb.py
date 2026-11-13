#!/usr/bin/env python3
"""Validator helper (never shipped to the agent): take a finished submission in /task/output and
remove NaCl pairs / round the reported box, to probe the verifier's tolerances from both sides."""
import argparse
import json
from pathlib import Path

import numpy as np
import parmed as pmd

ap = argparse.ArgumentParser()
ap.add_argument("--output", default="/task/output")
ap.add_argument("--salt", type=float, required=True, help="SOP concentration, mol/L")
ap.add_argument("--max-deficit", type=float, default=None,
                help="remove pairs as long as the concentration stays within this deficit")
ap.add_argument("--min-deficit", type=float, default=None,
                help="remove pairs until the concentration is at least this far below the SOP value")
ap.add_argument("--round-box", type=int, default=None)
ap.add_argument("--shift-box", type=float, default=0.0)
a = ap.parse_args()
out = Path(a.output)

parm = pmd.amber.AmberParm(str(out / "system.prmtop"))
na = [r for r in parm.residues if r.name in ("Na+", "NA")]
cl = [r for r in parm.residues if r.name in ("Cl-", "CL")]
n_wat = sum(1 for r in parm.residues if r.name in ("WAT", "HOH"))
pairs = min(len(na), len(cl))


def conc(p):
    return p / n_wat * 55.5


k = 0
if a.max_deficit is not None:
    while pairs - k - 1 > 0 and a.salt - conc(pairs - k - 1) <= a.max_deficit:
        k += 1
if a.min_deficit is not None:
    while pairs - k > 0 and a.salt - conc(pairs - k) < a.min_deficit:
        k += 1
drop = [at.idx for r in na[len(na) - k:] + cl[len(cl) - k:] for at in r.atoms] if k else []
keep = np.ones(len(parm.atoms), dtype=bool)
keep[drop] = False
if k:
    parm.strip("@" + ",".join(str(i + 1) for i in sorted(drop)))
    parm.write_parm(str(out / "system.prmtop"))
for name in ("system.inpcrd", "equilibrated.rst7"):
    old = pmd.amber.Rst7(str(out / name))
    xyz = np.asarray(old.coordinates).reshape(-1, 3)[keep]
    new = pmd.amber.Rst7(natom=int(keep.sum()))
    new.coordinates = xyz
    if old.hasvels:
        new.vels = np.asarray(old.vels).reshape(-1, 3)[keep]
    new.box = list(old.box)
    new.time = 0.0
    new.write(str(out / name))
rep = json.loads((out / "report.json").read_text())
rep["n_na"] -= k
rep["n_cl"] -= k
if a.round_box is not None:
    rep["box_angstrom"] = [round(v, a.round_box) for v in rep["box_angstrom"]]
if a.shift_box:
    rep["box_angstrom"] = [v + a.shift_box for v in rep["box_angstrom"]]
(out / "report.json").write_text(json.dumps(rep) + "\n")
print(json.dumps({"removed_pairs": k, "pairs_left": pairs - k, "conc": round(conc(pairs - k), 4)}))
