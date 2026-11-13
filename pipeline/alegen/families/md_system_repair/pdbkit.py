"""PDB handling for the md_system_repair generator: read a complete RCSB entry (the hidden truth),
plant defects, and write the defective file the agent sees. Pure python + numpy (host side)."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

# heavy atoms of the 20 standard amino acids (used only to assert that the source entry is complete)
HEAVY = {
    "ALA": "N CA C O CB", "ARG": "N CA C O CB CG CD NE CZ NH1 NH2", "ASN": "N CA C O CB CG OD1 ND2",
    "ASP": "N CA C O CB CG OD1 OD2", "CYS": "N CA C O CB SG", "GLN": "N CA C O CB CG CD OE1 NE2",
    "GLU": "N CA C O CB CG CD OE1 OE2", "GLY": "N CA C O", "HIS": "N CA C O CB CG ND1 CD2 CE1 NE2",
    "ILE": "N CA C O CB CG1 CG2 CD1", "LEU": "N CA C O CB CG CD1 CD2", "LYS": "N CA C O CB CG CD CE NZ",
    "MET": "N CA C O CB CG SD CE", "PHE": "N CA C O CB CG CD1 CD2 CE1 CE2 CZ", "PRO": "N CA C O CB CG CD",
    "SER": "N CA C O CB OG", "THR": "N CA C O CB OG1 CG2", "TRP": "N CA C O CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2",
    "TYR": "N CA C O CB CG CD1 CD2 CE1 CE2 CZ OH", "VAL": "N CA C O CB CG1 CG2",
}
HEAVY = {k: v.split() for k, v in HEAVY.items()}
CHARGE = {"ASP": -1, "GLU": -1, "LYS": 1, "ARG": 1}
GAMMA = {"SER": "OG", "THR": "OG1", "ILE": "CG1", "LEU": "CG", "ASN": "CG", "ASP": "CG", "GLN": "CG",
         "GLU": "CG", "LYS": "CG", "MET": "CG", "ARG": "CG", "PHE": "CG", "TYR": "CG", "HIS": "CG",
         "TRP": "CG", "CYS": "SG", "VAL": "CG1", "PRO": "CG"}

# idealised hetero groups (Angstrom), centred near the origin
HET_TEMPLATES = {
    "GOL": [("C1", "C", (-1.249, 0.521, 0.250)), ("O1", "O", (-2.400, -0.296, 0.050)),
            ("C2", "C", (0.000, -0.285, -0.050)), ("O2", "O", (0.020, -1.420, 0.800)),
            ("C3", "C", (1.249, 0.560, 0.180)), ("O3", "O", (2.420, -0.210, -0.090))],
    "SO4": [("S", "S", (0.0, 0.0, 0.0)), ("O1", "O", (0.854, 0.854, 0.854)), ("O2", "O", (-0.854, -0.854, 0.854)),
            ("O3", "O", (-0.854, 0.854, -0.854)), ("O4", "O", (0.854, -0.854, -0.854))],
    "EDO": [("C1", "C", (-0.650, 0.380, 0.0)), ("O1", "O", (-1.830, -0.400, 0.0)),
            ("C2", "C", (0.650, -0.380, 0.0)), ("O2", "O", (1.830, 0.400, 0.0))],
    "CL": [("CL", "CL", (0.0, 0.0, 0.0))],
    "NA": [("NA", "NA", (0.0, 0.0, 0.0))],
}
HETNAM = {"GOL": "GLYCEROL", "SO4": "SULFATE ION", "EDO": "1,2-ETHANEDIOL", "CL": "CHLORIDE ION",
          "NA": "SODIUM ION", "MSE": "SELENOMETHIONINE"}


@dataclass
class Atom:
    name: str
    xyz: np.ndarray
    elem: str
    b: float
    occ: float = 1.0
    alt: str = " "


@dataclass
class Res:
    name: str
    num: int                       # numbering of the source entry
    atoms: list = field(default_factory=list)
    label: str = ""                # residue number + insertion code as written in structure.pdb
    het: bool = False

    def get(self, name):
        for a in self.atoms:
            if a.name == name:
                return a
        return None


def read_entry(path) -> dict:
    """Protein chain A (heavy atoms, first conformer), crystal waters, SSBOND pairs, CRYST1."""
    residues: list[Res] = []
    waters = []
    ss = []
    cryst = None
    cur = None
    for line in open(path):
        rec = line[:6]
        if rec == "CRYST1":
            cryst = line.rstrip("\n")
        elif rec == "SSBOND":
            ss.append((int(line[17:21]), int(line[31:35])))
        elif rec == "ENDMDL":
            break
        elif rec in ("ATOM  ", "HETATM"):
            elem = line[76:78].strip().upper() or line[12:14].strip().upper()
            if elem in ("H", "D"):
                continue
            if line[16] not in (" ", "A"):
                continue
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            resn, chain, num = line[17:20].strip(), line[21], int(line[22:26])
            if rec == "HETATM":
                if resn == "HOH":
                    waters.append((xyz, float(line[60:66])))
                continue
            if chain != "A":
                continue
            if cur is None or cur.num != num:
                cur = Res(resn, num)
                residues.append(cur)
            cur.atoms.append(Atom(line[12:16].strip(), xyz, elem, float(line[60:66])))
    for i, r in enumerate(residues):
        want = set(HEAVY[r.name]) | ({"OXT"} if i == len(residues) - 1 else set())
        have = {a.name for a in r.atoms}
        if have != want:
            raise ValueError(f"source entry residue {r.name}{r.num} is not complete: {sorted(have ^ want)}")
    nums = [r.num for r in residues]
    if nums != list(range(nums[0], nums[0] + len(nums))):
        raise ValueError("source entry numbering is not contiguous")
    idx = {r.num: i for i, r in enumerate(residues)}
    pairs = sorted(tuple(sorted((idx[a], idx[b]))) for a, b in ss)
    # geometry cross-check of the deposited SSBOND records (truth by two routes)
    geo = []
    cys = [i for i, r in enumerate(residues) if r.name == "CYS"]
    for x in range(len(cys)):
        for y in range(x + 1, len(cys)):
            d = np.linalg.norm(residues[cys[x]].get("SG").xyz - residues[cys[y]].get("SG").xyz)
            if d < 2.5:
                geo.append((cys[x], cys[y]))
    if sorted(geo) != pairs:
        raise ValueError(f"SSBOND records {pairs} disagree with SG-SG geometry {geo}")
    return {"residues": residues, "waters": waters, "ss": pairs, "cryst": cryst}


def dihedral(p0, p1, p2, p3) -> float:
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1n, v), w), np.dot(v, w))))


def rotate_about(points, origin, axis, deg):
    k = axis / np.linalg.norm(axis)
    t = np.radians(deg)
    p = points - origin
    return origin + p * np.cos(t) + np.cross(k, p) * np.sin(t) + np.outer(p @ k, k) * (1 - np.cos(t))


def random_rotation(rng) -> np.ndarray:
    q = rng.normal(size=4)
    q /= np.linalg.norm(q)
    a, b, c, d = q
    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a-b*b+c*c-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a-b*b-c*c+d*d]])


def _aname(name: str, elem: str) -> str:
    if len(name) >= 4 or len(elem) == 2:
        return f"{name:<4}"
    return f" {name:<3}"


def atom_line(rec, serial, name, alt, resn, chain, num, icode, xyz, occ, b, elem) -> str:
    return (f"{rec:<6}{serial:>5} {_aname(name, elem)}{alt}{resn:>3} {chain}{num:>4}{icode or ' '}   "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{occ:6.2f}{b:6.2f}          {elem:>2}  ")


def split_label(label: str):
    if label and label[-1].isalpha():
        return int(label[:-1]), label[-1]
    return int(label), " "
