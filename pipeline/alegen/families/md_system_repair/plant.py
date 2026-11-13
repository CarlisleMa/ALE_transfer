"""Defect planting: complete structure (truth) -> defective PDB text + everything the verifier needs."""
from __future__ import annotations

import copy

import numpy as np

from .pdbkit import (CHARGE, GAMMA, HET_TEMPLATES, HETNAM, Atom, Res, atom_line, dihedral,
                     random_rotation, rotate_about)

TRUNCATABLE = {"LYS": ["CB", "CG"], "ARG": ["CB", "CG"], "GLU": ["CB", "CG"], "GLN": ["CB", "CG"],
               "ASP": ["CB"], "ASN": ["CB"], "LEU": ["CB", "CG"], "TYR": ["CB"], "PHE": ["CB"]}
ALT_OK = ("SER", "THR", "ILE", "LEU", "ASN", "ASP", "GLN", "GLU", "LYS", "MET", "ARG")


def _all_xyz(residues, skip=()):
    pts, owner = [], []
    for i, r in enumerate(residues):
        if i in skip:
            continue
        for a in r.atoms:
            pts.append(a.xyz)
            owner.append(i)
    return np.array(pts), np.array(owner)


def plant(entry: dict, rng: np.random.Generator, level: int) -> dict:
    full = entry["residues"]                  # complete truth
    n = len(full)
    ss = entry["ss"]
    cys = {i for p in ss for i in p}
    mets = [i for i, r in enumerate(full) if r.name == "MET"]
    plan = {"gap": None, "mse": None, "trunc": {}, "alt": {}, "offset": 0, "ins": None, "chain": "A"}

    # ---------------------------------------------------------------- choose what to plant
    used: set[int] = set()
    if level >= 3:
        length = int(rng.integers(4, 7))
        cands = []
        for s in range(6, n - 6 - length):
            seg = range(s - 1, s + length + 1)                    # segment plus both flanks
            if any(i in cys for i in seg) or any(full[i].name in ("MET", "PRO") for i in seg):
                continue
            d = np.linalg.norm(full[s - 1].get("C").xyz - full[s + length].get("N").xyz)
            if d > 6.0:
                cands.append(s)
        s = int(cands[int(rng.integers(0, len(cands)))])
        plan["gap"] = (s, s + length - 1)
        used |= set(range(s - 2, s + length + 2))
    if level >= 2:
        if mets:
            free = [m for m in mets if m not in used]
            plan["mse"] = int(free[int(rng.integers(0, len(free)))])
            used.add(plan["mse"])
        pool = [i for i in range(2, n - 2) if i not in used and full[i].name in TRUNCATABLE]
        k = 3 if n < 80 else 4
        for i in sorted(int(x) for x in rng.choice(pool, size=k, replace=False)):
            keep = TRUNCATABLE[full[i].name]
            keep = keep[:1] if (len(keep) == 1 or rng.random() < 0.5) else keep
            plan["trunc"][i] = ["N", "CA", "C", "O"] + keep
            used.add(i)
    if level >= 3:
        pts_all, owner = _all_xyz(full)
        pool = [i for i in range(2, n - 2) if i not in used and full[i].name in ALT_OK]
        order = [int(x) for x in rng.permutation(pool)]
        want = 3
        for i in order:
            if len(plan["alt"]) >= want:
                break
            r = full[i]
            side = [a for a in r.atoms if a.name not in ("N", "CA", "C", "O", "CB", "OXT")]
            ca, cb = r.get("CA").xyz, r.get("CB").xyz
            others = pts_all[owner != i]
            own_bb = np.array([r.get(x).xyz for x in ("N", "C", "O")])
            for deg in (120.0, -120.0) if rng.random() < 0.5 else (-120.0, 120.0):
                new = rotate_about(np.array([a.xyz for a in side]), cb, cb - ca, deg)
                d_o = np.min(np.linalg.norm(new[:, None, :] - others[None, :, :], axis=2))
                d_b = np.min(np.linalg.norm(new[:, None, :] - own_bb[None, :, :], axis=2))
                if d_o >= 3.0 and d_b >= 2.5:
                    plan["alt"][i] = {"decoy": {a.name: new[j] for j, a in enumerate(side)}, "deg": deg}
                    break
        if len(plan["alt"]) < 2:
            raise RuntimeError("could not place alternate conformers without clashes")
        alt_ids = sorted(plan["alt"])
        # at least one residue lists the minor (decoy) conformer first, at least one the major first
        modes = ["decoy_first", "truth_first"] + ["decoy_first" if rng.random() < 0.5 else "truth_first"
                                                  for _ in alt_ids[2:]]
        modes = [modes[int(j)] for j in rng.permutation(len(alt_ids))]
        for i, m in zip(alt_ids, modes):
            major = round(float(rng.uniform(0.58, 0.72)), 2)
            plan["alt"][i].update(mode=m, occ_truth=major, occ_decoy=round(1.0 - major, 2))
            used.add(i)
        plan["offset"] = int(rng.choice([17, 100, 201, 300]))
        plan["chain"] = str(rng.choice(list("BCEPX")))
        gap = plan["gap"]
        cands = []
        for c in sorted(cys):                      # put an insertion code on a disulfide cysteine
            for p in (c - 1, c - 2):
                trio = (p, p + 1, p + 2)
                if p < 2 or p + 2 > n - 3:
                    continue
                if any(gap[0] - 1 <= t <= gap[1] + 1 for t in trio):
                    continue
                cands.append(p)
        plan["ins"] = int(cands[int(rng.integers(0, len(cands)))])

    # ---------------------------------------------------------------- labels as written in the file
    def label(i: int) -> str:
        num = i + full[0].num + plan["offset"]
        p = plan["ins"]
        if p is None or i <= p:
            return str(num)
        if i == p + 1:
            return f"{p + full[0].num + plan['offset']}A"
        if i == p + 2:
            return f"{p + full[0].num + plan['offset']}B"
        return str(num - 2)

    gap = plan["gap"]
    in_gap = (lambda i: gap is not None and gap[0] <= i <= gap[1])
    kept = [i for i in range(n) if not in_gap(i)]

    # ---------------------------------------------------------------- hetero groups and waters
    chain = plan["chain"]
    het_lines_src = []      # (resname, [(name, elem, xyz)], b)
    xtal_sites = []         # positions of crystal waters / ions the SOP says to discard
    if level >= 2:
        prot = np.array([a.xyz for i in kept for a in full[i].atoms])
        waters = [w for w in entry["waters"] if np.min(np.linalg.norm(prot - w[0], axis=1)) >= 2.4]
        wpos = np.array([w[0] for w in waters])
        dmin = np.array([np.min(np.linalg.norm(prot - w, axis=1)) for w in wpos])
        centre = prot.mean(axis=0)
        placed = []          # all additive atom positions so far
        groups = ["GOL", "SO4", "SO4" if rng.random() < 0.5 else "EDO", "CL", "NA"]
        dead = set()
        for g in groups:
            tmpl = HET_TEMPLATES[g]
            ok = False
            for _ in range(400):
                wi = int(rng.integers(0, len(wpos)))
                if wi in dead or dmin[wi] < (2.9 if len(tmpl) == 1 else 3.2):
                    continue
                out = wpos[wi] - centre
                out /= np.linalg.norm(out)
                shift = 0.0 if len(tmpl) == 1 else float(rng.uniform(1.5, 2.5))
                rot = random_rotation(rng)
                xyz = np.array([t[2] for t in tmpl]) @ rot.T + wpos[wi] + shift * out
                if np.min(np.linalg.norm(prot[None, :, :] - xyz[:, None, :], axis=2)) < 3.1:
                    continue
                if placed and np.min(np.linalg.norm(np.array(placed)[None, :, :] - xyz[:, None, :], axis=2)) < 3.2:
                    continue
                ok = True
                break
            if not ok:
                raise RuntimeError(f"could not place additive {g}")
            placed += list(xyz)
            near = np.where(np.min(np.linalg.norm(wpos[:, None, :] - xyz[None, :, :], axis=2), axis=1) < 2.6)[0]
            dead |= set(int(x) for x in near) | {wi}
            het_lines_src.append((g, [(t[0], t[1], xyz[j]) for j, t in enumerate(tmpl)],
                                  float(rng.uniform(25, 45))))
            if g in ("CL", "NA"):
                xtal_sites.append(xyz[0])
        for wi, w in enumerate(waters):
            if wi in dead:
                continue
            het_lines_src.append(("HOH", [("O", "O", w[0])], w[1]))
            xtal_sites.append(w[0])

    # ---------------------------------------------------------------- write the defective file
    seqres_names = [("MSE" if i == plan["mse"] else r.name) for i, r in enumerate(full)]
    L = []
    L.append("HEADER    UNKNOWN FUNCTION                                    XXXX")
    L.append("TITLE     REFINED MODEL, CRYSTAL FORM I (IN-HOUSE DATA SET)")
    L.append("REMARK   2 RESOLUTION.    1.50 ANGSTROMS.")
    L.append("REMARK   3 PROGRAM     : REFMAC")
    if gap is not None:
        L += ["REMARK 465", "REMARK 465 MISSING RESIDUES",
              "REMARK 465 THE FOLLOWING RESIDUES WERE NOT LOCATED IN THE",
              "REMARK 465 EXPERIMENT. (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN",
              "REMARK 465 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)",
              "REMARK 465", "REMARK 465   M RES C SSSEQI"]
        for i in range(gap[0], gap[1] + 1):
            lab = label(i)
            num, ic = (lab[:-1], lab[-1]) if lab[-1].isalpha() else (lab, " ")
            L.append(f"REMARK 465     {full[i].name} {chain} {int(num):>5}{ic}")
    if plan["trunc"]:
        L += ["REMARK 470", "REMARK 470 MISSING ATOM",
              "REMARK 470 THE FOLLOWING RESIDUES HAVE MISSING ATOMS (M=MODEL NUMBER;",
              "REMARK 470 RES=RESIDUE NAME; C=CHAIN IDENTIFIER; SSEQ=SEQUENCE NUMBER;",
              "REMARK 470 I=INSERTION CODE):", "REMARK 470   M RES CSSEQI  ATOMS"]
        for i, keep in sorted(plan["trunc"].items()):
            lab = label(i)
            num, ic = (lab[:-1], lab[-1]) if lab[-1].isalpha() else (lab, " ")
            miss = [a.name for a in full[i].atoms if a.name not in keep]
            L.append(f"REMARK 470     {full[i].name} {chain}{int(num):>4}{ic}   " + " ".join(f"{m:<4}" for m in miss))
    for k in range(0, n, 13):
        L.append(f"SEQRES {k // 13 + 1:>3} {chain} {n:>4}  " + " ".join(seqres_names[k:k + 13]))
    if plan["mse"] is not None:
        lab = label(plan["mse"])
        num, ic = (lab[:-1], lab[-1]) if lab[-1].isalpha() else (lab, " ")
        L.append(f"MODRES XXXX MSE {chain} {int(num):>4}{ic} MET  SELENOMETHIONINE")
    hets = sorted({h[0] for h in het_lines_src if h[0] != "HOH"} | ({"MSE"} if plan["mse"] is not None else set()))
    for h in hets:
        L.append(f"HETNAM     {h} {HETNAM[h]}")
    if level == 1:
        for k, (a, b) in enumerate(ss):
            L.append(f"SSBOND {k + 1:>3} CYS {chain} {int(label(a)):>4}    CYS {chain} {int(label(b)):>4}")
    if entry["cryst"]:
        L.append(entry["cryst"][:70].rstrip())

    serial = 0
    last = None
    for i in kept:
        r = full[i]
        lab = label(i)
        num, ic = (int(lab[:-1]), lab[-1]) if lab[-1].isalpha() else (int(lab), " ")
        is_mse = (i == plan["mse"])
        resn = "MSE" if is_mse else r.name
        rec = "HETATM" if is_mse else "ATOM"
        keep = plan["trunc"].get(i)
        alt = plan["alt"].get(i)
        for a in r.atoms:
            if keep is not None and a.name not in keep:
                continue
            name, elem = a.name, a.elem
            if is_mse and name == "SD":
                name, elem = "SE", "SE"
            if alt is not None and a.name in alt["decoy"]:
                first = ("A", alt["decoy"][a.name], alt["occ_decoy"]) if alt["mode"] == "decoy_first" else \
                        ("A", a.xyz, alt["occ_truth"])
                second = ("B", a.xyz, alt["occ_truth"]) if alt["mode"] == "decoy_first" else \
                         ("B", alt["decoy"][a.name], alt["occ_decoy"])
                for al, xyz, occ in (first, second):
                    serial += 1
                    L.append(atom_line(rec, serial, name, al, resn, chain, num, ic, xyz, occ, a.b, elem))
            else:
                serial += 1
                L.append(atom_line(rec, serial, name, " ", resn, chain, num, ic, a.xyz, 1.0, a.b, elem))
        last = (resn, num, ic)
    serial += 1
    L.append(f"TER   {serial:>5}      {last[0]:>3} {chain}{last[1]:>4}{last[2]}")
    hnum = (int(label(n - 1).rstrip("AB")) // 100 + 2) * 100
    for g, atoms, b in het_lines_src:
        hnum += 1
        for name, elem, xyz in atoms:
            serial += 1
            L.append(atom_line("HETATM", serial, name, " ", g, chain, hnum, " ", xyz, 1.0, b, elem))
    L.append("END")
    pdb_text = "\n".join(x.rstrip() for x in L) + "\n"

    # ---------------------------------------------------------------- truth for the verifier
    new_index = {i: k for k, i in enumerate(kept)}
    residues = []
    for i in kept:
        r = full[i]
        residues.append({"name": r.name, "label": label(i),
                         "atoms": {a.name: [round(float(v), 3) for v in a.xyz] for a in r.atoms},
                         "resolved": (plan["trunc"].get(i) or [a.name for a in r.atoms])})
    fragments = [[0, len(kept) - 1]] if gap is None else [[0, new_index[gap[0] - 1]],
                                                          [new_index[gap[1] + 1], len(kept) - 1]]
    net = sum(CHARGE.get(full[i].name, 0) for i in kept)        # charged true termini cancel (+1, -1)
    watch = []
    for i, alt in sorted(plan["alt"].items()):
        r = full[i]
        g = GAMMA[r.name]
        chi = dihedral(r.get("N").xyz, r.get("CA").xyz, r.get("CB").xyz, r.get(g).xyz)
        watch.append({"idx": new_index[i], "gamma": g, "chi1": round(chi, 2),
                      "decoy_chi1": round(dihedral(r.get("N").xyz, r.get("CA").xyz, r.get("CB").xyz,
                                                   alt["decoy"][g]), 2),
                      "mode": alt["mode"]})
    truth = {
        "level": level,
        "residues": residues,
        "fragments": fragments,
        "disulfides": [[new_index[a], new_index[b]] for a, b in ss],
        "disulfide_labels": [[label(a), label(b)] for a, b in ss],
        "net_charge": int(net),
        "chi1_watch": watch,
        "xtal_sites": [[round(float(v), 3) for v in p] for p in xtal_sites],
        "plan": {"gap_labels": [label(i) for i in range(gap[0], gap[1] + 1)] if gap else [],
                 "gap_names": [full[i].name for i in range(gap[0], gap[1] + 1)] if gap else [],
                 "mse": label(plan["mse"]) if plan["mse"] is not None else None,
                 "truncated": {label(i): v for i, v in plan["trunc"].items()},
                 "altloc": {label(i): v["mode"] for i, v in plan["alt"].items()},
                 "offset": plan["offset"], "insertion_after": label(plan["ins"]) if plan["ins"] is not None else None,
                 "chain": chain, "additives": [h[0] for h in het_lines_src if h[0] != "HOH"],
                 "n_xtal_water": sum(1 for h in het_lines_src if h[0] == "HOH")},
    }
    # clean model of the truth (what a perfect repair would feed to the builder), used for the
    # planted-truth probe: resolved residues, true conformers, MET, no hetero groups.
    clean = copy.deepcopy([full[i] for i in kept])
    for k, i in enumerate(kept):
        clean[k].label = label(i)
    return {"pdb": pdb_text, "truth": truth, "clean": clean, "plan": plan}
