#!/usr/bin/env python3
"""Independent second solution for level 4 (no Biopython, no scipy): own fixed-column PDB reader,
chains identified by shared 6-mers, positions mapped through exact matching blocks (difflib),
surface distances by a bucketed lattice search in numpy (sasd_np.py).

The non-default flags replay plausible-but-wrong analyses for the validator's attack suite.
"""
import argparse
import csv
import difflib
import json
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sasd_np  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--input", default="/task/input")
ap.add_argument("--output", default="/task/output")
ap.add_argument("--distance", choices=["surface", "straight"], default="surface",
                help="straight: CA-CA straight line with the cutoff of the notes")
ap.add_argument("--copies", choices=["any", "within", "last"], default="any",
                help="any: every assignment of the two positions to chains; within: two positions of one protein "
                     "are looked up in one chain only; last: one chain per protein (the level-3 logic)")
ap.add_argument("--numbering", choices=["blocks", "raw", "offset"], default="blocks")
ap.add_argument("--chains", choices=["sequence", "first"], default="sequence")
ap.add_argument("--select", choices=["total", "inter", "confidence"], default="total")
ap.add_argument("--prune", action="store_true", help="discard cross-links that no model satisfies")
ap.add_argument("--pooled", action="store_true", help="do not separate inter- and intra-protein cross-links")
a = ap.parse_args()
inp, out = Path(a.input), Path(a.output)
out.mkdir(parents=True, exist_ok=True)
THREE = dict(ALA="A", ARG="R", ASN="N", ASP="D", CYS="C", GLN="Q", GLU="E", GLY="G", HIS="H", ILE="I",
             LEU="L", LYS="K", MET="M", PHE="F", PRO="P", SER="S", THR="T", TRP="W", TYR="Y", VAL="V")

full, acc = {}, None
for line in (inp / "sequences.fasta").read_text().splitlines():
    if line.startswith(">"):
        acc = line[1:].split()[0]
        full[acc] = ""
    elif acc:
        full[acc] += line.strip()


def number_after(text, key):
    return float(text.split(key)[1].split()[0])


cutoff = radius = tether = None
for line in (inp / "xlms_notes.md").read_text().splitlines():
    if "within range if" in line and "is at most" in line:
        cutoff = number_after(line, "is at most")
    if "Blocked points" in line and "is less than" in line:
        radius = number_after(line, "is less than")
    if "Tethers" in line and "may be at most" in line:
        tether = number_after(line, "may be at most")
assert cutoff and radius and tether, "distance criterion not found in xlms_notes.md"
LIMIT = cutoff + float(os.environ.get("XL_LIMIT_EXTRA", "3.5"))   # followed a little beyond the cutoff for distances.tsv

links = [(r["xl_id"], r["protein_a"], int(r["position_a"]), r["protein_b"], int(r["position_b"]))
         for r in csv.DictReader(open(inp / "crosslinks.csv"))]


def read_pdb(path):
    """-> ({chain: [(resnum, aa, CA xyz)]} in file order, all atom coordinates)"""
    chains, atoms = {}, []
    for line in open(path):
        if not line.startswith("ATOM"):
            continue
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        atoms.append(xyz)
        if line[12:16].strip() == "CA":
            chains.setdefault(line[21], []).append((int(line[22:26]), THREE.get(line[17:20], "X"), xyz))
    return chains, np.array(atoms)


def kmers(s, k=6):
    return {s[i:i + k] for i in range(len(s) - k + 1)}


full_k = {x: kmers(s) for x, s in full.items()}


def place(chain):
    """-> (accession, {full-length position: residue number})"""
    seq = "".join(r[1] for r in chain)
    x = max(full, key=lambda y: len(kmers(seq) & full_k[y]))
    nums = [r[0] for r in chain]
    if a.numbering == "raw":
        return x, {n: n for n in nums}
    if a.numbering == "offset":
        start = full[x].find(seq[:10])
        return x, ({} if start < 0 else {n + (start + 1) - nums[0]: n for n in nums})
    mp = {}
    for i, j, n in difflib.SequenceMatcher(None, seq, full[x], autojunk=False).get_matching_blocks():
        for k in range(n):
            mp[j + k + 1] = nums[i + k]
    assert len(mp) == len(seq), "model chain is not a set of exact pieces of the full-length sequence"
    return x, mp


files = sorted((inp / "models").glob("model_*.pdb"))
names = [p.stem for p in files]
first_layout = None
dist = {}        # (model, xl_id) -> smallest distance over the assignments that were tried, None = not evaluable
for p in files:
    chains, atoms = read_pdb(p)
    lay = {}     # accession -> [(chain id, {position: residue number})] in file order
    for cid, chain in chains.items():
        x, mp = place(chain)
        lay.setdefault(x, []).append((cid, mp))
    if first_layout is None:
        first_layout = lay
    if a.chains == "first":
        lay = first_layout
    ca = {(cid, r[0]): np.array(r[2]) for cid, chain in chains.items() for r in chain}

    def atom(x, pos, which):
        cid, mp = which
        return (cid, mp[pos]) if pos in mp and (cid, mp[pos]) in ca else None

    todo = {}    # xl_id -> list of (anchor key, anchor key)
    for xl, pa, xa, pb, xb in links:
        ca_list, cb_list = lay.get(pa, []), lay.get(pb, [])
        if a.copies == "last":
            ca_list, cb_list = ca_list[-1:], cb_list[-1:]
        cand = []
        for wa in ca_list:
            for wb in cb_list:
                same_chain = wa[0] == wb[0]
                if a.copies == "within" and pa == pb and not same_chain:
                    continue
                if a.copies == "any" and pa == pb and xa == xb and same_chain:
                    continue                                   # a residue is not cross-linked to itself
                ka, kb = atom(pa, xa, wa), atom(pb, xb, wb)
                if ka is not None and kb is not None:
                    cand.append((ka, kb))
        todo[xl] = cand
    pairs = {tuple(sorted(c)) for cand in todo.values() for c in cand if c[0] != c[1]}
    d = {}
    if a.distance == "straight":
        for k1, k2 in pairs:
            d[(k1, k2)] = float(np.sqrt(((ca[k1] - ca[k2]) ** 2).sum()))
    elif pairs:
        lat = sasd_np.Lattice(atoms, radius)
        tp = {}
        left = set(pairs)
        while left:                                            # sources: a greedy vertex cover of the needed pairs
            deg = {}
            for k1, k2 in left:
                deg[k1] = deg.get(k1, 0) + 1
                deg[k2] = deg.get(k2, 0) + 1
            src = max(sorted(deg), key=lambda k: deg[k])
            for k in [src] + [q for pr in left if src in pr for q in pr if q != src]:
                if k not in tp:
                    tp[k] = lat.tether_points(ca[k], tether)
            field = lat.field(tp[src][0], tp[src][1], LIMIT) if len(tp[src][0]) else None
            for pr in [pr for pr in left if src in pr]:
                other = pr[0] if pr[1] == src else pr[1]
                v = np.inf
                if field is not None and len(tp[other][0]):
                    v = float((field[tp[other][0]] + tp[other][1]).min())
                d[pr] = v if v <= LIMIT else np.inf
                left.discard(pr)
    for xl, cand in todo.items():
        if not cand:
            dist[(p.stem, xl)] = None
        else:
            dist[(p.stem, xl)] = min(0.0 if c[0] == c[1] else d[tuple(sorted(c))] for c in cand)

sat = {k: (v is not None and v <= cutoff) for k, v in dist.items()}
if a.prune:
    for xl, *_ in links:
        if not any(sat[(n, xl)] for n in names):
            for n in names:
                dist[(n, xl)] = None
rows = []
for n in names:
    c = dict(inter_satisfied=0, inter_total=0, intra_satisfied=0, intra_total=0)
    for xl, pa, xa, pb, xb in links:
        if dist[(n, xl)] is None:
            continue
        cls = "inter" if pa != pb else "intra"
        c[cls + "_total"] += 1
        c[cls + "_satisfied"] += int(sat[(n, xl)])
    c["total_satisfied"] = c["inter_satisfied"] + c["intra_satisfied"]
    if a.pooled:
        c.update(inter_satisfied=c["total_satisfied"], inter_total=c["inter_total"] + c["intra_total"],
                 intra_satisfied=0, intra_total=0)
    rows.append((n, c))
cols = ["inter_satisfied", "inter_total", "intra_satisfied", "intra_total", "total_satisfied"]
with open(out / "restraints.tsv", "w") as fh:
    fh.write("model\t" + "\t".join(cols) + "\n")
    for n, c in rows:
        fh.write(n + "\t" + "\t".join(str(c[k]) for k in cols) + "\n")
with open(out / "link_status.tsv", "w") as fh:
    fh.write("model\txl_id\tstatus\n")
    for n in names:
        for xl, *_ in links:
            st = "not_evaluable" if dist[(n, xl)] is None else ("satisfied" if sat[(n, xl)] else "violated")
            fh.write(f"{n}\t{xl}\t{st}\n")
with open(out / "distances.tsv", "w") as fh:                   # diagnostic, not a required output
    fh.write("model\txl_id\tdistance_A\n")
    for n in names:
        for xl, *_ in links:
            v = dist[(n, xl)]
            fh.write(f"{n}\t{xl}\t{'NA' if v is None else ('inf' if not np.isfinite(v) else f'{v:.3f}')}\n")

by = dict(rows)
if a.select == "confidence" and (inp / "ranking.json").exists():
    conf = json.loads((inp / "ranking.json").read_text())["confidence"]
    ranking = sorted(names, key=lambda n: -conf[n])
elif a.select == "inter":
    ranking = sorted(names, key=lambda n: (-by[n]["inter_satisfied"], n))
else:
    def excess(n):                                             # tie-break: how far the violated links are out
        v = [dist[(n, xl)] for xl, *_ in links if dist[(n, xl)] is not None and not sat[(n, xl)]]
        return sum(min(x, LIMIT) - cutoff for x in v)
    ranking = sorted(names, key=lambda n: (-by[n]["total_satisfied"], excess(n), n))
json.dump({"best_model": ranking[0], "ranking": ranking}, open(out / "selection.json", "w"), indent=1)

seen = []
for xl, pa, xa, pb, xb in links:
    for k in ((pa, xa), (pb, xb)):
        if k not in seen:
            seen.append(k)
model01 = read_pdb(files[0])[0]
with open(out / "mapping.tsv", "w") as fh:
    fh.write("accession\tposition\tmodel_chain\tmodel_resnum\n")
    for x, pos in seen:
        which = first_layout.get(x, [])
        if a.copies == "last":
            which = which[-1:]
        hit = [(cid, mp[pos]) for cid, mp in which if pos in mp and any(r[0] == mp[pos] for r in model01[cid])]
        if hit:
            fh.write(f"{x}\t{pos}\t{','.join(sorted(h[0] for h in hit))}\t{hit[0][1]}\n")
        else:
            fh.write(f"{x}\t{pos}\tNA\tNA\n")
