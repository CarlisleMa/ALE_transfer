#!/usr/bin/env python3
"""Reference analysis for level 4.

Steps: read the hand-over notes (cutoff, blocking radius, tether length); parse every model with
Biopython; identify the protein of every chain and map full-length positions onto residue numbers
by local alignment; build the lattice of the notes (k-d tree for the blocked points, sparse graph
over the open points with one extra node per lysine CA that is joined to its tether points);
single-source Dijkstra (scipy) from every lysine that starts a needed pair; a cross-link is
satisfied if any allowed assignment of its two positions to chains is within range; counts, per-link
status, ranking by the total number of satisfied cross-links, mapping table.
"""
import argparse
import csv
import json
import os
import re
from pathlib import Path

import numpy as np
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra
from scipy.spatial import cKDTree

ap = argparse.ArgumentParser()
ap.add_argument("--input", default="/task/input")
ap.add_argument("--output", default="/task/output")
a = ap.parse_args()
inp, out = Path(a.input), Path(a.output)
out.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------------- inputs
full = {rec.id: str(rec.seq) for rec in SeqIO.parse(str(inp / "sequences.fasta"), "fasta")}
with open(inp / "crosslinks.csv") as fh:
    links = [(r["xl_id"], r["protein_a"], int(r["position_a"]), r["protein_b"], int(r["position_b"]))
             for r in csv.DictReader(fh)]
notes = (inp / "xlms_notes.md").read_text()
num = r"([0-9]+(?:\.[0-9]+)?)\s*A\b"
m_cut = re.search(r"within range if its surface distance is at most\s+" + num, notes)
m_rad = re.search(r"non-hydrogen atom of the model is less than\s+" + num, notes)
m_tet = re.search(r"segments may be at most\s+" + num, notes)
if not (m_cut and m_rad and m_tet):
    raise SystemExit("distance criterion not found in xlms_notes.md")
CUTOFF, RADIUS, TETHER = float(m_cut.group(1)), float(m_rad.group(1)), float(m_tet.group(1))
# distances are followed a little beyond the cutoff for the diagnostic table
LIMIT = CUTOFF + float(os.environ.get("XL_LIMIT_EXTRA", "3.5"))
model_files = sorted((inp / "models").glob("model_*.pdb"))
names = [p.stem for p in model_files]

parser = PDBParser(QUIET=True)
aligner = PairwiseAligner()
aligner.mode = "local"
aligner.match_score = 1.0
aligner.mismatch_score = -6.0
aligner.open_gap_score = -4.0
aligner.extend_gap_score = 0.0


def read_model(path):
    """-> ([(chain id, accession, {position: residue number}, {residue number: CA})], all heavy atoms)"""
    chains, atoms = [], []
    for ch in parser.get_structure(path.stem, str(path))[0]:
        res = [r for r in ch if r.id[0] == " " and "CA" in r]
        seq = "".join(seq1(r.get_resname()) for r in res)
        nums = [r.id[1] for r in res]
        acc = max(full, key=lambda x: aligner.score(full[x], seq))
        aln = aligner.align(full[acc], seq)[0]
        mp = {}
        for (t0, t1), (q0, q1) in zip(*aln.aligned):
            for k in range(t1 - t0):
                if full[acc][t0 + k] != seq[q0 + k]:
                    raise SystemExit(f"model sequence does not match {acc} at position {t0 + k + 1}")
                mp[t0 + k + 1] = nums[q0 + k]
        if len(mp) != len(seq):
            raise SystemExit(f"could not place every model residue on {acc}")
        chains.append((ch.id, acc, mp, {r.id[1]: r["CA"].coord.astype(np.float64) for r in res}))
        atoms += [at.coord for r in ch for at in r if at.element != "H"]
    return chains, np.array(atoms, dtype=np.float64)


OFFS = np.array([(dx, dy, dz) for dx in (-1, 0, 1) for dy in (-1, 0, 1) for dz in (-1, 0, 1) if (dx, dy, dz) != (0, 0, 0)])
STEP = np.sqrt((OFFS ** 2).sum(1).astype(float))


def surface_distances(atoms, anchors, pairs):
    """Surface distance for the needed anchor pairs. anchors: {key: CA xyz}; pairs: set of (key, key)."""
    pad = int(np.ceil(RADIUS)) + 2
    lo = np.floor(atoms.min(0)).astype(int) - pad
    hi = np.ceil(atoms.max(0)).astype(int) + pad
    shape = hi - lo + 1
    grid = np.stack(np.meshgrid(*[np.arange(lo[i], hi[i] + 1) for i in range(3)], indexing="ij"), axis=-1)
    pts = grid.reshape(-1, 3)
    nearest, _ = cKDTree(atoms).query(pts.astype(float), k=1, workers=2)
    is_open = (nearest >= RADIUS).reshape(shape)
    node = np.full(tuple(shape), -1, dtype=np.int64)
    n_open = int(is_open.sum())
    node[is_open] = np.arange(n_open)
    open_pts = pts[is_open.ravel()]
    # neighbour table of the open points (-1 where the neighbour is blocked or outside the box)
    padded = np.full(tuple(shape + 2), -1, dtype=np.int64)
    padded[1:-1, 1:-1, 1:-1] = node
    ijk = open_pts - lo + 1
    nbr = np.stack([padded[ijk[:, 0] + o[0], ijk[:, 1] + o[1], ijk[:, 2] + o[2]] for o in OFFS], axis=1)
    valid = nbr >= 0
    keys = sorted(anchors)
    kidx = {k: n_open + i for i, k in enumerate(keys)}
    tree = cKDTree(open_pts.astype(float))
    teth = {}
    for k in keys:
        idx = np.array(tree.query_ball_point(anchors[k], TETHER), dtype=np.int64)
        length = np.linalg.norm(open_pts[idx] - anchors[k], axis=1) if len(idx) else np.zeros(0)
        keep = length <= TETHER
        teth[k] = (idx[keep], length[keep])
    # directed graph: lattice steps in both directions, plus anchor -> its tether points (one way only,
    # so that no path can pass through an anchor)
    counts = np.concatenate([valid.sum(1), [len(teth[k][0]) for k in keys]])
    indptr = np.concatenate([[0], np.cumsum(counts)])
    indices = np.concatenate([nbr[valid]] + [teth[k][0] for k in keys])
    data = np.concatenate([np.broadcast_to(STEP, nbr.shape)[valid]] + [teth[k][1] for k in keys])
    n = n_open + len(keys)
    graph = csr_matrix((data, indices, indptr), shape=(n, n))
    result = {}
    left = set(pairs)
    while left:
        deg = {}
        for k1, k2 in left:
            deg[k1] = deg.get(k1, 0) + 1
            deg[k2] = deg.get(k2, 0) + 1
        src = max(sorted(deg), key=lambda k: deg[k])
        d = dijkstra(graph, directed=True, indices=kidx[src], limit=LIMIT)
        for pr in [pr for pr in left if src in pr]:
            other = pr[0] if pr[1] == src else pr[1]
            idx, length = teth[other]
            v = float((d[idx] + length).min()) if len(idx) else np.inf
            result[pr] = v if v <= LIMIT else np.inf
            left.discard(pr)
    return result


# ---------------------------------------------------------------------------------- evaluation
dist = {}                      # (model, xl_id) -> smallest surface distance over the assignments, None = not evaluable
layout01 = None
for path in model_files:
    chains, atoms = read_model(path)
    if layout01 is None:
        layout01 = chains
    by_acc = {}
    for cid, acc, mp, ca in chains:
        by_acc.setdefault(acc, []).append((cid, mp, ca))
    todo, anchors = {}, {}
    for xl, pa, xa, pb, xb in links:
        cand = []
        for cid_a, mp_a, ca_a in by_acc.get(pa, []):
            for cid_b, mp_b, ca_b in by_acc.get(pb, []):
                if xa not in mp_a or xb not in mp_b:
                    continue
                ka, kb = (cid_a, mp_a[xa]), (cid_b, mp_b[xb])
                if ka == kb:                       # the two ends of a cross-link are two different residues
                    continue
                anchors[ka], anchors[kb] = ca_a[mp_a[xa]], ca_b[mp_b[xb]]
                cand.append(tuple(sorted((ka, kb))))
        todo[xl] = cand
    d = surface_distances(atoms, anchors, {c for cand in todo.values() for c in cand})
    for xl, cand in todo.items():
        dist[(path.stem, xl)] = min(d[c] for c in cand) if cand else None

# ---------------------------------------------------------------------------------- outputs
status = {k: ("not_evaluable" if v is None else ("satisfied" if v <= CUTOFF else "violated")) for k, v in dist.items()}
table = {}
for n in names:
    row = {"inter_satisfied": 0, "inter_total": 0, "intra_satisfied": 0, "intra_total": 0}
    for xl, pa, xa, pb, xb in links:
        st = status[(n, xl)]
        if st == "not_evaluable":
            continue
        cls = "inter" if pa != pb else "intra"
        row[cls + "_total"] += 1
        row[cls + "_satisfied"] += int(st == "satisfied")
    table[n] = row
cols = ["inter_satisfied", "inter_total", "intra_satisfied", "intra_total"]
with open(out / "restraints.tsv", "w") as fh:
    fh.write("model\t" + "\t".join(cols) + "\n")
    for n in names:
        fh.write(n + "\t" + "\t".join(str(table[n][c]) for c in cols) + "\n")
with open(out / "link_status.tsv", "w") as fh:
    fh.write("model\txl_id\tstatus\n")
    for n in names:
        for xl, *_ in links:
            fh.write(f"{n}\t{xl}\t{status[(n, xl)]}\n")
with open(out / "distances.tsv", "w") as fh:       # diagnostic, not a required output
    fh.write("model\txl_id\tdistance_A\n")
    for n in names:
        for xl, *_ in links:
            v = dist[(n, xl)]
            fh.write(f"{n}\t{xl}\t{'NA' if v is None else ('inf' if not np.isfinite(v) else f'{v:.3f}')}\n")
ranking = sorted(names, key=lambda n: (-(table[n]["inter_satisfied"] + table[n]["intra_satisfied"]), n))
(out / "selection.json").write_text(json.dumps({"best_model": ranking[0], "ranking": ranking}, indent=1) + "\n")

positions = sorted({(p, x) for _, pa, xa, pb, xb in links for p, x in ((pa, xa), (pb, xb))})
with open(out / "mapping.tsv", "w") as fh:
    fh.write("accession\tposition\tmodel_chain\tmodel_resnum\n")
    for acc, pos in positions:
        hit = [(cid, mp[pos]) for cid, acc2, mp, ca in layout01 if acc2 == acc and pos in mp and mp[pos] in ca]
        if hit:
            fh.write(f"{acc}\t{pos}\t{','.join(sorted(h[0] for h in hit))}\t{hit[0][1]}\n")
        else:
            fh.write(f"{acc}\t{pos}\tNA\tNA\n")
