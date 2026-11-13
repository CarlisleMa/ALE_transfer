#!/usr/bin/env python3
"""Reference analysis: map cross-linked positions onto every model (Biopython alignment of each
model chain against the full-length sequences), measure CA-CA distances, count satisfied
inter-/intra-protein cross-links per model and rank the models by inter-protein support.

The non-default flags replay plausible-but-wrong analyses for the validator's attack suite; the
oracle (solve.sh) uses the defaults.
"""
import argparse
import csv
import json
import re
from pathlib import Path

import numpy as np
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

ap = argparse.ArgumentParser()
ap.add_argument("--input", default="/task/input")
ap.add_argument("--output", default="/task/output")
ap.add_argument("--numbering", choices=["align", "raw", "offset"], default="align",
                help="align: sequence alignment per chain; raw: position = model residue number; "
                     "offset: one constant shift per protein")
ap.add_argument("--chains", choices=["sequence", "first"], default="sequence",
                help="sequence: identify the chains of every model by sequence; first: reuse model_01's letters")
ap.add_argument("--select", choices=["crosslinks", "confidence"], default="crosslinks")
ap.add_argument("--prune", action="store_true", help="discard cross-links that no model satisfies")
ap.add_argument("--pooled", action="store_true", help="do not separate inter- and intra-protein cross-links")
a = ap.parse_args()
inp, out = Path(a.input), Path(a.output)
out.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------------- inputs
full = {rec.id: str(rec.seq) for rec in SeqIO.parse(str(inp / "sequences.fasta"), "fasta")}
with open(inp / "crosslinks.csv") as fh:
    links = [(r["xl_id"], r["protein_a"], int(r["position_a"]), r["protein_b"], int(r["position_b"]))
             for r in csv.DictReader(fh)]
notes = (inp / "xlms_notes.md").read_text()
m = re.search(r"CA-CA distance[^.]*?at most\s+([0-9]+(?:\.[0-9]+)?)\s*A", notes)
if not m:
    raise SystemExit("distance criterion not found in xlms_notes.md")
cutoff = float(m.group(1))
model_files = sorted((inp / "models").glob("model_*.pdb"))
names = [p.stem for p in model_files]

parser = PDBParser(QUIET=True)


def read_model(path):
    """chain id -> (sequence, residue numbers, {residue number: CA coordinate})"""
    chains = {}
    for ch in parser.get_structure(path.stem, str(path))[0]:
        res = [r for r in ch if r.id[0] == " " and "CA" in r]
        chains[ch.id] = ("".join(seq1(r.get_resname()) for r in res), [r.id[1] for r in res],
                         {r.id[1]: r["CA"].coord.astype(float) for r in res})
    return chains


aligner = PairwiseAligner()
aligner.mode = "local"
aligner.match_score = 1.0
aligner.mismatch_score = -6.0
aligner.open_gap_score = -4.0
aligner.extend_gap_score = 0.0


def which_protein(chain_seq):
    return max(full, key=lambda acc: aligner.score(full[acc], chain_seq))


def position_map(acc, chain):
    """full-length position (1-based) -> residue number in this chain, for the chosen numbering mode."""
    seq, nums, _ = chain
    if a.numbering == "raw":
        return {n: n for n in nums}
    if a.numbering == "offset":
        start = full[acc].find(seq[:10])
        if start < 0:
            return {}
        shift = (start + 1) - nums[0]
        return {n + shift: n for n in nums}
    aln = aligner.align(full[acc], seq)[0]
    mp = {}
    for (t0, t1), (q0, q1) in zip(*aln.aligned):
        for k in range(t1 - t0):
            if full[acc][t0 + k] != seq[q0 + k]:
                raise SystemExit(f"model sequence does not match {acc} at position {t0 + k + 1}")
            mp[t0 + k + 1] = nums[q0 + k]
    if len(mp) != len(seq):
        raise SystemExit(f"could not place every model residue on {acc}")
    return mp


models = {n: read_model(p) for n, p in zip(names, model_files)}
first = names[0]
layout = {}      # model -> acc -> (chain id or None, {position: residue number})
ref = {}
for cid, chain in models[first].items():
    acc = which_protein(chain[0])
    ref[acc] = (cid, position_map(acc, chain))
for n in names:
    if a.chains == "first":
        layout[n] = ref
        continue
    lay = {}
    for cid, chain in models[n].items():
        acc = which_protein(chain[0])
        lay[acc] = (cid, position_map(acc, chain))
    layout[n] = lay


def ca(n, acc, pos):
    cid, mp = layout[n].get(acc, (None, {}))
    if cid not in models[n] or pos not in mp:
        return None
    return models[n][cid][2].get(mp[pos])


# ---------------------------------------------------------------------------------- distances
D = np.full((len(links), len(names)), np.nan)
for i, (_, pa, xa, pb, xb) in enumerate(links):
    for j, n in enumerate(names):
        ca1, ca2 = ca(n, pa, xa), ca(n, pb, xb)
        if ca1 is not None and ca2 is not None:
            D[i, j] = np.linalg.norm(ca1 - ca2)
evaluable = ~np.isnan(D)
sat = evaluable & (np.nan_to_num(D, nan=1e9) <= cutoff)
if a.prune:
    keep = sat.any(axis=1)
    evaluable, sat = evaluable & keep[:, None], sat & keep[:, None]
inter = np.array([pa != pb for _, pa, _, pb, _ in links])

# ---------------------------------------------------------------------------------- outputs
table = []
for j, n in enumerate(names):
    row = {"model": n,
           "inter_satisfied": int(sat[inter, j].sum()), "inter_total": int(evaluable[inter, j].sum()),
           "intra_satisfied": int(sat[~inter, j].sum()), "intra_total": int(evaluable[~inter, j].sum())}
    excess = np.where(evaluable[inter, j], np.clip(np.nan_to_num(D[inter, j], nan=0.0) - cutoff, 0, None), 0.0)
    row["inter_violation_sum_A"] = round(float(excess.sum()), 2)
    if a.pooled:
        row.update(inter_satisfied=int(sat[:, j].sum()), inter_total=int(evaluable[:, j].sum()),
                   intra_satisfied=0, intra_total=0)
    table.append(row)
cols = ["model", "inter_satisfied", "inter_total", "intra_satisfied", "intra_total", "inter_violation_sum_A"]
with open(out / "restraints.tsv", "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for row in table:
        fh.write("\t".join(str(row[c]) for c in cols) + "\n")

if a.select == "confidence":
    rk = inp / "ranking.json"
    if rk.exists():
        conf = json.loads(rk.read_text())["confidence"]
        ranking = sorted(names, key=lambda n: -conf[n])
    else:
        ranking = list(names)
else:
    by = {r["model"]: r for r in table}
    ranking = sorted(names, key=lambda n: (-by[n]["inter_satisfied"], by[n]["inter_violation_sum_A"], n))
(out / "selection.json").write_text(json.dumps({"best_model": ranking[0], "ranking": ranking}, indent=1) + "\n")

positions = sorted({(p, x) for _, pa, xa, pb, xb in links for p, x in ((pa, xa), (pb, xb))})
with open(out / "mapping.tsv", "w") as fh:
    fh.write("accession\tposition\tmodel_chain\tmodel_resnum\n")
    for acc, pos in positions:
        cid, mp = ref.get(acc, (None, {}))
        if cid is not None and pos in mp and mp[pos] in models[first][cid][2]:
            fh.write(f"{acc}\t{pos}\t{cid}\t{mp[pos]}\n")
        else:
            fh.write(f"{acc}\t{pos}\tNA\tNA\n")
