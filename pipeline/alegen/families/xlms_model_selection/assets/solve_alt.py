#!/usr/bin/env python3
"""Independent second solution (no Biopython): own fixed-column PDB reader, chains identified by
shared 6-mers, positions mapped through exact matching blocks (difflib) with a lysine check, and
the ranking by satisfied fraction with the mean inter-protein distance as tie-break. Used to show
that the verifier accepts a different but correct implementation."""
import csv
import difflib
import json
import sys
from pathlib import Path

import numpy as np

inp = Path(sys.argv[1] if len(sys.argv) > 2 else "/task/input")
out = Path(sys.argv[2] if len(sys.argv) > 2 else "/task/output")
out.mkdir(parents=True, exist_ok=True)
THREE = dict(ALA="A", ARG="R", ASN="N", ASP="D", CYS="C", GLN="Q", GLU="E", GLY="G", HIS="H", ILE="I",
             LEU="L", LYS="K", MET="M", PHE="F", PRO="P", SER="S", THR="T", TRP="W", TYR="Y", VAL="V")

# full-length sequences
full, acc = {}, None
for line in (inp / "sequences.fasta").read_text().splitlines():
    if line.startswith(">"):
        acc = line[1:].split()[0]
        full[acc] = ""
    elif acc:
        full[acc] += line.strip()

# distance criterion: the number that follows "at most" in the hand-over notes
cutoff = None
for line in (inp / "xlms_notes.md").read_text().splitlines():
    if "CA-CA" in line and "at most" in line:
        cutoff = float(line.split("at most")[1].split()[0])
assert cutoff, "no distance criterion found"

links = list(csv.DictReader(open(inp / "crosslinks.csv")))


def read_pdb(path):
    chains = {}
    for line in open(path):
        if line.startswith("ATOM") and line[12:16].strip() == "CA" and line[16] in " A":
            chains.setdefault(line[21], []).append(
                (int(line[22:26]), THREE.get(line[17:20], "X"),
                 (float(line[30:38]), float(line[38:46]), float(line[46:54]))))
    return chains


def kmers(s, k=6):
    return {s[i:i + k] for i in range(len(s) - k + 1)}


full_k = {a: kmers(s) for a, s in full.items()}


def place(chain):
    """Returns (accession, {full position: (residue number, CA)})."""
    seq = "".join(r[1] for r in chain)
    a = max(full, key=lambda x: len(kmers(seq) & full_k[x]))
    blocks = difflib.SequenceMatcher(None, seq, full[a], autojunk=False).get_matching_blocks()
    mp = {}
    for i, j, n in blocks:
        for k in range(n):
            mp[j + k + 1] = (chain[i + k][0], np.array(chain[i + k][2]))
    assert len(mp) == len(seq), "model chain is not a set of exact pieces of the full-length sequence"
    return a, mp


models = sorted((inp / "models").glob("model_*.pdb"))
names = [p.stem for p in models]
rows, score, first_layout = [], {}, None
for p in models:
    lay = {}
    for cid, chain in read_pdb(p).items():
        a, mp = place(chain)
        assert a not in lay, "two chains of the same protein"
        lay[a] = (cid, mp)
    if first_layout is None:
        first_layout = lay
    c = dict(inter_satisfied=0, inter_total=0, intra_satisfied=0, intra_total=0)
    inter_d = []
    for x in links:
        pa, pb = x["protein_a"], x["protein_b"]
        ra, rb = lay[pa][1].get(int(x["position_a"])), lay[pb][1].get(int(x["position_b"]))
        if ra is None or rb is None:
            continue
        assert full[pa][int(x["position_a"]) - 1] == "K" and full[pb][int(x["position_b"]) - 1] == "K"
        d = float(np.sqrt(((ra[1] - rb[1]) ** 2).sum()))
        cls = "inter" if pa != pb else "intra"
        c[cls + "_total"] += 1
        c[cls + "_satisfied"] += d <= cutoff
        if cls == "inter":
            inter_d.append(d)
    rows.append((p.stem, c))
    score[p.stem] = (-c["inter_satisfied"] / max(c["inter_total"], 1), float(np.mean(inter_d)) if inter_d else 1e9)

with open(out / "restraints.tsv", "w") as fh:
    fh.write("model\tinter_satisfied\tinter_total\tintra_satisfied\tintra_total\n")
    for n, c in rows:
        fh.write(f"{n}\t{c['inter_satisfied']}\t{c['inter_total']}\t{c['intra_satisfied']}\t{c['intra_total']}\n")
ranking = sorted(names, key=lambda n: score[n])
json.dump({"best_model": ranking[0], "ranking": ranking}, open(out / "selection.json", "w"), indent=1)
seen = []
for x in links:
    for k in ((x["protein_a"], int(x["position_a"])), (x["protein_b"], int(x["position_b"]))):
        if k not in seen:
            seen.append(k)
with open(out / "mapping.tsv", "w") as fh:
    fh.write("accession\tposition\tmodel_chain\tmodel_resnum\n")
    for a, pos in seen:
        cid, mp = first_layout[a]
        hit = mp.get(pos)
        fh.write(f"{a}\t{pos}\t{cid if hit else 'NA'}\t{hit[0] if hit else 'NA'}\n")
