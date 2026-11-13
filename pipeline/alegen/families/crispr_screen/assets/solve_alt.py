#!/usr/bin/env python3
"""Independent second solution (no cutadapt, no MAGeCK): anchor on the scaffold in either
orientation, 1-mismatch-tolerant spacer lookup, control-guide size factors, and an empirical null
built from random sets of non-targeting guides. Used to show the verifier accepts a different but
correct method."""
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd

inp, out = Path("/task/input"), Path("/task/output")
out.mkdir(parents=True, exist_ok=True)
ANCHOR = "GTTTTAGAGCTA"
COMP = str.maketrans("ACGTN", "TGCAN")

lib = pd.read_csv(inp / "library.tsv", sep="\t")
samples = pd.read_csv(inp / "samples.tsv", sep="\t")
sids = list(samples["sample_id"])
exact = dict(zip(lib["sequence"], lib["sgRNA"]))
near = {}
for seq, gid in exact.items():
    for i in range(20):
        for b in "ACGT":
            if b != seq[i]:
                v = seq[:i] + b + seq[i + 1:]
                near[v] = None if v in near and near[v] != gid else gid   # ambiguous -> None


def spacer(read):
    for r in (read, read.translate(COMP)[::-1]):
        k = r.find(ANCHOR)
        if k >= 20:
            return r[k - 20:k]
    return None


tab = {g: [0] * len(sids) for g in lib["sgRNA"]}
total = {}
for j, r in enumerate(samples.itertuples()):
    n = 0
    with gzip.open(inp / "fastq" / r.fastq_file, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 != 1:
                continue
            n += 1
            sp = spacer(line.strip())
            if sp is None:
                continue
            g = exact.get(sp) or near.get(sp)
            if g:
                tab[g][j] += 1
    total[r.sample_id] = n

counts = lib[["sgRNA", "gene"]].copy()
for j, s in enumerate(sids):
    counts[s] = [tab[g][j] for g in counts["sgRNA"]]
counts.to_csv(out / "counts.tsv", sep="\t", index=False)

ntc = (counts["gene"] == "NON_TARGETING").values
C = counts[sids].astype(float).values
ref = np.exp(np.log(C[ntc] + 1).mean(axis=1))                  # geometric mean per control guide
sf = np.median((C[ntc] + 1) / ref[:, None], axis=0)            # median-of-ratios on controls only
N = C / sf
t_idx = [i for i, r in enumerate(samples.itertuples()) if r.timepoint != "T0"]
c_idx = [i for i, r in enumerate(samples.itertuples()) if r.timepoint == "T0"]
lfc = np.log2((N[:, t_idx].mean(axis=1) + 1.0) / (N[:, c_idx].mean(axis=1) + 1.0))
lfc -= np.median(lfc[ntc])
counts["lfc"] = lfc

rng = np.random.default_rng(0)
null_guides = lfc[ntc]
rows = []
for gene, sub in counts[~ntc].groupby("gene"):
    k = len(sub)
    obs = float(np.median(sub["lfc"]))
    draws = np.median(rng.choice(null_guides, size=(20000, k), replace=True), axis=1)
    p = (np.sum(draws <= obs) + 1) / (len(draws) + 1)
    rows.append((gene, obs, p))
res = pd.DataFrame(rows, columns=["gene", "lfc", "pvalue"]).sort_values("pvalue").reset_index(drop=True)
m = len(res)
adj = res["pvalue"].values * m / (np.arange(m) + 1)
res["fdr"] = np.minimum(np.minimum.accumulate(adj[::-1])[::-1], 1.0)
res["is_hit"] = (res["fdr"] < 0.05).astype(int)
res.to_csv(out / "gene_summary.tsv", sep="\t", index=False)
(out / "qc.json").write_text(json.dumps(
    {"samples": {s: {"total_reads": int(total[s]), "mapped_reads": int(counts[s].sum())} for s in sids}}, indent=1))
