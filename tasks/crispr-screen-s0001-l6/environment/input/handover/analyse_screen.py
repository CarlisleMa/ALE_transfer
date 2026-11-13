#!/usr/bin/env python3
"""analyse_screen.py - gene table and hit list for the VAL-2 drug screen (K.B.)

usage: analyse_screen.py COUNTS.tsv SAMPLES.tsv OUT_GENES.tsv OUT_HITS.tsv
"""
import sys

import numpy as np

TOP = 40                             # genes that go into the hit list
counts_path, sheet_path, out_genes, out_hits = sys.argv[1:5]

rows = [line.rstrip("\n").split("\t") for line in open(counts_path) if line.strip()]
cols, rows = rows[0][2:], rows[1:]
gene = np.array([r[1] for r in rows])
C = np.array([[float(v) for v in r[2:]] for r in rows])
arm = {}
for line in list(open(sheet_path))[1:]:
    f = line.rstrip("\n").split("\t")
    if len(f) > 1:
        arm[f[0]] = f[1]

keep = [j for j in range(len(cols)) if C[:, j].sum() > 0]        # a library without reads is left out
t0 = [j for j in keep if arm[cols[j]] == "T0"]
veh = [j for j in keep if arm[cols[j]] == "vehicle"]
drug = [j for j in keep if arm[cols[j]] == "drug"]
ntc = gene == "NON_TARGETING"

# size factors: median of ratios to the geometric mean, all sgRNAs
size = np.ones(len(cols))
with np.errstate(divide="ignore", invalid="ignore"):
    logc = np.log(np.where(C[:, keep] > 0, C[:, keep], np.nan))
    size[keep] = np.exp(np.nanmedian(logc - np.nanmean(logc, axis=1, keepdims=True), axis=0))
N = C / size


def lfc(a, b):
    return np.log2((N[:, a].mean(axis=1) + 0.5) / (N[:, b].mean(axis=1) + 0.5))


fit = lfc(veh, t0)
drg = lfc(drug, veh)

# gene level: median over the sgRNAs; z against the spread of the non-targeting sgRNAs
null_sd = np.std(drg[ntc]) / 2.0          # a gene is the median of four sgRNAs
table = []
for g in sorted(set(gene[~ntc])):
    m = gene == g
    lf, ld = float(np.median(fit[m])), float(np.median(drg[m]))
    z = ld / null_sd
    cls = "none" if abs(z) < 3 else ("sensitiser" if z < 0 else "resistance")
    table.append((g, lf, ld, z, cls))

with open(out_genes, "w") as out:
    out.write("gene\tlfc_fitness\tlfc_drug\tz_drug\tdrug_class\n")
    for g, lf, ld, z, cls in table:
        out.write("%s\t%.4f\t%.4f\t%.2f\t%s\n" % (g, lf, ld, z, cls))

with open(out_hits, "w") as out:
    out.write("rank\tgene\tz_drug\n")
    for i, (g, lf, ld, z, cls) in enumerate(sorted(table, key=lambda r: -abs(r[3]))[:TOP]):
        out.write("%d\t%s\t%.2f\n" % (i + 1, g, z))
