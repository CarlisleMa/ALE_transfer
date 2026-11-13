#!/usr/bin/env python3
"""Independent second solution. Shares no tool with the oracle:

  quantification  two disjoint 25-mers per read looked up in a k-mer -> gene table (numpy), reads whose
                  seeds point to different genes are dropped; no salmon, no EM
  identity check  each sample is typed from XIST and chrY reads; samples whose type contradicts the
                  sample sheet are EXCLUDED (the oracle exchanges their records instead)
  statistics      limma-voom re-implemented with numpy/scipy/statsmodels: TMM factors, log-CPM,
                  lowess mean-variance trend -> precision weights -> weighted least squares on
                  ~ batch + condition, empirical-Bayes moderated t, Benjamini-Hochberg
"""
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import special, stats
from statsmodels.nonparametric.smoothers_lowess import lowess

import argparse

ap = argparse.ArgumentParser()
ap.add_argument("--input", default="/task/input")
ap.add_argument("--output", default="/task/output")
ap.add_argument("--counts", default=None, help="skip quantification and read this gene x sample table")
a = ap.parse_args()
inp, out = Path(a.input), Path(a.output)
out.mkdir(parents=True, exist_ok=True)
K = 25
Y_GENES = {"RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY"}

samples = pd.read_csv(inp / "samples.tsv", sep="\t", dtype=str)
t2g = pd.read_csv(inp / "reference" / "tx2gene.tsv", sep="\t", dtype=str)
genes = list(dict.fromkeys(t2g["gene_id"]))
gidx = {g: i for i, g in enumerate(genes)}
tx_gene = dict(zip(t2g["transcript_id"], t2g["gene_id"]))
name_of = dict(zip(t2g["gene_id"], t2g["gene_name"]))
sids = list(samples["sample_id"])

CODE = np.full(256, 4, dtype=np.uint8)
for i, b in enumerate(b"ACGT"):
    CODE[b] = i
POW = (np.uint64(4) ** np.arange(K - 1, -1, -1, dtype=np.uint64))


def kmer_codes(arr2d):
    """(n, K) uint8 base codes -> uint64 k-mer codes; rows containing N get the sentinel 2^63."""
    bad = (arr2d > 3).any(axis=1)
    v = (arr2d.astype(np.uint64) * POW).sum(axis=1)
    v[bad] = np.uint64(1) << np.uint64(63)
    return v


# ---- k-mer -> gene table (both strands) ------------------------------------------------------
keys, vals = [], []
name, seq = None, []


def flush():
    if name is None:
        return
    g = gidx[tx_gene.get(name, tx_gene.get(name.split(".")[0]))]
    s = CODE[np.frombuffer("".join(seq).encode(), dtype=np.uint8)]
    for strand in (s, (3 - s[::-1]).astype(np.uint8)):
        strand = np.where(strand > 3, 4, strand).astype(np.uint8)
        win = np.lib.stride_tricks.sliding_window_view(strand, K)
        keys.append(kmer_codes(win))
        vals.append(np.full(len(win), g, dtype=np.int32))


if a.counts:      # development switch: statistics only, on a supplied count table
    keys, vals = [np.zeros(1, dtype=np.uint64)], [np.zeros(1, dtype=np.int32)]
else:
    with gzip.open(inp / "reference" / "transcripts.fa.gz", "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip().upper())
    flush()
keys, vals = np.concatenate(keys), np.concatenate(vals)
order = np.argsort(keys, kind="stable")
keys, vals = keys[order], vals[order]
first = np.concatenate([[True], keys[1:] != keys[:-1]])
grp = np.cumsum(first) - 1
lo = np.minimum.reduceat(vals, np.flatnonzero(first))
hi = np.maximum.reduceat(vals, np.flatnonzero(first))
ukeys = keys[first]
ugene = np.where(lo == hi, lo, -1).astype(np.int32)      # -1: k-mer shared between genes


def lookup(codes):
    pos = np.searchsorted(ukeys, codes)
    pos[pos >= len(ukeys)] = 0
    hit = ukeys[pos] == codes
    return np.where(hit, ugene[pos], -2)                  # -2: not in the reference


counts = np.zeros((len(genes), len(sids)), dtype=np.int64)
if a.counts:
    counts = pd.read_csv(a.counts, sep="\t", index_col=0)[sids].reindex(genes).fillna(0).round().astype(np.int64).values
for j, r in enumerate(samples.itertuples() if not a.counts else []):
    with gzip.open(inp / "fastq" / r.fastq_file, "rb") as fh:
        reads = fh.read().split(b"\n")[1::4]
    L = min(len(x) for x in reads)
    arr = CODE[np.frombuffer(b"".join(x[:L] for x in reads), dtype=np.uint8).reshape(len(reads), L)]
    g1 = lookup(kmer_codes(arr[:, :K]))
    g2 = lookup(kmer_codes(arr[:, L - K:]))
    gene = np.where(g1 >= 0, g1, g2)
    ok = (gene >= 0) & ~((g1 >= 0) & (g2 >= 0) & (g1 != g2))
    counts[:, j] = np.bincount(gene[ok], minlength=len(genes))

cdf = pd.DataFrame(counts, index=pd.Index(genes, name="gene_id"), columns=sids)
cdf.to_csv(out / "gene_counts.tsv", sep="\t")

# ---- identity check: exclude what contradicts the sheet --------------------------------------
sym = np.array([name_of[g] for g in genes])
xist = counts[sym == "XIST"].sum(axis=0) + 0.5
chry = counts[np.isin(sym, list(Y_GENES))].sum(axis=0) + 0.5
inferred = np.where(np.log2(chry / xist) > 0, "M", "F")
recorded = samples["donor_sex"].values
keep = inferred == recorded
qc = pd.DataFrame({"sample_id": sids, "inferred_sex": inferred,
                   "condition_used": np.where(keep, samples["condition"].values, "excluded")})
qc.to_csv(out / "sample_qc.tsv", sep="\t", index=False)

# ---- voom-style linear model -------------------------------------------------------------------
C = counts[:, keep].astype(float)
cond = (samples["condition"].values[keep] == "treated").astype(float)
batch = samples["library_prep_batch"].values[keep]
cols = [np.ones(len(cond))]
for lev in sorted(set(batch))[1:]:
    cols.append((batch == lev).astype(float))
cols.append(cond)
X = np.column_stack(cols)
n, p = X.shape
tested = (C >= 10).sum(axis=1) >= min(int(cond.sum()), int(n - cond.sum()))
Y = C[tested]


def tmm_factors(Y):
    N = Y.sum(axis=0)
    uq = np.quantile(Y / N, 0.75, axis=0)
    ref = int(np.argmin(np.abs(uq - uq.mean())))
    f = np.ones(Y.shape[1])
    for j in range(Y.shape[1]):
        a, b = Y[:, j], Y[:, ref]
        ok = (a > 0) & (b > 0)
        M = np.log2((a[ok] / N[j]) / (b[ok] / N[ref]))
        A = 0.5 * np.log2((a[ok] / N[j]) * (b[ok] / N[ref]))
        w = 1.0 / ((N[j] - a[ok]) / (N[j] * a[ok]) + (N[ref] - b[ok]) / (N[ref] * b[ok]))
        km = (M > np.quantile(M, 0.3)) & (M < np.quantile(M, 0.7)) & (A > np.quantile(A, 0.05)) & (A < np.quantile(A, 0.95))
        f[j] = 2.0 ** (np.sum(w[km] * M[km]) / np.sum(w[km])) if km.any() else 1.0
    return f / np.exp(np.mean(np.log(f)))


lib = Y.sum(axis=0) * tmm_factors(Y)
E = np.log2((Y + 0.5) / (lib + 1.0) * 1e6)
XtXi = np.linalg.inv(X.T @ X)
B = E @ X @ XtXi
res = E - B @ X.T
s = np.sqrt((res ** 2).sum(axis=1) / (n - p))
sx = E.mean(axis=1) + np.mean(np.log2(lib + 1.0)) - np.log2(1e6)
trend = lowess(np.sqrt(s), sx, frac=0.5, return_sorted=True)
tx_, ty_ = np.unique(trend[:, 0], return_index=True)
fit_cnt = (B @ X.T) + np.log2(lib + 1.0) - np.log2(1e6)
W = 1.0 / np.interp(fit_cnt, tx_, trend[ty_, 1]) ** 4

beta = np.zeros(len(Y))
s2 = np.zeros(len(Y))
unsc = np.zeros(len(Y))
for i in range(len(Y)):
    w = W[i]
    A_ = np.linalg.inv((X * w[:, None]).T @ X)
    b = A_ @ (X * w[:, None]).T @ E[i]
    r = E[i] - X @ b
    beta[i], s2[i], unsc[i] = b[-1], np.sum(w * r * r) / (n - p), A_[-1, -1]

df = n - p
z = np.log(s2)
e = z - special.digamma(df / 2) + np.log(df / 2)
evar = np.var(e, ddof=1) - special.polygamma(1, df / 2)
if evar > 0:
    x = 0.5 + 1.0 / evar                              # trigamma inverse (Newton, as in limma)
    for _ in range(50):
        tri = special.polygamma(1, x)
        dif = tri * (1 - tri / evar) / special.polygamma(2, x)
        x += dif
        if -dif / x < 1e-8:
            break
    d0 = 2 * x
    s02 = np.exp(e.mean() + special.digamma(d0 / 2) - np.log(d0 / 2))
    post = (d0 * s02 + df * s2) / (d0 + df)
    df_tot = min(d0 + df, len(Y) * df)
else:
    post, df_tot = np.full(len(Y), np.exp(e.mean())), len(Y) * df
t = beta / np.sqrt(post * unsc)
pval = 2 * stats.t.sf(np.abs(t), df_tot)
o = np.argsort(pval)
adj = pval[o] * len(pval) / (np.arange(len(pval)) + 1)
padj = np.empty(len(pval))
padj[o] = np.minimum(np.minimum.accumulate(adj[::-1])[::-1], 1.0)

de = pd.DataFrame({"gene_id": genes, "gene_name": [name_of[g] for g in genes], "log2fc": 0.0,
                   "pvalue": np.nan, "padj": np.nan})
ti = np.flatnonzero(tested)
de.loc[ti, "log2fc"] = np.round(beta, 4)
de.loc[ti, "pvalue"] = pval
de.loc[ti, "padj"] = padj
de["significant"] = (de["padj"] < 0.05).astype(int)
de.to_csv(out / "de_results.tsv", sep="\t", index=False, na_rep="NA")
print(f"kept={int(keep.sum())} tested={int(tested.sum())} significant={int(de['significant'].sum())}")
