#!/usr/bin/env python3
"""Independent second solution for level 4. Shares no quantifier, no statistics and no audit logic
with the oracle:

  quantification  two disjoint 25-mers per read looked up in a k-mer -> gene table (numpy); no salmon
  identity        SOP criterion on XIST / chrY reads; the mixed library and every library whose sex
                  contradicts the sample sheet are EXCLUDED (the oracle exchanges the two records)
  statistics      limma-voom re-implemented with numpy/scipy/statsmodels (see the level-3 solve_alt.py)
  audit           the shipped matrix is compared with the k-mer counts (profile distance and depth),
                  not with a salmon re-run; the scripts are inspected statically (shlex / ast / regex);
                  the shipped DE script is executed on the k-mer counts of the accepted libraries to
                  see what it does, and on the shipped matrix to see whether the draft is its output
  erratum         significance in the voom analysis
"""
import ast
import gzip
import re
import shlex
import subprocess
import sys
import tempfile
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

# ---- identity control of the SOP: exclude the mixture and whatever contradicts the sheet ---------
sym = np.array([name_of[g] for g in genes])
xist = counts[sym == "XIST"].sum(axis=0).astype(float)
chry = counts[np.isin(sym, list(Y_GENES))].sum(axis=0).astype(float)
share = np.minimum(xist, chry) / np.maximum(xist + chry, 1e-9)
single = share < 0.02
inferred = np.where(~single, "mixed", np.where(xist > chry, "F", "M"))
recorded = samples["donor_sex"].values
keep = single & (inferred == recorded)
qc = pd.DataFrame({"sample_id": sids, "inferred_sex": inferred,
                   "condition_used": np.where(keep, samples["condition"].values, "excluded"),
                   "minor_marker_share": np.round(share, 4)})
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

# ---- audit of the previous analysis ----------------------------------------------------------------
prev = inp / "previous_analysis"
ITEMS = ["quant_library_type", "merge_read_counts", "matrix_sample_files", "matrix_all_reads",
         "de_session_term", "de_contrast_direction", "de_adjusted_pvalues", "draft_reproducible"]
verdict, why = {}, {}
if prev.is_dir() and not a.counts:
    tmp = Path(tempfile.mkdtemp(prefix="alt_audit_"))
    # library type: what the quantifier was told
    lib = None
    for line in (prev / "scripts" / "01_quantify.sh").read_text().replace("\\\n", " ").splitlines():
        try:
            tok = shlex.split(line, comments=True)
        except ValueError:
            continue
        if "quant" in tok:
            for i, t in enumerate(tok[:-1]):
                if t in ("-l", "--libType"):
                    lib = tok[i + 1]
    verdict["quant_library_type"] = "ok" if lib in ("A", "U") else "defective"
    why["quant_library_type"] = f"salmon quant -l {lib}; libraries are unstranded single-end"

    M = pd.read_csv(prev / "gene_counts_merged.tsv", sep="\t", index_col=0)[sids].reindex(genes).fillna(0.0).values
    cs = M.sum(axis=0)
    # read counts or normalised abundances
    src02 = (prev / "scripts" / "02_merge_counts.py").read_text()
    used_cols = set(re.findall(r"\[\s*[\"'](NumReads|TPM|EffectiveLength|Length)[\"']\s*\]", src02))
    flat = cs.max() / max(cs.min(), 1e-9) < 1.001          # every column has the same total
    verdict["merge_read_counts"] = "ok" if used_cols == {"NumReads"} and not flat else "defective"
    why["merge_read_counts"] = f"quant.sf column(s) {sorted(used_cols)}; column sums {cs.min():.0f}-{cs.max():.0f}"

    # shipped matrix against the k-mer counts: profile distance after removing depth and gene scale
    expr = (counts.mean(axis=1) >= 20) & ~np.isin(sym, list(Y_GENES) + ["XIST"])

    def profile(mat):
        v = np.log2(mat[expr] + 1.0)
        v = v - np.median(v, axis=0, keepdims=True)
        v = v - v.mean(axis=1, keepdims=True)
        return v - np.median(v, axis=0, keepdims=True)

    pm, pc = profile(M), profile(counts.astype(float))
    D = np.array([[np.median(np.abs(pm[:, j] - pc[:, k])) for k in range(len(sids))] for j in range(len(sids))])
    best = D.argmin(axis=1)
    wrong = [f"{sids[j]} looks like {sids[best[j]]} ({D[j, best[j]]:.3f} vs {D[j, j]:.3f})"
             for j in range(len(sids)) if best[j] != j]
    verdict["matrix_sample_files"] = "defective" if wrong else "ok"
    why["matrix_sample_files"] = "; ".join(wrong)
    own = [j for j in range(len(sids)) if best[j] == j]
    if not flat:      # counts: a partial file shows as a depth deficit against the other columns
        ratio = cs / counts.sum(axis=0)
        ref = np.median(ratio[own])
        part = [f"{sids[j]} has {ratio[j] / ref:.2f} of the expected depth" for j in own if ratio[j] < 0.75 * ref]
    else:             # abundances carry no depth: fall back to the noise level of the column
        ref = np.median(D[own, own])
        part = [f"{sids[j]} is {D[j, j] / ref:.1f} x noisier than the other columns" for j in own if D[j, j] > 2.5 * ref]
    verdict["matrix_all_reads"] = "defective" if part else "ok"
    why["matrix_all_reads"] = "; ".join(part)

    # the DE script: read the model, then run it on clean data to see what it reports
    src03 = (prev / "scripts" / "03_de.py").read_text()
    formula = ""
    for node in ast.walk(ast.parse(src03)):
        if isinstance(node, ast.keyword) and node.arg == "design" and isinstance(node.value, ast.Constant):
            formula = str(node.value.value)
    terms = {t.strip() for t in formula.replace("~", "").split("+")} - {"", "1", "condition"}
    verdict["de_session_term"] = "ok" if terms else "defective"
    why["de_session_term"] = f"design {formula!r}"

    def run_de(matrix, sheet, dst):
        r = subprocess.run([sys.executable, str(prev / "scripts" / "03_de.py"), str(matrix), str(dst)]
                           + ([str(sheet)] if sheet else []), capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit("03_de.py failed:\n" + r.stderr[-2000:])
        return pd.read_csv(dst, sep="\t").set_index("gene_id").reindex(genes)

    samples[keep].to_csv(tmp / "samples_accepted.tsv", sep="\t", index=False)
    clean = run_de(out / "gene_counts.tsv", tmp / "samples_accepted.tsv", tmp / "de_clean.tsv")
    mine = de.set_index("gene_id")
    clear = mine.index[(mine["padj"] < 0.01) & (mine["log2fc"].abs() >= 1.0)]
    agree = float(np.mean(np.sign(clean.loc[clear, "log2fc"].fillna(0.0)) == np.sign(mine.loc[clear, "log2fc"])))
    verdict["de_contrast_direction"] = "ok" if agree > 0.5 else "defective"
    why["de_contrast_direction"] = f"sign agreement with the voom analysis on {len(clear)} clear responders: {agree:.2f}"
    t = clean.dropna(subset=["pvalue", "padj"])
    raw = bool(len(t)) and float(np.max(np.abs(t["pvalue"] - t["padj"]))) < 1e-12
    verdict["de_adjusted_pvalues"] = "defective" if raw else "ok"
    why["de_adjusted_pvalues"] = "padj is the raw p-value" if raw else "padj >= pvalue, not identical"

    redo = run_de(prev / "gene_counts_merged.tsv", None, tmp / "de_redo.tsv")
    draft = pd.read_csv(prev / "de_draft.tsv", sep="\t").set_index("gene_id").reindex(genes)
    dl = float(np.max(np.abs(redo["log2fc"].fillna(0.0) - draft["log2fc"].fillna(0.0))))
    verdict["draft_reproducible"] = "ok" if dl <= 0.002 else "defective"
    why["draft_reproducible"] = f"max |log2fc(re-run) - log2fc(draft)| = {dl:.4f}"

    pd.DataFrame({"item": ITEMS, "verdict": [verdict[i] for i in ITEMS],
                  "evidence": [why[i] for i in ITEMS]}).to_csv(out / "audit.tsv", sep="\t", index=False)

    # ---- erratum -----------------------------------------------------------------------------------
    mailed = pd.read_csv(prev / "mail_top_responders.tsv", sep="\t", dtype=str)
    sig = dict(zip(de["gene_id"], de["significant"]))
    pd.DataFrame({"gene_id": mailed["gene_id"],
                  "verdict": ["confirmed" if sig.get(g, 0) == 1 else "retracted" for g in mailed["gene_id"]]}
                 ).to_csv(out / "erratum.tsv", sep="\t", index=False)
    print("audit:", verdict)
