#!/usr/bin/env python3
"""Independent second solution, level 4 (no cutadapt, no MAGeCK, no scipy). Shares no step with the
reference pipeline:

  counts   pure python: scaffold anchor on either strand, exact then 1-mismatch spacer lookup
  arms     reference-gene signature per day-14 library (reference sensitisers minus reference
           resistance genes), no clustering
  calls    no replicate variance model: each sgRNA log ratio is standardised by the robust spread
           of the log ratios of sgRNAs with similar abundance (as drugZ does), Stouffer per gene
  audit    the inherited tables are explained numerically: each inherited count column is matched
           to (FASTQ file, read layout) by re-counting every file in every layout, and the inherited
           gene table is reproduced from the inherited count table under each combination of
           normalisation / baseline / sign. The inherited scripts are not executed.
"""
import gzip
import itertools
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

inp, out, work = Path("/task/input"), Path("/task/output"), Path("/tmp/l4_alt")
ho = inp / "handover"
out.mkdir(parents=True, exist_ok=True)
work.mkdir(parents=True, exist_ok=True)
ANCHOR = "GTTTTAGAGCTA"
COMP = str.maketrans("ACGTN", "TGCAN")
ITEMS = ["library_release", "sample_columns", "stagger_reads", "read2_library",
         "arm_labels", "control_normalisation", "drug_baseline", "contrast_sign"]

lib = pd.read_csv(inp / "library.tsv", sep="\t")
samples = pd.read_csv(inp / "samples.tsv", sep="\t")
ref = pd.read_csv(inp / "reference_genes.tsv", sep="\t")
sids = list(samples["sample_id"])
sheet = dict(zip(samples["sample_id"], samples["arm"]))
idx = {g: k for k, g in enumerate(lib["sgRNA"])}
exact = dict(zip(lib["sequence"], range(len(lib))))
near = {}
for seq, k in exact.items():
    for i in range(20):
        for b in "ACGT":
            if b != seq[i]:
                v = seq[:i] + b + seq[i + 1:]
                near[v] = None if v in near and near[v] != k else k   # ambiguous -> None

# ---------------------------------------------------------------- counting, all layouts in one pass
LAYOUTS = ["anchor_fwd", "anchor_rev", "fixed_fwd", "fixed_rev"]
lay = {s: {m: np.zeros(len(lib), dtype=int) for m in LAYOUTS} for s in sids}
mine = {s: np.zeros(len(lib), dtype=int) for s in sids}
total = {}
for r in samples.itertuples():
    n = 0
    L, M = lay[r.sample_id], mine[r.sample_id]
    with gzip.open(inp / "fastq" / r.fastq_file, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 != 1:
                continue
            n += 1
            read = line.strip()
            rc = read.translate(COMP)[::-1]
            k = read.find(ANCHOR)
            sp_f = read[k - 20:k] if k >= 20 else None
            k = rc.find(ANCHOR)
            sp_r = rc[k - 20:k] if k >= 20 else None
            hit = None
            for sp in (sp_f, sp_r):
                if sp is not None:
                    hit = exact.get(sp)
                    if hit is None:
                        hit = near.get(sp)
                    if hit is not None:
                        break
            if hit is not None:
                M[hit] += 1
            for name, sp in (("anchor_fwd", sp_f), ("anchor_rev", sp_r), ("fixed_fwd", read[24:44]), ("fixed_rev", rc[15:35])):
                e = exact.get(sp) if sp else None
                if e is not None:
                    L[name][e] += 1
    total[r.sample_id] = n

counts = lib[["sgRNA", "gene"]].copy()
for s in sids:
    counts[s] = mine[s]
counts.to_csv(out / "counts.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- arms from the reference genes
ntc = (lib["gene"] == "NON_TARGETING").values
C = counts[sids].astype(float).values
geo = np.exp(np.log(C[ntc] + 1).mean(axis=1))
sf = np.median((C[ntc] + 1) / geo[:, None], axis=0)             # median of ratios, controls only
N = C / sf
day14 = [j for j, s in enumerate(sids) if sheet[s] != "T0"]
X = np.log2(N[:, day14] + 8.0)
X = X - X.mean(axis=1, keepdims=True)
cls_of = lib["gene"].map(dict(zip(ref["gene"], ref["expected_class"]))).values
sig = X[cls_of == "sensitiser"].mean(axis=0) - X[cls_of == "resistance"].mean(axis=0)
arms = {s: "T0" for s in sids}
for j, v in zip(day14, sig):
    arms[sids[j]] = "drug" if v < 0 else "vehicle"
T0 = [j for j, s in enumerate(sids) if arms[s] == "T0"]
V = [j for j, s in enumerate(sids) if arms[s] == "vehicle"]
D = [j for j, s in enumerate(sids) if arms[s] == "drug"]
print("signature:", dict(zip([sids[j] for j in day14], np.round(sig, 2))), "arms:", arms)


# ---------------------------------------------------------------- effect sizes
def guide_lfc(Nm, a, b, pc):
    return np.log2((Nm[:, a].mean(axis=1) + pc) / (Nm[:, b].mean(axis=1) + pc))


genes = sorted(set(lib["gene"]) - {"NON_TARGETING"})
gene_arr = lib["gene"].values


def gene_median(v, garr=gene_arr, glist=genes):
    s = pd.Series(v).groupby(garr).median()
    return s.reindex(glist).values


lf = guide_lfc(N, V, T0, 1.0)
lf -= np.median(lf[ntc])
ld = guide_lfc(N, D, V, 1.0)
ld -= np.median(ld[ntc])
res = pd.DataFrame({"gene": genes, "lfc_fitness": gene_median(lf), "lfc_drug": gene_median(ld)})

# ---------------------------------------------------------------- calls
# Spread of the sgRNA log ratios among sgRNAs of similar abundance (robust, so the few modifiers do
# not inflate it); sgRNA z-scores are combined per gene and rescaled to the bulk of the genes.
ab = np.log2(N[:, D + V].mean(axis=1) + 1.0)
o = np.argsort(ab)
lo = ld[o]
sd = np.empty(len(o))
for i in range(len(o)):
    w = lo[max(0, i - 75):i + 76]
    sd[i] = 1.4826 * np.median(np.abs(w - np.median(w)))
z = np.empty(len(o))
z[o] = lo / sd
gz = pd.Series(z).groupby(gene_arr).agg(lambda v: v.sum() / math.sqrt(len(v))).reindex(genes).values
scale = 1.4826 * np.median(np.abs(gz - np.median(gz)))
gz = (gz - np.median(gz)) / max(scale, 1.0)
pv = np.array([math.erfc(abs(v) / math.sqrt(2)) for v in gz])     # two-sided normal p-value
order = np.argsort(pv)
adj = pv[order] * len(pv) / (np.arange(len(pv)) + 1)
fdr = np.empty(len(pv))
fdr[order] = np.minimum(np.minimum.accumulate(adj[::-1])[::-1], 1.0)
res["z_drug"], res["fdr"] = gz, fdr
res["drug_class"] = np.where(fdr < 0.05, np.where(gz < 0, "sensitiser", "resistance"), "none")
res.to_csv(out / "gene_summary.tsv", sep="\t", index=False)
called = dict(zip(res["gene"], res["drug_class"]))
exp = dict(zip(ref["gene"], ref["expected_class"]))
print("release criterion:", sum(called[x] == exp[x] for x in exp if exp[x] != "none"), "reference modifiers right,",
      sum(called[x] != "none" for x in exp if exp[x] == "none"), "reference none genes called")

# ---------------------------------------------------------------- audit
draft = pd.read_csv(ho / "counts_draft.tsv", sep="\t").set_index("sgRNA").reindex(lib["sgRNA"])
audit = {"library_release": "ok" if (draft["gene"].values == lib["gene"].values).all() else "defective"}
combos = {"anchor_both": ("anchor_fwd", "anchor_rev"), "anchor_fwd": ("anchor_fwd",),
          "fixed_both": ("fixed_fwd", "fixed_rev"), "fixed_fwd": ("fixed_fwd",)}
col_sum = draft[sids].sum(axis=0)
empty = [s for s in sids if col_sum[s] <= 0.01 * col_sum.median()]
source, layout = {}, {}
for s in sids:
    if s in empty:
        continue
    d = draft[s].values
    best = None
    for f in sids:
        for name, parts in combos.items():
            v = sum(lay[f][m] for m in parts)
            dist = np.abs(d - v).sum() / max(d.sum(), 1)
            if best is None or dist < best[0]:
                best = (dist, f, name)
    source[s], layout[s] = best[1], best[2]
print("inherited columns come from:", source, layout, "empty:", empty)
audit["sample_columns"] = "ok" if all(source[s] == s for s in source) else "defective"
audit["stagger_reads"] = "defective" if any(v.startswith("fixed") for v in layout.values()) else "ok"
audit["read2_library"] = "defective" if empty else "ok"
audit["arm_labels"] = "ok" if arms == sheet else "defective"

# reproduce the inherited gene table from the inherited count table (sheet arms, as that script does)
gd = pd.read_csv(ho / "gene_draft.tsv", sep="\t").set_index("gene")
dC = draft[sids].astype(float).values
dgene = draft["gene"].values
dntc = dgene == "NON_TARGETING"
keep = [j for j, s in enumerate(sids) if s not in empty]
grp = {k: [j for j in keep if sheet[sids[j]] == k] for k in ("T0", "vehicle", "drug")}
dgenes = list(gd.index)
best = None
for norm, base, sign in itertools.product(("control", "library"), ("vehicle", "T0"), (1, -1)):
    size = np.ones(len(sids))
    if norm == "control":
        tot = dC[dntc].sum(axis=0)
        size[keep] = tot[keep] / tot[keep].mean()
    else:
        with np.errstate(divide="ignore", invalid="ignore"):
            lg = np.log(np.where(dC[:, keep] > 0, dC[:, keep], np.nan))
            size[keep] = np.exp(np.nanmedian(lg - np.nanmean(lg, axis=1, keepdims=True), axis=0))
    Nd = dC / size
    f_ = guide_lfc(Nd, grp["vehicle"], grp["T0"], 0.5)
    d_ = sign * guide_lfc(Nd, grp["drug"], grp[base], 0.5)
    if norm == "control":
        f_, d_ = f_ - np.median(f_[dntc]), d_ - np.median(d_[dntc])
    err = np.median(np.abs(gene_median(f_, dgene, dgenes) - gd["lfc_fitness"].values)) + \
        np.median(np.abs(gene_median(d_, dgene, dgenes) - gd["lfc_drug"].values))
    if best is None or err < best[0]:
        best = (err, norm, base, sign)
print("inherited gene table reproduced with", best)
audit["control_normalisation"] = "ok" if best[1] == "control" else "defective"
audit["drug_baseline"] = "ok" if best[2] == "vehicle" else "defective"
audit["contrast_sign"] = "ok" if best[3] == 1 else "defective"
pd.DataFrame({"item": ITEMS, "verdict": [audit[i] for i in ITEMS]}).to_csv(out / "audit.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- erratum, qc
mailed = pd.read_csv(ho / "hits_mailed.tsv", sep="\t")["gene"]
pd.DataFrame({"gene": mailed, "verdict": ["confirmed" if called[x] != "none" else "retracted" for x in mailed]}).to_csv(
    out / "erratum.tsv", sep="\t", index=False)
(out / "qc.json").write_text(json.dumps(
    {"samples": {s: {"total_reads": int(total[s]), "mapped_reads": int(counts[s].sum()), "arm": arms[s]} for s in sids}}, indent=1))
