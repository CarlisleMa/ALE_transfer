#!/usr/bin/env python3
"""Reference analysis of the pooled CRISPR dropout screen (cutadapt -> mageck count -> control-based
statistics). The non-default flags replay plausible-but-wrong analyses for the validator's attack
suite; the oracle (solve.sh) uses the defaults.
"""
import argparse
import gzip
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

PREFIX_TAIL = "GGACGAAACACCG"            # 3' end of the U6 promoter (+G) directly upstream of the spacer
SCAFFOLD_HEAD = "GTTTTAGAGCTAGAAATAGC"   # 5' end of the sgRNA scaffold directly downstream

ap = argparse.ArgumentParser()
ap.add_argument("--input", default="/task/input")
ap.add_argument("--output", default="/task/output")
ap.add_argument("--work", default="/tmp/crispr_work")
ap.add_argument("--extract", choices=["cutadapt", "fixed", "mageck_auto"], default="cutadapt")
ap.add_argument("--no-revcomp", action="store_true")
ap.add_argument("--norm", choices=["control", "median"], default="control")
ap.add_argument("--test", choices=["control", "rra"], default="control")
a = ap.parse_args()

inp, out, work = Path(a.input), Path(a.output), Path(a.work)
out.mkdir(parents=True, exist_ok=True)
work.mkdir(parents=True, exist_ok=True)

lib = pd.read_csv(inp / "library.tsv", sep="\t")
samples = pd.read_csv(inp / "samples.tsv", sep="\t")
sids = list(samples["sample_id"])
fq = {r.sample_id: inp / "fastq" / r.fastq_file for r in samples.itertuples()}
treat = [r.sample_id for r in samples.itertuples() if r.timepoint != "T0"]
ctrl = [r.sample_id for r in samples.itertuples() if r.timepoint == "T0"]
ntc_ids = list(lib.loc[lib["gene"] == "NON_TARGETING", "sgRNA"])


def sh(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"command failed: {' '.join(map(str, cmd))}\n{r.stderr[-2000:]}")
    return r


def n_records(path):
    n = 0
    with gzip.open(path, "rt") as fh:
        for _ in fh:
            n += 1
    return n // 4


total_reads = {s: n_records(fq[s]) for s in sids}

# ----------------------------------------------------------------------------- guide counts
lib[["sgRNA", "sequence", "gene"]].to_csv(work / "library.csv", index=False, header=False)


def mageck_count(fastqs, trim):
    """One `mageck count` per sample; a sample MAGeCK cannot process is reported as all zeros
    (what an analyst gets after 'excluding the failed sample')."""
    cols = {}
    for s in sids:
        r = subprocess.run(["mageck", "count", "-l", str(work / "library.csv"), "-n", str(work / f"cnt_{s}"),
                            "--sample-label", s, "--fastq", str(fastqs[s]), "--trim-5", trim,
                            "--norm-method", "none"], capture_output=True, text=True)
        f = work / f"cnt_{s}.count.txt"
        if r.returncode == 0 and f.exists():
            t = pd.read_csv(f, sep="\t")
            cols[s] = dict(zip(t["sgRNA"], t[s]))
        else:
            cols[s] = {}
    df = pd.DataFrame({"sgRNA": lib["sgRNA"]})
    for s in sids:
        df[s] = df["sgRNA"].map(cols[s]).fillna(0).astype(int)
    return df


if a.extract == "cutadapt":
    trimmed = {}
    for s in sids:
        o = work / f"{s}.spacer.fastq.gz"
        cmd = ["cutadapt", "-j", "4", "-e", "0.1", "-g", f"{PREFIX_TAIL}...{SCAFFOLD_HEAD}",
               "--discard-untrimmed", "-m", "20", "-M", "20", "-o", str(o), str(fq[s])]
        if not a.no_revcomp:
            cmd.insert(1, "--revcomp")
        sh(cmd)
        trimmed[s] = o
    counts = mageck_count(trimmed, "0")
elif a.extract == "mageck_auto":   # raw reads, explicit list of stagger offsets, forward strand only
    counts = mageck_count(fq, ",".join(str(k) for k in range(24, 32)))
else:  # fixed offset, forward strand only
    seq2id = dict(zip(lib["sequence"], lib["sgRNA"]))
    tab = {g: [0] * len(sids) for g in lib["sgRNA"]}
    for j, s in enumerate(sids):
        with gzip.open(fq[s], "rt") as fh:
            for i, line in enumerate(fh):
                if i % 4 == 1:
                    g = seq2id.get(line[24:44])
                    if g:
                        tab[g][j] += 1
    counts = pd.DataFrame([[g, *v] for g, v in tab.items()], columns=["sgRNA", *sids])
counts["gene"] = counts["sgRNA"].map(dict(zip(lib["sgRNA"], lib["gene"])))

counts = lib[["sgRNA", "gene"]].merge(counts.drop(columns=["gene"]), on="sgRNA", how="left").fillna(0)
for s in sids:
    counts[s] = counts[s].astype(int)
counts[["sgRNA", "gene", *sids]].to_csv(out / "counts.tsv", sep="\t", index=False)

# ----------------------------------------------------------------------------- statistics
C = counts.set_index("sgRNA")[sids].astype(float)
is_ntc = counts["gene"].values == "NON_TARGETING"
alive = [s for s in sids if C[s].sum() > 0]
treat = [s for s in treat if s in alive]
ctrl = [s for s in ctrl if s in alive]
C = C[alive]
if a.norm == "control":
    tot = C[is_ntc].sum(axis=0)
    sf = tot / tot.mean()
else:  # median-of-ratios over the whole library (the usual default)
    logc = np.log(C.where(C > 0))
    geo = logc.mean(axis=1)
    sf = np.exp((logc.sub(geo, axis=0)).median(axis=0))
N = C / sf
glfc = np.log2((N[treat].mean(axis=1) + 0.5) / (N[ctrl].mean(axis=1) + 0.5))
if a.norm == "control":
    glfc = glfc - np.median(glfc[is_ntc])
counts["lfc"] = glfc.values
targ = counts[~is_ntc]
gene_lfc = targ.groupby("gene")["lfc"].median()

if a.test == "control":
    null = counts.loc[is_ntc, "lfc"].values
    pvals = targ.groupby("gene")["lfc"].apply(
        lambda v: stats.mannwhitneyu(v.values, null, alternative="less").pvalue)
    order = np.argsort(pvals.values)
    ranked = pvals.values[order] * len(pvals) / (np.arange(len(pvals)) + 1)
    fdr = np.minimum.accumulate(ranked[::-1])[::-1]
    fdr_s = pd.Series(np.minimum(fdr, 1.0), index=pvals.index[order])
else:
    C = C.copy()
    C.insert(0, "gene", counts["gene"].values)
    C.to_csv(work / "for_test.txt", sep="\t", index_label="sgRNA")
    cmd = ["mageck", "test", "-k", str(work / "for_test.txt"), "-t", ",".join(treat), "-c", ",".join(ctrl),
           "-n", str(work / "rra")]
    if a.norm == "control":
        (work / "ntc.txt").write_text("\n".join(ntc_ids) + "\n")
        cmd += ["--norm-method", "control", "--control-sgrna", str(work / "ntc.txt")]
    sh(cmd)
    g = pd.read_csv(work / "rra.gene_summary.txt", sep="\t").set_index("id")
    fdr_s = g["neg|fdr"]

summ = pd.DataFrame({"gene": gene_lfc.index, "lfc": gene_lfc.values})
summ["fdr"] = summ["gene"].map(fdr_s).fillna(1.0)
summ["is_hit"] = (summ["fdr"] < 0.05).astype(int)
summ.to_csv(out / "gene_summary.tsv", sep="\t", index=False)

qc = {"samples": {s: {"total_reads": int(total_reads[s]), "mapped_reads": int(counts[s].sum())} for s in sids}}
(out / "qc.json").write_text(json.dumps(qc, indent=1))
print(json.dumps({"hits": int(summ["is_hit"].sum()), "mapping_rate": {s: round(counts[s].sum() / total_reads[s], 3) for s in sids}}))
