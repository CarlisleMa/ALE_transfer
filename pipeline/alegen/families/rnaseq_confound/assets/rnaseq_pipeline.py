#!/usr/bin/env python3
"""Reference analysis: salmon -> gene-level counts -> sample identity check -> pyDESeq2.

Defaults are the oracle. The non-default flags replay plausible-but-wrong analyses for the
validator's attack suite (a careless analyst still delivers every file):
  --batch no          condition-only design
  --labels given      sample sheet used as delivered
  --labels cond_only  the two mismatched samples get each other's condition, nothing else
  --labels exclude    mismatched samples are left out (a second valid resolution)
"""
import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

Y_GENES = ["RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY"]

ap = argparse.ArgumentParser()
ap.add_argument("--input", default="/task/input")
ap.add_argument("--output", default="/task/output")
ap.add_argument("--work", default="/tmp/rnaseq_work")
ap.add_argument("--threads", type=int, default=2)
ap.add_argument("--labels", choices=["fix", "exclude", "given", "cond_only"], default="fix")
ap.add_argument("--batch", choices=["yes", "no"], default="yes")
ap.add_argument("--counts", default=None, help="skip quantification and read this gene x sample table")
a = ap.parse_args()
inp, out, work = Path(a.input), Path(a.output), Path(a.work)
out.mkdir(parents=True, exist_ok=True)
work.mkdir(parents=True, exist_ok=True)

samples = pd.read_csv(inp / "samples.tsv", sep="\t", dtype=str)
t2g = pd.read_csv(inp / "reference" / "tx2gene.tsv", sep="\t", dtype=str)
genes = list(dict.fromkeys(t2g["gene_id"]))
name_of = dict(zip(t2g["gene_id"], t2g["gene_name"]))
sids = list(samples["sample_id"])


def sh(cmd):
    r = subprocess.run([str(c) for c in cmd], capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"command failed: {' '.join(map(str, cmd))}\n{r.stderr[-3000:]}")
    return r


# ----------------------------------------------------------------------------- quantification
if a.counts:
    counts = pd.read_csv(a.counts, sep="\t", index_col=0)[sids].reindex(genes).fillna(0.0)
else:
    idx = work / "salmon_index"
    if not (idx / "versionInfo.json").exists():
        sh(["salmon", "index", "-t", inp / "reference" / "transcripts.fa.gz", "-i", idx, "-k", "31",
            "-p", a.threads])
    tx_gene = dict(zip(t2g["transcript_id"], t2g["gene_id"]))
    cols = {}
    for r in samples.itertuples():
        q = work / "quant" / r.sample_id
        if not (q / "quant.sf").exists():
            sh(["salmon", "quant", "-i", idx, "-l", "A", "-r", inp / "fastq" / r.fastq_file,
                "--validateMappings", "-p", a.threads, "-o", q])
        sf = pd.read_csv(q / "quant.sf", sep="\t")
        sf["gene"] = sf["Name"].map(lambda n: tx_gene.get(n, tx_gene.get(n.split(".")[0])))
        cols[r.sample_id] = sf.groupby("gene")["NumReads"].sum()
    counts = pd.DataFrame(cols).reindex(genes).fillna(0.0)[sids]
counts.index.name = "gene_id"
counts.round(3).to_csv(out / "gene_counts.tsv", sep="\t")

# ----------------------------------------------------------------------------- identity check
cpm = counts / counts.sum(axis=0) * 1e6
sym = pd.Series(name_of).reindex(counts.index)
y_cpm = cpm[sym.isin(Y_GENES).values].sum(axis=0)
x_cpm = cpm[(sym == "XIST").values].sum(axis=0)
inferred = {s: ("M" if y_cpm[s] > x_cpm[s] else "F") for s in sids}
meta = samples.set_index("sample_id").copy()
meta["inferred_sex"] = pd.Series(inferred)
mismatch = [s for s in sids if meta.loc[s, "donor_sex"] != inferred[s]]

used = meta[["condition", "library_prep_batch"]].copy()
used.columns = ["condition", "batch"]
excluded = []
if mismatch and a.labels in ("fix", "cond_only"):
    pair_ok = (len(mismatch) == 2
               and meta.loc[mismatch[0], "donor_sex"] == inferred[mismatch[1]]
               and meta.loc[mismatch[1], "donor_sex"] == inferred[mismatch[0]])
    if pair_ok:   # the two files carry each other's library: exchange the records
        s1, s2 = mismatch
        cols_x = ["condition", "batch"] if a.labels == "fix" else ["condition"]
        tmp = used.loc[s1, cols_x].copy()
        used.loc[s1, cols_x] = used.loc[s2, cols_x].values
        used.loc[s2, cols_x] = tmp.values
    else:         # identity cannot be reconstructed
        excluded = mismatch
elif mismatch and a.labels == "exclude":
    excluded = mismatch

qc = pd.DataFrame({"sample_id": sids,
                   "inferred_sex": [inferred[s] for s in sids],
                   "recorded_sex": [meta.loc[s, "donor_sex"] for s in sids],
                   "condition_used": ["excluded" if s in excluded else used.loc[s, "condition"] for s in sids],
                   "batch_used": ["excluded" if s in excluded else used.loc[s, "batch"] for s in sids],
                   "xist_cpm": [round(float(x_cpm[s]), 1) for s in sids],
                   "chrY_cpm": [round(float(y_cpm[s]), 1) for s in sids]})
qc.to_csv(out / "sample_qc.tsv", sep="\t", index=False)

# ----------------------------------------------------------------------------- differential expression
from pydeseq2.dds import DeseqDataSet   # noqa: E402
from pydeseq2.ds import DeseqStats      # noqa: E402

keep = [s for s in sids if s not in excluded]
md = used.loc[keep].copy()
md["batch"] = ["b" + str(b).replace("-", "_").replace(" ", "_") for b in md["batch"]]
mat = counts[keep].T.round().astype(int)
use_batch = a.batch == "yes" and md["batch"].nunique() > 1
design = "~batch + condition" if use_batch else "~condition"
dds = DeseqDataSet(counts=mat, metadata=md, design=design, quiet=True, n_cpus=a.threads)
dds.deseq2()
st = DeseqStats(dds, contrast=["condition", "treated", "vehicle"], alpha=0.05, quiet=True, n_cpus=a.threads)
st.summary()
res = st.results_df.reindex(genes)
de = pd.DataFrame({"gene_id": genes,
                   "gene_name": [name_of[g] for g in genes],
                   "base_mean": res["baseMean"].fillna(0.0).round(2).values,
                   "log2fc": res["log2FoldChange"].fillna(0.0).round(4).values,
                   "pvalue": res["pvalue"].values,
                   "padj": res["padj"].values})
de["significant"] = (de["padj"] < 0.05).astype(int)
de.to_csv(out / "de_results.tsv", sep="\t", index=False, na_rep="NA")
print(f"design={design} n={len(keep)} mismatch={mismatch} excluded={excluded} "
      f"significant={int(de['significant'].sum())}")
