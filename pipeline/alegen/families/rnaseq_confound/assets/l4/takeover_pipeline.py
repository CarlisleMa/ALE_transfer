#!/usr/bin/env python3
"""Level-4 reference: take over the RNA-seq project.

  1  fresh quantification of every FASTQ (salmon, automatic library type) -> gene_counts.tsv
  2  identity control of the SOP on the sex-chromosome control genes       -> sample_qc.tsv
  3  pyDESeq2, ~ session + condition, accepted libraries only               -> de_results.tsv
  4  audit of the previous analysis: the inherited scripts are re-run as shipped and their output is
     compared with the shipped files; the scripts' behaviour is compared with step 1-3 -> audit.tsv
  5  erratum for the mailed gene list                                       -> erratum.tsv

Defaults are the oracle. The other flags replay careless analyses for the validator's attack suite
(every one of them still delivers all files).
"""
import argparse
import gzip
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

Y_GENES = ["RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY"]
ITEMS = ["quant_library_type", "merge_read_counts", "matrix_sample_files", "matrix_all_reads",
         "de_session_term", "de_contrast_direction", "de_adjusted_pvalues", "draft_reproducible"]

ap = argparse.ArgumentParser()
ap.add_argument("--input", default="/task/input")
ap.add_argument("--output", default="/task/output")
ap.add_argument("--work", default="/tmp/takeover_work")
ap.add_argument("--threads", type=int, default=2)
ap.add_argument("--counts-from", choices=["own", "inherited"], default="own")
ap.add_argument("--libtype", default="A")
ap.add_argument("--quantity", choices=["NumReads", "TPM"], default="NumReads")
ap.add_argument("--dup", default=None, help="X:Y  column X is replaced by the quantification of sample Y")
ap.add_argument("--trunc", default=None, help="Z:N  column Z is quantified from the first N reads only")
ap.add_argument("--identity", choices=["sop", "binary"], default="sop")
ap.add_argument("--labels", choices=["fix", "exclude", "given"], default="fix")
ap.add_argument("--batch", choices=["yes", "no"], default="yes")
ap.add_argument("--contrast", choices=["normal", "reversed"], default="normal")
ap.add_argument("--padj", choices=["bh", "raw"], default="bh")
ap.add_argument("--audit", choices=["real", "all_ok", "all_defective"], default="real")
ap.add_argument("--audit-force", action="append", default=[], help="item=verdict, overrides the audit result")
ap.add_argument("--erratum", choices=["own", "all_confirmed", "all_retracted"], default="own")
ap.add_argument("--dev-counts", default=None, help="development: skip quantification, read this table")
a = ap.parse_args()
inp, out, work = Path(a.input), Path(a.output), Path(a.work)
prev = inp / "previous_analysis"
out.mkdir(parents=True, exist_ok=True)
work.mkdir(parents=True, exist_ok=True)

samples = pd.read_csv(inp / "samples.tsv", sep="\t", dtype=str)
t2g = pd.read_csv(inp / "reference" / "tx2gene.tsv", sep="\t", dtype=str)
genes = list(dict.fromkeys(t2g["gene_id"]))
name_of = dict(zip(t2g["gene_id"], t2g["gene_name"]))
tx_gene = dict(zip(t2g["transcript_id"], t2g["gene_id"]))
sids = list(samples["sample_id"])
fq_of = dict(zip(samples["sample_id"], samples["fastq_file"]))


def sh(cmd, **kw):
    r = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, **kw)
    if r.returncode != 0:
        sys.exit(f"command failed: {' '.join(map(str, cmd))}\n{r.stderr[-3000:]}")
    return r


def read_matrix(path):
    return pd.read_csv(path, sep="\t", index_col=0)[sids].reindex(genes).fillna(0.0)


def salmon_gene_column(fastq, qdir, libtype, quantity):
    if not (qdir / "quant.sf").exists():
        sh(["salmon", "quant", "-i", work / "salmon_index", "-l", libtype, "-r", fastq, "--validateMappings",
            "-p", a.threads, "-o", qdir])
    sf = pd.read_csv(qdir / "quant.sf", sep="\t")
    return sf.groupby(sf["Name"].map(tx_gene))[quantity].sum()


# ----------------------------------------------------------------------------- 1 quantification
if a.dev_counts:
    counts = read_matrix(a.dev_counts)
elif a.counts_from == "inherited":
    counts = read_matrix(prev / "gene_counts_merged.tsv")
else:
    if not (work / "salmon_index" / "versionInfo.json").exists():
        sh(["salmon", "index", "-t", inp / "reference" / "transcripts.fa.gz", "-i", work / "salmon_index",
            "-k", "31", "-p", a.threads])
    cols = {s: salmon_gene_column(inp / "fastq" / fq_of[s], work / "quant" / s, a.libtype, a.quantity) for s in sids}
    if a.dup:
        x, y = a.dup.split(":")
        cols[x] = cols[y]
    if a.trunc:
        z, n = a.trunc.split(":")
        part = work / f"{z}.first{n}.fastq.gz"
        with gzip.open(inp / "fastq" / fq_of[z], "rb") as src, gzip.open(part, "wb", compresslevel=1) as dst:
            for i, line in enumerate(src):
                if i >= 4 * int(n):
                    break
                dst.write(line)
        cols[z] = salmon_gene_column(part, work / "quant" / (z + "_partial"), a.libtype, a.quantity)
    counts = pd.DataFrame(cols).reindex(genes).fillna(0.0)[sids]
counts.index.name = "gene_id"
counts.round(3).to_csv(out / "gene_counts.tsv", sep="\t")

# ----------------------------------------------------------------------------- 2 identity control
sym = pd.Series(name_of).reindex(counts.index)
X = counts[(sym == "XIST").values].sum(axis=0)
Y = counts[sym.isin(Y_GENES).values].sum(axis=0)
minor = {s: float(min(X[s], Y[s]) / max(X[s] + Y[s], 1e-9)) for s in sids}
mixed = [s for s in sids if minor[s] >= 0.02] if a.identity == "sop" else []
inferred = {s: ("mixed" if s in mixed else ("F" if X[s] > Y[s] else "M")) for s in sids}
meta = samples.set_index("sample_id")
mismatch = [s for s in sids if s not in mixed and meta.loc[s, "donor_sex"] != inferred[s]]

used = meta[["condition", "library_prep_batch"]].copy()
used.columns = ["condition", "batch"]
excluded = list(mixed)
if mismatch and a.labels == "fix":
    pair_ok = (len(mismatch) == 2
               and meta.loc[mismatch[0], "donor_sex"] == inferred[mismatch[1]]
               and meta.loc[mismatch[1], "donor_sex"] == inferred[mismatch[0]])
    if pair_ok:   # the two files carry each other's library: exchange the records
        s1, s2 = mismatch
        tmp = used.loc[s1].copy()
        used.loc[s1] = used.loc[s2].values
        used.loc[s2] = tmp.values
    else:         # identity cannot be reconstructed
        excluded += mismatch
elif mismatch and a.labels == "exclude":
    excluded += mismatch

qc = pd.DataFrame({"sample_id": sids,
                   "inferred_sex": [inferred[s] for s in sids],
                   "recorded_sex": [meta.loc[s, "donor_sex"] for s in sids],
                   "condition_used": ["excluded" if s in excluded else used.loc[s, "condition"] for s in sids],
                   "session_used": ["excluded" if s in excluded else used.loc[s, "batch"] for s in sids],
                   "xist_reads": [round(float(X[s]), 1) for s in sids],
                   "chrY_reads": [round(float(Y[s]), 1) for s in sids],
                   "minor_marker_share": [round(minor[s], 4) for s in sids]})
qc.to_csv(out / "sample_qc.tsv", sep="\t", index=False)

# ----------------------------------------------------------------------------- 3 differential expression
from pydeseq2.dds import DeseqDataSet   # noqa: E402
from pydeseq2.ds import DeseqStats      # noqa: E402

keep = [s for s in sids if s not in excluded]
md = used.loc[keep].copy()
md["batch"] = ["s" + str(b).replace("-", "") for b in md["batch"]]
mat = counts[keep].T.round().astype(int)
use_batch = a.batch == "yes" and md["batch"].nunique() > 1
design = "~batch + condition" if use_batch else "~condition"
dds = DeseqDataSet(counts=mat, metadata=md, design=design, quiet=True, n_cpus=a.threads)
dds.deseq2()
contrast = ["condition", "treated", "vehicle"] if a.contrast == "normal" else ["condition", "vehicle", "treated"]
st = DeseqStats(dds, contrast=contrast, alpha=0.05, quiet=True, n_cpus=a.threads)
st.summary()
res = st.results_df.reindex(genes)
de = pd.DataFrame({"gene_id": genes,
                   "gene_name": [name_of[g] for g in genes],
                   "base_mean": res["baseMean"].fillna(0.0).round(2).values,
                   "log2fc": res["log2FoldChange"].fillna(0.0).round(4).values,
                   "pvalue": res["pvalue"].values,
                   "padj": (res["padj"] if a.padj == "bh" else res["pvalue"]).values})
de["significant"] = (de["padj"] < 0.05).astype(int)
de.to_csv(out / "de_results.tsv", sep="\t", index=False, na_rep="NA")

# ----------------------------------------------------------------------------- 4 audit
verdict = {i: "ok" for i in ITEMS}
evidence = {i: "" for i in ITEMS}


def flag(item, bad, why):
    verdict[item] = "defective" if bad else "ok"
    evidence[item] = why


if a.audit == "real" and not a.dev_counts:
    scripts = prev / "scripts"
    rw = work / "previous_rerun"
    if rw.exists():
        shutil.rmtree(rw)
    rw.mkdir(parents=True)
    # (a) quantification scripts exactly as shipped
    sh(["bash", scripts / "01_quantify.sh", rw])
    sh(["python3", scripts / "02_merge_counts.py", rw, rw / "matrix_rerun.tsv"])
    R = read_matrix(rw / "matrix_rerun.tsv")
    M = read_matrix(prev / "gene_counts_merged.tsv")
    # library type: does the shipped setting lose reads that the automatic setting keeps?
    rate = []
    for s in sids:
        j = json.loads((rw / "quant" / s / "aux_info" / "meta_info.json").read_text())
        rate.append(j["num_mapped"] / max(j["num_processed"], 1))
    own_rate = []
    for s in sids:
        f = work / "quant" / s / "aux_info" / "meta_info.json"
        if f.exists():
            j = json.loads(f.read_text())
            own_rate.append(j["num_mapped"] / max(j["num_processed"], 1))
    ratio = float(np.median(rate) / np.median(own_rate)) if own_rate else 1.0
    flag("quant_library_type", ratio < 0.9, f"reads assigned with the shipped setting / automatic setting = {ratio:.3f}")
    # read counts or normalised abundances: every column of a TPM table sums to the same constant
    tot = R.sum(axis=0)
    n_reads = np.array([json.loads((rw / "quant" / s / "aux_info" / "meta_info.json").read_text())["num_mapped"] for s in sids])
    is_counts = bool(np.allclose(tot.values, n_reads, rtol=0.02))
    flag("merge_read_counts", not is_counts,
         f"column sums of the re-run matrix {tot.min():.0f}-{tot.max():.0f}, reads assigned {n_reads.min()}-{n_reads.max()}")
    # shipped matrix against the re-run of the shipped scripts
    wrong_file, partial, notes = [], [], []
    for s in sids:
        if np.abs(M[s] - R[s]).sum() <= 1e-3 * max(R[s].sum(), 1.0):
            continue
        dist = {k: float(np.abs(M[s] - R[k]).sum() / max(R[k].sum(), 1.0)) for k in sids}
        best = min(dist, key=dist.get)
        if best != s and dist[best] <= 1e-3:
            wrong_file.append(s)
            notes.append(f"{s} equals the quantification of {best}")
            continue
        lm, lr = np.log2(M[s] + 1.0), np.log2(R + 1.0)
        cor = {k: float(np.corrcoef(lm, lr[k])[0, 1]) for k in sids}
        bestc = max(cor, key=cor.get)
        if bestc != s:
            wrong_file.append(s)
            notes.append(f"{s} resembles {bestc}")
        else:
            partial.append(s)
            notes.append(f"{s} holds {M[s].sum() / max(R[s].sum(), 1.0):.2f} of the re-run total")
    flag("matrix_sample_files", bool(wrong_file), "; ".join(n for n in notes if "equals" in n or "resembles" in n))
    flag("matrix_all_reads", bool(partial), "; ".join(n for n in notes if "holds" in n))
    # (b) DE script exactly as shipped, on the shipped matrix
    src = (scripts / "03_de.py").read_text()
    sh(["python3", scripts / "03_de.py", prev / "gene_counts_merged.tsv", rw / "de_rerun.tsv"])
    rr = pd.read_csv(rw / "de_rerun.tsv", sep="\t").set_index("gene_id").reindex(genes)
    dr = pd.read_csv(prev / "de_draft.tsv", sep="\t").set_index("gene_id").reindex(genes)
    dl = float(np.nanmax(np.abs(rr["log2fc"].fillna(0.0) - dr["log2fc"].fillna(0.0))))
    dp = float(np.nanmax(np.abs(rr["padj"].fillna(1.0) - dr["padj"].fillna(1.0))))
    flag("draft_reproducible", dl > 1e-3 or dp > 1e-3, f"max |log2fc difference| {dl:.4f}, max |padj difference| {dp:.4f}")
    m = re.search(r"design\s*=\s*[\"']([^\"']+)[\"']", src)
    formula = m.group(1) if m else ""
    terms = [t.strip() for t in formula.replace("~", "").split("+")]
    flag("de_session_term", not any(t not in ("", "condition", "1") for t in terms), f"design {formula!r}")
    # direction: sign of the shipped script's log2fc against the corrected analysis, on clear responders
    strong = de[(de["padj"] < 0.01) & (de["log2fc"].abs() >= 1.0)]["gene_id"] if a.contrast == "normal" else []
    if len(strong):
        agree = float(np.mean(np.sign(rr.loc[strong, "log2fc"].fillna(0.0).values)
                              == np.sign(de.set_index("gene_id").loc[strong, "log2fc"].values)))
    else:
        agree = 1.0
    flag("de_contrast_direction", agree < 0.5, f"sign agreement with the corrected analysis on {len(strong)} genes: {agree:.2f}")
    both = rr.dropna(subset=["pvalue", "padj"])
    same = bool(len(both)) and bool(np.allclose(both["pvalue"], both["padj"], rtol=1e-9, atol=0.0))
    flag("de_adjusted_pvalues", same, "padj equals pvalue for every tested gene" if same else "padj differs from pvalue")
elif a.audit == "all_defective":
    verdict = {i: "defective" for i in ITEMS}
for kv in a.audit_force:
    k, v = kv.split("=")
    verdict[k] = v
pd.DataFrame({"item": ITEMS, "verdict": [verdict[i] for i in ITEMS],
              "evidence": [evidence[i] for i in ITEMS]}).to_csv(out / "audit.tsv", sep="\t", index=False)

# ----------------------------------------------------------------------------- 5 erratum
if (prev / "mail_top_responders.tsv").exists():
    mailed = pd.read_csv(prev / "mail_top_responders.tsv", sep="\t", dtype=str)
    sig = dict(zip(de["gene_id"], de["significant"]))
    if a.erratum == "own":
        v = ["confirmed" if sig.get(g, 0) == 1 else "retracted" for g in mailed["gene_id"]]
    else:
        v = ["confirmed" if a.erratum == "all_confirmed" else "retracted"] * len(mailed)
    pd.DataFrame({"gene_id": mailed["gene_id"], "verdict": v}).to_csv(out / "erratum.tsv", sep="\t", index=False)
print(f"design={design} n={len(keep)} mixed={mixed} mismatch={mismatch} excluded={excluded} "
      f"significant={int(de['significant'].sum())} audit={verdict}")
