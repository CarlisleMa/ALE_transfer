#!/usr/bin/env python3
"""Reference analysis of the 16S V4 amplicon batch:

  cutadapt (both read orientations) -> vsearch mergepairs -> expected-error filter -> pooled
  dereplication -> UNOISE3 -> uchime_denovo -> all merged reads mapped back to the variants ->
  genus by best global hit -> reagent contaminants by enrichment in the extraction blanks.

The non-default flags replay plausible-but-wrong analyses for the validator's attack suite; the
oracle (solve.sh) uses the defaults.
"""
import argparse
import subprocess
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ap = argparse.ArgumentParser()
ap.add_argument("--input", default="/task/input")
ap.add_argument("--output", default="/task/output")
ap.add_argument("--work", default="/tmp/amplicon_work")
ap.add_argument("--threads", default="2")
ap.add_argument("--primers", choices=["both", "fwd_keep", "fwd_discard", "none", "fixed"], default="both")
ap.add_argument("--cluster", choices=["unoise", "otu97"], default="unoise")
ap.add_argument("--chimera", choices=["uchime", "uchime3", "off"], default="uchime")
ap.add_argument("--contam", choices=["enrichment", "none", "any_in_neg"], default="enrichment")
ap.add_argument("--counts", choices=["all", "exact"], default="all")
a = ap.parse_args()

inp, out, work = Path(a.input), Path(a.output), Path(a.work)
out.mkdir(parents=True, exist_ok=True)
work.mkdir(parents=True, exist_ok=True)
T = a.threads

samples = pd.read_csv(inp / "samples.tsv", sep="\t", dtype=str)
sids = list(samples["sample_id"])
prim = pd.read_csv(inp / "primers.tsv", sep="\t", dtype=str)
FWD = prim.loc[prim["strand"] == "forward", "sequence_5to3"].iloc[0]
REV = prim.loc[prim["strand"] == "reverse", "sequence_5to3"].iloc[0]


def sh(cmd):
    r = subprocess.run([str(c) for c in cmd], capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"command failed: {' '.join(map(str, cmd))}\n{r.stderr[-2000:]}")
    return r


def has_reads(p):
    p = Path(p)
    if not p.exists() or p.stat().st_size == 0:
        return False
    if p.suffix == ".gz":
        import gzip
        with gzip.open(p, "rb") as fh:
            return bool(fh.read(1))
    return True


def merge(f, r, o):
    if not (has_reads(f) and has_reads(r)):
        Path(o).write_text("")
        return
    sh(["vsearch", "--fastq_mergepairs", f, "--reverse", r, "--fastqout", o, "--threads", T, "--quiet"])


def read_fastq_seqs(p):
    seqs = []
    with open(p) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                seqs.append(line.strip())
    return seqs


# ------------------------------------------------------------------ 1. primers, orientation, merging
merged = {}
for row in samples.itertuples():
    s = row.sample_id
    r1, r2 = inp / "fastq" / row.fastq_r1, inp / "fastq" / row.fastq_r2
    w = work / s
    w.mkdir(exist_ok=True)
    m = w / "merged.fq"
    if a.primers == "both":
        sh(["cutadapt", "-j", T, "-e", "0.1", "-O", "12", "-g", FWD, "-G", REV, "--pair-filter=any",
            "--untrimmed-output", w / "u1.fq", "--untrimmed-paired-output", w / "u2.fq",
            "-o", w / "a1.fq", "-p", w / "a2.fq", r1, r2])
        merge(w / "a1.fq", w / "a2.fq", w / "ma.fq")
        if has_reads(w / "u1.fq"):
            sh(["cutadapt", "-j", T, "-e", "0.1", "-O", "12", "-g", REV, "-G", FWD, "--pair-filter=any",
                "--discard-untrimmed", "-o", w / "b1.fq", "-p", w / "b2.fq", w / "u1.fq", w / "u2.fq"])
            merge(w / "b2.fq", w / "b1.fq", w / "mb.fq")        # mate 2 is the 515F side
        else:
            (w / "mb.fq").write_text("")
        m.write_text((w / "ma.fq").read_text() + (w / "mb.fq").read_text())
    elif a.primers in ("fwd_keep", "fwd_discard"):
        cmd = ["cutadapt", "-j", T, "-e", "0.1", "-O", "12", "-g", FWD, "-G", REV,
               "-o", w / "a1.fq", "-p", w / "a2.fq", r1, r2]
        if a.primers == "fwd_discard":
            cmd[1:1] = ["--discard-untrimmed"]
        sh(cmd)
        merge(w / "a1.fq", w / "a2.fq", m)
    else:
        merge(r1, r2, w / "raw_merged.fq")
        if a.primers == "fixed" and has_reads(w / "raw_merged.fq"):
            sh(["vsearch", "--fastq_filter", w / "raw_merged.fq", "--fastq_stripleft", str(len(FWD)),
                "--fastq_stripright", str(len(REV)), "--fastqout", m, "--quiet"])
        else:
            m.write_text((w / "raw_merged.fq").read_text())
    merged[s] = m

# ------------------------------------------------------------------ 2. quality filter, pooled denoising
with open(work / "pooled.fa", "w") as pooled:
    for s in sids:
        if not has_reads(merged[s]):
            continue
        sh(["vsearch", "--fastq_filter", merged[s], "--fastq_maxee", "1.0", "--fastq_minlen", "150",
            "--fastaout", work / s / "filtered.fa", "--relabel", f"{s}.", "--quiet"])
        pooled.write((work / s / "filtered.fa").read_text())
sh(["vsearch", "--derep_fulllength", work / "pooled.fa", "--sizeout", "--minuniquesize", "2",
    "--output", work / "derep.fa", "--relabel", "U", "--quiet"])
if a.cluster == "unoise":
    sh(["vsearch", "--cluster_unoise", work / "derep.fa", "--minsize", "8", "--unoise_alpha", "2.0",
        "--centroids", work / "centroids.fa", "--sizein", "--sizeout", "--threads", T, "--quiet"])
else:
    sh(["vsearch", "--cluster_size", work / "derep.fa", "--id", "0.97", "--centroids", work / "centroids.fa",
        "--sizein", "--sizeout", "--threads", T, "--quiet"])
if a.chimera != "off":
    sh(["vsearch", "--uchime3_denovo" if a.chimera == "uchime3" else "--uchime_denovo", work / "centroids.fa",
        "--nonchimeras", work / "nonchim.fa", "--sizein", "--sizeout", "--quiet"])
    final = work / "nonchim.fa"
else:
    final = work / "centroids.fa"

asv = []          # (size, sequence)
name, chunks = None, []
for line in list(open(final)) + [">end;size=0"]:
    if line.startswith(">"):
        if name is not None:
            asv.append((int(name.split("size=")[1].split(";")[0]), "".join(chunks).upper()))
        name, chunks = line.strip(), []
    else:
        chunks.append(line.strip())
asv.sort(key=lambda x: (-x[0], x[1]))
ids = [f"ASV{n + 1:04d}" for n in range(len(asv))]
seqs = [s for _, s in asv]
(out / "asv.fasta").write_text("".join(f">{i}\n{s}\n" for i, s in zip(ids, seqs)))

# ------------------------------------------------------------------ 3. count table: all merged reads
exact = {s: i for i, s in enumerate(seqs)}
table = np.zeros((len(seqs), len(sids)), dtype=int)
with open(work / "nonexact.fa", "w") as fh:
    n_q = 0
    for j, s in enumerate(sids):
        uniq = Counter(read_fastq_seqs(merged[s])) if has_reads(merged[s]) else {}
        for seq, c in uniq.items():
            if seq in exact:
                table[exact[seq], j] += c
            elif len(seq) >= 32:
                fh.write(f">{j}_{n_q};size={c}\n{seq}\n")
                n_q += 1
if n_q and seqs and a.counts == "all":
    sh(["vsearch", "--usearch_global", work / "nonexact.fa", "--db", out / "asv.fasta", "--id", "0.97",
        "--strand", "plus", "--maxaccepts", "4", "--maxrejects", "8", "--top_hits_only",
        "--userout", work / "map.tsv", "--userfields", "query+target+id", "--threads", T, "--quiet"])
    seen = set()
    idx = {i: n for n, i in enumerate(ids)}
    for line in open(work / "map.tsv"):
        q, t, _ = line.rstrip("\n").split("\t")
        if q in seen:
            continue
        seen.add(q)
        j = int(q.split("_")[0])
        c = int(q.split("size=")[1].rstrip(";"))
        table[idx[t], j] += c
tab = pd.DataFrame(table, columns=sids)
tab.insert(0, "asv_id", ids)
tab.to_csv(out / "asv_table.tsv", sep="\t", index=False)

# ------------------------------------------------------------------ 4. genus: best global hit >= 97 %
sh(["vsearch", "--usearch_global", out / "asv.fasta", "--db", inp / "reference_db.fasta", "--id", "0.97",
    "--strand", "plus", "--maxaccepts", "0", "--maxrejects", "0", "--top_hits_only",
    "--userout", work / "tax.tsv", "--userfields", "query+target+id", "--threads", T, "--quiet"])
ref_genus = dict(pd.read_csv(inp / "reference_taxonomy.tsv", sep="\t", dtype=str).values)
hits = {}
for line in open(work / "tax.tsv"):
    q, t, _ = line.rstrip("\n").split("\t")
    hits.setdefault(q, set()).add(ref_genus[t])
genus = [next(iter(hits[i])) if len(hits.get(i, ())) == 1 else "Unclassified" for i in ids]
pd.DataFrame({"asv_id": ids, "genus": genus}).to_csv(out / "taxonomy.tsv", sep="\t", index=False)

# ------------------------------------------------------------------ 5. reagent contaminants
neg = [j for j, t in enumerate(samples["sample_type"]) if t == "negative_control"]
real = [j for j in range(len(sids)) if j not in neg]
flag = np.zeros(len(seqs), dtype=bool)
if a.contam == "any_in_neg":
    flag = table[:, neg].sum(axis=1) > 0
elif a.contam == "enrichment" and neg:
    # null: reads of an ASV in the blanks are spill-over from the other libraries of the run, i.e.
    # they occur at the ASV's pooled relative abundance among the non-blank libraries
    n_neg = table[:, neg].sum(axis=1)
    big_n = max(int(n_neg.sum()), 1)
    share_real = table[:, real].sum(axis=1) / max(int(table[:, real].sum()), 1)
    p = stats.poisson.sf(n_neg - 1, big_n * np.maximum(share_real, 1e-9))
    in_all = (table[:, neg] > 0).all(axis=1)
    flag = in_all & (n_neg / big_n >= 2.0 * share_real) & (p < 1e-6)
(out / "contaminants.txt").write_text("asv_id\n" + "".join(f"{i}\n" for i, f in zip(ids, flag) if f))
print(f"ASVs={len(ids)} assigned_pairs={table.sum(axis=0).tolist()} contaminants={int(flag.sum())}")
