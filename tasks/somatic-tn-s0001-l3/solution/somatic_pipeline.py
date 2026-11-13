#!/usr/bin/env python3
"""Reference solution: bwa mem -> fixmate/sort/markdup -> fragment-aware allele counting with pysam
-> pairing by germline genotype concordance -> tumour-vs-matched-normal calling with mapping-quality
and read-pair-orientation filters.

The switches reproduce what a less careful analysis would deliver (used by the validation battery):
  --pairing sheet     trust the sample sheet instead of verifying the pairs
  --min-mapq 0        count reads of any mapping quality
  --no-orientation    skip the read-pair-orientation check of the supporting reads
  --no-normal         tumour-only calling
  --strand-instead-of-orientation   require support on both alignment strands (classic strand-bias
                      filter) instead of in both read-pair orientations
  --spectrum-filter X in tumours that core_notes.md lists as FFPE, drop every C>T, G>A, G>T and C>A call
                      with an allele fraction below X (blanket damage filter; use with --no-orientation)
"""
import argparse
import csv
import gzip
import itertools
import json
import os
import subprocess
from collections import defaultdict
from pathlib import Path

import numpy as np
import pysam

IN = Path("/task/input")
OUT = Path("/task/output")
WORK = Path(os.environ.get("WORK", "/tmp/somatic_work"))
THREADS = 2
MIN_BQ = 20
END_MARGIN = 5          # bases at either end of an alignment that are not used as SNV evidence
DAMAGE_CLASSES = {("C", "T"), ("G", "A"), ("G", "T"), ("C", "A")}

LUT = np.full(256, 4, dtype=np.int64)
for _i, _b in enumerate(b"ACGT"):
    LUT[_b] = _i


def sh(cmd: str) -> None:
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")


def normalise(seq, pos0, ref, alt):
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos0 = ref[1:], alt[1:], pos0 + 1
    if len(ref) == len(alt):
        return pos0, ref, alt
    while pos0 > 0 and ref[-1] == alt[-1]:
        pos0 -= 1
        b = seq[pos0]
        ref, alt = b + ref[:-1], b + alt[:-1]
    return pos0, ref, alt


def align(sid: str, fq1: str, fq2: str) -> Path:
    bam = WORK / f"{sid}.bam"
    sh(f"set -o pipefail; bwa mem -t {THREADS} -R '@RG\\tID:{sid}\\tSM:{sid}\\tPL:ILLUMINA' {WORK}/ref.fa "
       f"{IN}/{fq1} {IN}/{fq2} 2> {WORK}/{sid}.bwa.log | samtools fixmate -m -u - - "
       f"| samtools sort -u -T {WORK}/{sid}.tmp - | samtools markdup - {bam} 2> {WORK}/{sid}.markdup.log "
       f"&& samtools index {bam}")
    return bam


def count(bam_path: Path, seq: str, min_mapq: int, by_strand: bool = False):
    """Per read-pair orientation (0 = F1R2, 1 = F2R1): base counts [2, 4, L], indel counts, spanning reads [L].
    by_strand=True splits by alignment strand instead (the plain strand-bias view)."""
    L = len(seq)
    cnt = np.zeros((2, 4, L), dtype=np.int32)
    indel = [defaultdict(int), defaultdict(int)]
    diff = np.zeros(L + 1, dtype=np.int32)
    with pysam.AlignmentFile(str(bam_path)) as bam:
        for r in bam.fetch(until_eof=True):
            if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate or r.mapping_quality < min_mapq:
                continue
            o = int(r.is_reverse) if by_strand else (0 if (r.is_read1 != r.is_reverse) else 1)
            start = r.reference_start
            lo = start
            # count every fragment once: where mates overlap, the reverse read yields to the forward read
            if (r.is_reverse and r.is_paired and not r.mate_is_unmapped and not r.mate_is_reverse
                    and r.next_reference_id == r.reference_id and r.next_reference_start <= start):
                lo = max(lo, r.next_reference_start + r.query_length)
            end = r.reference_end
            if lo >= end:
                continue
            diff[lo] += 1
            diff[end] -= 1
            codes = LUT[np.frombuffer(r.query_sequence.encode(), dtype=np.uint8)]
            quals = np.asarray(r.query_qualities)
            # mismatches in the last bases of an alignment are often an indel or clip drawn as substitutions
            q_lo, q_hi = r.query_alignment_start + END_MARGIN, r.query_alignment_end - END_MARGIN
            qpos, rpos = 0, start
            for op, ln in r.cigartuples:
                if op in (0, 7, 8):
                    a = max(rpos, lo)
                    if a < rpos + ln:
                        off = a - rpos
                        c = codes[qpos + off:qpos + ln]
                        qi = np.arange(qpos + off, qpos + ln)
                        ok = (quals[qpos + off:qpos + ln] >= MIN_BQ) & (c < 4) & (qi >= q_lo) & (qi < q_hi)
                        cnt[o][c[ok], np.arange(a, rpos + ln)[ok]] += 1
                    qpos += ln
                    rpos += ln
                elif op == 2:
                    if rpos - 1 >= lo:
                        indel[o][normalise(seq, rpos - 1, seq[rpos - 1:rpos + ln], seq[rpos - 1])] += 1
                    rpos += ln
                elif op == 1:
                    if rpos - 1 >= lo:
                        indel[o][normalise(seq, rpos - 1, seq[rpos - 1],
                                           seq[rpos - 1] + r.query_sequence[qpos:qpos + ln])] += 1
                    qpos += ln
                elif op == 4:
                    qpos += ln
                elif op == 3:
                    rpos += ln
    return {"cnt": cnt, "indel": indel, "span": np.cumsum(diff[:-1])}


def genotype_concordance(tum, nor, seq):
    """Fraction of sites variable among the normals at which tumour and normal show the same genotype class."""
    refc = LUT[np.frombuffer(seq.encode(), dtype=np.uint8)]
    L = len(seq)
    tot = {s: c["cnt"].sum(axis=0) for s, c in {**tum, **nor}.items()}      # [4, L]
    dp = {s: t.sum(axis=0) for s, t in tot.items()}
    sites = set()
    for s in nor:
        t = tot[s].astype(float).copy()
        t[refc, np.arange(L)] = 0
        alt = t.argmax(axis=0)
        frac = t.max(axis=0) / np.maximum(dp[s], 1)
        for x in np.nonzero((frac >= 0.25) & (dp[s] >= 20))[0]:
            sites.add((int(x), int(alt[x])))

    def cls(s, x, b):
        if dp[s][x] < 15:
            return None
        f = tot[s][b, x] / dp[s][x]
        return 0 if f < 0.12 else (2 if f > 0.88 else 1)

    conc = {}
    for t in tum:
        for n in nor:
            same = n_ok = 0
            for x, b in sites:
                a, c = cls(t, x, b), cls(n, x, b)
                if a is None or c is None:
                    continue
                n_ok += 1
                same += a == c
            conc[(t, n)] = (same / n_ok if n_ok else 0.0, n_ok)
    return conc


def call(t, n, seq, targets, a):
    """Somatic candidates of tumour counts `t` against normal counts `n` (None = tumour-only)."""
    L = len(seq)
    refc = LUT[np.frombuffer(seq.encode(), dtype=np.uint8)]
    inside = np.zeros(L, dtype=bool)
    for s, e in targets:
        inside[s:e] = True
    out = []

    def keep(t_alt, t_dp, o_alt, n_alt, n_dp):
        if t_alt < a.min_alt or t_dp < 20 or t_alt / t_dp < a.min_af:
            return False
        if n is not None and (n_dp < 10 or n_alt > max(1, 0.02 * n_dp)):
            return False
        if not a.no_orientation and (min(o_alt) < 2 or min(o_alt) / t_alt < 0.10):
            return False
        return True

    t_dp = t["cnt"].sum(axis=(0, 1))
    n_dp = n["cnt"].sum(axis=(0, 1)) if n is not None else None
    for b in range(4):
        t_alt = t["cnt"][:, b, :].sum(axis=0)
        cand = np.nonzero(inside & (refc != b) & (t_alt >= a.min_alt))[0]
        for x in cand:
            x = int(x)
            o_alt = (int(t["cnt"][0, b, x]), int(t["cnt"][1, b, x]))
            na = int(n["cnt"][:, b, x].sum()) if n is not None else 0
            nd = int(n_dp[x]) if n is not None else 0
            if keep(int(t_alt[x]), int(t_dp[x]), o_alt, na, nd):
                out.append((x, seq[x], "ACGT"[b], int(t_alt[x]) / int(t_dp[x]), int(t_alt[x]), int(t_dp[x]), o_alt, na, nd))
    for key in set(t["indel"][0]) | set(t["indel"][1]):
        x = key[0]
        if not inside[x]:
            continue
        o_alt = (t["indel"][0].get(key, 0), t["indel"][1].get(key, 0))
        ta, td = sum(o_alt), int(t["span"][x + 1])
        na = (n["indel"][0].get(key, 0) + n["indel"][1].get(key, 0)) if n is not None else 0
        nd = int(n["span"][x + 1]) if n is not None else 0
        if td and keep(ta, td, o_alt, na, nd):
            out.append((x, key[1], key[2], min(1.0, ta / td), ta, td, o_alt, na, nd))
    return sorted(out)


def write_vcf(path: Path, contig: str, L: int, tumour: str, normal: str, calls) -> None:
    lines = ["##fileformat=VCFv4.2", f"##contig=<ID={contig},length={L}>",
             '##INFO=<ID=AF,Number=1,Type=Float,Description="Fraction of tumour fragments carrying ALT">',
             '##INFO=<ID=TAD,Number=1,Type=Integer,Description="Tumour ALT fragments">',
             '##INFO=<ID=TDP,Number=1,Type=Integer,Description="Tumour fragments covering the site">',
             '##INFO=<ID=F1R2,Number=1,Type=Integer,Description="Tumour ALT fragments in F1R2 orientation">',
             '##INFO=<ID=F2R1,Number=1,Type=Integer,Description="Tumour ALT fragments in F2R1 orientation">',
             '##INFO=<ID=NAD,Number=1,Type=Integer,Description="Normal ALT fragments">',
             '##INFO=<ID=NDP,Number=1,Type=Integer,Description="Normal fragments covering the site">',
             f"##tumour_sample={tumour}", f"##normal_sample={normal}",
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    for x, ref, alt, af, ta, td, o_alt, na, nd in calls:
        lines.append(f"{contig}\t{x + 1}\t.\t{ref}\t{alt}\t.\tPASS\tAF={af:.4f};TAD={ta};TDP={td};"
                     f"F1R2={o_alt[0]};F2R1={o_alt[1]};NAD={na};NDP={nd}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as fh:
        fh.write("\n".join(lines) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairing", choices=["verify", "sheet"], default="verify")
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--no-orientation", action="store_true")
    ap.add_argument("--no-normal", action="store_true")
    ap.add_argument("--strand-instead-of-orientation", action="store_true")
    ap.add_argument("--spectrum-filter", type=float, default=0.0)
    ap.add_argument("--min-alt", type=int, default=5)
    ap.add_argument("--min-af", type=float, default=0.04)
    a = ap.parse_args()

    WORK.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    sh(f"cp {IN}/reference/ref.fa {WORK}/ref.fa && samtools faidx {WORK}/ref.fa && bwa index {WORK}/ref.fa 2> {WORK}/index.log")
    fa = pysam.FastaFile(str(WORK / "ref.fa"))
    contig = fa.references[0]
    seq = fa.fetch(contig).upper()
    targets = [(int(f[1]), int(f[2])) for f in csv.reader(open(IN / "targets.bed"), delimiter="\t") if len(f) >= 3]
    sheet = list(csv.DictReader(open(IN / "manifest.tsv"), delimiter="\t"))

    counts, tum, sheet_nor = {}, {}, {}
    for row in sheet:
        sid = row["sample_id"]
        bam = align(sid, row["fastq_1"], row["fastq_2"])
        counts[sid] = count(bam, seq, a.min_mapq, a.strand_instead_of_orientation)
        (tum if row["sample_type"] == "tumour" else sheet_nor)[row["patient"]] = sid
    patients = sorted(tum)
    normals = sorted(sheet_nor.values())

    conc = genotype_concordance({tum[p]: counts[tum[p]] for p in patients}, {n: counts[n] for n in normals}, seq)
    best = max(itertools.permutations(normals),
               key=lambda perm: sum(conc[(tum[p], n)][0] for p, n in zip(patients, perm)))
    verified = dict(zip(patients, best))
    pairing = verified if a.pairing == "verify" else dict(sheet_nor)
    qc = {"concordance": {f"{t}|{n}": [round(v[0], 4), v[1]] for (t, n), v in conc.items()},
          "sheet": sheet_nor, "verified": verified}
    (WORK / "pairing_qc.json").write_text(json.dumps(qc, indent=1))

    ffpe = set()
    for line in open(IN / "core_notes.md"):
        cells = [c.strip() for c in line.split("|")]
        if len(cells) > 3 and "FFPE" in cells[2]:
            ffpe.add(cells[1])
    summary = {}
    with open(OUT / "pairing.tsv", "w") as fh:
        fh.write("patient\ttumour_sample\tnormal_sample\n")
        for p in patients:
            fh.write(f"{p}\t{tum[p]}\t{pairing[p]}\n")
    for p in patients:
        calls = call(counts[tum[p]], None if a.no_normal else counts[pairing[p]], seq, targets, a)
        if a.spectrum_filter and tum[p] in ffpe:
            calls = [c for c in calls if not ((c[1], c[2]) in DAMAGE_CLASSES and c[3] < a.spectrum_filter)]
        write_vcf(OUT / "somatic" / f"{p}.vcf.gz", contig, len(seq), tum[p], pairing[p], calls)
        summary[p] = {"tumour_sample": tum[p], "normal_sample": pairing[p],
                      "n_snv": sum(1 for c in calls if len(c[1]) == 1 and len(c[2]) == 1),
                      "n_indel": sum(1 for c in calls if len(c[1]) != len(c[2]))}
    (OUT / "summary.json").write_text(json.dumps({"patients": summary}, indent=1) + "\n")


if __name__ == "__main__":
    main()
