#!/usr/bin/env python3
"""Independent second solution: no pysam. Reads `samtools view` text, walks CIGARs itself, counts
further mismatches by comparing the read with the reference (not the NM tag), and decides indel
support by left-normalising each read's gap instead of comparing haplotype strings."""
import subprocess
from collections import defaultdict
from pathlib import Path

inp, out = Path("/task/input"), Path("/task/output")
out.mkdir(parents=True, exist_ok=True)

lines = (inp / "reference" / "ref.fa").read_text().splitlines()
contig = lines[0][1:].split()[0]
ref = "".join(lines[1:]).upper()
cands = []
for line in open(inp / "candidates.vcf"):
    if not line.startswith("#"):
        f = line.rstrip("\n").split("\t")
        cands.append((f[2], int(f[1]) - 1, f[3], f[4]))


def parse_cigar(c):
    out_, num = [], ""
    for ch in c:
        if ch.isdigit():
            num += ch
        else:
            out_.append((ch, int(num)))
            num = ""
    return out_


def left_normalise(gap, d):
    while gap > 0 and ref[gap - 1] == ref[gap + d - 1]:
        gap -= 1
    return gap


def reads_at(bam, pos0):
    """Non-duplicate primary reads overlapping pos0 (1-based region query), first per read name."""
    r = subprocess.run(["samtools", "view", "-F", "0xD04", str(bam), f"{contig}:{pos0 + 1}-{pos0 + 1}"],
                       capture_output=True, text=True, check=True)
    seen = set()
    for ln in r.stdout.splitlines():
        f = ln.split("\t")
        if f[0] in seen:
            continue
        seen.add(f[0])
        yield int(f[1]), int(f[3]) - 1, int(f[4]), parse_cigar(f[5]), f[9], f[10]


def review(bam, pos0, r, alt):
    is_del, d = len(r) > len(alt), len(r) - len(alt)
    sup, depth, strands = [], 0, defaultdict(int)
    for flag, start, mapq, cigar, seq, qual in reads_at(bam, pos0):
        depth += 1
        strand = "r" if flag & 16 else "f"
        strands[strand] += 1
        rp, qp, mism, base_at, q_at, gaps = start, 0, 0, None, None, []
        for op, ln in cigar:
            if op in "M=X":
                for k in range(ln):
                    if seq[qp + k] != ref[rp + k]:
                        mism += 1
                    if rp + k == pos0:
                        base_at, q_at = qp + k, ord(qual[qp + k]) - 33
                rp += ln
                qp += ln
            elif op == "D":
                gaps.append((rp, ln, qp))
                rp += ln
            elif op in "IS":
                qp += ln
            elif op == "N":
                rp += ln
        n = len(seq)
        if is_del:
            for g, ln, q in gaps:
                if ln == d and left_normalise(g, ln) == pos0 + 1:
                    sup.append((99, mapq, strand, min(q, n - q), mism))
                    break
        elif base_at is not None and seq[base_at] == alt:
            sup.append((q_at, mapq, strand, min(base_at, n - 1 - base_at), mism - 1))
    return sup, depth, strands


rows = ["variant_id\tcall\ttag\ttumor_alt\ttumor_depth\ttumor_vaf\tnormal_alt\tnormal_depth"]
for vid, pos0, r, alt in cands:
    ts, td, st = review(inp / "bam" / "tumor.bam", pos0, r, alt)
    ns, nd, _ = review(inp / "bam" / "normal.bam", pos0, r, alt)
    n = len(ts)
    if len(ns) >= 2 and nd and len(ns) / nd >= 0.05:
        tag = "GERMLINE"
    elif n < 4:
        tag = "LOW_SUPPORT"
    elif 4 * sum(1 for s in ts if s[0] < 20) >= 3 * n:
        tag = "LOW_BASE_QUALITY"
    elif 4 * sum(1 for s in ts if s[1] < 20) >= 3 * n:
        tag = "LOW_MAPPING_QUALITY"
    elif sum(s[4] for s in ts) >= 4 * n:
        tag = "MULTI_MISMATCH"
    elif 10 * max(sum(1 for s in ts if s[2] == "f"), sum(1 for s in ts if s[2] == "r")) >= 9 * n \
            and st["f"] > 0 and st["r"] > 0:
        tag = "STRAND_BIAS"
    elif 4 * sum(1 for s in ts if s[3] < 5) >= 3 * n:
        tag = "READ_END"
    else:
        tag = "PASS"
    vaf = n / td if td else 0.0
    rows.append(f"{vid}\t{'somatic' if tag == 'PASS' else 'fail'}\t{tag}\t{n}\t{td}\t{vaf:.4f}\t{len(ns)}\t{nd}")
(out / "review.tsv").write_text("\n".join(rows) + "\n")
(out / "igv_session.xml").write_text(
    '<?xml version="1.0" encoding="UTF-8"?>\n<Session genome="/task/input/reference/ref.fa" version="8">\n'
    '<Resources>\n<Resource path="/task/input/bam/normal.bam"/>\n<Resource path="/task/input/bam/tumor.bam"/>\n'
    '<Resource path="/task/input/candidates.vcf"/>\n</Resources>\n</Session>\n')
