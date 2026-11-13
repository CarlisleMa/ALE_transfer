#!/usr/bin/env python3
"""Reference review (pysam). Evidence is collected per molecule, indel support is decided by
comparing the haplotype a read's gap produces with the haplotype of the candidate, and the SOP
tags are applied in their order of precedence. The non-default flags replay careless reviews for
the validator's attack suite."""
import argparse
from pathlib import Path

import pysam

ap = argparse.ArgumentParser()
ap.add_argument("--input", default="/task/input")
ap.add_argument("--output", default="/task/output")
ap.add_argument("--counts-only", action="store_true")
ap.add_argument("--keep-duplicates", action="store_true")
ap.add_argument("--exact-indel", action="store_true")
a = ap.parse_args()
inp, out = Path(a.input), Path(a.output)
out.mkdir(parents=True, exist_ok=True)

fa = pysam.FastaFile(str(inp / "reference" / "ref.fa"))
contig = fa.references[0]
ref = fa.fetch(contig)
cands = []
for line in open(inp / "candidates.vcf"):
    if line.startswith("#"):
        continue
    c, pos, vid, r, alt = line.split("\t")[:5]
    cands.append((vid, int(pos) - 1, r, alt))


def evidence(bam, pos0, r, alt):
    """Per-molecule evidence at one candidate: list of dicts for supporting molecules, depth, strands."""
    is_del = len(r) > len(alt)
    d = len(r) - len(alt)
    hap_alt = ref[:pos0 + 1] + ref[pos0 + 1 + d:] if is_del else None
    seen, support = set(), []
    strands = {"f": 0, "r": 0}
    for read in bam.fetch(contig, pos0, pos0 + 1):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.is_duplicate and not a.keep_duplicates:
            continue
        if read.query_name in seen:      # a molecule is counted once
            continue
        seen.add(read.query_name)
        strand = "r" if read.is_reverse else "f"
        strands[strand] += 1
        nm = read.get_tag("NM")
        if is_del:
            rp, qp, hit = read.reference_start, 0, None
            for op, ln in read.cigartuples:
                if op == 2:
                    if a.exact_indel:
                        ok = ln == d and rp == pos0 + 1
                    else:
                        ok = ln == d and ref[:rp] + ref[rp + ln:] == hap_alt
                    if ok:
                        hit = qp
                    rp += ln
                elif op in (0, 7, 8):
                    rp += ln
                    qp += ln
                elif op in (1, 4):
                    qp += ln
                elif op == 3:
                    rp += ln
            if hit is not None:
                support.append({"bq": 40, "mapq": read.mapping_quality, "strand": strand,
                                "end": min(hit, read.query_length - hit), "other_mm": nm - d})
        else:
            qpos = None
            for q, rpos in read.get_aligned_pairs(matches_only=True):
                if rpos == pos0:
                    qpos = q
                    break
            if qpos is None or read.query_sequence[qpos] != alt:
                continue
            support.append({"bq": read.query_qualities[qpos], "mapq": read.mapping_quality, "strand": strand,
                            "end": min(qpos, read.query_length - 1 - qpos), "other_mm": nm - 1})
    return support, len(seen), strands


def tag_for(t_sup, t_depth, t_strands, n_sup, n_depth):
    n = len(t_sup)
    if len(n_sup) >= 2 and n_depth and len(n_sup) / n_depth >= 0.05:
        return "GERMLINE"
    if n < 4:
        return "LOW_SUPPORT"
    if a.counts_only:
        return "PASS"
    if sum(s["bq"] < 20 for s in t_sup) >= 0.75 * n:
        return "LOW_BASE_QUALITY"
    if sum(s["mapq"] < 20 for s in t_sup) >= 0.75 * n:
        return "LOW_MAPPING_QUALITY"
    if sum(s["other_mm"] for s in t_sup) / n >= 4:
        return "MULTI_MISMATCH"
    f = sum(s["strand"] == "f" for s in t_sup)
    if max(f, n - f) >= 0.9 * n and t_strands["f"] > 0 and t_strands["r"] > 0:
        return "STRAND_BIAS"
    if sum(s["end"] < 5 for s in t_sup) >= 0.75 * n:
        return "READ_END"
    return "PASS"


tb = pysam.AlignmentFile(str(inp / "bam" / "tumor.bam"))
nb = pysam.AlignmentFile(str(inp / "bam" / "normal.bam"))
rows = ["variant_id\tcall\ttag\ttumor_alt\ttumor_depth\ttumor_vaf\tnormal_alt\tnormal_depth"]
for vid, pos0, r, alt in cands:
    ts, td, tst = evidence(tb, pos0, r, alt)
    ns, nd, _ = evidence(nb, pos0, r, alt)
    tag = tag_for(ts, td, tst, ns, nd)
    call = "somatic" if tag == "PASS" else "fail"
    vaf = len(ts) / td if td else 0.0
    rows.append(f"{vid}\t{call}\t{tag}\t{len(ts)}\t{td}\t{vaf:.4f}\t{len(ns)}\t{nd}")
(out / "review.tsv").write_text("\n".join(rows) + "\n")

first = cands[0]
(out / "igv_session.xml").write_text(
    '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'
    f'<Session genome="/task/input/reference/ref.fa" locus="{contig}:{max(1, first[1] - 100)}-{first[1] + 100}" version="8">\n'
    ' <Resources>\n'
    '  <Resource path="/task/input/bam/tumor.bam"/>\n'
    '  <Resource path="/task/input/bam/normal.bam"/>\n'
    '  <Resource path="/task/input/candidates.bed"/>\n'
    ' </Resources>\n'
    '</Session>\n')
print(f"reviewed {len(cands)} candidates: {sum(1 for r in rows[1:] if r.split(chr(9))[1] == 'somatic')} somatic")
