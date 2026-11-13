#!/usr/bin/env python3
"""Level-4 review logic on SAM text (standard library only).

The alignments of a level-4 task can only be read through IGV; this module works on what IGV
hands out as text: the SAM files written by the alignment track's "Export Alignments..." command.
It applies the review SOP (R-07 rev. 4) molecule by molecule. The non-default switches replay
careless reviews for the validator's attack suite.

Used by the reference solution (GUI route) inside the task image and by the generator on the
host (cross-check of the planted truth, careless variants for the battery).
"""
import argparse
import bisect
import re
from pathlib import Path

HEADER = "variant_id\tcall\ttag\ttumor_alt\ttumor_depth\ttumor_vaf\tnormal_alt\tnormal_depth"
_CIG = re.compile(r"(\d+)([MIDNSHP=X])")


class Read:
    __slots__ = ("name", "flag", "start", "end", "mapq", "cigar", "tlen", "seq", "qual")

    def __init__(self, f):
        self.name, self.flag, self.start, self.mapq = f[0], int(f[1]), int(f[3]) - 1, int(f[4])
        self.cigar = [(int(n), op) for n, op in _CIG.findall(f[5])]
        self.tlen, self.seq, self.qual = int(f[8]), f[9], f[10]
        self.end = self.start + sum(n for n, op in self.cigar if op in "MDN=X")


def load_sam(paths):
    """Reads of one sample from one or more SAM files (overlapping exports are merged)."""
    seen, reads = set(), []
    for p in paths:
        with open(p) as fh:
            for line in fh:
                if line.startswith("@") or not line.strip():
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 11 or f[5] == "*":
                    continue
                key = (f[0], int(f[1]) & 0xC0, f[3])
                if key in seen:
                    continue
                seen.add(key)
                reads.append(Read(f))
    reads.sort(key=lambda r: r.start)
    return reads


def load_reference(path):
    lines = Path(path).read_text().splitlines()
    return lines[0][1:].split()[0], "".join(lines[1:]).upper()


def load_candidates(path):
    out = []
    for line in open(path):
        if not line.startswith("#") and line.strip():
            f = line.rstrip("\n").split("\t")
            out.append((f[2], int(f[1]) - 1, f[3], f[4]))
    return out


def left_normalise(ref, gap, d):
    while gap > 0 and ref[gap - 1] == ref[gap + d - 1]:
        gap -= 1
    return gap


class Sample:
    def __init__(self, reads, ref, keep_duplicates=False, exact_indel=False):
        self.reads, self.ref = reads, ref
        self.starts = [r.start for r in reads]
        self.keep_duplicates, self.exact_indel = keep_duplicates, exact_indel
        self.span = max((r.end - r.start for r in reads), default=0)

    def evidence(self, pos0, r, alt):
        """Molecule-level evidence at one candidate.
        Returns ({molecule: attributes} for the supporting molecules, depth, covering molecules per strand)."""
        ref = self.ref
        is_del, d = len(r) > len(alt), len(r) - len(alt)
        lo = bisect.bisect_left(self.starts, pos0 - self.span)
        hi = bisect.bisect_right(self.starts, pos0)
        seen, support, strands = set(), {}, {"f": 0, "r": 0}
        for rd in self.reads[lo:hi]:
            if rd.flag & 0x904 or not (rd.start <= pos0 < rd.end):
                continue
            if rd.flag & 0x400 and not self.keep_duplicates:
                continue
            if rd.name in seen:                     # a molecule is counted once
                continue
            seen.add(rd.name)
            strand = "r" if rd.flag & 16 else "f"
            strands[strand] += 1
            rp, qp, mism, q_at, gaps = rd.start, 0, 0, None, []
            for n, op in rd.cigar:
                if op in "M=X":
                    seg = rd.seq[qp:qp + n]
                    rseg = ref[rp:rp + n]
                    mism += sum(1 for a, b in zip(seg, rseg) if a != b)
                    if rp <= pos0 < rp + n:
                        q_at = qp + pos0 - rp
                    rp += n
                    qp += n
                elif op == "D":
                    gaps.append((rp, n, qp))
                    rp += n
                elif op == "N":
                    rp += n
                elif op in "IS":
                    qp += n
            n_seq = len(rd.seq)
            if is_del:
                for g, n, q in gaps:
                    same = (g == pos0 + 1) if self.exact_indel else (left_normalise(ref, g, n) == pos0 + 1)
                    if n == d and same:
                        support[rd.name] = {"bq": 99, "mapq": rd.mapq, "strand": strand, "end": min(q, n_seq - q),
                                            "other_mm": mism, "isize": abs(rd.tlen)}
                        break
            elif q_at is not None and rd.seq[q_at] == alt:
                support[rd.name] = {"bq": ord(rd.qual[q_at]) - 33, "mapq": rd.mapq, "strand": strand,
                                    "end": min(q_at, n_seq - 1 - q_at), "other_mm": mism - 1,
                                    "isize": abs(rd.tlen)}
        return support, len(seen), strands


def fractions(t_sup, t_strands):
    """The quantities the SOP thresholds are applied to (also used by the generator's margin check)."""
    s = list(t_sup.values())
    n = len(s)
    if not n:
        return {"n": 0}
    f = sum(1 for x in s if x["strand"] == "f")
    return {"n": n,
            "low_bq": sum(1 for x in s if x["bq"] < 20) / n,
            "low_mapq": sum(1 for x in s if x["mapq"] < 20) / n,
            "other_mm": sum(x["other_mm"] for x in s) / n,
            "one_strand": max(f, n - f) / n,
            "both_strands_cover": t_strands["f"] > 0 and t_strands["r"] > 0,
            "read_end": sum(1 for x in s if x["end"] < 5) / n,
            "long_insert": sum(1 for x in s if x["isize"] > 1000) / n}


def tag_for(t_sup, t_strands, n_sup, n_depth, opt):
    n = len(t_sup)
    if opt.get("normal_any_alt"):
        if len(n_sup) >= 1:
            return "GERMLINE"
    elif len(n_sup) >= 2 and n_depth and len(n_sup) / n_depth >= 0.05:
        return "GERMLINE"
    if n < 4:
        return "LOW_SUPPORT"
    if opt.get("counts_only"):
        return "PASS"
    q = fractions(t_sup, t_strands)
    if q["low_bq"] >= 0.75:
        return "LOW_BASE_QUALITY"
    if q["low_mapq"] >= 0.75:
        return "LOW_MAPPING_QUALITY"
    if q["other_mm"] >= 4:
        return "MULTI_MISMATCH"
    if q["one_strand"] >= 0.9 and q["both_strands_cover"]:
        return "STRAND_BIAS"
    if q["read_end"] >= 0.75:
        return "READ_END"
    if q["long_insert"] >= 0.75 and not opt.get("no_pair_rule"):
        return "DISCORDANT_PAIRS"
    return "PASS"


def review(cands, ref, tumor, normal, opt=None):
    """cands: [(id, pos0, ref, alt)], tumor/normal: Sample. Returns (rows, diagnostics)."""
    opt = opt or {}
    res = []
    for vid, pos0, r, alt in cands:
        ts, td, tst = tumor.evidence(pos0, r, alt)
        ns, nd, _ = normal.evidence(pos0, r, alt)
        res.append({"id": vid, "pos0": pos0, "snv": len(r) == 1 and len(alt) == 1, "sup": ts, "depth": td,
                    "n_alt": len(ns), "n_depth": nd, "tag": tag_for(ts, tst, ns, nd, opt),
                    "frac": fractions(ts, tst), "shared": None})
    # multi-nucleotide events: adjacent candidates that both pass and sit on the same molecules
    for a, b in zip(res, res[1:]):
        if not (a["snv"] and b["snv"] and b["pos0"] == a["pos0"] + 1):
            continue
        union = set(a["sup"]) | set(b["sup"])
        shared = len(set(a["sup"]) & set(b["sup"])) / len(union) if union else 0.0
        a["shared"] = b["shared"] = shared
        if a["tag"] != "PASS" or b["tag"] != "PASS" or opt.get("no_mnv_rule"):
            continue
        if shared >= 0.9 or opt.get("merge_adjacent"):
            b["tag"] = "MNV_PARTNER"
    rows = [HEADER]
    for x in res:
        n, d = len(x["sup"]), x["depth"]
        rows.append(f"{x['id']}\t{'somatic' if x['tag'] == 'PASS' else 'fail'}\t{x['tag']}\t{n}\t{d}\t"
                    f"{(n / d if d else 0.0):.4f}\t{x['n_alt']}\t{x['n_depth']}")
    return rows, res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default="/task/input")
    ap.add_argument("--tumor", nargs="+", required=True, help="SAM file(s) exported from the tumour track")
    ap.add_argument("--normal", nargs="+", required=True, help="SAM file(s) exported from the normal track")
    ap.add_argument("--out", default="/task/output/review.tsv")
    for flag in ("counts-only", "keep-duplicates", "exact-indel", "no-mnv-rule", "merge-adjacent", "no-pair-rule",
                 "normal-any-alt"):
        ap.add_argument("--" + flag, action="store_true")
    a = ap.parse_args()
    inp = Path(a.input)
    _, ref = load_reference(inp / "reference" / "ref.fa")
    cands = load_candidates(inp / "candidates.vcf")
    t = Sample(load_sam(a.tumor), ref, a.keep_duplicates, a.exact_indel)
    n = Sample(load_sam(a.normal), ref, a.keep_duplicates, a.exact_indel)
    rows, res = review(cands, ref, t, n, vars(a))
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text("\n".join(rows) + "\n")
    print(f"reviewed {len(cands)} candidates: {sum(1 for x in res if x['tag'] == 'PASS')} somatic")


if __name__ == "__main__":
    main()
