"""Family: manual review (refinement) of somatic variant candidates in tumour/normal BAMs.

The classic GUI workflow of cancer genomics: a caller proposes candidates, an analyst opens IGV and
decides call by call whether the evidence is real, following the lab's review SOP. Here every
candidate's mechanism is planted while the reads are synthesised (BAMs are written directly, so the
generator knows exactly which read carries what), which makes the review outcome a fact rather
than an opinion.

GUI-capable: the image ships IGV plus a desktop MCP server (screenshot/click/type tools); the task
is graded on artefacts, so it can equally be solved from the terminal (ALE observed that agents do
exactly that).

Decision points (level 2 adds the artefact classes, level 3 the last two):
  artefact classes        strand, read end, base quality, mapping quality, multi-mismatch reads
  duplicate_only_support  all but one supporting read are flagged duplicates of one molecule
  indel_representation    reads draw the gap at the right end of a repeat, the candidate list is
                          left-aligned: position-exact evidence counting finds nothing

Level 4 (see `level4.py`): the alignments are shipped as CRAM 3.1 on a locked-down workstation
image in which IGV is the only reader (T), and the review needs per-molecule logic that allele
counts cannot give (P): multi-nucleotide events, support confined to discordant pairs, traces of
the variant in the normal, mixed artefacts.
"""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

import numpy as np

from ... import desktop, dockerlib
from ...spec import Attack, Probe, TaskBuild, Trap

ASSETS = Path(__file__).parent / "assets"
CACHE = Path(__file__).resolve().parents[4] / ".cache" / "variant_review"
REF_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz"
REF_SHA256 = "05f9d97d6fbfd08a44ca45b50837ca2ae9c471f35ba79dffec04d2cb5eaaf695"
L = 60_000
RL = 100
CONTIG = "chr22_roi"
FAIL_TAGS = ["GERMLINE", "LOW_SUPPORT", "LOW_BASE_QUALITY", "LOW_MAPPING_QUALITY", "MULTI_MISMATCH",
             "STRAND_BIAS", "READ_END"]


def _reference(seed: int) -> str:
    fa = CACHE / "chr22.hg38.fa.gz"
    if not fa.exists():
        import subprocess
        CACHE.mkdir(parents=True, exist_ok=True)
        subprocess.run(["curl", "-fsSL", "-o", str(fa), REF_URL], check=True)
    if hashlib.sha256(fa.read_bytes()).hexdigest() != REF_SHA256:
        raise RuntimeError("reference checksum mismatch")
    with gzip.open(fa, "rt") as fh:
        next(fh)
        seq = "".join(line.strip() for line in fh).upper()
    rng = np.random.default_rng(50_000 + seed)
    for _ in range(2000):
        s = int(rng.integers(18_000_000, 48_000_000))
        w = seq[s:s + L]
        gc = (w.count("G") + w.count("C")) / L
        if "N" not in w and 0.38 <= gc <= 0.58:
            return w
    raise RuntimeError("no suitable window")


def _homopolymers(ref: str, lo: int = 7, hi: int = 12) -> list[tuple[int, int]]:
    """(start, length) of homopolymer runs with lo <= length <= hi."""
    out, i = [], 0
    while i < len(ref):
        j = i
        while j < len(ref) and ref[j] == ref[i]:
            j += 1
        if lo <= j - i <= hi:
            out.append((i, j - i))
        i = j
    return out


class Frag:
    __slots__ = ("start", "flen", "hap", "r1_fwd", "dup", "mapq", "subs", "dele", "xmm", "name", "proper")

    def __init__(self, start, flen, hap, r1_fwd, name):
        self.start, self.flen, self.hap, self.r1_fwd, self.name = start, flen, hap, r1_fwd, name
        self.dup, self.mapq = False, 60
        self.proper = True  # level 4: pairs with an aberrant insert size are not flagged "proper pair"
        self.subs = {}      # ref pos -> (base, qual)
        self.dele = None    # (gap_start, length) drawn by the "aligner"
        self.xmm = 0        # extra mismatches per read (paralog-like reads)

    def fwd_span(self):
        return self.start, self.start + RL

    def rev_span(self):
        return self.start + self.flen - RL, self.start + self.flen

    def covering(self, x, margin=0):
        a, b = self.fwd_span()
        if a + margin <= x < b - margin:
            return "f"
        a, b = self.rev_span()
        if a + margin <= x < b - margin:
            return "r"
        return None


def _simulate(ref: str, cands: list[dict], depth: int, tumour: bool, rng, tag: str):
    """Return a list of Frag for one sample with all candidate mechanisms applied."""
    n = int(L * depth / (2 * RL))
    frags = []
    for i in range(n):
        flen = int(np.clip(rng.normal(340, 30), 232, 480))
        start = int(rng.integers(0, L - flen))
        frags.append(Frag(start, flen, int(rng.integers(0, 2)), bool(rng.integers(0, 2)), f"{tag}:{i}"))
    extra_id = [0]

    def extra(start, flen, r1_fwd=None):
        extra_id[0] += 1
        f = Frag(int(start), int(flen), int(rng.integers(0, 2)),
                 bool(rng.integers(0, 2)) if r1_fwd is None else r1_fwd, f"{tag}:x{extra_id[0]}")
        frags.append(f)
        return f

    def place(x, which, off):
        """Fragment whose `which` read ('f'/'r') has ref position x at read offset `off`."""
        flen = int(np.clip(rng.normal(340, 30), 232, 480))
        if which == "f":
            start = x - off
        else:
            start = x - off - (flen - RL)
        start = min(max(start, 0), L - flen)
        return extra(start, flen)

    for c in cands:
        x, alt, cls = c["pos0"], c["alt"], c["cls"]
        cov = [(f, f.covering(x)) for f in frags]
        cov = [(f, w) for f, w in cov if w]
        if cls in ("somatic_snv",):
            if tumour:
                for f, _ in cov:
                    if rng.random() < c["vaf"]:
                        f.subs[x] = (alt, 37)
            elif c.get("normal_trace"):     # level 4: one stray molecule in the normal, far below the SOP rule
                f, _ = cov[int(rng.integers(0, len(cov)))]
                f.subs[x] = (alt, 37)
        elif cls == "mnv":                  # level 4: two adjacent substitutions, planted from the left one
            if tumour and c["role"] == "L":
                x2, alt2 = x + 1, c["mate_alt"]
                both = [(f, w) for f, w in cov if f.covering(x2) == w]
                if c["phase"] == "cis":     # the same molecules carry both
                    for f, _ in both:
                        if rng.random() < c["vaf"]:
                            f.subs[x] = (alt, 37)
                            f.subs[x2] = (alt2, 37)
                else:                       # two disjoint sets of molecules, same size, same strand split:
                    fw = [f for f, w in both if w == "f"]       # allele counts cannot tell this from cis
                    rv = [f for f, w in both if w == "r"]
                    fw = [fw[i] for i in rng.permutation(len(fw))]
                    rv = [rv[i] for i in rng.permutation(len(rv))]
                    kf, kr = int(round(c["vaf"] * len(fw))), int(round(c["vaf"] * len(rv)))
                    for f in fw[:kf] + rv[:kr]:
                        f.subs[x] = (alt, 37)
                    for f in fw[kf:2 * kf] + rv[kr:2 * kr]:
                        f.subs[x2] = (alt2, 37)
        elif cls == "long_insert":          # level 4: support only from pairs with an aberrant insert size
            if tumour:
                for k in range(c["n_alt"]):
                    w = "f" if k % 2 == 0 else "r"
                    off = int(rng.integers(15, 85))
                    if k == c["n_alt"] - 1:     # one ordinary pair among them (mixed, still far from the rule)
                        f = place(x, w, off)
                    else:
                        flen = int(rng.integers(1200, 1900))
                        f = extra(x - off if w == "f" else x - off - (flen - RL), flen)
                        f.proper = False
                    f.subs[x] = (alt, 37)
        elif cls == "germline":
            for f, _ in cov:
                if f.hap == 1:
                    f.subs[x] = (alt, 37)
        elif cls == "low_support":
            if tumour:
                for f, _ in [cov[i] for i in rng.permutation(len(cov))[:c["n_alt"]]]:
                    f.subs[x] = (alt, 37)
        elif cls == "dup_only":
            if tumour:
                # 2-3 real molecules on both strands, each amplified into several flagged duplicates:
                # read-level counting sees a well-supported, balanced variant; molecule counting does not
                for k in range(c["n_orig"]):
                    f0 = place(x, "f" if k % 2 == 0 else "r", int(rng.integers(25, 75)))
                    f0.subs[x] = (alt, 37)
                    for _ in range(c["n_dup"]):
                        d = extra(f0.start, f0.flen, f0.r1_fwd)
                        d.subs, d.dup, d.hap = dict(f0.subs), True, f0.hap
        elif cls == "strand_bias":
            if tumour:
                w = c["strand"]
                for _ in range(c["n_alt"]):
                    f = place(x, w, int(rng.integers(15, 85)))
                    f.subs[x] = (alt, 37)
        elif cls == "read_end":
            if tumour:
                for k in range(c["n_alt"]):
                    w = "f" if k % 2 == 0 else "r"
                    off = int(rng.integers(0, 3)) if rng.random() < 0.5 else RL - 1 - int(rng.integers(0, 3))
                    f = place(x, w, off)
                    f.subs[x] = (alt, 37)
        elif cls == "low_mapq":
            if tumour:
                for k in range(c["n_alt"]):
                    f = place(x, "f" if k % 2 == 0 else "r", int(rng.integers(15, 85)))
                    f.subs[x] = (alt, 37)
                    f.mapq = int(rng.integers(0, 6))
        elif cls == "multi_mismatch":
            if tumour:
                for k in range(c["n_alt"]):
                    f = place(x, "f" if k % 2 == 0 else "r", int(rng.integers(15, 85)))
                    f.subs[x] = (alt, 37)
                    f.xmm = int(rng.integers(5, 8))
        elif cls == "low_baseq":
            if tumour:
                for k in range(c["n_alt"]):
                    f = place(x, "f" if k % 2 == 0 else "r", int(rng.integers(15, 85)))
                    f.subs[x] = (alt, int(rng.integers(4, 12)))
        elif cls in ("somatic_del", "somatic_del_shifted"):
            if tumour:
                d = c["del_len"]
                gap = c["gap_drawn0"]          # where the "aligner" draws the gap in the reads
                lo, hi = c["span0"]            # region that must lie inside the read with margin
                for f, w in cov:
                    a, b = f.fwd_span() if w == "f" else f.rev_span()
                    if a + 12 <= lo and hi + d + 12 <= b and rng.random() < c["vaf"]:
                        f.dele = (gap, d)
    return frags


def _records(ref: str, frags: list[Frag], germ: dict, protected: set, rng, sm: str):
    """Turn fragments into SAM records (tuples) for one sample."""
    bases = "ACGT"
    comp = str.maketrans("ACGT", "TGCA")
    out = []
    for f in frags:
        pair = []
        for which in ("f", "r"):
            a, b = f.fwd_span() if which == "f" else f.rev_span()
            seq = list(ref[a:b])
            qual = [37] * RL
            # germline background
            for x in range(a, b):
                g = germ.get(x)
                if g and (g[1] == 2 or f.hap == 1):
                    seq[x - a] = g[0]
            for x, (bs, q) in f.subs.items():
                if a <= x < b:
                    seq[x - a] = bs
                    qual[x - a] = q
            # paralog-like extra mismatches
            if f.xmm:
                for p in rng.choice(np.arange(RL), size=f.xmm, replace=False):
                    x = a + int(p)
                    if x in f.subs or x in protected:
                        continue
                    seq[p] = bases[(bases.index(seq[p]) + int(rng.integers(1, 4))) % 4]
            # sequencing errors (never at a candidate position)
            for p in np.nonzero(rng.random(RL) < 0.0015)[0]:
                x = a + int(p)
                if x in protected or x in f.subs:
                    continue
                seq[p] = bases[(bases.index(seq[p]) + int(rng.integers(1, 4))) % 4]
                qual[p] = int(rng.integers(6, 19)) if rng.random() < 0.9 else 37
            cigar = f"{RL}M"
            if f.dele and a < f.dele[0] and f.dele[0] + f.dele[1] < b:
                g0, d = f.dele
                left = g0 - a
                tail = ref[b:b + d]                       # keep the read length at RL
                seq = seq[:left] + seq[left + d:] + list(tail)
                qual = qual[:left] + qual[left + d:] + [37] * d
                cigar = f"{left}M{d}D{RL - left}M"
            nm = sum(1 for i, ch in enumerate(seq) if cigar == f"{RL}M" and ch != ref[a + i])
            if cigar != f"{RL}M":
                g0, d = f.dele
                left = g0 - a
                nm = d + sum(1 for i in range(left) if seq[i] != ref[a + i]) \
                    + sum(1 for i in range(left, RL) if seq[i] != ref[a + d + i])
            pair.append((which, a, cigar, "".join(seq), "".join(chr(q + 33) for q in qual), nm))
        (wf, af, cf, sf, qf, nf), (wr, ar, cr, sr, qr, nr) = pair
        dupf = 1024 if f.dup else 0
        first_f = f.r1_fwd
        pp = 2 if f.proper else 0
        flag_f = 1 + pp + 32 + (64 if first_f else 128) + dupf
        flag_r = 1 + pp + 16 + (128 if first_f else 64) + dupf
        tl = f.flen
        out.append((f.name, flag_f, af, f.mapq, cf, ar, tl, sf, qf, nf, sm))
        out.append((f.name, flag_r, ar, f.mapq, cr, af, -tl, sr, qr, nr, sm))
    out.sort(key=lambda r: (r[2], r[0]))
    return out


def _write_bam(path: Path, recs, sm: str):
    import pysam
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": CONTIG, "LN": L}],
              "RG": [{"ID": sm, "SM": sm, "PL": "ILLUMINA", "LB": f"lib_{sm}"}],
              "PG": [{"ID": "bwa", "PN": "bwa", "VN": "0.7.17-r1188", "CL": "bwa mem ref.fa R1.fq.gz R2.fq.gz"},
                     {"ID": "MarkDuplicates", "PN": "MarkDuplicates", "PP": "bwa"}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for name, flag, pos, mapq, cigar, mpos, tlen, seq, qual, nm, rg in recs:
            a = pysam.AlignedSegment(bam.header)
            a.query_name, a.flag, a.reference_id, a.reference_start = name, flag, 0, pos
            a.mapping_quality, a.cigarstring = mapq, cigar
            a.next_reference_id, a.next_reference_start, a.template_length = 0, mpos, tlen
            a.query_sequence = seq
            a.query_qualities = pysam.qualitystring_to_array(qual)
            a.set_tag("NM", nm)
            a.set_tag("RG", rg)
            bam.write(a)
    pysam.index(str(path))


INSTRUCTION = """\
You are the variant-review analyst of a cancer genomics laboratory. An automated caller produced
somatic variant candidates for one tumour/normal pair; before reporting, every candidate must be
reviewed against the aligned reads following the laboratory's review SOP.

Task directory:
- `/task`

Inputs (read-only) under `/task/input`:
- `bam/tumor.bam`, `bam/normal.bam` (+ `.bai`) - coordinate-sorted alignments, duplicates marked
- `reference/ref.fa` (+ `.fai`) - the reference the reads were aligned to
- `candidates.vcf` - the {n} candidates to review (`ID` column = variant id); `candidates.bed` holds
  the same loci for display
- `review_sop.md` - the laboratory's review SOP (evidence definition, failure tags, their order)

Available software (on `PATH`): `samtools`, Python 3.12 with `pysam`, numpy, pandas; IGV 2.19
(`/task/software/open_igv.sh` starts it on the virtual display with the reference, both BAMs and
the candidate track loaded; its batch port 60151 is enabled). If your harness can operate a GUI, the
MCP server `desktop` provides screenshot, click, type and key tools for that display. The review
may be done in IGV, from the terminal, or both. There is no internet access.

Your task: review every candidate according to `review_sop.md`.

Required outputs under `/task/output`:
- `review.tsv` - tab-separated, one row per candidate, columns
  `variant_id`, `call`, `tag`, `tumor_alt`, `tumor_depth`, `tumor_vaf`, `normal_alt`, `normal_depth`.
  `call` is `somatic` or `fail`; `tag` is `PASS` for a somatic call, otherwise exactly one failure
  tag from the SOP; the count columns are in molecules as defined in the SOP and
  `tumor_vaf = tumor_alt / tumor_depth`.
- `igv_session.xml` - an IGV session file that opens the reference, both BAMs and the candidate
  track, so that a colleague can re-open your review (paths as under `/task/input`).

Constraints:
- Do not modify anything under `/task/input`; write deliverables only under `/task/output`.
- Do not ask for confirmation. Execute directly.
"""

SOP = """\
# Somatic variant review SOP (laboratory document R-07)

Review each candidate in the tumour and in the matched normal. The reviewer records one call per
candidate: `somatic`, or `fail` with exactly one failure tag.

## What counts as evidence
- Evidence is counted in sequenced molecules. Reads marked as duplicates add no evidence and no
  depth. Secondary and supplementary alignments are ignored.
- Depth at a candidate = molecules whose read covers the candidate position.
- A molecule supports a single-nucleotide candidate if its read shows the alternate base at that
  position, whatever the base quality or mapping quality (those are judged by the tags below).
- For insertions and deletions, judge what the read sequence shows, not where the aligner drew the
  gap: a read supports the candidate if it carries a gap that produces the same sequence as the
  candidate allele.

## Failure tags, in order of precedence (record the first that applies)
1. `GERMLINE` - the alternate allele is seen in the normal in at least 2 molecules and at least 5%
   of the normal depth.
2. `LOW_SUPPORT` - fewer than 4 supporting molecules in the tumour.
3. `LOW_BASE_QUALITY` - at least three quarters of the supporting molecules show the alternate
   base with base quality below 20.
4. `LOW_MAPPING_QUALITY` - at least three quarters of the supporting molecules have mapping
   quality below 20.
5. `MULTI_MISMATCH` - the supporting reads carry, on average, four or more further mismatches to
   the reference within the same read.
6. `STRAND_BIAS` - at least 90% of the supporting reads are on one strand although reads on both
   strands cover the position.
7. `READ_END` - at least three quarters of the supporting reads show the alternate base within 5
   bases of either end of the read.

A candidate to which no tag applies is called `somatic` (tag `PASS`).
"""

JUDGE = """\
**Where the truth comes from.** The BAM files are synthesised read by read. For every candidate the
generator plants one mechanism (true clonal/subclonal somatic variant, germline variant, too few
molecules, duplicates of a single molecule, strand-specific artefact, read-end artefact, low base
quality, low mapping quality, paralog-like reads with many mismatches, or a somatic deletion in a
repeat) with parameters far from every SOP threshold, so the SOP outcome of each candidate is known
before any read exists. `reference/truth.json` stores the planted call and tag and, for true
somatic variants, the realised molecule-level allele fraction.

**Hard gates (score 0).** `review.tsv` missing/unparseable; a candidate missing or listed twice;
`call` not in {somatic, fail}; `tag` outside the SOP vocabulary; `PASS` used with `fail` or a
failure tag with `somatic`; non-numeric count columns.

**Components.**

| weight | component | metric | full credit | zero |
|---|---|---|---|---|
| 0.45 | calls | Matthews correlation of `call` with the planted truth | >= 0.96 | <= 0.70 |
| 0.30 | failure_tags | fraction of planted failures given the planted tag | >= 0.97 | <= 0.75 |
| 0.20 | allele_fraction | fraction of planted somatic variants called `somatic` with `tumor_vaf` within 0.05 of the realised molecule-level fraction | >= 0.95 | <= 0.70 |
| 0.05 | igv_session | well-formed IGV session referencing the reference, both BAMs and the candidate track | all four | none |

**Why these bands.** The truth is exact, so a correct review reproduces every call and tag; the
bands exist to give graded reward, not to absorb noise: with 48 candidates one wrong call gives
MCC ~0.957 and still (almost) full credit, while ignoring a whole class of four candidates gives
MCC <= 0.85. The allele-fraction tolerance (0.05 absolute) covers legitimate counting conventions
(depth at the anchor base versus across a repeat, reads versus molecules when no pairs overlap, as
here) but not position-exact indel counts, which report 0. The session file is worth only 0.05 because it can be checked for structure, not for use.

**What is deliberately not scored.** Whether IGV was actually opened; screenshots.
"""


def generate(seed: int, level: int, workdir: Path) -> TaskBuild:
    rng = np.random.default_rng(7_000 + seed)
    ref = _reference(seed)

    # ---- plan candidates ---------------------------------------------------------------
    plan = [("somatic_snv", 12), ("somatic_del", 2), ("germline", 5), ("low_support", 4)]
    if level >= 2:
        plan += [("strand_bias", 5), ("read_end", 4), ("low_mapq", 3), ("multi_mismatch", 3), ("low_baseq", 2)]
    else:
        plan += [("germline", 6), ("low_support", 5), ("somatic_snv", 6)]
    if level >= 3:
        plan += [("dup_only", 4), ("somatic_del_shifted", 4)]
    else:
        plan += [("somatic_snv", 4), ("low_support", 4)]
    if level >= 4:      # each mnv entry is one locus with two adjacent candidates
        plan += [("mnv_cis", 4), ("mnv_trans", 4), ("long_insert", 4)]
    classes = [c for c, k in plan for _ in range(k)]
    classes = [classes[i] for i in rng.permutation(len(classes))]
    n = len(classes)

    # homopolymers for the repeat-context deletions are chosen first; everything else avoids them
    homo = [h for h in _homopolymers(ref) if 3000 < h[0] < L - 3000]
    homo = [homo[i] for i in rng.permutation(len(homo))]
    chosen_h = []
    for h in homo:
        if all(abs(h[0] - u[0]) > 3000 for u in chosen_h):
            chosen_h.append(h)
        if len(chosen_h) == classes.count("somatic_del_shifted"):
            break
    if len(chosen_h) < classes.count("somatic_del_shifted"):
        raise RuntimeError("not enough homopolymers in this window; change the seed")
    slots = np.linspace(2500, L - 2500, n).astype(int)
    cands = []
    hi_iter = iter(chosen_h)
    for i, cls in enumerate(classes):
        x = int(slots[i] + rng.integers(-200, 200))
        c = {"id": "", "cls": cls}
        if cls == "somatic_del_shifted":
            # deletion of 1-2 bases inside a homopolymer: VCF is left-aligned, reads are right-aligned
            hs, hl = next(hi_iter)
            d = int(rng.integers(1, 3))
            c.update({"pos0": hs - 1, "ref": ref[hs - 1:hs + d], "alt": ref[hs - 1], "del_len": d,
                      "gap_drawn0": hs + hl - d, "span0": (hs - 1, hs + hl), "vaf": float(rng.uniform(0.22, 0.42))})
        else:
            for _ in range(200):     # keep >= 500 bp to every repeat-context deletion and placed candidate
                taken = [h[0] for h in chosen_h] + [q["pos0"] for q in cands if "pos0" in q]
                if all(abs(x - t) >= 500 for t in taken):
                    break
                x = int(slots[i] + rng.integers(-550, 550))
            else:
                if level < 4:
                    raise RuntimeError("could not place candidates; change the seed")
                for _ in range(3000):   # level 4 has more loci in the same window: accept 380 bp before giving up
                    x = int(slots[i] + rng.integers(-550, 550))
                    if all(abs(x - t) >= 380 for t in taken):
                        break
                else:
                    raise RuntimeError("could not place candidates; change the seed")
            if cls == "somatic_del":
                d = int(rng.integers(2, 5))
                while not (ref[x] != ref[x + d] and ref[x - 1] != ref[x + d - 1]):
                    x += 1          # unique representation: the gap cannot be shifted either way
                c.update({"pos0": x - 1, "ref": ref[x - 1:x + d], "alt": ref[x - 1], "del_len": d,
                          "gap_drawn0": x, "span0": (x - 1, x + d), "vaf": float(rng.uniform(0.2, 0.4))})
            else:
                r = ref[x]
                alt = "ACGT"[("ACGT".index(r) + int(rng.integers(1, 4))) % 4]
                c.update({"pos0": x, "ref": r, "alt": alt})
                if cls == "somatic_snv":
                    c["vaf"] = float(rng.uniform(0.25, 0.45)) if rng.random() < 0.6 else float(rng.uniform(0.09, 0.16))
                elif cls == "low_support":
                    c["n_alt"] = int(rng.integers(1, 3))
                elif cls == "dup_only":
                    c["n_orig"], c["n_dup"] = int(rng.integers(2, 4)), int(rng.integers(3, 5))
                elif cls == "strand_bias":
                    c["n_alt"], c["strand"] = int(rng.integers(9, 14)), ("f" if rng.random() < 0.5 else "r")
                elif cls in ("read_end", "low_mapq", "multi_mismatch", "low_baseq"):
                    c["n_alt"] = int(rng.integers(9, 14))
                elif cls == "long_insert":
                    c["n_alt"] = int(rng.integers(10, 14))
                elif cls in ("mnv_cis", "mnv_trans"):
                    r2 = ref[x + 1]
                    alt2 = "ACGT"[("ACGT".index(r2) + int(rng.integers(1, 4))) % 4]
                    cis = cls == "mnv_cis"
                    c.update({"cls": "mnv", "role": "L", "phase": "cis" if cis else "trans", "mate_alt": alt2,
                              "vaf": float(rng.uniform(0.25, 0.42)) if cis else float(rng.uniform(0.16, 0.24))})
                    cands.append(c)
                    c = {"id": "", "cls": "mnv", "role": "R", "phase": c["phase"], "pos0": x + 1, "ref": r2, "alt": alt2}
        cands.append(c)
    cands.sort(key=lambda c: c["pos0"])
    for i, c in enumerate(cands):
        c["id"] = f"cand_{i + 1:03d}"
    gaps = [b["pos0"] - a["pos0"] for a, b in zip(cands, cands[1:])
            if not (a["cls"] == "mnv" and a.get("role") == "L" and b["pos0"] == a["pos0"] + 1)]
    if min(gaps) < 300:
        raise RuntimeError("candidates too close; change the seed")
    if level >= 4:
        from . import level4
        level4.plan_extras(cands, rng)

    protected = set()
    for c in cands:
        lo = c["pos0"] - 2
        hi = (c["span0"][1] + 3) if "span0" in c else c["pos0"] + 3
        protected.update(range(lo, hi))
    # germline background SNPs shared by tumour and normal (not near candidates)
    germ = {}
    while len(germ) < 45:
        x = int(rng.integers(500, L - 500))
        if any(abs(x - c["pos0"]) < 150 for c in cands):
            continue
        germ[x] = ("ACGT"[("ACGT".index(ref[x]) + int(rng.integers(1, 4))) % 4], 2 if rng.random() < 0.3 else 1)

    # ---- reads ---------------------------------------------------------------------------
    bam_dir = workdir / "bam"
    bam_dir.mkdir(parents=True, exist_ok=True)
    realised = {}
    for sm, depth, tumour in (("tumor", 90, True), ("normal", 45, False)):
        frags = _simulate(ref, cands, depth, tumour, rng, sm[0].upper())
        if level >= 4 and tumour:
            level4.mix(frags, cands, rng, Frag, RL)
        # unrelated flagged duplicates for realism
        for f in [frags[i] for i in rng.permutation(len(frags))[:int(0.06 * len(frags))]]:
            d = Frag(f.start, f.flen, f.hap, f.r1_fwd, f.name + "d")
            d.subs, d.dele, d.xmm, d.mapq, d.dup = dict(f.subs), f.dele, f.xmm, f.mapq, True
            frags.append(d)
        for c in cands:   # realised molecule-level evidence (non-duplicate)
            x = c["pos0"] if "del_len" not in c else c["span0"][0]
            cov = [f for f in frags if not f.dup and f.covering(x)]
            if "del_len" in c:
                alt = sum(1 for f in cov if f.dele is not None)
            else:
                alt = sum(1 for f in cov if c["pos0"] in f.subs)
            realised[(c["id"], sm)] = (alt, len(cov))
        recs = _records(ref, frags, germ, protected, rng, sm)
        _write_bam(bam_dir / f"{sm}.bam", recs, sm)

    # ---- truth ---------------------------------------------------------------------------
    tag_of = {"somatic_snv": "PASS", "somatic_del": "PASS", "somatic_del_shifted": "PASS", "germline": "GERMLINE",
              "low_support": "LOW_SUPPORT", "dup_only": "LOW_SUPPORT", "strand_bias": "STRAND_BIAS",
              "read_end": "READ_END", "low_mapq": "LOW_MAPPING_QUALITY", "multi_mismatch": "MULTI_MISMATCH",
              "low_baseq": "LOW_BASE_QUALITY", "long_insert": "DISCORDANT_PAIRS"}
    truth = {"variants": {}}
    for c in cands:
        if c["cls"] == "mnv":   # level 4: the right-hand candidate of a same-molecule pair is merged into the left
            tag = "MNV_PARTNER" if (c["phase"] == "cis" and c["role"] == "R") else "PASS"
        else:
            tag = tag_of[c["cls"]]
        ta, td = realised[(c["id"], "tumor")]
        na, nd = realised[(c["id"], "normal")]
        if tag == "PASS" and ta < 6:           # keep true variants unambiguous
            raise RuntimeError(f"somatic candidate {c['id']} realised with only {ta} molecules; change the seed")
        truth["variants"][c["id"]] = {"call": "somatic" if tag == "PASS" else "fail", "tag": tag, "class": c["cls"],
                                      "tumor_alt": ta, "tumor_depth": td, "tumor_vaf": ta / td,
                                      "normal_alt": na, "normal_depth": nd}

    fasta = f">{CONTIG}\n" + "\n".join(ref[i:i + 60] for i in range(0, L, 60)) + "\n"
    (workdir / "ref.fa").write_text(fasta)
    import pysam
    pysam.faidx(str(workdir / "ref.fa"))
    vcf = ["##fileformat=VCFv4.2", f"##contig=<ID={CONTIG},length={L}>", "##source=candidate_caller_v3",
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    bed = []
    for c in cands:
        vcf.append(f"{CONTIG}\t{c['pos0'] + 1}\t{c['id']}\t{c['ref']}\t{c['alt']}\t.\t.\t.")
        bed.append(f"{CONTIG}\t{c['pos0']}\t{c['pos0'] + len(c['ref'])}\t{c['id']}")
    if level >= 4:
        return level4.build(seed=seed, level=level, workdir=workdir, ref=ref, cands=cands, truth=truth,
                            bam_dir=bam_dir, vcf="\n".join(vcf) + "\n", bed="\n".join(bed) + "\n",
                            ref_url=REF_URL, ref_sha256=REF_SHA256, contig=CONTIG)

    open_igv = """#!/bin/bash
# Start IGV on the virtual display with the reference, both BAMs and the candidate track.
export DISPLAY="${DISPLAY:-:99}"
bash /task/software/desktop/start_desktop.sh
mkdir -p "$HOME/igv"
if [ ! -s "$HOME/igv/prefs.properties" ]; then
  printf 'PROVISIONING_URL_DEFAULT=/task/software/igv/oauth_offline.json\\nIGV.genome.sequence.dir=/task/software/igv/genomes.tsv\\nCHECK_VERSION=false\\n' > "$HOME/igv/prefs.properties"
fi
nohup igv -g /task/input/reference/ref.fa /task/input/bam/tumor.bam,/task/input/bam/normal.bam,/task/input/candidates.bed "$@" > /tmp/igv.log 2>&1 &
echo "IGV starting (log: /tmp/igv.log); batch port 60151"
"""
    oauth = json.dumps({"client_id": "offline", "client_secret": "offline", "auth_uri": "http://127.0.0.1/auth",
                        "token_uri": "http://127.0.0.1/token", "project_id": "offline"})

    # ---- battery ---------------------------------------------------------------------------
    ids = [c["id"] for c in cands]
    tv = truth["variants"]

    def review(fn):
        rows = ["variant_id\tcall\ttag\ttumor_alt\ttumor_depth\ttumor_vaf\tnormal_alt\tnormal_depth"]
        for i in ids:
            call, tag, vaf = fn(i)
            t = tv[i]
            rows.append(f"{i}\t{call}\t{tag}\t{t['tumor_alt']}\t{t['tumor_depth']}\t{vaf:.4f}\t{t['normal_alt']}\t{t['normal_depth']}")
        return "\n".join(rows) + "\n"

    session = ('<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'
               f'<Session genome="/task/input/reference/ref.fa" locus="{CONTIG}:1-2000" version="8">\n <Resources>\n'
               '  <Resource path="/task/input/bam/tumor.bam"/>\n  <Resource path="/task/input/bam/normal.bam"/>\n'
               '  <Resource path="/task/input/candidates.bed"/>\n </Resources>\n</Session>\n')
    prng = np.random.default_rng(seed + 991)
    one_wrong = ids[int(prng.integers(0, n))]

    def flip(i):
        t = tv[i]
        if i == one_wrong:
            return ("fail", "LOW_SUPPORT", t["tumor_vaf"]) if t["call"] == "somatic" else ("somatic", "PASS", t["tumor_vaf"])
        return t["call"], t["tag"], t["tumor_vaf"] + float(prng.normal(0, 0.008))

    probes = [
        Probe("planted_truth", "review written directly from the planted truth", 0.999,
              files={"review.tsv": review(lambda i: (tv[i]["call"], tv[i]["tag"], tv[i]["tumor_vaf"])),
                     "igv_session.xml": session}),
        Probe("within_tolerance", "one candidate reviewed wrongly, allele fractions off by ~0.01", 0.97,
              files={"review.tsv": review(flip), "igv_session.xml": session}),
        Probe("independent_solution", "second implementation: samtools-view text parsing with its own CIGAR walker",
              0.99, script="SOLUTION=2 bash /solution/solve.sh"),
    ]
    helper = {"review.py": ASSETS / "review.py"}
    attacks = [
        Attack("trust_the_caller", "every candidate called somatic", 0.30,
               files={"review.tsv": review(lambda i: ("somatic", "PASS", tv[i]["tumor_vaf"])), "igv_session.xml": session}),
        Attack("fail_everything", "every candidate failed as LOW_SUPPORT", 0.30,
               files={"review.tsv": review(lambda i: ("fail", "LOW_SUPPORT", 0.0)), "igv_session.xml": session}),
        Attack("counts_only_review", "only GERMLINE and LOW_SUPPORT are checked (no read-level inspection)", 0.75,
               script="python3 /alegen_script/review.py --counts-only", script_files=helper),
    ]
    if level >= 3:
        attacks += [
            Attack("duplicates_counted", "duplicate-flagged reads counted as evidence", 0.75,
                   script="python3 /alegen_script/review.py --keep-duplicates", script_files=helper),
            Attack("position_exact_indels", "a deletion only counts when the gap is drawn at the candidate position", 0.75,
                   script="python3 /alegen_script/review.py --exact-indel", script_files=helper),
        ]

    traps = []
    if level >= 2:
        traps.append(Trap("artefact_classes", "Five artefact mechanisms look like real variants in a pile-up count "
                          "(strand, read end, base quality, mapping quality, multi-mismatch reads).",
                          "review_sop.md defines the tags; per-read inspection (IGV colouring by strand/MAPQ, or "
                          "pysam per-read attributes) shows each one", "Approach / incomplete workflow"))
    if level >= 3:
        traps += [
            Trap("duplicate_only_support", "8-15 supporting reads on both strands that are only 2-3 molecules plus their flagged duplicates.",
                 "flag 0x400 in the BAM; IGV hides duplicates by default; SOP counts molecules",
                 "Execution / implementation bug"),
            Trap("indel_representation", "Four true deletions in homopolymers are left-aligned in the candidate list "
                 "but right-aligned in the reads; position-exact counting finds zero support.",
                 "the gap is visible a few bases downstream in IGV; SOP says to judge the sequence, not the gap position",
                 "Understanding / domain knowledge gap"),
        ]

    return TaskBuild(
        task_id=f"variant-review-s{seed:04d}-l{level}", family="variant_review", seed=seed, level=level,
        title="Somatic variant review against tumour/normal alignments (IGV)",
        summary="Review somatic variant candidates read by read in tumour/normal BAMs following a laboratory SOP, "
                "using IGV and/or the terminal, and deliver calls, failure tags and an IGV session.",
        subdomain="3.2 Genomics & Sequence Analysis",
        software=["IGV (GUI)", "samtools", "pysam", "Python", "desktop MCP server"],
        instruction_md=INSTRUCTION.format(n=n),
        agent_dockerfile=dockerlib.agent_dockerfile(conda=["igv=2.19.8", "samtools=1.25", "pysam=0.24.1"],
                                                    pip=desktop.PIP, extra=desktop.APT),
        verifier_dockerfile=dockerlib.verifier_dockerfile(pip=[]),
        input_files={"bam/tumor.bam": bam_dir / "tumor.bam", "bam/tumor.bam.bai": bam_dir / "tumor.bam.bai",
                     "bam/normal.bam": bam_dir / "normal.bam", "bam/normal.bam.bai": bam_dir / "normal.bam.bai",
                     "reference/ref.fa": workdir / "ref.fa", "reference/ref.fa.fai": workdir / "ref.fa.fai",
                     "candidates.vcf": "\n".join(vcf) + "\n", "candidates.bed": "\n".join(bed) + "\n",
                     "review_sop.md": SOP},
        software_files={**desktop.FILES, "open_igv.sh": open_igv, "igv/oauth_offline.json": oauth,
                        "igv/genomes.tsv": "id\tname\turl\n"},
        solution_files={"solve.sh": ASSETS / "solve.sh", "review.py": ASSETS / "review.py",
                        "review_alt.py": ASSETS / "review_alt.py"},
        tests_files={"score.py": ASSETS / "score.py", "reference/truth.json": json.dumps(truth)},
        traps=traps, attacks=attacks, probes=probes,
        capabilities=["read-level evidence inspection (GUI or CLI)", "following a multi-rule SOP with precedence",
                      "molecule vs read counting", "indel representation", "GUI-native deliverable (IGV session)"],
        judge_explanation_md=JUDGE, pass_threshold=0.9,
        oracle_steps=6, oracle_runtime_min=1.0, needs_gui=True, cpus=4, memory_mb=8000, storage_mb=10000,
        mcp_servers=[desktop.MCP_SERVER],
        extra_metadata={"reference_source": REF_URL, "reference_sha256": REF_SHA256},
    )
