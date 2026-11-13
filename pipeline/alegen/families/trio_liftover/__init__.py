"""Family: rare-disease trio, cross-build variant prioritisation (hg19 call set -> GRCh38 -> SOP).

Truth-by-construction: every variant is planted on GRCh38 as a normalised (left-aligned, biallelic)
allele with planted trio genotypes; the provider VCF shown to the agent is *derived* from that truth
by mapping each allele back to hg19 through the UCSC chain (only inside ungapped, plus-strand blocks
whose flanking sequence is identical in both assemblies) and then degrading the representation
(b37 contig name, right-shifted / padded alleles, multi-allelic records). The candidate set is the
planted causal set; every other record fails at least one SOP rule by construction.

Tool-derived quantities are only used as *guards*: a planted allele is scored only if Picard
LiftoverVcf and bcftools +liftover both reproduce the planted GRCh38 key (otherwise gray zone), and
generation is retried if `bcftools csq` disagrees with a designed consequence of a causal/decoy
variant or if the SOP applied to csq consequences does not return exactly the planted causal set.

Levels: 1 = declared hg19 header, canonical representation, no inheritance decoys;
        2 = + undeclared build, b37 contig name, non-canonical representation;
        3 = + genotype-quality and phase decoys.
"""
from __future__ import annotations

import bisect
import gzip
import json
from pathlib import Path

import numpy as np

from ... import dockerlib
from ...spec import Attack, Probe, TaskBuild, Trap
from . import refdata

ASSETS = Path(__file__).parent / "assets"
COMP = str.maketrans("ACGT", "TGCA")
BASES = "ACGT"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: _AA[16 * i + 4 * j + k] for i, a in enumerate("TCAG") for j, b in enumerate("TCAG")
         for k, c in enumerate("TCAG")}
QUALIFYING = ("missense", "stop_gained", "frameshift", "splice_acceptor", "splice_donor")
MIN_DP, MIN_GQ, AB_LO, AB_HI, AF_CUT = 10, 20, 0.25, 0.75, 0.001
_REF: refdata.Ref | None = None
L4_ORACLE_STEPS = 24      # distinct tool invocations of assets/trio_pipeline_l4.py (counted in its docstring)


def _ref() -> refdata.Ref:
    global _REF
    if _REF is None:
        _REF = refdata.Ref()
    return _REF


# ------------------------------------------------------------------------------- allele algebra
def normalise(seq: str, pos: int, ref: str, alt: str):
    """Left-align and trim (Tan et al. 2015). pos is 0-based."""
    changed = True
    while changed:
        changed = False
        if ref and alt and ref[-1] == alt[-1]:
            ref, alt, changed = ref[:-1], alt[:-1], True
        if not ref or not alt:
            pos -= 1
            ref, alt, changed = seq[pos] + ref, seq[pos] + alt, True
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def shift_max(seq: str, pos: int, ref: str, alt: str, cap: int = 40) -> int:
    """How many bases a normalised indel can be moved to the right without changing the haplotype."""
    k = 0
    if len(ref) > len(alt):
        start, n = pos + 1, len(ref) - 1
        while k < cap and seq[start + k] == seq[start + n + k]:
            k += 1
    elif len(alt) > len(ref):
        s = alt[1:]
        while k < cap and seq[pos + 1 + k] == s[k % len(s)]:
            k += 1
    return k


def shifted(seq: str, pos: int, ref: str, alt: str, k: int):
    start = pos + 1 + k
    if len(ref) > len(alt):
        n = len(ref) - 1
        return start - 1, seq[start - 1:start + n], seq[start - 1]
    s = alt[1:]
    r = k % len(s)
    return start - 1, seq[start - 1], seq[start - 1] + s[r:] + s[:r]


def padded(seq: str, pos: int, ref: str, alt: str, n: int):
    tail = seq[pos + len(ref):pos + len(ref) + n]
    return pos, ref + tail, alt + tail


# ------------------------------------------------------------------------------- gene models
class GeneModel:
    def __init__(self, g: dict, seq: str):
        self.name, self.strand, self.span = g["name"], g["strand"], (g["start"] - 1, g["end"])
        self.segs = g["cds"]
        self.exons = g["exons"]
        asc = [p for s, e in self.segs for p in range(s, e)]
        self.asc = asc
        self.gpos = asc if self.strand == "+" else asc[::-1]
        cs = "".join(seq[p] for p in self.gpos)
        self.cds = cs if self.strand == "+" else cs.translate(COMP)
        n = len(self.cds) // 3
        prot = "".join(CODON.get(self.cds[3 * i:3 * i + 3], "X") for i in range(n))
        self.ncodon = n
        self.valid = (len(self.cds) % 3 == 0 and n >= 80 and prot[0] == "M" and prot[-1] == "*"
                      and "*" not in prot[:-1] and "X" not in prot)
        self._starts = [s for s, _ in self.segs]

    def seg(self, p: int):
        i = bisect.bisect_right(self._starts, p) - 1
        if i >= 0 and self.segs[i][0] <= p < self.segs[i][1]:
            return self.segs[i]
        return None

    def edge(self, p: int) -> int:
        s = self.seg(p)
        return -1 if s is None else min(p - s[0], s[1] - 1 - p)

    def protein_after(self, seq: str, pos: int, ref: str, alt: str) -> str:
        i = bisect.bisect_left(self.asc, pos)
        g = "".join(seq[p] for p in self.asc)
        g = g[:i] + alt + g[i + len(ref):]
        c = g if self.strand == "+" else g.translate(COMP)[::-1]
        return "".join(CODON.get(c[3 * k:3 * k + 3], "X") for k in range(len(c) // 3))


def design_snv(gm: GeneModel, seq: str, rng, cls: str):
    for _ in range(400):
        c = int(rng.integers(8, gm.ncodon - 8))
        gps = gm.gpos[3 * c:3 * c + 3]
        if any(gm.edge(p) < 6 for p in gps):
            continue
        codon = gm.cds[3 * c:3 * c + 3]
        opts = []
        for k in range(3):
            for b in BASES:
                if b == codon[k]:
                    continue
                naa = CODON[codon[:k] + b + codon[k + 1:]]
                got = "synonymous" if naa == CODON[codon] else ("stop_gained" if naa == "*" else "missense")
                if got == cls:
                    opts.append((k, b))
        if not opts:
            continue
        k, b = opts[int(rng.integers(len(opts)))]
        p = gps[k]
        return p, seq[p], (b if gm.strand == "+" else b.translate(COMP))
    return None


def design_frameshift(gm: GeneModel, seq: str, rng, need_shift: bool = False, inframe: bool = False):
    prot0 = "".join(CODON[gm.cds[3 * i:3 * i + 3]] for i in range(gm.ncodon))
    for _ in range(600):
        c = int(rng.integers(8, gm.ncodon - 14))
        p = gm.gpos[3 * c + int(rng.integers(0, 3))]
        if gm.edge(p) < 16:
            continue
        n = 3 if inframe else int(rng.integers(1, 3))
        if inframe or rng.random() < 0.6:
            ref, alt = seq[p:p + 1 + n], seq[p]
        else:
            ref, alt = seq[p], seq[p] + "".join(BASES[int(x)] for x in rng.integers(0, 4, size=n))
        pos, r, a = normalise(seq, p, ref, alt)
        k = shift_max(seq, pos, r, a)
        if (need_shift and k < 1) or k > 10:
            continue
        lo, hi = pos, pos + len(r) + k
        if gm.seg(lo) is None or gm.seg(lo) != gm.seg(hi) or gm.edge(lo) < 8 or gm.edge(hi) < 8:
            continue
        if not inframe:
            prot1 = gm.protein_after(seq, pos, r, a)
            d = next((i for i in range(min(len(prot0), len(prot1))) if prot0[i] != prot1[i]), None)
            if d is None or "*" in prot1[d:d + 3] or d < 5:
                continue
        return pos, r, a
    return None


def design_splice(gm: GeneModel, seq: str, rng):
    introns = [(gm.segs[i][1], gm.segs[i + 1][0]) for i in range(len(gm.segs) - 1)
               if gm.segs[i + 1][0] - gm.segs[i][1] >= 120]
    for _ in range(60):
        if not introns:
            return None
        a, b = introns[int(rng.integers(len(introns)))]
        donor = rng.random() < 0.5
        if gm.strand == "+":
            sites = [(a, "G"), (a + 1, "T")] if donor else [(b - 2, "A"), (b - 1, "G")]
        else:
            sites = [(b - 1, "C"), (b - 2, "A")] if donor else [(a, "C"), (a + 1, "T")]
        if any(seq[p] != x for p, x in sites):
            continue
        p = sites[int(rng.integers(2))][0]
        alt = [x for x in BASES if x != seq[p]][int(rng.integers(3))]
        return p, seq[p], alt, ("splice_donor" if donor else "splice_acceptor")
    return None


# ------------------------------------------------------------------------------- genotypes
def fmt_good(rng, g: int, clean: bool = False):
    dp = int(np.clip(rng.poisson(34), 16, 70))
    if g == 0:
        alt = 1 if (not clean and dp >= 28 and rng.random() < 0.03) else 0
        return dp - alt, alt, dp, min(99, 3 * (dp - alt))
    if g == 2:
        ref = 1 if (not clean and dp >= 28 and rng.random() < 0.03) else 0
        return ref, dp - ref, dp, min(99, 3 * (dp - ref))
    while True:
        alt = int(rng.binomial(dp, 0.5))
        if 0.34 <= alt / dp <= 0.66:
            return dp - alt, alt, dp, 99


def fmt_lowdp(rng, g: int):
    dp = int(rng.integers(3, 7))
    if g == 0:
        return dp, 0, dp, 3 * dp
    if g == 2:
        return 0, dp, dp, 3 * dp
    alt = int(rng.integers(1, dp))
    return dp - alt, alt, dp, int(rng.integers(5, 19))


def gt_str(g: int) -> str:
    return ("0/0", "0/1", "1/1")[g]


def apply_sop(variants, cons_of, panel, quality=True, phase=True):
    """The lab SOP as code, on planted data. variants: dicts with gt/fmt per role P, M, F.
    cons_of(v) -> list of (gene, class). Returns {index: (gene, class, model)}."""
    passing = []
    for i, v in enumerate(variants):
        if v["pos"] is None:
            continue
        if quality:
            if any(v["fmt"][r][2] < MIN_DP or v["fmt"][r][3] < MIN_GQ for r in "PMF"):
                continue
            if v["gt"]["P"] == 1:
                ref, alt = v["fmt"]["P"][0], v["fmt"]["P"][1]
                if ref + alt == 0 or not AB_LO <= alt / (ref + alt) <= AB_HI:
                    continue
        if v["af"] is not None and v["af"] >= AF_CUT:
            continue
        for gene, cls in cons_of(v):
            if cls in QUALIFYING and gene in panel:
                passing.append((i, gene, cls))
    out = {}
    by_gene = {}
    for i, gene, cls in passing:
        g = variants[i]["gt"]
        if (g["P"], g["M"], g["F"]) == (1, 0, 0):
            out[i] = (gene, cls, "de_novo")
        elif (g["P"], g["M"], g["F"]) == (2, 1, 1):
            out[i] = (gene, cls, "homozygous_recessive")
        elif g["P"] == 1:
            by_gene.setdefault(gene, []).append((i, cls))
    for gene, lst in by_gene.items():
        mat = [(i, c) for i, c in lst if (variants[i]["gt"]["M"], variants[i]["gt"]["F"]) == (1, 0)]
        pat = [(i, c) for i, c in lst if (variants[i]["gt"]["M"], variants[i]["gt"]["F"]) == (0, 1)]
        hit = (mat + pat) if (mat and pat) else ([] if phase else (lst if len(lst) >= 2 else []))
        for i, c in hit:
            out[i] = (gene, c, "compound_het")
    return out


# ------------------------------------------------------------------------------- level 4 helpers
_MISMATCH = None


def real_mismatches(R):
    """Positions inside aligned blocks of the main chain where the two shipped assemblies carry different
    bases: [(pos38, pos19, distance to the nearest block edge)], 0-based."""
    global _MISMATCH
    if _MISMATCH is None:
        a19 = np.frombuffer(R.hg19.encode(), dtype=np.uint8)
        a38 = np.frombuffer(R.hg38.encode(), dtype=np.uint8)
        out = []
        for t, q, size in R.blocks:
            if not any(q < b and q + size > a for a, b in R.windows38):
                continue
            x, y = a19[t:t + size], a38[q:q + size]
            for i in np.nonzero((x != y) & (x != 78) & (y != 78))[0]:
                out.append((q + int(i), t + int(i), int(min(i, size - 1 - i))))
        _MISMATCH = sorted(out)
    return _MISMATCH


def revcomp(x: str) -> str:
    return x.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def invert_in_chain(text: str, t_a: int, t_b: int, q_a: int, tsize: int, qsize: int, new_id: str):
    """Chain file text in which the plus/plus block that aligns target [t_a, t_b) to query [q_a, ...) is split
    and the interval is carried by an extra single-block chain on the minus strand of the query (UCSC
    convention: minus-strand query coordinates count from the end of the query sequence). None if no single
    block holds the interval with at least one aligned base on either side."""
    length = t_b - t_a
    out, done, lines, i = [], False, text.split("\n"), 0
    while i < len(lines):
        f = lines[i].split()
        out.append(lines[i])
        i += 1
        if f and f[0] == "chain" and f[4] == "+" and f[9] == "+":
            t, q = int(f[5]), int(f[10])
            while i < len(lines) and lines[i].split():
                g = lines[i].split()
                size = int(g[0])
                if not done and t < t_a and t_b < t + size and q + (t_a - t) == q_a:
                    out.append(f"{t_a - t} {length} {length}")
                    out.append(" ".join([str(t + size - t_b)] + g[1:]))
                    done = True
                else:
                    out.append(lines[i])
                if len(g) == 3:
                    t, q = t + size + int(g[1]), q + size + int(g[2])
                i += 1
    if not done:
        return None
    return ("\n".join(out).rstrip("\n") + "\n\n"
            + f"chain {50 * length} chr22 {tsize} + {t_a} {t_b} chr22 {qsize} - {qsize - (q_a + length)} {qsize - q_a} {new_id}\n"
            + f"{length}\n\n")


def snv_class(a0: str, a1: str) -> str:
    return "synonymous" if a1 == a0 else ("stop_gained" if a1 == "*" else "missense")


def codon_pairs(gm: GeneModel, usable):
    """All pairs of single-base changes inside one codon of a gene: tuples
    ((pos, alt, class alone), (pos, alt, class alone), class of the codon with both changes, aa0, aa1, aa2, aa12)
    with genomic 0-based positions and plus-strand alleles; both orders are listed."""
    out = []
    comp = (lambda b: b) if gm.strand == "+" else (lambda b: b.translate(COMP))
    for c in range(8, gm.ncodon - 8):
        gps = gm.gpos[3 * c:3 * c + 3]
        if abs(gps[2] - gps[0]) != 2 or any(gm.edge(p) < 6 or not usable(p) for p in gps):
            continue
        codon = gm.cds[3 * c:3 * c + 3]
        a0 = CODON[codon]
        for i in range(3):
            for j in range(3):
                if i == j:
                    continue
                for bi in BASES:
                    for bj in BASES:
                        if bi == codon[i] or bj == codon[j]:
                            continue
                        c1 = codon[:i] + bi + codon[i + 1:]
                        c2 = codon[:j] + bj + codon[j + 1:]
                        c12 = c1[:j] + bj + c1[j + 1:]
                        a1, a2, a12 = CODON[c1], CODON[c2], CODON[c12]
                        out.append(((gps[i], comp(bi), snv_class(a0, a1)), (gps[j], comp(bj), snv_class(a0, a2)),
                                    snv_class(a0, a12), a0, a1, a2, a12))
    return out


def haplotype_cons(variants, models, cds_at, base):
    """Rule 5 of the level-4 SOP as code on planted data. base(v) -> [(gene, class)] of a variant on its own.
    Returns cons(v): the same, but single-nucleotide variants that share a codon on one chromosome copy of
    the proband get the class of the codon with all of these changes applied."""
    by_codon = {}
    for i, v in enumerate(variants):
        if v["pos"] is None or len(v["ref"]) != 1 or len(v["alt"]) != 1 or v["gt"]["P"] == 0:
            continue
        hit = cds_at.get(v["pos"])
        if hit is None:
            continue
        g = v["gt"]
        copy = "b" if g["P"] == 2 else ("m" if g["M"] > 0 and g["F"] == 0 else ("f" if g["F"] > 0 and g["M"] == 0 else "u"))
        by_codon.setdefault((hit[0], hit[1] // 3), []).append((i, hit[1] % 3, copy))
    joint = {}
    for (gname, c), members in by_codon.items():
        if len(members) < 2:
            continue
        gm = models[gname]
        codon = gm.cds[3 * c:3 * c + 3]
        a0 = CODON[codon]

        def cls_with(idxs):
            new, seen = list(codon), set()
            for i, k, _ in members:
                if i in idxs:
                    if k in seen:
                        return None            # two alleles of one position cannot share a copy
                    seen.add(k)
                    alt = variants[i]["alt"]
                    new[k] = alt if gm.strand == "+" else alt.translate(COMP)
            a1 = CODON["".join(new)]
            if a0 == "*":
                return "stop_retained" if a1 == "*" else "stop_lost"
            return snv_class(a0, a1)

        homs = {i for i, _, cp in members if cp == "b"}
        for i, k, cp in members:
            if cp == "b":
                groups = [homs | {j for j, _, c2 in members if c2 == x} for x in "mf"]
                groups += [homs | {j} for j, _, c2 in members if c2 == "u"]
            elif cp == "u":
                groups = [homs | {i}]
            else:
                groups = [homs | {j for j, _, c2 in members if c2 == cp}]
            joint[i] = (gname, {None if len(grp) == 1 else cls_with(grp) for grp in groups})
    index = {id(v): i for i, v in enumerate(variants)}

    def cons(v):
        alone = base(v)
        j = joint.get(index[id(v)])
        if j is None:
            return alone
        gname, classes = j
        own = next((c for g, c in alone if g == gname), "other")
        classes = sorted({own if c is None else c for c in classes})
        rest = [(g, c) for g, c in alone if g != gname]
        return rest + [(gname, classes[0] if len(classes) == 1 else "&".join(classes))]
    return cons


INSTRUCTION = """\
You are the clinical bioinformatician of a rare-disease genetics laboratory. An external sequencing
provider delivered a joint variant call set for a family trio (an affected child and both parents).
The laboratory needs this call set on its own reference and the candidate causal variants for the
gene panel of the referral.

Task directory:
- `/task`

Inputs (read-only) under `/task/input`:
- `trio.vcf.gz` (+ `.tbi`) - the provider's joint call set of the three family members, chromosome 22
  only, FORMAT fields `GT:AD:DP:GQ`
- `pedigree.ped` - six-column PED file (family, individual, father, mother, sex, phenotype; phenotype
  2 = affected)
- `panel_genes.txt` - gene symbols of the referral's panel, one per line
- `lab_sop.md` - the laboratory's SOP for trio candidate review
- `resources/GRCh37.chr22.fa.gz`, `resources/GRCh38.chr22.fa.gz` (+ `.fai`, `.gzi`) - chromosome 22 of
  both human assemblies, bgzip-compressed (sequence outside the regions the laboratory works on is
  masked with N; coordinates are those of the full chromosome)
- `resources/hg19ToHg38.chr22.over.chain`, `resources/hg38ToHg19.chr22.over.chain` - UCSC liftOver
  chain files restricted to chromosome 22
- `resources/GRCh38.chr22.canonical.gff3` - the laboratory's transcript set on GRCh38 (Ensembl GFF3
  conventions, one transcript per gene)
- `resources/popfreq.GRCh38.tsv.gz` - population allele frequencies, columns `chrom`, `pos`, `ref`,
  `alt`, `af` (GRCh38, 1-based, one ALT allele per row, left-aligned minimal representation)

Available software (on `PATH`): `bcftools` (with its plugins, including `+liftover`), `samtools`,
`bgzip`, `tabix`, `picard`, Python 3.12 with pysam, numpy, scipy and pandas. There is no internet
access.

Your task:
1. Deliver the provider's call set on GRCh38 in the laboratory's canonical representation.
2. Report which input alleles cannot be placed on GRCh38.
3. Apply the laboratory SOP and report the candidate causal variants of the affected child.

Required outputs under `/task/output`:
- `trio.GRCh38.vcf.gz` plus its index (`.tbi` or `.csi`) - bgzip-compressed VCF on GRCh38 with the
  contig named `chr22`, sorted by position. One ALT allele per record; every record in left-aligned,
  minimal (parsimonious) representation with `REF` equal to the GRCh38 reference sequence. The same
  three samples under their original names, each with its genotype (`GT`) for every record. It holds
  every allele of the input that has a position on GRCh38 chromosome 22 according to the provided
  chain file, and nothing else.
- `unmapped.tsv` - tab-separated with header `chrom_in`, `pos_in`, `ref`, `alt`: one row per ALT
  allele of the input that has no position on GRCh38 chromosome 22 according to the provided chain
  file, with chromosome, position and alleles as written in the input file. Header only if there is
  none.
- `candidates.tsv` - tab-separated with header `chrom`, `pos`, `ref`, `alt`, `gene`, `consequence`,
  `model`: one row per candidate variant under the SOP. `chrom`, `pos`, `ref`, `alt` exactly as in
  your `trio.GRCh38.vcf.gz`; `gene` is the gene symbol as written in `panel_genes.txt`;
  `consequence` is one of the qualifying classes named in the SOP; `model` is one of `de_novo`,
  `compound_het`, `homozygous_recessive`. Header only if there is no candidate.

Constraints:
- Do not modify anything under `/task/input`; write deliverables only under `/task/output`.
- Do not ask for confirmation. Execute directly.
"""

SOP = """\
# SOP-GEN-017 (rev. 4) - Review of trio call sets for candidate causal variants

Scope: germline joint call sets of an affected child ("proband") and both parents, restricted to the
gene panel of the referral. All interpretation resources of the laboratory (transcript set,
population frequencies) are kept on GRCh38 with UCSC-style contig names, and results are issued on
GRCh38.

1. **Representation.** Variants are compared with laboratory resources and reported as one ALT
   allele per row, left-aligned and minimal.
2. **Usable genotypes.** A genotype is usable if `DP` >= {dp} and `GQ` >= {gq}. A variant is reviewed
   only if the genotypes of all three family members are usable.
3. **Allele balance.** A heterozygous genotype of the proband is accepted only if the fraction of
   reads supporting the ALT allele, `AD_alt / (AD_ref + AD_alt)`, lies between {ablo} and {abhi}
   (limits included).
4. **Population frequency.** A variant is rare if its allele frequency in the laboratory frequency
   table is below {af}. A variant that is absent from the table is treated as having frequency 0.
5. **Consequence.** Consequences are taken relative to the laboratory transcript set, each variant
   on its own (no merging of neighbouring variants into haplotype consequences). Qualifying
   classes: `missense`, `stop_gained`, `frameshift`, `splice_acceptor`, `splice_donor`. Every other
   consequence (synonymous, in-frame insertion/deletion, UTR, intronic, intergenic, ...) does not
   qualify.
6. **Panel.** Only variants with a qualifying consequence in a gene of the referral's panel are
   reviewed.
7. **Inheritance models.** A rare variant that passes rules 2-6 is reported if it fits one of:
   - `de_novo` - proband heterozygous, both parents homozygous reference;
   - `homozygous_recessive` - proband homozygous for the ALT allele, both parents heterozygous;
   - `compound_het` - the proband is heterozygous for two different variants in the same gene, one
     of them inherited from the mother only (mother heterozygous, father homozygous reference) and
     the other from the father only (father heterozygous, mother homozygous reference). Both
     variants are reported.
   Variants that fit none of these models are not reported.
"""

JUDGE = """\
**Where the truth comes from.** Every allele is planted on GRCh38 in left-aligned minimal form with
planted trio genotypes; the provider VCF is derived from it (mapped back to hg19 inside ungapped
plus-strand chain blocks with identical flanks, then renamed/right-shifted/padded/merged). The
candidate list is the planted causal set. `bcftools csq`, Picard LiftoverVcf and bcftools +liftover
are used only as guards at generation time: an allele is *scored* only if both lift-over tools
reproduce the planted key (hg19-only alleles: only if both tools reject them); everything else is a
gray zone that is neither rewarded nor penalised; generation is retried when csq disagrees with a
designed consequence of a causal/decoy variant, or when the SOP applied to csq consequences does not
return exactly the planted causal set.

**Hard gates (score 0).** Missing/unparseable file; VCF without exactly the three input samples;
no `.tbi`/`.csi` index; VCF larger than 10 MB compressed / 50 MB uncompressed; more than 200
candidate rows; `model` outside the three allowed values; candidate alleles that are not plain
nucleotide strings.

**Components.**

| weight | component | metric | full credit | zero |
|---|---|---|---|---|
| 0.15 | lifted_callset | F1 of literal (chrom, pos, ref, alt) keys against the scored planted GRCh38 alleles | >= 0.985 | <= 0.90 |
| 0.10 | genotypes | min over the three samples of the fraction of scored planted alleles present with the planted genotype | >= 0.985 | <= 0.90 |
| 0.10 | unmapped | F1 of (pos_in, ref, alt) against the planted hg19-only alleles | >= 0.95 | <= 0.50 |
| 0.45 | candidates | F1 of candidate keys against the planted causal set | 1.0 | <= 0.70 |
| 0.20 | candidate_annotation | fraction of planted causal variants reported with the right gene, consequence class and model | 1.0 | 0 |

**Why these tolerances.** Both reference solutions reproduce every scored key (F1 = 1.000, genotype
concordance 1.000), so the 0.985 edge only leaves room for a handful of legitimately debatable
records; a call set that was lifted but not brought into canonical representation reaches F1 of
about 0.77 (measured) and therefore gets nothing for that component. The candidate set is tiny (3-4
variants), so one false or one missing candidate already moves F1 to <= 0.89: a clinical short list
with a wrong entry is not a pass. Keys are compared literally because the instruction defines the
representation; chromosome naming is part of the key.

**What is deliberately not scored.** INFO/FORMAT fields other than GT, record order, extra
columns, the gray zone.
"""

INSTRUCTION_L4 = """\
You are the clinical bioinformatician of a rare-disease genetics laboratory. An external sequencing
provider delivered a joint variant call set for a family trio (an affected child and both parents).
The laboratory needs to know how the family differs from the laboratory's own reference assembly,
and the candidate causal variants for the gene panel of the referral.

Task directory:
- `/task`

Inputs (read-only) under `/task/input`:
- `trio.vcf.gz` (+ `.tbi`) - the provider's joint call set of the three family members, chromosome 22
  only, FORMAT fields `GT:AD:DP:GQ`. The provider states that the file is complete with respect to
  the reference sequence the provider worked with: at every position of that reference which is not
  covered by a record of the file, all three family members are homozygous for the base of the
  provider's reference.
- `pedigree.ped` - six-column PED file (family, individual, father, mother, sex, phenotype; phenotype
  2 = affected)
- `panel_genes.txt` - gene symbols of the referral's panel, one per line
- `lab_sop.md` - the laboratory's SOP for trio candidate review
- `resources/GRCh37.chr22.fa.gz`, `resources/GRCh38.chr22.fa.gz` (+ `.fai`, `.gzi`) - chromosome 22 of
  both human assemblies, bgzip-compressed, exactly as they are in use at the provider and at the
  laboratory (sequence outside the regions the laboratory works on is masked with N; coordinates are
  those of the full chromosome)
- `resources/hg19ToHg38.chr22.over.chain`, `resources/hg38ToHg19.chr22.over.chain` - UCSC liftOver
  chain files restricted to chromosome 22
- `resources/GRCh38.chr22.canonical.gff3` - the laboratory's transcript set on GRCh38 (Ensembl GFF3
  conventions, one transcript per gene)
- `resources/popfreq.GRCh38.tsv.gz` - population allele frequencies, columns `chrom`, `pos`, `ref`,
  `alt`, `af` (GRCh38, 1-based, one ALT allele per row, left-aligned minimal representation)

Available software (on `PATH`): `bcftools` (with its plugins, including `+liftover`), `samtools`,
`bgzip`, `tabix`, `picard`, Python 3.12 with pysam, numpy, scipy and pandas. There is no internet
access.

Your task:
1. Deliver the call set of the family on GRCh38 in the laboratory's canonical representation.
2. Report which input alleles lie in sequence that cannot be placed on GRCh38.
3. Apply the laboratory SOP and report the candidate causal variants of the affected child.

Required outputs under `/task/output`:
- `trio.GRCh38.vcf.gz` plus its index (`.tbi` or `.csi`) - bgzip-compressed VCF on GRCh38 with the
  contig named `chr22`, sorted by position. It describes how the three family members differ from
  the GRCh38 reference sequence: it holds every allele that differs from the GRCh38 reference and is
  carried by at least one family member, for all GRCh38 positions that the provided chain file
  aligns to the provider's reference (the aligned blocks of the chain; sequence that the chain does
  not align and N-masked sequence are out of scope), and nothing else. One ALT allele per record;
  every record in left-aligned, minimal (parsimonious) representation with `REF` equal to the GRCh38
  reference sequence; differences of single bases are written as single-base records. The same
  three samples under their original names. `GT` of a sample states how many of its two chromosome
  copies carry the ALT allele of that record (`0/0`, `0/1` or `1/1`); a copy that carries any other
  allele counts as `0` in that record.
- `unmapped.tsv` - tab-separated with header `chrom_in`, `pos_in`, `ref`, `alt`: one row per ALT
  allele of the input that lies at a position of the provider's reference which the provided chain
  file does not align to GRCh38 chromosome 22, with chromosome, position and alleles as written in
  the input file. An input allele at an aligned position is never listed here, whatever it turns
  out to be on GRCh38. Header only if there is none.
- `candidates.tsv` - tab-separated with header `chrom`, `pos`, `ref`, `alt`, `gene`, `consequence`,
  `model`: one row per candidate variant under the SOP. `chrom`, `pos`, `ref`, `alt` exactly as in
  your `trio.GRCh38.vcf.gz`; `gene` is the gene symbol as written in `panel_genes.txt`;
  `consequence` is the qualifying class that the SOP assigns to the variant; `model` is one of
  `de_novo`, `compound_het`, `homozygous_recessive`. Header only if there is no candidate.

Constraints:
- Do not modify anything under `/task/input`; write deliverables only under `/task/output`.
- Do not ask for confirmation. Execute directly.
"""

SOP_L4 = """\
# SOP-GEN-017 (rev. 5) - Review of trio call sets for candidate causal variants

Scope: germline joint call sets of an affected child ("proband") and both parents, restricted to the
gene panel of the referral. All interpretation resources of the laboratory (transcript set,
population frequencies) are kept on GRCh38 with UCSC-style contig names, and results are issued on
GRCh38. Every rule below refers to the call set as expressed against GRCh38: "ALT allele" is an
allele that differs from the GRCh38 reference, and "carries n copies" counts the chromosome copies
of a family member with that allele.

1. **Representation.** Variants are compared with laboratory resources and reported as one ALT
   allele per row, left-aligned and minimal.
2. **Usable genotypes.** A genotype is usable if `DP` >= {dp} and `GQ` >= {gq}. A variant is reviewed
   only if the genotypes of all three family members are usable. Genotypes that follow from a
   provider's statement of completeness (positions without a record in the delivered file) are
   usable.
3. **Allele balance.** If the proband carries one copy of the ALT allele, the variant is accepted
   only if the fraction of reads supporting the ALT allele, `AD_alt / (AD_ref + AD_alt)` with the
   read counts of the ALT allele and of the GRCh38 reference allele, lies between {ablo} and {abhi}
   (limits included).
4. **Population frequency.** A variant is rare if its allele frequency in the laboratory frequency
   table is below {af}. A variant that is absent from the table is treated as having frequency 0.
5. **Consequence.** Consequences are taken relative to the laboratory transcript set and describe
   what the chromosome copy of the proband that carries the variant encodes.
   - A variant is assessed on its own, with one exception: if one chromosome copy of the proband
     carries two or more single-nucleotide variants in the same codon, these are assessed together.
     Each of them then has the consequence of the codon with all of these changes applied, compared
     with the reference codon (`synonymous` if the amino acid is unchanged, `missense` if it is
     changed, `stop_gained` if a stop codon results). Single-nucleotide variants in one codon that
     lie on different copies are assessed separately, and so are all other combinations of variants.
   - Which copy carries a variant follows from the trio genotypes. A variant of which the proband
     carries two copies is on both copies. A variant of which the proband carries one copy is on
     the copy inherited from the mother if the mother carries the allele and the father does not,
     and on the copy inherited from the father if the father carries it and the mother does not. In
     all other cases (neither parent or both parents carry the allele) the copy cannot be assigned;
     such a variant is assessed together only with variants of which the proband carries two copies.
   - All variants of the call set that the proband carries count for this, whatever their
     frequency, genotype quality or own consequence.
   - A variant of which the proband carries two copies qualifies only if its consequence qualifies
     on both copies (if the two classes differ, report both joined by `&`).

   Qualifying classes: `missense`, `stop_gained`, `frameshift`, `splice_acceptor`, `splice_donor`.
   Every other consequence (synonymous, in-frame insertion/deletion, UTR, intronic, intergenic, ...)
   does not qualify.
6. **Panel.** Only variants with a qualifying consequence in a gene of the referral's panel are
   reviewed.
7. **Inheritance models.** A rare variant that passes rules 2-6 is reported if it fits one of:
   - `de_novo` - the proband carries one copy of the ALT allele, neither parent carries it;
   - `homozygous_recessive` - the proband carries two copies, each parent carries exactly one;
   - `compound_het` - the proband carries one copy each of two different variants in the same gene,
     one of them inherited from the mother only (the mother carries one copy, the father none) and
     the other from the father only (the father carries one copy, the mother none). Both variants
     are reported.
   Variants that fit none of these models are not reported.
"""

JUDGE_L4 = """\
**Where the truth comes from.** The family's alleles are planted on GRCh38 (left-aligned minimal,
planted number of ALT copies per member); the provider VCF is *derived* by expressing these
haplotypes against the provider's reference. For alleles at positions where both assemblies agree
this is the level-3 construction (mapped back inside ungapped chain blocks with identical flanks,
then renamed/right-shifted/padded/merged; scored only if Picard LiftoverVcf and bcftools +liftover
both reproduce the planted key and genotypes). New at level 4:

- *Reference differences.* Positions inside aligned blocks where the provider's reference base
  differs from GRCh38: the 173 differences that exist between the two real assemblies inside the
  shipped windows plus five planted into the shipped GRCh37 slice inside coding sequence. There the
  truth needs no tool: the planted allele pair of each member is written once against GRCh38
  (truth) and once against the provider's base (input). A site where all three members are
  homozygous for the provider's base has no input record and is a homozygous ALT allele on GRCh38;
  a provider record whose ALT is the GRCh38 base turns into the opposite allele with complementary
  genotypes; a site where everybody carries the GRCh38 base disappears; some sites carry a third
  allele.
- *Inverted block.* A 3-7 kb segment of the shipped GRCh37 slice is replaced by its reverse
  complement and both chain files carry it as a separate minus-strand chain. Alleles inside it are
  written on the other strand and left-aligned there; they are scored only if both guard tools
  bring them back to the planted GRCh38 allele.
- *Haplotype consequence.* The generator runs the SOP, including the per-copy codon rule, as code
  on the planted data and requires exactly the planted causal set; `bcftools csq -l` confirms every
  designed single-variant class.

**Hard gates (score 0).** Missing/unparseable file; VCF without exactly the three input samples;
no `.tbi`/`.csi` index; VCF larger than 10 MB compressed / 50 MB uncompressed; more than 200
candidate rows; `model` outside the three allowed values; candidate alleles that are not plain
nucleotide strings.

**Components.** In the three call-set components a true positive needs the literal (chrom, pos,
ref, alt) key *and* the planted number of ALT copies in all three samples; a submitted record that
is not a scored allele is a false positive of the category its position belongs to.

| weight | component | metric | full credit | zero |
|---|---|---|---|---|
| 0.08 | callset_standard | F1 over scored alleles at plus-strand positions where both assemblies agree | >= 0.985 | <= 0.90 |
| 0.12 | callset_reference_differences | F1 over alleles at positions where the provider's reference base differs from GRCh38 | >= 0.97 | <= 0.80 |
| 0.07 | callset_inverted_block | F1 over scored alleles inside the minus-strand block | >= 0.95 | <= 0.70 |
| 0.08 | unmapped | F1 of (pos_in, ref, alt) against the planted alleles in unaligned provider sequence | >= 0.95 | <= 0.50 |
| 0.45 | candidates | F1 of candidate keys against the planted causal set (13 variants) | 1.0 | <= 0.87 |
| 0.20 | candidate_annotation | fraction of planted causal variants reported with the right gene, consequence class and model | 1.0 | 0 |

Records in which no sample carries the ALT allele are ignored (they claim nothing). Keys are compared
literally because the instruction defines the representation.

**Why these tolerances** (measured on seeds 1 and 2). Both reference solutions reproduce every
scored allele of all three categories (F1 = 1.000) and the exact candidate list. The standard band
is the level-3 band. Reference differences: a solution that re-expresses the provider's records but
adds no site reaches F1 0.64-0.65, one that splits multi-allelic records before the lift-over 0.91,
one that rejects the records 0; the band 0.80-0.97 gives nothing to the first, a third of the credit
to the second and full credit only above both. Inverted block: handled gives 1.000, not handled 0;
the band only has to leave room for one or two of about 45 alleles. With 13 planted candidates one
false row gives F1 0.963 and one missing row 0.960; the 0.87-1.0 ramp turns both into a total below
the pass mark (0.87 / 0.85), so a clinical short list with a wrong entry still does not pass, while
two errors keep partial credit and the annotation component is linear in the number of correct rows.

**What is deliberately not scored.** INFO/FORMAT fields other than GT, record order, extra columns,
records without a carrier, and the gray zone: alleles in 1-60 bp chain gaps on which the two guard
tools disagree, and assembly differences at neighbouring bases (one record or two is a matter of
convention).
"""


def generate(seed: int, level: int, workdir: Path) -> TaskBuild:
    last = None
    for attempt in range(6):
        try:
            return _generate(seed, level, Path(workdir), attempt)
        except _Retry as e:          # a guard tripped: redraw with the next sub-seed
            last = e
    raise RuntimeError(f"trio_liftover: generation guards failed 6 times: {last}")


class _Retry(Exception):
    pass


def _generate(seed: int, level: int, workdir: Path, attempt: int) -> TaskBuild:
    R = _ref()
    rng = np.random.default_rng([4242, seed, level, attempt])
    s38, s19 = R.hg38, R.hg19
    L38 = len(s38)
    nonnorm = level >= 2
    decoys3 = level >= 3
    lvl4 = level >= 4

    models = {g["name"]: GeneModel(g, s38) for g in R.genes.values()}
    valid = sorted(n for n, m in models.items() if m.valid)
    near_exon = np.zeros(L38, dtype=bool)
    gene_cov = np.zeros(L38, dtype=np.int8)
    for m in models.values():
        gene_cov[max(0, m.span[0] - 200):m.span[1] + 200] += 1
        for a, b in m.exons:
            near_exon[max(0, a - 30):b + 30] = True
    occ = np.zeros(L38, dtype=bool)
    rd_zone = np.zeros(L38, dtype=bool)   # level 4: neighbourhood of sites where the assemblies differ
    inv_zone = np.zeros(L38, dtype=bool)  # level 4: GRCh38 interval that is inverted in the provider's reference
    inv = None
    variants: list[dict] = []

    def placeable(pos, ref, alt, margin=12):
        k = shift_max(s38, pos, ref, alt) if len(ref) != len(alt) else 0
        if k > 20:
            return None
        lo, hi = pos - margin, pos + len(ref) + k + margin
        if occ[lo:hi].any() or rd_zone[lo:hi].any():
            return None
        m = R.to19(pos)
        if m is None or m[1] < 100 + (hi - lo):
            return None
        if not R.flank_identical(pos, m[0], 50, len(ref) + k + 50):
            return None
        return m[0], k, lo, hi

    def forceable(pos, ref):
        """Can an allele be planted here although the spot is occupied by its record partner?"""
        m = R.to19(pos)
        return m is not None and m[1] >= 150 and R.flank_identical(pos, m[0], 50, len(ref) + 50) and not inv_zone[pos]

    def add(pos, ref, alt, kind, gt, gene=None, cls=None, af=None, fmt=None, force=False, clean=False):
        """Plant one normalised GRCh38 allele. Returns the variant dict or None if it cannot be placed."""
        if normalise(s38, pos, ref, alt) != (pos, ref, alt):
            return None
        pl = placeable(pos, ref, alt)
        if pl is None and not force:
            return None
        if pl is None:
            if not forceable(pos, ref):
                return None
            pl = (R.to19(pos)[0], 0, pos, pos)
        p19, k, lo, hi = pl
        occ[lo:hi] = True
        v = dict(pos=pos, ref=ref, alt=alt, kind=kind, gene=gene, cls=cls, af=af, pos19=p19, kmax=k,
                 gt=dict(zip("PMF", gt)), rep="norm", group=None)
        v["fmt"] = fmt or {r: fmt_good(rng, v["gt"][r], clean) for r in "PMF"}
        if lvl4 and inv_zone[pos]:
            v["inv"] = True
        variants.append(v)
        return v

    def rand_noncoding():
        while True:
            a, b = R.windows38[int(rng.integers(len(R.windows38)))]
            p = int(rng.integers(a + 2000, b - 2000))
            if not near_exon[p - 40:p + 40].any() and "N" not in s38[p - 40:p + 40]:
                return p

    def rand_indel(p):
        n = int(min(rng.geometric(0.45), 8))
        if rng.random() < 0.55:
            return normalise(s38, p, s38[p:p + 1 + n], s38[p])
        ins = "".join(BASES[int(x)] for x in rng.integers(0, 4, size=n)) if rng.random() < 0.5 else s38[p + 1:p + 1 + n]
        return normalise(s38, p, s38[p], s38[p] + ins)

    def trio_from_af(af):
        while True:
            m, f = int(rng.binomial(2, af)), int(rng.binomial(2, af))
            if m + f > 0:
                break
        c = int(rng.random() < m / 2) + int(rng.random() < f / 2)
        return c, m, f

    def rare_af():
        return None if rng.random() < 0.5 else float(f"{10 ** rng.uniform(-5.0, -3.25):.3g}")

    def common_af():
        return float(f"{10 ** rng.uniform(-2.2, -0.25):.3g}")

    def inherited():
        return (1, 1, 0) if rng.random() < 0.5 else (1, 0, 1)

    def bg_af_l4():
        return (None if rng.random() < 0.25 else float(f"{10 ** rng.uniform(-4.6, -3.2):.3g}")) \
            if rng.random() < 0.18 else common_af()

    def damaging(gm, kinds=("missense", "missense", "stop_gained", "frameshift", "splice")):
        k = kinds[int(rng.integers(len(kinds)))]
        if k == "frameshift":
            d = design_frameshift(gm, s38, rng)
            return None if d is None else (*d, "frameshift")
        if k == "splice":
            return design_splice(gm, s38, rng)
        d = design_snv(gm, s38, rng, k)
        return None if d is None else (*d, k)

    def must(v, what):
        if v is None:
            raise _Retry(f"could not place {what}")
        return v

    def plant_coding(gname, kind, gt, kinds=None, af=None, fmt=None, clean=True, tries=60, **kw):
        gm = models[gname]
        for _ in range(tries):
            d = damaging(gm, kinds) if kinds else damaging(gm)
            if d is None:
                continue
            pos, ref, alt, cls = d
            if gene_cov[pos] != 1:
                continue
            v = add(pos, ref, alt, kind, gt, gene=gname, cls=cls, af=af, fmt=fmt, clean=clean, **kw)
            if v is not None:
                return v
        return None

    # ---- roles of genes (all seed dependent) -------------------------------------------------
    order = [valid[i] for i in rng.permutation(len(valid))]
    n_panel = int(rng.integers(20, 27))
    if lvl4:
        n_panel += 6                      # more genes carry a planted causal/decoy structure
    panel = sorted(order[:n_panel])
    pool = list(order[:n_panel])          # panel genes that do not yet carry a causal/decoy variant
    np_pool = list(order[n_panel:])       # the same for non-panel genes
    with_hr = bool(rng.random() < 0.6) or lvl4
    multi_causal = "dn" if rng.random() < 0.5 else "ch_mis"

    def host(what, fn, genes=None):
        """Plant something in the first unused gene that can take it (fn is all-or-nothing)."""
        genes = pool if genes is None else genes
        for gname in list(genes):
            v = fn(gname)
            if v is not None:
                genes.remove(gname)
                return gname, v
        raise _Retry(f"could not place {what}")

    def shiftable_frameshift(gname, kind, gt, af):
        for _ in range(80):
            d = design_frameshift(models[gname], s38, rng, need_shift=True)
            if d and gene_cov[d[0]] == 1:
                v = add(*d, kind, gt, gene=gname, cls="frameshift", af=af, clean=True)
                if v:
                    return v
        return None

    def causal_in_multiallelic(gname, kind, gt):
        """Rare missense SNV that shares a record with a common in-frame deletion (deletion = 1st ALT)."""
        gm = models[gname]
        for _ in range(300):
            d = design_frameshift(gm, s38, rng, inframe=True)
            if d is None:
                continue
            dpos, dref, dalt = d
            if shift_max(s38, dpos, dref, dalt) != 0 or gene_cov[dpos] != 1:
                continue
            j = int(rng.integers(1, len(dref)))
            p = dpos + j
            ci = gm.gpos.index(p) // 3
            codon = gm.cds[3 * ci:3 * ci + 3]
            k = gm.gpos.index(p) % 3
            alts = []
            for b in BASES:
                if b != codon[k]:
                    naa = CODON[codon[:k] + b + codon[k + 1:]]
                    if naa != "*" and naa != CODON[codon]:
                        alts.append(b if gm.strand == "+" else b.translate(COMP))
            if not alts or placeable(dpos, dref, dalt) is None or not forceable(p, s38[p]):
                continue
            carrier = (0, 1, 0) if rng.random() < 0.5 else (0, 0, 1)
            if gt == (1, 1, 0):
                carrier = (0, 0, 1)
            elif gt == (1, 0, 1):
                carrier = (0, 1, 0)
            a = add(dpos, dref, dalt, "common_inframe_del", carrier, gene=gname, cls="inframe_deletion",
                    af=common_af())
            b = add(p, s38[p], alts[int(rng.integers(len(alts)))], kind, gt, gene=gname, cls="missense",
                    af=None, force=True, clean=True)
            if b is None:
                raise _Retry("multi-allelic causal SNV not placeable")
            a["group"] = b["group"] = ("delsnv", len(variants))
            # what a splitter makes of the shared record: one DP/GQ/ref-depth per sample
            for r in "PMF":
                src = b if b["gt"][r] else a
                ref_d, _, dp, gq = src["fmt"][r]
                a["fmt"][r] = (ref_d, a["fmt"][r][1] if a["gt"][r] else 0, dp, gq)
                b["fmt"][r] = (ref_d, b["fmt"][r][1] if b["gt"][r] else 0, dp, gq)
            return b
        return None

    # ---- level 4: one aligned block is inverted between the provider's reference and GRCh38 --------
    if lvl4:
        mm_inv = [m[0] for m in real_mismatches(R)]

        def choose_inversion(gname):
            gm = models[gname]
            for _ in range(40):
                a, b = gm.segs[int(rng.integers(len(gm.segs)))]
                if b - a < 150 or (gene_cov[a:b] != 1).any():
                    continue
                q_a, q_b = a - int(rng.integers(1500, 3000)), b + int(rng.integers(1500, 3000))
                ma, mb = R.to19(q_a), R.to19(q_b - 1)
                if ma is None or mb is None or mb[0] - ma[0] != q_b - 1 - q_a or min(ma[1], mb[1]) < 500:
                    continue
                t_a, t_b = ma[0], mb[0] + 1
                if s19[t_a:t_b] != s38[q_a:q_b] or "N" in s38[q_a:q_b] or any(q_a - 200 <= m < q_b + 200 for m in mm_inv):
                    continue
                fwd = invert_in_chain(R.fwd_chain_text, t_a, t_b, q_a, refdata.HG19_LEN, refdata.HG38_LEN, "900001")
                rev = invert_in_chain(R.rev_chain_text, q_a, q_b, t_a, refdata.HG38_LEN, refdata.HG19_LEN, "900002")
                if fwd is None or rev is None:
                    continue
                return dict(q_a=q_a, q_b=q_b, t_a=t_a, t_b=t_b, fwd=fwd, rev=rev, gene=gname, exon=(a, b))
            return None

        def in_exon(gname, kind, gt, exon, frameshift):
            gm = models[gname]
            for _ in range(800):
                d = design_frameshift(gm, s38, rng) if frameshift else design_snv(gm, s38, rng, "missense")
                if d is None or not exon[0] + 10 <= d[0] < exon[1] - 12 or gene_cov[d[0]] != 1:
                    continue
                v = add(d[0], d[1], d[2], kind, gt, gene=gname, cls="frameshift" if frameshift else "missense",
                        af=rare_af(), clean=True)
                if v is not None:
                    return v
            return None

        for gname in list(pool):
            inv = choose_inversion(gname)
            if inv is not None:
                pool.remove(gname)
                break
        if inv is None:
            raise _Retry("no panel gene with an exon that can sit in an inverted block")
        inv_zone[inv["q_a"]:inv["q_b"]] = True
        for edge in (inv["q_a"], inv["q_b"]):       # nothing straddles or touches the ends of the block
            occ[edge - 80:edge + 80] = True
        par = inherited()
        must(in_exon(inv["gene"], "causal_inv_ch", par, inv["exon"], True), "frameshift in the inverted block")
        must(in_exon(inv["gene"], "causal_inv_ch", (1, 0, 1) if par == (1, 1, 0) else (1, 1, 0), inv["exon"], False),
             "missense in the inverted block")

    # ---- causal set ---------------------------------------------------------------------------
    causal = []
    if multi_causal == "dn":
        _, v = host("de novo", lambda g: causal_in_multiallelic(g, "causal_dn", (1, 0, 0)))
    else:
        _, v = host("de novo", lambda g: plant_coding(g, "causal_dn", (1, 0, 0)))
    causal.append(v)
    lof_parent = (1, 1, 0) if rng.random() < 0.5 else (1, 0, 1)
    mis_parent = (1, 0, 1) if lof_parent == (1, 1, 0) else (1, 1, 0)
    # LoF partner: a frameshift that can be right-shifted
    g_ch, v = host("compound-het LoF", lambda g: shiftable_frameshift(g, "causal_ch", lof_parent, rare_af()))
    causal.append(v)
    if multi_causal == "ch_mis":
        causal.append(must(causal_in_multiallelic(g_ch, "causal_ch", mis_parent), "compound-het missense"))
    else:
        causal.append(must(plant_coding(g_ch, "causal_ch", mis_parent, kinds=("missense",), af=rare_af()),
                           "compound-het missense"))
    if with_hr:
        _, v = host("hom-rec", lambda g: plant_coding(g, "causal_hr", (2, 1, 1), af=rare_af()))
        causal.append(v)

    # ---- decoys: each fails exactly one SOP rule ---------------------------------------------
    host("synonymous de novo", lambda g: plant_coding(g, "decoy_syn_dn", (1, 0, 0), kinds=("synonymous",)))
    host("non-panel de novo", lambda g: plant_coding(g, "decoy_nonpanel_dn", (1, 0, 0)), np_pool)
    for gname in np_pool[:3]:                # true-looking recessive variants outside the panel
        plant_coding(gname, "decoy_nonpanel", (2, 1, 1), af=rare_af())
    np_rest = np_pool[3:]
    if nonnorm:
        # common frameshift in a recessive pattern, delivered right-shifted
        host("common shifted frameshift", lambda g: shiftable_frameshift(
            g, "decoy_common_shifted", (2, 1, 1), float(f"{rng.uniform(0.012, 0.06):.3g}")))

        def wrong_build(gname):
            """Record whose hg19 coordinate, read as GRCh38, is a missense position of this panel gene
            (same reference base in both assemblies); its true GRCh38 position is non-coding."""
            gm = models[gname]
            for _ in range(1500):
                d = design_snv(gm, s38, rng, "missense")
                if d is None:
                    continue
                e, _, alt = d
                m = R.to38(e)
                if (s19[e] != s38[e] or gene_cov[e] != 1 or m is None or m[1] < 200 or inv_zone[m[0]]
                        or near_exon[m[0] - 40:m[0] + 40].any() or alt == s38[m[0]]):
                    continue
                wb = add(m[0], s38[m[0]], alt, "decoy_wrong_build", (1, 0, 0), clean=True)
                if wb:
                    wb["wb_gene"] = gname
                    return wb
            return None

        host("wrong-build decoy", wrong_build)
    if decoys3:
        for _ in range(2):                    # two rare damaging hets from the SAME parent
            par = inherited()
            g, _ = host("same-parent pair", lambda g: plant_coding(g, "decoy_same_parent", par, af=rare_af()))
            must(plant_coding(g, "decoy_same_parent", par, kinds=("missense",), af=rare_af()), "same-parent partner")
        for i in range(2):                    # "de novo" with an unusable parental genotype
            _, v = host("low-quality parent", lambda g: plant_coding(g, "decoy_parent_quality", (1, 0, 0)))
            who = "M" if rng.random() < 0.5 else "F"
            if i == 0:
                v["fmt"][who] = fmt_lowdp(rng, 0)
            else:
                dp = int(rng.integers(12, 19))
                v["fmt"][who] = (dp - 2, 2, dp, int(rng.integers(5, 16)))
        for _ in range(2):                    # "de novo" with skewed allele balance in the proband
            _, v = host("allele-balance decoy", lambda g: plant_coding(g, "decoy_allele_balance", (1, 0, 0)))
            dp = int(rng.integers(40, 64))
            alt = max(3, int(round(dp * rng.uniform(0.09, 0.19))))
            v["fmt"]["P"] = (dp - alt, alt, dp, int(rng.integers(30, 80)))

    # ---- level 4: sites where the two assemblies differ, and codon pairs --------------------------
    rd_sites, rd_planted, cds_at = [], {}, {}
    if lvl4:
        for gname in valid:
            for i, p in enumerate(models[gname].gpos):
                cds_at[p] = None if p in cds_at else (gname, i)      # positions coding in two genes: unused
        cds_at = {p: h for p, h in cds_at.items() if h is not None}

        def site_depths(bases, clean=True):
            """Read depths of one site: {role: ({base: reads}, DP, GQ)} for true allele pairs."""
            out = {}
            for r in "PMF":
                x, y = bases[r]
                dp = int(np.clip(rng.poisson(34), 16, 70))
                if x == y:
                    out[r] = ({x: dp}, dp, min(99, 3 * dp))
                else:
                    while True:
                        n = int(rng.binomial(dp, 0.5))
                        if 0.34 <= n / dp <= 0.66:
                            break
                    out[r] = ({x: dp - n, y: n}, dp, 99)
            return out

        def rd_site(pos, b19, bases, kind, gene=None, classes=None, afs=None, plant=True):
            """A GRCh38 position whose base differs from the provider's reference base b19.
            bases = {role: (allele, allele)} are the true alleles of the family; the GRCh38 truth alleles
            are added to `variants`, the provider record is derived later. All-or-nothing when planting."""
            b38 = s38[pos]
            if plant:
                if b19 == b38 or inv_zone[pos] or placeable(pos, b38, b19) is None or occ[pos - 60:pos + 61].any():
                    return None
                rd_planted[pos] = b19
            rd_zone[pos - 60:pos + 61] = True
            occ[pos - 12:pos + 13] = True
            carried = {b for r in "PMF" for b in bases[r]}
            site = dict(pos=pos, pos19=R.to19(pos)[0], b38=b38, b19=b19, bases=bases, kind=kind,
                        depth=site_depths(bases) if carried != {b19} else None, variants=[])
            for alt in sorted(carried - {b38}):
                gt = {r: sum(1 for b in bases[r] if b == alt) for r in "PMF"}
                if site["depth"] is None:        # no provider record: usable by definition (SOP rule 2)
                    fmt = {r: (0, 40, 40, 99) for r in "PMF"}
                else:
                    fmt = {r: (site["depth"][r][0].get(b38, 0), site["depth"][r][0].get(alt, 0),
                               site["depth"][r][1], site["depth"][r][2]) for r in "PMF"}
                v = dict(pos=pos, ref=b38, alt=alt, kind=kind, gene=gene, cls=(classes or {}).get(alt),
                         af=(afs or {}).get(alt), pos19=site["pos19"], kmax=0, gt=gt, rep="norm", group=None,
                         fmt=fmt, rd=site)
                variants.append(v)
                site["variants"].append(v)
            rd_sites.append(site)
            return site

        def high_af():
            return float(f"{rng.uniform(0.35, 0.99):.3g}")

        def swap_site(gname, kind, gts):
            """The provider's reference base is a rare damaging allele of this gene; gts = copies of it."""
            gm = models[gname]
            for _ in range(60):
                cls = "stop_gained" if rng.random() < 0.2 else "missense"
                d = design_snv(gm, s38, rng, cls)
                if d is None or gene_cov[d[0]] != 1:
                    continue
                pos, ref, alt = d
                bases = {r: tuple([alt] * g + [ref] * (2 - g)) for r, g in zip("PMF", gts)}
                site = rd_site(pos, alt, bases, kind, gene=gname, classes={alt: cls}, afs={alt: rare_af()})
                if site is not None:
                    return site["variants"][0]
            return None

        def pick_pairs(gname, want):
            opts = [o for o in codon_pairs(models[gname], lambda q: gene_cov[q] == 1) if want(o)]
            return [opts[int(k)] for k in rng.permutation(len(opts))[:60]]

        def trans_pair(gname):
            """Two missense SNVs in one codon, one from each parent (together they would be silent)."""
            for (p1, a1, c1), (p2, a2, c2), joint, *_ in pick_pairs(
                    gname, lambda o: o[0][2] == "missense" and o[1][2] == "missense" and o[2] == "synonymous"):
                if placeable(p1, s38[p1], a1) is None or not forceable(p2, s38[p2]):
                    continue
                x = add(p1, s38[p1], a1, "causal_trans", (1, 1, 0), gene=gname, cls="missense", af=rare_af(), clean=True)
                y = add(p2, s38[p2], a2, "causal_trans", (1, 0, 1), gene=gname, cls="missense", af=rare_af(),
                        force=True, clean=True)
                if x is None or y is None:
                    raise _Retry("trans codon pair not placeable")
                return x
            return None

        def cis_rescue(gname):
            """Two qualifying SNVs in one codon on the same parental copy that together are silent, plus a rare
            damaging variant from the other parent: not a compound heterozygote."""
            par = inherited()
            other = (1, 0, 1) if par == (1, 1, 0) else (1, 1, 0)
            for (p1, a1, c1), (p2, a2, c2), joint, *_ in pick_pairs(
                    gname, lambda o: o[0][2] in QUALIFYING and o[1][2] in QUALIFYING and o[2] == "synonymous"):
                if placeable(p1, s38[p1], a1) is None or not forceable(p2, s38[p2]):
                    continue
                z = plant_coding(gname, "decoy_cis_rescue", other, af=rare_af())
                if z is None:
                    return None
                x = add(p1, s38[p1], a1, "decoy_cis_rescue", par, gene=gname, cls=c1, af=rare_af(), clean=True)
                y = add(p2, s38[p2], a2, "decoy_cis_rescue", par, gene=gname, cls=c2, af=rare_af(), force=True, clean=True)
                if x is None or y is None:
                    raise _Retry("cis codon pair not placeable")
                return x
            return None

        def hidden_pair(gname, kind, want):
            """De novo SNV X next to a base H of the provider's reference that differs from GRCh38 in the same
            codon; nobody in the family differs from the provider's reference at H (no record)."""
            for (px, ax, cx), (ph, ah, ch), joint, *_ in pick_pairs(gname, want):
                if not forceable(px, s38[px]):
                    continue
                site = rd_site(ph, ah, {r: (ah, ah) for r in "PMF"}, "rd_hidden_coding", gene=gname,
                               classes={ah: ch}, afs={ah: high_af()})
                if site is None:
                    continue
                x = add(px, s38[px], ax, kind, (1, 0, 0), gene=gname, cls=cx, force=True, clean=True)
                if x is None:
                    raise _Retry("de novo next to a hidden assembly difference not placeable")
                x["joint"] = joint
                return x
            return None

        host("trans codon pair", trans_pair)
        host("cis codon pair", cis_rescue)
        # X alone silent; with H: a stop codon or an amino acid that neither the reference nor H alone encodes
        host("hidden-difference de novo", lambda g: hidden_pair(
            g, "causal_hidden_dn", lambda o: o[0][2] == "synonymous" and o[1][2] == "missense"
            and o[2] in ("stop_gained", "missense") and o[6] not in (o[3], o[5])))
        # X alone qualifying; with H the codon encodes the reference amino acid again
        host("hidden-difference rescue", lambda g: hidden_pair(
            g, "decoy_hidden_rescue", lambda o: o[0][2] in QUALIFYING and o[1][2] == "missense" and o[2] == "synonymous"))
        host("swap hom-rec", lambda g: swap_site(g, "causal_hr_swap", (2, 1, 1)))
        z_par = inherited()

        def swap_comphet(gname):
            w = swap_site(gname, "causal_ch_swap", (1, 0, 1) if z_par == (1, 1, 0) else (1, 1, 0))
            if w is None:
                return None
            if plant_coding(gname, "causal_ch_swap", z_par, af=rare_af()) is None:
                raise _Retry("partner of the swap compound-het not placeable")
            return w

        host("swap compound het", swap_comphet)

        def tri_denovo(gname):
            """De novo allele at a position where the provider's reference base differs from GRCh38 and
            nobody carries the provider's base: the provider record reads REF=provider base, ALT=GRCh38 base
            and the new allele, parents 1/1, child 1/2. The provider's base itself would be a rare qualifying
            allele of the gene, so a genotype that wrongly keeps a copy of it shows up as a second de novo."""
            gm = models[gname]
            for _ in range(80):
                cls = "stop_gained" if rng.random() < 0.2 else "missense"
                d = design_snv(gm, s38, rng, cls)
                if d is None or gene_cov[d[0]] != 1:
                    continue
                pos, ref, alt = d
                hit = cds_at.get(pos)
                if hit is None:
                    continue
                codon, k = gm.cds[hit[1] - hit[1] % 3:hit[1] - hit[1] % 3 + 3], hit[1] % 3
                to_tx = (lambda b: b) if gm.strand == "+" else (lambda b: b.translate(COMP))
                others = [b for b in BASES if b not in (ref, alt)
                          and snv_class(CODON[codon], CODON[codon[:k] + to_tx(b) + codon[k + 1:]]) in QUALIFYING]
                if not others:
                    continue
                b19 = others[int(rng.integers(len(others)))]
                site = rd_site(pos, b19, {"P": (ref, alt), "M": (ref, ref), "F": (ref, ref)}, "causal_tri_dn",
                               gene=gname, classes={alt: cls})
                if site is not None:
                    return site["variants"][0]
            return None

        host("three-allele de novo", tri_denovo)

        # assembly differences that exist in the real sequences: most stay without a provider record
        mm = real_mismatches(R)
        mm_pos = [m[0] for m in mm]
        for q, t, edge in mm:
            crowded = any(abs(q - o) <= 12 for o in mm_pos if o != q)
            b38, b19 = s38[q], s19[t]
            u = rng.random()
            third = [b for b in BASES if b not in (b38, b19)][int(rng.integers(2))]
            # near exons every such site stays without a record: all three homozygous, so no SOP rule can
            # depend on how reads of a third allele are counted
            if crowded or edge < 100 or near_exon[q] or u < 0.5:
                rd_site(q, b19, {r: (b19, b19) for r in "PMF"}, "rd_hidden", afs={b19: high_af() if rng.random() < 0.75 else None},
                        plant=False)
                continue
            if u < 0.58:                      # everybody carries the GRCh38 base: nothing left on GRCh38
                rd_site(q, b19, {r: (b38, b38) for r in "PMF"}, "rd_fixed38", plant=False)
                continue
            f = float(rng.uniform(0.25, 0.85))
            while True:
                par = {r: tuple(b19 if rng.random() < f else b38 for _ in range(2)) for r in "MF"}
                if u >= 0.88:                 # a third allele segregates in the family
                    who = "M" if rng.random() < 0.5 else "F"
                    par[who] = (third, par[who][1])
                child = (par["M"][int(rng.integers(2))], par["F"][int(rng.integers(2))])
                bases = {"P": child, "M": par["M"], "F": par["F"]}
                carried = {b for r in "PMF" for b in bases[r]}
                if carried not in ({b19}, {b38}):
                    break
            rd_site(q, b19, bases, "rd_tri" if u >= 0.88 else "rd_swap",
                    afs={b19: float(f"{rng.uniform(0.05, 0.95):.3g}"), third: bg_af_l4()}, plant=False)
    spare_panel = pool                        # panel genes without a planted causal/decoy variant

    # ---- background ---------------------------------------------------------------------------
    for gname in spare_panel:                 # a single inherited rare damaging het per spare panel gene
        if rng.random() < 0.75:
            plant_coding(gname, "bg_single_het", inherited(), af=rare_af(), clean=False, tries=20)
    for gname in np_rest:
        for _ in range(int(rng.integers(0, 3))):
            plant_coding(gname, "bg_nonpanel", inherited(), af=rare_af(), clean=False, tries=10)
    for gname in valid:                       # common coding variation and rare synonymous variation
        for _ in range(int(rng.integers(1, 4))):
            af = common_af()
            plant_coding(gname, "bg_common_coding", trio_from_af(af), kinds=("missense", "synonymous", "synonymous"),
                         af=af, clean=False, tries=10)
        if rng.random() < 0.4:
            plant_coding(gname, "bg_rare_syn", inherited(), kinds=("synonymous",), af=rare_af(), clean=False, tries=10)

    n_snv, n_indel, n_multi, n_delsnv = (int(rng.integers(1750, 1950)), int(rng.integers(330, 380)),
                                         int(rng.integers(120, 140)), int(rng.integers(100, 120)))

    def bg_af():
        return (None if rng.random() < 0.25 else float(f"{10 ** rng.uniform(-4.6, -3.2):.3g}")) \
            if rng.random() < 0.18 else common_af()

    def bg_gt(af):
        return trio_from_af(af) if af is not None and af >= 0.005 else ((0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1))[
            int(rng.integers(4))]

    def noisy(v):
        for r in "PMF":
            if rng.random() < 0.02:
                v["fmt"][r] = fmt_lowdp(rng, v["gt"][r])

    done = 0
    while done < n_snv:
        p = rand_noncoding()
        af = bg_af()
        v = add(p, s38[p], [b for b in BASES if b != s38[p]][int(rng.integers(3))], "bg_noncoding", bg_gt(af), af=af)
        if v:
            noisy(v)
            done += 1
    done = 0
    while done < n_indel:
        af = bg_af()
        v = add(*rand_indel(rand_noncoding()), "bg_noncoding", bg_gt(af), af=af)
        if v:
            noisy(v)
            done += 1
    carriers = [((0, 1, 0), (1, 0, 1)), ((0, 0, 1), (1, 1, 0)), ((1, 1, 0), (0, 0, 1)), ((1, 0, 1), (0, 1, 0)),
                ((0, 1, 0), (0, 0, 1)), ((2, 1, 1), (0, 0, 0))]

    def fuse(a, b, tag):
        a["group"] = b["group"] = (tag, len(variants))
        for r in "PMF":
            src = b if b["gt"][r] else a
            ref_d, _, dp, gq = src["fmt"][r]
            a["fmt"][r] = (ref_d, a["fmt"][r][1] if a["gt"][r] else 0, dp, gq)
            b["fmt"][r] = (ref_d, b["fmt"][r][1] if b["gt"][r] else 0, dp, gq)

    done = 0
    while done < n_multi:                    # two SNVs at one position, never both in one sample
        p = rand_noncoding()
        alts = [b for b in BASES if b != s38[p]]
        i, j = (int(x) for x in rng.permutation(3)[:2])
        ga, gb = carriers[int(rng.integers(5))]
        if not forceable(p, s38[p]):
            continue
        a = add(p, s38[p], alts[i], "bg_noncoding", ga, af=common_af())
        if a is None:
            continue
        b = add(p, s38[p], alts[j], "bg_noncoding", gb, af=bg_af(), force=True)
        if b is None:
            raise _Retry("second SNV of a two-SNV site not placeable")
        fuse(a, b, "multisnv")
        done += 1
    done = 0
    while done < n_delsnv:                   # deletion + SNV inside the deleted span
        p = rand_noncoding()
        n = int(rng.integers(2, 6))
        d = normalise(s38, p, s38[p:p + 1 + n], s38[p])
        if shift_max(s38, *d) != 0 or len(d[1]) < 3:
            continue
        ga, gb = carriers[int(rng.integers(5))]
        q = d[0] + int(rng.integers(1, len(d[1])))
        if not forceable(q, s38[q]):
            continue
        a = add(*d, "bg_noncoding", ga, af=common_af())
        if a is None:
            continue
        b = add(q, s38[q], [x for x in BASES if x != s38[q]][int(rng.integers(3))], "bg_noncoding", gb,
                af=bg_af(), force=True)
        if b is None:
            raise _Retry("SNV inside a deletion record not placeable")
        fuse(a, b, "delsnv")
        done += 1

    if lvl4:                                 # background inside the inverted block
        n_inv, tries, done = int(rng.integers(32, 44)), 0, 0
        while done < n_inv and tries < 6000:
            tries += 1
            p = int(rng.integers(inv["q_a"] + 100, inv["q_b"] - 100))
            if near_exon[p - 40:p + 40].any():
                continue
            af = bg_af()
            if rng.random() < 0.6:
                v = add(p, s38[p], [b for b in BASES if b != s38[p]][int(rng.integers(3))], "bg_noncoding", bg_gt(af), af=af)
            else:
                v = add(*rand_indel(p), "bg_noncoding", bg_gt(af), af=af)
            if v:
                noisy(v)
                done += 1
        if done < 20:
            raise _Retry("too few alleles fit into the inverted block")

    # ---- alleles in hg19-only sequence --------------------------------------------------------
    only19 = []
    gaps = []
    for (t, q, s), (t2, q2, s2) in zip(R.blocks, R.blocks[1:]):
        if t2 > t + s and any(a <= t + s and t2 <= b for a, b in R.windows19):
            gaps.append((t + s, t2))
    big = max(gaps, key=lambda g: g[1] - g[0])
    small = [g for g in gaps if g != big]
    used19 = set()
    n_only = int(rng.integers(26, 34))
    while len(only19) < n_only:
        if len(only19) < n_only - 8:
            p = int(rng.integers(big[0] + 300, big[1] - 300))
        else:
            g = small[int(rng.integers(len(small)))]
            p = int(rng.integers(g[0], g[1]))
        if (R.is_covered19(p) or s19[p] == "N" or any(abs(p - u) < 12 for u in used19)
                or R.is_covered19(p - 1) != R.is_covered19(p + 1) and False):
            continue
        used19.add(p)
        gt = bg_gt(0.2 if rng.random() < 0.7 else None)
        v = dict(pos=None, ref=s19[p], alt=[b for b in BASES if b != s19[p]][int(rng.integers(3))], kind="hg19_only",
                 gene=None, cls=None, af=None, pos19=p, kmax=0, gt=dict(zip("PMF", gt)), rep="norm", group=None)
        v["fmt"] = {r: fmt_good(rng, v["gt"][r]) for r in "PMF"}
        variants.append(v)
        only19.append(v)

    # ---- planted truth on GRCh38: uniqueness, SOP guard on designed classes --------------------
    keys = [(v["pos"], v["ref"], v["alt"]) for v in variants if v["pos"] is not None]
    if len(set(keys)) != len(keys):
        raise _Retry("duplicate planted key")
    panel_set = set(panel)
    by_design = lambda v: [(v["gene"], v["cls"])] if v["gene"] else []      # noqa: E731
    designed = apply_sop(variants, haplotype_cons(variants, models, cds_at, by_design) if lvl4 else by_design, panel_set)
    causal_idx = {i for i, v in enumerate(variants) if v["kind"].startswith("causal")}
    if set(designed) != causal_idx:
        raise _Retry("SOP on designed classes does not return the planted causal set")
    for i, (g, c, m) in designed.items():
        if variants[i].get("joint", c) != c:
            raise RuntimeError("generator bug: joint codon class differs from the designed one")

    # ---- derive the provider VCF (hg19) -------------------------------------------------------
    fam = f"F{int(rng.integers(1000, 9999))}"
    nums = [int(x) for x in rng.permutation(np.arange(11, 99))[:3]]
    sid = {"P": f"{fam}-{nums[0]}", "M": f"{fam}-{nums[1]}", "F": f"{fam}-{nums[2]}"}
    col_roles = [("P", "M", "F")[i] for i in rng.permutation(3)]
    samples = [sid[r] for r in col_roles]
    chrom_in = "22" if nonnorm else "chr22"

    s19x = s19
    if lvl4:
        # the provider's reference: the shipped GRCh37 slice with the planted single-base differences and
        # with the reverse complement of one segment (the inverted block)
        parts, last = [], 0
        for t, b in sorted((R.to19(q)[0], b) for q, b in rd_planted.items()):
            if s19[t] != s38[R.to38(t)[0]] or b == s19[t] or inv["t_a"] <= t < inv["t_b"]:
                raise RuntimeError("generator bug: planted assembly difference at a wrong place")
            parts += [s19[last:t], b]
            last = t + 1
        s19x = "".join(parts + [s19[last:]])
        s19x = s19x[:inv["t_a"]] + revcomp(s19x[inv["t_a"]:inv["t_b"]]) + s19x[inv["t_b"]:]
        if revcomp(s19x[inv["t_a"]:inv["t_b"]]) != s38[inv["q_a"]:inv["q_b"]] or len(s19x) != len(s19):
            raise RuntimeError("generator bug: inverted block")
    recs = []      # (pos19_0based, ref, [alts], [variants in ALT order])
    groups = {}
    for v in variants:
        if v.get("rd") is not None:          # level 4: derived per site below
            continue
        if v.get("inv"):                     # level 4: the allele on the other strand, left-aligned there
            lo = inv["t_a"] + (inv["q_b"] - 1 - (v["pos"] + len(v["ref"]) - 1))
            p, r, a = normalise(s19x, lo, revcomp(v["ref"]), revcomp(v["alt"]))
            v["pos19"], v["id19"], v["rep"] = p, f"{p + 1}_{r}_{a}", "inverted"
            recs.append((p, r, [a], [v]))
            continue
        if v["group"] is not None and nonnorm:
            groups.setdefault(v["group"], []).append(v)
            continue
        p, r, a = v["pos19"], v["ref"], v["alt"]
        if nonnorm and len(r) != len(a) and (v["kind"] in ("causal_ch", "decoy_common_shifted") or rng.random() < 0.65):
            if v["kmax"] >= 1:
                k = int(rng.integers(1, v["kmax"] + 1))
                p, r, a = shifted(s19, p, r, a, k)
                v["rep"] = f"shift{k}"
            else:
                p, r, a = padded(s19, p, r, a, int(rng.integers(1, 4)))
                v["rep"] = "pad"
        elif nonnorm and len(r) == 1 and len(a) == 1 and v["kind"] == "bg_noncoding" and rng.random() < 0.02:
            p, r, a = padded(s19, p, r, a, int(rng.integers(1, 3)))
            v["rep"] = "pad"
        if normalise(s19, p, r, a) != (v["pos19"], v["ref"], v["alt"]):
            raise RuntimeError("generator bug: input representation is not equivalent to the planted allele")
        recs.append((p, r, [a], [v]))
    for (tag, _), (a, b) in groups.items():
        if tag == "multisnv":
            recs.append((a["pos19"], a["ref"], [a["alt"], b["alt"]], [a, b]))
        else:
            off = b["pos"] - a["pos"]
            recs.append((a["pos19"], a["ref"], [a["alt"], a["ref"][:off] + b["alt"] + a["ref"][off + 1:]], [a, b]))
        a["rep"] = b["rep"] = "multiallelic"
    recs.sort(key=lambda x: (x[0], x[1]))

    hdr = ["##fileformat=VCFv4.2", '##FILTER=<ID=PASS,Description="All filters passed">',
           "##source=JointGenotyping", f"##fileDate=2024{int(rng.integers(1, 13)):02d}{int(rng.integers(1, 28)):02d}"]
    if nonnorm:
        hdr.append(f"##contig=<ID=22,length={refdata.HG19_LEN}>")
    else:
        hdr += ["##reference=file:///ref/ucsc.hg19.fasta",
                f"##contig=<ID=chr22,length={refdata.HG19_LEN},assembly=hg19>"]
    hdr += ['##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">',
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">',
            '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples)]
    lines = []
    for p, r, alts, vs in recs:
        cols = []
        for role in col_roles:
            if len(vs) == 1:
                v = vs[0]
                rd, ad, dp, gq = v["fmt"][role]
                cols.append(f"{gt_str(v['gt'][role])}:{rd},{ad}:{dp}:{gq}")
            else:
                a, b = vs
                rd, _, dp, gq = a["fmt"][role]
                if b["gt"][role]:
                    g = "0/2" if b["gt"][role] == 1 else "2/2"
                    cols.append(f"{g}:{rd},0,{b['fmt'][role][1]}:{dp}:{gq}")
                else:
                    cols.append(f"{gt_str(a['gt'][role])}:{rd},{a['fmt'][role][1]},0:{dp}:{gq}")
        qual = int(rng.integers(120, 4000))
        lines.append(f"{chrom_in}\t{p + 1}\t.\t{r}\t{','.join(alts)}\t{qual}\tPASS\t.\tGT:AD:DP:GQ\t" + "\t".join(cols))
    gen = workdir / "gen"
    gen.mkdir(parents=True, exist_ok=True)
    if lvl4:
        # the family's alleles at sites where the assemblies differ, expressed against the provider's reference
        for site in rd_sites:
            if site["depth"] is None:        # everybody homozygous for the provider's base: no record
                continue
            alts = sorted({b for r in "PMF" for b in site["bases"][r]} - {site["b19"]})
            alts = [alts[int(k)] for k in rng.permutation(len(alts))]
            alleles = [site["b19"]] + alts
            cols = []
            for role in col_roles:
                d, dp, gq = site["depth"][role]
                gt = sorted(alleles.index(b) for b in site["bases"][role])
                cols.append(f"{gt[0]}/{gt[1]}:{','.join(str(d.get(b, 0)) for b in alleles)}:{dp}:{gq}")
            qual = int(rng.integers(120, 4000))
            lines.append(f"{chrom_in}\t{site['pos19'] + 1}\t.\t{site['b19']}\t{','.join(alts)}\t{qual}\tPASS\t.\t"
                         "GT:AD:DP:GQ\t" + "\t".join(cols))
        lines.sort(key=lambda ln: (int(ln.split("\t", 2)[1]), ln.split("\t", 4)[3]))
        with open(gen / "GRCh37.chr22.fa", "w") as fh:
            fh.write(">chr22\n")
            for i in range(0, len(s19x), 60):
                fh.write(s19x[i:i + 60] + "\n")
        (gen / "hg19ToHg38.chr22.over.chain").write_text(inv["fwd"])
    (gen / "trio.vcf").write_text("\n".join(hdr + lines) + "\n")
    mapped = sorted((v for v in variants if v["pos"] is not None), key=lambda v: (v["pos"], v["ref"], v["alt"]))
    (gen / "truth.sites.vcf").write_text(
        "##fileformat=VCFv4.2\n" + f"##contig=<ID=chr22,length={refdata.HG38_LEN}>\n"
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + "".join(f"chr22\t{v['pos'] + 1}\t.\t{v['ref']}\t{v['alt']}\t.\t.\t.\n" for v in mapped))
    (gen / "check.sh").write_text((ASSETS / ("gen_check_l4.sh" if lvl4 else "gen_check.sh")).read_text())

    # ---- guards inside the tool image: csq on the truth, two independent lift-overs of the input -
    tag = dockerlib.tool_image(conda=refdata.TOOLS)
    res = refdata.DERIVED / "resources"
    r = dockerlib.run_in_image(tag, f"RENAME={'1' if nonnorm else '0'} bash /work/check.sh",
                               mounts={gen: "/work", res: "/res"}, copy_out={"/out": gen / "out"}, cpus=2)
    out = gen / "out"
    if lvl4:
        (gen / "GRCh37.chr22.fa").unlink()      # 51 MB of plain text; the compressed copy is in out/
    if r.returncode != 0 or not (out / "done").exists():
        raise RuntimeError("generation check failed in the tool image:\n" + (r.stdout + r.stderr)[-3000:])

    def tsv(name):
        return [ln.split("\t") for ln in (out / name).read_text().splitlines() if ln]

    tn = [(int(a) - 1, b, c) for a, b, c in tsv("truth_norm.tsv")]
    if tn != [(v["pos"], v["ref"], v["alt"]) for v in mapped]:
        raise RuntimeError("generator bug: planted GRCh38 alleles are not in normal form according to bcftools norm")
    csq = {}
    for a, b, c, d in tsv("truth_csq.tsv"):
        cons = []
        for e in d.split(","):
            f = e.split("|")
            if len(f) >= 2 and f[0] not in (".", ""):
                cons.append((f[1], f[0].lstrip("*")))
        csq[(int(a) - 1, b, c)] = cons
    for v in variants:
        if v["pos"] is None:
            continue
        got = csq.get((v["pos"], v["ref"], v["alt"]), [])
        v["csq"] = got
        k = v["kind"]
        if k.startswith("causal") or k.startswith("decoy_same") or k in ("decoy_parent_quality", "decoy_allele_balance",
                                                                        "decoy_common_shifted", "decoy_cis_rescue",
                                                                        "decoy_hidden_rescue", "rd_hidden_coding"):
            if got != [(v["gene"], v["cls"])]:
                raise _Retry(f"csq disagrees with the designed consequence of {k}: {got} vs {v['gene']}/{v['cls']}")
        if k in ("decoy_syn_dn", "decoy_wrong_build") and any(c.split("&")[0] in QUALIFYING or
                                                              any(x in QUALIFYING for x in c.split("&")) for _, c in got):
            raise _Retry(f"csq gives a qualifying consequence for {k}")
    of_csq = lambda v: [(g, x) for g, c in v.get("csq", []) for x in c.split("&")]      # noqa: E731
    by_csq = apply_sop(variants, haplotype_cons(variants, models, cds_at, of_csq) if lvl4 else of_csq, panel_set)
    if set(by_csq) != causal_idx:
        raise _Retry("SOP on csq consequences does not return the planted causal set")

    ids = {ln[0] for ln in tsv("split_ids.txt")}
    lifted = {}
    for name in ("picard.tsv", "plugin.tsv"):
        d = {}
        for f in tsv(name):
            d.setdefault(f[0], []).append((int(f[1]) - 1, f[2], f[3], tuple(f[4:7])))
        lifted[name] = d
    out_samples = [ln[0] for ln in tsv("samples.txt")]
    gray_keys, n_gray = set(), 0
    mm_all = [m[0] for m in real_mismatches(R)] if lvl4 else []
    gray_rd = {q for q in mm_all if any(0 < abs(q - o) <= 2 for o in mm_all)}   # neighbouring differences: one
    for v in variants:                                                           # record or two is a convention
        if v.get("rd") is not None:          # exact by construction; never looked up in the tool outputs
            v["scored"] = v["pos"] not in gray_rd
            continue
        vid = v.get("id19") or f"{v['pos19'] + 1}_{v['ref']}_{v['alt']}"
        if vid not in ids:
            raise RuntimeError(f"generator bug: {vid} not found after split/left-align on hg19")
        hits = [lifted[n].get(vid, []) for n in ("picard.tsv", "plugin.tsv")]
        if v["pos"] is None:
            v["scored"] = not hits[0] and not hits[1]
        else:
            want = (v["pos"], v["ref"], v["alt"])
            exp_gt = tuple(gt_str(v["gt"][r]) for r in col_roles)
            v["scored"] = all(len(h) == 1 and h[0][:3] == want and
                              tuple(g.replace("|", "/") for g in h[0][3]) == exp_gt for h in hits)
        if not v["scored"]:
            n_gray += 1
            if v["kind"].startswith(("causal", "decoy")):
                raise _Retry(f"{v['kind']} is not reproduced by both lift-over tools")
            for h in hits:
                gray_keys.update((x[0] + 1, x[1], x[2]) for x in h)
            if v["pos"] is not None:
                gray_keys.add((v["pos"] + 1, v["ref"], v["alt"]))
    if out_samples != samples:
        raise RuntimeError("sample order changed in the tool image")
    if n_gray > 0.02 * len(variants):
        raise _Retry(f"gray zone too large: {n_gray}")

    # ---- shipped inputs -----------------------------------------------------------------------
    af_rows = [(v["pos"] + 1, v["ref"], v["alt"], v["af"]) for v in mapped if v["af"] is not None]
    n_extra = int(rng.integers(2500, 3200))
    while n_extra:
        p = rand_noncoding() if rng.random() < 0.8 else int(rng.integers(*R.windows38[int(rng.integers(5))]))
        if occ[p - 3:p + 3].any() or s38[p] == "N":
            continue
        occ[p] = True
        af_rows.append((p + 1, s38[p], [b for b in BASES if b != s38[p]][int(rng.integers(3))],
                        float(f"{10 ** rng.uniform(-5.2, -0.3):.3g}")))
        n_extra -= 1
    af_rows.sort()
    popfreq = "chrom\tpos\tref\talt\taf\n" + "".join(f"chr22\t{p}\t{r}\t{a}\t{af:.6g}\n" for p, r, a, af in af_rows)
    ped = ("#family\tindividual\tfather\tmother\tsex\tphenotype\n"
           + "".join(sorted([f"{fam}\t{sid['P']}\t{sid['F']}\t{sid['M']}\t{int(rng.integers(1, 3))}\t2\n",
                             f"{fam}\t{sid['F']}\t0\t0\t1\t1\n", f"{fam}\t{sid['M']}\t0\t0\t2\t1\n"])))
    input_files = {
        "trio.vcf.gz": out / "trio.vcf.gz", "trio.vcf.gz.tbi": out / "trio.vcf.gz.tbi",
        "pedigree.ped": ped, "panel_genes.txt": "\n".join(panel) + "\n",
        "lab_sop.md": SOP.format(dp=MIN_DP, gq=MIN_GQ, ablo=AB_LO, abhi=AB_HI, af=AF_CUT),
        "resources/popfreq.GRCh38.tsv.gz": gzip.compress(popfreq.encode(), mtime=0),
    }
    for f in sorted(res.iterdir()):
        input_files[f"resources/{f.name}"] = f
    if lvl4:
        for name in ("GRCh37.chr22.fa.gz", "GRCh37.chr22.fa.gz.fai", "GRCh37.chr22.fa.gz.gzi"):
            input_files[f"resources/{name}"] = out / name
        input_files["lab_sop.md"] = SOP_L4.format(dp=MIN_DP, gq=MIN_GQ, ablo=AB_LO, abhi=AB_HI, af=AF_CUT)
        input_files["resources/hg19ToHg38.chr22.over.chain"] = inv["fwd"]
        input_files["resources/hg38ToHg19.chr22.over.chain"] = inv["rev"]

    # ---- truth --------------------------------------------------------------------------------
    sites = [dict(pos=v["pos"] + 1, ref=v["ref"], alt=v["alt"], scored=bool(v["scored"]), kind=v["kind"],
                  gt={sid[r]: gt_str(v["gt"][r]) for r in "PMF"}) for v in mapped]
    cands = [dict(pos=variants[i]["pos"] + 1, ref=variants[i]["ref"], alt=variants[i]["alt"], gene=g, consequence=c,
                  model=m) for i, (g, c, m) in sorted(designed.items())]
    truth = dict(samples=samples, proband=sid["P"], sites=sites, gray_keys=sorted(gray_keys),
                 unmapped=[dict(pos=v["pos19"] + 1, ref=v["ref"], alt=v["alt"], scored=bool(v["scored"])) for v in only19],
                 candidates=cands, level=level)
    if lvl4:
        for st, v in zip(sites, mapped):
            st["cat"] = "inverted" if v.get("inv") else ("refdiff" if v.get("rd") is not None else "standard")
        truth["inverted_interval"] = [inv["q_a"] + 1, inv["q_b"]]
        truth["refdiff_positions"] = sorted({q + 1 for q in mm_all} | {q + 1 for q in rd_planted})
        truth["gray_positions"] = sorted(q + 1 for q in gray_rd)
    stats = dict(attempt=attempt, n_input_records=len(lines), n_alleles=len(variants), n_mapped=len(mapped),
                 n_hg19_only=len(only19), n_gray=n_gray, n_panel=len(panel), with_hom_rec=with_hr,
                 multiallelic_causal=multi_causal,
                 n_noncanonical=sum(v["rep"] != "norm" for v in variants),
                 kinds={k: sum(v["kind"] == k for v in variants) for k in sorted({v["kind"] for v in variants})})
    if lvl4:
        stats["level4"] = dict(
            planted_reference_bases=len(rd_planted), real_assembly_differences=len(mm_all),
            sites_without_record=sum(st["depth"] is None for st in rd_sites),
            sites_with_record=sum(st["depth"] is not None for st in rd_sites),
            refdiff_alleles_scored=sum(1 for v in mapped if v.get("rd") is not None and v["scored"]),
            inverted_block_gene=inv["gene"], inverted_block_grch38=[inv["q_a"] + 1, inv["q_b"]],
            inverted_block_provider=[inv["t_a"] + 1, inv["t_b"]], inverted_block_length=inv["q_b"] - inv["q_a"],
            inverted_alleles_scored=sum(1 for v in mapped if v.get("inv") and v["scored"]),
            inverted_indels_scored=sum(1 for v in mapped if v.get("inv") and v["scored"] and len(v["ref"]) != len(v["alt"])),
            n_candidates=len(cands))
    (gen / "stats.json").write_text(json.dumps(stats, indent=1))
    (gen / "truth.json").write_text(json.dumps(truth))

    # ---- probes / attacks ---------------------------------------------------------------------
    def vcf_text(site_list, drop=(), phased=False, extra=()):
        h = ["##fileformat=VCFv4.2", f"##contig=<ID=chr22,length={refdata.HG38_LEN},assembly=GRCh38>",
             '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples)]
        body = []
        for j, s in enumerate(site_list):
            if j in drop:
                continue
            gts = [s["gt"][x].replace("/", "|") if phased else s["gt"][x] for x in samples]
            body.append((s["pos"], f"chr22\t{s['pos']}\t.\t{s['ref']}\t{s['alt']}\t.\tPASS\t.\tGT\t" + "\t".join(gts)))
        body += list(extra)
        body.sort(key=lambda x: x[0])
        return gzip.compress(("\n".join(h + [b for _, b in body]) + "\n").encode(), mtime=0)

    def unm_text(rows, chrom=chrom_in):
        return "chrom_in\tpos_in\tref\talt\n" + "".join(f"{chrom}\t{u['pos']}\t{u['ref']}\t{u['alt']}\n" for u in rows)

    def cand_text(rows, suffix="", extra_cols=False):
        h = "chrom\tpos\tref\talt\tgene\tconsequence\tmodel" + ("\tnote" if extra_cols else "") + "\n"
        return h + "".join(f"chr22\t{c['pos']}\t{c['ref']}\t{c['alt']}\t{c['gene']}\t{c['consequence']}{suffix}\t"
                           f"{c['model']}" + ("\treviewed" if extra_cols else "") + "\n" for c in rows)

    fake_idx = gzip.compress(b"TBI\x01" + b"\x00" * 32, mtime=0)
    prng = np.random.default_rng(seed + 77)
    scored_idx = [j for j, s in enumerate(sites) if s["scored"]]
    drop = set(int(x) for x in prng.choice(scored_idx, size=max(1, len(scored_idx) // 150), replace=False))
    perfect = {"trio.GRCh38.vcf.gz": vcf_text(sites), "trio.GRCh38.vcf.gz.tbi": fake_idx,
               "unmapped.tsv": unm_text(truth["unmapped"]), "candidates.tsv": cand_text(cands)}
    decoy_rows = [dict(pos=v["pos"] + 1, ref=v["ref"], alt=v["alt"], gene=v["gene"], consequence=v["cls"], model="de_novo")
                  for v in variants if v["kind"] in ("decoy_parent_quality", "decoy_allele_balance", "decoy_nonpanel_dn")]
    probes = [
        Probe("planted_truth", "submission written directly from the planted truth", 0.999, files=perfect),
        Probe("within_tolerance", "0.7% of the scored records missing, one hg19-only allele missing, phased GT, "
              "UCSC contig name in unmapped.tsv, *_variant consequence names, shuffled rows, extra column", 0.99,
              files={"trio.GRCh38.vcf.gz": vcf_text(sites, drop=drop, phased=True),
                     "trio.GRCh38.vcf.gz.csi": fake_idx,
                     "unmapped.tsv": unm_text(truth["unmapped"][1:], chrom="chr22"),
                     "candidates.tsv": cand_text(cands[::-1], suffix="_variant", extra_cols=True)}),
        Probe("independent_solution", "second implementation: own chain lift-over, own left-alignment and own "
              "consequence prediction in Python (no Picard, no bcftools csq)", 0.99,
              script="SOLUTION=2 bash /solution/solve.sh"),
    ]
    pipe = {"trio_pipeline.py": ASSETS / "trio_pipeline.py"}
    attacks = [
        Attack("constant_answer", "valid empty call set, no unmapped alleles, no candidates", 0.05,
               files={"trio.GRCh38.vcf.gz": vcf_text([]), "trio.GRCh38.vcf.gz.tbi": fake_idx,
                      "unmapped.tsv": unm_text([]), "candidates.tsv": cand_text([])}),
        Attack("every_panel_variant", "perfect call set, but every allele seen in the proband is listed as de novo", 0.40,
               files={**perfect, "candidates.tsv": cand_text(
                   [dict(pos=s["pos"], ref=s["ref"], alt=s["alt"], gene=panel[0], consequence="missense", model="de_novo")
                    for s in sites if s["gt"][sid["P"]] != "0/0"][:190])}),
        Attack("one_false_candidate", "planted truth plus one extra candidate row (just beyond tolerance)",
               0.89 if lvl4 else 0.88, files={**perfect, "candidates.tsv": cand_text(cands + decoy_rows[:1])}),
        Attack("one_missing_candidate", "planted truth with one causal variant missing", 0.89 if lvl4 else 0.80,
               files={**perfect, "candidates.tsv": cand_text(cands[1:])}),
    ]
    if nonnorm:
        attacks += [
            Attack("naive_no_liftover", "contig renamed, coordinates left on the provider's assembly and annotated as "
                   "if they were GRCh38", 0.30, script="python3 /alegen_script/trio_pipeline.py --no-lift", script_files=pipe),
            Attack("naive_no_rename", "lift-over run on the provider's contig name as delivered: every record is "
                   "rejected and reported as unmapped", 0.30,
                   script="python3 /alegen_script/trio_pipeline.py --no-rename", script_files=pipe),
            Attack("naive_no_normalise", "lifted correctly, but multi-allelic records are not split and alleles are not "
                   "left-aligned before frequency lookup and reporting", 0.75,
                   script="python3 /alegen_script/trio_pipeline.py --no-norm", script_files=pipe),
        ]
    if decoys3:
        attacks += [
            Attack("naive_ignore_quality", "correct lift-over and representation, SOP genotype-quality and allele-balance "
                   "rules not applied", 0.75, script="python3 /alegen_script/trio_pipeline.py --ignore-quality",
                   script_files=pipe),
            Attack("naive_phase_blind", "correct lift-over and representation, any two rare damaging heterozygous "
                   "variants in a gene called compound_het", 0.75,
                   script="python3 /alegen_script/trio_pipeline.py --phase-blind", script_files=pipe),
            Attack("naive_default_pipeline", "no lift-over, no normalisation, no quality rules, phase-blind", 0.30,
                   script="python3 /alegen_script/trio_pipeline.py --no-lift --no-norm --ignore-quality --phase-blind",
                   script_files=pipe),
        ]
    if lvl4:
        l4 = {"trio_pipeline_l4.py": ASSETS / "trio_pipeline_l4.py"}
        run4 = "python3 /alegen_script/trio_pipeline_l4.py "
        attacks = [a for a in attacks if a.name not in ("naive_no_normalise", "naive_ignore_quality", "naive_phase_blind")]
        hidden_keys = {(v["pos"] + 1, v["ref"], v["alt"]) for v in variants if v.get("rd") is not None
                       and v["rd"]["depth"] is None}

        def provider_row(i, c):
            """The causal variant as the provider's file shows it (position and alleles on the provider's reference)."""
            v = variants[i]
            if v.get("rd") is not None:
                st = v["rd"]
                key = (st["pos19"] + 1, st["b19"], st["b38"] if v["alt"] == st["b19"] else v["alt"])
            elif v.get("inv"):
                q, r_, a_ = v["id19"].split("_")
                key = (int(q), r_, a_)
            else:
                key = (v["pos19"] + 1, v["ref"], v["alt"])
            return dict(c, pos=key[0], ref=key[1], alt=key[2])
        attacks += [
            Attack("l3_reference", "the level-3 reference method on the level-4 inputs: Picard LiftoverVcf (records "
                   "whose REF is not the GRCh38 base are rejected and reported as unmapped), no sites added, one "
                   "consequence per variant", 0.75, script="python3 /alegen_script/trio_pipeline.py", script_files=pipe),
            Attack("ignore_reference_differences", "P_reference_differences ignored, everything else as in the "
                   "reference solution: records at positions where the assemblies differ are rejected, no sites added",
                   0.75, script=run4 + "--liftover picard --no-hidden", script_files=l4),
            Attack("records_only_no_added_sites", "the lift-over tool re-expresses the provider's records at positions "
                   "where the assemblies differ, but positions without a record are not considered", 0.75,
                   script=run4 + "--no-hidden", script_files=l4),
            Attack("added_sites_missing_from_truth", "planted truth without the alleles at positions that have no "
                   "record in the provider's file (candidates as planted)", 0.89,
                   files={**perfect, "trio.GRCh38.vcf.gz": vcf_text(
                       [st for st in sites if (st["pos"], st["ref"], st["alt"]) not in hidden_keys])}),
            Attack("split_before_liftover", "multi-allelic records are split against the provider's reference before "
                   "the lift-over (the other ALT then counts as the provider's reference allele, which turns into a "
                   "wrong ALT copy where the reference base changes); everything else as in the reference solution",
                   0.89, script=run4 + "--split-first", script_files=l4),
            Attack("inverted_block_unmapped", "P_minus_strand_block ignored: only plus-strand chains are used, the "
                   "records of the inverted block are reported as unmapped; everything else as in the reference "
                   "solution", 0.75, script=run4 + "--drop-inverted unmapped", script_files=l4),
            Attack("inverted_block_lost", "positions of the inverted block are lifted without reverse complementing "
                   "the alleles; the records then contradict the GRCh38 reference and are dropped silently", 0.75,
                   script=run4 + "--drop-inverted silent", script_files=l4),
            Attack("per_variant_consequence", "P_haplotype_consequence ignored, everything else as in the reference "
                   "solution: every variant gets the consequence it has on its own", 0.75,
                   script=run4 + "--per-variant", script_files=l4),
            Attack("codon_merge_ignoring_phase", "variants of the proband in one codon are always assessed together, "
                   "whichever chromosome copy carries them", 0.75, script=run4 + "--merge-copies", script_files=l4),
            Attack("provider_genotypes", "candidates taken from the provider's coordinates: inheritance patterns are "
                   "read from the genotypes as written against the provider's reference", 0.75,
                   script=run4 + "--provider-genotypes", script_files=l4),
            Attack("candidates_at_provider_coordinates", "planted call set and planted candidate list, but the candidate "
                   "rows carry position and alleles as in the provider's file", 0.40,
                   files={**perfect, "candidates.tsv": cand_text(
                       [provider_row(i, c) for (i, _), c in zip(sorted(designed.items()), cands)])}),
            Attack("naive_no_normalise", "level-4 reference solution without splitting / left-alignment", 0.75,
                   script=run4 + "--no-norm", script_files=l4),
            Attack("naive_ignore_quality", "level-4 reference solution without the SOP genotype-quality and "
                   "allele-balance rules", 0.75, script=run4 + "--ignore-quality", script_files=l4),
            Attack("naive_phase_blind", "level-4 reference solution, any two rare damaging heterozygous variants in "
                   "a gene called compound_het", 0.75, script=run4 + "--phase-blind", script_files=l4),
        ]

    traps = []
    if nonnorm:
        traps += [
            Trap("build_mismatch", "The provider call set is on hg19/GRCh37 but nothing says so; annotating it against "
                 "the GRCh38 resources puts every variant in the wrong place (and one record then looks like a de novo "
                 "missense variant of a panel gene).",
                 "the ##contig length 51304566 is the hg19 length of chr22 (GRCh38: 50818468); REF alleles match "
                 "GRCh37.chr22.fa.gz at 100% and GRCh38.chr22.fa.gz at ~10-12%; bcftools csq/norm stop with 'the fasta "
                 "reference does not match the VCF REF allele'; lab_sop.md says all lab resources are on GRCh38",
                 "Understanding / unverified assumption about inputs"),
            Trap("contig_naming", "The VCF names the contig '22' while FASTA, chain, GFF3 and frequency table use 'chr22'; "
                 "Picard rejects every record (NoTarget), bcftools norm stops with 'sequence 22 not found'.",
                 "tool error/100% reject rate on the first run; bcftools view -h | grep contig",
                 "Execution / tool-interface mismatch"),
            Trap("representation", "A quarter of the alleles are right-shifted, padded or share a multi-allelic record "
                 "(a causal SNV is the 2nd ALT of a record whose 1st ALT is a common deletion; a common frameshift is "
                 "delivered right-shifted so a literal lookup in the frequency table misses it).",
                 "bcftools norm --check-ref/-f reports realigned and split records; multi-allelic ALT columns are "
                 "visible in the VCF; the frequency table and the SOP state the canonical representation",
                 "Approach / missing normalisation step"),
        ]
    if decoys3:
        traps += [
            Trap("genotype_quality_decoys", "Four apparent de novo variants in panel genes rest on an unusable parental "
                 "genotype (low DP/GQ) or on a proband allele balance of 9-19%.",
                 "rules 2 and 3 of lab_sop.md; the FORMAT fields DP, GQ and AD of the three samples",
                 "Verification / skipping stated acceptance criteria"),
            Trap("same_parent_pairs", "Two panel genes carry two rare damaging heterozygous variants that both come from "
                 "the same parent (in cis): not a compound heterozygote.",
                 "rule 7 of lab_sop.md; parental genotypes in the VCF and the pedigree",
                 "Understanding / domain knowledge gap"),
        ]
    if lvl4:
        n_hidden = stats["level4"]["sites_without_record"]
        traps += [
            Trap("P_reference_differences", f"Inside aligned chain blocks the provider's reference and GRCh38 differ at "
                 f"{len(truth['refdiff_positions'])} single bases. A provider record whose ALT is the GRCh38 base is the "
                 f"opposite allele with complementary genotypes on GRCh38 (a default Picard run rejects it); {n_hidden} "
                 "of these positions have no provider record at all because the family is homozygous for the provider's "
                 "base, which is a homozygous ALT allele on GRCh38 that no lift-over tool can produce from the input; "
                 "some sites carry a third allele. A homozygous-recessive variant, a compound-het partner and a de novo "
                 "variant (third allele) sit on such sites.",
                 "the instruction defines the call set by meaning (how the family differs from GRCh38, given that the "
                 "provider's file is complete against the provider's reference); Picard reports MismatchedRefAllele "
                 "rejects although the positions lift; bcftools +liftover reports swapped/added reference alleles; "
                 "comparing the two FASTA files through the chain shows the differing bases",
                 "Approach / prior decoupled from truth (tool output taken as the deliverable)"),
            Trap("P_minus_strand_block", f"A {inv['q_b'] - inv['q_a']}-bp segment around an exon of a panel gene is "
                 "inverted between the provider's reference and GRCh38: the chain files carry it as a separate "
                 f"minus-strand chain. It holds {stats['level4']['inverted_alleles_scored']} scored alleles "
                 f"({stats['level4']['inverted_indels_scored']} indels, which have to be reverse complemented, "
                 "re-anchored and left-aligned again) including a compound-het pair (frameshift + missense).",
                 "the header line of the extra chain has '-' in the query strand column; Picard LiftoverVcf and "
                 "bcftools +liftover handle it (INFO/FLIP); a position-only lift-over yields REF alleles that "
                 "contradict GRCh38 for every record of the block",
                 "Execution / coordinate conversion without strand"),
            Trap("P_haplotype_consequence", "SOP rule 5 assesses single-nucleotide variants that share a codon on one "
                 "chromosome copy of the proband together. Planted: two missense variants in one codon on different "
                 "copies (a true compound heterozygote; merged they would be silent), two qualifying variants in one "
                 "codon on the same copy that are silent together (their gene must not be reported), and de novo "
                 "variants next to a base where the provider's reference differs from GRCh38 (silent alone and "
                 "damaging on the real haplotype, or the reverse).",
                 "rule 5 of lab_sop.md states the definition; neighbouring records in the same codon are visible in "
                 "any annotation table; bcftools csq without -l reports compound consequences for phased genotypes",
                 "Approach / per-variant logic where the stated rule is per haplotype"),
            Trap("C_dependent_stages", "Lift-over -> call set against GRCh38 -> phase assignment -> consequence -> "
                 "inheritance model: a position without a provider record that is missed in the first stage changes "
                 "the consequence of a de novo variant in the same codon and thereby the candidate list (one causal "
                 "variant lost, one false candidate gained).",
                 "the candidate list is derived from the delivered call set; rule 5 says that every variant the "
                 "proband carries counts", "Planning / early error surfaces only in the final deliverable"),
            Trap("Q_gradient", f"{len(cands)} scored candidates in {len({c['gene'] for c in cands})} genes instead "
                 "of 3-4 in 2-3 genes, so that the candidate "
                 "components give partial credit per element. This is a scoring change, not an in-data control in "
                 "the sense of AUTHORING.md's Q element.", "n/a (no decision point)", "n/a"),
        ]
        return TaskBuild(
            task_id=f"trio-liftover-s{seed:04d}-l{level}", family="trio_liftover", seed=seed, level=level,
            title="Rare-disease trio: family call set against GRCh38 and haplotype-aware SOP candidate review",
            summary="Express a provider's joint trio call set against GRCh38 (including the positions where the two "
                    "reference assemblies differ), report the alleles in unaligned sequence, and return the candidate "
                    "causal variants of the proband under a lab SOP with per-copy codon consequences.",
            subdomain="3.2 Genomics & Sequence Analysis",
            software=["bcftools", "bcftools +liftover", "bcftools csq", "samtools/htslib", "Picard", "Python/pysam"],
            instruction_md=INSTRUCTION_L4,
            agent_dockerfile=dockerlib.agent_dockerfile(conda=refdata.TOOLS),
            verifier_dockerfile=dockerlib.verifier_dockerfile(pip=[]),
            input_files=input_files,
            solution_files={"solve.sh": ASSETS / "solve_l4.sh", "trio_pipeline_l4.py": ASSETS / "trio_pipeline_l4.py",
                            "solve_alt_l4.py": ASSETS / "solve_alt_l4.py"},
            tests_files={"score.py": ASSETS / "score_l4.py", "reference/truth.json": json.dumps(truth)},
            traps=traps, attacks=attacks, probes=probes,
            capabilities=["verifying the genome build of delivered data", "cross-build lift-over with reject handling",
                          "re-expressing genotypes where two reference assemblies differ",
                          "variant normalisation before database lookup", "trio phasing by transmission",
                          "haplotype-aware consequence prediction", "trio inheritance reasoning under a written SOP",
                          "exact output contracts"],
            judge_explanation_md=JUDGE_L4, pass_threshold=0.9,
            oracle_steps=L4_ORACLE_STEPS, oracle_runtime_min=0.3, cpus=4, memory_mb=8000, storage_mb=10000,
            extra_metadata={"n_input_records": stats["n_input_records"], "n_alleles": stats["n_alleles"],
                            "n_gray_zone": stats["n_gray"], "generation_attempt": attempt, **stats["level4"]},
        )

    return TaskBuild(
        task_id=f"trio-liftover-s{seed:04d}-l{level}", family="trio_liftover", seed=seed, level=level,
        title="Rare-disease trio: provider call set to GRCh38 and SOP candidate review",
        summary="Bring a provider's joint trio call set onto GRCh38 in canonical representation, report the alleles "
                "that cannot be placed, and return the candidate causal variants of the proband under the lab SOP.",
        subdomain="3.2 Genomics & Sequence Analysis",
        software=["bcftools", "bcftools +liftover", "samtools/htslib", "Picard", "Python/pysam"],
        instruction_md=INSTRUCTION,
        agent_dockerfile=dockerlib.agent_dockerfile(conda=refdata.TOOLS),
        verifier_dockerfile=dockerlib.verifier_dockerfile(pip=[]),
        input_files=input_files,
        solution_files={"solve.sh": ASSETS / "solve.sh", "trio_pipeline.py": ASSETS / "trio_pipeline.py",
                        "solve_alt.py": ASSETS / "solve_alt.py"},
        tests_files={"score.py": ASSETS / "score.py", "reference/truth.json": json.dumps(truth)},
        traps=traps, attacks=attacks, probes=probes,
        capabilities=["verifying the genome build of delivered data", "cross-build lift-over with reject handling",
                      "variant normalisation before database lookup", "trio inheritance reasoning under a written SOP",
                      "exact output contracts"],
        judge_explanation_md=JUDGE, pass_threshold=0.9,
        oracle_steps=14, oracle_runtime_min=0.2, cpus=4, memory_mb=8000, storage_mb=10000,
        extra_metadata={"n_input_records": stats["n_input_records"], "n_alleles": stats["n_alleles"],
                        "n_gray_zone": stats["n_gray"], "generation_attempt": attempt},
    )
