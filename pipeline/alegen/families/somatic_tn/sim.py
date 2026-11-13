"""Forward simulator for the somatic_tn family (host side, numpy only).

Everything that is scored is decided here before a read exists:
  * the panel reference (a slice of hg38 chr22, optionally with one exact intra-slice duplication),
  * the germline genotype of every patient at a set of population SNP sites,
  * the somatic SNVs / indels of every tumour and their allele fractions,
  * library artefacts that exist in one read-pair orientation only,
  * which fragment carries which allele.
Reads are then written as paired FASTQ. No aligner or caller is run to obtain the truth.
"""
from __future__ import annotations

import bisect
import gzip
import hashlib
import shutil
import subprocess
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[4]
CACHE = REPO / ".cache" / "somatic_tn"
SRC_CACHE = REPO / ".cache" / "variant_review" / "chr22.hg38.fa.gz"
REF_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz"
REF_SHA256 = "05f9d97d6fbfd08a44ca45b50837ca2ae9c471f35ba79dffec04d2cb5eaaf695"

CONTIG = "chr22_panel"
RL = 100                      # read length (2 x 100)
L0 = 120_000                  # length of the real slice
N_TARGETS = 12
EDGE = 200                    # planted sites keep this distance from target edges (full depth there)
K = 35                        # k-mer size of the uniqueness mask
MIN_CANDIDATES = 15000        # eligible positions a reference window must offer
DUP_MARGIN = 470              # dup-interior variants keep this distance from the copy edges (> max insert)
INS_MEAN, INS_SD, INS_MIN, INS_MAX = 270, 45, 150, 450
QUAL_P = [0.955, 0.035, 0.010]          # fraction of bases with quality F (37), : (25), , (11)
QUAL_ERR = [0.0004, 0.004, 0.04]        # substitution error rate per quality class
QUAL_CH = np.frombuffer(b"F:,", dtype=np.uint8)
DUP_P = [0.92, 0.065, 0.015]            # a fragment is sequenced 1 / 2 / 3 times (PCR duplicates)

# scoring floor (see docs/families/somatic_tn.md): a planted somatic variant is scored only if the
# simulated tumour library holds enough evidence for it
MIN_ALT_FRAG = 8
MIN_AF = 0.07
MIN_MINOR_FRAG = 3
MIN_MINOR_FRAC = 0.25

COMP = str.maketrans("ACGT", "TGCA")
BASES = "ACGT"


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


# --------------------------------------------------------------------------- reference

def load_chr22() -> str:
    """chr22 of hg38 with UCSC soft-masking preserved (lower case = repeat)."""
    fa = CACHE / "chr22.hg38.fa.gz"
    if not fa.exists():
        CACHE.mkdir(parents=True, exist_ok=True)
        if SRC_CACHE.exists():
            shutil.copyfile(SRC_CACHE, fa)
        else:
            subprocess.run(["curl", "-fsSL", "-o", str(fa), REF_URL], check=True)
    if hashlib.sha256(fa.read_bytes()).hexdigest() != REF_SHA256:
        raise RuntimeError("reference checksum mismatch")
    with gzip.open(fa, "rt") as fh:
        next(fh)
        return "".join(line.strip() for line in fh)


def _pick_targets(mask: np.ndarray, rng, long_len: int | None):
    """12 non-overlapping capture targets in low-repeat sequence. Returns sorted [(s, e)] or None."""
    L = len(mask)
    csum = np.concatenate([[0], np.cumsum(mask)])
    lens = [int(rng.integers(2600, 4400)) for _ in range(N_TARGETS)]
    if long_len:
        lens[0] = long_len
    out: list[tuple[int, int]] = []
    for ln in lens:
        ok = False
        for _ in range(400):
            s = int(rng.integers(600, L - ln - 600))
            e = s + ln
            if (csum[e] - csum[s]) / ln > 0.35:
                continue
            if any(s < b + 1500 and a < e + 1500 for a, b in out):
                continue
            out.append((s, e))
            ok = True
            break
        if not ok:
            return None
    return out


def build_reference(seed: int, with_dup: bool) -> dict:
    """Pick the slice, the targets and (optionally) plant one exact duplication."""
    chrom = load_chr22()
    rng = np.random.default_rng(70_000 + seed)
    dup_len = int(rng.integers(3400, 4000))
    long_len = dup_len + int(rng.integers(700, 1000))
    for _ in range(400):
        s0 = int(rng.integers(18_000_000, 48_000_000))
        w = chrom[s0:s0 + L0]
        up = w.upper()
        if "N" in up:
            continue
        gc = (up.count("G") + up.count("C")) / L0
        if not 0.36 <= gc <= 0.60:
            continue
        mask = np.frombuffer(w.encode(), dtype=np.uint8) >= 97      # lower case
        targets = _pick_targets(mask, rng, long_len)
        if targets is None:
            continue
        src = targets[0]                      # the long target holds the duplicated segment
        dup = None
        if with_dup:
            ds = src[0] + int(rng.integers(300, src[1] - src[0] - dup_len - 300))
            st = sorted(targets)
            gaps = [(st[i][1], st[i + 1][0]) for i in range(len(st) - 1)]
            gaps = [g for g in gaps if g[1] - g[0] >= 2400]
            if not gaps:
                continue
            g = gaps[int(rng.integers(0, len(gaps)))]
            y = (g[0] + g[1]) // 2
            seg = up[ds:ds + dup_len]
            up = up[:y] + seg + up[y:]
            mask = np.concatenate([mask[:y], mask[ds:ds + dup_len], mask[y:]])
            targets = [((a + dup_len, b + dup_len) if a >= y else (a, b)) for a, b in targets]
            ds2 = ds + dup_len if ds >= y else ds
            targets.append((y - 350, y + dup_len + 350))
            dup = {"a": [ds2, ds2 + dup_len], "b": [y, y + dup_len], "len": dup_len}
        targets = sorted(targets)
        ref = {"seq": up, "mask": mask, "targets": targets, "dup": dup,
               "source": {"chrom": "chr22", "start0": s0, "length": L0, "assembly": "hg38"}}
        ref["cand"] = candidate_sites(ref)
        if len(ref["cand"]) < MIN_CANDIDATES:
            continue                           # too repeat-rich to place every planted site with spacing
        return ref
    raise RuntimeError("no suitable reference window")


def unique_mask(seq: str) -> np.ndarray:
    """True where no K-mer overlapping the position occurs twice in the reference (either strand)."""
    n = len(seq) - K + 1
    cnt: dict[str, int] = {}
    keys = []
    for i in range(n):
        k = seq[i:i + K]
        r = rc(k)
        if r < k:
            k = r
        keys.append(k)
        cnt[k] = cnt.get(k, 0) + 1
    rep = np.zeros(len(seq) + 1, dtype=np.int32)
    for i, k in enumerate(keys):
        if cnt[k] > 1:
            rep[i] += 1
            rep[i + K] -= 1
    return np.cumsum(rep[:-1]) == 0


def left_normalise(seq: str, pos0: int, ref: str, alt: str) -> tuple[int, str, str]:
    """Parsimonious, left-aligned representation (same convention as `bcftools norm`). pos0 is 0-based."""
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


# --------------------------------------------------------------------------- planting

# relative weight of the six substitution classes for true somatic SNVs (pyrimidine reference)
SPECTRUM = {("C", "T"): 0.36, ("C", "A"): 0.16, ("C", "G"): 0.10,
            ("T", "C"): 0.20, ("T", "A"): 0.09, ("T", "G"): 0.09}
PAIR = {"A": "T", "C": "G", "G": "C", "T": "A"}


def _draw_alt(refb: str, rng) -> str:
    pyr = refb if refb in "CT" else PAIR[refb]
    opts = [(k, v) for k, v in SPECTRUM.items() if k[0] == pyr]
    p = np.array([v for _, v in opts])
    alt = opts[int(rng.choice(len(opts), p=p / p.sum()))][0][1]
    return alt if refb in "CT" else PAIR[alt]


def candidate_sites(ref: dict) -> np.ndarray:
    """Positions eligible for a planted (non-duplicated-segment) site: inside a target with full depth,
    +-60 bp free of soft-masked repeat and of any K-mer that occurs twice in the reference, >= 150 bp
    away from the duplicated segment, not next to a homopolymer of 5 or more."""
    seq, mask = ref["seq"], ref["mask"]
    uniq = unique_mask(seq)
    bad = np.concatenate([[0], np.cumsum(~uniq | mask)])
    ok = np.zeros(len(seq), dtype=bool)
    for s, e in ref["targets"]:
        ok[s + EDGE:e - EDGE] = True
    if ref["dup"]:
        for a, b in (ref["dup"]["a"], ref["dup"]["b"]):
            ok[max(0, a - 150):b + 150] = False
    cand = []
    for x in np.nonzero(ok)[0]:
        x = int(x)
        w = seq[x - 6:x + 7]
        if bad[x + 61] - bad[x - 60] == 0 and not any(b * 5 in w for b in BASES):
            cand.append(x)
    return np.array(cand)


class Planter:
    """Hands out candidate positions that keep a minimum distance from every site already planted."""

    def __init__(self, ref: dict, rng):
        self.seq, self.rng = ref["seq"], rng
        self.cand = ref["cand"]
        self.taken: list[int] = []

    def free(self, x: int, dist: int) -> bool:
        i = bisect.bisect_left(self.taken, x)
        return not ((i < len(self.taken) and self.taken[i] - x < dist) or (i > 0 and x - self.taken[i - 1] < dist))

    def take(self, x: int) -> None:
        bisect.insort(self.taken, x)

    def site(self, dist: int = 40, pred=None, pool=None) -> int:
        pool = self.cand if pool is None else pool
        for _ in range(20000):
            x = int(pool[int(self.rng.integers(0, len(pool)))])
            if self.free(x, dist) and (pred is None or pred(x)):
                self.take(x)
                return x
        raise RuntimeError("could not place a site")


def plant(ref: dict, seed: int, level: int) -> dict:
    """Decide genotypes, somatic variants and artefacts for three patients."""
    rng = np.random.default_rng(71_000 + seed)
    seq = ref["seq"]
    pl = Planter(ref, rng)
    patients = ["P1", "P2", "P3"]

    # ---- germline: population SNP sites with an allele frequency, genotypes drawn per haplotype
    n_pop = 110
    pop = []
    for _ in range(n_pop):
        x = pl.site(25)
        alt = _draw_alt(seq[x], rng)
        pop.append((x, alt, float(np.clip(rng.beta(0.7, 1.3), 0.06, 0.92))))
    while True:
        geno = {p: {} for p in patients}
        for x, alt, f in pop:
            for p in patients:
                h = (int(rng.random() < f), int(rng.random() < f))
                if h != (0, 0):
                    geno[p][x] = (alt, h)
        diff = [sum(1 for x, _, _ in pop if geno[a].get(x, (0, (0, 0)))[1] != geno[b].get(x, (0, (0, 0)))[1])
                for a, b in (("P1", "P2"), ("P1", "P3"), ("P2", "P3"))]
        if min(diff) >= 30 and min(len(geno[p]) for p in patients) >= 30:
            break

    # ---- orientation damage (level 3): which tumours, which chemistry
    n_dam = int(rng.integers(1, 3))
    damaged = sorted(rng.permutation(patients)[:n_dam].tolist()) if level >= 3 else []
    # reference base -> (alt base, read-pair orientation that shows it); 0 = F1R2, 1 = F2R1
    damage = {"deamination": {"C": ("T", 0), "G": ("A", 1)}, "oxidation": {"G": ("T", 0), "C": ("A", 1)}}

    dup = ref["dup"]
    out = {"patients": {}, "pop_sites": [(x, a) for x, a, _ in pop], "damage_mode": "deamination+oxidation" if damaged else None,
           "damaged": damaged}
    for p in patients:
        purity = float(rng.uniform(0.58, 0.90))
        ccf = {"clonal": 1.0, "sub1": float(rng.uniform(0.45, 0.65)), "sub2": float(rng.uniform(0.20, 0.34))}
        som = []
        plan = ([("clonal", "snv")] * int(rng.integers(13, 17)) + [("sub1", "snv")] * int(rng.integers(6, 9))
                + [("sub2", "snv")] * int(rng.integers(5, 8)) + [("low", "snv")] * 3)
        n_del, n_ins = int(rng.integers(3, 5)), int(rng.integers(2, 4))
        plan += [(str(rng.choice(["clonal", "sub1"])), "del") for _ in range(n_del)]
        plan += [(str(rng.choice(["clonal", "sub1"])), "ins") for _ in range(n_ins)]
        n_forced = 0
        for cls, kind in plan:
            vaf = float(rng.uniform(0.02, 0.04)) if cls == "low" else purity * ccf[cls] / 2
            if kind == "snv" and p in damaged and cls == "sub1" and n_forced < 5:
                # five true variants of the first subclone of a damaged tumour are of the damage substitution
                # classes: a blanket "drop these classes below an allele-fraction cut-off" then costs recall
                rule = damage["deamination" if n_forced % 2 == 0 else "oxidation"]
                x = pl.site(30, pred=lambda q: seq[q] in rule)
                som.append({"pos0": x, "ref": seq[x], "alt": rule[seq[x]][0], "kind": "snv", "cls": cls, "vaf": vaf})
                n_forced += 1
            elif kind == "snv":
                x = pl.site(30)
                som.append({"pos0": x, "ref": seq[x], "alt": _draw_alt(seq[x], rng), "kind": "snv",
                            "cls": cls, "vaf": vaf})
            elif kind == "del":
                d = int(rng.choice([1, 1, 2, 3, 4, 6, 8]))
                x = pl.site(60, pred=lambda q: seq[q] != seq[q + d] and seq[q + 1] != seq[q + d + 1])
                som.append({"pos0": x, "ref": seq[x:x + d + 1], "alt": seq[x], "kind": "del", "cls": cls, "vaf": vaf})
            else:
                k = int(rng.choice([1, 1, 2, 3, 5]))
                x = pl.site(60)
                while True:
                    ins = "".join(BASES[i] for i in rng.integers(0, 4, size=k))
                    if ins[-1] != seq[x] and ins[0] != seq[x + 1] and len(set(ins)) >= min(2, k):
                        break
                som.append({"pos0": x, "ref": seq[x], "alt": seq[x] + ins, "kind": "ins", "cls": cls, "vaf": vaf})
        for v in som:   # representation check: the planted form is already the normalised one
            assert left_normalise(seq, v["pos0"], v["ref"], v["alt"]) == (v["pos0"], v["ref"], v["alt"]), v

        # ---- clonal SNVs in the interior of the duplicated segment (either copy)
        if dup:
            lo, hi = DUP_MARGIN, dup["len"] - DUP_MARGIN
            offs: list[int] = []
            n_dup = int(rng.integers(9, 12))
            while len(offs) < n_dup:
                o = int(rng.integers(lo, hi))
                if all(abs(o - q) >= 40 for q in offs) and pl.free(dup["a"][0] + o, 40) and pl.free(dup["b"][0] + o, 40):
                    offs.append(o)
            for o in offs:
                copy = "a" if rng.random() < 0.5 else "b"
                other = "b" if copy == "a" else "a"
                x, m = dup[copy][0] + o, dup[other][0] + o
                pl.take(x)
                pl.take(m)
                som.append({"pos0": x, "ref": seq[x], "alt": _draw_alt(seq[x], rng), "kind": "snv", "cls": "dup",
                            "vaf": purity / 2, "mirror0": m})

        # ---- artefact sites of a damaged tumour library: both chemistries, so that the artefacts are not
        # one substitution class, and allele fractions that reach into the range of clonal variants
        arts = []
        if p in damaged:
            n_art = int(rng.integers(14, 19))
            n_deam = n_art // 2 + int(rng.integers(-2, 3))
            for chem in ["deamination"] * n_deam + ["oxidation"] * (n_art - n_deam):
                rule = damage[chem]
                x = pl.site(30, pred=lambda q: seq[q] in rule)
                alt, orient = rule[seq[x]]
                arts.append({"pos0": x, "ref": seq[x], "alt": alt, "orient": orient, "chemistry": chem,
                             "rate": float(rng.uniform(0.16, 0.56))})
        out["patients"][p] = {"purity": purity, "ccf": ccf, "germline": geno[p], "somatic": som, "artefacts": arts}
    return out


# --------------------------------------------------------------------------- reads

def _fragments(ref: dict, depth: float, eff: list[float], rng):
    L = len(ref["seq"])
    for (ts, te), ef in zip(ref["targets"], eff):
        n = int(round(depth * ef * (te - ts + 100) / (2 * RL)))
        mids = rng.integers(ts - 50, te + 50, size=n)
        lens = np.clip(rng.normal(INS_MEAN, INS_SD, size=n), INS_MIN, INS_MAX).astype(int)
        haps = rng.integers(0, 2, size=n)
        ori = rng.integers(0, 2, size=n)
        dupn = rng.choice(3, size=n, p=DUP_P)
        for i in range(n):
            s = max(0, int(mids[i]) - int(lens[i]) // 2)
            e = min(L, s + int(lens[i]))
            if e - s < INS_MIN:
                s = max(0, e - INS_MIN)
            yield s, e, int(haps[i]), int(ori[i]), int(dupn[i])


def simulate_sample(ref: dict, pat: dict, tumour: bool, depth: float, eff: list[float], rng):
    """Returns (r1 list, r2 list, stats) where stats[site index] = realised fragment counts."""
    seq = ref["seq"]
    sites = []   # (pos0, kind, payload)
    for x, (alt, h) in pat["germline"].items():
        sites.append((x, "g", (alt, h)))
    if tumour:
        for j, v in enumerate(pat["somatic"]):
            sites.append((v["pos0"], "s", j))
        for j, a in enumerate(pat["artefacts"]):
            sites.append((a["pos0"], "a", j))
    sites.sort()
    spos = [s[0] for s in sites]
    n_som, n_art = len(pat["somatic"]), len(pat["artefacts"])
    som_alt = np.zeros((n_som, 2), dtype=int)
    som_cov = np.zeros(n_som, dtype=int)
    art_alt = np.zeros(n_art, dtype=int)
    art_cov = np.zeros(n_art, dtype=int)
    r1s: list[str] = []
    r2s: list[str] = []
    n_frag = 0
    for s, e, hap, ori, dupn in _fragments(ref, depth, eff, rng):
        n_frag += 1
        lo, hi = bisect.bisect_left(spos, s), bisect.bisect_left(spos, e)
        if lo == hi:
            fs = seq[s:e]
        else:
            fl = list(seq[s:e])
            indels = []
            for x, kind, pay in sites[lo:hi]:
                covered = x < s + RL or x >= e - RL
                if kind == "g":
                    if pay[1][hap]:
                        fl[x - s] = pay[0]
                elif kind == "s":
                    v = pat["somatic"][pay]
                    carries = rng.random() < v["vaf"]
                    if v["kind"] != "snv" and x + len(v["ref"]) + 2 > e:
                        carries = False                     # deletion would run off the fragment
                    if covered:
                        som_cov[pay] += 1
                        if carries:
                            som_alt[pay, ori] += 1
                    if carries:
                        if v["kind"] == "snv":
                            fl[x - s] = v["alt"]
                        else:
                            indels.append(v)
                else:
                    a = pat["artefacts"][pay]
                    hit = ori == a["orient"] and rng.random() < a["rate"]
                    if covered:
                        art_cov[pay] += 1
                        if hit:
                            art_alt[pay] += 1
                    if hit:
                        fl[x - s] = a["alt"]
            for v in sorted(indels, key=lambda q: -q["pos0"]):
                i = v["pos0"] - s
                fl[i:i + len(v["ref"])] = list(v["alt"])
            fs = "".join(fl)
        left, right = fs[:RL], rc(fs[-RL:])
        a, b = (left, right) if ori == 0 else (right, left)
        for _ in range(1 + dupn):
            r1s.append(a)
            r2s.append(b)
    stats = {"som_alt": som_alt, "som_cov": som_cov, "art_alt": art_alt, "art_cov": art_cov, "fragments": n_frag}
    return r1s, r2s, stats


_LUT = np.zeros(256, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    _LUT[_b] = _i
_BASE_ARR = np.frombuffer(b"ACGT", dtype=np.uint8)


def _sequence(reads: list[str], rng):
    n = len(reads)
    arr = np.frombuffer("".join(reads).encode(), dtype=np.uint8).reshape(n, RL).copy()
    q = rng.choice(3, size=(n, RL), p=QUAL_P)
    err = rng.random((n, RL)) < np.array(QUAL_ERR)[q]
    k = int(err.sum())
    arr[err] = _BASE_ARR[(_LUT[arr[err]] + rng.integers(1, 4, size=k)) % 4]
    return arr, QUAL_CH[q]


def write_fastq_pair(p1: Path, p2: Path, r1s: list[str], r2s: list[str], index: str, rng) -> int:
    n = len(r1s)
    order = rng.permutation(n)
    r1s = [r1s[i] for i in order]
    r2s = [r2s[i] for i in order]
    tiles = np.sort(rng.integers(0, 78, size=n))
    tile_ids = [f"{1 + t // 39}{1 + (t % 39) // 13}{(t % 13) + 1:02d}" for t in range(78)]
    xs = rng.integers(1000, 32000, size=n)
    ys = np.zeros(n, dtype=int)
    for t in np.unique(tiles):
        m = tiles == t
        ys[m] = np.sort(rng.integers(1000, 36000, size=int(m.sum())))
    for path, reads, mate in ((p1, r1s, 1), (p2, r2s, 2)):
        arr, qual = _sequence(reads, rng)
        with open(path, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=6, mtime=0) as fh:
            buf = []
            for i in range(n):
                buf.append(b"@A01789:214:HKF3WDSX7:2:%s:%d:%d %d:N:0:%s\n" % (
                    tile_ids[tiles[i]].encode(), xs[i], ys[i], mate, index.encode()))
                buf.append(arr[i].tobytes())
                buf.append(b"\n+\n")
                buf.append(qual[i].tobytes())
                buf.append(b"\n")
                if len(buf) >= 50000:
                    fh.write(b"".join(buf))
                    buf = []
            fh.write(b"".join(buf))
    return n


def fasta(name: str, seq: str, width: int = 60) -> tuple[str, str]:
    body = "\n".join(seq[i:i + width] for i in range(0, len(seq), width))
    fa = f">{name}\n{body}\n"
    fai = f"{name}\t{len(seq)}\t{len(name) + 2}\t{width}\t{width + 1}\n"
    return fa, fai
