"""Forward simulator for the rnaseq_confound family (host side, numpy only).

Order of events (truth by construction):
  1. pick a gene subset from real Ensembl transcripts (sequence realism only; no real expression),
  2. decide every effect (treatment log2FC, library-prep batch log2FC, donor sex) and the design,
  3. draw negative-binomial fragment counts per gene and library  -> hidden truth,
  4. turn every fragment into a single-end read sampled from the gene's transcripts.
Nothing in here runs an analysis tool.
"""
from __future__ import annotations

import gzip
import hashlib
import re
import subprocess
from pathlib import Path

import numpy as np

CACHE = Path(__file__).resolve().parents[4] / ".cache" / "rnaseq_confound"
ENSEMBL_RELEASE = 113
SOURCES = {
    "Homo_sapiens.GRCh38.cdna.all.fa.gz": (
        "https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz",
        "34f848e2dd9c2a4e30d6ff2c7918a3e06c51fe1716c1e955e71e0c36ce28d5ad"),
    "Homo_sapiens.GRCh38.ncrna.fa.gz": (
        "https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz",
        "3a5e436d438d17af848492011c361848b59aab85f169fbe631d540fa07aa57b5"),
}
Y_GENES = ["RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY"]
XIST_TX = "ENST00000429829"          # Ensembl canonical XIST transcript (19 kb)
ELIGIBLE = "eligible_v2.fa.gz"       # derived cache: <= 3 transcripts per eligible gene + sex genes
READLEN = 50
K = 31
COMP = str.maketrans("ACGT", "TGCA")


# --------------------------------------------------------------------------- reference material

def _sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _fetch(name: str) -> Path:
    url, sha = SOURCES[name]
    p = CACHE / name
    if not p.exists():
        CACHE.mkdir(parents=True, exist_ok=True)
        tmp = p.with_suffix(".part")
        subprocess.run(["curl", "-fsSL", "-o", str(tmp), url], check=True)
        tmp.rename(p)
    if _sha256(p) != sha:
        raise RuntimeError(f"{p} does not match the recorded sha256; delete it and retry")
    return p


def _fasta(path: Path):
    name, seq = None, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line[0] == ">":
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].strip(), []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def _build_eligible() -> Path:
    """Reduce the two Ensembl FASTA files to the transcripts a task may use.

    Autosomal/X protein-coding genes: protein-coding transcripts on primary chromosomes, 800-5000 nt,
    A/C/G/T only; the three longest per gene are kept (longest first). Sex genes: the longest
    protein-coding transcript <= 5000 nt of each Y gene, and the canonical XIST transcript.
    """
    out = CACHE / ELIGIBLE
    if out.exists():
        return out
    genes: dict[str, dict] = {}
    sex: dict[str, tuple] = {}
    for h, s in _fasta(_fetch("Homo_sapiens.GRCh38.cdna.all.fa.gz")):
        f = h.split()
        m = re.search(r"gene_symbol:(\S+)", h)
        if not m or "transcript_biotype:protein_coding" not in f or "gene_biotype:protein_coding" not in f:
            continue
        sym = m.group(1)
        loc = re.match(r"chromosome:GRCh38:(\d+|X|Y):", f[2])
        if not loc or set(s) - set("ACGT"):
            continue
        gid = re.search(r"gene:(\S+)", h).group(1).split(".")[0]
        tid = f[0].split(".")[0]
        if loc.group(1) == "Y":
            if sym in Y_GENES and len(s) <= 5000 and (sym not in sex or len(s) > len(sex[sym][3])):
                sex[sym] = (tid, gid, "Y", s)
            continue
        if not 800 <= len(s) <= 5000 or not re.fullmatch(r"[A-Za-z0-9-]+", sym):
            continue
        genes.setdefault(gid, {"sym": sym, "chrom": loc.group(1), "tx": []})["tx"].append((tid, s))
    for h, s in _fasta(_fetch("Homo_sapiens.GRCh38.ncrna.fa.gz")):
        if h.startswith(XIST_TX + "."):
            sex["XIST"] = (XIST_TX, re.search(r"gene:(\S+)", h).group(1).split(".")[0], "X", s)
    if set(sex) != set(Y_GENES) | {"XIST"}:
        raise RuntimeError(f"sex genes missing from the Ensembl files: {sorted(sex)}")
    tmp = out.with_suffix(".part")
    with gzip.open(tmp, "wt", compresslevel=6) as fh:
        for sym in ["XIST"] + Y_GENES:
            tid, gid, chrom, s = sex[sym]
            fh.write(f">{tid}|{gid}|{sym}|{chrom}|sex\n{s}\n")
        for gid in sorted(genes):
            g = genes[gid]
            seen, kept = set(), []
            for tid, s in sorted(g["tx"], key=lambda t: (-len(t[1]), t[0])):
                if s not in seen:
                    seen.add(s)
                    kept.append((tid, s))
            for tid, s in kept[:3]:
                fh.write(f">{tid}|{gid}|{g['sym']}|{g['chrom']}|auto\n{s}\n")
    tmp.rename(out)
    return out


_ELIGIBLE_CACHE: list = []


def load_eligible() -> tuple[list[dict], list[dict]]:
    """-> (sex_genes, other_genes); each gene = {gene_id, symbol, chrom, tx: [(tid, seq), ...]}."""
    if _ELIGIBLE_CACHE:
        return _ELIGIBLE_CACHE[0]
    sex: dict[str, dict] = {}
    oth: dict[str, dict] = {}
    for h, s in _fasta(_build_eligible()):
        tid, gid, sym, chrom, kind = h.split("|")
        d = sex if kind == "sex" else oth
        d.setdefault(gid, {"gene_id": gid, "symbol": sym, "chrom": chrom, "tx": []})["tx"].append((tid, s))
    _ELIGIBLE_CACHE.append((list(sex.values()), list(oth.values())))
    return _ELIGIBLE_CACHE[0]


_CODE = np.full(256, 0, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    _CODE[_b] = _i
_POW = np.uint64(4) ** np.arange(K - 1, -1, -1, dtype=np.uint64)


def _kmers(seq: str) -> set[int]:
    """Exact 2-bit codes of the K-mers of seq and of its reverse complement (strand-agnostic)."""
    c = _CODE[np.frombuffer(seq.encode(), dtype=np.uint8)]
    out: set[int] = set()
    for strand in (c, (3 - c[::-1]).astype(np.uint8)):
        win = np.lib.stride_tricks.sliding_window_view(strand, K)
        out.update((win.astype(np.uint64) * _POW).sum(axis=1).tolist())
    return out


def select_genes(rng: np.random.Generator, n_other: int, frac_multi: float = 0.2,
                 max_shared: float = 0.01) -> list[dict]:
    """Seeded gene subset. A candidate is rejected when more than `max_shared` of the 31-mers of
    any of its transcripts already occur (either strand) in an accepted gene, so that a read can
    almost always be assigned to one gene and "fragments per gene" is a well-defined truth."""
    sex, other = load_eligible()
    chosen = []
    seen: set[int] = set()
    for g in sex:
        chosen.append(dict(g, kind="sex"))
        for _, s in g["tx"]:
            seen |= _kmers(s)
    for i in rng.permutation(len(other)):
        g = other[int(i)]
        n_tx = min(len(g["tx"]), int(rng.integers(2, 4))) if rng.random() < frac_multi else 1
        tx = g["tx"][:n_tx]
        ks = [_kmers(s) for _, s in tx]
        if any(len(k & seen) > max_shared * len(k) for k in ks):
            continue
        for k in ks:
            seen |= k
        chosen.append(dict(g, tx=tx, kind="auto"))
        if len(chosen) == n_other + len(sex):
            break
    order = rng.permutation(len(chosen))
    return [chosen[int(i)] for i in order]


# --------------------------------------------------------------------------- design and counts

def make_design(rng: np.random.Generator, level: int, n_per_group: int = 6, major: int = 5) -> dict:
    """Libraries (true identities) and the sample sheet rows they are filed under.

    level 1: one library-prep batch.  level >= 2: two batches, `major`/`n-major` treated vs the
    mirror image for vehicle.  level 3: the FASTQ files of one treated and one vehicle library of
    opposite donor sex, taken from the two well-populated cells, are filed under each other's id.
    """
    n = 2 * n_per_group
    cond = np.array(["treated"] * n_per_group + ["vehicle"] * n_per_group)
    batch = np.zeros(n, dtype=int)
    if level >= 2:
        flip = int(rng.integers(0, 2))
        batch[:major] = flip                       # treated majority batch
        batch[major:n_per_group] = 1 - flip
        batch[n_per_group:n_per_group + major] = 1 - flip   # vehicle majority batch
        batch[n_per_group + major:] = flip
    sex = np.empty(n, dtype="<U1")
    for lo in (0, n_per_group):
        half = n_per_group // 2
        s = np.array(["F"] * half + ["M"] * (n_per_group - half))
        while True:
            s = s[rng.permutation(n_per_group)]
            if level < 2 or len(set(s[:major])) == 2:   # both sexes inside the majority cell
                break
        sex[lo:lo + n_per_group] = s
    ids = [f"RS{k:02d}" for k in rng.permutation(np.arange(1, 100))[:n]]
    donors = [f"D{k:03d}" for k in rng.permutation(np.arange(100, 1000))[:n]]
    order = rng.permutation(n)
    libs = [{"lib": int(i), "condition": str(cond[i]), "batch": int(batch[i]), "sex": str(sex[i]),
             "donor": donors[i], "sheet_id": ids[i]} for i in order]
    libs.sort(key=lambda d: d["sheet_id"])
    swapped: list[str] = []
    if level >= 3:
        t = [d for d in libs if d["condition"] == "treated" and d["batch"] == batch[0]]
        a = t[int(rng.integers(0, len(t)))]
        v = [d for d in libs if d["condition"] == "vehicle" and d["batch"] == batch[n_per_group]
             and d["sex"] != a["sex"]]
        b = v[int(rng.integers(0, len(v)))]
        swapped = sorted([a["sheet_id"], b["sheet_id"]])
    # file_id = the sample id whose FASTQ file holds this library
    for d in libs:
        d["file_id"] = d["sheet_id"]
    if swapped:
        a = next(d for d in libs if d["sheet_id"] == swapped[0])
        b = next(d for d in libs if d["sheet_id"] == swapped[1])
        a["file_id"], b["file_id"] = b["sheet_id"], a["sheet_id"]
    return {"libs": libs, "swapped": swapped}


def plant_effects(rng: np.random.Generator, genes: list[dict], level: int, libs: list[dict],
                  target_reads: int, frac_de: float = 0.10, frac_batch: float = 0.30,
                  de_min_reads: float = 0.0, lfc_range: tuple = (1.0, 3.0),
                  batch_range: tuple = (0.5, 2.0), disp_median: float = 0.03) -> dict:
    """Decide every effect. `base` is the expected read count of a gene in a vehicle library of the
    reference batch at size factor 1; it is rescaled so that the mean expected library size over
    the design equals `target_reads`. Treatment effects are planted only on genes whose expected
    count is at least `de_min_reads` (a responder nobody can detect at this depth would make the
    sensitivity score a lottery)."""
    G = len(genes)
    is_sex = np.array([g["kind"] == "sex" for g in genes])
    base = np.exp(rng.normal(0.0, 0.85, size=G))
    base = np.clip(base / np.median(base), 0.18, 12.0)
    base[is_sex] = rng.uniform(0.8, 2.5, size=int(is_sex.sum()))
    base = base / base.sum() * target_reads
    auto = np.flatnonzero(~is_sex)
    ok = auto[base[auto] >= de_min_reads]
    n_de = int(round(frac_de * len(auto)))
    de = rng.permutation(ok)[:n_de]
    lfc = np.zeros(G)
    mag = rng.uniform(lfc_range[0], lfc_range[1], size=n_de)
    sign = np.where(np.arange(n_de) % 2 == 0, 1.0, -1.0)[rng.permutation(n_de)]
    lfc[de] = np.round(mag * sign, 3)
    blfc = np.zeros(G)
    if level >= 2:
        de_set = set(de.tolist())
        for pool in (de, np.array([i for i in auto if i not in de_set])):
            nb = int(round(frac_batch * len(pool)))
            pick = rng.permutation(pool)[:nb]
            bm = rng.uniform(batch_range[0], batch_range[1], size=nb)
            bs = np.where(np.arange(nb) % 2 == 0, 1.0, -1.0)[rng.permutation(nb)]
            blfc[pick] = np.round(bm * bs, 3)
    disp = np.clip(np.exp(rng.normal(np.log(disp_median), 0.4, size=G)), disp_median / 3, disp_median * 3.3)
    # rescale to the read budget (up-regulation adds more reads than down-regulation removes)
    tot = np.mean([np.sum(base * 2.0 ** (lfc * (d["condition"] == "treated") + blfc * d["batch"])) for d in libs])
    base = base * target_reads / tot
    return {"base": base, "lfc": lfc, "batch_lfc": blfc, "disp": disp, "is_sex": is_sex}


def draw_counts(rng: np.random.Generator, genes: list[dict], eff: dict, libs: list[dict]) -> np.ndarray:
    """Negative-binomial fragment counts, genes x libraries (columns in the order of `libs`)."""
    G = len(genes)
    C = np.zeros((G, len(libs)), dtype=np.int64)
    for j, d in enumerate(libs):
        size = rng.uniform(0.85, 1.15)
        mu = eff["base"] * size * 2.0 ** (eff["lfc"] * (d["condition"] == "treated")
                                          + eff["batch_lfc"] * d["batch"])
        for i, g in enumerate(genes):
            if g["kind"] == "sex":
                if g["symbol"] == "XIST":
                    mu[i] = mu[i] if d["sex"] == "F" else 0.4
                else:
                    mu[i] = mu[i] if d["sex"] == "M" else 0.0
        r = 1.0 / eff["disp"]
        lam = rng.gamma(shape=r, scale=mu / r)
        C[:, j] = rng.poisson(lam)
    return C


# --------------------------------------------------------------------------- reads

class Transcriptome:
    def __init__(self, genes: list[dict], rng: np.random.Generator):
        code = np.full(256, 255, dtype=np.uint8)
        for k, b in enumerate(b"ACGT"):
            code[b] = k
        self.tx_ids, self.tx_gene, lens, seqs, self.iso_frac = [], [], [], [], []
        for gi, g in enumerate(genes):
            n = len(g["tx"])
            w = np.array([1.0]) if n == 1 else rng.dirichlet(np.full(n, 2.0))
            self.iso_frac.append(np.sort(w)[::-1])
            for tid, s in g["tx"]:
                self.tx_ids.append(tid)
                self.tx_gene.append(gi)
                lens.append(len(s))
                seqs.append(code[np.frombuffer(s.encode(), dtype=np.uint8)])
        self.lens = np.array(lens)
        self.offsets = np.concatenate([[0], np.cumsum(self.lens)[:-1]])
        self.seq = np.concatenate(seqs)
        self.tx_gene = np.array(self.tx_gene)
        self.gene_tx = [np.flatnonzero(self.tx_gene == gi) for gi in range(len(genes))]

    def reads(self, gene_counts: np.ndarray, rng: np.random.Generator, err: float) -> np.ndarray:
        """(n, READLEN) uint8 array of ASCII bases: one read per fragment, random strand."""
        per_tx = np.zeros(len(self.lens), dtype=np.int64)
        for gi, c in enumerate(gene_counts):
            tx = self.gene_tx[gi]
            per_tx[tx] = c if len(tx) == 1 else rng.multinomial(int(c), self.iso_frac[gi])
        t = np.repeat(np.arange(len(self.lens)), per_tx)
        n = len(t)
        start = self.offsets[t] + np.floor(rng.random(n) * (self.lens[t] - READLEN + 1)).astype(np.int64)
        r = self.seq[start[:, None] + np.arange(READLEN)[None, :]]
        rev = rng.random(n) < 0.5
        r[rev] = 3 - r[rev][:, ::-1]
        mask = rng.random((n, READLEN)) < err
        r[mask] = (r[mask] + rng.integers(1, 4, size=int(mask.sum()))) % 4
        r = r[rng.permutation(n)]
        return np.frombuffer(b"ACGT", dtype=np.uint8)[r]


def write_fastq(path: Path, reads: np.ndarray) -> None:
    qual = b"I" * READLEN
    raw = reads.tobytes()
    parts = [b"@%d\n%s\n+\n%s\n" % (i + 1, raw[i * READLEN:(i + 1) * READLEN], qual)
             for i in range(reads.shape[0])]
    with open(path, "wb") as fh:
        # mtime=0 and no file name in the header -> byte-identical regeneration
        with gzip.GzipFile(filename="", mode="wb", fileobj=fh, compresslevel=9, mtime=0) as gz:
            gz.write(b"".join(parts))
