"""Seed-independent reference material for the trio_liftover family (downloaded once, cached).

Real sequence is kept only inside a few gene-dense windows of chr22; everything else is N-masked so
that coordinates stay real while the shipped files stay small. The same windows are used for every
seed; what varies with the seed is the panel, the causal genes and every variant.
"""
from __future__ import annotations

import gzip
import hashlib
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
CACHE = ROOT / ".cache" / "trio_liftover"
DL = CACHE / "dl"
DERIVED = CACHE / "derived_v1"

SOURCES = {
    "hg19.chr22.fa.gz": ("https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr22.fa.gz",
                         "4a5a0049b23232b163e575bf506b8af86e7d563fa9990311fb9fa3346d6865ec"),
    "hg38.chr22.fa.gz": ("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz",
                         "05f9d97d6fbfd08a44ca45b50837ca2ae9c471f35ba79dffec04d2cb5eaaf695"),
    "hg19ToHg38.over.chain.gz": ("https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
                                 "5c0598e500ceb5a78c73086929e8ef993aec309bcafb595139b53d440b125a1d"),
    "hg38ToHg19.over.chain.gz": ("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz",
                                 "14a712e8e147d9fc8e9d87d51977b46f6f8ddb93efbe5d0843d86b6205f587b1"),
    "Homo_sapiens.GRCh38.113.chromosome.22.gff3.gz": (
        "https://ftp.ensembl.org/pub/release-113/gff3/homo_sapiens/Homo_sapiens.GRCh38.113.chromosome.22.gff3.gz",
        "20b5eff9aa90c936d5cc9c6fa514c6a9678cfebe0aa3b3882b84337284ecc102"),
}

# 0-based half-open GRCh38 windows (5 x 600 kb); chosen for protein-coding gene density, for covering
# four different hg19->hg38 offsets (+12477, -354289, -342187/-396448, -396004, -438429) and for
# containing hg19 segments that have no GRCh38 alignment.
WINDOWS_HG38 = [(18_900_000, 19_500_000), (21_400_000, 22_000_000), (23_600_000, 24_200_000),
                (41_550_000, 42_150_000), (50_050_000, 50_650_000)]
MAIN_CHAIN_ID = "23"
HG19_LEN, HG38_LEN = 51_304_566, 50_818_468
TOOLS = ["bcftools=1.21", "samtools=1.21", "htslib=1.21", "bcftools-liftover-plugin=1.21", "picard-slim=3.3.0",
         "pysam=0.23.0"]


def _sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def ensure_downloads() -> None:
    DL.mkdir(parents=True, exist_ok=True)
    for name, (url, sha) in SOURCES.items():
        p = DL / name
        if not p.exists():
            subprocess.run(["curl", "-fsSL", "-o", str(p), url], check=True)
        if _sha256(p) != sha:
            raise RuntimeError(f"sha256 mismatch for cached download {p}")


def read_fa_gz(p: Path) -> str:
    with gzip.open(p, "rt") as fh:
        fh.readline()
        return "".join(line.strip() for line in fh).upper()


def parse_chains(p: Path, t_name: str = "chr22", q_name: str = "chr22"):
    """Return (list of chains, kept text). A chain = dict(header fields, blocks=[(t0, q0, size)])
    with t = source assembly, q = destination assembly (0-based, q on the q strand)."""
    chains, text, cur, keep = [], [], None, False
    with gzip.open(p, "rt") as fh:
        for line in fh:
            f = line.split()
            if not f:
                if keep:
                    text.append("\n")
                continue
            if f[0] == "chain":
                keep = f[2] == t_name and f[7] == q_name
                if keep:
                    cur = dict(id=f[12], tstrand=f[4], qstrand=f[9], qsize=int(f[8]), blocks=[])
                    cur["_t"], cur["_q"] = int(f[5]), int(f[10])
                    chains.append(cur)
                    text.append(line)
                continue
            if keep:
                text.append(line)
                size = int(f[0])
                cur["blocks"].append((cur["_t"], cur["_q"], size))
                if len(f) == 3:
                    cur["_t"] += size + int(f[1])
                    cur["_q"] += size + int(f[2])
    return chains, "".join(text)


def parse_genes(p: Path, windows) -> dict:
    """Protein-coding genes with an Ensembl-canonical mRNA lying fully inside a window."""
    genes, tx2gene = {}, {}
    with gzip.open(p, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            typ, s, e, strand, attr = f[2], int(f[3]), int(f[4]), f[6], f[8]
            if typ == "gene" and "biotype=protein_coding" in attr:
                gid = re.search(r"ID=gene:(\w+)", attr).group(1)
                nm = re.search(r"Name=([^;]+)", attr)
                if nm and any(a + 5000 <= s - 1 and e <= b - 5000 for a, b in windows):
                    genes[gid] = dict(id=gid, name=nm.group(1), strand=strand, start=s, end=e, exons=[], cds=[],
                                      lines=[f])
            elif typ == "mRNA" and "Ensembl_canonical" in attr and "biotype=protein_coding" in attr:
                gid = re.search(r"Parent=gene:(\w+)", attr).group(1)
                if gid in genes:
                    tid = re.search(r"ID=transcript:(\w+)", attr).group(1)
                    tx2gene[tid] = gid
                    genes[gid]["tx"] = tid
                    genes[gid]["lines"].append(f)
            elif typ in ("exon", "CDS", "five_prime_UTR", "three_prime_UTR"):
                m = re.search(r"Parent=transcript:(\w+)", attr)
                if m and m.group(1) in tx2gene:
                    g = genes[tx2gene[m.group(1)]]
                    g["lines"].append(f)
                    if typ == "exon":
                        g["exons"].append((s - 1, e))
                    elif typ == "CDS":
                        g["cds"].append((s - 1, e))
    out = {}
    for gid, g in genes.items():
        if "tx" in g and g["cds"]:
            g["exons"].sort()
            g["cds"].sort()
            out[gid] = g
    return out


class Ref:
    """Loaded reference material (masked sequences, chain blocks, gene models)."""

    def __init__(self):
        ensure_downloads()
        DERIVED.mkdir(parents=True, exist_ok=True)
        fwd_chains, fwd_text = parse_chains(DL / "hg19ToHg38.over.chain.gz")
        _, rev_text = parse_chains(DL / "hg38ToHg19.over.chain.gz")
        main = [c for c in fwd_chains if c["id"] == MAIN_CHAIN_ID][0]
        assert main["qstrand"] == "+" and main["tstrand"] == "+"
        self.blocks = main["blocks"]                      # (hg19_0, hg38_0, size), sorted in both
        # every hg19 interval covered by ANY shipped chain (for choosing hg19-only positions)
        cov = sorted((t, t + s) for c in fwd_chains for (t, q, s) in c["blocks"])
        self.covered = cov
        self.windows38 = list(WINDOWS_HG38)
        self.windows19 = []
        for a, b in self.windows38:
            lo = min(t + max(0, a - q) for t, q, s in self.blocks if q < b and q + s > a)
            hi = max(t + min(s, b - q) for t, q, s in self.blocks if q < b and q + s > a)
            self.windows19.append((lo, hi))
        hg19 = read_fa_gz(DL / "hg19.chr22.fa.gz")
        hg38 = read_fa_gz(DL / "hg38.chr22.fa.gz")
        assert len(hg19) == HG19_LEN and len(hg38) == HG38_LEN
        self.hg19 = self._mask(hg19, self.windows19)
        self.hg38 = self._mask(hg38, self.windows38)
        self.genes = parse_genes(DL / "Homo_sapiens.GRCh38.113.chromosome.22.gff3.gz", self.windows38)
        self.fwd_chain_text, self.rev_chain_text = fwd_text, rev_text
        self._starts38 = [q for _, q, _ in self.blocks]
        self._starts19 = [t for t, _, _ in self.blocks]

    @staticmethod
    def _mask(seq: str, windows) -> str:
        parts, last = [], 0
        for a, b in sorted(windows):
            parts.append("N" * (a - last))
            parts.append(seq[a:b])
            last = b
        parts.append("N" * (len(seq) - last))
        return "".join(parts)

    # ---- coordinate conversion through the main chain -----------------------------------------
    def to19(self, p38: int):
        """hg38 0-based -> (hg19 0-based, distance to the nearest block edge) or None."""
        import bisect
        i = bisect.bisect_right(self._starts38, p38) - 1
        if i < 0:
            return None
        t, q, s = self.blocks[i]
        if p38 >= q + s:
            return None
        return t + (p38 - q), min(p38 - q, q + s - 1 - p38)

    def to38(self, p19: int):
        import bisect
        i = bisect.bisect_right(self._starts19, p19) - 1
        if i < 0:
            return None
        t, q, s = self.blocks[i]
        if p19 >= t + s:
            return None
        return q + (p19 - t), min(p19 - t, t + s - 1 - p19)

    def is_covered19(self, p19: int) -> bool:
        import bisect
        i = bisect.bisect_right(self.covered, (p19, 1 << 60)) - 1
        return i >= 0 and self.covered[i][0] <= p19 < self.covered[i][1]

    def flank_identical(self, p38: int, p19: int, left: int, right: int) -> bool:
        a, b = self.hg38[p38 - left:p38 + right], self.hg19[p19 - left:p19 + right]
        return a == b and "N" not in a

    # ---- shipped resource files (plain text; compressed/indexed later inside the tool image) ---
    def gff3_text(self) -> str:
        out = ["##gff-version 3", f"##sequence-region   chr22 1 {HG38_LEN}",
               "#!genome-build GRCh38.p14", "#!note one canonical protein-coding transcript per gene (lab transcript set)"]
        for g in sorted(self.genes.values(), key=lambda g: g["start"]):
            out.append("###")
            for f in g["lines"]:
                out.append("\t".join(["chr22"] + f[1:]))
        return "\n".join(out) + "\n"

    def write_plain(self, d: Path) -> None:
        d.mkdir(parents=True, exist_ok=True)
        for name, seq in (("GRCh37.chr22.fa", self.hg19), ("GRCh38.chr22.fa", self.hg38)):
            with open(d / name, "w") as fh:
                fh.write(">chr22\n")
                for i in range(0, len(seq), 60):
                    fh.write(seq[i:i + 60] + "\n")
        (d / "hg19ToHg38.chr22.over.chain").write_text(self.fwd_chain_text)
        (d / "hg38ToHg19.chr22.over.chain").write_text(self.rev_chain_text)
        (d / "GRCh38.chr22.canonical.gff3").write_text(self.gff3_text())
