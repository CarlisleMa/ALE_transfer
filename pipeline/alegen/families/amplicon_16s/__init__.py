"""Family: 16S rRNA V4 amplicon reads -> exact sequence variants, abundances, taxonomy, contaminants.

Truth-by-construction: the generator picks real V4 amplicons (NCBI RefSeq targeted loci, bacteria
16S; cut out with 515F/806R), decides which of them are specimen strains, mock strains and reagent
contaminants and how many read pairs each contributes to each sample, and only then emits reads.
The reference is the planted sequence set and the per-sample read-origin counts, never the output
of a denoiser.

Decision points (level 2 adds the first three, level 3 adds the last two):
  primer_spacers          0-3 nt heterogeneity spacers in front of degenerate primers
  mixed_orientation       ligation library: about half of the pairs start with the reverse primer
  chimeras                5-8 % of the read pairs are PCR recombinants of two strains of the sample
  strain_variants         strain pairs whose amplicons differ at 1-2 positions and whose ratio
                          changes from specimen to specimen (97 % OTUs merge them)
  reagent_contaminants    reagent taxa in every library, dominating the blanks and the low-biomass
                          specimens, next to cross-talk from real specimens into the blanks
Level 7 (level 3 plus M_control_swap, the reads of one blank and of one specimen filed under each other's
name with the sample sheet unchanged) is in l7.py.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import re
import subprocess
from pathlib import Path

import numpy as np

from ... import dockerlib
from ...spec import Attack, Probe, TaskBuild, Trap

ASSETS = Path(__file__).parent / "assets"
CACHE = Path(__file__).resolve().parents[4] / ".cache" / "amplicon_16s"
REFSEQ_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz"
REFSEQ_SHA256 = "c5dc51b877cfed87b526b5cb9425f05d59ce96805ce4d4eb2015d65d75436e7e"

FWD, REV = "GTGYCAGCMGCCGCGGTAA", "GGACTACNVGGGTWTCTAAT"      # 515F (Parada) / 806R (Apprill)
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "M": "AC", "K": "GT", "S": "CG",
         "W": "AT", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
COMP = str.maketrans("ACGTRYMKSWBDHVN", "TGCAYRKMSWVHDBN")
READLEN = 250
K = 8
CONDA = ["vsearch=2.30.1", "cutadapt=5.2"]
PASS = 0.9


def _rc(s: str) -> str:
    return s.translate(COMP)[::-1]


def _expand(primer: str) -> list[str]:
    out = [""]
    for c in primer:
        out = [o + b for o in out for b in IUPAC[c]]
    return out


def _regex(primer: str) -> re.Pattern:
    return re.compile("".join(b if len(IUPAC[b]) == 1 else f"[{IUPAC[b]}]" for b in primer))


# ----------------------------------------------------------------------------- reference pool

def _load_pool() -> list[tuple[str, str]]:
    """Unique V4 amplicons (primer-free, 251-255 nt, ACGT only) that belong to exactly one genus,
    as a list of (sequence, genus) sorted by sequence. Cached as v4_pool.tsv."""
    CACHE.mkdir(parents=True, exist_ok=True)
    pool_f = CACHE / "v4_pool.tsv"
    if pool_f.exists():
        return [tuple(l.split("\t")) for l in pool_f.read_text().splitlines()]
    fn = CACHE / "bacteria.16SrRNA.fna.gz"
    if not fn.exists():
        subprocess.run(["curl", "-fsSL", "-o", str(fn), REFSEQ_URL], check=True)
    digest = hashlib.sha256(fn.read_bytes()).hexdigest()
    if digest != REFSEQ_SHA256:
        raise RuntimeError(f"{fn} sha256 {digest} differs from the pinned release {REFSEQ_SHA256}")
    rf, rr = _regex(FWD), _regex(_rc(REV))
    by_seq: dict[str, set[str]] = {}
    name, chunks = None, []

    def flush():
        if name is None:
            return
        s = "".join(chunks)
        m = rf.search(s)
        m2 = rr.search(s, m.end()) if m else None
        if not m2:
            return
        amp = s[m.end():m2.start()]
        toks = name.split()
        if not 251 <= len(amp) <= 255 or set(amp) - set("ACGT") or len(toks) < 3:
            return
        if not re.fullmatch(r"[A-Z][a-z]+", toks[1]) or toks[1] == "Candidatus":
            return
        by_seq.setdefault(amp, set()).add(toks[1])

    with gzip.open(fn, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                name, chunks = line[1:].strip(), []
            else:
                chunks.append(line.strip().upper())
    flush()
    pool = sorted((s, next(iter(g))) for s, g in by_seq.items() if len(g) == 1)
    pool_f.write_text("".join(f"{s}\t{g}\n" for s, g in pool))
    return pool


def _close_pairs(pool: list[tuple[str, str]]) -> list[tuple[int, int, int]]:
    """Index pairs of same-genus, same-length amplicons at Hamming distance 1 or 2."""
    f = CACHE / "v4_close_pairs.json"
    if f.exists():
        return [tuple(x) for x in json.loads(f.read_text())]
    bylen: dict[int, list[int]] = {}
    for i, (s, _) in enumerate(pool):
        bylen.setdefault(len(s), []).append(i)
    out = set()
    for L, idx in sorted(bylen.items()):
        cuts = [0, L // 3, 2 * L // 3, L]
        for k in range(3):
            buckets: dict[str, list[int]] = {}
            for i in idx:
                buckets.setdefault(pool[i][0][cuts[k]:cuts[k + 1]], []).append(i)
            for v in buckets.values():
                if not 2 <= len(v) <= 300:
                    continue
                for x in range(len(v)):
                    for y in range(x + 1, len(v)):
                        a, b = v[x], v[y]
                        if pool[a][1] != pool[b][1]:
                            continue
                        d = sum(1 for p, q in zip(pool[a][0], pool[b][0]) if p != q)
                        if 1 <= d <= 2:
                            out.add((a, b, d))
    res = sorted(out)
    f.write_text(json.dumps(res))
    return res


def _kset(s: str) -> frozenset:
    return frozenset(s[i:i + K] for i in range(len(s) - K + 1))


def _maybe_within(ka: frozenset, kb: frozenset, edits: int) -> bool:
    """False only if the edit distance is provably > edits (each edit destroys at most K k-mers)."""
    return len(ka & kb) >= len(ka) - K * edits


def _is_bimera(q: str, others: list[str]) -> bool:
    best_p = best_s = 0
    for o in others:
        n = min(len(q), len(o))
        p = 0
        while p < n and q[p] == o[p]:
            p += 1
        s = 0
        while s < n and q[-1 - s] == o[-1 - s]:
            s += 1
        best_p, best_s = max(best_p, p), max(best_s, s)
    return best_p + best_s >= len(q)


# ----------------------------------------------------------------------------- read simulation

def _chimera(a: str, b: str, rng: np.random.Generator, cache: dict) -> str | None:
    """Recombinant: 5' part of a, 3' part of b, switching inside a shared 12-mer."""
    key = (a, b)
    if key not in cache:
        pos_b: dict[str, list[int]] = {}
        for j in range(len(b) - 11):
            pos_b.setdefault(b[j:j + 12], []).append(j)
        anchors = []
        for i in range(15, len(a) - 27):
            for j in pos_b.get(a[i:i + 12], ()):
                if abs(i - j) <= 6:
                    anchors.append((i, j))
        cache[key] = anchors
    anchors = cache[key]
    if not anchors:
        return None
    i, j = anchors[int(rng.integers(0, len(anchors)))]
    return a[:i] + b[j:]


def _write_pairs(path1: Path, path2: Path, templates: list[str], rev_frac: float, sp_f: list[str],
                 sp_r: list[str], rng: np.random.Generator, tag: str) -> None:
    n = len(templates)
    fvar, rvar = _expand(FWD), _expand(_rc(REV))
    fi = rng.integers(0, len(fvar), size=n)
    ri = rng.integers(0, len(rvar), size=n)
    sfi = rng.integers(0, len(sp_f), size=n)
    sri = rng.integers(0, len(sp_r), size=n)
    flip = rng.random(n) < rev_frac
    r1, r2 = [], []
    for k, t in enumerate(templates):
        top = sp_f[sfi[k]] + fvar[fi[k]] + t + rvar[ri[k]] + _rc(sp_r[sri[k]])
        a, b = top[:READLEN], _rc(top)[:READLEN]
        if flip[k]:
            a, b = b, a
        r1.append(a)
        r2.append(b)
    pos = np.arange(READLEN) / (READLEN - 1)
    lut = np.zeros(256, dtype=np.int8)
    for i, c in enumerate(b"ACGT"):
        lut[c] = i
    bases = np.frombuffer(b"ACGT", dtype=np.uint8)
    qchar = np.array([45, 58, 70], dtype=np.uint8)            # Q12 '-', Q25 ':', Q37 'F'
    perr = np.array([10 ** -1.2, 10 ** -2.5, 10 ** -3.7])
    tiles = rng.integers(1101, 1120, size=n)
    xs = rng.integers(1000, 28000, size=n)
    ys = rng.integers(1000, 28000, size=n)
    for path, reads, mate, scale in ((path1, r1, 1, 1.0), (path2, r2, 2, 1.6)):
        arr = np.frombuffer("".join(reads).encode(), dtype=np.uint8).reshape(n, READLEN)
        p12 = 0.002 + 0.10 * scale * pos ** 4
        p25 = 0.03 + 0.30 * scale * pos ** 2
        u = rng.random((n, READLEN))
        q = np.where(u < p12, 0, np.where(u < p12 + p25, 1, 2))
        err = rng.random((n, READLEN)) < perr[q]
        code = lut[arr]
        code = np.where(err, (code + rng.integers(1, 4, size=(n, READLEN))) % 4, code)
        seqs = bases[code]
        quals = qchar[q]
        # mtime=0: the same seed must give byte-identical files (gzip stores a timestamp otherwise)
        with gzip.GzipFile(filename=str(path), mode="wb", compresslevel=6, mtime=0) as fh:
            buf = []
            for k in range(n):
                buf.append(b"@M06142:%s:%d:%d:%d %d:N:0:1\n" % (tag.encode(), tiles[k], xs[k], ys[k], mate))
                buf.append(seqs[k].tobytes())
                buf.append(b"\n+\n")
                buf.append(quals[k].tobytes())
                buf.append(b"\n")
            fh.write(b"".join(buf))


# ----------------------------------------------------------------------------- texts

INSTRUCTION = """\
You are the bioinformatician of a microbiome sequencing core. A study batch of 16S rRNA gene (V4
region) amplicon libraries has come off the sequencer and the investigator needs the sequence
variants, their abundances, their taxonomy and the list of reagent contaminants.

Task directory:
- `/task`

Inputs (read-only) under `/task/input`:
- `fastq/<sample_id>_R1.fastq.gz`, `fastq/<sample_id>_R2.fastq.gz` - demultiplexed paired-end reads
  of {n_samples} libraries
- `samples.tsv` - sample sheet (`sample_id`, `sample_type` = study / negative_control / mock,
  `specimen`, `qpcr_16s_copies_per_ul`, `fastq_r1`, `fastq_r2`)
- `primers.tsv` - the amplification primers (IUPAC notation)
- `reference_db.fasta`, `reference_taxonomy.tsv` - the core's 16S V4 reference set (primer-free
  amplicon region) and the genus of every reference sequence
- `mock_community.tsv` - nominal composition of the mock community standard
- `library_prep_notes.md` - wet-lab notes for this batch
- `sop.md` - the core's definitions of the deliverables

Available software (on `PATH`): `vsearch`, `cutadapt`, Python 3.12 with numpy, scipy and pandas.
There is no internet access.

Your task:
1. Determine the amplicon sequence variants (ASVs) present in the batch, as defined in `sop.md`.
2. Count, for every sample, the read pairs that derive from each ASV.
3. Assign a genus to every ASV with the rule in `sop.md`.
4. Identify the ASVs that are reagent contaminants.

Required outputs under `/task/output`:
- `asv.fasta` - one record per ASV; the header line is `>` followed by an ASV id of your choice
  (unique, no whitespace); the sequence is the part of the amplicon that lies between the two
  primers (primers excluded), written 5'->3' in the direction that runs from the 515F primer to
  the 806R primer, using only A, C, G, T.
- `asv_table.tsv` - tab-separated, header `asv_id` followed by one column per `sample_id` of
  `samples.tsv` (all {n_samples} samples). One row per ASV of `asv.fasta`. Each value is the integer
  number of read pairs of that sample that derive from that ASV; read pairs that carry sequencing
  errors still count for the ASV they derive from, read pairs that derive from no ASV are not
  counted. Reagent-contaminant ASVs keep their rows and counts.
- `taxonomy.tsv` - tab-separated with at least the columns `asv_id`, `genus`; one row per ASV.
  `genus` is a genus name from `reference_taxonomy.tsv` or `Unclassified`.
- `contaminants.txt` - first line `asv_id`, then one ASV id per line for every ASV that is a reagent
  contaminant (only the first line if there is none). The investigator computes the community
  profile of each specimen from `asv_table.tsv` after removing the ASVs listed here.

Constraints:
- Do not modify anything under `/task/input`; write deliverables only under `/task/output`.
- Every sample of the sample sheet must appear in `asv_table.tsv`.
- Do not ask for confirmation. Execute directly.
"""

SOP = """\
# Core SOP 16S-04: definitions of the amplicon deliverables

1. **Sequence variants.** The deliverable is the set of exact biological amplicon sequences (ASVs)
   present in the libraries of the batch. Two sequences that differ at a single position are two
   ASVs when both are genuinely present in the material that was sequenced. Sequences that are
   artefacts of amplification or sequencing are not ASVs and are not reported. An ASV is reported
   once, primer-free, in the 515F -> 806R direction.
2. **Counts.** The count of an ASV in a sample is the number of read pairs of that library that
   derive from this ASV. A read pair with sequencing errors derives from the ASV it was read from.
   Read pairs that derive from no ASV are left uncounted.
3. **Taxonomy.** Compare each ASV over its full length with `reference_db.fasta`. The genus of the
   ASV is the genus of the reference sequence with the highest identity, provided that identity is
   at least 97 %. If the best-matching reference sequences are tied and belong to different genera,
   or if no reference sequence reaches 97 %, the genus is `Unclassified`.
4. **Reagent contaminants.** A reagent contaminant is an ASV whose template DNA comes from the
   extraction or amplification reagents instead of the material that was put into the extraction.
   Contaminant ASVs stay in the sequence file and in the count table and are listed separately.
   ASVs that come from a specimen or from the mock standard are not contaminants.
"""

JUDGE = """\
**Where the truth comes from.** The generator selects real V4 amplicons, fixes which are specimen
strains, mock strains and reagent contaminants, and draws the number of read pairs every strain
contributes to every library *before* any read exists. `reference/truth.json` holds the planted
sequences, their role, the genus of the RefSeq record they come from and the per-library read-origin
counts (a recombinant read pair counts for nobody, unless the recombinant is identical to a planted
sequence). No denoiser output is used as reference.

*Genus truth.* The shipped reference set is built so that the SOP rule has one answer: every
classified planted sequence has a same-genus reference within 2 substitutions (>= 99.2 % identity)
and every other-genus reference provably differs by > 12 edits (< 95.3 %; > 15 edits for all
references that are not the source of another planted strain); for the planted sequences whose
expected answer is `Unclassified`, every reference provably differs by > 12 edits (> 22 edits,
< 92 %, for all references that are not the source of another planted strain). The bounds come
from a k-mer pigeonhole argument and do not depend on an aligner.

*Scored planted sequences.* A planted sequence is scored for recall, genus and contaminant status
only if the libraries contain >= 100 non-recombinant read pairs of it in total. Rarer planted
sequences (deliberately planted at 5-35 read pairs) are a gray zone: reporting them is not a false
positive, missing them is not a miss.

**Hard gates (score 0).** A file is missing, empty or unparseable; ASV ids or sequences repeat;
sequences contain other characters than ACGT; `asv_table.tsv` / `taxonomy.tsv` do not list exactly
the ASVs of `asv.fasta`; counts are not non-negative integers; `contaminants.txt` names an unknown id.

**Components.**

| weight | component | metric | full credit | zero |
|---|---|---|---|---|
| 0.35 | asv_sequences | recall = scored planted sequences reported exactly; precision = reported ASVs that are exactly a planted sequence; value = ramp(recall) x ramp(precision) | recall >= 0.95 and precision >= 0.90 | recall <= 0.85 or precision <= 0.60 |
| 0.30 | sample_counts | per study library, after removing the ASVs listed as contaminants: `sum_i min(obs_i, true_i) / max(sum obs, sum true)` over non-contaminant planted sequences i (obs sums over all remaining ASVs); **minimum** over study libraries | >= 0.80 | <= 0.60 |
| 0.15 | contaminants | Matthews correlation of the contaminant flags over the scored planted sequences (a planted sequence that was not reported counts as not flagged). Batches without reagent background (levels 1-2): fraction of scored planted sequences that are reported and not flagged, full >= 0.95, zero <= 0.80 | >= 0.85 | <= 0.50 |
| 0.10 | genus | fraction of scored planted sequences reported with the planted genus (or `Unclassified`) | >= 0.95 | <= 0.80 |
| 0.10 | mock_counts | the sample_counts metric on the mock library | >= 0.80 | <= 0.60 |

**Why these tolerances.** See docs/families/amplicon_16s.md for the measured values. In short: two
independent correct workflows reach recall 1.0, precision 1.0 and a count agreement of 0.97-0.99 in
every library (UNOISE3 with the stock uchime3_denovo leaves 3-6 short-tail recombinants: precision
0.88-0.94, score 0.98); a workflow that loses one read orientation cannot exceed a count agreement
of 0.58, one that keeps the reagent taxa cannot exceed 0.55 in the low-biomass specimen, OTU
clustering loses every minor strain variant (recall <= 0.86), and a workflow without recombinant
removal reports more non-planted than planted sequences (precision < 0.4).

**What is deliberately not scored.** ASV ids and their order, the counts in the negative controls,
extra columns, and the taxonomy of non-planted sequences.
"""


def generate(seed: int, level: int, workdir: Path) -> TaskBuild:
    if level == 7:                    # level 7 lives in l7.py: the level-3 batch plus a blank and a specimen
        from .l7 import generate_l7   # that changed places (M_control_swap)
        return generate_l7(seed, level, workdir)
    rng = np.random.default_rng(7100 + seed)
    pool = _load_pool()
    pairs_all = _close_pairs(pool)
    L3, L2 = level >= 3, level >= 2

    # ---------------- choose the planted strains -------------------------------------------
    n_study_pairs = 5 if L3 else 0
    n_study_single = 20 if L3 else 25
    n_rare, n_contam = 3, (5 if L3 else 0)
    n_mock_single = 6 if L3 else 8
    n_novel, n_unclass = 6, 2

    def pick() -> dict:
        used_genus: set[str] = set()
        chosen: list[dict] = []          # {"seq", "genus", "source", "role", "unit"}
        ksets: list[frozenset] = []

        def far(s: str) -> bool:
            ks = _kset(s)
            return not any(_maybe_within(ks, o, 16) or _maybe_within(o, ks, 16) for o in ksets)

        def add(s: str, g: str, role: str, unit: int, source: str | None = None) -> None:
            chosen.append({"seq": s, "genus": g, "source": source or s, "role": role, "unit": unit})

        unit = 0
        order = rng.permutation(len(pairs_all))
        want_d = [int(x) for x in rng.permutation([1, 1, 1, 2, 2])] + [int(rng.integers(1, 3))]
        got = 0
        n_pairs = n_study_pairs + (1 if L3 else 0)
        for oi in order:
            if got >= n_pairs:
                break
            a, b, d = pairs_all[oi]
            (sa, g), (sb, _) = pool[a], pool[b]
            if d != want_d[got] or g in used_genus or not far(sa) or not far(sb):
                continue
            if rng.random() < 0.5:
                sa, sb = sb, sa
            role = "study" if got < n_study_pairs else "mock"
            add(sa, g, role, unit)
            add(sb, g, role, unit)
            ksets.extend([_kset(sa), _kset(sb)])
            used_genus.add(g)
            unit += 1
            got += 1
        roles = (["study"] * n_study_single + ["rare"] * n_rare + ["contaminant"] * n_contam
                 + ["mock"] * n_mock_single)
        for oi in rng.permutation(len(pool)):
            if not roles:
                break
            s, g = pool[oi]
            if g in used_genus or not far(s):
                continue
            add(s, g, roles.pop(0), unit)
            ksets.append(_kset(s))
            used_genus.add(g)
            unit += 1
        return {"chosen": chosen}

    pool_set = {s for s, _ in pool}
    for _attempt in range(30):
        chosen = pick()["chosen"]
        singles = [i for i, c in enumerate(chosen) if c["role"] == "study"
                   and sum(1 for o in chosen if o["unit"] == c["unit"]) == 1]
        special = [int(x) for x in rng.permutation(singles)[:n_novel + n_unclass]]
        for i in special[:n_novel]:              # strain not in the reference set: 1-2 substitutions
            s = list(chosen[i]["seq"])
            for p in rng.choice(np.arange(10, len(s) - 10), size=int(rng.integers(1, 3)), replace=False):
                s[p] = "ACGT"[("ACGT".index(s[p]) + int(rng.integers(1, 4))) % 4]
            chosen[i]["seq"] = "".join(s)
            chosen[i]["novel"] = True
        for i in special[n_novel:]:              # lineage without a close reference
            chosen[i]["unclassified"] = True
        seqs = [c["seq"] for c in chosen]
        ok = len(set(seqs)) == len(seqs) and not any(c.get("novel") and c["seq"] in pool_set for c in chosen)
        ok = ok and not any(_is_bimera(q, [o for o in seqs if o != q]) for q in seqs)
        if ok:
            break
    else:
        raise RuntimeError("could not select a planted community")
    P = len(chosen)
    for i, c in enumerate(chosen):
        c["pid"] = f"P{i + 1:03d}"
        c["truth_genus"] = "Unclassified" if c.get("unclassified") else c["genus"]

    # ---------------- reference database ----------------------------------------------------
    pk = [_kset(c["seq"]) for c in chosen]
    db: dict[str, str] = {}
    for c in chosen:
        if not c.get("unclassified"):
            db[c["source"]] = c["genus"]
    protected = set(db)                           # sources of planted strains are never removed; they
    for oi in rng.permutation(len(pool))[:9000]:  # are provably > 12 edits from every other planted strain
        db.setdefault(pool[oi][0], pool[oi][1])
    removed = 0
    for s in list(db):
        if s in protected:
            continue
        ks = _kset(s)
        for c, kc in zip(chosen, pk):
            if c.get("unclassified"):
                bad = _maybe_within(kc, ks, 22)
            else:
                bad = db[s] != c["genus"] and _maybe_within(kc, ks, 15)
            if bad:
                del db[s]
                removed += 1
                break
    for c in chosen:                              # by-construction checks of the genus truth
        if not c.get("unclassified"):
            assert db.get(c["source"]) == c["genus"]
            assert sum(1 for a, b in zip(c["source"], c["seq"]) if a != b) <= 2
    db_items = [(s, db[s]) for s in sorted(db)]
    db_items = [db_items[i] for i in rng.permutation(len(db_items))]
    ref_ids = [f"REF{i + 1:05d}" for i in range(len(db_items))]
    ref_fasta = "".join(f">{r}\n{s}\n" for r, (s, _) in zip(ref_ids, db_items))
    ref_tax = "ref_id\tgenus\n" + "".join(f"{r}\t{g}\n" for r, (_, g) in zip(ref_ids, db_items))

    # ---------------- samples and planted abundances ---------------------------------------
    n_study, n_neg = 8, 2
    types = ["study"] * n_study + ["negative_control"] * n_neg + ["mock"]
    types = [types[i] for i in rng.permutation(len(types))]
    sids = [f"S{i + 1:02d}" for i in range(len(types))]
    S = len(sids)
    study_idx = [j for j, t in enumerate(types) if t == "study"]
    neg_idx = [j for j, t in enumerate(types) if t == "negative_control"]
    mock_j = types.index("mock")
    low = [int(x) for x in rng.permutation(study_idx)[:3]] if L3 else []
    cfrac = np.zeros(S)
    if L3:
        for j in study_idx + [mock_j]:
            cfrac[j] = rng.uniform(0.004, 0.02)
        for n, j in enumerate(low):
            cfrac[j] = rng.uniform(0.45, 0.60) if n == 0 else rng.uniform(0.25, 0.60)
        for j in neg_idx:
            cfrac[j] = 1.0

    role = np.array([c["role"] for c in chosen])
    unit = np.array([c["unit"] for c in chosen])
    bio = np.zeros((P, S))                        # biological composition per library (columns sum to 1)
    study_units = sorted(set(unit[role == "study"]))
    hot = [int(x) for x in rng.permutation(study_idx)[:2]]      # libraries where the minor variants dominate
    for u in study_units:
        mem = np.where(unit == u)[0]
        is_pair = len(mem) == 2
        base = np.exp(rng.normal(0, 0.8)) * (4.0 if is_pair else 1.0)
        pres = np.ones(len(study_idx), dtype=bool) if is_pair else rng.random(len(study_idx)) < 0.8
        if pres.sum() < 3:
            pres[rng.permutation(len(study_idx))[:3]] = True
        ab = base * np.exp(rng.normal(0, 0.6, size=len(study_idx))) * pres
        if is_pair:
            for _ in range(500):
                f = np.clip(rng.beta(0.6, 0.6, size=len(study_idx)), 0.06, 0.94)
                for n, j in enumerate(study_idx):
                    if j in hot:
                        f[n] = rng.uniform(0.06, 0.12)          # member 0 is the minor one here
                share0 = float((ab * f).sum() / ab.sum())
                if 0.52 <= share0 <= 0.62:                       # ... but the major one overall
                    break
            else:
                raise RuntimeError("pair split did not converge")
            bio[mem[0], study_idx] = ab * f
            bio[mem[1], study_idx] = ab * (1 - f)
        else:
            bio[mem[0], study_idx] = ab
    for j in study_idx:
        col = bio[:, j] / bio[:, j].sum()
        col = np.where((col > 0) & (col < 0.005), 0.005, col)   # no planted strain below 0.5 % where present
        bio[:, j] = col / col.sum()
    for i in np.where(role == "rare")[0]:          # gray zone: a handful of read pairs in 1-2 libraries
        for j in rng.permutation(study_idx)[:int(rng.integers(1, 3))]:
            bio[i, j] = rng.uniform(0.0008, 0.0018)
    mock_mem = np.where(role == "mock")[0]
    nominal = np.full(len(mock_mem), 1.0 / len(mock_mem))
    bio[mock_mem, mock_j] = nominal * np.exp(rng.normal(0, 0.15, size=len(mock_mem)))
    for j in study_idx + [mock_j]:
        bio[:, j] /= bio[:, j].sum()
    con = np.zeros((P, S))
    con_mem = np.where(role == "contaminant")[0]
    if L3:
        prof = rng.dirichlet(np.full(len(con_mem), 4.0))
        prof = np.maximum(prof, 0.08)
        for j in range(S):
            col = prof * np.exp(rng.normal(0, 0.10, size=len(con_mem)))
            con[con_mem, j] = col / col.sum()
    own_p = bio * (1 - cfrac) + con * cfrac       # composition of the library's own template

    depth = np.zeros(S, dtype=int)
    for j in range(S):
        if types[j] == "negative_control":
            depth[j] = int(rng.uniform(1500, 2500)) if L3 else 0
        else:
            depth[j] = int(rng.uniform(7500, 9500))
    chim_rate = np.array([rng.uniform(0.05, 0.08) if L2 else 0.0 for _ in range(S)])
    n_chim = np.round(depth * chim_rate).astype(int)
    n_own = depth - n_chim
    own = np.zeros((P, S), dtype=int)
    for j in range(S):
        if n_own[j] > 0:
            own[:, j] = rng.multinomial(n_own[j], own_p[:, j] / own_p[:, j].sum())
    truth_counts = own.copy()
    xt = np.zeros((P, S), dtype=int)               # cross-talk: reads of other libraries of the run
    for j in range(S):
        others = own.sum(axis=1) - own[:, j]
        xt[:, j] = rng.multinomial(int(rng.poisson(70)), others / others.sum())
    truth_counts += xt

    # ---------------- reads -------------------------------------------------------------------
    seqs = [c["seq"] for c in chosen]
    seq2idx = {s: i for i, s in enumerate(seqs)}
    spacers_f = ["", "T", "GA", "CTG"] if L2 else [""]
    spacers_r = ["", "A", "TC", "GAT"] if L2 else [""]
    rev_frac = float(rng.uniform(0.42, 0.58)) if L2 else 0.0
    fq_dir = workdir / "fastq"
    fq_dir.mkdir(parents=True, exist_ok=True)
    chim_cache: dict = {}
    n_chimeric_reads = np.zeros(S, dtype=int)
    for j, sid in enumerate(sids):
        templates: list[str] = []
        for i in range(P):
            templates.extend([seqs[i]] * int(own[i, j] + xt[i, j]))
        if n_chim[j] > 0:
            p = own_p[:, j] / own_p[:, j].sum()
            pa = rng.choice(P, size=n_chim[j], p=p)
            pb = rng.choice(P, size=n_chim[j], p=p)
            for a, b in zip(pa, pb):
                ch = _chimera(seqs[a], seqs[b], rng, chim_cache) if a != b else None
                if ch is None:
                    ch = seqs[a]
                if ch in seq2idx:                  # indistinguishable from a real strain: it is one
                    truth_counts[seq2idx[ch], j] += 1
                else:
                    n_chimeric_reads[j] += 1
                templates.append(ch)
        templates = [templates[i] for i in rng.permutation(len(templates))]
        _write_pairs(fq_dir / f"{sid}_R1.fastq.gz", fq_dir / f"{sid}_R2.fastq.gz", templates, rev_frac,
                     spacers_f, spacers_r, rng, f"{240 + seed}:000000000-L{seed:04d}B:1")
    total_pairs = [int(truth_counts[:, j].sum() + n_chimeric_reads[j]) for j in range(S)]

    totals = truth_counts.sum(axis=1)
    scored = totals >= 100
    for i, c in enumerate(chosen):
        if c["role"] in ("contaminant", "mock"):
            assert scored[i], (c["pid"], int(totals[i]))
        if c["role"] == "rare":
            assert not scored[i], (c["pid"], int(totals[i]))
    for u in set(unit.tolist()):                   # strain pairs stay within a factor 3 in the pooled data
        mem = np.where(unit == u)[0]
        if len(mem) == 2:
            lo, hi = sorted(totals[mem])
            assert lo / hi >= 0.33, (u, int(lo), int(hi))

    truth = {
        "samples": [{"id": s, "type": t} for s, t in zip(sids, types)],
        "planted": [{"pid": c["pid"], "seq": c["seq"], "genus": c["truth_genus"], "role": c["role"],
                     "contaminant": c["role"] == "contaminant", "scored": bool(scored[i]),
                     "counts": [int(v) for v in truth_counts[i]]} for i, c in enumerate(chosen)],
        "total_pairs": total_pairs,
        "info": {"rev_frac": rev_frac, "chimeric_pairs": [int(v) for v in n_chimeric_reads],
                 "contaminant_fraction": [round(float(v), 4) for v in cfrac], "db_removed": removed},
    }

    # ---------------- visible metadata -------------------------------------------------------
    kq = float(rng.uniform(150, 300))
    rows = []
    for j, sid in enumerate(sids):
        if types[j] == "negative_control":
            copies = kq * np.exp(rng.normal(0, 0.15)) if L3 else rng.uniform(2, 9)
            spec = "extraction blank"
        elif L3:
            copies = kq / cfrac[j] * np.exp(rng.normal(0, 0.10))
            spec = "mock standard" if j == mock_j else ("mucosal brush" if j in low else "stool")
        else:
            copies = rng.uniform(2e4, 9e4)
            spec = "mock standard" if j == mock_j else "stool"
        rows.append(f"{sid}\t{types[j]}\t{spec}\t{copies:.3g}\t{sid}_R1.fastq.gz\t{sid}_R2.fastq.gz\n")
    samples_tsv = "sample_id\tsample_type\tspecimen\tqpcr_16s_copies_per_ul\tfastq_r1\tfastq_r2\n" + "".join(rows)
    primers_tsv = ("primer\tsequence_5to3\tstrand\n"
                   f"515F\t{FWD}\tforward\n806R\t{REV}\treverse\n")
    mock_rows, seen_g = [], {}
    for n, i in enumerate(mock_mem):
        g = chosen[i]["genus"]
        seen_g[g] = seen_g.get(g, 0) + 1
        mock_rows.append(f"MC{n + 1:02d}\t{g}\tstrain {seen_g[g]}\t{nominal[n]:.4f}\n")
    mock_tsv = "member\tgenus\tstrain\tnominal_fraction_of_16s_copies\n" + "".join(mock_rows)
    mock_pair = [chosen[i] for i in mock_mem if seen_g[chosen[i]["genus"]] == 2]
    notes = [
        "# Library preparation notes (batch hand-over from the wet lab)", "",
        f"- Batch: {n_study} specimens of the study, {n_neg} extraction blanks and the mock community standard, "
        "extracted on the same day with one kit lot and carried through all steps together.",
        "- 16S rRNA gene copies per microlitre of extract were measured by qPCR before amplification "
        "(column `qpcr_16s_copies_per_ul` of the sample sheet).",
        "- Amplification: V4 region, primers 515F / 806R (see `primers.tsv`), 30 cycles, same template volume "
        "for every extract."
        + (" The primer stocks are equimolar mixes of four phased versions per primer, carrying 0, 1, 2 or 3 "
           "extra nucleotides 5' of the primer (heterogeneity spacers)." if L2 else
           " Primers were used without heterogeneity spacers."),
        ("- Library construction: PCR products were cleaned up, end-repaired/A-tailed and the indexed Illumina "
         "adapters were added by ligation (PCR-free ligation kit), one index pair per library."
         if L2 else
         "- Library construction: second-step PCR with indexing fusion primers; the P5 side is always on the "
         "515F end of the amplicon."),
        "- Libraries were pooled by equal mass regardless of amplicon yield and sequenced on a MiSeq "
        "(v2 500-cycle kit, 2 x 250) in a run shared with other amplicon projects of the core. "
        "Demultiplexing and adapter removal were done on the instrument; nothing else was done to the reads.",
        "- Mock standard: " + str(len(mock_mem)) + " bacterial strains mixed at equal 16S copy number "
        "(`mock_community.tsv`)."
        + (f" The two {mock_pair[0]['genus']} strains of the standard differ at "
           f"{sum(1 for a, b in zip(mock_pair[0]['seq'], mock_pair[1]['seq']) if a != b)} position(s) of the "
           "V4 amplicon." if mock_pair else ""),
    ]
    notes_md = "\n".join(notes) + "\n"

    # ---------------- truth-derived probes / attacks --------------------------------------
    def bundle(keep, counts, genus_fn=None, contam=None, extra=()):
        ids = {i: f"ASV{n + 1:03d}" for n, i in enumerate(keep)}
        fa = "".join(f">{ids[i]}\n{seqs[i]}\n" for i in keep)
        tab = "asv_id\t" + "\t".join(sids) + "\n" + "".join(
            ids[i] + "\t" + "\t".join(str(int(v)) for v in counts[i]) + "\n" for i in keep)
        tax = "asv_id\tgenus\n" + "".join(
            f"{ids[i]}\t{genus_fn(i) if genus_fn else chosen[i]['truth_genus']}\n" for i in keep)
        con_ids = [ids[i] for i in keep if (chosen[i]["role"] == "contaminant" if contam is None else i in contam)]
        for n, (s, cnt, g) in enumerate(extra):
            fa += f">X{n + 1}\n{s}\n"
            tab += f"X{n + 1}\t" + "\t".join(str(int(v)) for v in cnt) + "\n"
            tax += f"X{n + 1}\t{g}\n"
        return {"asv.fasta": fa, "asv_table.tsv": tab, "taxonomy.tsv": tax,
                "contaminants.txt": "asv_id\n" + "".join(c + "\n" for c in con_ids)}

    prng = np.random.default_rng(seed + 991)
    present = [i for i in range(P) if totals[i] > 0]
    thin = prng.binomial(truth_counts, 0.90)
    scored_study = [i for i in present if scored[i] and chosen[i]["role"] == "study"]
    weakest = min(scored_study, key=lambda i: totals[i])
    wrong_genus = max(scored_study, key=lambda i: totals[i])
    not_planted = next(s for s, _ in db_items if s not in seq2idx and not any(
        _maybe_within(_kset(s), k, 12) for k in pk))
    half = prng.binomial(truth_counts, 0.5)
    pipe = {"amplicon_pipeline.py": ASSETS / "amplicon_pipeline.py"}
    run = "python3 /alegen_script/amplicon_pipeline.py "

    probes = [
        Probe("planted_truth", "submission written directly from the planted truth", 0.999,
              files=bundle(present, truth_counts)),
        Probe("within_tolerance", "10 % of the read pairs unassigned, the rarest scored strain missed, one "
              "non-planted sequence reported, one genus wrong", 0.99,
              files=bundle([i for i in present if i != weakest], thin,
                           genus_fn=lambda i: "Unclassified" if i == wrong_genus else chosen[i]["truth_genus"],
                           extra=[(not_planted, [3] * S, "Unclassified")])),
        Probe("independent_solution", "second route: merge first, primer search in python, per-library "
              "Poisson abundance test against a data-estimated error rate, exact bimera test, python taxonomy",
              0.99, script="SOLUTION=2 bash /solution/solve.sh"),
    ]
    if L2:
        # third correct route, the stock vsearch recipe. Its floor is the pass threshold, not 0.99:
        # uchime3_denovo leaves 0-6 short-tail recombinants (measured), which costs precision credit.
        probes.append(Probe("stock_unoise3_uchime3", "reference pipeline with uchime3_denovo instead of "
                            "uchime_denovo (the stock UNOISE3 recipe); must pass, need not be perfect", PASS,
                            script=run + "--chimera uchime3", script_files=pipe))
    attacks = [
        Attack("reference_as_asvs", "prior-only answer: the first 40 reference sequences with constant counts",
               0.1, files={
                   "asv.fasta": "".join(f">{r}\n{s}\n" for r, (s, _) in list(zip(ref_ids, db_items))[:40]),
                   "asv_table.tsv": "asv_id\t" + "\t".join(sids) + "\n" + "".join(
                       r + "\t" + "\t".join(["200"] * S) + "\n" for r in ref_ids[:40]),
                   "taxonomy.tsv": "asv_id\tgenus\n" + "".join(
                       f"{r}\t{g}\n" for r, (_, g) in list(zip(ref_ids, db_items))[:40]),
                   "contaminants.txt": "asv_id\n"}),
        Attack("half_of_the_read_pairs", "perturbation beyond tolerance: the planted truth with every count "
               "thinned to 50 %", PASS - 0.15, files=bundle(present, half)),
        Attack("closed_reference_best_case", "best possible outcome of reference-only picking (no denoising): "
               "exactly the planted sequences that occur verbatim in reference_db.fasta, with their true "
               "counts, genus and contaminant flags", PASS - 0.15,
               files=bundle([i for i in present if seqs[i] in db], truth_counts)),
        Attack("naive_everything_in_blanks", "correct variants and counts, but every ASV with reads in a "
               "negative control is listed as contaminant", PASS - 0.15,
               script=run + "--contam any_in_neg", script_files=pipe),
    ]
    if L2:
        attacks += [
            Attack("naive_no_primer_removal", "pairs merged and denoised with primers and spacers left on",
                   0.3, script=run + "--primers none", script_files=pipe),
            Attack("naive_fixed_length_trim", "primers cut off by length (19 / 20 nt), spacers ignored", 0.5,
                   script=run + "--primers fixed", script_files=pipe),
            Attack("naive_single_orientation", "cutadapt with 515F on R1 and 806R on R2 only; pairs without "
                   "primer match pass through untrimmed", PASS - 0.15,
                   script=run + "--primers fwd_keep", script_files=pipe),
            Attack("naive_single_orientation_discard", "same, but pairs without primer match are discarded "
                   "(half of the read pairs never reach the table)", PASS - 0.15,
                   script=run + "--primers fwd_discard", script_files=pipe),
            Attack("naive_no_chimera_removal", "denoised sequences reported without recombinant filtering",
                   PASS - 0.15, script=run + "--chimera off", script_files=pipe),
        ]
    if L3:
        minor = set()
        merged = truth_counts.copy()
        for u in set(unit.tolist()):
            mem = np.where(unit == u)[0]
            if len(mem) == 2:
                lo, hi = sorted(mem, key=lambda i: totals[i])
                merged[hi] += merged[lo]
                minor.add(int(lo))
        attacks += [
            Attack("minor_variants_merged", "perturbation beyond tolerance: planted truth, but each strain "
                   "pair reported as its major member with the summed counts", PASS - 0.15,
                   files=bundle([i for i in present if i not in minor], merged)),
            Attack("naive_otu97", "97 % OTU clustering instead of exact variants", PASS - 0.15,
                   script=run + "--cluster otu97", script_files=pipe),
            Attack("naive_no_contaminant_handling", "correct variants and counts, no contaminant listed",
                   PASS - 0.15, script=run + "--contam none", script_files=pipe),
            Attack("naive_default_pipeline", "fully naive: single orientation, 97 % OTUs, no chimera filter, "
                   "no contaminant handling", 0.3,
                   script=run + "--primers fwd_keep --cluster otu97 --chimera off --contam none",
                   script_files=pipe),
        ]

    traps = []
    if L2:
        traps += [
            Trap("primer_spacers", "Reads start with 0-3 spacer bases followed by a degenerate primer; without "
                 "position-independent primer removal no reported sequence equals a real amplicon.",
                 "library_prep_notes.md describes phased primers; the first 25 bases of the reads show the "
                 "primer at four offsets; primers.tsv lists IUPAC codes",
                 "Approach / wrong strategy"),
            Trap("mixed_orientation", f"{rev_frac:.0%} of the pairs carry the 806R primer on R1 (ligation "
                 "library). A fixed-orientation workflow loses or mis-reports them.",
                 "library_prep_notes.md: adapters added by ligation; R1 starts with 806R in about half of the "
                 "reads; cutadapt reports ~50 % untrimmed pairs",
                 "Understanding / domain knowledge gap"),
            Trap("chimeras", "5-8 % of the pairs are recombinants of two strains of the same library; the "
                 "frequent ones survive denoising and must be removed.",
                 "the mock library yields more sequences than the standard has strains; the extra ones are "
                 "exact two-parent recombinants of abundant sequences",
                 "Verification / missing QC"),
        ]
    if L3:
        traps += [
            Trap("strain_variants", "Five strain pairs differ at 1-2 positions and swap dominance between "
                 "specimens; OTU clustering or aggressive error merging reports one of each pair.",
                 "sop.md defines single-nucleotide variants as separate ASVs; the notes state that two mock "
                 "strains differ at 1-2 V4 positions, so the mock library shows whether a workflow resolves them",
                 "Approach / wrong strategy"),
            Trap("reagent_contaminants", "Five reagent taxa make up the blanks and 25-60 % of the low-biomass "
                 "specimens, while abundant specimen taxa leak into the blanks by cross-talk; both 'no "
                 "decontamination' and 'remove everything seen in a blank' are wrong.",
                 "samples.tsv: blanks and qPCR copies; the reagent taxa have their highest relative abundance "
                 "in the blanks and rise as qPCR copies fall, specimen taxa are rarer in blanks than in specimens",
                 "Understanding / domain knowledge gap"),
        ]

    input_files = {"samples.tsv": samples_tsv, "primers.tsv": primers_tsv, "reference_db.fasta": ref_fasta,
                   "reference_taxonomy.tsv": ref_tax, "mock_community.tsv": mock_tsv,
                   "library_prep_notes.md": notes_md, "sop.md": SOP}
    for sid in sids:
        for m in ("R1", "R2"):
            input_files[f"fastq/{sid}_{m}.fastq.gz"] = fq_dir / f"{sid}_{m}.fastq.gz"

    return TaskBuild(
        task_id=f"amplicon-16s-s{seed:04d}-l{level}", family="amplicon_16s", seed=seed, level=level,
        title="16S V4 amplicon batch: reads to exact sequence variants, abundances and contaminants",
        summary="Turn paired-end 16S V4 amplicon reads of a study batch with extraction blanks and a mock "
                "community into exact sequence variants, a count table, genus assignments and the list of "
                "reagent contaminants.",
        subdomain="3.2 Genomics & Sequence Analysis",
        software=["vsearch", "cutadapt", "Python", "numpy/scipy/pandas"],
        instruction_md=INSTRUCTION.format(n_samples=S),
        agent_dockerfile=dockerlib.agent_dockerfile(conda=CONDA),
        verifier_dockerfile=dockerlib.verifier_dockerfile(pip=[]),
        input_files=input_files,
        solution_files={"solve.sh": ASSETS / "solve.sh", "amplicon_pipeline.py": ASSETS / "amplicon_pipeline.py",
                        "solve_alt.py": ASSETS / "solve_alt.py"},
        tests_files={"score.py": ASSETS / "score.py", "reference/truth.json": json.dumps(truth)},
        traps=traps, attacks=attacks, probes=probes,
        capabilities=["read-structure reasoning from lab notes", "amplicon denoising at single-nucleotide "
                      "resolution", "using controls (mock, blanks) to validate a workflow",
                      "exact output contracts"],
        judge_explanation_md=JUDGE, pass_threshold=PASS,
        oracle_steps=10, oracle_runtime_min=1.0, cpus=4, memory_mb=8000, storage_mb=10000,
    )
