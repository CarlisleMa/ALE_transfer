#!/usr/bin/env python3
"""Independent second solution, level 7. It shares only the read merger with the reference pipeline and
differs where the answer is decided:

  order of operations   pairs are merged first (primers still on); primers are then located on the
                        merged read in either strand by an IUPAC search written here (no cutadapt)
  denoising             per library, not pooled: a sequence is accepted when its abundance cannot be
                        explained (Poisson test) as substitution errors of a more abundant accepted
                        sequence, with the error rate estimated from the library itself (no UNOISE)
  chimeras              exact two-parent (prefix + suffix) test against more abundant sequences of the
                        same library (no uchime)
  counts                every merged read goes to its nearest variant by Hamming distance (numpy)
  taxonomy              k-mer screen + banded edit distance written here (no usearch_global)
  sample identity       mass balance against the qPCR column, without the sample types and without any
                        comparison of profiles: ASVs whose absolute abundance (relative abundance x qPCR
                        copies) is the same in most libraries are found by a robust spread; a library whose
                        reads do not belong to its row of the sample sheet is off that common amount by the
                        ratio of the two qPCR values. Two libraries that are off in opposite directions and
                        fall into place when their rows are exchanged have changed places.
                        (The reference pipeline compares profiles: Bray-Curtis and richness, no qPCR.)
  contaminants          constancy of the absolute abundance (relative abundance x qPCR copies) between
                        blanks and specimens, instead of enrichment in the blanks
"""
import json
import re
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

inp, out, work = Path("/task/input"), Path("/task/output"), Path("/tmp/amplicon_alt")
out.mkdir(parents=True, exist_ok=True)
work.mkdir(parents=True, exist_ok=True)
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "M": "AC", "K": "GT", "S": "CG",
         "W": "AT", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
COMP = str.maketrans("ACGTN", "TGCAN")

samples = pd.read_csv(inp / "samples.tsv", sep="\t", dtype=str)
sids = list(samples["sample_id"])
prim = pd.read_csv(inp / "primers.tsv", sep="\t", dtype=str)
FWD = prim.loc[prim["strand"] == "forward", "sequence_5to3"].iloc[0]
REV = prim.loc[prim["strand"] == "reverse", "sequence_5to3"].iloc[0]


def rc(s):
    return s.translate(COMP)[::-1]


class Primer:
    def __init__(self, seq):
        self.sets = [IUPAC[c] for c in seq]
        self.n = len(seq)
        self.rx = re.compile(".{0,6}?" + "".join(c if len(IUPAC[c]) == 1 else f"[{IUPAC[c]}]" for c in seq))

    def end(self, read):
        """Position right after the primer near the 5' end of read (<= 2 mismatches), or -1."""
        m = self.rx.match(read)
        if m:
            return m.end()
        for off in range(0, 7):
            win = read[off:off + self.n]
            if len(win) == self.n and sum(1 for c, st in zip(win, self.sets) if c not in st) <= 2:
                return off + self.n
        return -1


PF, PR = Primer(FWD), Primer(REV)
QERR = np.array([10 ** (-(q - 33) / 10) if q >= 33 else 1.0 for q in range(256)])


def amplicon(seq, qual):
    """Primer-free amplicon in the 515F->806R direction plus its expected errors, or None."""
    for s, q in ((seq, qual), (rc(seq), qual[::-1])):
        a = PF.end(s)
        if a < 0:
            continue
        b = PR.end(rc(s))
        if b < 0:
            continue
        core = s[a:len(s) - b]
        if len(core) < 100:
            return None
        ee = float(QERR[np.frombuffer(q[a:len(s) - b].encode(), dtype=np.uint8)].sum())
        return core, ee
    return None


# ---------------------------------------------------------------- 1. merge, then orient and trim
all_reads, good_reads = {}, {}
for row in samples.itertuples():
    s = row.sample_id
    m = work / f"{s}.merged.fq"
    r = subprocess.run(["vsearch", "--fastq_mergepairs", str(inp / "fastq" / row.fastq_r1), "--reverse",
                        str(inp / "fastq" / row.fastq_r2), "--fastqout", str(m), "--threads", "2", "--quiet"],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(r.stderr[-2000:])
    every, good = Counter(), Counter()
    with open(m) as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            seq = fh.readline().strip()
            fh.readline()
            qual = fh.readline().strip()
            res = amplicon(seq, qual)
            if res is None or "N" in res[0]:
                continue
            every[res[0]] += 1
            if res[1] <= 1.0:
                good[res[0]] += 1
    all_reads[s], good_reads[s] = every, good


def enc(s):
    return np.frombuffer(s.encode(), dtype=np.uint8)


# ---------------------------------------------------------------- 2. per-library denoising + bimeras
ALPHA, MIN_ABUND, SAFETY = 1e-8, 4, 10.0
votes = defaultdict(lambda: [0, 0])          # sequence -> [libraries accepted, libraries bimeric]
for s in sids:
    uniq = sorted(good_reads[s].items(), key=lambda kv: (-kv[1], kv[0]))
    if not uniq:
        continue
    top, n_top = uniq[0]
    e_top = enc(top)
    one_off = 0
    for seq, c in uniq[1:]:
        if len(seq) == len(top) and c <= 0.05 * n_top and int((enc(seq) != e_top).sum()) == 1:
            one_off += c
    # rate of one specific substitution at one specific position, with a safety factor for hot spots
    rate = SAFETY * min(max(one_off / n_top, 0.02), 0.5) / (3 * len(top))
    centers = []                             # (sequence, count, encoded)
    for seq, c in uniq:
        if c < MIN_ABUND:
            break
        e = enc(seq)
        explained = False
        for cs, cc, ce in centers:
            if len(cs) != len(seq):
                continue
            d = int((ce != e).sum())
            if d <= 3 and stats.poisson.sf(c - 1, cc * rate ** d) >= ALPHA:
                explained = True
                break
        if not explained:
            centers.append((seq, c, e))
    for n, (seq, c, _) in enumerate(centers):
        best_p, arg_p, best_s, arg_s = 0, -1, 0, -1
        for k, (ps, pc, _) in enumerate(centers[:n]):
            if pc < 2 * c:
                continue
            lim = min(len(seq), len(ps))
            p = 0
            while p < lim and seq[p] == ps[p]:
                p += 1
            q = 0
            while q < lim and seq[-1 - q] == ps[-1 - q]:
                q += 1
            if p > best_p:
                best_p, arg_p = p, k
            if q > best_s:
                best_s, arg_s = q, k
        bimera = best_p + best_s >= len(seq) and arg_p != arg_s
        votes[seq][1 if bimera else 0] += 1

asv = [s for s, (acc, bim) in votes.items() if acc > 0]
pooled = Counter()
for s in sids:
    for seq in asv:
        pooled[seq] += good_reads[s].get(seq, 0)
asv.sort(key=lambda q: (-pooled[q], q))
ids = [f"SV{n + 1:03d}" for n in range(len(asv))]
(out / "asv.fasta").write_text("".join(f">{i}\n{q}\n" for i, q in zip(ids, asv)))

# ---------------------------------------------------------------- 3. counts: nearest variant
index = {q: n for n, q in enumerate(asv)}
by_len = defaultdict(list)
for n, q in enumerate(asv):
    by_len[len(q)].append(n)
mats = {L: np.stack([enc(asv[n]) for n in v]) for L, v in by_len.items()}
table = np.zeros((len(asv), len(sids)), dtype=int)
for j, s in enumerate(sids):
    for seq, c in all_reads[s].items():
        n = index.get(seq)
        if n is None and len(seq) in mats:
            d = (mats[len(seq)] != enc(seq)).sum(axis=1)
            k = int(d.argmin())
            if d[k] <= 0.03 * len(seq):
                n = by_len[len(seq)][k]
        if n is not None:
            table[n, j] += c
tab = pd.DataFrame(table, columns=sids)
tab.insert(0, "asv_id", ids)
tab.to_csv(out / "asv_table.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- 4. genus: best reference >= 97 %
ref = {}
name = None
for line in open(inp / "reference_db.fasta"):
    if line.startswith(">"):
        name = line[1:].split()[0]
        ref[name] = ""
    else:
        ref[name] += line.strip().upper()
ref_genus = dict(pd.read_csv(inp / "reference_taxonomy.tsv", sep="\t", dtype=str).values)
ref_names = list(ref)
exact_ref = {v: k for k, v in ref.items()}
KM = 8
inv = defaultdict(list)
for n, r in enumerate(ref_names):
    sq = ref[r]
    for km in {sq[i:i + KM] for i in range(len(sq) - KM + 1)}:
        inv[km].append(n)


def edit_distance(x, y, band=14):
    """Banded Levenshtein distance; returns band + 1 when the distance exceeds the band."""
    if abs(len(x) - len(y)) > band:
        return band + 1
    big = band + 1
    prev = {j: j for j in range(0, band + 1)}
    for i in range(1, len(x) + 1):
        cur = {}
        lo, hi = max(0, i - band), min(len(y), i + band)
        xi = x[i - 1]
        for j in range(lo, hi + 1):
            if j == 0:
                cur[j] = i
                continue
            v = prev.get(j - 1, big) + (xi != y[j - 1])
            w = prev.get(j, big) + 1
            if w < v:
                v = w
            w = cur.get(j - 1, big) + 1
            if w < v:
                v = w
            cur[j] = v
        if min(cur.values()) > band:
            return big
        prev = cur
    return min(prev.get(len(y), big), big)


genus = []
for q in asv:
    if q in exact_ref:
        genus.append(ref_genus[exact_ref[q]])
        continue
    tally = Counter()
    for km in {q[i:i + KM] for i in range(len(q) - KM + 1)}:
        tally.update(inv.get(km, ()))
    best, best_g = 0.0, set()
    for n, _ in tally.most_common(25):
        t = ref[ref_names[n]]
        ident = 1.0 - edit_distance(q, t) / max(len(q), len(t))
        if ident > best + 1e-9:
            best, best_g = ident, {ref_genus[ref_names[n]]}
        elif abs(ident - best) <= 1e-9:
            best_g.add(ref_genus[ref_names[n]])
    genus.append(next(iter(best_g)) if best >= 0.97 and len(best_g) == 1 else "Unclassified")
pd.DataFrame({"asv_id": ids, "genus": genus}).to_csv(out / "taxonomy.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- 5. do the reads belong to their row of the sheet?
# DNA that comes with the reagents is the same absolute amount in every extraction: its relative
# abundance times the 16S copies of the extract is constant, blanks included. ASVs of that kind are
# recognised by exactly that property (robust spread over the libraries of log10 absolute abundance),
# without the sample types. For a library whose FASTQ files hold another extract, the qPCR value in its
# row is not the one of its reads, and all those ASVs are off by the ratio of the two values.
copies = np.array(samples["qpcr_16s_copies_per_ul"].astype(float))
types = list(samples["sample_type"])
S = len(sids)
tot = np.maximum(table.sum(axis=0), 1)
SPREAD_MAX, OFF_MIN, BACK_MAX = 0.30, 1.0, 0.5            # decades
report = {"constant_amount_asvs": 0}
if len(asv) and np.all(copies > 0):
    logab = np.log10((table + 0.5) / tot * copies)
    centre = np.median(logab, axis=1)
    spread = 1.4826 * np.median(np.abs(logab - centre[:, None]), axis=1)
    const = np.where(((table > 0).sum(axis=1) >= int(np.ceil(0.8 * S))) & (spread <= SPREAD_MAX))[0]
    report["constant_amount_asvs"] = int(len(const))
    if len(const) >= 2:
        off = np.median(logab[const] - centre[const, None], axis=0)
        report["decades_off_the_common_amount"] = {s: round(float(v), 2) for s, v in zip(sids, off)}
        odd = [j for j in range(S) if abs(off[j]) >= OFF_MIN]
        if len(odd) == 2 and off[odd[0]] * off[odd[1]] < 0:
            j1, j2 = odd
            shift = float(np.log10(copies[j2] / copies[j1]))
            back = (float(off[j1] + shift), float(off[j2] - shift))
            if max(abs(back[0]), abs(back[1])) < BACK_MAX:
                report["exchanged"] = {"libraries": [sids[j1], sids[j2]],
                                       "decades_off_after_exchange": [round(back[0], 2), round(back[1], 2)]}
                copies[[j1, j2]] = copies[[j2, j1]]
                types[j1], types[j2] = types[j2], types[j1]
print("sheet check: " + json.dumps(report))

# ---------------------------------------------------------------- 6. reagent contaminants
# With every library in the row that belongs to it: reagent DNA has the same absolute abundance in
# blanks and specimens; DNA from a specimen is, in absolute terms, (nearly) absent from the blanks.
neg = [j for j, t in enumerate(types) if t == "negative_control"]
real = [j for j in range(len(sids)) if j not in neg]
absab = table / tot * copies
flag = np.zeros(len(asv), dtype=bool)
if neg:
    in_all = (table[:, neg] >= 2).all(axis=1)
    flag = in_all & (absab[:, neg].mean(axis=1) >= 0.3 * absab[:, real].mean(axis=1))
(out / "contaminants.txt").write_text("asv_id\n" + "".join(f"{i}\n" for i, f in zip(ids, flag) if f))
print(f"variants={len(asv)} assigned={table.sum(axis=0).tolist()} contaminants={int(flag.sum())}")
