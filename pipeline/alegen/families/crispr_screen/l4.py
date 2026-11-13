"""Level 4 of the pooled CRISPR screen family: taking over a drug-modifier screen (H + P + Q).

Level 3 was solved by a small model in two minutes because the lab notes named every decision
point. Level 4 keeps the three level-3 decision points (stagger, Read-2 library, hit-enriched
library) and adds what the pilot showed to be hard:

  H  inherited project   a predecessor's scripts, count table, gene table and mailed hit list; a
                         seeded subset of an 8-item checklist is defective. Deliverables: the
                         corrected results, `audit.tsv`, `erratum.tsv`.
  P  prior != truth      (a) the sample sheet has one vehicle and one drug library exchanged; the
                         sgRNA profiles show it. (b) "depleted in the drug arm compared with T0"
                         (the standard dropout contrast) does not identify sensitisers, because
                         fitness and drug effects are planted independently in a hit-enriched library.
  Q  quality bar         drug-arm libraries come from few cells, so noise grows steeply as sgRNAs
                         drop out. Calls made against the spread of the non-targeting sgRNAs are
                         plausible and wrong; reference genes of known behaviour let the analyst
                         measure this, and the SOP states the release criterion on them.

Everything is planted before a read exists (`plan`). The inherited files are produced by running
the predecessor's scripts (seeded mutations of a correct stdlib/numpy pipeline) on the task data.
"""
from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

from ... import dockerlib
from ...spec import Attack, Probe, TaskBuild, Trap
from . import ASSETS, PREFIX, SCAFFOLD, _guides, _rc, _write_fastq

N_FIT, N_NEU, N_NTC, PER = 190, 40, 80, 4
N_SENS, N_RES = 18, 18
N_REF_MOD_EACH, N_REF_INERT = 5, 30
N_V1_PAIRS = 10                      # release-1 annotation: 10 modifiers + 10 non-modifiers mislabelled
COV = {"T0": 500.0, "vehicle": 150.0, "drug": 25.0}   # cells per sgRNA that reach the PCR
JITTER, BIO = 0.06, 0.05
DEPTH = (105, 125)                   # reads per sgRNA (8 FASTQ files have to stay <= 25 MB)
ZMIN = 6.0                           # a modifier is planted only where a sound analysis has power
MODERATE, STRONG = (0.7, 1.4), (1.5, 3.0)
CLASS = {-1: "sensitiser", 0: "none", 1: "resistance"}

ITEMS = ["library_release", "sample_columns", "stagger_reads", "read2_library",
         "arm_labels", "control_normalisation", "drug_baseline", "contrast_sign"]
SEEDABLE = [i for i in ITEMS if i != "arm_labels"]
RECOMPUTE_ONLY = ["library_release", "sample_columns"]   # no trace in the scripts


# ------------------------------------------------------------------------------------------
# numpy ports of the reference statistics: used only to check that an instance is solvable
# (arms recoverable, release criterion reachable). They never define the truth.
# ------------------------------------------------------------------------------------------
def _phi(x: float) -> float:
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))


def _bh(p: np.ndarray) -> np.ndarray:
    o = np.argsort(p)
    adj = p[o] * len(p) / (np.arange(len(p)) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(len(p))
    out[o] = np.minimum(adj, 1.0)
    return out


def _mcc3(t, o) -> float:
    labs = sorted(set(t) | set(o))
    n = len(t)
    c = sum(1 for a, b in zip(t, o) if a == b)
    tk = [sum(1 for a in t if a == k) for k in labs]
    pk = [sum(1 for b in o if b == k) for k in labs]
    num = c * n - sum(a * b for a, b in zip(tk, pk))
    den = math.sqrt((n * n - sum(a * a for a in pk)) * (n * n - sum(a * a for a in tk)))
    return num / den if den else 0.0


def _norm_ntc(C, ntc):
    tot = C[ntc].sum(axis=0).astype(float)
    return C / (tot / tot.mean())


def _arms_from_data(C, ntc, day14: list[int]) -> list[int]:
    """sign of the first principal component of centred log abundance over the day-14 libraries"""
    X = np.log2(_norm_ntc(C, ntc)[:, day14] + 8.0)
    X = X - X.mean(axis=1, keepdims=True)
    _, _, vt = np.linalg.svd(X, full_matrices=False)
    return [int(v > 0) for v in vt[0]]


def _meanvar_calls(C, ntc, gene_idx, n_genes, D, V):
    N = _norm_ntc(C, ntc)
    L = np.log2(N + 4.0)
    d = L[:, D].mean(axis=1) - L[:, V].mean(axis=1)
    d -= np.median(d[ntc])
    ss = ((L[:, D] - L[:, D].mean(axis=1, keepdims=True)) ** 2).sum(axis=1) + \
         ((L[:, V] - L[:, V].mean(axis=1, keepdims=True)) ** 2).sum(axis=1)
    s2 = ss / (len(D) + len(V) - 2)
    o = np.argsort(L[:, D + V].mean(axis=1))
    sm = np.empty(len(o))
    for i in range(len(o)):
        sm[i] = s2[o][max(0, i - 50):i + 51].mean()
    s2s = np.empty(len(o))
    s2s[o] = sm
    z = d / np.sqrt(s2s * (1 / len(D) + 1 / len(V)))
    gz = np.array([z[gene_idx == g].sum() / math.sqrt((gene_idx == g).sum()) for g in range(n_genes)])
    mad = np.median(np.abs(gz - np.median(gz))) * 1.4826
    gz = (gz - np.median(gz)) / max(mad, 1.0)
    q = _bh(np.array([2 * (1 - _phi(abs(v))) for v in gz]))
    return np.where(q < 0.05, np.sign(gz), 0).astype(int)


# ------------------------------------------------------------------------------------------
# planted truth
# ------------------------------------------------------------------------------------------
def plan(seed: int, sub: int = 0) -> dict:
    rng = np.random.default_rng([40000 + seed, sub])
    n_genes = N_FIT + N_NEU
    gene_ids = sorted(f"TGT{i:04d}" for i in rng.permutation(np.arange(1, 9000))[:n_genes])
    fit_mask = np.zeros(n_genes, dtype=bool)
    fit_mask[rng.permutation(n_genes)[:N_FIT]] = True
    f = np.where(fit_mask, -rng.uniform(1.5, 4.0, size=n_genes), 0.0)
    eff = rng.uniform(0.75, 1.0, size=(n_genes, PER))
    for g in range(n_genes):
        if rng.random() < 0.45:                       # at most one inactive sgRNA per gene
            eff[g, int(rng.integers(0, PER))] = rng.uniform(0.0, 0.08)

    G = n_genes * PER + N_NTC
    gene_idx = np.concatenate([np.repeat(np.arange(n_genes), PER), -np.ones(N_NTC, dtype=int)])
    effv = np.concatenate([eff.ravel(), np.zeros(N_NTC)])
    seqs = _guides(rng, G)
    lib_order = rng.permutation(G)                    # row order of the library file
    gene_idx, effv = gene_idx[lib_order], effv[lib_order]
    seqs = [seqs[i] for i in rng.permutation(G)]
    sg_ids = [f"VAL2_{i + 1:04d}" for i in range(G)]  # neutral ids: the id does not name the gene
    ntc = gene_idx < 0

    base = rng.lognormal(0.0, 0.45, size=G)
    fl = np.where(ntc, 0.0, f[gene_idx] * effv)       # expected sgRNA fitness LFC
    ab14 = base * 2.0 ** fl
    share = ab14 / ab14.sum()
    depth_mean = float(np.mean(DEPTH))
    sb2 = JITTER ** 2 + (BIO * math.log(2)) ** 2

    def expected_z(g: int, d: float) -> float:
        z, idx = 0.0, np.where(gene_idx == g)[0]
        for i in idx:
            dd = d * effv[i]
            rv, cv = depth_mean * G * share[i], COV["vehicle"] * G * share[i]
            rd, cd = rv * 2.0 ** dd, COV["drug"] * G * share[i] * 2.0 ** dd
            var = ((1 / rv + 1 / cv + sb2) / 3 + (1 / rd + 1 / cd + sb2) / 3) / math.log(2) ** 2
            z += dd / math.sqrt(var)
        return abs(z) / math.sqrt(len(idx))

    # drug effects: drawn over all genes, independent of fitness; only where detectable
    d = np.zeros(n_genes)
    cand = list(rng.permutation(n_genes))
    for sign, n in ((-1, N_SENS), (1, N_RES)):
        got = 0
        while got < n:
            g = int(cand.pop(0))
            mag = rng.uniform(*MODERATE) if rng.random() < 0.5 else rng.uniform(*STRONG)
            if expected_z(g, sign * mag) < ZMIN:
                mag = rng.uniform(*STRONG)
                if expected_z(g, sign * mag) < ZMIN:
                    continue
            d[g] = sign * mag
            got += 1
    cls = np.sign(d).astype(int)
    dl = np.where(ntc, 0.0, d[gene_idx] * effv)

    # reference genes (the in-data control of Q)
    sens = [g for g in rng.permutation(n_genes) if cls[g] < 0]
    res = [g for g in rng.permutation(n_genes) if cls[g] > 0]
    inert = [g for g in rng.permutation(n_genes) if cls[g] == 0]
    ref = sorted(int(g) for g in sens[:N_REF_MOD_EACH] + res[:N_REF_MOD_EACH] + inert[:N_REF_INERT])
    is_ref = np.zeros(n_genes, dtype=bool)
    is_ref[ref] = True

    # release 1 of the library annotation: the gene labels of 20 genes were shifted by one
    hid_mod = [int(g) for g in rng.permutation(n_genes) if cls[g] != 0 and not is_ref[g]][:N_V1_PAIRS]
    hid_non = [int(g) for g in rng.permutation(n_genes) if cls[g] == 0 and not is_ref[g]][:N_V1_PAIRS]
    ring = [g for pair in zip(hid_mod, hid_non) for g in pair]
    v1_label = np.arange(n_genes)
    for i, g in enumerate(ring):
        v1_label[g] = ring[(i + 1) % len(ring)]       # sgRNAs of g carry the next gene's name in r1

    # libraries
    first = int(rng.integers(3, 10))                  # ids cross from one to two digits
    sids = [f"LIB{first + i}" for i in range(8)]
    sheet_arm = ["T0", "T0", "vehicle", "vehicle", "vehicle", "drug", "drug", "drug"]
    true_arm = list(sheet_arm)
    sw_v, sw_d = 2 + int(rng.integers(0, 3)), 5 + int(rng.integers(0, 3))
    true_arm[sw_v], true_arm[sw_d] = "drug", "vehicle"   # exchanged before the sheet was written
    rc_j = int(rng.integers(0, 8))

    counts = np.zeros((G, 8), dtype=int)
    n_junk = []
    for j in range(8):
        a = true_arm[j]
        ab = base * rng.lognormal(0.0, JITTER, size=G)
        if a != "T0":
            ab = ab * 2.0 ** (fl + rng.normal(0.0, BIO, size=G))
        if a == "drug":
            ab = ab * 2.0 ** dl
        cells = rng.poisson(COV[a] * G * ab / ab.sum()).astype(float)
        depth = int(G * rng.uniform(*DEPTH))
        counts[:, j] = rng.multinomial(depth, cells / cells.sum())
        n_junk.append(int(depth * rng.uniform(0.03, 0.06)))

    lat_f = np.array([np.median(fl[gene_idx == g]) for g in range(n_genes)])
    lat_d = np.array([np.median(dl[gene_idx == g]) for g in range(n_genes)])

    # which checklist items are defective in the inherited project
    n_seeded = int(rng.integers(2, 5))
    chosen = [RECOMPUTE_ONLY[int(rng.integers(0, 2))]]
    rest = [i for i in SEEDABLE if i not in chosen]
    chosen += [rest[i] for i in rng.permutation(len(rest))[:n_seeded - 1]]
    defects = sorted(set(chosen) | {"arm_labels"}, key=ITEMS.index)

    return dict(seed=seed, sub=sub, n_genes=n_genes, G=G, gene_ids=gene_ids, gene_idx=gene_idx, ntc=ntc,
                sg_ids=sg_ids, seqs=seqs, f=f, d=d, cls=cls, lat_f=lat_f, lat_d=lat_d, is_ref=is_ref,
                v1_label=v1_label, v1_ring=ring, sids=sids, sheet_arm=sheet_arm, true_arm=true_arm,
                swapped=(sw_v, sw_d), rc_j=rc_j, counts=counts, n_junk=n_junk, defects=defects)


def solvable(p: dict) -> tuple[bool, dict]:
    """The stated acceptance criteria must be reachable from the data of this instance."""
    C, ntc = p["counts"], p["ntc"]
    day14 = list(range(2, 8))
    split = _arms_from_data(C, ntc, day14)
    truth = [int(p["true_arm"][j] == "drug") for j in day14]
    arms_ok = split == truth or split == [1 - t for t in truth]
    D = [j for j in day14 if p["true_arm"][j] == "drug"]
    V = [j for j in day14 if p["true_arm"][j] == "vehicle"]
    calls = _meanvar_calls(C, ntc, p["gene_idx"], p["n_genes"], D, V)
    ref, cls = p["is_ref"], p["cls"]
    ref_mod_ok = int(((calls == cls) & (cls != 0) & ref).sum())
    ref_inert_called = int(((calls != 0) & (cls == 0) & ref).sum())
    mcc = _mcc3([int(x) for x in cls[~ref]], [int(x) for x in calls[~ref]])
    info = {"arms_recoverable": arms_ok, "ref_modifiers_right": ref_mod_ok,
            "ref_inert_called": ref_inert_called, "hidden_mcc_numpy_port": round(mcc, 4)}
    return arms_ok and ref_mod_ok >= 2 * N_REF_MOD_EACH - 1 and ref_inert_called <= 1 and mcc >= 0.93, info


# ------------------------------------------------------------------------------------------
# the predecessor's scripts (text templates; the seeded mutations are the alternatives)
# ------------------------------------------------------------------------------------------
COUNT_HEAD = '''#!/usr/bin/env python3
"""count_guides.py - sgRNA read counts for the VAL-2 screens (K.B.)

usage: count_guides.py LIBRARY.tsv SAMPLES.tsv OUT.tsv FASTQ [FASTQ ...]

One FASTQ per row of the sample sheet, given in the same order as the rows.
Only exact matches to a library spacer are counted.
"""
import gzip
import sys

COMP = str.maketrans("ACGTN", "TGCAN")
SCAFFOLD = "GTTTTAGAGCTA"


def revcomp(s):
    return s.translate(COMP)[::-1]

'''
COUNT_CAND = {
    (True, True): '''
def candidates(read):
    # the spacer is the 20 nt in front of the scaffold; run B delivered Read 2, so look on both strands
    for r in (read, revcomp(read)):
        k = r.find(SCAFFOLD)
        if k >= 20:
            yield r[k - 20:k]
''',
    (False, True): '''
def candidates(read):
    yield read[24:44]             # Read 1: the spacer follows the 24 nt of the U6 primer site
    yield revcomp(read)[15:35]    # run B delivered Read 2: same spacer seen from the other end
''',
    (True, False): '''
def candidates(read):
    # the spacer is the 20 nt in front of the scaffold
    k = read.find(SCAFFOLD)
    if k >= 20:
        yield read[k - 20:k]
''',
    (False, False): '''
def candidates(read):
    yield read[24:44]             # the spacer follows the 24 nt of the U6 primer site
''',
}
COUNT_TAIL = '''

def main():
    lib_path, sheet_path, out_path = sys.argv[1:4]
    fastqs = sys.argv[4:]
    guides, by_seq = [], {}
    with open(lib_path) as fh:
        head = fh.readline().rstrip("\\n").split("\\t")
        i_id, i_seq, i_gene = head.index("sgRNA"), head.index("sequence"), head.index("gene")
        for line in fh:
            f = line.rstrip("\\n").split("\\t")
            by_seq[f[i_seq]] = len(guides)
            guides.append((f[i_id], f[i_gene]))
    with open(sheet_path) as fh:
        fh.readline()
        samples = [line.split("\\t")[0] for line in fh if line.strip()]
    if len(fastqs) != len(samples):
        sys.exit("need one FASTQ per sample sheet row")
    counts = [[0] * len(samples) for _ in guides]
    for j, path in enumerate(fastqs):
        with gzip.open(path, "rt") as fh:
            for n, line in enumerate(fh):
                if n % 4 != 1:
                    continue
                for spacer in candidates(line.rstrip()):
                    k = by_seq.get(spacer)
                    if k is not None:
                        counts[k][j] += 1
                        break
    with open(out_path, "w") as out:
        out.write("sgRNA\\tgene\\t" + "\\t".join(samples) + "\\n")
        for (sg, gene), row in zip(guides, counts):
            out.write(sg + "\\t" + gene + "\\t" + "\\t".join(map(str, row)) + "\\n")


if __name__ == "__main__":
    main()
'''

ANALYSE_HEAD = '''#!/usr/bin/env python3
"""analyse_screen.py - gene table and hit list for the VAL-2 drug screen (K.B.)

usage: analyse_screen.py COUNTS.tsv SAMPLES.tsv OUT_GENES.tsv OUT_HITS.tsv
"""
import sys

import numpy as np

TOP = __TOP__                             # genes that go into the hit list
counts_path, sheet_path, out_genes, out_hits = sys.argv[1:5]

rows = [line.rstrip("\\n").split("\\t") for line in open(counts_path) if line.strip()]
cols, rows = rows[0][2:], rows[1:]
gene = np.array([r[1] for r in rows])
C = np.array([[float(v) for v in r[2:]] for r in rows])
arm = {}
for line in list(open(sheet_path))[1:]:
    f = line.rstrip("\\n").split("\\t")
    if len(f) > 1:
        arm[f[0]] = f[1]

keep = [j for j in range(len(cols)) if C[:, j].sum() > 0]        # a library without reads is left out
t0 = [j for j in keep if arm[cols[j]] == "T0"]
veh = [j for j in keep if arm[cols[j]] == "vehicle"]
drug = [j for j in keep if arm[cols[j]] == "drug"]
ntc = gene == "NON_TARGETING"
'''
ANALYSE_NORM = {
    True: '''
# size factors from the non-targeting sgRNAs
size = np.ones(len(cols))
tot = C[ntc].sum(axis=0)
size[keep] = tot[keep] / tot[keep].mean()
N = C / size
''',
    False: '''
# size factors: median of ratios to the geometric mean, all sgRNAs
size = np.ones(len(cols))
with np.errstate(divide="ignore", invalid="ignore"):
    logc = np.log(np.where(C[:, keep] > 0, C[:, keep], np.nan))
    size[keep] = np.exp(np.nanmedian(logc - np.nanmean(logc, axis=1, keepdims=True), axis=0))
N = C / size
''',
}
ANALYSE_LFC = '''

def lfc(a, b):
    return np.log2((N[:, a].mean(axis=1) + 0.5) / (N[:, b].mean(axis=1) + 0.5))


fit = lfc(veh, t0)
drg = __DRUG__
'''
ANALYSE_CENTRE = {
    True: '''fit = fit - np.median(fit[ntc])
drg = drg - np.median(drg[ntc])
''',
    False: "",
}
ANALYSE_DRUG = {   # (baseline is vehicle, sign is drug-minus-reference)
    (True, True): "lfc(drug, veh)",
    (False, True): "lfc(drug, t0)",
    (True, False): "lfc(veh, drug)",
    (False, False): "lfc(t0, drug)",
}
ANALYSE_TAIL = '''
# gene level: median over the sgRNAs; z against the spread of the non-targeting sgRNAs
null_sd = np.std(drg[ntc]) / 2.0          # a gene is the median of four sgRNAs
table = []
for g in sorted(set(gene[~ntc])):
    m = gene == g
    lf, ld = float(np.median(fit[m])), float(np.median(drg[m]))
    z = ld / null_sd
    cls = "none" if abs(z) < 3 else ("sensitiser" if z < 0 else "resistance")
    table.append((g, lf, ld, z, cls))

with open(out_genes, "w") as out:
    out.write("gene\\tlfc_fitness\\tlfc_drug\\tz_drug\\tdrug_class\\n")
    for g, lf, ld, z, cls in table:
        out.write("%s\\t%.4f\\t%.4f\\t%.2f\\t%s\\n" % (g, lf, ld, z, cls))

with open(out_hits, "w") as out:
    out.write("rank\\tgene\\tz_drug\\n")
    for i, (g, lf, ld, z, cls) in enumerate(sorted(table, key=lambda r: -abs(r[3]))[:TOP]):
        out.write("%d\\t%s\\t%.2f\\n" % (i + 1, g, z))
'''


def count_script(defects: list[str]) -> str:
    return COUNT_HEAD + COUNT_CAND[("stagger_reads" not in defects, "read2_library" not in defects)].lstrip("\n") + COUNT_TAIL


def analyse_script(defects: list[str], top: int) -> str:
    norm_ok = "control_normalisation" not in defects
    drug = ANALYSE_DRUG[("drug_baseline" not in defects, "contrast_sign" not in defects)]
    return (ANALYSE_HEAD.replace("__TOP__", str(top)) + ANALYSE_NORM[norm_ok] + ANALYSE_LFC.replace("__DRUG__", drug)
            + ANALYSE_CENTRE[norm_ok] + ANALYSE_TAIL)


# ------------------------------------------------------------------------------------------
# texts
# ------------------------------------------------------------------------------------------
INSTRUCTION = """\
You are the computational biologist of a functional-genomics core facility. You are taking over a
pooled CRISPR-Cas9 knockout drug-modifier screen (secondary validation library VAL-2) from a
colleague who has left. He had started the analysis and mailed a draft hit list to the PI. You
deliver the final results, an audit of the inherited project and an erratum for the mailed list.

Task directory:
- `/task`

Inputs (read-only) under `/task/input`:
- `fastq/<sample_id>.fastq.gz` - one single-end FASTQ per library ({n_samples} libraries)
- `samples.tsv` - sample sheet as it was submitted with the libraries (`sample_id`, `arm`,
  `replicate`, `fastq_file`, `sequencing_run`); `arm` is `T0`, `vehicle` or `drug`
- `library.tsv` - current release (r2) of the sgRNA library: `sgRNA`, `sequence`, `gene`.
  Non-targeting control sgRNAs have `gene = NON_TARGETING`. Gene identities were blinded by the client.
- `reference_genes.tsv` - library genes whose behaviour under this compound is established
  (`gene`, `expected_class`)
- `screen_design.md` - wet-lab notes on the library, the treatment arms, PCR and sequencing
- `screen_sop.md` - the facility SOP for drug-modifier screens: definitions, acceptance criteria,
  the hand-over audit checklist
- `handover/` - the colleague's project folder: `HANDOVER.md`, `count_guides.py`,
  `analyse_screen.py`, `counts_draft.tsv`, `gene_draft.tsv`, `hits_mailed.tsv`

Available software (on `PATH`): `cutadapt`, `mageck` (MAGeCK), Python 3.12 with numpy, scipy and
pandas. There is no internet access.

Your task:
1. Count how many reads carry each library sgRNA in every library.
2. For every targeted gene, estimate the fitness effect of the knockout and its effect on the
   response to the drug, and classify the gene as sensitiser, resistance gene or neither.
3. Audit the inherited project against the checklist in section 6 of the SOP.
4. Write the erratum for the hit list that was mailed.

Required outputs under `/task/output`:
- `counts.tsv` - tab-separated, header `sgRNA`, `gene`, then one column per `sample_id` in the
  order of `samples.tsv`. One row per sgRNA of `library.tsv` (non-targeting sgRNAs included), `gene`
  as in `library.tsv`. Values are raw integer read counts of the FASTQ file listed for that sample.
- `gene_summary.tsv` - tab-separated with at least the columns `gene`, `lfc_fitness`, `lfc_drug`,
  `drug_class`; one row per targeted gene (no row for `NON_TARGETING`). Extra columns are allowed.
  - `lfc_fitness`: log2 fold change of the gene, vehicle-treated libraries of day 14 versus T0,
    defined as the median over the gene's sgRNAs of the sgRNA-level log2 ratio of depth-normalised
    abundance (mean of the replicates of one group over mean of the replicates of the other). The
    scale must be such that a knockout without any fitness effect has `lfc_fitness` = 0.
  - `lfc_drug`: the same quantity for drug-treated versus vehicle-treated libraries of day 14;
    negative means less abundant under drug. A knockout that does not change the drug response
    has `lfc_drug` = 0.
  - `drug_class`: `sensitiser` if the gene's sgRNAs are significantly less abundant under drug than
    under vehicle, `resistance` if they are significantly more abundant, otherwise `none`; false
    discovery rate 5%. The calls must satisfy the release criterion of SOP section 5.
  - Libraries are grouped by the treatment they actually received (SOP section 3).
- `audit.tsv` - tab-separated, columns `item`, `verdict`; one row per checklist item of SOP
  section 6, `verdict` is `ok` or `defective`.
- `erratum.tsv` - tab-separated, columns `gene`, `verdict`; one row per gene of
  `handover/hits_mailed.tsv`. `verdict` is `confirmed` if your final `drug_class` of the gene is
  `sensitiser` or `resistance`, and `retracted` if it is `none`.
- `qc.json` - `{{"samples": {{"<sample_id>": {{"total_reads": <int>, "mapped_reads": <int>, "arm": "<T0|vehicle|drug>"}}, ...}}}}`
  where `total_reads` is the number of reads in the FASTQ, `mapped_reads` the number of reads
  assigned to a library sgRNA and `arm` the group in which you analysed the library.

Constraints:
- Do not modify anything under `/task/input`; write deliverables only under `/task/output`.
- Every library in the sample sheet belongs to the experiment and must be reported and used.
- Do not ask for confirmation. Execute directly.
"""

SOP = """\
# SOP FG-017: pooled CRISPR drug-modifier screens, analysis and release (rev. 4)

## 1. Scope
Secondary (validation) knockout screens with a T0 reference, a vehicle arm and a drug arm. This
SOP states what a released analysis has to satisfy. It does not prescribe software.

## 2. Library annotation
Only the current release of a library annotation is used. For VAL-2 this is release 2
(`library.tsv`). Release 2 corrected the gene assignment of {n_v1} sgRNAs ({n_v1_genes} genes) that
were mislabelled in release 1 (a sub-pool of the oligo order sheet was shifted against its gene
column); sgRNA ids and sequences did not change. Release 1 is withdrawn.

## 3. Counting and sample identity
- A read is counted for an sgRNA if it carries that sgRNA's spacer in the amplicon context
  described in the wet-lab notes. Every library of the sample sheet is counted and reported.
- The sample sheet is what the wet lab wrote down. Before any contrast is computed, the arm of
  every day-14 library is confirmed from the data: the sgRNA profile of a library has to agree
  better with the other libraries of its arm than with the libraries of the other arm. A library
  that does not is analysed in the arm its profile supports. The arm in which each library was
  analysed is recorded in `qc.json`.

## 4. Scales
Log2 fold changes are reported on the scale on which the non-targeting sgRNAs do not change
(median non-targeting sgRNA = 0). The fitness effect is vehicle day 14 against T0. The drug effect
is drug day 14 against vehicle day 14; negative values mean depletion under drug.

## 5. Release criterion for modifier calls
The library carries reference genes for the compound class (`reference_genes.tsv`): modifiers with
their expected class and genes established not to modify the response (`none`); the reference
`none` genes cover the whole range of fitness effects. A call set is released only if
- at least {min_ref_mod} of the {n_ref_mod} reference modifiers are called with their expected class, and
- at most {max_ref_inert} of the {n_ref_inert} reference `none` genes is called (either class).
The reference genes are reported like every other gene.

## 6. Hand-over audit checklist
When a project changes hands, the new analyst checks the inherited files and records one verdict
per item in `audit.tsv`: `ok` if the statement holds, `defective` if it does not. Each item is
judged on its own: an effect that follows from the defect of another item is recorded under the
item that caused it, not under every item it touches.

| item | statement |
|---|---|
| `library_release` | The sgRNA-to-gene annotation in the inherited count table is that of the current library release. |
| `sample_columns` | No sample column of the inherited count table was counted from a FASTQ file other than the one the sample sheet lists for that sample. (A column without any counts is the subject of `read2_library`, not of this item.) |
| `stagger_reads` | The inherited count table contains the reads of all stagger variants of the amplicon, not only of one read layout. |
| `read2_library` | The library that was delivered as Read 2 has been counted: no sample column of the inherited count table is empty. |
| `arm_labels` | The inherited analysis script takes the arm of each library from the sample sheet; every library is thereby analysed in the arm of the treatment it actually received. |
| `control_normalisation` | The inherited analysis script reports log2 fold changes on the scale of section 4 (non-targeting sgRNAs at 0). |
| `drug_baseline` | The drug effect in the inherited analysis script compares the drug arm with the vehicle arm of day 14, not with T0. |
| `contrast_sign` | In the inherited analysis script a negative drug effect means depletion under drug. |

The inherited gene table and hit list are the output of the inherited analysis script run on the
inherited count table with the sample sheet.

## 7. Erratum
If results have left the facility, an erratum lists every gene of the mailed list with one verdict:
`confirmed` if the gene is a modifier (either class) in the released analysis, `retracted` if not.
"""

HANDOVER = """\
# Hand-over, VAL-2 drug screen (K.B.)

Sorry for the state of this, my last week was short.

- `count_guides.py` - my counting script. `counts_draft.tsv` is its output for the eight libraries.
- `analyse_screen.py` - gene table and hit list from a count table and the sample sheet.
  `gene_draft.tsv` and `hits_mailed.tsv` are its output for `counts_draft.tsv`.
- `hits_mailed.tsv` went to the PI by mail on Friday as "top {top} modifier candidates". Nobody has
  ordered reagents on it yet as far as I know.

Reads assigned per library in `counts_draft.tsv`:

{assigned}

Not done: the SOP checks, the reference genes, anything on the re-sequenced library beyond counting it
the way the script does. I used the dropout workflow we always use and adapted it for two arms.
"""

JUDGE = """\
**Where the truth comes from.** The generator fixes, before any read exists: each gene's fitness
effect and, independently, its drug effect (sensitiser, resistance or none); each sgRNA's
efficiency; the treatment every library really received (one vehicle and one drug library are
exchanged on the sample sheet); how many reads each sgRNA contributes to each library; which
checklist items are defective in the inherited project. `reference/truth.json` stores these planted
values. The inherited files are the output of the predecessor's scripts (seeded mutations of a
correct pipeline) run on the task data; the audit truth is the list of mutations, the erratum
truth is the planted class of every mailed gene. No tool output is used as reference. A numpy port
of the reference statistics is run at generation time only to confirm that the stated acceptance
criteria can be met (arms recoverable, SOP section 5 reachable); it never defines a truth value.
A drug effect is planted only where a sound analysis has power (expected gene z >= 6), so no scored
gene is a borderline case.

**Hard gates (score 0).** Missing/unparseable file; `counts.tsv` not listing every library sgRNA
exactly once with non-negative integers; `gene_summary.tsv` not listing every targeted gene exactly
once; `drug_class` outside {sensitiser, resistance, none}; `audit.tsv` not listing every checklist
item exactly once with `ok`/`defective`; `erratum.tsv` not listing every mailed gene exactly once
with `confirmed`/`retracted`; non-finite numbers.

**Components.**

| weight | component | metric | full credit | zero |
|---|---|---|---|---|
| 0.15 | guide_counts | per library `sum(min(obs,true)) / max(sum obs, sum true)`, then the **minimum** over libraries | >= 0.90 | <= 0.50 |
| 0.20 | gene_lfc | the largest of three median absolute errors against the planted latent values: `lfc_fitness` (all genes), `lfc_drug` (planted non-modifiers), `lfc_drug` (planted modifiers) | <= 0.25 | >= 0.60 |
| 0.30 | drug_class | three-class Matthews correlation with the planted class, over the genes that are **not** reference genes | >= 0.90 | <= 0.75 |
| 0.20 | handover_audit | Matthews correlation of `defective` verdicts with the planted defects (8 items) | 1.0 | <= 0.50 |
| 0.15 | erratum | balanced accuracy (mean of sensitivity and specificity) of the verdicts against "mailed gene is a planted modifier" | >= 0.85 | <= 0.60 |

`qc.json` is required and format-gated (an entry per library); its content is written to the
verifier details (`qc_arms_right`, `qc_read_numbers_right`) but carries no weight: the arms are
scored through `lfc_drug` and `drug_class`, the read numbers through `guide_counts`.

**Why these tolerances.** Counts and the LFC band are those of level 3 (exact-match counting keeps
about 96% of the reads; correct LFC estimates have median errors of 0.08-0.2, the wrong
normalisation, the wrong baseline, the sample-sheet arms or a reversed contrast give 0.7-2.4). The
modifier subset is scored separately because the median over all genes is blind to a reversed or
diluted drug contrast (85% of the genes have a drug effect of 0). For `drug_class` the stated
requirement is a false discovery rate of 5% (one or two false calls among ~26 scored modifiers,
MCC >= 0.96). Full credit starts at MCC 0.90, i.e. about five false calls or four misses, three
times the nominal error; two different abundance-aware tests reach 0.96-1.0 and MAGeCK RRA on the
right contrast 0.89-1.0 (it is a must-pass probe with floor 0.95). Calls against the
spread of the non-targeting sgRNAs give 20-55 false calls (MCC 0.55-0.8) and the level-3 rank test
both over- and under-calls (0.65-0.9); the band is put between the two groups. Reference genes are
excluded from this component because their class is given to the agent. Audit verdicts are facts,
so only a complete audit earns full credit; one wrong item costs about 0.09 of the total. The erratum
uses balanced accuracy because a mailed list can contain as few as four true modifiers (the generator
guarantees at least four of each verdict), where MCC moves by 0.1 per verdict; "all confirmed" and
"all retracted" both score exactly 0.5, i.e. zero credit. Every component is therefore blind to its
trivial answers ("all ok", "all defective", "all retracted", "nothing called").

**What is deliberately not scored.** P-values, the method used, the reference genes' own classes
(they are the agent's control), free-text explanations.
"""


def _fmt_library(p: dict, label: np.ndarray | None = None) -> str:
    out = ["sgRNA\tsequence\tgene\n"]
    for i in range(p["G"]):
        g = p["gene_idx"][i]
        if g < 0:
            name = "NON_TARGETING"
        else:
            name = p["gene_ids"][int(label[g]) if label is not None else int(g)]
        out.append(f"{p['sg_ids'][i]}\t{p['seqs'][i]}\t{name}\n")
    return "".join(out)


def _run(cmd: list[str]) -> None:
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError("predecessor script failed: " + " ".join(cmd) + "\n" + r.stderr[-2000:])


def generate(seed: int, level: int, workdir: Path) -> TaskBuild:
    # ---- planted truth; re-draw (deterministically) until the stated criteria are reachable ----
    for sub in range(12):
        p = plan(seed, sub)
        ok, solv = solvable(p)
        if ok:
            break
    else:
        raise RuntimeError(f"crispr_screen level 4 seed {seed}: no solvable instance in 12 draws")
    G, n_genes, sids, gene_ids = p["G"], p["n_genes"], p["sids"], p["gene_ids"]
    counts, defects = p["counts"], p["defects"]
    workdir.mkdir(parents=True, exist_ok=True)

    # ---- reads -------------------------------------------------------------------------------
    rng = np.random.default_rng([41000 + seed, p["sub"]])
    fq_dir = workdir / "fastq"
    fq_dir.mkdir(parents=True, exist_ok=True)
    fwd = [PREFIX + s + SCAFFOLD for s in p["seqs"]]
    rev = [_rc(PREFIX + s + SCAFFOLD[:40]) for s in p["seqs"]]
    total_reads = []
    for j, sid in enumerate(sids):
        origin = np.concatenate([np.repeat(np.arange(G), counts[:, j]), -np.ones(p["n_junk"][j], dtype=int)])
        origin = origin[rng.permutation(len(origin))]
        reverse = j == p["rc_j"]
        _write_fastq(fq_dir / f"{sid}.fastq.gz", origin, fwd, rev, reverse, True, rng, "B" if reverse else "A")
        total_reads.append(int(len(origin)))
    guide_reads = [int(v) for v in counts.sum(axis=0)]

    library_tsv = _fmt_library(p)
    library_r1_tsv = _fmt_library(p, p["v1_label"])
    samples_tsv = "sample_id\tarm\treplicate\tfastq_file\tsequencing_run\n"
    rep = {"T0": 0, "vehicle": 0, "drug": 0}
    for j, sid in enumerate(sids):
        a = p["sheet_arm"][j]
        rep[a] += 1
        samples_tsv += f"{sid}\t{a}\t{rep[a]}\t{sid}.fastq.gz\t{'B' if j == p['rc_j'] else 'A'}\n"
    ref_tsv = "gene\texpected_class\n" + "".join(
        f"{gene_ids[g]}\t{CLASS[int(p['cls'][g])]}\n" for g in range(n_genes) if p["is_ref"][g])

    # ---- the inherited project: run the predecessor's (mutated) scripts on the task data -----
    ho = workdir / "handover"
    ho.mkdir(parents=True, exist_ok=True)
    (workdir / "library.tsv").write_text(library_tsv)
    (workdir / "library_r1.tsv").write_text(library_r1_tsv)
    (workdir / "samples.tsv").write_text(samples_tsv)
    (ho / "count_guides.py").write_text(count_script(defects))
    files_sheet = [str(fq_dir / f"{s}.fastq.gz") for s in sids]
    files_used = sorted(files_sheet) if "sample_columns" in defects else files_sheet   # shell glob order
    lib_used = workdir / ("library_r1.tsv" if "library_release" in defects else "library.tsv")
    _run([sys.executable, str(ho / "count_guides.py"), str(lib_used), str(workdir / "samples.tsv"),
          str(ho / "counts_draft.tsv")] + files_used)
    truth_cls = {gene_ids[g]: int(p["cls"][g]) for g in range(n_genes)}
    for top in (40, 60, 80, 120):
        (ho / "analyse_screen.py").write_text(analyse_script(defects, top))
        _run([sys.executable, str(ho / "analyse_screen.py"), str(ho / "counts_draft.tsv"),
              str(workdir / "samples.tsv"), str(ho / "gene_draft.tsv"), str(ho / "hits_mailed.tsv")])
        mailed = [line.split("\t")[1] for line in (ho / "hits_mailed.tsv").read_text().splitlines()[1:]]
        n_conf = sum(1 for g in mailed if truth_cls[g] != 0)
        if n_conf >= 4 and len(mailed) - n_conf >= 4:
            break
    else:
        raise RuntimeError(f"seed {seed}: mailed list has a single verdict class ({n_conf}/{len(mailed)})")
    draft_rows = (ho / "counts_draft.tsv").read_text().splitlines()
    draft_tot = [sum(int(r.split("\t")[2 + j]) for r in draft_rows[1:]) for j in range(8)]
    assigned = "\n".join(f"- {sid}: {draft_tot[j]}" for j, sid in enumerate(sids))
    handover_md = HANDOVER.format(top=top, assigned=assigned)

    rc_sample = sids[p["rc_j"]]
    notes = "\n".join([
        "# Screen design notes (wet lab)", "",
        "- Cells: Cas9-stable clone (identity blinded by the client), lentiviral infection at MOI ~0.3, "
        "500x library coverage until the arms were split.",
        f"- Library: VAL-2 secondary validation library in lentiGuide-Puro. {N_FIT} genes carried forward as "
        f"candidate fitness genes from the primary genome-wide dropout screen and {N_NEU} genes picked as putative "
        f"neutral; {PER} sgRNAs per gene, plus {N_NTC} non-targeting sgRNAs. The reference genes for the compound "
        "class are part of these genes.",
        "- Arms: T0 = day 3 after puromycin selection (2 replicates). The culture was then split into a vehicle "
        "arm (DMSO, 3 replicates) and a drug arm (compound at IC70, 3 replicates), both harvested 14 days later.",
        f"- Harvest: vehicle arm about {int(COV['vehicle'] * 1000 / 1000) * 1000} cells per replicate went into "
        f"the gDNA prep (about {int(COV['vehicle'])} cells per sgRNA). The drug arm was limited by the treatment: all "
        f"surviving cells were taken, about {int(COV['drug']) * 1000} per replicate (about {int(COV['drug'])} per sgRNA).",
        "- gDNA PCR: one-step PCR with the P5 stagger primer mix (8 forward primers, 0-7 nt stagger) and one "
        "P7 index primer per library.",
        f"- Amplicon, sense strand 5'->3': ...{PREFIX} [20-nt spacer] {SCAFFOLD[:30]}... (sgRNA scaffold).",
        "- Sequencing run A: NextSeq 550, 75-cycle single read.",
        f"- {rc_sample} failed QC on run A. The core re-sequenced the same PCR product on run B (NovaSeq, "
        "paired-end lane shared with other users); only the Read 2 file was delivered and that is the file "
        "provided here.",
        "- The day-14 gDNA tubes were re-racked when the -20 freezer was replaced; the sample sheet was "
        "filled in from the new rack positions.",
    ]) + "\n"
    sop_md = SOP.format(n_v1=2 * N_V1_PAIRS * PER, n_v1_genes=2 * N_V1_PAIRS, n_ref_mod=2 * N_REF_MOD_EACH,
                        min_ref_mod=2 * N_REF_MOD_EACH - 1, n_ref_inert=N_REF_INERT, max_ref_inert=1)

    # ---- hidden reference --------------------------------------------------------------------
    audit_truth = {it: ("defective" if it in defects else "ok") for it in ITEMS}
    erratum_truth = {g: ("confirmed" if truth_cls[g] != 0 else "retracted") for g in mailed}
    genes_truth = {gene_ids[g]: {"lfc_fitness": float(p["lat_f"][g]), "lfc_drug": float(p["lat_d"][g]),
                                 "class": CLASS[int(p["cls"][g])], "reference": bool(p["is_ref"][g])}
                   for g in range(n_genes)}
    lib_gene = ["NON_TARGETING" if g < 0 else gene_ids[int(g)] for g in p["gene_idx"]]
    truth = {"samples": sids, "guides": dict(zip(p["sg_ids"], lib_gene)),
             "true_counts": {sg: [int(v) for v in counts[i]] for i, sg in enumerate(p["sg_ids"])},
             "genes": genes_truth, "arms": dict(zip(sids, p["true_arm"])), "audit": audit_truth,
             "erratum": erratum_truth, "total_reads": total_reads, "guide_reads": guide_reads,
             "provenance": {"sub_draw": p["sub"], "solvability_check": solv, "defects": defects,
                            "swapped_libraries": [sids[p["swapped"][0]], sids[p["swapped"][1]]],
                            "read2_library": rc_sample, "mailed_top": top}}

    # ---- truth-derived probes / attacks ------------------------------------------------------
    def tsv_counts(mat):
        return "sgRNA\tgene\t" + "\t".join(sids) + "\n" + "".join(
            f"{sg}\t{lib_gene[i]}\t" + "\t".join(str(int(v)) for v in mat[i]) + "\n" for i, sg in enumerate(p["sg_ids"]))

    def tsv_genes(lf=lambda g: genes_truth[g]["lfc_fitness"], ld=lambda g: genes_truth[g]["lfc_drug"],
                  cl=lambda g: genes_truth[g]["class"]):
        return "gene\tlfc_fitness\tlfc_drug\tdrug_class\n" + "".join(
            f"{g}\t{lf(g):.4f}\t{ld(g):.4f}\t{cl(g)}\n" for g in gene_ids)

    def tsv_audit(fn=lambda it: audit_truth[it]):
        return "item\tverdict\n" + "".join(f"{it}\t{fn(it)}\n" for it in ITEMS)

    def tsv_erratum(fn=lambda g: erratum_truth[g]):
        return "gene\tverdict\n" + "".join(f"{g}\t{fn(g)}\n" for g in mailed)

    def qc_json(mapped, arms=None):
        arms = arms or p["true_arm"]
        return json.dumps({"samples": {s: {"total_reads": total_reads[j], "mapped_reads": int(mapped[j]),
                                           "arm": arms[j]} for j, s in enumerate(sids)}})

    def truth_files(**over):
        files = {"counts.tsv": tsv_counts(counts), "gene_summary.tsv": tsv_genes(), "audit.tsv": tsv_audit(),
                 "erratum.tsv": tsv_erratum(), "qc.json": qc_json(guide_reads)}
        files.update(over)
        return files

    prng = np.random.default_rng(seed + 4077)
    thinned = prng.binomial(counts, 0.955)
    noise_f = {g: float(prng.normal(0, 0.08)) for g in gene_ids}
    noise_d = {g: float(prng.normal(0, 0.08)) for g in gene_ids}
    hid_non = [g for g in gene_ids if genes_truth[g]["class"] == "none" and not genes_truth[g]["reference"]]
    hid_mod = [g for g in gene_ids if genes_truth[g]["class"] != "none" and not genes_truth[g]["reference"]]
    fp2, fn1 = set(hid_non[:2]), set(hid_mod[:1])
    fp10, fn8 = set(hid_non[:10]), set(hid_mod[:8])

    def cls_err(fp, fn):
        def f(g):
            if g in fp:
                return "sensitiser"
            if g in fn:
                return "none"
            return genes_truth[g]["class"]
        return f

    wrong_one = [g for g in mailed if g not in fp2 and g not in fn1][0]
    flip = {"confirmed": "retracted", "retracted": "confirmed"}
    pipeline_py = {"crispr_l4_pipeline.py": ASSETS / "crispr_l4_pipeline.py", "library_r1.tsv": library_r1_tsv}
    run = "python3 /alegen_script/crispr_l4_pipeline.py "

    probes = [
        Probe("planted_truth", "submission written directly from the planted truth", 0.999, files=truth_files()),
        Probe("within_tolerance", "4.5% of reads lost, LFC noise sd 0.08 on both scales, two false and one missed "
              "modifier among the scored genes, one wrong erratum verdict", 0.99,
              files=truth_files(**{
                  "counts.tsv": tsv_counts(thinned),
                  "gene_summary.tsv": tsv_genes(lambda g: genes_truth[g]["lfc_fitness"] + noise_f[g],
                                                lambda g: genes_truth[g]["lfc_drug"] + noise_d[g], cls_err(fp2, fn1)),
                  "erratum.tsv": tsv_erratum(lambda g: flip[erratum_truth[g]] if g == wrong_one else erratum_truth[g]),
                  "qc.json": qc_json(thinned.sum(axis=0))})),
        Probe("independent_solution", "second implementation without cutadapt, MAGeCK or scipy: pure-python "
              "counting, arms from the reference-gene signature, calls from the spread of sgRNAs of similar abundance "
              "(no replicate variance model), the inherited tables reproduced numerically instead of re-running the "
              "inherited scripts", 0.99, script="SOLUTION=2 bash /solution/solve.sh"),
        Probe("mageck_rra_calls", "third route to the calls: the reference pipeline with the calling step replaced "
              "by MAGeCK RRA (mageck test, control normalisation) on the drug-versus-vehicle contrast with the arms "
              "from the data. Floor 0.95 instead of 0.99: RRA estimates its variance from the vehicle arm and the drug "
              "arm is noisier, so its FDR is slightly anti-conservative (0-6 false calls measured, MCC 0.89-1.0)",
              0.95, script=run + "--mut call_mageck", script_files=pipeline_py),
    ]
    attacks = [
        Attack("l3_reference", "level-3 reference method on the level-4 inputs: correct counts, control "
               "normalisation, every contrast against T0 with the sample-sheet arms, sensitiser = level-3 hit in the "
               "drug arm that is no hit in the vehicle arm, no audit (all ok), nothing retracted", 0.75,
               script=run + "--l3", script_files=pipeline_py),
        Attack("trust_inherited", "ignores H: finishes from the inherited count table with the reference "
               "statistics, audit all ok", 0.75, script=run + "--mut inherited_counts", script_files=pipeline_py),
        Attack("trust_sample_sheet", "ignores P(a): everything as the reference, but the arms are taken from the "
               "sample sheet", 0.75, script=run + "--mut sheet_arms", script_files=pipeline_py),
        Attack("drug_vs_t0", "ignores P(b): everything as the reference, but the drug effect is the drug arm "
               "against T0", 0.75, script=run + "--mut baseline_t0", script_files=pipeline_py),
        Attack("calls_vs_nontargeting_spread", "ignores Q: everything as the reference, but genes are called "
               "against the spread of the non-targeting sgRNAs (the level-3 second solution's null)", 0.75,
               script=run + "--mut call_ntc_z", script_files=pipeline_py),
        Attack("calls_l3_rank_test", "ignores Q, second variant: the level-3 reference test (Mann-Whitney "
               "against the non-targeting sgRNAs), two-sided on the right contrast", 0.89,
               script=run + "--mut call_mwu", script_files=pipeline_py),
        Attack("unfixed_library_release", "one defect operator left in: gene annotation of release 1", 0.75,
               script=run + "--mut library_r1", script_files=pipeline_py),
        Attack("unfixed_sample_columns", "one defect operator left in: FASTQ files taken in lexicographic order", 0.75,
               script=run + "--mut file_order", script_files=pipeline_py),
        Attack("unfixed_stagger_reads", "one defect operator left in (level-3 trap): fixed-offset parsing", 0.86,
               script=run + "--mut fixed_offset", script_files=pipeline_py),
        Attack("unfixed_read2_library", "one defect operator left in (level-3 trap): the Read-2 library is not "
               "rescued, reported as zeros and left out of the statistics", 0.86,
               script=run + "--mut no_revcomp", script_files=pipeline_py),
        Attack("unfixed_control_normalisation", "one defect operator left in (level-3 trap): median-of-ratios "
               "normalisation over the whole library", 0.81, script=run + "--mut median_norm", script_files=pipeline_py),
        Attack("unfixed_contrast_sign", "one defect operator left in: vehicle against drug", 0.75,
               script=run + "--mut reversed", script_files=pipeline_py),
        Attack("naive_default_pipeline", "fully careless: forward strand only (the Read-2 library is lost and left "
               "out), sample-sheet arms, median-of-ratios normalisation, calls against the non-targeting spread, no "
               "audit (all ok)", 0.50,
               script=run + "--mut no_revcomp --mut sheet_arms --mut median_norm --mut call_ntc_z --mut no_audit",
               script_files=pipeline_py),
        Attack("audit_all_defective", "correct results, every checklist item reported defective", 0.81,
               files=truth_files(**{"audit.tsv": tsv_audit(lambda it: "defective")})),
        Attack("audit_all_ok", "correct results, every checklist item reported ok", 0.81,
               files=truth_files(**{"audit.tsv": tsv_audit(lambda it: "ok")})),
        Attack("erratum_nothing_retracted", "correct results, every mailed gene confirmed", 0.86,
               files=truth_files(**{"erratum.tsv": tsv_erratum(lambda g: "confirmed")})),
        Attack("prior_only", "no analysis: counts all 1, every gene a sensitiser with lfc -3 / -1, audit all "
               "defective, everything retracted", 0.05,
               files={"counts.tsv": tsv_counts(np.ones_like(counts)),
                      "gene_summary.tsv": tsv_genes(lambda g: -3.0, lambda g: -1.0, lambda g: "sensitiser"),
                      "audit.tsv": tsv_audit(lambda it: "defective"),
                      "erratum.tsv": tsv_erratum(lambda g: "retracted"),
                      "qc.json": qc_json([1] * 8, p["sheet_arm"])}),
        Attack("lfc_drug_shifted_beyond_tolerance", "planted truth, every lfc_drug shifted by +0.7", 0.81,
               files=truth_files(**{"gene_summary.tsv": tsv_genes(ld=lambda g: genes_truth[g]["lfc_drug"] + 0.7)})),
        Attack("calls_beyond_tolerance", "planted truth with ten false and eight missed modifiers among the "
               "scored genes", 0.75,
               files=truth_files(**{"gene_summary.tsv": tsv_genes(cl=cls_err(fp10, fn8))})),
    ]

    traps = [
        Trap("staggered_amplicon", "The spacer starts 24-31 nt into the read; a fixed-offset parser assigns ~1/8 of reads.",
             "screen_design.md names the P5 stagger mix; per-library mapping rate collapses with a fixed offset",
             "Approach / wrong strategy"),
        Trap("read2_only_sample", f"{rc_sample} is the reverse complement of the amplicon; forward search maps ~0% of it.",
             "screen_design.md says only Read 2 was delivered; read names end in ' 2:N:0'; mapping rate ~0 for one library",
             "Understanding / domain knowledge gap"),
        Trap("hit_enriched_library", "83% of genes are true fitness genes. Median/total normalisation shifts every "
             "lfc_fitness by > 1 log2 unit.",
             "screen_design.md describes a validation library of carried-forward hits; the SOP states the scale "
             "(non-targeting sgRNAs at 0)", "Understanding / domain knowledge gap"),
        Trap("H_inherited_project", "The predecessor's count table, gene table and mailed list carry "
             f"{len(defects)} of 8 checklist defects ({', '.join(defects)}); two leave no trace in the scripts.",
             "SOP section 6 lists the checklist; every item is decided by recounting / re-running the inherited "
             "scripts and comparing with the inherited tables", "Verification / unverified reuse of prior work"),
        Trap("P_sample_sheet_arms", f"{sids[p['swapped'][0]]} (sheet: vehicle) received drug and "
             f"{sids[p['swapped'][1]]} (sheet: drug) received vehicle. With the sheet arms the drug contrast is "
             "diluted to one third.",
             "SOP section 3 states the criterion (a library must agree better with its own arm); library-library "
             "correlation of centred log abundance, or the reference modifiers, separate the arms 3:3; notes mention "
             "the re-racked tubes", "Verification / trusting metadata over data"),
        Trap("P_drug_vs_t0", "Fitness and drug effects are planted independently in a hit-enriched library: "
             "'depleted in the drug arm compared with T0' returns the fitness genes, not the sensitisers.",
             "instruction and SOP section 4 define the drug effect against vehicle; reference `none` genes with "
             "strong fitness effects are called by the T0 contrast", "Approach / wrong strategy"),
        Trap("Q_modifier_calls", "Drug-arm libraries come from ~25 cells per sgRNA, so the noise of depleted "
             "sgRNAs is several times that of the non-targeting sgRNAs: a null taken from the non-targeting "
             "sgRNAs calls dozens of fitness genes, the level-3 rank test both over- and under-calls.",
             "SOP section 5 states the release criterion on reference_genes.tsv (30 `none` genes over the whole "
             "fitness range, 10 modifiers); notes give the cell numbers at harvest; replicate scatter grows with "
             "depletion", "Verification / no self-check of result quality"),
    ]

    return TaskBuild(
        task_id=f"crispr-screen-s{seed:04d}-l{level}", family="crispr_screen", seed=seed, level=level,
        title="Pooled CRISPR drug-modifier screen: taking over an unfinished analysis",
        summary="Take over a pooled CRISPR drug-modifier screen from a predecessor: recount staggered amplicon "
                "reads, establish which libraries received drug, call sensitisers and resistance genes to a stated "
                "release criterion, audit the inherited analysis and write the erratum for a mailed hit list.",
        subdomain="3.2 Genomics & Sequence Analysis",
        software=["cutadapt", "MAGeCK", "Python", "numpy/scipy/pandas"],
        instruction_md=INSTRUCTION.format(n_samples=8),
        agent_dockerfile=dockerlib.agent_dockerfile(conda=["mageck=0.5.9.5", "cutadapt=5.2"]),
        verifier_dockerfile=dockerlib.verifier_dockerfile(pip=[]),
        input_files={**{f"fastq/{s}.fastq.gz": fq_dir / f"{s}.fastq.gz" for s in sids},
                     "library.tsv": library_tsv, "samples.tsv": samples_tsv, "reference_genes.tsv": ref_tsv,
                     "screen_design.md": notes, "screen_sop.md": sop_md,
                     "handover/HANDOVER.md": handover_md,
                     "handover/count_guides.py": ho / "count_guides.py",
                     "handover/analyse_screen.py": ho / "analyse_screen.py",
                     "handover/counts_draft.tsv": ho / "counts_draft.tsv",
                     "handover/gene_draft.tsv": ho / "gene_draft.tsv",
                     "handover/hits_mailed.tsv": ho / "hits_mailed.tsv"},
        solution_files={"solve.sh": ASSETS / "solve_l4.sh", "crispr_l4_pipeline.py": ASSETS / "crispr_l4_pipeline.py",
                        "solve_l4_alt.py": ASSETS / "solve_l4_alt.py"},
        tests_files={"score.py": ASSETS / "score_l4.py", "reference/truth.json": json.dumps(truth)},
        traps=traps, attacks=attacks, probes=probes,
        capabilities=["read-structure reasoning from lab notes", "per-sample QC before analysis",
                      "auditing inherited work by recomputation", "trusting data over metadata",
                      "self-checking calls against in-data controls", "choosing statistics that match the noise",
                      "exact output contracts"],
        judge_explanation_md=JUDGE, pass_threshold=0.9,
        oracle_steps=19, oracle_runtime_min=1.0, cpus=4, memory_mb=8000, storage_mb=10000,
    )
