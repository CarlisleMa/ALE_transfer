#!/usr/bin/env python3
"""Reference solution, level 4 (taking over a drug-modifier CRISPR screen).

  1. counts          cutadapt (linked adapter, both strands) + mageck count, current library release
  2. sample identity arms of the day-14 libraries from the data (first principal component of the
                     centred log abundances), named by the majority of the sample-sheet labels and
                     cross-checked with the reference modifiers
  3. statistics      size factors from the non-targeting sgRNAs; lfc_fitness = vehicle vs T0,
                     lfc_drug = drug vs vehicle; calls from an abundance-aware z test (replicate
                     variance smoothed over abundance, Stouffer over the gene's sgRNAs, empirical-null
                     rescaling, BH 5%); self-check against the reference genes (SOP section 5)
  4. audit           inherited count table against a re-run of the inherited count script and against
                     the own counts; inherited analysis script re-run on the own counts with the
                     corrected arms and compared with the own estimates
  5. erratum         every mailed gene against the final classes

`--mut` and `--l3` replay careless analyses for the validator's attack suite (one decision ignored
at a time; the audit item of the ignored decision is then reported `ok`), and `--mut call_mageck`
swaps the calling step for MAGeCK RRA (a must-pass probe). The oracle uses neither.
"""
import argparse
import gzip
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

PREFIX_TAIL = "GGACGAAACACCG"            # 3' end of the U6 promoter (+G) directly upstream of the spacer
SCAFFOLD_HEAD = "GTTTTAGAGCTAGAAATAGC"   # 5' end of the sgRNA scaffold directly downstream
ITEMS = ["library_release", "sample_columns", "stagger_reads", "read2_library",
         "arm_labels", "control_normalisation", "drug_baseline", "contrast_sign"]
MUT_ITEM = {"library_r1": "library_release", "file_order": "sample_columns", "fixed_offset": "stagger_reads",
            "no_revcomp": "read2_library", "sheet_arms": "arm_labels", "median_norm": "control_normalisation",
            "baseline_t0": "drug_baseline", "reversed": "contrast_sign"}

ap = argparse.ArgumentParser()
ap.add_argument("--input", default="/task/input")
ap.add_argument("--output", default="/task/output")
ap.add_argument("--work", default="/tmp/crispr_l4_work")
ap.add_argument("--mut", action="append", default=[],
                choices=list(MUT_ITEM) + ["inherited_counts", "call_ntc_z", "call_mwu", "call_mageck", "no_audit"])
ap.add_argument("--l3", action="store_true", help="level-3 reference method, nothing else")
a = ap.parse_args()
mut = set(a.mut)

inp, out, work = Path(a.input), Path(a.output), Path(a.work)
ho = inp / "handover"
out.mkdir(parents=True, exist_ok=True)
work.mkdir(parents=True, exist_ok=True)

lib = pd.read_csv(inp / "library.tsv", sep="\t")
samples = pd.read_csv(inp / "samples.tsv", sep="\t")
sids = list(samples["sample_id"])
fq = {r.sample_id: inp / "fastq" / r.fastq_file for r in samples.itertuples()}
sheet_arm = dict(zip(samples["sample_id"], samples["arm"]))
run_of = dict(zip(samples["sample_id"], samples["sequencing_run"]))
ref = pd.read_csv(inp / "reference_genes.tsv", sep="\t")
COMP = str.maketrans("ACGTN", "TGCAN")


def sh(cmd):
    r = subprocess.run([str(c) for c in cmd], capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"command failed: {' '.join(map(str, cmd))}\n{r.stderr[-2000:]}")
    return r


def n_records(path):
    n = 0
    with gzip.open(path, "rt") as fh:
        for _ in fh:
            n += 1
    return n // 4


total_reads = {s: n_records(fq[s]) for s in sids}

# ============================================================================== 1. guide counts
lib[["sgRNA", "sequence", "gene"]].to_csv(work / "library.csv", index=False, header=False)


def count_cutadapt(revcomp=True):
    """per-file counts (dict sample -> Series indexed by sgRNA); a file MAGeCK cannot use gives zeros"""
    res = {}
    for s in sids:
        o = work / f"{s}.{'rc' if revcomp else 'fw'}.spacer.fastq.gz"
        cmd = ["cutadapt", "-j", "2", "-e", "0.1", "-g", f"{PREFIX_TAIL}...{SCAFFOLD_HEAD}",
               "--discard-untrimmed", "-m", "20", "-M", "20", "-o", o, fq[s]]
        if revcomp:
            cmd.insert(1, "--revcomp")
        sh(cmd)
        pre = work / f"cnt_{s}_{'rc' if revcomp else 'fw'}"
        r = subprocess.run(["mageck", "count", "-l", str(work / "library.csv"), "-n", str(pre), "--sample-label", s,
                            "--fastq", str(o), "--trim-5", "0", "--norm-method", "none"], capture_output=True, text=True)
        f = Path(str(pre) + ".count.txt")
        if r.returncode == 0 and f.exists():
            t = pd.read_csv(f, sep="\t")
            res[s] = lib["sgRNA"].map(dict(zip(t["sgRNA"], t[s]))).fillna(0).astype(int).values
        else:
            res[s] = np.zeros(len(lib), dtype=int)
    return pd.DataFrame(res, index=lib["sgRNA"])


def count_fixed():
    seq2k = {s: k for k, s in enumerate(lib["sequence"])}
    res = {}
    for s in sids:
        v = np.zeros(len(lib), dtype=int)
        with gzip.open(fq[s], "rt") as fh:
            for i, line in enumerate(fh):
                if i % 4 == 1:
                    k = seq2k.get(line[24:44])
                    if k is None:
                        k = seq2k.get(line.rstrip().translate(COMP)[::-1][15:35])
                    if k is not None:
                        v[k] += 1
        res[s] = v
    return pd.DataFrame(res, index=lib["sgRNA"])


clean = count_cutadapt(revcomp=True)                      # counts[sgRNA, sample], the deliverable
gene_clean = pd.Series(lib["gene"].values, index=lib["sgRNA"])
draft = pd.read_csv(ho / "counts_draft.tsv", sep="\t").set_index("sgRNA").reindex(lib["sgRNA"])

final, gene_map = clean, gene_clean
if "no_revcomp" in mut:
    final = count_cutadapt(revcomp=False)
if "fixed_offset" in mut:
    final = count_fixed()
if "file_order" in mut:                                   # files taken in lexicographic order
    order = sorted(sids, key=lambda s: str(fq[s]))
    final = pd.DataFrame({s: final[order[j]].values for j, s in enumerate(sids)}, index=final.index)
if "inherited_counts" in mut:
    final = draft[sids].fillna(0).astype(int)
    gene_map = draft["gene"]
if "library_r1" in mut:
    r1 = pd.read_csv(Path(__file__).parent / "library_r1.tsv", sep="\t")
    gene_map = pd.Series(r1["gene"].values, index=r1["sgRNA"]).reindex(lib["sgRNA"])


# ============================================================================== 2. sample identity
def size_factors(C, ntc, median_norm=False):
    if median_norm:                                       # median of ratios over the whole library
        logc = np.log(C.where(C > 0))
        return np.exp(logc.sub(logc.mean(axis=1), axis=0).median(axis=0))
    tot = C[ntc].sum(axis=0)
    return tot / tot.mean()


def arms_from_data(C, ntc):
    """T0 as on the sheet; the day-14 libraries split by the sign of their first principal component"""
    alive = [s for s in sids if C[s].sum() > 0]
    day14 = [s for s in alive if sheet_arm[s] != "T0"]
    N = C[day14] / size_factors(C[day14], ntc)
    X = np.log2(N + 8.0)
    X = X.sub(X.mean(axis=1), axis=0)
    _, _, vt = np.linalg.svd(X.values, full_matrices=False)
    side = dict(zip(day14, vt[0] > 0))
    # name the two groups: majority of the sheet labels ...
    votes = sum(1 if (sheet_arm[s] == "drug") == side[s] else -1 for s in day14)
    drug_side = votes >= 0
    arms = {s: ("T0" if sheet_arm[s] == "T0" else ("drug" if side.get(s, False) == drug_side else "vehicle"))
            for s in sids}
    return arms, X


def reference_signature(X, gene_of):
    """reference sensitisers minus reference resistance genes, per day-14 library (negative = drug)"""
    exp = dict(zip(ref["gene"], ref["expected_class"]))
    g = gene_of.reindex(X.index).map(exp)
    return (X[(g == "sensitiser").values].mean(axis=0) - X[(g == "resistance").values].mean(axis=0)).round(2).to_dict()


# ============================================================================== 3. statistics
def bh(p):
    p = np.asarray(p, dtype=float)
    o = np.argsort(p)
    adj = p[o] * len(p) / (np.arange(len(p)) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    res = np.empty(len(p))
    res[o] = np.minimum(adj, 1.0)
    return res


def analyse(C, gene_of, arms, median_norm=False, baseline="vehicle", reverse=False, caller="meanvar"):
    alive = [s for s in sids if C[s].sum() > 0]
    C = C[alive].astype(float)
    ntc = (gene_of.reindex(C.index) == "NON_TARGETING").values
    grp = {k: [s for s in alive if arms[s] == k] for k in ("T0", "vehicle", "drug")}
    N = C / size_factors(C, ntc, median_norm)

    def lfc(x, y):
        v = np.log2((N[x].mean(axis=1) + 0.5) / (N[y].mean(axis=1) + 0.5))
        return v if median_norm else v - np.median(v[ntc])

    fit = lfc(grp["vehicle"], grp["T0"])
    A, B = grp["drug"], grp[baseline]
    if reverse:
        A, B = B, A
    drg = lfc(A, B)
    genes = gene_of.reindex(C.index).values
    tab = pd.DataFrame({"gene": genes, "fit": fit.values, "drg": drg.values})[~ntc]
    summ = tab.groupby("gene").agg(lfc_fitness=("fit", "median"), lfc_drug=("drg", "median"))

    if caller == "meanvar":
        L = np.log2(N + 4.0)
        d = L[A].mean(axis=1) - L[B].mean(axis=1)
        d = d - np.median(d[ntc])
        ss = ((L[A].sub(L[A].mean(axis=1), axis=0)) ** 2).sum(axis=1) + ((L[B].sub(L[B].mean(axis=1), axis=0)) ** 2).sum(axis=1)
        s2 = ss / max(len(A) + len(B) - 2, 1)
        o = np.argsort(L[A + B].mean(axis=1).values)          # variance as a function of abundance
        sm = pd.Series(s2.values[o]).rolling(101, center=True, min_periods=25).mean().values
        s2s = np.empty(len(o))
        s2s[o] = sm
        z = d.values / np.sqrt(s2s * (1 / len(A) + 1 / len(B)))
        zt = pd.DataFrame({"gene": genes, "z": z})[~ntc]
        gz = zt.groupby("gene")["z"].agg(lambda v: v.sum() / np.sqrt(len(v)))
        mad = stats.median_abs_deviation(gz.values, scale="normal")
        gz = (gz - gz.median()) / max(mad, 1.0)               # empirical null: most genes are no modifiers
        p = 2 * stats.norm.sf(np.abs(gz.values))
        sign = np.sign(gz.values)
    elif caller == "ntc_z":                                   # gene median against medians of 4 random controls
        rng = np.random.default_rng(0)
        null = np.median(rng.choice(drg.values[ntc], size=(20000, 4), replace=True), axis=1)
        gz = summ["lfc_drug"] / null.std()
        p = 2 * stats.norm.sf(np.abs(gz.values))
        sign = np.sign(gz.values)
    else:                                                     # level-3 reference test, two-sided
        null = drg.values[ntc]
        res = tab.groupby("gene")["drg"].apply(lambda v: stats.mannwhitneyu(v.values, null, alternative="two-sided").pvalue)
        p = res.reindex(summ.index).values
        sign = np.sign(summ["lfc_drug"].values)
    summ["fdr"] = bh(p)
    summ["drug_class"] = np.where(summ["fdr"] < 0.05, np.where(sign < 0, "sensitiser", "resistance"), "none")
    return summ


def l3_hits(C, gene_of, treat, ctrl):
    """the level-3 reference statistics for one 'timepoint versus T0' contrast"""
    C = C.astype(float)
    ntc = (gene_of.reindex(C.index) == "NON_TARGETING").values
    N = C / size_factors(C, ntc)
    g = np.log2((N[treat].mean(axis=1) + 0.5) / (N[ctrl].mean(axis=1) + 0.5))
    g = g - np.median(g[ntc])
    tab = pd.DataFrame({"gene": gene_of.reindex(C.index).values, "lfc": g.values})[~ntc]
    null = g.values[ntc]
    med = tab.groupby("gene")["lfc"].median()
    pv = tab.groupby("gene")["lfc"].apply(lambda v: stats.mannwhitneyu(v.values, null, alternative="less").pvalue)
    return med, pd.Series(bh(pv.values) < 0.05, index=pv.index)


ntc_clean = (gene_clean == "NON_TARGETING").values
arms_clean, Xc = arms_from_data(clean, ntc_clean)
own = analyse(clean, gene_clean, arms_clean)
print("arms from the data:", arms_clean, file=sys.stderr)
print("reference signature (sensitisers - resistance, per day-14 library):", reference_signature(Xc, gene_clean), file=sys.stderr)

if a.l3:
    arms = dict(sheet_arm)
    t0 = [s for s in sids if arms[s] == "T0"]
    lf, hit_v = l3_hits(clean, gene_clean, [s for s in sids if arms[s] == "vehicle"], t0)
    ld, hit_d = l3_hits(clean, gene_clean, [s for s in sids if arms[s] == "drug"], t0)
    summ = pd.DataFrame({"lfc_fitness": lf, "lfc_drug": ld})
    summ["drug_class"] = np.where(hit_d & ~hit_v, "sensitiser", "none")
else:
    ntc_f = (gene_map.reindex(final.index) == "NON_TARGETING").values
    if "sheet_arms" in mut:
        arms = dict(sheet_arm)
    else:
        arms, _ = arms_from_data(final, ntc_f)
    summ = analyse(final, gene_map, arms, median_norm="median_norm" in mut,
                   baseline="T0" if "baseline_t0" in mut else "vehicle", reverse="reversed" in mut,
                   caller="ntc_z" if "call_ntc_z" in mut else ("mwu" if "call_mwu" in mut else "meanvar"))
if "call_mageck" in mut:
    # third way to the calls (a must-pass probe of the validator): MAGeCK RRA on the same contrast
    tab = final.copy()
    tab.insert(0, "gene", gene_map.reindex(final.index).values)
    tab.to_csv(work / "for_rra.txt", sep="\t", index_label="sgRNA")
    (work / "ntc.txt").write_text("\n".join(lib.loc[lib["gene"] == "NON_TARGETING", "sgRNA"]) + "\n")
    sh(["mageck", "test", "-k", work / "for_rra.txt", "-t", ",".join(s for s in sids if arms[s] == "drug"),
        "-c", ",".join(s for s in sids if arms[s] == "vehicle"), "--norm-method", "control",
        "--control-sgrna", work / "ntc.txt", "-n", work / "rra"])
    rra = pd.read_csv(work / "rra.gene_summary.txt", sep="\t").set_index("id").reindex(summ.index)
    summ["fdr"] = np.minimum(rra["neg|fdr"], rra["pos|fdr"]).values
    summ["drug_class"] = np.where(rra["neg|fdr"].values < 0.05, "sensitiser",
                                  np.where(rra["pos|fdr"].values < 0.05, "resistance", "none"))
all_genes = sorted(set(lib.loc[lib["gene"] != "NON_TARGETING", "gene"]))
summ = summ.reindex(all_genes)
summ["drug_class"] = summ["drug_class"].fillna("none")
summ = summ.fillna(0.0)

# self-check against the reference genes (SOP section 5)
exp = dict(zip(ref["gene"], ref["expected_class"]))
mods = [g for g in exp if exp[g] != "none"]
inert = [g for g in exp if exp[g] == "none"]
n_mod_ok = sum(summ.loc[g, "drug_class"] == exp[g] for g in mods)
n_inert_called = sum(summ.loc[g, "drug_class"] != "none" for g in inert)
print(f"release criterion: {n_mod_ok}/{len(mods)} reference modifiers right, "
      f"{n_inert_called}/{len(inert)} reference none genes called", file=sys.stderr)

# ============================================================================== 4. audit
audit = {}
audit["library_release"] = "ok" if (draft["gene"].values == lib["gene"].values).all() else "defective"

# re-run the inherited count script the way its docstring says (sheet order, current library)
rerun_f = work / "rerun_counts.tsv"
sh([sys.executable, ho / "count_guides.py", inp / "library.tsv", inp / "samples.tsv", rerun_f] + [fq[s] for s in sids])
rerun = pd.read_csv(rerun_f, sep="\t").set_index("sgRNA").reindex(lib["sgRNA"])
audit["sample_columns"] = "ok" if all((draft[s].values == rerun[s].values).all() for s in sids) else "defective"
run_a = [s for s in sids if run_of[s] == "A"]
ratio = {s: rerun[s].sum() / max(clean[s].sum(), 1) for s in run_a}
audit["stagger_reads"] = "ok" if min(ratio.values()) > 0.5 else "defective"
col_tot = draft[sids].sum(axis=0)
audit["read2_library"] = "ok" if col_tot.min() > 0.01 * col_tot.median() else "defective"
audit["arm_labels"] = "ok" if arms_clean == sheet_arm else "defective"

# re-run the inherited analysis script on the own counts with the corrected arms
fixed_sheet = samples.copy()
fixed_sheet["arm"] = [arms_clean[s] for s in sids]
fixed_sheet.to_csv(work / "samples_fixed.tsv", sep="\t", index=False)
cc = clean.copy()
cc.insert(0, "gene", gene_clean.values)
cc.to_csv(work / "counts_clean.tsv", sep="\t", index_label="sgRNA")
sh([sys.executable, ho / "analyse_screen.py", work / "counts_clean.tsv", work / "samples_fixed.tsv",
    work / "their_genes.tsv", work / "their_hits.tsv"])
their = pd.read_csv(work / "their_genes.tsv", sep="\t").set_index("gene").reindex(own.index)
shift = float((their["lfc_fitness"] - own["lfc_fitness"]).median())
audit["control_normalisation"] = "ok" if abs(shift) < 0.3 else "defective"
vs_t0 = analyse(clean, gene_clean, arms_clean, baseline="T0")["lfc_drug"]
cands = {("vehicle", 1): own["lfc_drug"], ("T0", 1): vs_t0, ("vehicle", -1): -own["lfc_drug"], ("T0", -1): -vs_t0}


def spread(x):                                               # shift-invariant distance
    return float((x - x.median()).abs().median())


best = min(cands, key=lambda k: spread(their["lfc_drug"] - cands[k]))
audit["drug_baseline"] = "ok" if best[0] == "vehicle" else "defective"
audit["contrast_sign"] = "ok" if best[1] == 1 else "defective"
print("audit:", audit, "| normalisation shift", round(shift, 2), "| inherited drug effect looks like", best, file=sys.stderr)

if a.l3 or "inherited_counts" in mut or "no_audit" in mut:
    audit = {it: "ok" for it in ITEMS}
for m in mut:
    if m in MUT_ITEM:
        audit[MUT_ITEM[m]] = "ok"

# ============================================================================== 5. deliverables
cnt = final.copy()
cnt.insert(0, "gene", lib["gene"].values)
cnt.to_csv(out / "counts.tsv", sep="\t", index_label="sgRNA")
summ[["lfc_fitness", "lfc_drug", "drug_class"] + [c for c in summ.columns if c == "fdr"]].to_csv(
    out / "gene_summary.tsv", sep="\t", index_label="gene")
pd.DataFrame({"item": ITEMS, "verdict": [audit[it] for it in ITEMS]}).to_csv(out / "audit.tsv", sep="\t", index=False)
mailed = pd.read_csv(ho / "hits_mailed.tsv", sep="\t")["gene"]
verdict = ["confirmed" if (a.l3 or summ.loc[g, "drug_class"] != "none") else "retracted" for g in mailed]
pd.DataFrame({"gene": mailed, "verdict": verdict}).to_csv(out / "erratum.tsv", sep="\t", index=False)
qc = {"samples": {s: {"total_reads": int(total_reads[s]), "mapped_reads": int(final[s].sum()), "arm": arms[s]} for s in sids}}
(out / "qc.json").write_text(json.dumps(qc, indent=1))
print(json.dumps({"calls": summ["drug_class"].value_counts().to_dict(),
                  "mapping_rate": {s: round(final[s].sum() / total_reads[s], 3) for s in sids}}))
