"""Verifier for the pooled CRISPR dropout screen family.

Truth is planted by the generator (see reference/truth.json):
  * true_counts[sgRNA][sample]  - number of reads the simulator emitted for that guide
  * gene lfc                    - latent median guide log2 fold change (0 for neutral genes)
  * fitness genes               - the genes given a non-zero effect
Nothing here is derived from a tool run, so the reference cannot be "wrong because the reference
pipeline was wrong"; tolerances only have to cover sampling noise and legitimate estimator choices.
"""
import math
import sys

sys.path.insert(0, "/tests")
from _scoring import CandidateError, ramp, run, to_float  # noqa: E402


def mcc(tp, fp, tn, fn):
    den = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return (tp * tn - fp * fn) / den if den else 0.0


def main(ctx):
    truth = ctx.ref_json("truth.json")
    samples = truth["samples"]
    guides = truth["guides"]            # sgRNA -> gene
    tcounts = truth["true_counts"]      # sgRNA -> [counts per sample]
    genes = truth["genes"]              # gene -> {"lfc": float, "fitness": 0/1}

    # ---------------- counts.tsv -------------------------------------------------------
    rows = ctx.table("counts.tsv", ["sgRNA", "gene"] + samples, max_rows=len(guides) + 10)
    seen = {}
    for r in rows:
        g = r["sgRNA"]
        ctx.gate(g in guides, f"counts.tsv: sgRNA id {g[:30]!r} is in the library")
        ctx.gate(g not in seen, f"counts.tsv: sgRNA {g[:30]!r} listed once")
        vals = []
        for s in samples:
            v = to_float(r[s], f"counts.tsv {g} {s}")
            if v < 0 or abs(v - round(v)) > 1e-6:
                raise CandidateError("counts.tsv: counts must be non-negative integers")
            vals.append(int(round(v)))
        seen[g] = vals
    ctx.gate(len(seen) == len(guides), "counts.tsv lists every library sgRNA exactly once")

    per_sample = []
    for j, s in enumerate(samples):
        inter = sum(min(seen[g][j], tcounts[g][j]) for g in guides)
        tot_c = sum(seen[g][j] for g in guides)
        tot_t = sum(tcounts[g][j] for g in guides)
        agree = inter / max(tot_c, tot_t, 1)
        per_sample.append(agree)
        ctx.info[f"count_agreement_{s}"] = round(agree, 4)
    # every sample has to be right: a rescued-or-not sample is not averaged away
    ctx.component("guide_counts", ramp(min(per_sample), 0.50, 0.90), 0.35,
                  "min over samples of sum(min(obs,true))/max(sum obs,sum true); full credit >= 0.90")

    # ---------------- gene_summary.tsv -------------------------------------------------
    grows = ctx.table("gene_summary.tsv", ["gene", "lfc", "is_hit"], max_rows=len(genes) + 10)
    gs = {}
    for r in grows:
        g = r["gene"]
        ctx.gate(g in genes, f"gene_summary.tsv: gene id {g[:30]!r} is in the library")
        ctx.gate(g not in gs, f"gene_summary.tsv: gene {g[:30]!r} listed once")
        hit = r["is_hit"].strip()
        if hit not in ("0", "1"):
            raise CandidateError("gene_summary.tsv: is_hit must be 0 or 1")
        gs[g] = (to_float(r["lfc"], f"lfc of {g}"), int(hit))
    ctx.gate(len(gs) == len(genes), "gene_summary.tsv lists every targeted gene exactly once")

    errs = sorted(abs(gs[g][0] - genes[g]["lfc"]) for g in genes)
    med_err = errs[len(errs) // 2]
    ctx.info["lfc_median_abs_error"] = round(med_err, 4)
    ctx.component("gene_lfc", ramp(med_err, 0.60, 0.20), 0.30,
                  "median |lfc - planted lfc| over genes; full credit <= 0.20, none >= 0.60")

    tp = sum(1 for g in genes if gs[g][1] == 1 and genes[g]["fitness"] == 1)
    fp = sum(1 for g in genes if gs[g][1] == 1 and genes[g]["fitness"] == 0)
    fn = sum(1 for g in genes if gs[g][1] == 0 and genes[g]["fitness"] == 1)
    tn = sum(1 for g in genes if gs[g][1] == 0 and genes[g]["fitness"] == 0)
    m = mcc(tp, fp, tn, fn)
    ctx.info.update({"hits_tp": tp, "hits_fp": fp, "hits_fn": fn, "hits_tn": tn, "hits_mcc": round(m, 4)})
    # MCC, not F1: in a hit-enriched library "call everything a hit" has F1 > 0.9 but MCC = 0
    ctx.component("hit_calls", ramp(m, 0.50, 0.90), 0.30,
                  "Matthews correlation of is_hit vs planted fitness genes; full credit >= 0.90")

    # ---------------- qc.json ----------------------------------------------------------
    qc = ctx.json("qc.json")
    ctx.gate(isinstance(qc, dict) and isinstance(qc.get("samples"), dict), "qc.json has a samples object")
    ok = 0
    for j, s in enumerate(samples):
        e = qc["samples"].get(s)
        if not isinstance(e, dict):
            continue
        try:
            tot = to_float(e.get("total_reads"), "total_reads")
            mapped = to_float(e.get("mapped_reads"), "mapped_reads")
        except CandidateError:
            continue
        t_tot, t_map = truth["total_reads"][j], truth["guide_reads"][j]
        if tot == t_tot and 0.90 * t_map <= mapped <= 1.005 * t_map:
            ok += 1
    ctx.component("qc_report", ok / len(samples), 0.05,
                  "per sample: total_reads exact, mapped_reads within [0.90, 1.005] x reads carrying a library guide")


run(main, pass_threshold=0.9)
