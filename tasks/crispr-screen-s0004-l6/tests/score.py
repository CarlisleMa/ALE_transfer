"""Verifier for level 4 of the pooled CRISPR screen family (taking over a drug-modifier screen).

Everything in reference/truth.json was planted by the generator before any read existed:
  * true_counts[sgRNA][library]   reads the simulator emitted for that guide
  * genes[gene]                   latent median-guide lfc_fitness / lfc_drug, planted class, reference flag
  * arms[library]                 the treatment the library really received (the sample sheet has two exchanged)
  * audit[item]                   which checklist items were mutated in the predecessor's pipeline
  * erratum[gene]                 "confirmed" iff the mailed gene is a planted modifier
No value comes from a tool run, so tolerances only cover sampling noise and estimator choices.
"""
import math
import sys

sys.path.insert(0, "/tests")
from _scoring import CandidateError, ramp, run, to_float  # noqa: E402

CLASSES = {"sensitiser": "sensitiser", "sensitizer": "sensitiser", "resistance": "resistance", "none": "none"}


def mcc_multi(truth, pred):
    """Multi-class Matthews correlation (Gorodkin's R_K); equals the usual MCC for two classes."""
    labels = sorted(set(truth) | set(pred))
    n = len(truth)
    c = sum(1 for a, b in zip(truth, pred) if a == b)
    tk = [sum(1 for a in truth if a == k) for k in labels]
    pk = [sum(1 for b in pred if b == k) for k in labels]
    num = c * n - sum(a * b for a, b in zip(tk, pk))
    den = math.sqrt((n * n - sum(a * a for a in pk)) * (n * n - sum(a * a for a in tk)))
    return num / den if den else 0.0


def median(xs):
    xs = sorted(xs)
    return xs[len(xs) // 2] if xs else 0.0


def main(ctx):
    truth = ctx.ref_json("truth.json")
    samples = truth["samples"]
    guides = truth["guides"]            # sgRNA -> gene
    tcounts = truth["true_counts"]      # sgRNA -> [counts per library]
    genes = truth["genes"]              # gene -> {lfc_fitness, lfc_drug, class, reference}

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
    ctx.component("guide_counts", ramp(min(per_sample), 0.50, 0.90), 0.1875,
                  "min over libraries of sum(min(obs,true))/max(sum obs,sum true); full credit >= 0.90")

    # ---------------- gene_summary.tsv -------------------------------------------------
    grows = ctx.table("gene_summary.tsv", ["gene", "lfc_fitness", "lfc_drug", "drug_class"], max_rows=len(genes) + 10)
    gs = {}
    for r in grows:
        g = r["gene"]
        ctx.gate(g in genes, f"gene_summary.tsv: gene id {g[:30]!r} is in the library")
        ctx.gate(g not in gs, f"gene_summary.tsv: gene {g[:30]!r} listed once")
        cls = CLASSES.get(r["drug_class"].strip().lower())
        if cls is None:
            raise CandidateError("gene_summary.tsv: drug_class must be sensitiser, resistance or none")
        gs[g] = (to_float(r["lfc_fitness"], f"lfc_fitness of {g}"), to_float(r["lfc_drug"], f"lfc_drug of {g}"), cls)
    ctx.gate(len(gs) == len(genes), "gene_summary.tsv lists every targeted gene exactly once")

    e_fit = median([abs(gs[g][0] - genes[g]["lfc_fitness"]) for g in genes])
    e_null = median([abs(gs[g][1] - genes[g]["lfc_drug"]) for g in genes if genes[g]["class"] == "none"])
    e_mod = median([abs(gs[g][1] - genes[g]["lfc_drug"]) for g in genes if genes[g]["class"] != "none"])
    worst = max(e_fit, e_null, e_mod)
    ctx.info.update({"lfc_fitness_median_abs_error": round(e_fit, 4),
                     "lfc_drug_median_abs_error_non_modifiers": round(e_null, 4),
                     "lfc_drug_median_abs_error_modifiers": round(e_mod, 4)})
    ctx.component("gene_lfc", ramp(worst, 0.60, 0.25), 0.25,
                  "largest of three median |error| (lfc_fitness; lfc_drug of non-modifiers; lfc_drug of modifiers); "
                  "full credit <= 0.25, none >= 0.60")

    scored = [g for g in sorted(genes) if not genes[g]["reference"]]
    t_cls = [genes[g]["class"] for g in scored]
    o_cls = [gs[g][2] for g in scored]
    m = mcc_multi(t_cls, o_cls)
    ctx.info.update({
        "drug_class_mcc": round(m, 4),
        "drug_class_false_calls": sum(1 for a, b in zip(t_cls, o_cls) if a == "none" and b != "none"),
        "drug_class_missed_or_wrong_direction": sum(1 for a, b in zip(t_cls, o_cls) if a != "none" and b != a),
        "drug_class_scored_modifiers": sum(1 for a in t_cls if a != "none"),
        "reference_modifiers_right": sum(1 for g in genes if genes[g]["reference"] and genes[g]["class"] != "none"
                                         and gs[g][2] == genes[g]["class"]),
        "reference_none_called": sum(1 for g in genes if genes[g]["reference"] and genes[g]["class"] == "none"
                                     and gs[g][2] != "none")})
    ctx.component("drug_class", ramp(m, 0.75, 0.90), 0.375,
                  "three-class Matthews correlation with the planted class over the non-reference genes; "
                  "full credit >= 0.90, none <= 0.75")

    # ---------------- erratum.tsv ------------------------------------------------------
    et = truth["erratum"]
    erows = ctx.table("erratum.tsv", ["gene", "verdict"], max_rows=len(et) + 10)
    ev = {}
    for r in erows:
        g = r["gene"].strip()
        ctx.gate(g in et, f"erratum.tsv: gene {g[:30]!r} is on the mailed list")
        ctx.gate(g not in ev, f"erratum.tsv: gene {g[:30]!r} listed once")
        v = r["verdict"].strip().lower()
        if v not in ("confirmed", "retracted"):
            raise CandidateError("erratum.tsv: verdict must be confirmed or retracted")
        ev[g] = v
    ctx.gate(len(ev) == len(et), "erratum.tsv lists every mailed gene exactly once")
    eg = sorted(et)
    n_conf = sum(1 for g in eg if et[g] == "confirmed")
    n_retr = len(eg) - n_conf
    sens = sum(1 for g in eg if et[g] == "confirmed" and ev[g] == "confirmed") / max(n_conf, 1)
    spec = sum(1 for g in eg if et[g] == "retracted" and ev[g] == "retracted") / max(n_retr, 1)
    bal = (sens + spec) / 2
    ctx.info.update({"erratum_balanced_accuracy": round(bal, 4), "erratum_wrong": sum(1 for g in eg if et[g] != ev[g]),
                     "erratum_genes": len(eg), "erratum_true_confirmed": n_conf})
    # balanced accuracy, not MCC: a mailed list can hold as few as four true modifiers, where MCC
    # swings by 0.1 per verdict; "all confirmed" and "all retracted" both give exactly 0.5
    ctx.component("erratum", ramp(bal, 0.60, 0.85), 0.1875,
                  "balanced accuracy of confirmed/retracted against 'mailed gene is a planted modifier'; "
                  "full credit >= 0.85, none <= 0.60")

    # ---------------- qc.json (required, format-gated; content reported, not scored) ----
    qc = ctx.json("qc.json")
    ctx.gate(isinstance(qc, dict) and isinstance(qc.get("samples"), dict), "qc.json has a samples object")
    ok = 0
    arms_right = 0
    for j, s in enumerate(samples):
        e = qc["samples"].get(s)
        ctx.gate(isinstance(e, dict), f"qc.json has an entry for library {s}")
        arm_ok = str(e.get("arm", "")).strip().lower() == truth["arms"][s].lower()
        arms_right += int(arm_ok)
        try:
            tot = to_float(e.get("total_reads"), "total_reads")
            mapped = to_float(e.get("mapped_reads"), "mapped_reads")
        except CandidateError:
            continue
        t_tot, t_map = truth["total_reads"][j], truth["guide_reads"][j]
        if tot == t_tot and 0.90 * t_map <= mapped <= 1.005 * t_map:
            ok += 1
    # the arms and the read numbers are already scored through lfc_drug / drug_class and guide_counts
    ctx.info["qc_arms_right"] = arms_right
    ctx.info["qc_read_numbers_right"] = ok


run(main, pass_threshold=0.9)
