"""Verifier for level 4 of the XL-MS model selection family.

Truth (reference/truth.json):
  * every model is a rigid-body variant of a hidden native (heterodimer 2 moved relative to
    heterodimer 1), so `near_native` and the interface RMSD are exact properties of the transform;
  * the status of every cross-link in every model follows the definition written in the task's
    xlms_notes.md (surface distance on a lattice, any assignment of the two positions to the copies).
    Surface distances are computed, not planted: the generator uses only cross-links whose surface
    distance is at least 3 A away from the cutoff in every model, so that a correct implementation
    cannot land on the other side, and the release gate requires two independent implementations
    (scipy Dijkstra, numpy bucket search) to reproduce every status;
  * the mapping of cross-linked positions onto model_01 is the planted construct layout.
"""
import math
import re
import sys

sys.path.insert(0, "/tests")
from _scoring import CandidateError, ramp, run, to_float  # noqa: E402

COUNT_COLS = ["inter_satisfied", "inter_total", "intra_satisfied", "intra_total"]
NA = {"", "na", "n/a", "nan", "none", "null", "-", "."}
STATUS = {"satisfied": 0, "violated": 1, "not_evaluable": 2}


def model_name(s):
    s = str(s).strip().replace("\\", "/").split("/")[-1]
    return s[:-4] if s.lower().endswith(".pdb") else s


def as_int(s, what):
    v = to_float(s, what)
    if abs(v - round(v)) > 1e-9:
        raise CandidateError(f"{what} must be an integer")
    return int(round(v))


def main(ctx):
    truth = ctx.ref_json("truth.json")
    models = truth["models"]
    names = list(models)
    status = truth["status"]
    xl_ids = [x["xl_id"] for x in truth["crosslinks"]]

    # ---------------- restraints.tsv ---------------------------------------------------
    rows = ctx.table("restraints.tsv", ["model"] + COUNT_COLS, max_rows=len(names) + 20)
    got = {}
    for r in rows:
        m = model_name(r["model"])
        ctx.gate(m in models, f"restraints.tsv: {m[:30]!r} is one of the models")
        ctx.gate(m not in got, f"restraints.tsv: {m[:30]!r} listed once")
        got[m] = [as_int(r[c], f"restraints.tsv {m} {c}") for c in COUNT_COLS]
    ctx.gate(len(got) == len(names), "restraints.tsv lists every model exactly once")

    # ---------------- link_status.tsv --------------------------------------------------
    lrows = ctx.table("link_status.tsv", ["model", "xl_id", "status"], max_rows=len(names) * len(xl_ids) + 50)
    calls = {}
    for r in lrows:
        m, x = model_name(r["model"]), (r["xl_id"] or "").strip()
        st = (r["status"] or "").strip().lower().replace(" ", "_").replace("-", "_")
        ctx.gate(m in models and x in status[names[0]], f"link_status.tsv: {m[:30]!r} / {x[:20]!r} is a model and a cross-link")
        ctx.gate((m, x) not in calls, f"link_status.tsv: {m[:30]} / {x[:20]} listed once")
        ctx.gate(st in STATUS, f"link_status.tsv: status {st[:20]!r} is satisfied, violated or not_evaluable")
        calls[(m, x)] = st
    ctx.gate(len(calls) == len(names) * len(xl_ids), "link_status.tsv has one row for every model and cross-link")

    # ---------------- selection.json ---------------------------------------------------
    sel = ctx.json("selection.json")
    ctx.gate(isinstance(sel, dict) and isinstance(sel.get("best_model"), str)
             and isinstance(sel.get("ranking"), list), "selection.json has best_model (string) and ranking (list)")
    ctx.gate(len(sel["ranking"]) <= len(names) + 5 and all(isinstance(x, str) for x in sel["ranking"]),
             "selection.json: ranking is a list of model names")
    best = model_name(sel["best_model"])
    ranking = [model_name(x) for x in sel["ranking"]]
    ctx.gate(best in models, "selection.json: best_model is one of the models")
    ctx.gate(sorted(ranking) == sorted(names), "selection.json: ranking lists every model exactly once")
    ctx.gate(ranking[0] == best, "selection.json: the first element of ranking is best_model")

    # ---------------- mapping.tsv ------------------------------------------------------
    mrows = ctx.table("mapping.tsv", ["accession", "position", "model_chain", "model_resnum"], max_rows=5000)
    mp = {}
    for r in mrows:
        key = (r["accession"], as_int(r["position"], "mapping.tsv position"))
        ctx.gate(key not in mp, f"mapping.tsv: position {key[0][:20]}:{key[1]} listed once")
        mp[key] = (r["model_chain"] or "", r["model_resnum"] or "")

    # ---------------- components -------------------------------------------------------
    exact = sum(1 for n in names if got[n] == [models[n][c] for c in COUNT_COLS])
    frac = exact / len(names)
    ctx.info["models_with_exact_counts"] = f"{exact}/{len(names)}"
    ctx.component("restraint_counts", 1.0 if exact == len(names) else 0.6 * frac, 0.30,
                  "fraction of models whose four counts are all exact; 1.0 only if every model is exact, "
                  "otherwise 0.6 x fraction")

    # three-class Matthews correlation (Gorodkin's R_K) over all model x cross-link statuses
    conf = [[0] * 3 for _ in range(3)]
    for n in names:
        for x in xl_ids:
            conf[STATUS[status[n][x]]][STATUS[calls[(n, x)]]] += 1
    s = sum(map(sum, conf))
    c = sum(conf[k][k] for k in range(3))
    t = [sum(conf[k]) for k in range(3)]
    p = [sum(conf[j][k] for j in range(3)) for k in range(3)]
    den = math.sqrt((s * s - sum(v * v for v in p)) * (s * s - sum(v * v for v in t)))
    mcc = (c * s - sum(a * b for a, b in zip(p, t))) / den if den > 0 else 0.0
    ctx.info["link_status_wrong"] = s - c
    ctx.info["link_status_mcc"] = round(mcc, 4)
    ctx.component("link_status", ramp(mcc, 0.95, 0.995), 0.15,
                  "three-class Matthews correlation of the per-model per-cross-link statuses; full credit "
                  ">= 0.995, none <= 0.95")

    ctx.info["best_model"] = best
    ctx.info["best_model_irmsd"] = models[best]["irmsd"]
    ctx.info["best_model_dockq"] = models[best]["dockq"]
    ctx.component("best_model", 1.0 if models[best]["near_native"] else 0.0, 0.30,
                  "best_model is near-native (hidden interface RMSD < 2 A)")

    pos = {m: i for i, m in enumerate(ranking)}
    near = [n for n in names if models[n]["near_native"]]
    wrong = [n for n in names if not models[n]["near_native"]]
    auc = sum(1 for a in near for b in wrong if pos[a] < pos[b]) / (len(near) * len(wrong))
    ctx.info["ranking_auc"] = round(auc, 4)
    ctx.component("ranking", ramp(auc, 0.70, 0.95), 0.10,
                  "ROC AUC of the submitted order, near-native versus wrong models; full credit >= 0.95, none <= 0.70")

    ok = 0
    for tr in truth["mapping"]:
        cand = mp.get((tr["accession"], tr["position"]))
        if cand is None:
            continue
        chain, num = cand[0].strip(), cand[1].strip()
        if tr["chains"] is None:
            ok += int(chain.lower() in NA and num.lower() in NA)
        else:
            try:
                given = sorted(x for x in re.split(r"[,;/|\s]+", chain) if x)
                ok += int(given == sorted(tr["chains"]) and float(num) == float(tr["resnum"]))
            except ValueError:
                pass
    mfrac = ok / len(truth["mapping"])
    ctx.info["mapping_correct"] = f"{ok}/{len(truth['mapping'])}"
    ctx.component("residue_mapping", 1.0 if ok == len(truth["mapping"]) else 0.6 * mfrac, 0.15,
                  "fraction of cross-linked positions with the exact set of chains and residue number in model_01 "
                  "(NA when the residue is not in the model); 1.0 only if every position is right, "
                  "otherwise 0.6 x fraction")


run(main, pass_threshold=0.9)
