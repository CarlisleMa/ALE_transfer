"""Verifier for the XL-MS model selection family.

Truth is planted by the generator (reference/truth.json):
  * every model is a rigid-body variant of a hidden native, so `near_native` (interface RMSD < 2 A)
    and the interface RMSD / DockQ are exact properties of the transform that produced the model;
  * the satisfied / evaluable counts were computed with the planted residue mapping on the
    coordinates as written to the model files, and no evaluated distance is within 0.5 A of the
    cutoff, so the counts do not depend on how a correct implementation reads the files;
  * the mapping of cross-linked positions onto model_01 is the planted construct layout.
Nothing is derived from the reference solution.
"""
import sys

sys.path.insert(0, "/tests")
from _scoring import CandidateError, ramp, run, to_float  # noqa: E402

COUNT_COLS = ["inter_satisfied", "inter_total", "intra_satisfied", "intra_total"]
NA = {"", "na", "n/a", "nan", "none", "null", "-", "."}


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

    # ---------------- restraints.tsv ---------------------------------------------------
    rows = ctx.table("restraints.tsv", ["model"] + COUNT_COLS, max_rows=len(names) + 20)
    got = {}
    for r in rows:
        m = model_name(r["model"])
        ctx.gate(m in models, f"restraints.tsv: {m[:30]!r} is one of the models")
        ctx.gate(m not in got, f"restraints.tsv: {m[:30]!r} listed once")
        got[m] = [as_int(r[c], f"restraints.tsv {m} {c}") for c in COUNT_COLS]
    ctx.gate(len(got) == len(names), "restraints.tsv lists every model exactly once")

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

    ctx.info["best_model"] = best
    ctx.info["best_model_irmsd"] = models[best]["irmsd"]
    ctx.info["best_model_dockq"] = models[best]["dockq"]
    ctx.component("best_model", 1.0 if models[best]["near_native"] else 0.0, 0.35,
                  "best_model is near-native (hidden interface RMSD < 2 A)")

    pos = {m: i for i, m in enumerate(ranking)}
    near = [n for n in names if models[n]["near_native"]]
    wrong = [n for n in names if not models[n]["near_native"]]
    auc = sum(1 for a in near for b in wrong if pos[a] < pos[b]) / (len(near) * len(wrong))
    ctx.info["ranking_auc"] = round(auc, 4)
    ctx.component("ranking", ramp(auc, 0.70, 0.95), 0.20,
                  "ROC AUC of the submitted order, near-native versus wrong models; full credit >= 0.95, none <= 0.70")

    ok = 0
    for t in truth["mapping"]:
        c = mp.get((t["accession"], t["position"]))
        if c is None:
            continue
        chain, num = c[0].strip(), c[1].strip()
        if t["chain"] is None:
            ok += int(chain.lower() in NA and num.lower() in NA)
        else:
            try:
                ok += int(chain == t["chain"] and float(num) == float(t["resnum"]))
            except ValueError:
                pass
    mfrac = ok / len(truth["mapping"])
    ctx.info["mapping_correct"] = f"{ok}/{len(truth['mapping'])}"
    ctx.component("residue_mapping", 1.0 if ok == len(truth["mapping"]) else 0.6 * mfrac, 0.15,
                  "fraction of cross-linked positions with the exact chain and residue number in model_01 "
                  "(NA when the residue is not in the model); 1.0 only if every position is right, "
                  "otherwise 0.6 x fraction")


run(main, pass_threshold=0.9)
