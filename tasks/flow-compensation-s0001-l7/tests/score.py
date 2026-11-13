"""Verifier of level 4 of the flow-cytometry compensation / gating family.

Truth is planted by the generator (reference/truth.json):
  * spillover[sample]  - the matrix the simulator mixed the fluorochromes of that sample with (sample
                         reagent lots, detector settings of the sample's acquisition block)
  * lineage populations - per sample, number of events counted from the hidden class label
  * FMO-gated populations - the SOP's rule applied to planted quantities (signals compensated with the
                         planted matrix, parent = hidden label, threshold = 99th percentile of the
                         label-defined parent in the FMO tube of the same block)
Nothing is derived from a compensation or gating run of the reference pipeline.
"""
import csv
import io
import re
import sys

sys.path.insert(0, "/tests")
from _scoring import CandidateError, ramp, run, to_float  # noqa: E402


def norm(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def read_matrices(ctx, samples, fluors, dets):
    """spillover.csv, long layout: sample, fluorochrome, <8 detectors> -> {sample: 8x8 list}."""
    txt = ctx.text("spillover.csv", max_bytes=2_000_000)
    try:
        rows = [r for r in csv.reader(io.StringIO(txt)) if any(c.strip() for c in r)]
    except csv.Error as e:
        raise CandidateError(f"spillover.csv is not valid CSV: {e}")
    ctx.gate(2 <= len(rows) <= 2000, "spillover.csv has a header and at most 1999 rows")
    head = [norm(c) for c in rows[0]]
    need = ["sample", "fluorochrome"] + [norm(d) for d in dets]
    if any(head.count(c) != 1 for c in need):
        raise CandidateError("gate failed: spillover.csv has the columns sample, fluorochrome and the eight "
                             "detectors of panel.tsv (each once)")
    col = {c: head.index(c) for c in need}
    fl = {norm(f): i for i, f in enumerate(fluors)}
    sid = {s: s for s in samples}
    got = {}
    for r in rows[1:]:
        if len(r) < len(head):
            raise CandidateError("spillover.csv: a row is shorter than the header")
        s, f = r[col["sample"]].strip(), norm(r[col["fluorochrome"]])
        if s not in sid or f not in fl:
            continue                                   # rows for other tubes are ignored
        ctx.gate((s, f) not in got, f"spillover.csv: {s}/{fluors[fl[f]]} listed once")
        got[(s, f)] = [to_float(r[col[norm(d)]], f"spillover.csv entry {s}/{f}/{d}") for d in dets]
    ctx.gate(len(got) == len(samples) * len(fluors),
             "spillover.csv holds one row for every sample x fluorochrome")
    return {s: [got[(s, norm(f))] for f in fluors] for s in samples}


def main(ctx):
    truth = ctx.ref_json("truth.json")
    fluors, dets = truth["fluorochromes"], truth["detectors"]
    samples = truth["samples"]
    pops = truth["populations"]
    tol = truth["tolerance"]
    lineage, fmo_gated = truth["lineage"], truth["fmo_gated"]

    # ---------------- spillover.csv ----------------------------------------------------
    M = read_matrices(ctx, samples, fluors, dets)
    n = len(fluors)
    max_err, where = 0.0, ""
    for s in samples:
        S = truth["spillover"][s]
        for i in range(n):
            for j in range(n):
                e = abs(M[s][i][j] - S[i][j])
                if e > max_err:
                    max_err, where = e, f"{s}: {fluors[i]} -> {dets[j]} ({M[s][i][j]:.4f} vs {S[i][j]:.4f})"
    ctx.info["spillover_max_abs_error"] = round(max_err, 5)
    ctx.info["spillover_worst_entry"] = where
    ctx.component("spillover_matrices", ramp(max_err, 0.05, 0.02), 0.20,
                  "max |entry - planted entry| over the 8x8 matrix of every sample; full credit <= 0.02, "
                  "none >= 0.05")

    # ---------------- frequencies.tsv --------------------------------------------------
    rows = ctx.table("frequencies.tsv", ["sample", "population", "parent", "count", "percent_of_parent"],
                     max_rows=5000)
    pop_names = {norm(p["name"]): p["name"] for p in pops}
    parent_of = {p["name"]: p["parent"] for p in pops}
    got = {}
    for r in rows:
        s, p = r["sample"], pop_names.get(norm(r["population"]))
        if s not in samples or p is None:
            continue                      # extra rows (other populations) are ignored
        ctx.gate((s, p) not in got, f"frequencies.tsv: {s}/{p} listed once")
        if norm(r["parent"]) != norm(parent_of[p]):
            raise CandidateError(f"gate failed: frequencies.tsv: parent of {p} is {parent_of[p]}")
        c = to_float(r["count"], f"count of {s}/{p}")
        if c < 0 or abs(c - round(c)) > 1e-6:
            raise CandidateError("gate failed: frequencies.tsv: count is a non-negative integer")
        pct = to_float(r["percent_of_parent"], f"percent_of_parent of {s}/{p}")
        if not 0.0 <= pct <= 100.0 + 1e-9:
            raise CandidateError("gate failed: frequencies.tsv: percent_of_parent lies in [0, 100]")
        got[(s, p)] = (int(round(c)), pct)
    ctx.gate(len(got) == len(samples) * len(pops),
             "frequencies.tsv holds one row for every sample x SOP population")

    def block(names, abs_tol, rel_tol, tag):
        per_sample, n_ok, worst = [], 0, 0.0
        for s in samples:
            ok = 0
            for p in names:
                t = truth["percent_of_parent"][s][p]
                band = max(abs_tol, rel_tol * t)
                dev = abs(got[(s, p)][1] - t) / band
                worst = max(worst, dev)
                if dev <= 1.0:
                    ok += 1
                else:
                    ctx.info.setdefault("out_of_tolerance", []).append(
                        f"{s}/{p}: {got[(s, p)][1]:.2f} vs {t:.2f}")
            n_ok += ok
            per_sample.append(ok / len(names))
        ctx.info[f"{tag}_fraction_in_tolerance_per_sample"] = [round(v, 4) for v in per_sample]
        ctx.info[f"{tag}_worst_error_over_tolerance"] = round(worst, 3)
        overall = n_ok / (len(samples) * len(names))
        ctx.info[f"{tag}_fraction_in_tolerance_overall"] = round(overall, 4)
        return min(per_sample), overall

    lin_min, lin_all = block(lineage, tol["percent_abs"], tol["percent_rel"], "lineage")
    act_min, act_all = block(fmo_gated, tol["fmo_abs"], tol["fmo_rel"], "fmo_gated")
    if "out_of_tolerance" in ctx.info:
        ctx.info["out_of_tolerance"] = ctx.info["out_of_tolerance"][:40]
    # every sample has to be right: a wholly wrong sample is not averaged away
    ctx.component("lineage_every_sample", ramp(lin_min, 10.0 / 12.0, 1.0), 0.25,
                  "min over samples of the fraction of lineage populations with percent_of_parent within "
                  "tolerance; full credit only if all are, none if two or more of the twelve are out")
    ctx.component("lineage_overall", ramp(lin_all, 0.5, 1.0), 0.10,
                  "fraction of all sample x lineage-population items within tolerance")

    good, worst_c = 0, 0.0
    for s in samples:
        t = truth["counts"][s]["lymphocytes"]
        dev = abs(got[(s, "lymphocytes")][0] - t) / max(t, 1)
        worst_c = max(worst_c, dev)
        good += dev <= tol["count_rel"]
    ctx.info["lymphocyte_count_worst_rel_error"] = round(worst_c, 4)
    ctx.component("lymphocyte_counts", ramp(good / len(samples), 0.5, 1.0), 0.10,
                  "fraction of samples whose live single lymphocyte count is within the relative tolerance")

    ctx.component("activation_every_sample", ramp(act_min, 0.5, 1.0), 0.25,
                  "min over samples of the fraction of the FMO-gated populations within tolerance; full "
                  "credit only if both are in every sample")
    ctx.component("activation_overall", act_all, 0.10,
                  "fraction of all sample x FMO-gated-population items within tolerance")


run(main, pass_threshold=0.9)
