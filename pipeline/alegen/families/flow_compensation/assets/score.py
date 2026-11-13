"""Verifier for the flow-cytometry compensation / gating family.

Truth is planted by the generator (reference/truth.json):
  * spillover            - the matrix the simulator mixed the fluorochromes with
  * counts / percent     - per sample, number of events of every SOP population counted from the
                           hidden class label of each simulated event
Nothing is derived from a compensation or gating run, so tolerances only have to cover estimator
noise of the controls and legitimate differences between correct gating strategies.
"""
import csv
import io
import re
import sys

sys.path.insert(0, "/tests")
from _scoring import CandidateError, ramp, run, to_float  # noqa: E402


def norm(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def read_matrix(ctx, fluors, dets):
    txt = ctx.text("spillover.csv", max_bytes=1_000_000)
    try:
        rows = [r for r in csv.reader(io.StringIO(txt)) if any(c.strip() for c in r)]
    except csv.Error as e:
        raise CandidateError(f"spillover.csv is not valid CSV: {e}")
    ctx.gate(2 <= len(rows) <= 60, "spillover.csv has a header and at most 59 rows")
    head = [norm(c) for c in rows[0][1:]]
    labels = [norm(r[0]) for r in rows[1:]]
    nf, nd = [norm(f) for f in fluors], [norm(d) for d in dets]
    if set(nd) <= set(head) and set(nf) <= set(labels):
        transposed = False
    elif set(nf) <= set(head) and set(nd) <= set(labels):
        transposed = True          # detectors as rows: same matrix, other layout
    else:
        raise CandidateError("gate failed: spillover.csv is labelled with the fluorochromes (rows) and "
                             "detectors (columns) of panel.tsv")
    ctx.gate(len(set(head)) == len(head) and len(set(labels)) == len(labels),
             "spillover.csv: row and column labels are unique")
    cells = {}
    for r, lab in zip(rows[1:], labels):
        for j, h in enumerate(head):
            if j + 1 >= len(r):
                raise CandidateError(f"spillover.csv: row {lab!r} is shorter than the header")
            cells[(lab, h)] = r[j + 1]
    M = []
    for f in nf:
        row = []
        for d in nd:
            raw = cells[(d, f)] if transposed else cells[(f, d)]
            row.append(to_float(raw, f"spillover.csv entry {f}/{d}"))
        M.append(row)
    return M


def main(ctx):
    truth = ctx.ref_json("truth.json")
    fluors, dets = truth["fluorochromes"], truth["detectors"]
    samples = truth["samples"]
    pops = truth["populations"]
    tol = truth["tolerance"]

    # ---------------- spillover.csv ----------------------------------------------------
    M = read_matrix(ctx, fluors, dets)
    S = truth["spillover"]
    n = len(fluors)
    max_err = max(abs(M[i][j] - S[i][j]) for i in range(n) for j in range(n))
    ctx.info["spillover_max_abs_error"] = round(max_err, 5)
    ctx.component("spillover_matrix", ramp(max_err, 0.05, 0.02), 0.25,
                  "max |entry - planted entry| over the 8x8 matrix; full credit <= 0.02, none >= 0.05")

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

    per_sample, n_ok, worst = [], 0, 0.0
    for s in samples:
        ok = 0
        for p in parent_of:
            t = truth["percent_of_parent"][s][p]
            band = max(tol["percent_abs"], tol["percent_rel"] * t)
            dev = abs(got[(s, p)][1] - t) / band
            worst = max(worst, dev)
            if dev <= 1.0:
                ok += 1
            else:
                ctx.info.setdefault("out_of_tolerance", []).append(
                    f"{s}/{p}: {got[(s, p)][1]:.2f} vs {t:.2f}")
        n_ok += ok
        per_sample.append(ok / len(parent_of))
        ctx.info[f"fraction_in_tolerance_{s}"] = round(ok / len(parent_of), 4)
    ctx.info["worst_error_over_tolerance"] = round(worst, 3)
    if "out_of_tolerance" in ctx.info:
        ctx.info["out_of_tolerance"] = ctx.info["out_of_tolerance"][:40]
    # every sample has to be right: a wholly wrong sample is not averaged away
    ctx.component("frequencies_every_sample", ramp(min(per_sample), 0.75, 1.0), 0.40,
                  "min over samples of the fraction of populations with percent_of_parent within "
                  "max(abs, rel) tolerance; full credit only if all are")
    overall = n_ok / (len(samples) * len(parent_of))
    ctx.info["fraction_in_tolerance_overall"] = round(overall, 4)
    ctx.component("frequencies_overall", ramp(overall, 0.5, 1.0), 0.20,
                  "fraction of all sample x population items within tolerance")

    good, worst_c = 0, 0.0
    for s in samples:
        t = truth["counts"][s]["lymphocytes"]
        dev = abs(got[(s, "lymphocytes")][0] - t) / max(t, 1)
        worst_c = max(worst_c, dev)
        good += dev <= tol["count_rel"]
    ctx.info["lymphocyte_count_worst_rel_error"] = round(worst_c, 4)
    ctx.component("lymphocyte_counts", ramp(good / len(samples), 0.5, 1.0), 0.15,
                  "fraction of samples whose live single lymphocyte count is within the relative tolerance")


run(main, pass_threshold=0.9)
