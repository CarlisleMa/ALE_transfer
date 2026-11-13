"""Verifier for the 16S V4 amplicon family.

Truth is planted by the generator (reference/truth.json):
  * planted[].seq        exact amplicon (primer-free, 515F -> 806R) of every strain put into the libraries
  * planted[].counts     read pairs emitted from that sequence per library (recombinant pairs count
                         for nobody unless the recombinant equals a planted sequence)
  * planted[].genus      genus of the source record, or Unclassified where the generator removed every
                         close reference from the shipped database
  * planted[].contaminant  the strain was added as reagent background
  * planted[].scored     False for strains planted below what denoising can resolve (< 100 read pairs
                         in the whole batch): reporting them is not an error, missing them is not a miss
Nothing is derived from a tool run.
"""
import math
import sys

sys.path.insert(0, "/tests")
from _scoring import CandidateError, ramp, run, to_float  # noqa: E402

MAX_ASV = 5000


def mcc(tp, fp, tn, fn):
    den = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return (tp * tn - fp * fn) / den if den else 0.0


def parse_fasta(txt):
    recs, name, chunks = [], None, []
    for line in txt.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                recs.append((name, "".join(chunks)))
            toks = line[1:].split()
            if not toks:
                raise CandidateError("asv.fasta: record without an id")
            name, chunks = toks[0], []
            if len(recs) > MAX_ASV:
                raise CandidateError(f"asv.fasta has more than {MAX_ASV} records")
        else:
            if name is None:
                raise CandidateError("asv.fasta does not start with a '>' header")
            chunks.append(line.upper())
    if name is not None:
        recs.append((name, "".join(chunks)))
    return recs


def main(ctx):
    truth = ctx.ref_json("truth.json")
    samples = [s["id"] for s in truth["samples"]]
    stype = {s["id"]: s["type"] for s in truth["samples"]}
    planted = truth["planted"]
    by_seq = {p["seq"]: p for p in planted}

    # ---------------- asv.fasta --------------------------------------------------------
    recs = parse_fasta(ctx.text("asv.fasta", max_bytes=10_000_000))
    ctx.gate(1 <= len(recs) <= MAX_ASV, f"asv.fasta holds between 1 and {MAX_ASV} records")
    ids = [r[0] for r in recs]
    ctx.gate(len(set(ids)) == len(ids), "asv.fasta: ASV ids are unique")
    for i, s in recs:
        if not (20 <= len(s) <= 1000) or set(s) - set("ACGT"):
            raise CandidateError(f"asv.fasta: sequence of {i[:30]!r} must be 20-1000 nt of A, C, G, T")
    ctx.gate(len({s for _, s in recs}) == len(recs), "asv.fasta: every sequence is reported once")
    idset = set(ids)

    # ---------------- asv_table.tsv ----------------------------------------------------
    rows = ctx.table("asv_table.tsv", ["asv_id"] + samples, max_rows=MAX_ASV + 10)
    counts = {}
    for r in rows:
        i = r["asv_id"] or ""
        ctx.gate(i in idset, f"asv_table.tsv: id {i[:30]!r} is an ASV of asv.fasta")
        ctx.gate(i not in counts, f"asv_table.tsv: {i[:30]!r} listed once")
        vals = []
        for s in samples:
            v = to_float(r[s], f"asv_table.tsv {i[:30]} {s}")
            if v < 0 or abs(v - round(v)) > 1e-6:
                raise CandidateError("asv_table.tsv: counts must be non-negative integers")
            vals.append(int(round(v)))
        counts[i] = vals
    ctx.gate(len(counts) == len(ids), "asv_table.tsv lists every ASV of asv.fasta exactly once")

    # ---------------- taxonomy.tsv -----------------------------------------------------
    rows = ctx.table("taxonomy.tsv", ["asv_id", "genus"], max_rows=MAX_ASV + 10)
    genus = {}
    for r in rows:
        i = r["asv_id"] or ""
        ctx.gate(i in idset, f"taxonomy.tsv: id {i[:30]!r} is an ASV of asv.fasta")
        ctx.gate(i not in genus, f"taxonomy.tsv: {i[:30]!r} listed once")
        genus[i] = (r["genus"] or "").strip()
    ctx.gate(len(genus) == len(ids), "taxonomy.tsv lists every ASV of asv.fasta exactly once")

    # ---------------- contaminants.txt -------------------------------------------------
    lines = [l.strip() for l in ctx.text("contaminants.txt", max_bytes=1_000_000).splitlines() if l.strip()]
    if lines and lines[0].lower() == "asv_id":
        lines = lines[1:]
    flagged = set()
    for l in lines:
        ctx.gate(l in idset, f"contaminants.txt: id {l[:30]!r} is an ASV of asv.fasta")
        flagged.add(l)

    # ---------------- sequences ----------------------------------------------------------
    asv_of = {}                       # pid -> submitted ASV id
    n_match = 0
    for i, s in recs:
        p = by_seq.get(s)
        if p is not None:
            asv_of[p["pid"]] = i
            n_match += 1
    scored = [p for p in planted if p["scored"]]
    recall = sum(1 for p in scored if p["pid"] in asv_of) / len(scored)
    precision = n_match / len(recs)
    ctx.info.update({"n_asv": len(recs), "n_scored_planted": len(scored), "seq_recall": round(recall, 4),
                     "seq_precision": round(precision, 4)})
    ctx.component("asv_sequences", ramp(recall, 0.85, 0.95) * ramp(precision, 0.60, 0.90), 0.35,
                  "recall of scored planted sequences (full >= 0.95, none <= 0.85) x precision of reported "
                  "sequences (full >= 0.90, none <= 0.60), exact matches only")

    # ---------------- counts ---------------------------------------------------------------
    def agreement(j):
        obs_total = sum(counts[i][j] for i in ids if i not in flagged)
        true_total = sum(p["counts"][j] for p in planted if not p["contaminant"])
        inter = 0
        for p in planted:
            i = asv_of.get(p["pid"])
            if p["contaminant"] or i is None or i in flagged:
                continue
            inter += min(counts[i][j], p["counts"][j])
        return inter / max(obs_total, true_total, 1)

    per_lib = {s: agreement(j) for j, s in enumerate(samples)}
    study = [per_lib[s] for s in samples if stype[s] == "study"]
    mock = [per_lib[s] for s in samples if stype[s] == "mock"]
    ctx.info["count_agreement"] = {s: round(v, 4) for s, v in per_lib.items() if stype[s] != "negative_control"}
    ctx.info["count_agreement_min_study"] = round(min(study), 4)
    ctx.info["count_agreement_mock"] = round(min(mock), 4)
    # the minimum, not the mean: one specimen with a wrong profile is a wrong deliverable
    ctx.component("sample_counts", ramp(min(study), 0.60, 0.80), 0.30,
                  "min over study libraries of sum(min(obs,true))/max(sum obs,sum true) after removing the "
                  "listed contaminants; full >= 0.80, none <= 0.60")

    # ---------------- contaminant flags -----------------------------------------------------
    tp = fp = tn = fn = 0
    for p in scored:
        pred = asv_of.get(p["pid"]) in flagged if p["pid"] in asv_of else False
        if p["contaminant"]:
            tp, fn = tp + pred, fn + (not pred)
        else:
            fp, tn = fp + pred, tn + (not pred)
    ctx.info.update({"contam_tp": tp, "contam_fp": fp, "contam_fn": fn})
    if tp + fn == 0:
        # Batch without reagent background: the right answer is "none", which costs nothing to write.
        # Credit is therefore given per scored planted strain that is reported AND left unflagged, so
        # an answer that recovered nothing earns nothing here (same convention as the genus component).
        clean = sum(1 for p in scored if p["pid"] in asv_of and asv_of[p["pid"]] not in flagged) / len(scored)
        ctx.info["contam_unflagged_fraction"] = round(clean, 4)
        ctx.component("contaminants", ramp(clean, 0.80, 0.95), 0.15,
                      "no reagent background in this batch: fraction of scored planted sequences that are "
                      "reported and not listed as contaminant; full >= 0.95, none <= 0.80")
    else:
        m = mcc(tp, fp, tn, fn)
        ctx.info["contam_mcc"] = round(m, 4)
        ctx.component("contaminants", ramp(m, 0.50, 0.85), 0.15,
                      "Matthews correlation of the contaminant flags over scored planted sequences; full >= 0.85")

    # ---------------- genus -----------------------------------------------------------------
    ok = sum(1 for p in scored if p["pid"] in asv_of
             and genus[asv_of[p["pid"]]].lower() == p["genus"].lower())
    acc = ok / len(scored)
    ctx.info["genus_accuracy"] = round(acc, 4)
    ctx.component("genus", ramp(acc, 0.80, 0.95), 0.10,
                  "fraction of scored planted sequences reported with the planted genus; full >= 0.95")

    # ---------------- mock library ----------------------------------------------------------
    ctx.component("mock_counts", ramp(min(mock), 0.60, 0.80), 0.10,
                  "count agreement of the mock library; full >= 0.80, none <= 0.60")


run(main, pass_threshold=0.9)
