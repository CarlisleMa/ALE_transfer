"""Verifier for the somatic variant review family.

Truth (reference/truth.json) is the mechanism planted for every candidate while the reads were
synthesised, plus the realised molecule-level allele fraction of the true variants. Nothing is
derived from a review tool's output.
"""
import math
import sys
import xml.etree.ElementTree as ET

sys.path.insert(0, "/tests")
from _scoring import CandidateError, ramp, run, to_float  # noqa: E402

FAIL_TAGS = {"GERMLINE", "LOW_SUPPORT", "LOW_BASE_QUALITY", "LOW_MAPPING_QUALITY", "MULTI_MISMATCH",
             "STRAND_BIAS", "READ_END"}


def mcc(tp, fp, tn, fn):
    den = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return (tp * tn - fp * fn) / den if den else 0.0


def session_score(ctx):
    """Structure-only check of the IGV session (cannot prove IGV was used; weight 0.05)."""
    try:
        txt = ctx.text("igv_session.xml", max_bytes=1_000_000)
    except CandidateError:
        return 0.0
    if "<!DOCTYPE" in txt or "<!ENTITY" in txt:      # no entity expansion games
        return 0.0
    try:
        root = ET.fromstring(txt.encode())
    except ET.ParseError:
        return 0.0
    if root.tag != "Session":
        return 0.0
    paths = [r.get("path", "") for r in root.iter("Resource")]
    names = {p.replace("\\", "/").rsplit("/", 1)[-1] for p in paths}
    ok = [
        "ref.fa" in root.get("genome", ""),
        "tumor.bam" in names,
        "normal.bam" in names,
        bool({"candidates.bed", "candidates.vcf"} & names),
    ]
    return sum(ok) / len(ok)


def main(ctx):
    truth = ctx.ref_json("truth.json")["variants"]
    rows = ctx.table("review.tsv", ["variant_id", "call", "tag", "tumor_alt", "tumor_depth", "tumor_vaf",
                                    "normal_alt", "normal_depth"], max_rows=len(truth) + 10)
    sub = {}
    for r in rows:
        v = r["variant_id"]
        ctx.gate(v in truth, f"review.tsv: variant id {v[:30]!r} is a candidate")
        ctx.gate(v not in sub, f"review.tsv: {v[:30]!r} listed once")
        call, tag = r["call"], r["tag"]
        if call not in ("somatic", "fail"):
            raise CandidateError("review.tsv: call must be 'somatic' or 'fail'")
        if tag != "PASS" and tag not in FAIL_TAGS:
            raise CandidateError(f"review.tsv: unknown tag {tag[:30]!r}")
        if (call == "somatic") != (tag == "PASS"):
            raise CandidateError("review.tsv: tag PASS goes with call somatic, a failure tag with call fail")
        for col in ("tumor_alt", "tumor_depth", "normal_alt", "normal_depth"):
            to_float(r[col], col)
        vaf = to_float(r["tumor_vaf"], "tumor_vaf")
        if not 0.0 <= vaf <= 1.0:
            raise CandidateError("review.tsv: tumor_vaf must be within [0, 1]")
        sub[v] = (call, tag, vaf)
    ctx.gate(len(sub) == len(truth), "review.tsv lists every candidate exactly once")

    tp = sum(1 for v, t in truth.items() if t["call"] == "somatic" and sub[v][0] == "somatic")
    fn = sum(1 for v, t in truth.items() if t["call"] == "somatic" and sub[v][0] == "fail")
    fp = sum(1 for v, t in truth.items() if t["call"] == "fail" and sub[v][0] == "somatic")
    tn = sum(1 for v, t in truth.items() if t["call"] == "fail" and sub[v][0] == "fail")
    m = mcc(tp, fp, tn, fn)
    ctx.info.update({"tp": tp, "fp": fp, "tn": tn, "fn": fn, "mcc": round(m, 4)})
    ctx.component("calls", ramp(m, 0.70, 0.96), 0.45, "Matthews correlation of call vs planted truth")

    fails = [v for v, t in truth.items() if t["call"] == "fail"]
    tag_acc = sum(1 for v in fails if sub[v][1] == truth[v]["tag"]) / len(fails)
    ctx.info["failure_tag_accuracy"] = round(tag_acc, 4)
    ctx.component("failure_tags", ramp(tag_acc, 0.75, 0.97), 0.30,
                  "fraction of planted failures that carry the planted tag")

    som = [v for v, t in truth.items() if t["call"] == "somatic"]
    good = sum(1 for v in som if sub[v][0] == "somatic" and abs(sub[v][2] - truth[v]["tumor_vaf"]) <= 0.05)
    frac = good / len(som)
    ctx.info["somatic_with_correct_vaf"] = round(frac, 4)
    ctx.component("allele_fraction", ramp(frac, 0.70, 0.95), 0.20,
                  "fraction of planted somatic variants called somatic with tumor_vaf within 0.05 of realised")

    ctx.component("igv_session", session_score(ctx), 0.05,
                  "IGV session references the reference, both BAMs and the candidate track")


run(main, pass_threshold=0.9)
