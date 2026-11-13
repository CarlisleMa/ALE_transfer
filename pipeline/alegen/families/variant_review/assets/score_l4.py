"""Verifier for level 4 of the somatic variant review family.

Truth (reference/truth.json) is the mechanism planted for every candidate while the reads were
synthesised, plus the realised molecule-level allele fraction of the true variants. Nothing is
derived from a review tool's output. The session and the images are checked for what a file can
show (written in IGV's session format with the SOP's settings; PNG of sufficient size that shows
data and is not filed twice); together they carry 0.15 and cannot lift a failing review.
"""
import hashlib
import math
import re
import struct
import sys
import xml.etree.ElementTree as ET
import zlib

sys.path.insert(0, "/tests")
from _scoring import CandidateError, ramp, run, to_float  # noqa: E402

FAIL_TAGS = {"GERMLINE", "LOW_SUPPORT", "LOW_BASE_QUALITY", "LOW_MAPPING_QUALITY", "MULTI_MISMATCH",
             "STRAND_BIAS", "READ_END", "DISCORDANT_PAIRS", "MNV_PARTNER"}
W_CALLS, W_TAGS, W_VAF, W_SESSION, W_SNAP = 0.40, 0.30, 0.15, 0.07, 0.08
PASS_MARK = 0.9


def mcc(tp, fp, tn, fn):
    den = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return (tp * tn - fp * fn) / den if den else 0.0


def _base(path):
    return (path or "").replace("\\", "/").rsplit("/", 1)[-1]


def session_score(ctx, first):
    """IGV session: written in IGV's format, reference, tracks, colouring, position."""
    try:
        txt = ctx.text("review_session.xml", max_bytes=2_000_000)
    except CandidateError:
        return 0.0, {"present": False}
    if "<!DOCTYPE" in txt or "<!ENTITY" in txt:      # no entity expansion games
        return 0.0, {"present": True, "doctype": True}
    try:
        root = ET.fromstring(txt.encode())
    except ET.ParseError:
        return 0.0, {"present": True, "xml": False}
    if root.tag != "Session":
        return 0.0, {"present": True, "root": root.tag[:30]}
    # IGV writes one <Track> element per track inside <Panel> elements, with its Java class
    tracks = {}
    for panel in root.iter("Panel"):
        for tr in panel.iter("Track"):
            if tr.get("clazz", "").endswith("AlignmentTrack"):
                tracks[_base(tr.get("id"))] = tr
    names = {_base(r.get("path")) for r in root.iter("Resource")}
    feature_ids = {_base(tr.get("id")) for panel in root.iter("Panel") for tr in panel.iter("Track")}
    checks = {"igv_format": "tumor.cram" in tracks and "normal.cram" in tracks}
    checks["reference"] = _base(root.get("genome")) == "ref.fa"
    checks["tracks"] = ({"tumor.cram", "normal.cram"} <= names
                        and bool({"candidates.bed", "candidates.vcf"} & names & feature_ids))

    def by_strand(tr):
        ro = tr.find("RenderOptions") if tr is not None else None
        return ro is not None and ro.get("colorOption") == "READ_STRAND"
    checks["coloured_by_read_strand"] = by_strand(tracks.get("tumor.cram")) and by_strand(tracks.get("normal.cram"))
    m = re.match(r"^(.+):([\d,]+)-([\d,]+)$", (root.get("locus") or "").strip())
    ok = False
    if m and m.group(1) == first["contig"]:
        a, b = int(m.group(2).replace(",", "")), int(m.group(3).replace(",", ""))
        # IGV writes the start of the range 0-based: allow that and a 201-bp window
        ok = a <= first["pos"] <= b and b - a <= 210
    checks["on_first_candidate"] = ok
    if not checks["igv_format"]:
        return 0.0, checks
    return sum(checks.values()) / len(checks), checks


def png_digest(raw):
    """(ok, digest) for one image: PNG, >= 800x300, shows data. Pixel data are inflated with a cap."""
    if len(raw) < 100 or raw[:8] != b"\x89PNG\r\n\x1a\n" or raw[12:16] != b"IHDR":
        return False, None
    w, h, depth, ctype, _, _, interlace = struct.unpack(">IIBBBBB", raw[16:29])
    if not (800 <= w <= 10000 and 300 <= h <= 20000) or depth != 8 or ctype not in (2, 6) or interlace:
        return False, None
    stride = w * (3 if ctype == 2 else 4) + 1
    idat, off = [], 8
    while off + 8 <= len(raw):
        ln, typ = struct.unpack(">I4s", raw[off:off + 8])
        if typ == b"IDAT":
            idat.append(raw[off + 8:off + 8 + ln])
        off += 12 + ln
    d = zlib.decompressobj()
    try:
        data = d.decompress(b"".join(idat), stride * h + 1)
    except zlib.error:
        return False, None
    if len(data) != stride * h:
        return False, None
    rows = {data[i:i + stride] for i in range(0, len(data), stride)}
    if len(rows) < 40:                                # an empty panel or a flat image has a handful of row patterns
        return False, None
    return True, hashlib.sha256(data).hexdigest()


def snapshot_score(ctx, ids):
    digests, seen = {}, {}
    for v in ids:
        try:
            raw = ctx.path(f"snapshots/{v}.png", max_bytes=30_000_000).read_bytes()
        except CandidateError:
            continue
        ok, dg = png_digest(raw)
        if ok:
            digests[v] = dg
            seen[dg] = seen.get(dg, 0) + 1
    good = sum(1 for v, dg in digests.items() if seen[dg] == 1)
    return good / len(ids), {"valid": len(digests), "valid_and_unique": good, "candidates": len(ids)}


def main(ctx):
    ref = ctx.ref_json("truth.json")
    truth = ref["variants"]
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
    fails = [v for v, t in truth.items() if t["call"] == "fail"]
    tag_acc = sum(1 for v in fails if sub[v][1] == truth[v]["tag"]) / len(fails)
    som = [v for v, t in truth.items() if t["call"] == "somatic"]
    frac = sum(1 for v in som if sub[v][0] == "somatic" and abs(sub[v][2] - truth[v]["tumor_vaf"]) <= 0.05) / len(som)
    c_calls, c_tags, c_vaf = ramp(m, 0.86, 0.965), ramp(tag_acc, 0.85, 0.97), ramp(frac, 0.80, 0.96)
    review = (W_CALLS * c_calls + W_TAGS * c_tags + W_VAF * c_vaf) / (W_CALLS + W_TAGS + W_VAF)

    s_val, s_info = session_score(ctx, ref["first_candidate"])
    p_val, p_info = snapshot_score(ctx, list(truth))
    # documentation cannot lift a failing review over the pass mark
    docs = W_SESSION * s_val + W_SNAP * p_val
    total = (W_CALLS + W_TAGS + W_VAF) * review + docs
    scale = 1.0
    if review < PASS_MARK and total > PASS_MARK - 0.01 and docs > 0:
        scale = max(0.0, (PASS_MARK - 0.01 - (W_CALLS + W_TAGS + W_VAF) * review) / docs)

    ctx.info.update({"tp": tp, "fp": fp, "tn": tn, "fn": fn, "mcc": round(m, 4),
                     "failure_tag_accuracy": round(tag_acc, 4), "somatic_with_correct_vaf": round(frac, 4),
                     "review_fraction": round(review, 4), "session": s_info, "snapshots": p_info,
                     "documentation_scale": round(scale, 4),
                     "wrong": sorted(v for v, t in truth.items() if sub[v][1] != t["tag"])[:80]})
    ctx.component("calls", c_calls, W_CALLS, "Matthews correlation of call vs planted truth")
    ctx.component("failure_tags", c_tags, W_TAGS, "fraction of planted failures that carry the planted tag")
    ctx.component("allele_fraction", c_vaf, W_VAF,
                  "fraction of planted somatic variants called somatic with tumor_vaf within 0.05 of realised")
    ctx.component("igv_session", s_val * scale, W_SESSION,
                  "IGV-format session: reference, tracks, colouring by read strand, view on the first candidate")
    ctx.component("snapshots", p_val * scale, W_SNAP, "candidates with their own valid IGV image")


run(main, pass_threshold=PASS_MARK)
