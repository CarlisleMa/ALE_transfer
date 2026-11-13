"""Verifier for the trio lift-over / candidate review family, level 4.

Truth is planted by the generator (reference/truth.json):
  * sites       - every allele by which the family differs from GRCh38 (1-based pos, left-aligned minimal
                  REF/ALT) with the planted number of ALT copies of each sample. `cat` is `standard` for
                  alleles at positions where both assemblies carry the same base, `refdiff` for
                  alleles at positions where the provider's reference base differs from GRCh38 (with or
                  without a record in the provider's file) and `inverted` for alleles inside the block
                  that is aligned on the minus strand. `scored` is false for the few standard
                  alleles that the two lift-over tools used as guards at generation time did not both
                  reproduce, and for assembly differences at neighbouring bases.
  * gray_keys / gray_positions - that gray zone (neither rewarded nor penalised)
  * refdiff_positions - every aligned GRCh38 position whose base differs from the provider's reference
  * inverted_interval - first and last GRCh38 position of the minus-strand block
  * unmapped    - alleles planted in hg19 sequence without GRCh38 alignment (provider coordinates)
  * candidates  - the planted causal variants with gene, consequence class and inheritance model
Nothing is recomputed from the submission with a domain tool; the VCF is parsed as text after a
streaming, size-capped decompression (a real submission is < 1 MB uncompressed).
"""
import re
import sys
import zlib

sys.path.insert(0, "/tests")
from _scoring import CandidateError, prf, ramp, run  # noqa: E402

MAX_RECORDS = 60000
MODELS = {"de_novo", "compound_het", "homozygous_recessive"}
ALLELE = re.compile(r"^[ACGTN]+$")


def gunzip_capped(raw: bytes, cap: int, what: str) -> str:
    """Decompress (multi-member) gzip/BGZF without ever holding more than `cap` bytes of output."""
    out, total, data = [], 0, raw
    while data:
        d = zlib.decompressobj(31)
        buf = data
        try:
            while buf and not d.eof:
                chunk = d.decompress(buf, 1 << 20)
                total += len(chunk)
                if total > cap:
                    raise CandidateError(f"{what} decompresses to more than {cap} bytes")
                out.append(chunk)
                buf = d.unconsumed_tail
        except zlib.error:
            raise CandidateError(f"{what} is not valid gzip/BGZF")
        if not d.eof:
            raise CandidateError(f"{what} is a truncated gzip/BGZF file")
        data = d.unused_data
    try:
        return b"".join(out).decode("utf-8")
    except UnicodeDecodeError:
        raise CandidateError(f"{what} is not UTF-8 text")


def bare(c: str) -> str:
    c = c.strip()
    return c[3:] if c.lower().startswith("chr") else c


def gt_key(s: str):
    """Number of chromosome copies that carry the ALT allele of the record (None = no call)."""
    parts = re.split(r"[/|]", s.strip())
    if not parts or all(p == "." for p in parts) or any(p != "." and not p.isdigit() for p in parts):
        return None
    return sum(1 for p in parts if p == "1")


def cons_tokens(s: str) -> set:
    out = set()
    for t in re.split(r"[&,;/\s]+", s.strip().lower()):
        t = t.lstrip("*")
        if t.endswith("_variant"):
            t = t[:-len("_variant")]
        if t:
            out.add(t)
    return out


def to_pos(s, what):
    s = str(s).strip()
    if not s.isdigit() or len(s) > 10:
        raise CandidateError(f"{what}: position is not a positive integer: {s[:20]!r}")
    return int(s)


def main(ctx):
    truth = ctx.ref_json("truth.json")
    samples = truth["samples"]

    # ---------------- trio.GRCh38.vcf.gz ------------------------------------------------
    vcf_path = ctx.path("trio.GRCh38.vcf.gz", max_bytes=10_000_000)
    raw = vcf_path.read_bytes()
    if raw[:2] != b"\x1f\x8b":
        raise CandidateError("trio.GRCh38.vcf.gz is not gzip/BGZF compressed")
    text = gunzip_capped(raw, 50_000_000, "trio.GRCh38.vcf.gz")
    has_index = False
    for ext in (".tbi", ".csi"):
        try:
            ctx.path("trio.GRCh38.vcf.gz" + ext, max_bytes=20_000_000)
            has_index = True
        except CandidateError:
            pass
    ctx.gate(has_index, "trio.GRCh38.vcf.gz has a non-empty .tbi or .csi index next to it")
    cols, sub, n_rec, n_dup = None, {}, 0, 0
    for ln in text.split("\n"):
        if not ln or ln.startswith("##"):
            continue
        f = ln.rstrip("\r").split("\t")
        if ln.startswith("#CHROM"):
            cols = f
            ctx.gate(len(f) == 12 and sorted(f[9:]) == sorted(samples),
                     "VCF header lists exactly the three input samples under their original names")
            continue
        if cols is None:
            raise CandidateError("trio.GRCh38.vcf.gz: record before the #CHROM header line")
        n_rec += 1
        if n_rec > MAX_RECORDS:
            raise CandidateError(f"trio.GRCh38.vcf.gz has more than {MAX_RECORDS} records")
        if len(f) != len(cols):
            raise CandidateError("trio.GRCh38.vcf.gz: record with a wrong number of columns")
        key = (f[0], to_pos(f[1], "trio.GRCh38.vcf.gz"), f[3].upper(), f[4].upper())
        fmt = f[8].split(":")
        gts = {}
        if "GT" in fmt:
            gi = fmt.index("GT")
            for name, val in zip(cols[9:], f[9:]):
                v = val.split(":")
                gts[name] = gt_key(v[gi]) if gi < len(v) else None
        if key in sub:
            n_dup += 1
        else:
            sub[key] = gts
    ctx.gate(cols is not None, "trio.GRCh38.vcf.gz has a #CHROM header line")

    scored = {("chr22", s["pos"], s["ref"], s["alt"]): s for s in truth["sites"] if s["scored"]}
    gray = {("chr22", p, r, a) for p, r, a in truth["gray_keys"]}
    gray |= {("chr22", s["pos"], s["ref"], s["alt"]) for s in truth["sites"] if not s["scored"]}
    gray_pos = set(truth["gray_positions"])
    rd_pos = set(truth["refdiff_positions"])
    inv_lo, inv_hi = truth["inverted_interval"]
    cats = ("standard", "refdiff", "inverted")
    tp, fp, fn = dict.fromkeys(cats, 0), dict.fromkeys(cats, 0), dict.fromkeys(cats, 0)
    fp["standard"] += n_dup
    n_noncarrier = 0
    for k, gts in sub.items():
        counts = [gts.get(name) for name in samples]
        if all(c == 0 for c in counts):
            n_noncarrier += 1                    # nobody carries the ALT allele: the record claims nothing
            continue
        span = range(k[1], k[1] + len(k[2]))
        if k in scored:
            cat = scored[k]["cat"]
            if counts == [gt_key(scored[k]["gt"][name]) for name in samples]:
                tp[cat] += 1
            else:
                fp[cat] += 1                     # right allele, wrong genotypes: a wrong record
        elif k in gray or (k[0] == "chr22" and any(q in gray_pos for q in span)):
            continue
        else:
            on_chr = k[0] == "chr22"
            fp["refdiff" if on_chr and any(q in rd_pos for q in span) else
               ("inverted" if on_chr and inv_lo <= k[1] <= inv_hi else "standard")] += 1
    f1 = {}
    for cat in cats:
        fn[cat] = sum(1 for s in scored.values() if s["cat"] == cat) - tp[cat]
        f1[cat] = prf(tp[cat], fp[cat], fn[cat])[2]
        ctx.info.update({f"{cat}_tp": tp[cat], f"{cat}_fp": fp[cat], f"{cat}_fn": fn[cat], f"{cat}_f1": round(f1[cat], 4)})
    ctx.info["records_without_carrier_ignored"] = n_noncarrier
    ctx.component("callset_standard", ramp(f1["standard"], 0.90, 0.985), 0.08,
                  "F1 over alleles at positions where both assemblies agree (key and ALT copies of all three samples); "
                  "full >= 0.985, none <= 0.90")
    ctx.component("callset_reference_differences", ramp(f1["refdiff"], 0.80, 0.97), 0.12,
                  "F1 over alleles at positions where the provider's reference base differs from GRCh38; "
                  "full >= 0.97, none <= 0.80")
    ctx.component("callset_inverted_block", ramp(f1["inverted"], 0.70, 0.95), 0.07,
                  "F1 over alleles inside the block that the chain aligns on the minus strand; full >= 0.95, none <= 0.70")

    # ---------------- unmapped.tsv ------------------------------------------------------
    urows = ctx.table("unmapped.tsv", ["chrom_in", "pos_in", "ref", "alt"], max_rows=MAX_RECORDS)
    usub = set()
    for r in urows:
        usub.add((bare(r["chrom_in"]), to_pos(r["pos_in"], "unmapped.tsv"), r["ref"].upper(), r["alt"].upper()))
    u_scored = {("22", u["pos"], u["ref"], u["alt"]) for u in truth["unmapped"] if u["scored"]}
    u_gray = {("22", u["pos"], u["ref"], u["alt"]) for u in truth["unmapped"] if not u["scored"]}
    utp = len(usub & u_scored)
    ufp = len(usub - u_scored - u_gray)
    ufn = len(u_scored - usub)
    if u_scored:
        _, _, uf1 = prf(utp, ufp, ufn)
    else:
        uf1 = 1.0 if ufp == 0 else 0.0
    ctx.info.update({"unmapped_tp": utp, "unmapped_fp": ufp, "unmapped_fn": ufn, "unmapped_f1": round(uf1, 4)})
    ctx.component("unmapped", ramp(uf1, 0.50, 0.95), 0.08,
                  "F1 of (pos_in,ref,alt) vs planted hg19-only alleles; full >= 0.95, none <= 0.50")

    # ---------------- candidates.tsv ----------------------------------------------------
    crows = ctx.table("candidates.tsv", ["chrom", "pos", "ref", "alt", "gene", "consequence", "model"], max_rows=201)
    ctx.gate(len(crows) <= 200, "candidates.tsv has at most 200 rows")
    csub, cdup = {}, 0
    for r in crows:
        model = (r["model"] or "").strip().lower()
        if model not in MODELS:
            raise CandidateError(f"candidates.tsv: model must be one of {sorted(MODELS)}")
        ref, alt = (r["ref"] or "").upper(), (r["alt"] or "").upper()
        if not ALLELE.match(ref) or not ALLELE.match(alt):
            raise CandidateError("candidates.tsv: ref/alt must be plain nucleotide strings")
        key = (r["chrom"], to_pos(r["pos"], "candidates.tsv"), ref, alt)
        if key in csub:
            cdup += 1
        else:
            csub[key] = ((r["gene"] or "").strip(), cons_tokens(r["consequence"] or ""), model)
    want = {("chr22", c["pos"], c["ref"], c["alt"]): c for c in truth["candidates"]}
    ctp = sum(1 for k in want if k in csub)
    cfp = sum(1 for k in csub if k not in want) + cdup
    cfn = len(want) - ctp
    _, _, cf1 = prf(ctp, cfp, cfn)
    ctx.info.update({"cand_tp": ctp, "cand_fp": cfp, "cand_fn": cfn, "cand_f1": round(cf1, 4)})
    ctx.component("candidates", ramp(cf1, 0.87, 1.0), 0.45,
                  "F1 of candidate keys vs the planted causal set; full credit only for the exact set, none <= 0.87")
    good = 0
    for k, c in want.items():
        if k in csub:
            gene, toks, model = csub[k]
            if gene == c["gene"] and c["consequence"] in toks and model == c["model"]:
                good += 1
    ctx.info["cand_annotated_ok"] = good
    ctx.component("candidate_annotation", good / len(want), 0.20,
                  "fraction of planted causal variants reported with the planted gene, consequence class and model")


run(main, pass_threshold=0.9)
