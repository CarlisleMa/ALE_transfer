"""Verifier for the trio lift-over / candidate review family.

Truth is planted by the generator (reference/truth.json):
  * sites       - every allele planted on GRCh38 (1-based pos, left-aligned minimal REF/ALT) with the
                  planted genotype of each sample; `scored` is false for the few alleles that the two
                  lift-over tools used as guards at generation time did not both reproduce
  * gray_keys   - GRCh38 keys of that gray zone (neither rewarded nor penalised)
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
    parts = re.split(r"[/|]", s.strip())
    if not parts or any(not p.isdigit() for p in parts):
        return None
    return tuple(sorted(int(p) for p in parts))


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
    tp = sum(1 for k in scored if k in sub)
    fp = sum(1 for k in sub if k not in scored and k not in gray) + n_dup
    fn = len(scored) - tp
    _, _, f1 = prf(tp, fp, fn)
    ctx.info.update({"lifted_tp": tp, "lifted_fp": fp, "lifted_fn": fn, "lifted_f1": round(f1, 4)})
    ctx.component("lifted_callset", ramp(f1, 0.90, 0.985), 0.15,
                  "F1 of literal (chrom,pos,ref,alt) keys vs scored planted GRCh38 alleles; full >= 0.985, none <= 0.90")

    conc = []
    for name in samples:
        ok = 0
        for k, s in scored.items():
            g = sub.get(k, {}).get(name) if k in sub else None
            if g is not None and g == gt_key(s["gt"][name]):
                ok += 1
        conc.append(ok / max(len(scored), 1))
        ctx.info[f"gt_concordance_{name}"] = round(conc[-1], 4)
    ctx.component("genotypes", ramp(min(conc), 0.90, 0.985), 0.10,
                  "min over samples of the fraction of scored planted alleles present with the planted genotype")

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
    ctx.component("unmapped", ramp(uf1, 0.50, 0.95), 0.10,
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
    ctx.component("candidates", ramp(cf1, 0.70, 1.0), 0.45,
                  "F1 of candidate keys vs the planted causal set; full credit only for the exact set, none <= 0.70")
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
