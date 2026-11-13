"""Verifier for the somatic_tn family (stdlib only).

Truth is planted by the generator (reference/truth.json):
  * the library of every patient's tumour and of the normal from the same individual,
  * `scored`: planted somatic variants with enough simulated evidence (see the task card),
  * `gray`:   planted somatic variants below that floor or in the interior of the duplicated
              segment - calling them is neither rewarded nor penalised.
Every other PASS call inside the targets is a false positive. Nothing here is a tool output.
"""
import sys

sys.path.insert(0, "/tests")
from _scoring import CandidateError, ramp, run  # noqa: E402

MAX_RECORDS = 20000


def read_fasta(path):
    seq = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq).upper()


def normalise(seq, pos0, ref, alt):
    """Parsimonious left-aligned form (the `bcftools norm` convention); pos0 is 0-based."""
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos0 = ref[1:], alt[1:], pos0 + 1
    if len(ref) == len(alt):
        return pos0, ref, alt
    while pos0 > 0 and ref[-1] == alt[-1]:
        pos0 -= 1
        b = seq[pos0]
        ref, alt = b + ref[:-1], b + alt[:-1]
    return pos0, ref, alt


def parse_vcf(text, rel, contig, seq):
    """-> {(pos0, ref, alt): af or None} for PASS records, after normalisation."""
    calls = {}
    n = 0
    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        n += 1
        if n > MAX_RECORDS:
            raise CandidateError(f"{rel} has more than {MAX_RECORDS} records")
        f = line.rstrip("\r\n").split("\t")
        if len(f) < 8:
            raise CandidateError(f"{rel}: a record has fewer than 8 tab-separated columns")
        chrom, pos, _id, ref, alts, _qual, flt, info = f[:8]
        if flt.strip() not in ("PASS", "."):
            continue
        if chrom != contig:
            raise CandidateError(f"{rel}: contig {chrom[:30]!r} is not the panel reference contig")
        try:
            pos0 = int(pos) - 1
        except ValueError:
            raise CandidateError(f"{rel}: POS is not an integer")
        ref = ref.strip().upper()
        if pos0 < 0 or not ref or pos0 + len(ref) > len(seq):
            raise CandidateError(f"{rel}: record at POS {pos[:12]} lies outside the reference")
        if seq[pos0:pos0 + len(ref)] != ref:
            raise CandidateError(f"{rel}: REF at POS {pos0 + 1} does not match reference/ref.fa")
        afs = []
        for item in info.split(";"):
            if item.startswith("AF="):
                afs = item[3:].split(",")
        for i, alt in enumerate(alts.split(",")):
            alt = alt.strip().upper()
            if not alt or set(alt) - set("ACGT") or alt == ref:
                continue                                   # symbolic / missing alleles are not calls
            key = normalise(seq, pos0, ref, alt)
            af = None
            if afs:
                try:
                    af = float(afs[min(i, len(afs) - 1)])
                    if not (0.0 <= af <= 1.0):
                        af = None
                except ValueError:
                    af = None
            if key not in calls or calls[key] is None:
                calls[key] = af
    return calls


def main(ctx):
    truth = ctx.ref_json("truth.json")
    seq = read_fasta(ctx.reference / "ref.fa")
    contig = truth["contig"]
    targets = truth["targets"]
    pats = truth["patients"]

    # ---------------- pairing.tsv ------------------------------------------------------
    rows = ctx.table("pairing.tsv", ["patient", "tumour_sample", "normal_sample"], max_rows=1000)
    pairing = {}
    for r in rows:
        p = str(r.get("patient") or "")
        ctx.gate(p in pats, f"pairing.tsv: patient {p[:30]!r} is in the manifest")
        ctx.gate(p not in pairing, f"pairing.tsv: patient {p[:30]!r} listed once")
        pairing[p] = (str(r.get("tumour_sample") or ""), str(r.get("normal_sample") or ""))
    ctx.gate(len(pairing) == len(pats), "pairing.tsv lists every patient exactly once")
    wrong = [p for p in pats if pairing[p] != (pats[p]["tumour_sample"], pats[p]["normal_sample"])]
    ctx.info["pairing_wrong"] = wrong

    # ---------------- somatic/<patient>.vcf.gz -----------------------------------------
    def in_targets(pos0):
        return any(s <= pos0 < e for s, e in targets)

    precs, recs, af_meds, counts = [], [], [], {}
    for p in sorted(pats):
        rel = f"somatic/{p}.vcf.gz"
        calls = parse_vcf(ctx.text(rel, max_bytes=3_000_000), rel, contig, seq)
        counts[p] = (sum(1 for k in calls if len(k[1]) == 1 and len(k[2]) == 1),
                     sum(1 for k in calls if len(k[1]) != len(k[2])))
        scored = {(v["pos"] - 1, v["ref"], v["alt"]): v["af"] for v in pats[p]["scored"]}
        gray = {(v["pos"] - 1, v["ref"], v["alt"]) for v in pats[p]["gray"]}
        evaluated = {k: af for k, af in calls.items() if in_targets(k[0])}
        tp = [k for k in evaluated if k in scored]
        fp = [k for k in evaluated if k not in scored and k not in gray]
        prec = len(tp) / (len(tp) + len(fp)) if tp else 0.0
        rec = len(tp) / len(scored)
        errs = sorted(abs(evaluated[k] - scored[k]) if evaluated[k] is not None else 1.0 for k in tp)
        med = errs[len(errs) // 2] if errs else 1.0
        precs.append(prec)
        recs.append(rec)
        af_meds.append(med)
        ctx.info[p] = {"tp": len(tp), "fp": len(fp), "fn": len(scored) - len(tp),
                       "gray_called": sum(1 for k in evaluated if k in gray),
                       "outside_targets": len(calls) - len(evaluated),
                       "precision": round(prec, 4), "recall": round(rec, 4), "af_median_abs_error": round(med, 4),
                       "false_positives": sorted(f"{k[0] + 1}:{k[1]}>{k[2]}" for k in fp)[:60],
                       "missed": sorted(f"{k[0] + 1}:{k[1]}>{k[2]}" for k in scored if k not in evaluated)[:60]}
    min_prec, min_rec, worst_af = min(precs), min(recs), max(af_meds)
    ctx.info.update({"min_precision": round(min_prec, 4), "min_recall": round(min_rec, 4),
                     "worst_af_median_abs_error": round(worst_af, 4)})
    base = ramp(min_rec, 0.2, 0.5)      # a deliverable without real calls earns nothing

    ctx.component("pairing", (0.0 if wrong else 1.0) * base, 0.20,
                  "every patient's (tumour, normal) pair equals the planted one; x ramp(min recall, 0.2, 0.5)")
    ctx.component("somatic_precision", ramp(min_prec, 0.85, 0.95) * base, 0.30,
                  "min over patients of TP/(TP+FP); full credit >= 0.95, none <= 0.85")
    ctx.component("somatic_recall", ramp(min_rec, 0.85, 0.95) * ramp(min_prec, 0.50, 0.85), 0.30,
                  "min over patients of TP/scored; full credit >= 0.95, none <= 0.85; x ramp(min precision, 0.50, 0.85)")
    ctx.component("allele_fraction", ramp(worst_af, 0.08, 0.03) * base, 0.10,
                  "worst patient's median |AF - realised fragment fraction| over true positives; full <= 0.03, none >= 0.08")

    # ---------------- summary.json -----------------------------------------------------
    summ = ctx.json("summary.json")
    ctx.gate(isinstance(summ, dict) and isinstance(summ.get("patients"), dict), "summary.json has a patients object")
    ok = 0
    for p in pats:
        e = summ["patients"].get(p)
        if not isinstance(e, dict):
            continue
        try:
            good = (int(e.get("n_snv")) == counts[p][0] and int(e.get("n_indel")) == counts[p][1]
                    and not isinstance(e.get("n_snv"), bool)
                    and str(e.get("tumour_sample")) == pairing[p][0] and str(e.get("normal_sample")) == pairing[p][1])
        except (TypeError, ValueError, OverflowError):
            good = False
        ok += 1 if good else 0
    ctx.component("summary", ok / len(pats) * base, 0.10,
                  "per patient: n_snv / n_indel equal the PASS records of the VCF and the sample ids repeat pairing.tsv")


run(main, pass_threshold=0.9)
