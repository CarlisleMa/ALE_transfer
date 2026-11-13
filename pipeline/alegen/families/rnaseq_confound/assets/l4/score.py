"""Verifier for level 4 of the rnaseq_confound family (taking over an RNA-seq project).

reference/truth.json holds only planted quantities (fixed by the generator before a read existed):
  true_counts[gene][sample]  reads the simulator emitted for that gene into that FASTQ
  lfc[gene]                  planted treatment log2 fold change (0 = does not respond)
  sex / condition [sample]   donor sex and treatment of the library inside the sample's FASTQ
  swapped                    samples whose FASTQ holds another sample's library
  mixed                      samples whose FASTQ holds a mixture of two donors' libraries
  unscored_genes             sex-chromosome control genes (DE status depends on a modelling choice)
  audit[item]                planted state of every checklist item of the SOP (ok / defective)
  mailed[gene]               planted status of every mailed gene (confirmed = responds, else retracted)
"""
import math
import sys

sys.path.insert(0, "/tests")
from _scoring import CandidateError, ramp, run, to_float  # noqa: E402

# bands of the level-3 components are unchanged (docs/families/rnaseq_confound.md)
COUNT_AGREE_FULL, COUNT_AGREE_ZERO = 0.90, 0.50
FDP_FULL, FDP_ZERO = 0.12, 0.30
SENS_FULL, SENS_ZERO = 0.85, 0.50
LFC_FULL, LFC_ZERO = 0.33, 0.55
# level-4 components (measured values: "Level 4" section of the design notes)
AUDIT_FULL, AUDIT_ZERO = 1.0, 0.70      # one wrong verdict gives MCC 0.745-0.775 on 8 items with 3-5 defects
ERRATUM_FULL, ERRATUM_ZERO = 0.80, 0.40
UNRELEASED_FACTOR = 0.5                 # DE credit when the sample table shows an unacceptable model
NA = {"", "na", "nan", "null", "none"}
PADJ_ROUNDING = 0.005

W = {"gene_counts": 0.08, "sample_sex": 0.02, "sample_labels": 0.07, "de_false_discoveries": 0.19,
     "de_sensitivity": 0.16, "de_log2fc": 0.13, "audit": 0.20, "erratum": 0.15}


def cell(row, key):
    v = row.get(key)
    return v.strip() if isinstance(v, str) else ""


def mcc(pairs):
    """Matthews correlation of (truth_positive, called_positive) pairs; 0 when undefined."""
    tp = sum(1 for t, c in pairs if t and c)
    tn = sum(1 for t, c in pairs if not t and not c)
    fp = sum(1 for t, c in pairs if not t and c)
    fn = sum(1 for t, c in pairs if t and not c)
    den = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return ((tp * tn - fp * fn) / den if den else 0.0), (tp, fp, fn, tn)


def main(ctx):
    truth = ctx.ref_json("truth.json")
    samples = truth["samples"]
    genes = truth["genes"]
    tcounts = truth["true_counts"]
    lfc = truth["lfc"]
    unscored = set(truth["unscored_genes"])
    swapped = set(truth["swapped"])
    mixed = set(truth["mixed"])
    gset = set(genes)

    # ---------------- gene_counts.tsv --------------------------------------------------
    rows = ctx.table("gene_counts.tsv", ["gene_id"] + samples, max_rows=len(genes) + 10)
    seen = {}
    for r in rows:
        g = cell(r, "gene_id")
        ctx.gate(g in gset, f"gene_counts.tsv: gene id {g[:30]!r} is in tx2gene.tsv")
        ctx.gate(g not in seen, f"gene_counts.tsv: gene {g[:30]!r} listed once")
        vals = []
        for s in samples:
            v = to_float(cell(r, s), f"gene_counts.tsv {g} {s}")
            if v < 0:
                raise CandidateError("gene_counts.tsv: counts must be non-negative")
            vals.append(v)
        seen[g] = vals
    ctx.gate(len(seen) == len(genes), "gene_counts.tsv lists every gene exactly once")
    per_sample = []
    for j, s in enumerate(samples):
        inter = sum(min(seen[g][j], tcounts[g][j]) for g in genes)
        tot_c = sum(seen[g][j] for g in genes)
        tot_t = sum(tcounts[g][j] for g in genes)
        per_sample.append(inter / max(tot_c, tot_t, 1))
        ctx.info[f"count_agreement_{s}"] = round(per_sample[-1], 4)
    ctx.info["count_agreement_min"] = round(min(per_sample), 4)

    # A column has to hold the reads of its own FASTQ file: it must agree with the planted counts of
    # that file better than with those of any other file. (The quantification of another sample's
    # file agrees 0.70-0.72 with the truth of the column it sits in, which would otherwise still
    # earn half of the count credit.)
    def agreement(j, k):
        inter = sum(min(seen[g][j], tcounts[g][k]) for g in genes)
        return inter / max(sum(seen[g][j] for g in genes), sum(tcounts[g][k] for g in genes), 1)

    foreign = [s for j, s in enumerate(samples)
               if max(range(len(samples)), key=lambda k: agreement(j, k)) != j]
    ctx.info["count_columns_matching_another_file"] = foreign
    counts_credit = ctx.component(
        "gene_counts", 0.0 if foreign else ramp(min(per_sample), COUNT_AGREE_ZERO, COUNT_AGREE_FULL),
        W["gene_counts"],
        "min over samples of sum(min(obs,true))/max(sum obs,sum true); full credit >= 0.90; none when a "
        "column agrees better with the planted counts of another FASTQ file than with its own")
    quantified = 1.0 if counts_credit > 0 else 0.0

    # ---------------- sample_qc.tsv ----------------------------------------------------
    qrows = ctx.table("sample_qc.tsv", ["sample_id", "inferred_sex", "condition_used"], max_rows=len(samples) + 10)
    qc = {}
    for r in qrows:
        s = cell(r, "sample_id")
        ctx.gate(s in truth["sex"], f"sample_qc.tsv: sample id {s[:30]!r} is in samples.tsv")
        ctx.gate(s not in qc, f"sample_qc.tsv: sample {s[:30]!r} listed once")
        sex = cell(r, "inferred_sex").upper()
        sex = {"FEMALE": "F", "MALE": "M"}.get(sex, sex)
        cond = cell(r, "condition_used").lower()
        if cond not in ("treated", "vehicle", "excluded"):
            raise CandidateError("sample_qc.tsv: condition_used must be treated, vehicle or excluded")
        qc[s] = (sex, cond)
    ctx.gate(len(qc) == len(samples), "sample_qc.tsv lists every sample exactly once")
    sex_ok = sum(1 for s in samples if qc[s][0] == ("MIXED" if s in mixed else truth["sex"][s]))
    ctx.info["sex_correct"] = f"{sex_ok}/{len(samples)}"
    ctx.component("sample_sex", quantified if sex_ok == len(samples) else 0.0, W["sample_sex"],
                  "inferred_sex is the sex of the library inside the FASTQ, `mixed` for the library that fails "
                  "the identity control (credited only together with gene_counts)")
    lab_ok = 0
    for s in samples:
        if s in mixed:                      # two donors in one tube: nothing but exclusion is valid
            ok = qc[s][1] == "excluded"
        else:
            ok = qc[s][1] == truth["condition"][s] or (s in swapped and qc[s][1] == "excluded")
        lab_ok += ok
    labels_valid = lab_ok == len(samples)
    ctx.info["labels_valid"] = f"{lab_ok}/{len(samples)}"
    ctx.component("sample_labels", quantified if labels_valid else 0.0, W["sample_labels"],
                  "condition_used: mixed library excluded; every other library under its true condition "
                  "(the two mixed-up libraries may instead be excluded); nobody else dropped")
    # SOP section 2: a table from a model with an unaccepted or wrongly labelled library is not released
    released = 1.0 if labels_valid else UNRELEASED_FACTOR
    ctx.info["de_release_factor"] = released

    # ---------------- de_results.tsv ---------------------------------------------------
    drows = ctx.table("de_results.tsv", ["gene_id", "log2fc", "padj", "significant"], max_rows=len(genes) + 10)
    de = {}
    for r in drows:
        g = cell(r, "gene_id")
        ctx.gate(g in gset, f"de_results.tsv: gene id {g[:30]!r} is in tx2gene.tsv")
        ctx.gate(g not in de, f"de_results.tsv: gene {g[:30]!r} listed once")
        sig = cell(r, "significant")
        if sig not in ("0", "1"):
            raise CandidateError("de_results.tsv: significant must be 0 or 1")
        padj = cell(r, "padj")
        if padj.lower() in NA:
            if sig == "1":
                raise CandidateError("de_results.tsv: a gene without padj cannot be significant")
        else:
            p = to_float(padj, f"padj of {g}")
            if not 0.0 <= p <= 1.0:
                raise CandidateError("de_results.tsv: padj must lie in [0, 1]")
            if (sig == "1" and p >= 0.05 + PADJ_ROUNDING) or (sig == "0" and p < 0.05 - PADJ_ROUNDING):
                raise CandidateError("de_results.tsv: significant must be 1 exactly when padj < 0.05")
        de[g] = (to_float(cell(r, "log2fc"), f"log2fc of {g}"), int(sig))
    ctx.gate(len(de) == len(genes), "de_results.tsv lists every gene exactly once")

    scored = [g for g in genes if g not in unscored]
    planted = [g for g in scored if lfc[g] != 0.0]
    tp = sum(1 for g in planted if de[g][1] == 1)
    fp = sum(1 for g in scored if lfc[g] == 0.0 and de[g][1] == 1)
    wrong_sign = sum(1 for g in planted if de[g][1] == 1 and de[g][0] * lfc[g] <= 0)
    tp -= wrong_sign
    fp += wrong_sign
    sens = tp / len(planted)
    fdp = fp / (tp + fp) if tp + fp else 1.0
    errs = sorted(abs(de[g][0] - lfc[g]) for g in planted)
    med_err = errs[len(errs) // 2]
    ctx.info.update({"de_planted": len(planted), "de_tp": tp, "de_fp": fp, "de_sensitivity": round(sens, 4),
                     "de_fdp": round(fdp, 4), "lfc_median_abs_error": round(med_err, 4)})
    ctx.component("de_false_discoveries", released * ramp(fdp, FDP_ZERO, FDP_FULL) * ramp(sens, 0.0, 0.5),
                  W["de_false_discoveries"],
                  "false-discovery proportion among significant genes vs planted truth; full <= 0.12, none >= 0.30; "
                  "scaled down when fewer than half of the planted genes are found; halved when the sample table "
                  "is not valid")
    ctx.component("de_sensitivity", released * ramp(sens, SENS_ZERO, SENS_FULL) * ramp(fdp, 0.60, 0.30),
                  W["de_sensitivity"],
                  "fraction of planted DE genes called with the right sign; full >= 0.85, none <= 0.50; scaled "
                  "down when FDP exceeds 0.30; halved when the sample table is not valid")
    ctx.component("de_log2fc", released * ramp(med_err, LFC_ZERO, LFC_FULL), W["de_log2fc"],
                  "median |log2fc - planted log2FC| over planted DE genes; full <= 0.33, none >= 0.55; halved "
                  "when the sample table is not valid")

    # ---------------- audit.tsv --------------------------------------------------------
    items = truth["audit"]
    arows = ctx.table("audit.tsv", ["item", "verdict"], max_rows=len(items) + 10)
    got = {}
    for r in arows:
        it = cell(r, "item")
        ctx.gate(it in items, f"audit.tsv: item {it[:40]!r} is a checklist item of the SOP")
        ctx.gate(it not in got, f"audit.tsv: item {it[:40]!r} listed once")
        v = cell(r, "verdict").lower()
        if v not in ("ok", "defective"):
            raise CandidateError("audit.tsv: verdict must be ok or defective")
        got[it] = v
    ctx.gate(len(got) == len(items), "audit.tsv gives a verdict for every checklist item")
    a_mcc, (atp, afp, afn, atn) = mcc([(items[i] == "defective", got[i] == "defective") for i in items])
    ctx.info.update({"audit_mcc": round(a_mcc, 4), "audit_missed": afn, "audit_false_alarms": afp,
                     "audit_planted_defects": atp + afn})
    ctx.component("audit", ramp(a_mcc, AUDIT_ZERO, AUDIT_FULL), W["audit"],
                  "MCC of the verdicts against the planted defects; full credit only for a complete audit, "
                  "none at MCC <= 0.70 (two or more wrong verdicts)")

    # ---------------- erratum.tsv ------------------------------------------------------
    mailed = truth["mailed"]
    erows = ctx.table("erratum.tsv", ["gene_id", "verdict"], max_rows=len(mailed) + 10)
    eg = {}
    for r in erows:
        g = cell(r, "gene_id")
        ctx.gate(g in mailed, f"erratum.tsv: gene {g[:30]!r} is on the mailed list")
        ctx.gate(g not in eg, f"erratum.tsv: gene {g[:30]!r} listed once")
        v = cell(r, "verdict").lower()
        if v not in ("confirmed", "retracted"):
            raise CandidateError("erratum.tsv: verdict must be confirmed or retracted")
        eg[g] = v
    ctx.gate(len(eg) == len(mailed), "erratum.tsv lists every mailed gene exactly once")
    e_mcc, (etp, efp, efn, etn) = mcc([(mailed[g] == "confirmed", eg[g] == "confirmed") for g in mailed])
    ctx.info.update({"erratum_mcc": round(e_mcc, 4), "erratum_true_responders": etp + efn,
                     "erratum_wrongly_confirmed": efp, "erratum_wrongly_retracted": efn})
    ctx.component("erratum", ramp(e_mcc, ERRATUM_ZERO, ERRATUM_FULL), W["erratum"],
                  "MCC of confirmed/retracted against the planted responder status of the mailed genes; "
                  "full >= 0.80, none <= 0.40")


run(main, pass_threshold=0.9)
