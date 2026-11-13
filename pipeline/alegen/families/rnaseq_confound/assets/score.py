"""Verifier for the rnaseq_confound family (bulk RNA-seq, FASTQ -> differential expression).

Everything in reference/truth.json was fixed by the simulator before a read existed:
  true_counts[gene][sample]  number of reads the simulator emitted for that gene into that FASTQ
  lfc[gene]                  planted treatment log2 fold change (0 = not differentially expressed)
  sex[sample]                donor sex of the library that is actually inside the sample's FASTQ
  condition[sample]          treatment of the library that is actually inside the sample's FASTQ
  swapped                    samples whose FASTQ holds another sample's library
  unscored_genes             sex-specific genes (their DE status depends on a modelling choice)
No analysis tool produced any of it, so tolerances only have to cover sampling noise and
legitimate method choices (see the family design notes for the measured values).
"""
import sys

sys.path.insert(0, "/tests")
from _scoring import CandidateError, ramp, run, to_float  # noqa: E402

# measured bands (docs/families/rnaseq_confound.md, "Scoring" section)
COUNT_AGREE_FULL, COUNT_AGREE_ZERO = 0.90, 0.50
FDP_FULL, FDP_ZERO = 0.12, 0.30
SENS_FULL, SENS_ZERO = 0.85, 0.50
LFC_FULL, LFC_ZERO = 0.33, 0.55
NA = {"", "na", "nan", "null", "none"}
PADJ_ROUNDING = 0.005     # a padj written with 2-3 decimals may sit on the other side of 0.05


def cell(row, key):
    """Field of a parsed row as a stripped string ('' when the row is short)."""
    v = row.get(key)
    return v.strip() if isinstance(v, str) else ""


def main(ctx):
    truth = ctx.ref_json("truth.json")
    samples = truth["samples"]
    genes = truth["genes"]
    tcounts = truth["true_counts"]
    lfc = truth["lfc"]
    unscored = set(truth["unscored_genes"])
    swapped = set(truth["swapped"])
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
    counts_credit = ctx.component(
        "gene_counts", ramp(min(per_sample), COUNT_AGREE_ZERO, COUNT_AGREE_FULL), 0.15,
        "min over samples of sum(min(obs,true))/max(sum obs,sum true); full credit >= 0.90")
    # the identity check has to come from the data: without a usable quantification the QC table
    # can only be a copy of the sample sheet, which must not earn credit
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
    sex_ok = sum(1 for s in samples if qc[s][0] == truth["sex"][s])
    ctx.info["sex_correct"] = f"{sex_ok}/{len(samples)}"
    ctx.component("sample_sex", quantified if sex_ok == len(samples) else 0.0, 0.04,
                  "inferred_sex equals the sex of the library inside the FASTQ for every sample "
                  "(credited only together with gene_counts)")
    # a sample whose FASTQ holds another library may be relabelled or excluded; no other sample may be dropped
    lab_ok = 0
    for s in samples:
        ok = qc[s][1] == truth["condition"][s] or (s in swapped and qc[s][1] == "excluded")
        lab_ok += ok
    ctx.info["labels_valid"] = f"{lab_ok}/{len(samples)}"
    # weight 0.11 (not 0.10): a wrong sample table with otherwise perfect results must end below 0.9
    ctx.component("sample_labels", quantified if lab_ok == len(samples) else 0.0, 0.11,
                  "condition_used is the true condition of every library; mixed-up samples may instead be "
                  "excluded (credited only together with gene_counts)")

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
    tp -= wrong_sign          # significant in the wrong direction is a false discovery, not a hit
    fp += wrong_sign
    sens = tp / len(planted)
    fdp = fp / (tp + fp) if tp + fp else 1.0
    errs = sorted(abs(de[g][0] - lfc[g]) for g in planted)
    med_err = errs[len(errs) // 2]
    ctx.info.update({"de_planted": len(planted), "de_tp": tp, "de_fp": fp, "de_sensitivity": round(sens, 4),
                     "de_fdp": round(fdp, 4), "lfc_median_abs_error": round(med_err, 4)})
    # each rate is guarded against the trivial answer that maximises it (nothing / everything significant)
    ctx.component("de_false_discoveries", ramp(fdp, FDP_ZERO, FDP_FULL) * ramp(sens, 0.0, 0.5), 0.30,
                  "false-discovery proportion among significant genes vs planted truth; full credit <= 0.12, "
                  "none >= 0.30; scaled down when fewer than half of the planted genes are found")
    ctx.component("de_sensitivity", ramp(sens, SENS_ZERO, SENS_FULL) * ramp(fdp, 0.60, 0.30), 0.25,
                  "fraction of planted DE genes called significant with the right sign; full credit >= 0.85, "
                  "none <= 0.50; scaled down when the false-discovery proportion exceeds 0.30")
    ctx.component("de_log2fc", ramp(med_err, LFC_ZERO, LFC_FULL), 0.15,
                  "median |log2fc - planted log2FC| over planted DE genes; full credit <= 0.33, none >= 0.55")


run(main, pass_threshold=0.9)
