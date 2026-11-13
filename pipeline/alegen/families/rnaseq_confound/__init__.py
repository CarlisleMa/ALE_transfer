"""Family: bulk RNA-seq, FASTQ -> gene counts -> sample identity check -> differential expression.

Truth by construction: the simulator fixes the treatment log2 fold change of every gene, the
library-prep batch effect, the donor sex and the number of reads per gene and library, and only
then writes reads sampled from real Ensembl transcript sequences. Which genes respond is random,
so real gene names carry no information about the answer.

Decision points (level 2 adds the first, level 3 the second):
  batch_confound   two library-prep batches, 5+1 treated vs 1+5 vehicle; 30% of genes carry a batch
                   effect. A condition-only model reports batch genes as treatment genes.
  sample_swap      the FASTQ files of one treated and one vehicle library (opposite donor sex,
                   different batches) are filed under each other's sample id.
"""
from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np

from ...spec import Attack, Probe, TaskBuild, Trap
from ... import dockerlib
from . import image, sim

ASSETS = Path(__file__).parent / "assets"
N_GENES = 800
TARGET_READS = 93_000         # mean expected reads per library; 19.7 gzip bytes per read -> 21.5-23.5 MB of FASTQ
ERR = 0.003                   # per-base substitution rate
PASS = 0.9

INSTRUCTION = """\
You are the bioinformatician of a pharmacology group. The wet lab treated donor-derived primary
cells with a candidate compound or with vehicle, prepared poly-A RNA-seq libraries and now needs
gene-level counts, a per-sample QC table and the list of genes that respond to the compound.

Task directory:
- `/task`

Inputs (read-only) under `/task/input`:
- `fastq/<sample_id>.fastq.gz` - one single-end, unstranded FASTQ per sample ({n_samples} samples)
- `samples.tsv` - sample sheet as filled in by the wet lab (`sample_id`, `condition` = `treated` or
  `vehicle`, `donor_id`, `donor_sex` = `F` or `M`, `library_prep_batch`, `fastq_file`)
- `reference/transcripts.fa.gz` - transcript sequences of the gene panel that was sequenced
- `reference/tx2gene.tsv` - `transcript_id`, `gene_id`, `gene_name` for every transcript
- `lab_notes.md` - wet-lab notes for this experiment

Available software (on `PATH`): `salmon`, Python 3.12 with numpy, scipy, pandas, statsmodels,
scikit-learn and pydeseq2. There is no internet access and no R installation.

Your task:
1. Quantify every sample against the provided transcripts and summarise to genes.
2. Check every sample and decide how it enters the statistical analysis.
3. Test every gene for differential expression between compound-treated and vehicle-treated cells.

Required outputs under `/task/output`:
- `gene_counts.tsv` - tab-separated, header `gene_id` followed by one column per `sample_id`
  (all {n_samples} samples, named as in `samples.tsv`). One row per `gene_id` of `tx2gene.tsv`.
  Values are the estimated number of reads of that FASTQ file that originate from the gene
  (non-negative, not length- or depth-normalised; fractional values are allowed).
- `sample_qc.tsv` - tab-separated with at least the columns `sample_id`, `inferred_sex`,
  `condition_used`; one row per sample of `samples.tsv`.
  - `inferred_sex`: `F` or `M`, the sex of the donor of the library contained in that sample's
    FASTQ file, as determined from the sequencing data.
  - `condition_used`: `treated` or `vehicle` - the condition under which the reads of that FASTQ
    file entered your differential-expression model - or `excluded` if you left the sample out.
- `de_results.tsv` - tab-separated with at least the columns `gene_id`, `log2fc`, `padj`,
  `significant`; one row per `gene_id` of `tx2gene.tsv`.
  - `log2fc`: estimated log2 fold change of the gene's expression caused by the compound (treated
    over vehicle), on the scale on which a gene that does not respond has `log2fc` = 0. Must be a
    finite number for every gene (use 0 for a gene that cannot be estimated).
  - `padj`: multiple-testing-adjusted p-value of the test "the compound does not change the
    expression of this gene"; `NA` for a gene that was not tested.
  - `significant`: `1` if the gene is differentially expressed at a false-discovery rate of 5%
    (`padj` < 0.05, no fold-change cutoff), otherwise `0`.

Constraints:
- Do not modify anything under `/task/input`; write deliverables only under `/task/output`.
- Extra columns are allowed in every output file.
- Do not ask for confirmation. Execute directly.
"""

JUDGE = """\
**Where the truth comes from.** Before any read exists the generator fixes (a) which genes respond
to the compound and their log2 fold change, (b) the library-prep batch effect of every gene, (c) the
sex and treatment of every library and which FASTQ file it is written to, and (d) the number of
reads each gene contributes to each library (negative binomial around the planted mean).
`reference/truth.json` stores exactly these planted values. No analysis tool is run to obtain them.

**Hard gates (score 0).** Missing/unparseable file; a gene or sample missing or duplicated in any
table; negative counts; `condition_used` outside {treated, vehicle, excluded}; `significant` not in
{0,1} or inconsistent with `padj` (`1` requires `padj` < 0.05 and vice versa, with 0.005 slack for
rounded values; `NA` requires `0`); non-finite `log2fc`.

**Components.**

| weight | component | metric | full credit | zero |
|---|---|---|---|---|
| 0.15 | gene_counts | per sample `sum(min(obs,true)) / max(sum obs, sum true)`, **minimum** over samples | >= 0.90 | <= 0.50 |
| 0.04 | sample_sex | `inferred_sex` equals the planted sex of the library in the FASTQ, all samples; credited only if gene_counts > 0 | all | any wrong |
| 0.11 | sample_labels | `condition_used` equals the planted condition for every sample; the mixed-up samples may instead be `excluded`; no other sample may be dropped; credited only if gene_counts > 0 | all | any wrong |
| 0.30 | de_false_discoveries | false-discovery proportion among `significant` genes against the planted DE genes (a planted gene called with the wrong sign counts as false); multiplied by ramp(sensitivity, 0, 0.5) | <= 0.12 | >= 0.30 |
| 0.25 | de_sensitivity | fraction of planted DE genes called with the right sign; multiplied by ramp(FDP, 0.60, 0.30) | >= 0.85 | <= 0.50 |
| 0.15 | de_log2fc | median absolute error of `log2fc` over planted DE genes | <= 0.33 | >= 0.55 |

The six sex-specific genes are not scored in the DE components (whether they are "differentially
expressed" depends on whether donor sex is modelled, a legitimate analyst choice).

**Why these tolerances.** See `docs/families/rnaseq_confound.md` for the measured values of the
oracle (salmon + pyDESeq2, records exchanged), the independent solution (k-mer counter + voom-style
linear model, samples excluded), real DESeq2 / limma-voom / edgeR on the same counts, and every
trap-ignoring analysis on two seeds. The bands sit between the two groups with margin on both sides.

**What is deliberately not scored.** The numeric value of `padj` (method dependent), `log2fc` of
genes without a planted effect, the DE status of sex-specific genes, the batch the analyst used.
"""


def _lab_notes(level: int, libs: list[dict], batch_names: dict[int, str]) -> str:
    n = len(libs)
    lines = [
        "# Wet-lab notes: compound vs vehicle, primary cells",
        "",
        "- Cells: primary cells from {} unrelated donors (one culture per donor), passage 2.".format(n),
        "- Treatment: compound at 1 uM or vehicle (0.1% DMSO) for 24 h; one condition per donor.",
        "- RNA: column purification with on-column DNase; RIN >= 8.9 for all samples.",
        "- Gene panel: the transcripts in `reference/` (targeted panel plus sex-chromosome control genes).",
    ]
    if level >= 2:
        d0, d1 = batch_names[0], batch_names[1]
        n0 = sum(1 for d in libs if d["batch"] == 0)
        lines += [
            f"- Library prep: poly-A selection, unstranded. Two prep sessions because the first kit lot ran out: "
            f"{d0} ({n0} libraries, kit lot A, operator 1) and {d1} ({n - n0} libraries, kit lot B, operator 2, "
            f"different PCR cycler). The session of each library is in the sample sheet.",
        ]
    else:
        lines += [f"- Library prep: poly-A selection, unstranded; all {n} libraries in one session "
                  f"({batch_names[0]}), same kit lot and operator."]
    lines += [
        "- Sequencing: all libraries pooled equimolar on one lane, 50-cycle single read. The core delivered "
        "adapter-trimmed reads; read names and base qualities were stripped to placeholders by their export tool.",
    ]
    if level >= 3:
        lines += [
            "- Incident (recorded for completeness): the final libraries were stored at -20 C in a cardboard box. "
            "When the box was taken out for pooling, condensation had made the hand-written tube labels "
            "partly illegible. All tubes were re-labelled from the box map before pooling and index assignment. "
            "Library-prep session and all other sample-sheet fields were copied from the original records.",
        ]
    else:
        lines += ["- No deviations from the protocol were recorded."]
    return "\n".join(lines) + "\n"


def build_instance(seed: int, level: int, workdir: Path | None, write_reads: bool = True,
                   target_reads: int = TARGET_READS, major: int = 5, **effect_kw) -> dict:
    """Everything except the TaskBuild wrapping (also used by the calibration harness)."""
    # level 4 is a separate instance (other genes, effects, design), not the level-3 task of the same seed
    rng = np.random.default_rng((52_000 if level <= 3 else 64_000) + seed)
    genes = sim.select_genes(rng, N_GENES - 6)
    design = sim.make_design(rng, level, major=major)
    libs = design["libs"]
    eff = sim.plant_effects(rng, genes, level, libs, target_reads, **effect_kw)
    counts = sim.draw_counts(rng, genes, eff, libs)          # genes x libs
    mixture = None
    if level >= 4:                                           # own generator: levels 1-3 keep their random stream
        from . import l4
        rng4 = np.random.default_rng(91_000 + seed)
        mixture = l4.plant_mixture(rng4, genes, libs, design["swapped"], counts)
        counts = l4.apply_mixture(rng4, counts, libs, mixture)

    # the sample sheet is what the wet lab believes; the FASTQ named <file_id> holds the library
    dates = ["2024-03-05", "2024-03-19"] if rng.random() < 0.5 else ["2024-05-14", "2024-06-04"]
    batch_names = {0: dates[0], 1: dates[1]}
    sheet = sorted(libs, key=lambda d: d["sheet_id"])
    sids = [d["sheet_id"] for d in sheet]
    by_file = {d["file_id"]: d for d in libs}                # sample id -> library inside its FASTQ
    lib_col_of_file = {sid: libs.index(by_file[sid]) for sid in sids}
    file_counts = np.column_stack([counts[:, lib_col_of_file[s]] for s in sids])   # genes x samples (files)

    gene_ids = [g["gene_id"] for g in genes]
    samples_tsv = "sample_id\tcondition\tdonor_id\tdonor_sex\tlibrary_prep_batch\tfastq_file\n" + "".join(
        f"{d['sheet_id']}\t{d['condition']}\t{d['donor']}\t{d['sex']}\t{batch_names[d['batch']]}\t"
        f"{d['sheet_id']}.fastq.gz\n" for d in sheet)
    tx2gene = "transcript_id\tgene_id\tgene_name\n" + "".join(
        f"{tid}\t{g['gene_id']}\t{g['symbol']}\n" for g in genes for tid, _ in g["tx"])
    truth = {
        "samples": sids, "genes": gene_ids,
        "true_counts": {g: [int(v) for v in file_counts[i]] for i, g in enumerate(gene_ids)},
        "lfc": {g: float(eff["lfc"][i]) for i, g in enumerate(gene_ids)},
        "batch_lfc": {g: float(eff["batch_lfc"][i]) for i, g in enumerate(gene_ids)},
        "unscored_genes": [g["gene_id"] for g in genes if g["kind"] == "sex"],
        "sex": {s: by_file[s]["sex"] for s in sids},
        "condition": {s: by_file[s]["condition"] for s in sids},
        "batch": {s: batch_names[by_file[s]["batch"]] for s in sids},
        "swapped": design["swapped"],
    }
    inst = {"genes": genes, "libs": libs, "sheet": sheet, "sids": sids, "eff": eff, "truth": truth,
            "file_counts": file_counts, "samples_tsv": samples_tsv, "tx2gene": tx2gene,
            "lab_notes": _lab_notes(level, libs, batch_names), "batch_names": batch_names, "fastq": {}}
    if level >= 4:
        inst["mixture"] = mixture
        inst["rng4"] = rng4
        truth["mixed"] = [mixture["file_id"]]
        truth["mixture"] = {k: mixture[k] for k in ("file_id", "contaminant_file_id", "fraction")}
    if workdir is None:
        return inst
    workdir = Path(workdir)
    (workdir / "fastq").mkdir(parents=True, exist_ok=True)
    fa = workdir / "transcripts.fa.gz"
    with open(fa, "wb") as fh, gzip.GzipFile(filename="", mode="wb", fileobj=fh, compresslevel=9, mtime=0) as gz:
        for g in genes:
            for tid, s in g["tx"]:
                gz.write(f">{tid}\n".encode())
                gz.write(("\n".join(s[i:i + 80] for i in range(0, len(s), 80)) + "\n").encode())
    inst["transcripts"] = fa
    if write_reads:
        tr = sim.Transcriptome(genes, rng)
        for j, s in enumerate(sids):
            p = workdir / "fastq" / f"{s}.fastq.gz"
            sim.write_fastq(p, tr.reads(file_counts[:, j], rng, ERR))
            inst["fastq"][s] = p
    return inst


def _tsv_counts(gene_ids, sids, mat) -> str:
    return "gene_id\t" + "\t".join(sids) + "\n" + "".join(
        f"{g}\t" + "\t".join(str(int(v)) for v in mat[i]) + "\n" for i, g in enumerate(gene_ids))


def _tsv_qc(sids, sex, cond) -> str:
    return "sample_id\tinferred_sex\tcondition_used\n" + "".join(f"{s}\t{sex[s]}\t{cond[s]}\n" for s in sids)


def _tsv_de(gene_ids, lfc, sig) -> str:
    return "gene_id\tlog2fc\tpadj\tsignificant\n" + "".join(
        f"{g}\t{lfc[g]:.4f}\t{'0.001' if sig[g] else '0.5'}\t{int(sig[g])}\n" for g in gene_ids)


def generate(seed: int, level: int, workdir: Path) -> TaskBuild:
    if level >= 4:
        from . import l4
        return l4.generate(seed, level, workdir)
    inst = build_instance(seed, level, workdir)
    truth, sids, gene_ids = inst["truth"], inst["sids"], inst["truth"]["genes"]
    fc = inst["file_counts"]
    lfc = truth["lfc"]
    unscored = set(truth["unscored_genes"])
    planted = [g for g in gene_ids if lfc[g] != 0 and g not in unscored]
    null = [g for g in gene_ids if lfc[g] == 0 and g not in unscored]
    sheet_cond = {d["sheet_id"]: d["condition"] for d in inst["sheet"]}
    sheet_sex = {d["sheet_id"]: d["sex"] for d in inst["sheet"]}

    prng = np.random.default_rng(seed + 977)
    sig_true = {g: lfc[g] != 0 and g not in unscored for g in gene_ids}
    # inside tolerance: 6% of reads lost, 4 false positives (FDP ~0.05), 8% of planted genes missed,
    # log2fc noise sd 0.2 (median abs error ~0.13)
    thinned = prng.binomial(fc, 0.94)
    fp_in = set(prng.permutation(null)[:4].tolist())
    miss_in = set(prng.permutation(planted)[:int(0.08 * len(planted))].tolist())
    sig_in = {g: (sig_true[g] and g not in miss_in) or g in fp_in for g in gene_ids}
    lfc_in = {g: lfc[g] + float(prng.normal(0, 0.2)) for g in gene_ids}
    # just beyond tolerance: 40% false discoveries, everything else perfect
    n_fp_out = int(round(len(planted) * 0.40 / 0.60))
    fp_out = set(prng.permutation(null)[:n_fp_out].tolist())
    sig_out = {g: sig_true[g] or g in fp_out for g in gene_ids}
    lfc_out = {g: (lfc[g] if g not in fp_out else 0.8) for g in gene_ids}
    # just beyond tolerance: every planted log2fc attenuated by one third (what a 2/12 label error does)
    lfc_att = {g: lfc[g] * (2.0 / 3.0) for g in gene_ids}

    truth_files = {"gene_counts.tsv": _tsv_counts(gene_ids, sids, fc),
                   "sample_qc.tsv": _tsv_qc(sids, truth["sex"], truth["condition"]),
                   "de_results.tsv": _tsv_de(gene_ids, lfc, sig_true)}
    pipe = {"rnaseq_pipeline.py": ASSETS / "rnaseq_pipeline.py"}
    probes = [
        Probe("planted_truth", "submission written directly from the planted truth", 0.999, files=truth_files),
        Probe("within_tolerance", "6% of reads lost, 4 false discoveries, 8% of planted genes missed, "
              "log2fc noise sd 0.2", 0.99,
              files={"gene_counts.tsv": _tsv_counts(gene_ids, sids, thinned),
                     "sample_qc.tsv": truth_files["sample_qc.tsv"],
                     "de_results.tsv": _tsv_de(gene_ids, lfc_in, sig_in)}),
        Probe("independent_solution", "second tool chain: k-mer read counter (no salmon), voom-style weighted "
              "linear model (no DESeq2), mismatching samples excluded instead of relabelled",
              0.99, script="SOLUTION=2 bash /solution/solve.sh"),
    ]
    if level >= 3:
        probes.append(Probe("oracle_exclude_resolution", "oracle tool chain, but the two mixed-up samples are "
                            "excluded instead of relabelled", 0.99,
                            script="python3 /alegen_script/rnaseq_pipeline.py --labels exclude", script_files=pipe))
    ones = np.ones_like(fc)
    attacks = [
        Attack("sheet_copied_no_analysis", "prior-only answer: QC table copied from the sample sheet, counts all 1, "
               "no gene significant", 0.10,
               files={"gene_counts.tsv": _tsv_counts(gene_ids, sids, ones),
                      "sample_qc.tsv": _tsv_qc(sids, sheet_sex, sheet_cond),
                      "de_results.tsv": _tsv_de(gene_ids, {g: 0.0 for g in gene_ids}, {g: False for g in gene_ids})}),
        Attack("everything_significant", "true counts and QC, but every gene called significant with log2fc 1", 0.45,
               files={"gene_counts.tsv": truth_files["gene_counts.tsv"],
                      "sample_qc.tsv": truth_files["sample_qc.tsv"],
                      "de_results.tsv": _tsv_de(gene_ids, {g: 1.0 for g in gene_ids}, {g: True for g in gene_ids})}),
        Attack("fdp_beyond_tolerance", "perfect counts, QC and log2fc, but 40% of the significant genes are false",
               PASS - 0.15,
               files={"gene_counts.tsv": truth_files["gene_counts.tsv"],
                      "sample_qc.tsv": truth_files["sample_qc.tsv"],
                      "de_results.tsv": _tsv_de(gene_ids, lfc_out, sig_out)}),
        Attack("lfc_attenuated_beyond_tolerance", "perfect counts, QC and calls, but every log2fc shrunk by a third",
               0.86,
               files={"gene_counts.tsv": truth_files["gene_counts.tsv"],
                      "sample_qc.tsv": truth_files["sample_qc.tsv"],
                      "de_results.tsv": _tsv_de(gene_ids, lfc_att, sig_true)}),
    ]
    traps = []
    if level >= 2:
        attacks.append(Attack(
            "naive_no_batch_term", "correct quantification and sample handling, condition-only design", PASS - 0.15,
            script="python3 /alegen_script/rnaseq_pipeline.py --batch no", script_files=pipe))
        traps.append(Trap(
            "batch_confound",
            "Two library-prep sessions, 5 treated + 1 vehicle in one and 1 treated + 5 vehicle in the other; 30% of "
            "genes differ between sessions by 0.5-2 log2 units. A condition-only model reports batch genes as "
            "treatment genes and biases log2fc by two thirds of the batch effect.",
            "samples.tsv has a library_prep_batch column that is almost aliased with condition; lab_notes.md gives "
            "two prep sessions with different kit lot/operator/cycler; PCA or sample clustering of the counts "
            "separates the sessions, not the conditions",
            "Approach / wrong strategy"))
    if level >= 3:
        attacks += [
            Attack("naive_labels_as_given", "sex inferred correctly but the sample sheet is used as delivered "
                   "(batch term included)", PASS - 0.15,
                   script="python3 /alegen_script/rnaseq_pipeline.py --labels given", script_files=pipe),
            Attack("naive_default_pipeline", "sample sheet as delivered and condition-only design", 0.60,
                   script="python3 /alegen_script/rnaseq_pipeline.py --labels given --batch no", script_files=pipe),
        ]
        a, b = truth["swapped"]
        traps.append(Trap(
            "sample_swap",
            f"The FASTQ files of {a} and {b} hold each other's library (one treated, one vehicle, opposite donor "
            "sex, different prep sessions). Used as labelled, two of twelve samples carry the wrong condition and "
            "batch: effects are attenuated by a third and many true genes lose significance.",
            "the required inferred_sex column makes the analyst type every sample from XIST and chrY genes "
            "(named in tx2gene.tsv); exactly two samples contradict donor_sex, and they are a reciprocal pair; "
            "lab_notes.md records that all library tubes were re-labelled from a box map before pooling; after "
            "exchanging the two records the samples cluster with their recorded prep session",
            "Understanding / unverified assumption about inputs"))

    fastq_files = {f"fastq/{s}.fastq.gz": inst["fastq"][s] for s in sids}
    return TaskBuild(
        task_id=f"rnaseq-confound-s{seed:04d}-l{level}", family="rnaseq_confound", seed=seed, level=level,
        title="Bulk RNA-seq: FASTQ to differential expression under a confounded design",
        summary="Quantify single-end RNA-seq reads of twelve donor samples against a transcript panel, verify "
                "sample identity, and deliver gene counts, sample QC and treated-vs-vehicle differential expression.",
        subdomain="3.2 Genomics & Sequence Analysis",
        software=["salmon", "pyDESeq2", "statsmodels", "Python", "numpy/scipy/pandas"],
        instruction_md=INSTRUCTION.format(n_samples=len(sids)),
        agent_dockerfile=image.agent_dockerfile(),
        verifier_dockerfile=dockerlib.verifier_dockerfile(pip=[]),
        input_files={**fastq_files, "samples.tsv": inst["samples_tsv"],
                     "reference/transcripts.fa.gz": inst["transcripts"],
                     "reference/tx2gene.tsv": inst["tx2gene"], "lab_notes.md": inst["lab_notes"]},
        solution_files={"solve.sh": ASSETS / "solve.sh", "rnaseq_pipeline.py": ASSETS / "rnaseq_pipeline.py",
                        "solve_alt.py": ASSETS / "solve_alt.py"},
        tests_files={"score.py": ASSETS / "score.py", "reference/truth.json": json.dumps(truth)},
        traps=traps, attacks=attacks, probes=probes,
        capabilities=["experimental-design reasoning (aliased batch)", "sample identity QC before statistics",
                      "long tool chain from raw reads to statistics", "exact output contracts"],
        judge_explanation_md=JUDGE, pass_threshold=PASS,
        oracle_steps=15, oracle_runtime_min=0.5, cpus=4, memory_mb=8000, storage_mb=10000,
    )
