"""Level 4 of rnaseq_confound: the level-3 experiment as a project that is taken over (H + Q).

H  inherited project   `/task/input/previous_analysis/` holds the predecessor's scripts, merged count
                       matrix, draft DE table, mailed gene list and hand-over note. They are produced
                       at generation time by running *mutated copies of the reference pipeline* on the
                       task's own FASTQ files inside the family's tool image. A seeded subset (3-5) of
                       the eight checklist items of the SOP is defective; the others are fine, some in
                       an unusual-looking way. Deliverables: the corrected result, `audit.tsv`,
                       `erratum.tsv`.
Q  identity control    one library is a mixture of two donors (male library with 25-45% of a female
                       library of the other condition and session). Its sex call still agrees with
                       the sample sheet; only the SOP's stated criterion on the sex-chromosome
                       control genes (minor marker share < 0.02) identifies it. It cannot be
                       relabelled and has to be excluded.

Everything scored is planted: the defect set, which genes respond, which library is mixed. The tool
run only manufactures the *inputs* (what a predecessor would have left behind).
"""
from __future__ import annotations

import gzip
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from ... import dockerlib
from ...spec import Attack, Probe, TaskBuild, Trap
from . import image, sim

ASSETS = Path(__file__).parent / "assets"
A4 = ASSETS / "l4"
GENTOOLS_TAG = "alegen-rnaseq-confound-gentools:local"
PASS = 0.9
MINOR_SHARE_LIMIT = 0.02      # SOP acceptance criterion (stated to the agent)
N_MAILED = 40

# checklist of the SOP, in SOP order: (item id, kind) ; kind "code" = visible in the scripts,
# "recompute" = the scripts are fine, the shipped file is not what they produce
ITEMS = [
    ("quant_library_type", "code"),
    ("merge_read_counts", "code"),
    ("matrix_sample_files", "recompute"),
    ("matrix_all_reads", "recompute"),
    ("de_session_term", "code"),
    ("de_contrast_direction", "code"),
    ("de_adjusted_pvalues", "code"),
    ("draft_reproducible", "recompute"),
]
ITEM_IDS = [i for i, _ in ITEMS]


# --------------------------------------------------------------------------- Q: the mixed library

def _cell(d: dict) -> tuple:
    return d["condition"], d["batch"]


def plant_mixture(rng, genes: list[dict], libs: list[dict], swapped: list[str], counts: np.ndarray) -> dict:
    """Pick the mixed library c (male, 5-sample cell, not part of the file mix-up), its contaminant d
    (female, other condition, other session, 5-sample cell, not part of the mix-up) and the fraction.

    The fraction is 0.25-0.40, raised when necessary so that the expected minor marker share of the
    mixture is at least 3.5 x the SOP limit; it never exceeds 0.45, where the Y genes still outnumber
    XIST and a plain male/female call agrees with the sample sheet."""
    cells: dict[tuple, list] = {}
    for d in libs:
        cells.setdefault(_cell(d), []).append(d)
    major = [k for k, v in cells.items() if len(v) == max(len(x) for x in cells.values())]
    plain = [d for d in libs if d["sheet_id"] not in swapped]
    cand_c = [d for d in plain if d["sex"] == "M" and _cell(d) in major]
    c = cand_c[int(rng.integers(0, len(cand_c)))]
    cand_d = [d for d in plain if d["sex"] == "F" and _cell(d) in major and d["condition"] != c["condition"]
              and d["batch"] != c["batch"]]
    d = cand_d[int(rng.integers(0, len(cand_d)))]
    ic, id_ = libs.index(c), libs.index(d)
    sym = [g["symbol"] for g in genes]
    xi = sym.index("XIST")
    yi = [i for i, s in enumerate(sym) if s in sim.Y_GENES]
    x_d, y_c = float(counts[xi, id_]), float(counts[yi, ic].sum())
    f = float(rng.uniform(0.25, 0.40))
    target = 3.5 * MINOR_SHARE_LIMIT
    # share = f*x_d / (f*x_d + (1-f)*y_c) >= target
    f_min = target * y_c / (x_d * (1 - target) + target * y_c)
    f = round(min(max(f, f_min + 0.02), 0.45), 3)
    return {"lib": ic, "contaminant_lib": id_, "fraction": f, "file_id": c["file_id"],
            "contaminant_file_id": d["file_id"], "sheet_id": c["sheet_id"]}


def apply_mixture(rng, counts: np.ndarray, libs: list[dict], mix: dict) -> np.ndarray:
    """Reads of the mixed library: a (1-f) subsample of its own fragments plus an f subsample of the
    contaminating library's fragments (same molecules, independently sequenced)."""
    out = counts.copy()
    f = mix["fraction"]
    out[:, mix["lib"]] = rng.binomial(counts[:, mix["lib"]], 1.0 - f) + rng.binomial(counts[:, mix["contaminant_lib"]], f)
    return out


def marker_shares(genes: list[dict], sids: list[str], file_counts: np.ndarray) -> dict:
    sym = [g["symbol"] for g in genes]
    x = file_counts[sym.index("XIST")].astype(float)
    y = file_counts[[i for i, s in enumerate(sym) if s in sim.Y_GENES]].sum(axis=0).astype(float)
    return {s: (float(min(x[j], y[j]) / max(x[j] + y[j], 1.0)), float(x[j]), float(y[j])) for j, s in enumerate(sids)}


# --------------------------------------------------------------------------- H: the predecessor's files

Q01 = r'''#!/bin/bash
# Quantification of all libraries of the sample sheet with salmon.
# usage: bash 01_quantify.sh <work_dir>
# (PROJECT_INPUT overrides the input location, default /task/input)
set -euo pipefail
IN=${PROJECT_INPUT:-/task/input}
WORK=${1:?usage: 01_quantify.sh <work_dir>}
mkdir -p "$WORK/quant"
salmon index -t "$IN/reference/transcripts.fa.gz" -i "$WORK/index" -k 31 -p 2 > "$WORK/index.log" 2>&1
tail -n +2 "$IN/samples.tsv" | cut -f1,6 | while IFS=$'\t' read -r sid fq; do
  salmon quant -i "$WORK/index" -l @LIBTYPE@ -r "$IN/fastq/$fq" --validateMappings -p 2 \
    -o "$WORK/quant/$sid" > "$WORK/quant/$sid.log" 2>&1 < /dev/null
done
echo "quantified $(ls "$WORK/quant" | grep -vc '\.log$') libraries into $WORK/quant"
'''

Q02 = r'''#!/usr/bin/env python3
"""Gene-level matrix from the salmon output of 01_quantify.sh.
usage: python3 02_merge_counts.py <work_dir> <out.tsv>
"""
import os
import sys

import pandas as pd

IN = os.environ.get("PROJECT_INPUT", "/task/input")
work, out = sys.argv[1], sys.argv[2]
samples = pd.read_csv(f"{IN}/samples.tsv", sep="\t", dtype=str)
t2g = pd.read_csv(f"{IN}/reference/tx2gene.tsv", sep="\t", dtype=str)
gene_of = dict(zip(t2g["transcript_id"], t2g["gene_id"]))
cols = {}
for sid in samples["sample_id"]:
    sf = pd.read_csv(f"{work}/quant/{sid}/quant.sf", sep="\t")
    cols[sid] = sf.groupby(sf["Name"].map(gene_of))["@QUANTITY@"].sum()
matrix = pd.DataFrame(cols).reindex(list(dict.fromkeys(t2g["gene_id"]))).fillna(0.0)
matrix.index.name = "gene_id"
matrix.round(3).to_csv(out, sep="\t")
'''

Q03 = r'''#!/usr/bin/env python3
"""Differential expression, compound vs vehicle (DESeq2 model via pydeseq2).
usage: python3 03_de.py <gene_counts.tsv> <out.tsv> [samples.tsv]
"""
import os
import sys

import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
@EXTRA_IMPORT@
IN = os.environ.get("PROJECT_INPUT", "/task/input")
counts = pd.read_csv(sys.argv[1], sep="\t", index_col=0)
sheet = sys.argv[3] if len(sys.argv) > 3 else f"{IN}/samples.tsv"
samples = pd.read_csv(sheet, sep="\t", dtype=str).set_index("sample_id")
names = pd.read_csv(f"{IN}/reference/tx2gene.tsv", sep="\t", dtype=str).drop_duplicates("gene_id")
names = names.set_index("gene_id")["gene_name"]

meta = pd.DataFrame({"condition": samples["condition"],
                     "session": ["s" + d.replace("-", "") for d in samples["library_prep_batch"]]})
mat = counts[meta.index].T.round().astype(int)
dds = DeseqDataSet(counts=mat, metadata=meta, design="@DESIGN@", quiet=True, n_cpus=2)
dds.deseq2()
st = DeseqStats(dds, contrast=@CONTRAST@, alpha=0.05, quiet=True, n_cpus=2)
st.summary()
res = st.results_df.reindex(counts.index)
@PADJ_BLOCK@
table = pd.DataFrame({"gene_id": counts.index,
                      "gene_name": names.reindex(counts.index).values,
                      "base_mean": res["baseMean"].round(2).values,
                      "log2fc": (@LFC@).round(4).values,
                      "pvalue": res["pvalue"].values,
                      "padj": padj.values})
table.to_csv(sys.argv[2], sep="\t", index=False, na_rep="NA", float_format="%.6g")
'''

PADJ_DESEQ = 'padj = res["padj"]'
PADJ_MANUAL = ('tested = res["pvalue"].notna()\n'
               'padj = pd.Series(float("nan"), index=res.index)\n'
               'padj[tested] = multipletests(res.loc[tested, "pvalue"], method="fdr_bh")[1]')
PADJ_RAW = 'padj = res["pvalue"]'

HANDOVER = """\
# Hand-over note: compound vs vehicle RNA-seq

Written on my last day. Sorry for the state of things.

## What is in this folder
- `scripts/01_quantify.sh` - salmon index and quantification of the {n} FASTQ files
- `scripts/02_merge_counts.py` - gene-level matrix from the salmon output
- `scripts/03_de.py` - DESeq2 model (pydeseq2), compound vs vehicle
- `gene_counts_merged.tsv` - gene-level matrix (output of 01 + 02)
- `de_draft.tsv` - draft differential-expression table (output of 03)
- `mail_top_responders.tsv` - the {n_mailed} genes I mailed to the wet lab on {mail_date} as "top responder genes"

## Notes
- The salmon output directories were on the scratch disk, which has been purged since. Only the
  merged matrix is left.
- The mailed list came from my first differential-expression run (smallest adjusted p-values first,
  sex-chromosome control genes taken out). The table of that run was overwritten; `de_draft.tsv` is
  from a later run.
- Samples were used as listed in the sample sheet. I did not get to the identity control of the SOP.
- The wet lab has started to order qPCR primers for the mailed genes, so they need to know soon
  which of them hold up.
- All scripts take their paths as arguments, see the usage line at the top of each.
"""

SOP = """\
# SOP-BI-07: Taking over an RNA-seq differential-expression project

## 1. Principle
The analyst who releases a result is responsible for every number in it. Files left by a previous
analyst (scripts, intermediate tables, drafts) are not raw data. They may be reused only after they
have been verified against the raw data; anything may be recomputed instead. A script is judged by
what it does and a file by what it contains, not by comments, file names or the hand-over note.

## 2. Identity control of libraries
The gene panel contains sex-chromosome control genes: `XIST` (transcribed in female cells only) and
the Y-linked genes `RPS4Y1`, `DDX3Y`, `KDM5D`, `UTY`, `EIF1AY` (transcribed in male cells only).

For every FASTQ file, from the reads of that file that are assigned to genes:
- X = number of reads assigned to `XIST`
- Y = number of reads assigned to the five Y-linked genes together
- minor marker share = min(X, Y) / (X + Y)

Acceptance criterion: a library is accepted as coming from a single donor if its minor marker share
is below {limit}. Its sex is then `F` if X > Y and `M` otherwise.

A library with a minor marker share of {limit} or more contains material of more than one donor. Such
a library has no identity: it cannot be assigned to any record of the sample sheet, it is reported
with `inferred_sex` = `mixed`, and it does not enter the statistical model (`condition_used` =
`excluded`).

Passing this control shows that a library comes from one donor. It does not show that this donor is
the one named in the sample sheet.

A differential-expression table is released only if every library in the model passed the identity
control and entered the model under the condition its cells actually received. A library that
passed the control and whose data agree with its sample-sheet record is never left out of the model.

## 3. Audit of the previous analysis
Every item of the checklist below gets a verdict in `audit.tsv`: `ok` if the requirement is met,
`defective` if it is not.

- Judge every item on its own. A step that is carried out correctly on a faulty input is `ok`; the
  fault belongs to the item that caused it.
- The items concern the files in `previous_analysis/` as shipped (scripts and tables), not your own
  analysis.
- The checklist is complete: nothing else about the previous analysis is rated. In particular, the
  identity control of section 2 was not part of the previous analysis and is not an audit item.

| item | requirement |
|---|---|
| `quant_library_type` | The library type that `scripts/01_quantify.sh` passes to the quantifier fits the protocol of these libraries, so that no reads are discarded because of their orientation. |
| `merge_read_counts` | `scripts/02_merge_counts.py` builds the gene-level matrix from estimated read counts (not from length- or depth-normalised abundances), summed over all transcripts of a gene. |
| `matrix_sample_files` | Every column of `gene_counts_merged.tsv` was obtained from the FASTQ file that `samples.tsv` lists for the sample in the column header. |
| `matrix_all_reads` | Every column of `gene_counts_merged.tsv` was obtained from all reads of a FASTQ file, not from a part of a file. |
| `de_session_term` | The model in `scripts/03_de.py` contains the library-prep session in addition to the treatment. |
| `de_contrast_direction` | The `log2fc` that `scripts/03_de.py` writes is treated relative to vehicle: positive when expression is higher in treated cells. |
| `de_adjusted_pvalues` | The `padj` that `scripts/03_de.py` writes is the p-value adjusted for multiple testing over the tested genes. |
| `draft_reproducible` | `de_draft.tsv` is the output of `scripts/03_de.py` for `gene_counts_merged.tsv` and `/task/input/samples.tsv` as shipped (the same numbers, apart from rounding in the last printed digit). |

## 4. Erratum for results that have already left the group
For every gene of a list that was communicated outside the group, `erratum.tsv` states whether the
released analysis supports it: `confirmed` if the gene is differentially expressed in the released
table (`significant` = 1), `retracted` otherwise. Direction and size of an effect are taken from the
released table, not from the erratum.
"""

INSTRUCTION = """\
You are the bioinformatician of a pharmacology group and take over an RNA-seq project from a
colleague who has left. The wet lab treated donor-derived primary cells with a candidate compound or
with vehicle and prepared poly-A RNA-seq libraries. Your predecessor quantified the libraries, ran a
differential-expression analysis and mailed a list of "top responder genes" to the wet lab. You
release the final result: gene-level counts, a per-sample QC table, the genes that respond to the
compound, an audit of the previous analysis and an erratum for the mailed list.

Task directory:
- `/task`

Inputs (read-only) under `/task/input`:
- `fastq/<sample_id>.fastq.gz` - one single-end, unstranded FASTQ per sample ({n_samples} samples)
- `samples.tsv` - sample sheet as filled in by the wet lab (`sample_id`, `condition` = `treated` or
  `vehicle`, `donor_id`, `donor_sex` = `F` or `M`, `library_prep_batch`, `fastq_file`)
- `reference/transcripts.fa.gz` - transcript sequences of the gene panel that was sequenced
- `reference/tx2gene.tsv` - `transcript_id`, `gene_id`, `gene_name` for every transcript
- `lab_notes.md` - wet-lab notes for this experiment
- `SOP_takeover.md` - the group's SOP for taking over an analysis (identity control of libraries,
  audit checklist, erratum). It is binding for this task.
- `previous_analysis/` - what your predecessor left: `HANDOVER.md`, `scripts/01_quantify.sh`,
  `scripts/02_merge_counts.py`, `scripts/03_de.py`, `gene_counts_merged.tsv`, `de_draft.tsv`,
  `mail_top_responders.tsv`

Available software (on `PATH`): `salmon`, Python 3.12 with numpy, scipy, pandas, statsmodels,
scikit-learn and pydeseq2. There is no internet access and no R installation.

Your task:
1. Quantify every sample against the provided transcripts and summarise to genes.
2. Check every sample (SOP section 2) and decide how it enters the statistical analysis.
3. Test every gene for differential expression between compound-treated and vehicle-treated cells.
4. Audit the previous analysis against the checklist of SOP section 3.
5. Tell the wet lab which of the mailed genes hold up (SOP section 4).

Required outputs under `/task/output`:
- `gene_counts.tsv` - tab-separated, header `gene_id` followed by one column per `sample_id`
  (all {n_samples} samples, named as in `samples.tsv`). One row per `gene_id` of `tx2gene.tsv`.
  Values are the estimated number of reads of that FASTQ file that originate from the gene
  (non-negative, not length- or depth-normalised; fractional values are allowed).
- `sample_qc.tsv` - tab-separated with at least the columns `sample_id`, `inferred_sex`,
  `condition_used`; one row per sample of `samples.tsv`.
  - `inferred_sex`: `F` or `M`, the sex of the donor of the library contained in that sample's
    FASTQ file, as determined from the sequencing data; `mixed` for a library that fails the
    identity control of the SOP.
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
- `audit.tsv` - tab-separated with at least the columns `item`, `verdict`; one row per checklist item
  of SOP section 3, `item` exactly as written there; `verdict` = `ok` or `defective`.
- `erratum.tsv` - tab-separated with at least the columns `gene_id`, `verdict`; one row per gene of
  `previous_analysis/mail_top_responders.tsv`; `verdict` = `confirmed` if the gene is
  differentially expressed in your final analysis (`significant` = 1 in your `de_results.tsv`),
  otherwise `retracted`.

Constraints:
- Do not modify anything under `/task/input`; write deliverables only under `/task/output`.
- Extra columns are allowed in every output file.
- Do not ask for confirmation. Execute directly.
"""

JUDGE = """\
**Where the truth comes from.** Everything scored is planted. Before any read exists the generator
fixes which genes respond and by how much, the library-prep batch effects, sex and treatment of
every library, which two FASTQ files hold each other's library, which library is a mixture of two
donors, and which checklist items of the SOP are defective in the predecessor's files.
`reference/truth.json` stores exactly these planted values. The predecessor's files are produced by
running mutated copies of the reference pipeline on the task's own reads; they are inputs, never
the reference. The erratum truth is the planted responder status of the mailed genes.

**Hard gates (score 0).** As at level 3 (missing/unparseable file, gene or sample missing or
duplicated, negative counts, invalid `condition_used`, `significant` inconsistent with `padj`,
non-finite `log2fc`), plus: `audit.tsv` must give `ok`/`defective` for every checklist item exactly
once; `erratum.tsv` must give `confirmed`/`retracted` for every mailed gene exactly once.

**Components.**

| weight | component | metric | full credit | zero |
|---|---|---|---|---|
| 0.08 | gene_counts | per sample `sum(min(obs,true)) / max(sum obs, sum true)`, minimum over samples; 0 when a column agrees better with the planted counts of another FASTQ file than with those of its own | >= 0.90 | <= 0.50 |
| 0.02 | sample_sex | `inferred_sex` is the planted sex for every library, `mixed` for the mixed one; only if gene_counts > 0 | all | any wrong |
| 0.07 | sample_labels | mixed library `excluded`; every other library under its planted condition (the two mixed-up ones may be `excluded`); nobody else dropped; only if gene_counts > 0 | all | any wrong |
| 0.19 | de_false_discoveries | as level 3 | <= 0.12 | >= 0.30 |
| 0.16 | de_sensitivity | as level 3 | >= 0.85 | <= 0.50 |
| 0.13 | de_log2fc | as level 3 | <= 0.33 | >= 0.55 |
| 0.20 | audit | Matthews correlation of the verdicts against the planted defects | 1.0 (complete audit) | <= 0.70 |
| 0.15 | erratum | Matthews correlation of `confirmed`/`retracted` against the planted responder status of the mailed genes | >= 0.80 | <= 0.40 |

The three DE components are multiplied by 0.5 when the sample table is not valid (SOP section 2: a
table from a model that contains the mixed library, or a library under the wrong condition, is not
released). One wrong audit verdict gives MCC 0.745-0.775 on eight items with three to five defects,
i.e. 0.15-0.25 of the audit credit and a task score of at most 0.85 with everything else perfect: a
single missed defect or false alarm ends below the pass mark. Two wrong verdicts give no audit credit.

**Why these tolerances.** See the "Level 4" section of `docs/families/rnaseq_confound.md` for the
measured values of the oracle, the independent solution, the level-3 reference method and every
careless variant on two seeds.

**What is deliberately not scored.** The numeric value of `padj`, `log2fc` of genes without a planted
effect, the DE status of the sex-chromosome control genes, the evidence column of the audit, which
of the two valid resolutions of the file mix-up was chosen.
"""


def _render_scripts(plan: dict) -> dict:
    d = plan["defects"]
    lib = plan["libtype_bad"] if "quant_library_type" in d else plan["libtype_ok"]
    q01 = Q01.replace("@LIBTYPE@", lib)
    q02 = Q02.replace("@QUANTITY@", "TPM" if "merge_read_counts" in d else "NumReads")
    if "de_contrast_direction" in d:
        contrast, lfc = '["condition", "vehicle", "treated"]', 'res["log2FoldChange"]'
    elif plan["contrast_style"] == "negated":
        contrast, lfc = '["condition", "vehicle", "treated"]', '-res["log2FoldChange"]'
    else:
        contrast, lfc = '["condition", "treated", "vehicle"]', 'res["log2FoldChange"]'
    if "de_adjusted_pvalues" in d:
        padj, extra = PADJ_RAW, ""
    elif plan["padj_style"] == "manual_bh":
        padj, extra = PADJ_MANUAL, "from statsmodels.stats.multitest import multipletests\n"
    else:
        padj, extra = PADJ_DESEQ, ""

    def q03(design, contrast, lfc, padj, extra):
        return (Q03.replace("@DESIGN@", design).replace("@CONTRAST@", contrast).replace("@LFC@", lfc)
                .replace("@PADJ_BLOCK@", padj).replace("@EXTRA_IMPORT@\n", extra))

    return {"01_quantify.sh": q01, "02_merge_counts.py": q02,
            "03_de.py": q03("~condition" if "de_session_term" in d else "~session + condition", contrast, lfc, padj, extra),
            # the predecessor's first run: condition-only model, everything else plain (not shipped)
            "03_first_pass.py": q03("~condition", '["condition", "treated", "vehicle"]', 'res["log2FoldChange"]',
                                    PADJ_DESEQ, "")}


RECOMPUTE = [i for i, k in ITEMS if k == "recompute"]
# code-level items that can be planted as defective. `de_adjusted_pvalues` is never planted: raw
# p-values written as padj gave an FDP of 0.08-0.23 on the level-4 instances of seeds 1-3 (planted
# counts) and 0.20 on the reads of seed 1, where the analysis that keeps the defect scored 0.913
# before the audit miss is counted: not a decision point. The item stays in the checklist in two
# sound variants (pyDESeq2's padj, hand-written Benjamini-Hochberg).
CODE_MUTABLE = ["quant_library_type", "merge_read_counts", "de_session_term", "de_contrast_direction"]


def plan_defects(rng, inst: dict) -> dict:
    """Seeded choice of the defective checklist items (3-5 of 8: two of the three that only a
    recomputation shows, one to three of the code-level ones) and of the samples that carry the
    file-level defects."""
    sheet = inst["sheet"]
    special = set(inst["truth"]["swapped"]) | {inst["mixture"]["file_id"]}
    plain = [d for d in sheet if d["sheet_id"] not in special]
    cells: dict[tuple, list] = {}
    for d in sheet:
        cells.setdefault(_cell(d), []).append(d)
    major = sorted(k for k, v in cells.items() if len(v) == max(len(x) for x in cells.values()))
    # exactly two of the three recompute items are defective: rating all three `defective` without
    # recomputing is one false alarm, which is below the pass mark
    rec = [RECOMPUTE[int(i)] for i in rng.permutation(len(RECOMPUTE))[:2]]
    # one to three code-level defects. A partial file cannot be told from a complete one by depth
    # once the matrix holds TPM, so the TPM defect is only planted without the partial file (and
    # then with probability 0.75, otherwise it would be rare)
    k_c = int(rng.integers(1, 4))
    code = []
    if "matrix_all_reads" not in rec and rng.random() < 0.75:
        code.append("merge_read_counts")
    others = [c for c in CODE_MUTABLE if c != "merge_read_counts"]
    code += [others[int(i)] for i in rng.permutation(len(others))[:k_c - len(code)]]
    defects = set(rec) | set(code)
    # column X holds the quantification of Y's FASTQ: other condition and session, same recorded sex
    first = major[int(rng.integers(0, len(major)))]
    other = [k for k in major if k != first][0]
    pairs = [(x, y) for x in plain if _cell(x) == first for y in plain if _cell(y) == other and y["sex"] == x["sex"]]
    if not pairs:
        pairs = [(x, y) for x in plain if _cell(x) == first for y in plain if _cell(y) == other]
    x, y = pairs[int(rng.integers(0, len(pairs)))]
    rest = [d for d in plain if d is not x and d is not y]
    z = rest[int(rng.integers(0, len(rest)))]
    frac = float(rng.uniform(0.28, 0.42))
    rest_t = [d for d in rest if d is not z and d["condition"] == "treated"]
    rest_v = [d for d in rest if d is not z and d["condition"] == "vehicle"]
    p = rest_t[int(rng.integers(0, len(rest_t)))]
    q = rest_v[int(rng.integers(0, len(rest_v)))]
    j = inst["sids"].index(z["sheet_id"])
    n_reads = int(inst["file_counts"][:, j].sum())
    return {
        "defects": sorted(defects, key=ITEM_IDS.index),
        "libtype_ok": "A" if rng.random() < 0.5 else "U",
        "libtype_bad": "SF" if rng.random() < 0.5 else "SR",
        "contrast_style": "plain" if rng.random() < 0.5 else "negated",
        "padj_style": "deseq" if rng.random() < 0.5 else "manual_bh",
        "dup": [x["sheet_id"], y["sheet_id"]], "trunc": [z["sheet_id"], int(round(frac * n_reads))],
        "stale": [p["sheet_id"], q["sheet_id"]],
    }


GEN_SCRIPT = r'''
set -euo pipefail
export HOME=/tmp
S=/stage/previous_analysis/scripts
mkdir -p /out
bash $S/01_quantify.sh /tmp/w > /out/01.log
python3 $S/02_merge_counts.py /tmp/w /out/matrix_script.tsv
if [ -n "${TRUNC_SAMPLE:-}" ]; then
  mkdir -p /tmp/in2/fastq
  ln -s /task/input/reference /tmp/in2/reference
  head -n 1 /task/input/samples.tsv > /tmp/in2/samples.tsv
  awk -F'\t' -v s="$TRUNC_SAMPLE" '$1==s' /task/input/samples.tsv >> /tmp/in2/samples.tsv
  FQ=$(awk -F'\t' -v s="$TRUNC_SAMPLE" '$1==s {print $6}' /task/input/samples.tsv)
  python3 - "$FQ" "$TRUNC_READS" <<'PY'
import gzip, sys
fq, n = sys.argv[1], int(sys.argv[2])
with gzip.open("/task/input/fastq/" + fq, "rb") as src, gzip.open("/tmp/in2/fastq/" + fq, "wb", compresslevel=1) as dst:
    for i, line in enumerate(src):
        if i >= 4 * n:
            break
        dst.write(line)
PY
  PROJECT_INPUT=/tmp/in2 bash $S/01_quantify.sh /tmp/w2 > /out/01_partial.log
  PROJECT_INPUT=/tmp/in2 python3 $S/02_merge_counts.py /tmp/w2 /out/matrix_partial.tsv
fi
python3 /stage/gen/assemble_previous.py
'''


def ensure_gentools() -> str:
    """Tool image for the generation-time run of the predecessor's scripts. It is built from the same
    Dockerfile text as the task's agent image (empty input), so the conda layer is shared with every
    task image of the family and the run uses exactly the software the agent will have."""
    if subprocess.run(["docker", "image", "inspect", GENTOOLS_TAG], capture_output=True).returncode == 0:
        return GENTOOLS_TAG
    ctx = Path(tempfile.mkdtemp(prefix="alegen-rnaseq-gentools-"))
    try:
        for sub in ("input", "software"):
            (ctx / sub).mkdir()
            (ctx / sub / ".keep").write_text("")
        (ctx / "Dockerfile").write_text(image.agent_dockerfile())
        dockerlib.build_image(ctx, GENTOOLS_TAG)
    finally:
        shutil.rmtree(ctx, ignore_errors=True)
    return GENTOOLS_TAG


def build_previous(inst: dict, plan: dict, workdir: Path) -> dict:
    """Run the predecessor's (mutated) scripts on the task's reads and return the files left behind."""
    workdir = Path(workdir)
    stage = workdir / "l4_stage"
    if stage.exists():
        shutil.rmtree(stage)
    scripts = _render_scripts(plan)
    (stage / "previous_analysis" / "scripts").mkdir(parents=True)
    (stage / "gen").mkdir()
    (stage / "reference").mkdir()
    for name in ("01_quantify.sh", "02_merge_counts.py", "03_de.py"):
        (stage / "previous_analysis" / "scripts" / name).write_text(scripts[name])
    (stage / "gen" / "03_first_pass.py").write_text(scripts["03_first_pass.py"])
    shutil.copyfile(A4 / "gen" / "assemble_previous.py", stage / "gen" / "assemble_previous.py")
    d = plan["defects"]
    (stage / "gen" / "plan.json").write_text(json.dumps({
        "dup": plan["dup"] if "matrix_sample_files" in d else None,
        "trunc": plan["trunc"] if "matrix_all_reads" in d else None,
        "stale": plan["stale"] if "draft_reproducible" in d else None,
        "control_genes": ["XIST"] + sim.Y_GENES, "n_mailed": N_MAILED}))
    (stage / "samples.tsv").write_text(inst["samples_tsv"])
    (stage / "reference" / "tx2gene.tsv").write_text(inst["tx2gene"])
    shutil.copyfile(inst["transcripts"], stage / "reference" / "transcripts.fa.gz")
    env = {}
    if "matrix_all_reads" in d:
        env = {"TRUNC_SAMPLE": plan["trunc"][0], "TRUNC_READS": str(plan["trunc"][1])}
    out = workdir / "l4_out"
    r = dockerlib.run_in_image(
        ensure_gentools(), GEN_SCRIPT,
        mounts={stage: "/stage", workdir / "fastq": "/task/input/fastq", stage / "reference": "/task/input/reference",
                stage / "samples.tsv": "/task/input/samples.tsv"},
        copy_out={"/out": out}, env=env, cpus=2, timeout=1800)
    if r.returncode != 0 or not (out / "diag.json").exists():
        raise RuntimeError("generation of the previous analysis failed:\n" + (r.stdout or "")[-1500:] + (r.stderr or "")[-3000:])
    diag = json.loads((out / "diag.json").read_text())
    if not diag["rerun_identical"]:
        raise RuntimeError("the shipped DE script is not reproducible run to run; `draft_reproducible` would be ambiguous")
    if "draft_reproducible" in d and not (diag["stale_max_abs_lfc_diff"] >= 0.5 and diag["stale_genes_diff_gt_0.1"] >= 50):
        raise RuntimeError(f"stale draft too close to the reproducible one: {diag}")
    prev = out / "prev"
    files = {f"previous_analysis/scripts/{n}": scripts[n] for n in ("01_quantify.sh", "02_merge_counts.py", "03_de.py")}
    for n in ("gene_counts_merged.tsv", "de_draft.tsv"):
        files[f"previous_analysis/{n}"] = (prev / n).read_text()
    return {"files": files, "diag": diag}


def compose_mailed(rng, ranked: list[str], lfc: dict) -> list[str]:
    """Mailed list = the best-ranked genes of the predecessor's first run, composed so that 16-24 of
    the 40 are planted responders. (The plain top 40 of that run held only 4-7 non-responders on three
    seeds, which makes a correlation score a lottery; see the design notes.)"""
    n_true = int(rng.integers(16, 25))
    resp = [g for g in ranked if lfc[g] != 0.0][:n_true]
    rest = [g for g in ranked if lfc[g] == 0.0][:N_MAILED - n_true]
    chosen = set(resp) | set(rest)
    return [g for g in ranked if g in chosen]


def _tsv(header: list[str], rows: list[list]) -> str:
    return "\t".join(header) + "\n" + "".join("\t".join(str(v) for v in r) + "\n" for r in rows)


def generate(seed: int, level: int, workdir: Path) -> TaskBuild:
    from . import _lab_notes, _tsv_counts, _tsv_de, build_instance   # noqa: F401  (level-3 helpers)
    workdir = Path(workdir)
    inst = build_instance(seed, level, workdir)
    rng4 = inst["rng4"]
    truth, sids, gene_ids = inst["truth"], inst["sids"], inst["truth"]["genes"]
    fc = inst["file_counts"]
    lfc = truth["lfc"]
    unscored = set(truth["unscored_genes"])
    mix = inst["mixture"]
    mixed_id = mix["file_id"]

    # Q: margins of the stated identity criterion on the planted counts
    shares = marker_shares(inst["genes"], sids, fc)
    worst_pure = max(v[0] for s, v in shares.items() if s != mixed_id)
    if not (shares[mixed_id][0] >= 2.0 * MINOR_SHARE_LIMIT and worst_pure <= 0.5 * MINOR_SHARE_LIMIT):
        raise RuntimeError(f"seed {seed}: identity criterion without margin (mixed {shares[mixed_id][0]:.4f}, "
                           f"purest-failing {worst_pure:.4f}); use another seed")

    # H: planted defects, then the predecessor's files
    plan = plan_defects(rng4, inst)
    prev = build_previous(inst, plan, workdir)
    name_of = {g["gene_id"]: g["symbol"] for g in inst["genes"]}
    ranked = [ln.split("\t")[0] for ln in (workdir / "l4_out" / "first_pass_ranked.tsv").read_text().splitlines()[1:]]
    mailed = compose_mailed(rng4, ranked, lfc)
    mail_date = "2 July" if inst["batch_names"][0].startswith("2024-05") else "9 April"
    prev_files = dict(prev["files"])
    prev_files["previous_analysis/mail_top_responders.tsv"] = _tsv(
        ["rank", "gene_id", "gene_name"], [[i + 1, g, name_of[g]] for i, g in enumerate(mailed)])
    prev_files["previous_analysis/HANDOVER.md"] = HANDOVER.format(n=len(sids), n_mailed=len(mailed), mail_date=mail_date)
    truth["audit"] = {i: ("defective" if i in plan["defects"] else "ok") for i in ITEM_IDS}
    truth["mailed"] = {g: ("confirmed" if lfc[g] != 0.0 else "retracted") for g in mailed}
    truth["plan"] = plan
    n_resp = sum(v == "confirmed" for v in truth["mailed"].values())

    # ------------------------------------------------------------------ truth-derived submissions
    planted = [g for g in gene_ids if lfc[g] != 0 and g not in unscored]
    null = [g for g in gene_ids if lfc[g] == 0 and g not in unscored]
    sheet_cond = {d["sheet_id"]: d["condition"] for d in inst["sheet"]}
    sheet_sex = {d["sheet_id"]: d["sex"] for d in inst["sheet"]}
    prng = np.random.default_rng(seed + 4977)
    sig_true = {g: lfc[g] != 0 and g not in unscored for g in gene_ids}
    thinned = prng.binomial(fc, 0.94)
    fp_in = set(prng.permutation(null)[:4].tolist())
    miss_in = set(prng.permutation(planted)[:int(0.08 * len(planted))].tolist())
    sig_in = {g: (sig_true[g] and g not in miss_in) or g in fp_in for g in gene_ids}
    lfc_in = {g: lfc[g] + float(prng.normal(0, 0.2)) for g in gene_ids}
    n_fp_out = int(round(len(planted) * 0.40 / 0.60))
    fp_out = set(prng.permutation(null)[:n_fp_out].tolist())
    sig_out = {g: sig_true[g] or g in fp_out for g in gene_ids}
    lfc_out = {g: (lfc[g] if g not in fp_out else 0.8) for g in gene_ids}
    lfc_att = {g: lfc[g] * (2.0 / 3.0) for g in gene_ids}

    def qc_tsv(sex: dict, cond: dict) -> str:
        return _tsv(["sample_id", "inferred_sex", "condition_used"], [[s, sex[s], cond[s]] for s in sids])

    def audit_tsv(v: dict) -> str:
        return _tsv(["item", "verdict"], [[i, v[i]] for i in ITEM_IDS])

    def erratum_tsv(v: dict) -> str:
        return _tsv(["gene_id", "verdict"], [[g, v[g]] for g in mailed])

    true_sex = {s: ("mixed" if s == mixed_id else truth["sex"][s]) for s in sids}
    true_cond = {s: ("excluded" if s == mixed_id else truth["condition"][s]) for s in sids}
    all_ok = {i: "ok" for i in ITEM_IDS}
    flip = [g for g in mailed if truth["mailed"][g] == "confirmed"][-1:] + [g for g in mailed if truth["mailed"][g] == "retracted"][:1]
    err_in = {g: (v if g not in flip else ("retracted" if v == "confirmed" else "confirmed")) for g, v in truth["mailed"].items()}
    T = {"gene_counts.tsv": _tsv_counts(gene_ids, sids, fc),
         "sample_qc.tsv": qc_tsv(true_sex, true_cond),
         "de_results.tsv": _tsv_de(gene_ids, lfc, sig_true),
         "audit.tsv": audit_tsv(truth["audit"]),
         "erratum.tsv": erratum_tsv(truth["mailed"])}

    def but(**kw) -> dict:
        return {**T, **{k.replace("__", "."): v for k, v in kw.items()}}

    one_miss = dict(truth["audit"], **{plan["defects"][-1]: "ok"})
    false_alarm_item = [i for i in ITEM_IDS if i not in plan["defects"]][0]
    one_alarm = dict(truth["audit"], **{false_alarm_item: "defective"})
    by_reading = {i: (truth["audit"][i] if dict(ITEMS)[i] == "code" else "ok") for i in ITEM_IDS}

    pipe = {"takeover_pipeline.py": A4 / "takeover_pipeline.py"}
    run = "python3 /alegen_script/takeover_pipeline.py"
    probes = [
        Probe("planted_truth", "submission written directly from the planted truth", 0.999, files=T),
        Probe("within_tolerance", "6% of reads lost, 4 false discoveries, 8% of planted genes missed, log2fc noise "
              "sd 0.2, two erratum verdicts wrong (one each way); audit exact", 0.99,
              files=but(gene_counts__tsv=_tsv_counts(gene_ids, sids, thinned),
                        de_results__tsv=_tsv_de(gene_ids, lfc_in, sig_in), erratum__tsv=erratum_tsv(err_in))),
        Probe("independent_solution", "second tool chain: k-mer read counter (no salmon), voom-style weighted linear "
              "model (no DESeq2), mismatching and mixed libraries excluded, audit by profile comparison and static "
              "inspection instead of a salmon re-run", 0.99, script="SOLUTION=2 bash /solution/solve.sh"),
        Probe("oracle_exclude_resolution", "oracle tool chain, but the two mixed-up samples are excluded instead of "
              "relabelled", 0.99, script=run + " --labels exclude", script_files=pipe),
    ]
    ones = np.ones_like(fc)
    careless_extras = (
        "import csv\n"
        "items = " + repr(ITEM_IDS) + "\n"
        "sig = {r['gene_id']: r['significant'] for r in csv.DictReader(open('/task/output/de_results.tsv'), delimiter='\\t')}\n"
        "open('/task/output/audit.tsv', 'w').write('item\\tverdict\\n' + ''.join(f'{i}\\tok\\n' for i in items))\n"
        "rows = list(csv.DictReader(open('/task/input/previous_analysis/mail_top_responders.tsv'), delimiter='\\t'))\n"
        "open('/task/output/erratum.tsv', 'w').write('gene_id\\tverdict\\n' + ''.join(\n"
        "    f\"{r['gene_id']}\\t{'confirmed' if sig.get(r['gene_id']) == '1' else 'retracted'}\\n\" for r in rows))\n")
    low = PASS - 0.15
    attacks = [
        Attack("l3_reference", "the level-3 reference method on the level-4 inputs: fresh salmon quantification, "
               "male/female call, mix-up resolved, session term; no SOP identity criterion; the predecessor's work is "
               "taken as fine (audit all ok); erratum from its own DE table", low,
               script="python3 /alegen_script/rnaseq_pipeline.py && python3 /alegen_script/careless_extras.py",
               script_files={"rnaseq_pipeline.py": ASSETS / "rnaseq_pipeline.py", "careless_extras.py": careless_extras}),
        Attack("trust_inherited", "H ignored: the analysis is finished from the inherited count matrix (identity "
               "control, mix-up and session term handled correctly), nothing is re-checked, audit all ok", low,
               script=run + " --counts-from inherited --audit all_ok", script_files=pipe),
        Attack("q_identity_not_applied", "Q ignored: everything right (fresh quantification, complete audit, erratum) "
               "but libraries are only typed male/female, so the mixed library stays in the model", low,
               script=run + " --identity binary", script_files=pipe),
        Attack("sheet_copied_no_analysis", "prior-only answer: QC table copied from the sample sheet, counts all 1, "
               "no gene significant, audit all ok, every mailed gene confirmed", 0.10,
               files={"gene_counts.tsv": _tsv_counts(gene_ids, sids, ones), "sample_qc.tsv": qc_tsv(sheet_sex, sheet_cond),
                      "de_results.tsv": _tsv_de(gene_ids, {g: 0.0 for g in gene_ids}, {g: False for g in gene_ids}),
                      "audit.tsv": audit_tsv(all_ok), "erratum.tsv": erratum_tsv({g: "confirmed" for g in mailed})}),
        Attack("everything_significant", "true counts, QC and audit, but every gene called significant with log2fc 1 "
               "and every mailed gene confirmed", 0.45,
               files=but(de_results__tsv=_tsv_de(gene_ids, {g: 1.0 for g in gene_ids}, {g: True for g in gene_ids}),
                         erratum__tsv=erratum_tsv({g: "confirmed" for g in mailed}))),
        # single-component failures with everything else perfect: they cannot lose more than the weight
        # of the component, so their ceilings sit below the pass mark, not at pass - 0.15
        Attack("audit_all_ok", "everything perfect, but every checklist item rated ok", 0.85,
               files=but(audit__tsv=audit_tsv(all_ok))),
        Attack("audit_all_defective", "everything perfect, but every checklist item rated defective", 0.85,
               files=but(audit__tsv=audit_tsv({i: "defective" for i in ITEM_IDS}))),
        Attack("audit_by_reading_only", "everything perfect, scripts audited correctly by reading them, but nothing "
               "recomputed: the items that only a recomputation can show are rated ok", 0.85,
               files=but(audit__tsv=audit_tsv(by_reading))),
        Attack("audit_one_defect_missed", f"everything perfect, but the planted defect `{plan['defects'][-1]}` is rated ok",
               0.89, files=but(audit__tsv=audit_tsv(one_miss))),
        Attack("audit_one_false_alarm", f"everything perfect, but the sound item `{false_alarm_item}` is rated defective",
               0.89, files=but(audit__tsv=audit_tsv(one_alarm))),
        Attack("erratum_all_confirmed", "everything perfect, but the mailed list is confirmed wholesale", 0.89,
               files=but(erratum__tsv=erratum_tsv({g: "confirmed" for g in mailed}))),
        Attack("fdp_beyond_tolerance", "everything perfect, but 40% of the significant genes are false", 0.85,
               files=but(de_results__tsv=_tsv_de(gene_ids, lfc_out, sig_out))),
        Attack("lfc_attenuated_beyond_tolerance", "everything perfect, but every log2fc shrunk by a third", 0.89,
               files=but(de_results__tsv=_tsv_de(gene_ids, lfc_att, sig_true))),
        Attack("naive_labels_as_given", "level-3 trap ignored: SOP identity criterion applied, fresh quantification, "
               "complete audit, but the sample sheet is used as delivered", low,
               script=run + " --labels given", script_files=pipe),
    ]
    # one variant per planted defect that can propagate: the analyst's own pipeline repeats the defect
    # and the item is rated ok; everything else is right
    replay = {
        "quant_library_type": f"--libtype {plan['libtype_bad']}",
        "merge_read_counts": "--quantity TPM",
        "matrix_sample_files": f"--dup {plan['dup'][0]}:{plan['dup'][1]}",
        "matrix_all_reads": f"--trunc {plan['trunc'][0]}:{plan['trunc'][1]}",
        "de_session_term": "--batch no",
        "de_contrast_direction": "--contrast reversed",
    }
    for item in plan["defects"]:
        if item in replay:
            attacks.append(Attack(
                f"defect_kept_{item}", f"planted defect `{item}` neither reported nor repaired: the final analysis "
                "repeats it, everything else is right", low,
                script=f"{run} {replay[item]} --audit-force {item}=ok", script_files=pipe))
    a, b = truth["swapped"]
    traps = [
        Trap("batch_confound",
             "Two library-prep sessions, 5 treated + 1 vehicle in one and 1 treated + 5 vehicle in the other; 30% of "
             "genes differ between sessions by 0.5-2 log2 units. A condition-only model reports batch genes as "
             "treatment genes.",
             "samples.tsv has a library_prep_batch column that is almost aliased with condition; lab_notes.md gives "
             "two prep sessions; the SOP checklist asks whether the predecessor's model contains the session",
             "Approach / wrong strategy"),
        Trap("sample_swap",
             f"The FASTQ files of {a} and {b} hold each other's library (one treated, one vehicle, opposite donor sex, "
             "different prep sessions).",
             "the required inferred_sex column and SOP section 2 make the analyst type every library from XIST and "
             "chrY genes; exactly two libraries contradict donor_sex and they are a reciprocal pair; lab_notes.md "
             "records that all tubes were re-labelled from a box map",
             "Understanding / unverified assumption about inputs"),
        Trap("H_inherited_project",
             "previous_analysis/ holds the predecessor's scripts, merged matrix, draft table and mailed gene list. "
             f"Planted defects ({len(plan['defects'])} of {len(ITEM_IDS)} checklist items): {', '.join(plan['defects'])}. "
             "The other items are sound; `quant_library_type` uses "
             f"`-l {plan['libtype_ok'] if 'quant_library_type' not in plan['defects'] else plan['libtype_bad']}`, the "
             f"contrast is written in the `{plan['contrast_style']}` style and padj in the `{plan['padj_style']}` style. "
             f"{N_MAILED - n_resp} of the {N_MAILED} mailed genes are not planted responders. Reusing the inherited "
             "matrix or draft, or rating the files without recomputing, fails.",
             "SOP section 1 forbids unverified reuse and section 3 lists the eight requirements; code-level defects "
             "are visible in the three short scripts; file-level defects appear when the scripts are re-run on the "
             "FASTQ files and compared with the shipped matrix / draft (identical columns, a column with a third of "
             "the reads, a draft that the script does not reproduce)",
             "Verification / trusting inherited intermediate results"),
        Trap("Q_identity_control_mixed_library",
             f"The library in {mixed_id}.fastq.gz is a mixture: {1 - mix['fraction']:.0%} of the recorded male "
             f"donor's library and {mix['fraction']:.0%} of a female library of the other condition and session "
             f"({mix['contaminant_file_id']}). Its Y-gene reads still outnumber XIST, so a male/female call agrees "
             f"with the sample sheet; its minor marker share is {shares[mixed_id][0]:.3f} against the SOP limit of "
             f"{MINOR_SHARE_LIMIT} (largest value among the other libraries: {worst_pure:.4f}). It cannot be "
             "relabelled and has to be excluded.",
             "SOP section 2 states the acceptance criterion on the sex-chromosome control genes of the panel; "
             "computing X and Y per library, as the SOP defines them, shows one male library with XIST reads at a "
             "third to a half of the female level",
             "Verification / acceptance criterion on an in-data control not evaluated"),
    ]
    fastq_files = {f"fastq/{s}.fastq.gz": inst["fastq"][s] for s in sids}
    return TaskBuild(
        task_id=f"rnaseq-confound-s{seed:04d}-l{level}", family="rnaseq_confound", seed=seed, level=level,
        title="Bulk RNA-seq: taking over a half-finished differential-expression project",
        summary="Take over an RNA-seq project from a predecessor: re-quantify twelve donor samples, apply the SOP's "
                "identity control, deliver gene counts, sample QC and treated-vs-vehicle differential expression, "
                "audit the inherited scripts and tables against a checklist and issue an erratum for a gene list "
                "that was already mailed.",
        subdomain="3.2 Genomics & Sequence Analysis",
        software=["salmon", "pyDESeq2", "statsmodels", "Python", "numpy/scipy/pandas"],
        instruction_md=INSTRUCTION.format(n_samples=len(sids)),
        agent_dockerfile=image.agent_dockerfile(),
        verifier_dockerfile=dockerlib.verifier_dockerfile(pip=[]),
        input_files={**fastq_files, "samples.tsv": inst["samples_tsv"],
                     "reference/transcripts.fa.gz": inst["transcripts"],
                     "reference/tx2gene.tsv": inst["tx2gene"], "lab_notes.md": inst["lab_notes"],
                     "SOP_takeover.md": SOP.format(limit=MINOR_SHARE_LIMIT), **prev_files},
        solution_files={"solve.sh": A4 / "solve.sh", "takeover_pipeline.py": A4 / "takeover_pipeline.py",
                        "solve_alt_l4.py": A4 / "solve_alt_l4.py"},
        tests_files={"score.py": A4 / "score.py", "reference/truth.json": json.dumps(truth)},
        traps=traps, attacks=attacks, probes=probes,
        capabilities=["auditing inherited code and intermediate files by recomputation",
                      "applying a stated acceptance criterion to an in-data control",
                      "experimental-design reasoning (aliased batch)", "sample identity QC before statistics",
                      "long tool chain from raw reads to statistics", "exact output contracts"],
        judge_explanation_md=JUDGE, pass_threshold=PASS,
        oracle_steps=29, oracle_runtime_min=1.5, cpus=4, memory_mb=8000, storage_mb=10000,
    )
