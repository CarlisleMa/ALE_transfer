#!/usr/bin/env python3
"""Generation-time helper (runs inside the tool image, never shipped to the agent).

Turns the output of the predecessor's scripts into the files the predecessor left behind:
applies the planted "stale file" defects to the merged matrix, runs the shipped DE script for the
draft table (on a superseded sample sheet when the draft is planted as stale) and derives the
ranking of a first-pass run (condition-only model, sample sheet as delivered) from which the
generator composes the mailed gene list.
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

plan = json.loads(Path("/stage/gen/plan.json").read_text())
scripts = Path("/stage/previous_analysis/scripts")
out = Path("/out/prev")
out.mkdir(parents=True, exist_ok=True)
diag = {}

M = pd.read_csv("/out/matrix_script.tsv", sep="\t", index_col=0)
clean = M.copy()
if plan["dup"]:
    x, y = plan["dup"]
    M[x] = clean[y]
if plan["trunc"]:
    z = plan["trunc"][0]
    part = pd.read_csv("/out/matrix_partial.tsv", sep="\t", index_col=0)
    M[z] = part[z].reindex(M.index).fillna(0.0)
    diag["partial_total"] = float(M[z].sum())
    diag["full_total"] = float(clean[z].sum())
M.round(3).to_csv(out / "gene_counts_merged.tsv", sep="\t")


def run(script, matrix, dst, sheet=None):
    cmd = [sys.executable, str(script), str(matrix), str(dst)] + ([str(sheet)] if sheet else [])
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"{script} failed:\n{r.stderr[-3000:]}")
    return pd.read_csv(dst, sep="\t").set_index("gene_id")


honest = run(scripts / "03_de.py", out / "gene_counts_merged.tsv", "/out/de_honest.tsv")
again = run(scripts / "03_de.py", out / "gene_counts_merged.tsv", "/out/de_honest2.tsv")
diag["rerun_identical"] = bool(Path("/out/de_honest.tsv").read_bytes() == Path("/out/de_honest2.tsv").read_bytes())
if plan["stale"]:
    p, q = plan["stale"]
    sheet = pd.read_csv("/task/input/samples.tsv", sep="\t", dtype=str)
    cp = sheet.loc[sheet["sample_id"] == p, "condition"].values[0]
    cq = sheet.loc[sheet["sample_id"] == q, "condition"].values[0]
    sheet.loc[sheet["sample_id"] == p, "condition"] = cq
    sheet.loc[sheet["sample_id"] == q, "condition"] = cp
    sheet.to_csv("/tmp/samples_superseded.tsv", sep="\t", index=False)
    draft = run(scripts / "03_de.py", out / "gene_counts_merged.tsv", out / "de_draft.tsv", "/tmp/samples_superseded.tsv")
    d = (draft["log2fc"].fillna(0.0) - honest["log2fc"].fillna(0.0)).abs()
    diag["stale_max_abs_lfc_diff"] = float(d.max())
    diag["stale_genes_diff_gt_0.1"] = int((d > 0.1).sum())
else:
    (out / "de_draft.tsv").write_bytes(Path("/out/de_honest.tsv").read_bytes())

first = run("/stage/gen/03_first_pass.py", out / "gene_counts_merged.tsv", "/out/de_first_pass.tsv")
first = first[~first["gene_name"].isin(plan["control_genes"])]
first = first.assign(_p=first["padj"].fillna(1.0), _q=first["pvalue"].fillna(1.0)).sort_values(["_p", "_q"], kind="stable")
# ranking of the first run; the generator composes the mailed list from it on the host
first[["gene_name", "log2fc", "pvalue", "padj"]].to_csv("/out/first_pass_ranked.tsv", sep="\t")
diag["first_pass_significant"] = int((first["padj"] < 0.05).sum())
Path("/out/diag.json").write_text(json.dumps(diag))
