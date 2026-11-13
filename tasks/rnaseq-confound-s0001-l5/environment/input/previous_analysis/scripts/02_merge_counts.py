#!/usr/bin/env python3
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
    cols[sid] = sf.groupby(sf["Name"].map(gene_of))["NumReads"].sum()
matrix = pd.DataFrame(cols).reindex(list(dict.fromkeys(t2g["gene_id"]))).fillna(0.0)
matrix.index.name = "gene_id"
matrix.round(3).to_csv(out, sep="\t")
