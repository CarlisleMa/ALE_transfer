#!/usr/bin/env python3
"""Differential expression, compound vs vehicle (DESeq2 model via pydeseq2).
usage: python3 03_de.py <gene_counts.tsv> <out.tsv> [samples.tsv]
"""
import os
import sys

import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
IN = os.environ.get("PROJECT_INPUT", "/task/input")
counts = pd.read_csv(sys.argv[1], sep="\t", index_col=0)
sheet = sys.argv[3] if len(sys.argv) > 3 else f"{IN}/samples.tsv"
samples = pd.read_csv(sheet, sep="\t", dtype=str).set_index("sample_id")
names = pd.read_csv(f"{IN}/reference/tx2gene.tsv", sep="\t", dtype=str).drop_duplicates("gene_id")
names = names.set_index("gene_id")["gene_name"]

meta = pd.DataFrame({"condition": samples["condition"],
                     "session": ["s" + d.replace("-", "") for d in samples["library_prep_batch"]]})
mat = counts[meta.index].T.round().astype(int)
dds = DeseqDataSet(counts=mat, metadata=meta, design="~condition", quiet=True, n_cpus=2)
dds.deseq2()
st = DeseqStats(dds, contrast=["condition", "vehicle", "treated"], alpha=0.05, quiet=True, n_cpus=2)
st.summary()
res = st.results_df.reindex(counts.index)
padj = res["padj"]
table = pd.DataFrame({"gene_id": counts.index,
                      "gene_name": names.reindex(counts.index).values,
                      "base_mean": res["baseMean"].round(2).values,
                      "log2fc": (res["log2FoldChange"]).round(4).values,
                      "pvalue": res["pvalue"].values,
                      "padj": padj.values})
table.to_csv(sys.argv[2], sep="\t", index=False, na_rep="NA", float_format="%.6g")
