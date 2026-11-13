"""Development helper: condition-only pyDESeq2 with the sample sheet as delivered, genes ranked by padj."""
import sys
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
counts = pd.read_csv(sys.argv[1], sep="\t", index_col=0)
samples = pd.read_csv(sys.argv[2], sep="\t", dtype=str).set_index("sample_id")
meta = pd.DataFrame({"condition": samples["condition"]})
dds = DeseqDataSet(counts=counts[meta.index].T.round().astype(int), metadata=meta, design="~condition", quiet=True, n_cpus=2)
dds.deseq2()
st = DeseqStats(dds, contrast=["condition", "treated", "vehicle"], alpha=0.05, quiet=True, n_cpus=2)
st.summary()
r = st.results_df.sort_values(["padj", "pvalue"])
r.index.name = "gene_id"
r[["log2FoldChange", "pvalue", "padj"]].to_csv(sys.argv[3], sep="\t")
