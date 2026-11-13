#!/usr/bin/env Rscript
# Development-only cross-check (not shipped in the task): the real Bioconductor methods on a
# gene x sample count table. usage: de_methods.R counts.tsv meta.tsv <deseq2|voom|edger> <yes|no> out.tsv
suppressMessages({library(edgeR); library(limma); library(DESeq2)})
a <- commandArgs(trailingOnly = TRUE)
counts <- read.delim(a[1], row.names = 1, check.names = FALSE)
meta <- read.delim(a[2], row.names = 1, check.names = FALSE, colClasses = "character")
method <- a[3]; use_batch <- a[4] == "yes"
counts <- round(as.matrix(counts[, rownames(meta)]))
meta$condition <- factor(meta$condition, levels = c("vehicle", "treated"))
meta$batch <- factor(paste0("b", meta$batch))
use_batch <- use_batch && nlevels(meta$batch) > 1
design <- if (use_batch) model.matrix(~ batch + condition, meta) else model.matrix(~ condition, meta)
if (method == "deseq2") {
  f <- if (use_batch) ~ batch + condition else ~ condition
  dds <- suppressMessages(DESeq(DESeqDataSetFromMatrix(counts, meta, f), quiet = TRUE))
  r <- results(dds, contrast = c("condition", "treated", "vehicle"), alpha = 0.05)
  res <- data.frame(gene_id = rownames(r), log2fc = r$log2FoldChange, padj = r$padj)
} else {
  y <- DGEList(counts)
  keep <- filterByExpr(y, design)
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  if (method == "voom") {
    fit <- eBayes(lmFit(voom(y, design), design))
    tt <- topTable(fit, coef = "conditiontreated", number = Inf, sort.by = "none")
    res <- data.frame(gene_id = rownames(tt), log2fc = tt$logFC, padj = tt$adj.P.Val)
  } else {
    y <- estimateDisp(y, design)
    tt <- topTags(glmQLFTest(glmQLFit(y, design), coef = "conditiontreated"), n = Inf, sort.by = "none")$table
    res <- data.frame(gene_id = rownames(tt), log2fc = tt$logFC, padj = tt$FDR)
  }
}
all <- data.frame(gene_id = rownames(counts))
res <- merge(all, res, by = "gene_id", all.x = TRUE)
res$log2fc[is.na(res$log2fc)] <- 0
res$significant <- as.integer(!is.na(res$padj) & res$padj < 0.05)
write.table(res, a[5], sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
