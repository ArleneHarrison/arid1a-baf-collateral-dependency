#!/usr/bin/env Rscript

# Reconstruct the ARID1A-low versus ARID1A-high tertile comparison from a raw
# gene-by-sample count matrix. The middle tertile is excluded by design.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript src/09_tcga_deseq2.R counts.csv all_results.csv significant.csv")
}

suppressPackageStartupMessages(library(DESeq2))

counts_df <- read.csv(args[[1]], check.names = FALSE, row.names = 1)
counts_mat <- round(as.matrix(counts_df))
storage.mode(counts_mat) <- "integer"
counts_mat <- counts_mat[rowSums(counts_mat) >= 10, , drop = FALSE]
if (!"ARID1A" %in% rownames(counts_mat)) stop("ARID1A is absent from the count matrix")

coldata <- data.frame(row.names = colnames(counts_mat), intercept = factor("all"))
dds0 <- DESeqDataSetFromMatrix(countData = counts_mat, colData = coldata, design = ~1)
dds0 <- estimateSizeFactors(dds0)
arid1a <- counts(dds0, normalized = TRUE)["ARID1A", ]
cutpoints <- quantile(arid1a, probs = c(1/3, 2/3), names = FALSE)
group <- rep(NA_character_, length(arid1a))
group[arid1a <= cutpoints[[1]]] <- "low"
group[arid1a >= cutpoints[[2]]] <- "high"
keep <- !is.na(group)

coldata <- data.frame(row.names = colnames(counts_mat)[keep], group = factor(group[keep], levels = c("high", "low")))
dds <- DESeqDataSetFromMatrix(countData = counts_mat[, keep, drop = FALSE], colData = coldata, design = ~group)
dds <- DESeq(dds, quiet = TRUE)
result <- as.data.frame(results(dds, contrast = c("group", "low", "high")))
result$gene <- rownames(result)
result <- result[, c("gene", setdiff(colnames(result), "gene"))]
write.csv(result, args[[2]], row.names = FALSE)
significant <- subset(result, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(significant, args[[3]], row.names = FALSE)
message(sprintf("Wrote %d genes; %d passed padj < 0.05 and |log2FC| > 1", nrow(result), nrow(significant)))

