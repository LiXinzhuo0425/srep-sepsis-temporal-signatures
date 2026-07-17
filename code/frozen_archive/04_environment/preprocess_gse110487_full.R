#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: preprocess_gse110487_full.R <grouped_counts.tsv> <output_vst.tsv>")
}

project_lib <- "/Users/felix/Documents/New project/longitudinal_stage3/03_00_09_environment_lock/R_libs"
.libPaths(c(project_lib, .libPaths()))
suppressPackageStartupMessages(library(DESeq2))

counts <- read.delim(args[[1]], row.names = 1, check.names = FALSE)
counts <- as.matrix(counts)
storage.mode(counts) <- "integer"
if (any(counts < 0L) || anyNA(counts)) stop("Counts must be nonnegative integers without missing values")
if (anyDuplicated(rownames(counts))) stop("Grouped count row names are not unique")

keep <- rowSums(counts) > 0L
counts <- counts[keep, , drop = FALSE]
coldata <- data.frame(row.names = colnames(counts), intercept = factor(rep("all", ncol(counts))))
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ 1)
dds <- estimateSizeFactors(dds)
vst_obj <- varianceStabilizingTransformation(dds, blind = TRUE)
out <- assay(vst_obj)

write.table(cbind(gene = rownames(out), as.data.frame(out, check.names = FALSE)), args[[2]],
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("DESeq2=%s genes=%d samples=%d\n",
            as.character(packageVersion("DESeq2")), nrow(out), ncol(out)))
