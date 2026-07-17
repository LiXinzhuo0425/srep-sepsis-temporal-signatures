#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: preprocess_gse110487.R <grouped_counts.tsv> <required_genes.txt> <output_vst.tsv>")
}

project_lib <- "/Users/felix/Documents/New project/longitudinal_stage3/03_00_09_environment_lock/R_libs"
.libPaths(c(project_lib, .libPaths()))
suppressPackageStartupMessages(library(DESeq2))

counts_path <- args[[1]]
required_path <- args[[2]]
output_path <- args[[3]]

counts <- read.delim(counts_path, row.names = 1, check.names = FALSE)
counts <- as.matrix(counts)
storage.mode(counts) <- "integer"
if (any(counts < 0L) || anyNA(counts)) stop("Counts must be nonnegative integers without missing values")
if (anyDuplicated(rownames(counts))) stop("Grouped count row names are not unique")

coldata <- data.frame(row.names = colnames(counts), intercept = factor(rep("all", ncol(counts))))
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ 1)
dds <- estimateSizeFactors(dds)
vst_obj <- varianceStabilizingTransformation(dds, blind = TRUE)
vst_matrix <- assay(vst_obj)

required <- scan(required_path, what = character(), quiet = TRUE)
missing <- setdiff(required, rownames(vst_matrix))
if (length(missing) > 0) stop(paste("Missing required genes after grouping:", paste(missing, collapse = ",")))
out <- vst_matrix[required, , drop = FALSE]
write.table(cbind(gene = rownames(out), as.data.frame(out, check.names = FALSE)), output_path,
            sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("DESeq2=%s rows=%d samples=%d required=%d\n",
            as.character(packageVersion("DESeq2")), nrow(counts), ncol(counts), nrow(out)))
