# build_deg_matrices.R
# Builds the DEG matrices of the cepip example data (pvm3_deg_NeoC.rds and
# pvm3_deg_NonNeoC.rds) from the SEA-AD Part 8 export of the PVM_3 arm.
#
# Usage: Rscript build_deg_matrices.R <R_export_CP10K/PVM3_all folder> <output folder>
#
# Input per region set (neocortex, nonneocortex), in deg_support/:
#   norm.mtx (cells x genes, ln(CP10K + 1)), norm_barcodes.tsv, norm_features.tsv
# Output: a sparse genes x cells matrix (dgCMatrix) with gene symbols as row
# names and barcodes as column names. Genes expressed in fewer than 40 cells are
# dropped. The smallest DEG group of the paper has 970 cells, and a gene needs
# 5% of a group's cells (at least 49) to be tested, so no tested gene is lost.

suppressPackageStartupMessages(library(Matrix))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript build_deg_matrices.R <PVM3_all folder> <output folder>")
dir.create(args[2], showWarnings = FALSE, recursive = TRUE)
for (set in c(neocortex = "NeoC", nonneocortex = "NonNeoC")) {
    sub <- names(which(c(neocortex = "NeoC", nonneocortex = "NonNeoC") == set))
    d <- file.path(args[1], sub, "deg_support")
    X <- readMM(file.path(d, "norm.mtx"))
    bc <- readLines(file.path(d, "norm_barcodes.tsv"), warn = FALSE)
    ft <- read.table(file.path(d, "norm_features.tsv"), sep = "\t", header = FALSE,
                     stringsAsFactors = FALSE)[, 2]
    X <- t(as(X, "CsparseMatrix"))
    rownames(X) <- make.unique(ft)
    colnames(X) <- bc
    keep <- Matrix::rowSums(X > 0) >= 40
    X <- X[keep, , drop = FALSE]
    out <- file.path(args[2], paste0("pvm3_deg_", set, ".rds"))
    saveRDS(X, out, compress = "xz")
    cat(sprintf("%s: %d genes x %d cells, %d non-zero values -> %s (%.0f MB)\n", set, nrow(X),
                ncol(X), length(X@x), out, file.size(out) / 1e6))
}
