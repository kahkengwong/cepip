# reproduce_cmpb_pca.R
# Runs cepip on the prostate cancer (PCa) dataset of the published CEP-IP paper
# (Wong, Computer Methods and Programs in Biomedicine 2026) and compares the
# results with the published values (inst/reference/cmpb_pca_reference.csv,
# taken from Supplementary Tables 10 and 15 of that paper).
#
# Analysis: TRPM4 against the mean of seven ribosomal genes, in the cells of
# clusters 6, 9, 11, 14 and 19 of each of the seven tumor samples. Each sample
# is analysed on its own, without donors, with the published IP split on the
# TRPM4 axis (ip_axis = "x"). These are the settings of the published analysis.
#
# Data: GSE185344_Seurat_processed.RData (9.5 GB) from the Hugging Face dataset
# kahkengwong/CEP-IP_Framework. Reading it needs the SeuratObject package. The
# GAM data are built as in Part_3.03_GAM_REML_PRSS_Setup.r of the GitHub
# repository kahkengwong/CEP-IP_Framework:
#   - cells of the five clusters of the integrated object whose names contain
#     the sample id, in the column order of that object
#   - TRPM4 and the mean of the seven ribosomal genes from the "data" layer of
#     the RNA assay of the original object
#   - log2(value + 1) of both
# They are saved as pca_trpm4_ribo_cells.csv in the output folder, and a later
# run reads that file instead of the .RData.
#
# Usage (from a shell):
#   Rscript reproduce_cmpb_pca.R <output folder> [workers] [data folder] [ST10 .xlsx]
# With Supplementary Table 10 (.xlsx, from the Zenodo record
# 10.5281/zenodo.19360154), the EP, TRPM4, Ribo and fitted value of every cell
# are also compared.
#
# Output (all in the output folder):
#   PCa_TRPM4_Ribo/             tables and plots of the seven samples (cepip_save)
#   pca_trpm4_ribo_cells.csv    the GAM data of the seven samples
#   cmpb_pca_check.xlsx         the comparison with the published values
#   reproduction_check_log.txt  the summary and the verdict

suppressPackageStartupMessages(library(cepip))
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else "cepip_cmpb_pca"
workers <- if (length(args) >= 2) as.integer(args[2]) else 1L
data_dir <- if (length(args) >= 3) args[3] else "."
st10 <- if (length(args) >= 4) args[4] else NULL
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time()

SAMPLES <- c("HYW_4701_Tumor", "HYW_4847_Tumor", "HYW_4880_Tumor", "HYW_4881_Tumor",
             "HYW_5386_Tumor", "HYW_5742_Tumor", "HYW_5755_Tumor")
CLUSTERS <- c(6, 9, 11, 14, 19)
RIBO <- c("RPL10", "RPL27", "RPL28", "RPS2", "RPS8", "RPS12", "RPS26")
cells_csv <- file.path(out_dir, "pca_trpm4_ribo_cells.csv")

# ---- GAM data from the Seurat object ----
# The integrated object holds assay classes of the Seurat package. Its slots are
# read directly, so only SeuratObject is needed.
build_cells <- function(rdata) {
    if (!requireNamespace("SeuratObject", quietly = TRUE))
        stop("Reading the Seurat object needs SeuratObject: install.packages(\"SeuratObject\")",
             call. = FALSE)
    suppressPackageStartupMessages(requireNamespace("SeuratObject"))
    env <- new.env()
    cat("Loading", rdata, "\n")
    load(rdata, envir = env)
    integ <- if (is.list(env$prostate_results)) env$prostate_results$seurat_obj else NULL
    if (is.null(integ)) integ <- env$prostate_ca_seurat_integrated
    orig <- env$prostate_ca_seurat
    if (is.null(integ) || is.null(orig))
        stop("prostate_results$seurat_obj and prostate_ca_seurat were not found.", call. = FALSE)
    a <- orig@assays[["RNA"]]
    m <- tryCatch(SeuratObject::LayerData(a, layer = "data"), error = function(e) NULL)
    if (is.null(m)) m <- attr(a, "data")
    genes <- RIBO[RIBO %in% rownames(m)]
    integ_cells <- rownames(integ@meta.data)
    id <- integ@active.ident[integ_cells]
    cl_cells <- integ_cells[as.character(id) %in% as.character(CLUSTERS)]
    rows <- lapply(SAMPLES, function(s) {
        sel <- intersect(cl_cells, grep(s, integ_cells, value = TRUE))
        sub <- m[c("TRPM4", genes), sel, drop = FALSE]
        d <- log2(data.frame(TRPM4 = as.numeric(sub["TRPM4", ]),
                             Expression = Matrix::colMeans(sub[genes, , drop = FALSE])) + 1)
        data.frame(sample = s, cell = sel, TRPM4 = d$TRPM4, Expression = d$Expression,
                   stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
}
if (file.exists(cells_csv)) {
    cat("Reading", cells_csv, "\n")
    cells <- utils::read.csv(cells_csv, stringsAsFactors = FALSE)
} else {
    cells <- build_cells(file.path(data_dir, "GSE185344_Seurat_processed.RData"))
    utils::write.csv(cells, cells_csv, row.names = FALSE)
}

# ---- cepip on each sample ----
got <- list()
put <- function(set, step, item, value)
    got[[length(got) + 1L]] <<- data.frame(region_set = set, step = step, item = item,
                                           package = as.numeric(value), stringsAsFactors = FALSE)
log_lines <- character(0)
results <- list()
for (s in SAMPLES) {
    d <- cells[cells$sample == s, ]
    log_lines <- c(log_lines, paste("==", s), utils::capture.output(
        res <- cepip(d$TRPM4, d$Expression, id = d$cell, ip_axis = "x", workers = workers),
        type = "message"))
    results[[s]] <- res
    ip <- res$ip
    put(s, "01_fit", "cells", nrow(d))
    put(s, "01_fit", "EDF_total", sum(res$fit$model$edf))
    put(s, "01_fit", "DE", res$fit$de)
    put(s, "02_ep", "TREP_cells", sum(res$ep$is_TREP))
    put(s, "03_ip", "IP_x", ip$ip_x)
    put(s, "03_ip", "IP_CI_lower", ip$ip_ci[1])
    put(s, "03_ip", "IP_CI_upper", ip$ip_ci[2])
    put(s, "03_ip", "CI_width", ip$iprs$ci_width)
    put(s, "03_ip", "TRPM4_range", ip$iprs$x_range)
    put(s, "03_ip", "neg_pos_crossings", ip$n_crossings)
    put(s, "03_ip", "test_neg_log10_p", -log10(ip$test$p))
    for (cc in c("C1", "C2", "C3", "C4", "C5", "IPRS")) put(s, "04_iprs", cc, ip$iprs[[cc]])
    put(s, "04_iprs", "Rx_amplitude", ip$iprs$amplitude)
    put(s, "04_iprs", "Rx_total_variation", ip$iprs$total_variation)
}
cat(log_lines, sep = "\n")
cepip_save(results, file.path(out_dir, "PCa_TRPM4_Ribo"), log = log_lines)

# ---- Comparison with the published values ----
ref <- utils::read.csv(system.file("reference", "cmpb_pca_reference.csv", package = "cepip"),
                       stringsAsFactors = FALSE)
pk <- do.call(rbind, got)
cmp <- merge(ref, pk, by = c("region_set", "step", "item"), all.x = TRUE)
cmp$difference <- cmp$package - cmp$value
cmp$result <- ifelse(is.na(cmp$package), "NOT RUN",
                     ifelse(abs(cmp$difference) <= cmp$tolerance, "PASS", "FAIL"))
cmp <- cmp[order(cmp$region_set, cmp$step, cmp$item),
           c("arm", "region_set", "step", "item", "value", "package", "difference",
             "tolerance", "result")]
sheets <- list(published_values = cmp)

# Optional cell-level check against Supplementary Table 10. Both sides are sorted
# by EP. EP is compared as a relative difference, because cells with a value
# close to the sample mean have very large negative EP values.
if (!is.null(st10)) {
    ep_pub <- openxlsx::read.xlsx(st10, 1, startRow = 3, check.names = FALSE)
    names(ep_pub)[c(1, 3, 4, 6, 7, 11)] <- c("pt", "EP", "TREP", "TRPM4", "Ribo", "fitted")
    cl <- do.call(rbind, lapply(seq_along(SAMPLES), function(i) {
        e <- results[[SAMPLES[i]]]$ep
        e <- e[order(-e$EP), ]
        p <- ep_pub[ep_pub$pt == paste0("Pt.", i), ]
        p <- p[order(-as.numeric(p$EP)), ]
        same <- nrow(p) == nrow(e)
        data.frame(region_set = SAMPLES[i], cells = nrow(e), cells_published = nrow(p),
                   max_rel_diff_EP = if (same) max(abs(e$EP - as.numeric(p$EP)) /
                                                   pmax(1, abs(as.numeric(p$EP)))) else NA,
                   max_diff_TRPM4 = if (same) max(abs(e$x - as.numeric(p$TRPM4))) else NA,
                   max_diff_Ribo = if (same) max(abs(e$y - as.numeric(p$Ribo))) else NA,
                   max_diff_fitted = if (same) max(abs(e$fitted - as.numeric(p$fitted))) else NA,
                   same_TREP_cells = if (same) all(e$is_TREP == (p$TREP == "YES")) else FALSE,
                   stringsAsFactors = FALSE)
    }))
    cl$result <- ifelse(cl$cells == cl$cells_published & cl$max_rel_diff_EP < 1e-8 &
                            cl$max_diff_TRPM4 < 1e-10 & cl$max_diff_Ribo < 1e-10 &
                            cl$max_diff_fitted < 1e-8 & cl$same_TREP_cells, "PASS", "FAIL")
    sheets$cell_level_ST10 <- cl
}

wb <- openxlsx::createWorkbook()
bold <- openxlsx::createStyle(textDecoration = "bold")
for (nm in names(sheets)) {
    openxlsx::addWorksheet(wb, nm)
    openxlsx::writeData(wb, nm, sheets[[nm]], headerStyle = bold)
    openxlsx::freezePane(wb, nm, firstRow = TRUE)
}
openxlsx::saveWorkbook(wb, file.path(out_dir, "cmpb_pca_check.xlsx"), overwrite = TRUE)

fail <- any(cmp$result != "PASS") || (!is.null(sheets$cell_level_ST10) &&
                                      any(sheets$cell_level_ST10$result != "PASS"))
summ <- as.data.frame.matrix(table(cmp$step, cmp$result))
lines <- c(sprintf("cepip %s reproduction check of the CMPB PCa analysis, %s, %.1f min",
                   utils::packageVersion("cepip"), format(Sys.time()),
                   as.numeric(difftime(Sys.time(), t0, units = "mins"))),
           sprintf("mgcv %s | %s", utils::packageVersion("mgcv"), R.version.string), "",
           "Published values (7 tumor samples):",
           utils::capture.output(print(summ)), "",
           "Values outside their tolerance:",
           utils::capture.output(print(cmp[cmp$result != "PASS",
                                           c("region_set", "item", "value", "package")],
                                       row.names = FALSE)))
if (!is.null(sheets$cell_level_ST10))
    lines <- c(lines, "", "Cell-level check against Supplementary Table 10:",
               utils::capture.output(print(sheets$cell_level_ST10, digits = 3, row.names = FALSE)))
lines <- c(lines, "", sprintf("VERDICT: %s", if (fail) "FAIL" else "PASS"))
writeLines(lines, file.path(out_dir, "reproduction_check_log.txt"))
cat(lines, sep = "\n")
