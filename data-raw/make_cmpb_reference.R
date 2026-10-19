# Builds inst/reference/cmpb_pca_reference.csv from the published Supplementary
# Tables of the CMPB CEP-IP paper (Zenodo record 10.5281/zenodo.19360154):
#   Supplementary Table 10 (Deviance_Components sheet): cells, total EDF and DE
#   of the TRPM4-Ribo GAM of each PCa tumor sample.
#   Supplementary Table 15 (TRPM4 sheet): IP, CI, TREP cells, stat test, IPRS
#   and its components of each PCa tumor sample.
# Usage: Rscript data-raw/make_cmpb_reference.R <ST10 .xlsx> <ST15 .xlsx>
args <- commandArgs(trailingOnly = TRUE)
st10 <- args[1]; st15 <- args[2]
samples <- c("HYW_4701_Tumor", "HYW_4847_Tumor", "HYW_4880_Tumor", "HYW_4881_Tumor",
             "HYW_5386_Tumor", "HYW_5742_Tumor", "HYW_5755_Tumor")
rows <- list()
add <- function(set, step, item, value, tol)
    rows[[length(rows) + 1L]] <<- data.frame(arm = "PCa_TRPM4_Ribo", region_set = set,
                                            step = step, item = item, value = as.numeric(value),
                                            tolerance = tol, stringsAsFactors = FALSE)

dc <- openxlsx::read.xlsx(st10, "Deviance_Components", startRow = 3, sep.names = " ")
dc <- dc[dc[["Gene set"]] == "Ribo" & dc[["Group"]] == "PCa", ]
stopifnot(nrow(dc) == 7L, identical(dc$Sample, sprintf("Pt.%d_Tumor", 1:7)))

x <- openxlsx::read.xlsx(st15, "TRPM4", colNames = FALSE)
h <- unlist(x[3, ])
v <- x[4:10, ]; names(v) <- h
stopifnot(identical(gsub(" ", "_", v$Sample), samples))
num <- function(col) as.numeric(v[[col]])

for (i in seq_along(samples)) {
    s <- samples[i]
    add(s, "01_fit", "cells", dc[["Number of cells"]][i], 0)
    add(s, "01_fit", "EDF_total", dc[["EDF (total)"]][i], 1e-6)
    add(s, "01_fit", "DE", dc[["DE (mgcv)"]][i], 1e-9)
    add(s, "02_ep", "TREP_cells", num("TREP cells (n)")[i], 0)
    add(s, "03_ip", "IP_x", num("IP Value")[i], 1e-6)
    add(s, "03_ip", "IP_CI_lower", num("CI Lower 95")[i], 1e-6)
    add(s, "03_ip", "IP_CI_upper", num("CI Upper 95")[i], 1e-6)
    add(s, "03_ip", "CI_width", num("CI Width")[i], 1e-6)
    add(s, "03_ip", "TRPM4_range", num("TRPM4 Range")[i], 1e-9)
    add(s, "03_ip", "neg_pos_crossings", num("Neg-Pos Crossings")[i], 0)
    # p values reach 1e-58, so the check is on -log10(p)
    add(s, "03_ip", "test_neg_log10_p", -log10(num("p-value")[i]), 1e-5)
    add(s, "04_iprs", "C1", num("C1 Precision")[i], 1e-6)
    add(s, "04_iprs", "C2", num("C2 Statistical")[i], 1e-6)
    add(s, "04_iprs", "C3", num("C3 Abruptness")[i], 1e-6)
    add(s, "04_iprs", "C4", num("C4 GAM Ribbon")[i], 1e-6)
    add(s, "04_iprs", "C5", num("C5 Smoothness")[i], 1e-6)
    add(s, "04_iprs", "Rx_amplitude", num("Rx Amplitude")[i], 1e-6)
    add(s, "04_iprs", "Rx_total_variation", num("Rx Total Variance")[i], 1e-6)
    add(s, "04_iprs", "IPRS", num("IPRS")[i], 1e-6)
}
ref <- do.call(rbind, rows)
utils::write.csv(ref, file.path("inst", "reference", "cmpb_pca_reference.csv"), row.names = FALSE)
cat(nrow(ref), "published values written\n")
