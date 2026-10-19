# reproduce_seaad.R
# Runs cepip on the SEA-AD example data and compares every result with the
# published values (inst/reference/seaad_reference.csv).
#
# Arms:
#   PVM_3        MT-CO3 against the mean of its eight mitochondrial partners,
#                all Micro-PVM_3-SEAAD cells with MT-CO3 > 0 (the paper's arm)
#   HKG8_control ACTB against the mean of eight housekeeping genes, on the
#                paper's housekeeping-gene control cells
#   HKG_control  MT-CO3 against the mean of the nine housekeeping genes, on the
#                same cells as PVM_3 (same-cell control; GAM step only)
#
# Usage (from a shell):
#   Rscript reproduce_seaad.R <output folder> [workers] [data folder]
# Without a data folder, the files are downloaded with cepip_example().
# The DEG step runs only when the DEG matrices are available.
# Each arm gets one folder with its tables, plots and console log. The check
# against the published values is written to reproduction_check.csv and
# reproduction_check_log.txt in the output folder.

suppressPackageStartupMessages(library(cepip))
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else "cepip_seaad"
workers <- if (length(args) >= 2) as.integer(args[2]) else 1L
data_dir <- if (length(args) >= 3) args[3] else NULL
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

load_file <- function(which, file) {
    if (!is.null(data_dir)) {
        p <- file.path(data_dir, file)
        if (!file.exists(p)) return(NULL)
        if (grepl("rds$", p)) return(readRDS(p))
        return(utils::read.csv(p, check.names = FALSE, stringsAsFactors = FALSE))
    }
    tryCatch(cepip_example(which), error = function(e) NULL)
}
pvm3 <- load_file("pvm3", "pvm3_cells.csv")
hkg8 <- load_file("hkg8", "hkg8_cells.csv")
donors <- load_file("donors", "donors.csv")
if (is.null(pvm3) || is.null(hkg8) || is.null(donors))
    stop("The example data could not be loaded.")

MT8 <- c("MT-ATP6", "MT-CO2", "MT-ND3", "MT-ND4", "MT-CYB", "MT-ND1", "MT-ND2", "MT-CO1")
HK8 <- c("EEF1A1", "RPL13A", "B2M", "YWHAZ", "FTL", "PSAP", "RAB7A", "HNRNPA2B1")
HK9 <- c("ACTB", HK8)
REGIONS <- list(NeoC = c("DFC", "MTG", "AnG", "STG", "V1C", "ITG"),
                NonNeoC = c("HIP", "MEC", "LEC", "FI"))
SETS <- c("NeoC", "NonNeoC")
OUTCOMES <- c("cognitive_status", "MMSE", "CASI", "CERAD", "Braak", "Thal", "ADNC", "CAA",
              "APOE_e4_dose", "age_at_death", "sex", "PMI", "RIN", "brain_pH",
              "fresh_brain_weight")
COV4 <- c("age_at_death", "APOE_e4_dose", "sex", "PMI")
LABEL <- c(cognitive_status = "Cognitive status", MMSE = "MMSE score", CASI = "CASI score",
           CERAD = "CERAD score", Braak = "Braak stage", Thal = "Thal phase", ADNC = "ADNC rank",
           CAA = "CAA score", APOE_e4_dose = "APOE e4 dose")

got <- list()
put <- function(arm, set, step, item, value)
    got[[length(got) + 1L]] <<- data.frame(arm = arm, region_set = set, step = step, item = item,
                                           package = as.numeric(value), stringsAsFactors = FALSE)

# Per-region and pooled GAM fits (the paper's Part 9 table).
region_fits <- function(d, x, y, arm, set) {
    for (r in c(REGIONS[[set]], "Pooled")) {
        dd <- if (r == "Pooled") d else d[d$region == r, ]
        if (nrow(dd) < 10) next
        f <- cepip_fit(dd[[x]], y[match(dd$barcode, d$barcode)], donor = dd$donor,
                       workers = workers)
        put(arm, set, "01_fit", paste0("k_", r), f$k)
        put(arm, set, "01_fit", paste0("DE_", r), f$de)
    }
}

record <- function(res, arm, set) {
    ip <- res$ip
    v <- c(IP_Value = ip$ip_x, IP_CI_Lower = ip$ip_ci[1], IP_CI_Upper = ip$ip_ci[2],
           IP_Y_Value = ip$ip_y, IP_Y_CI_Lower = ip$ip_y_ci[1], IP_Y_CI_Upper = ip$ip_y_ci[2],
           unlist(ip$iprs[c("C1", "C2", "C3", "C4", "C5", "IPRS")]))
    for (nm in names(v)) put(arm, set, "02_ip", nm, v[[nm]])
    if (!is.null(res$skcs)) {
        put(arm, set, "03_skcs", "n_core", res$skcs$n_core)
        put(arm, set, "03_skcs", "n_TREP", sum(res$skcs$ep$is_TREP))
    }
    if (!is.null(res$cmrf)) {
        put(arm, set, "04_band", "delta_base", res$band$delta_base)
        put(arm, set, "04_band", "IP_u", res$band$IP_u)
        put(arm, set, "04_band", "n_candidates", sum(res$band$cells$is_candidate))
        for (cc in names(res$cmrf$counts)) put(arm, set, "05_cmrf", paste0("n_", cc), res$cmrf$counts[[cc]])
        put(arm, set, "05_cmrf", "n_removed", res$cmrf$gate$n_removed)
        put(arm, set, "05_cmrf", "p_final", res$cmrf$gate$p_final)
        put(arm, set, "05_cmrf", "dRMSE_final", res$cmrf$gate$effect_final)
        put(arm, set, "05_cmrf", "BV", res$cmrf$gate$bv)
    }
    if (!is.null(res$validate)) {
        put(arm, set, "06_validate", "k", res$validate$fit$k)
        put(arm, set, "06_validate", "DE", res$validate$de)
        put(arm, set, "06_validate", "n_TREP", res$validate$n_trep)
        put(arm, set, "06_validate", "agreement_pct", 100 * res$validate$agreement)
    }
}

donor_and_ols <- function(res, arm, set) {
    nt <- stats::setNames(donors[[paste0("n_PVM3_cells_", set)]], donors$donor)
    nt <- nt[!is.na(nt)]
    dp <- cepip_donor(res, n_total = nt)
    for (i in seq_len(nrow(dp))) put(arm, set, "07_donor", paste0("pct_PostIP_", dp$donor[i]),
                                     dp$pct_PostIP[i])
    tab <- merge(dp, donors, by = "donor")
    as <- cepip_associate(tab, "pct_PostIP", OUTCOMES, COV4)
    for (i in which(as$outcome %in% names(LABEL))) for (v in c("beta_std", "p_value", "FDR", "delta_R2"))
        put(arm, set, "09_associate", paste0(v, "_", LABEL[[as$outcome[i]]]), as[[v]][i])
    list(donor = dp, associate = as, outcomes = tab)
}

run_arm <- function(arm, cells, x, partners) {
    log_file <- file.path(tempdir(), paste0(arm, "_log.txt"))
    zz <- file(log_file, open = "wt")
    sink(zz); sink(zz, type = "message")
    on.exit({ sink(type = "message"); sink(); close(zz) }, add = TRUE)
    results <- list(); extras <- list()
    for (set in SETS) {
        d <- cells[cells$region_set == set, ]
        y <- cepip_composite(d, partners)
        cat("\n==", arm, set, ":", nrow(d), "cells\n")
        region_fits(d, x, y, arm, set)
        res <- cepip(d[[x]], y, donor = d$donor, id = d$barcode, cmrf = TRUE, ip_k = 10,
                     y_spine = d$composite_part8, workers = workers)
        record(res, arm, set)
        ex <- donor_and_ols(res, arm, set)
        deg_m <- if (arm == "PVM_3") load_file(paste0("deg_", set), paste0("pvm3_deg_", set, ".rds"))
        if (!is.null(deg_m)) {
            ex$deg <- cepip_deg(res, deg_m, verbose = TRUE)
            for (ct in names(ex$deg)) {
                info <- attr(ex$deg, "info")[[ct]]
                r <- ex$deg[[ct]]
                put(arm, set, "08_deg", paste0("n_cells_group_1_", ct), info$n_fg)
                put(arm, set, "08_deg", paste0("n_cells_group_2_", ct), info$n_bg)
                if (is.null(r)) next
                put(arm, set, "08_deg", paste0("n_paired_donors_", ct), info$n_donors_paired)
                put(arm, set, "08_deg", paste0("n_genes_tested_", ct), nrow(r))
                put(arm, set, "08_deg", paste0("n_significant_", ct), sum(r$significant))
                put(arm, set, "08_deg", paste0("n_significant_up_", ct), sum(r$significant & r$avg_log2FC > 0))
                put(arm, set, "08_deg", paste0("n_significant_down_", ct), sum(r$significant & r$avg_log2FC < 0))
            }
        }
        results[[set]] <- res; extras[[set]] <- ex
    }
    sink(type = "message"); sink(); close(zz); on.exit()
    cepip_save(results, file.path(out_dir, arm), extras = extras, log = readLines(log_file))
    invisible(results)
}

t0 <- Sys.time()
run_arm("PVM_3", pvm3, "MT-CO3", MT8)
run_arm("HKG8_control", hkg8, "ACTB", HK8)

# Same-cell control: GAM fits only (the paper reports this comparison at the GAM step).
hk_results <- list()
for (set in SETS) {
    d <- pvm3[pvm3$region_set == set, ]
    y <- cepip_composite(d, HK9)
    region_fits(d, "MT-CO3", y, "HKG_control", set)
    hk_results[[set]] <- cepip(d[["MT-CO3"]], y, donor = d$donor, id = d$barcode,
                               workers = workers, verbose = FALSE)
}
cepip_save(hk_results, file.path(out_dir, "HKG_control"))

# Comparison with the published values.
ref <- utils::read.csv(system.file("reference", "seaad_reference.csv", package = "cepip"),
                       stringsAsFactors = FALSE)
ref$arm[ref$arm == "HKG8"] <- "HKG8_control"
pk <- do.call(rbind, got)
cmp <- merge(ref, pk, by = c("arm", "region_set", "step", "item"), all.x = TRUE)
cmp$difference <- cmp$package - cmp$value
cmp$result <- ifelse(is.na(cmp$package), "NOT RUN",
                     ifelse(abs(cmp$difference) <= cmp$tolerance, "PASS", "FAIL"))
cmp <- cmp[order(cmp$arm, cmp$region_set, cmp$step, cmp$item), ]
utils::write.csv(cmp, file.path(out_dir, "reproduction_check.csv"), row.names = FALSE)
summ <- as.data.frame.matrix(table(paste(cmp$arm, cmp$step), cmp$result))
lines <- c(sprintf("cepip %s reproduction check, %s, %.0f min", utils::packageVersion("cepip"),
                   format(Sys.time()), as.numeric(difftime(Sys.time(), t0, units = "mins"))),
           utils::capture.output(print(summ)), "",
           "Values outside their tolerance:",
           utils::capture.output(print(cmp[cmp$result == "FAIL",
                                           c("arm", "region_set", "item", "value", "package")])),
           "", sprintf("VERDICT: %s", if (any(cmp$result == "FAIL")) "FAIL" else "PASS"))
writeLines(lines, file.path(out_dir, "reproduction_check_log.txt"))
cat(lines, sep = "\n")
