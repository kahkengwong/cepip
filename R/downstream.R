#' Cell groups for downstream analyses
#'
#' Returns the four groups used by [cepip_donor()] and [cepip_deg()]:
#' Module Post-IP (`PostIP`), Module Pre-IP (`PreIP`), non-Module Post-IP
#' (`NMPost`) and non-Module Pre-IP (`NMPre`). With CMRF, the Module is
#' OriRef plus Recovered; without CMRF (the published CEP-IP), the TREP cells
#' take the place of the Module. The IP side follows the IP axis of the
#' analysis (by default y < IPy is Pre-IP).
#'
#' @param x A `cepip` object, a `cepip_cmrf` object, or a `cepip_ip` object.
#' @return A data frame with the cell id, donor, module flag, IP side and group.
#' @examples
#' d <- cepip_toy
#' res <- cepip(d$x, d$y, donor = d$donor, id = d$id, verbose = FALSE)
#' table(cepip_groups(res)$group)
#' @export
cepip_groups <- function(x) {
    if (inherits(x, "cepip")) x <- if (!is.null(x$cmrf)) x$cmrf else x$ip
    if (inherits(x, "cepip_cmrf")) {
        cl <- x$cells
        mod <- cl$in_module
        side <- cl$ip_side_all
        ids <- cl$id; donor <- cl$donor
    } else if (inherits(x, "cepip_ip")) {
        if (!isTRUE(x$has_ip)) stop("No IP was found, so there are no Pre-IP and Post-IP groups.",
                                    call. = FALSE)
        ids <- x$ep$id; donor <- x$ep$donor
        mod <- x$ep$is_TREP
        side <- x$labels
    } else stop("x must be a cepip, cepip_cmrf or cepip_ip object.", call. = FALSE)
    grp <- ifelse(mod, ifelse(side == "Post-IP", "PostIP", "PreIP"),
                  ifelse(side == "Post-IP", "NMPost", "NMPre"))
    data.frame(id = ids, donor = donor, in_module = mod, ip_side = side, group = grp,
               stringsAsFactors = FALSE)
}

#' Per-donor share of cells in each group
#'
#' Counts each donor's cells per group and divides by the donor's total cell
#' count. The total should include every cell of the donor in the population,
#' also cells left out of the analysis (for example cells with GOI = 0), so
#' `n_total` can be given separately. Donors in `n_total` without any cell in
#' a group get 0.
#'
#' @param x A `cepip`, `cepip_cmrf` or `cepip_ip` object.
#' @param n_total Optional named vector of total cells per donor. By default
#'   the number of analysed cells per donor.
#' @return A data frame with one row per donor: total cells, counts and
#'   percentages of PostIP, PreIP, NMPost and NMPre, and the Module percentage.
#' @examples
#' d <- cepip_toy
#' res <- cepip(d$x, d$y, donor = d$donor, id = d$id, verbose = FALSE)
#' head(cepip_donor(res))
#' @export
cepip_donor <- function(x, n_total = NULL) {
    g <- cepip_groups(x)
    if (all(is.na(g$donor))) stop("The analysis has no donor labels.", call. = FALSE)
    if (is.null(n_total)) n_total <- table(g$donor)
    tot <- stats::setNames(as.numeric(n_total), names(n_total))
    miss <- setdiff(unique(g$donor), names(tot))
    if (length(miss)) stop(length(miss), " donor(s) of the analysis are missing from n_total.",
                           call. = FALSE)
    out <- data.frame(donor = names(tot), n_total = unname(tot), stringsAsFactors = FALSE)
    for (grp in c("PostIP", "PreIP", "NMPost", "NMPre")) {
        cnt <- table(factor(g$donor[g$group == grp], levels = names(tot)))
        out[[paste0("n_", grp)]] <- as.numeric(cnt)
        out[[paste0("pct_", grp)]] <- 100 * as.numeric(cnt) / out$n_total
    }
    out$pct_module <- out$pct_PostIP + out$pct_PreIP
    out[order(out$donor), , drop = FALSE]
}

#' Differential expression between cell groups, paired by donor
#'
#' Compares two of the four groups of [cepip_groups()] (foreground against
#' background). These four contrasts are built in:
#' `PostIP_vs_PreIP`, `NMPost_vs_NMPre`, `PostIP_vs_NMPost` and
#' `PreIP_vs_NMPre`. Other pairs of the four groups can be given as
#' `"A_vs_B"`.
#'
#' Donor unit (default):
#' 1. Genes expressed (value > 0) in at least `min_pct` of the foreground or
#'    of the background cells are kept.
#' 2. Donors with at least `min_cells_donor` cells in both groups are paired;
#'    at least `min_donors` paired donors are needed.
#' 3. For each donor and gene, log2FC = log2(mean(expm1(x_fg)) + 1) -
#'    log2(mean(expm1(x_bg)) + 1), on log1p-normalized values.
#' 4. Genes that change in fewer than 5% of the paired donors are dropped. Each
#'    gene is tested by a Wilcoxon signed-rank test of its non-zero donor
#'    log2FCs against 0 (at least 3 donors); BH is applied across genes.
#'    avg_log2FC is the mean over all paired donors and r_rb the matched-pairs
#'    rank-biserial correlation.
#'
#' The cell unit compares the cells directly with a Wilcoxon rank-sum test.
#' A gene is significant when its BH-adjusted p < `fdr` and |avg_log2FC| >=
#' `lfc`. Which groups to compare, and the thresholds, depend on the question
#' and the data; the defaults are those of the SEA-AD analysis.
#'
#' @param x A `cepip`, `cepip_cmrf` or `cepip_ip` object.
#' @param expr Log1p-normalized expression, genes x cells, with cell ids as
#'   column names (the ids used in the analysis) and gene names as row names.
#' @param contrast One or more contrasts.
#' @param unit `"donor"` or `"cell"`.
#' @param min_pct,min_cells_donor,min_donors,fdr,lfc See Details.
#' @param min_cells Minimum number of cells per group.
#' @param verbose Print progress messages.
#' @return A named list with one data frame per contrast (`NULL` when a
#'   contrast cannot be tested), with the attribute `info` giving the cell and
#'   donor counts.
#' @examples
#' d <- cepip_toy
#' res <- cepip(d$x, d$y, donor = d$donor, id = d$id, verbose = FALSE)
#' expr <- t(as.matrix(d[, c("x", paste0("P", 1:8))]))
#' colnames(expr) <- d$id
#' deg <- cepip_deg(res, expr, contrast = "PostIP_vs_PreIP")
#' deg$PostIP_vs_PreIP
#' @export
cepip_deg <- function(x, expr, contrast = c("PostIP_vs_PreIP", "NMPost_vs_NMPre",
                                            "PostIP_vs_NMPost", "PreIP_vs_NMPre"),
                      unit = c("donor", "cell"), min_pct = 0.05, min_cells_donor = 3L,
                      min_donors = 5L, fdr = 0.05, lfc = 0.25, min_cells = 10L,
                      verbose = FALSE) {
    unit <- match.arg(unit)
    g <- cepip_groups(x)
    if (is.null(colnames(expr))) stop("expr needs cell ids as column names.", call. = FALSE)
    if (!inherits(expr, "Matrix")) expr <- Matrix::Matrix(as.matrix(expr), sparse = TRUE)
    expr <- methods::as(expr, "CsparseMatrix")
    miss <- sum(!g$id %in% colnames(expr))
    if (miss) message(miss, " cell(s) of the analysis are not in expr and are skipped.")
    out <- list(); info <- list()
    for (ct in contrast) {
        parts <- strsplit(ct, "_vs_", fixed = TRUE)[[1]]
        if (length(parts) != 2L || !all(parts %in% c("PostIP", "PreIP", "NMPost", "NMPre")))
            stop("Unknown contrast: ", ct, call. = FALSE)
        fg <- intersect(g$id[g$group == parts[1]], colnames(expr))
        bg <- intersect(g$id[g$group == parts[2]], colnames(expr))
        fg <- setdiff(fg, bg)
        info[[ct]] <- list(n_fg = length(fg), n_bg = length(bg))
        if (length(fg) < min_cells || length(bg) < min_cells) {
            msg(verbose, "  %s: skipped (fewer than %d cells in a group)", ct, min_cells)
            out[ct] <- list(NULL); next
        }
        cells <- c(fg, bg)
        is_fg <- c(rep(TRUE, length(fg)), rep(FALSE, length(bg)))
        mat <- expr[, cells, drop = FALSE]
        dv <- g$donor[match(cells, g$id)]
        pct_fg <- Matrix::rowMeans(mat[, is_fg, drop = FALSE] > 0)
        pct_bg <- Matrix::rowMeans(mat[, !is_fg, drop = FALSE] > 0)
        keep <- (pct_fg >= min_pct) | (pct_bg >= min_pct)
        mat <- mat[keep, , drop = FALSE]
        fc_all <- log2fc(mat[, is_fg, drop = FALSE], mat[, !is_fg, drop = FALSE])
        if (unit == "donor") {
            pb <- donor_deg(mat, is_fg, dv, min_cells_donor, min_donors)
            info[[ct]] <- c(info[[ct]], pb$info)
            if (is.null(pb$res)) { out[ct] <- list(NULL); next }
            res <- pb$res
            res$avg_log2FC <- round(res$avg_log2FC, 5)
        } else {
            res <- cell_deg(mat, is_fg)
            res$avg_log2FC <- round(as.numeric(fc_all[res$gene]), 5)
            res$n_donors <- NA_integer_
        }
        res$p_val_adj <- signif(stats::p.adjust(res$p_val, method = "BH"), 5)
        res$p_val <- signif(res$p_val, 5)
        res$pct_fg <- round(as.numeric(pct_fg[res$gene]), 4)
        res$pct_bg <- round(as.numeric(pct_bg[res$gene]), 4)
        res$significant <- !is.na(res$p_val_adj) & res$p_val_adj < fdr & abs(res$avg_log2FC) >= lfc
        res$direction <- ifelse(res$avg_log2FC > 0, paste0("up_in_", parts[1]),
                                paste0("up_in_", parts[2]))
        res$contrast <- ct
        res$unit <- unit
        res <- res[order(res$p_val_adj, -abs(res$avg_log2FC)), ]
        rownames(res) <- NULL
        msg(verbose, "  %s: %d of %d genes significant", ct, sum(res$significant), nrow(res))
        out[[ct]] <- res
    }
    attr(out, "info") <- info
    attr(out, "settings") <- list(unit = unit, min_pct = min_pct, min_cells_donor = min_cells_donor,
                                  min_donors = min_donors, fdr = fdr, lfc = lfc)
    class(out) <- "cepip_deg"
    out
}

log2fc <- function(fg, bg) {
    fg@x <- expm1(fg@x); bg@x <- expm1(bg@x)
    log2(Matrix::rowMeans(fg) + 1) - log2(Matrix::rowMeans(bg) + 1)
}

donor_deg <- function(mat, is_fg, donor, min_cells, min_donors, min_donor_pct = 0.05) {
    d <- as.character(donor)
    ok <- !is.na(d) & nzchar(d)
    nf <- table(d[ok & is_fg]); nb <- table(d[ok & !is_fg])
    all_d <- union(names(nf), names(nb))
    cf <- stats::setNames(as.integer(nf[all_d]), all_d); cf[is.na(cf)] <- 0L
    cb <- stats::setNames(as.integer(nb[all_d]), all_d); cb[is.na(cb)] <- 0L
    paired <- all_d[cf >= min_cells & cb >= min_cells]
    info <- list(n_donors_paired = length(paired), n_donors_total = length(all_d))
    if (length(paired) < min_donors) return(list(res = NULL, info = info))
    fc <- vapply(paired, function(dd) {
        log2fc(mat[, which(ok & is_fg & d == dd), drop = FALSE],
               mat[, which(ok & !is_fg & d == dd), drop = FALSE])
    }, numeric(nrow(mat)))
    if (is.null(dim(fc))) fc <- matrix(fc, nrow = nrow(mat))
    rownames(fc) <- rownames(mat)
    fc <- fc[rowMeans(fc != 0) >= min_donor_pct, , drop = FALSE]
    if (!nrow(fc)) return(list(res = NULL, info = info))
    st <- t(apply(fc, 1L, function(v) {
        v <- v[is.finite(v)]
        nz <- v[v != 0]
        if (length(nz) < 3L) return(c(NA_real_, NA_real_, NA_real_))
        p <- tryCatch(suppressWarnings(stats::wilcox.test(nz, mu = 0, exact = FALSE,
                                                          correct = TRUE))$p.value,
                      error = function(e) NA_real_)
        r <- rank(abs(nz))
        wp <- sum(r[nz > 0]); wm <- sum(r[nz < 0])
        c(p, if ((wp + wm) > 0) (wp - wm) / (wp + wm) else NA_real_, mean(v))
    }))
    list(res = data.frame(gene = rownames(fc), avg_log2FC = st[, 3], p_val = st[, 1],
                          r_rb = st[, 2], n_donors = ncol(fc), stringsAsFactors = FALSE),
         info = info, fc = fc)
}

# Wilcoxon rank-sum test per gene on the cells (normal approximation with
# continuity and tie correction, two-sided).
cell_deg <- function(mat, is_fg) {
    n1 <- sum(is_fg); n2 <- sum(!is_fg); n <- n1 + n2
    res <- t(apply(as.matrix(mat), 1L, function(v) {
        r <- rank(v)
        w <- sum(r[is_fg]) - n1 * (n1 + 1) / 2
        tt <- table(v)
        sig <- sqrt(n1 * n2 / 12 * ((n + 1) - sum(tt^3 - tt) / (n * (n - 1))))
        z <- w - n1 * n2 / 2
        z <- (z - sign(z) * 0.5) / sig
        c(2 * stats::pnorm(-abs(z)), w / (n1 * n2))
    }))
    data.frame(gene = rownames(mat), p_val = res[, 1], auc = res[, 2],
               r_rb = 2 * (res[, 2] - 0.5), stringsAsFactors = FALSE)
}

#' Donor-level association with outcomes
#'
#' For each outcome, fits `outcome ~ predictor + covariates` and
#' `outcome ~ covariates` by ordinary least squares on the donors with complete
#' data, after z-scoring every variable. It reports the standardized beta of the
#' predictor, its p value, and delta R2, the R2 of the full model minus that of
#' the covariate-only model. An outcome is never used as its own covariate.
#' The BH false discovery rate is applied across all outcomes passed.
#'
#' In the SEA-AD analysis the predictor was the Module Post-IP percentage per
#' donor, the covariates were age at death, sex, PMI and APOE e4 dose, and the
#' BH family was 15 outcomes (nine clinical and pathology outcomes plus six
#' control variables). Outcomes and covariates depend on the data available.
#'
#' @param data Donor table (one row per donor).
#' @param predictor Column name of the predictor.
#' @param outcomes Column names of the outcomes.
#' @param covariates Column names of the covariates.
#' @param min_donors Minimum number of complete donors for a model.
#' @return A data frame with one row per outcome.
#' @examples
#' d <- cepip_toy
#' res <- cepip(d$x, d$y, donor = d$donor, id = d$id, verbose = FALSE)
#' dt <- cepip_donor(res)
#' set.seed(1)
#' dt$age <- rnorm(nrow(dt), 80, 5)
#' dt$score <- 0.05 * dt$pct_PostIP + rnorm(nrow(dt))
#' cepip_associate(dt, "pct_PostIP", outcomes = "score", covariates = "age",
#'                 min_donors = 10)
#' @export
cepip_associate <- function(data, predictor, outcomes, covariates = character(0),
                            min_donors = 20L) {
    miss <- setdiff(c(predictor, outcomes, covariates), names(data))
    if (length(miss)) stop("Columns not found: ", paste(miss, collapse = ", "), call. = FALSE)
    zs <- function(v) {
        v <- as.numeric(v)
        s <- stats::sd(v)
        if (!is.finite(s) || s < 1e-12) v - mean(v) else (v - mean(v)) / s
    }
    rows <- lapply(outcomes, function(oc) {
        cov <- setdiff(covariates, oc)
        vars <- c(oc, predictor, cov)
        d <- data[stats::complete.cases(data[, vars, drop = FALSE]), vars, drop = FALSE]
        empty <- data.frame(outcome = oc, n = nrow(d), covariates = paste(cov, collapse = " + "),
                            beta_std = NA_real_, se = NA_real_, ci_low = NA_real_,
                            ci_high = NA_real_, p_value = NA_real_, R2_full = NA_real_,
                            R2_covariates = NA_real_, delta_R2 = NA_real_,
                            stringsAsFactors = FALSE)
        if (nrow(d) < min_donors) return(empty)
        z <- as.data.frame(lapply(d, zs))
        names(z) <- c("y", "x", if (length(cov)) paste0("c", seq_along(cov)))
        cn <- if (length(cov)) paste(paste0("c", seq_along(cov)), collapse = " + ") else "1"
        mf <- stats::lm(stats::as.formula(paste("y ~ x", if (length(cov)) paste("+", cn) else "")),
                        data = z)
        mc <- stats::lm(stats::as.formula(paste("y ~", cn)), data = z)
        co <- summary(mf)$coefficients
        ci <- stats::confint(mf)["x", ]
        empty$beta_std <- co["x", 1]; empty$se <- co["x", 2]
        empty$ci_low <- ci[1]; empty$ci_high <- ci[2]; empty$p_value <- co["x", 4]
        empty$R2_full <- summary(mf)$r.squared
        empty$R2_covariates <- summary(mc)$r.squared
        empty$delta_R2 <- empty$R2_full - empty$R2_covariates
        empty
    })
    out <- do.call(rbind, rows)
    out$FDR <- stats::p.adjust(out$p_value, method = "BH")
    out$predictor <- predictor
    out
}
