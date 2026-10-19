#' CMRF tolerance band (Stage 1)
#'
#' Builds the CMRF Stage 1 band in copula space. This step is fast, so it can
#' be run several times to choose `kappa` and `bw` (see [plot_band()]) before
#' running the slower Stage 2 gate in [cepip_cmrf()]. It is optional:
#' [cepip_cmrf()] builds the band itself with `kappa = 0.5` and `bw = 0.1`
#' when no band is given.
#'
#' * Copula coordinates over all cells: u = rank(x) / (n + 1) and
#'   v = rank(y) / (n + 1), with average ranks for ties.
#' * IP_u = (number of cells with x <= IPx) / (n + 1) and IP_v = (number of cells
#'   with y <= IPy) / (n + 1).
#' * The spine g(u) is a GAM, v ~ s(u, bs = "tp", k) (with the donor random
#'   intercept and inverse-donor weights when donors are known), fitted by
#'   REML (`select = TRUE`, `gamma = 1.5`) on the SKCS-core TREP cells, with k
#'   chosen by PRSS as in [cepip_fit()]. The donor term is excluded when
#'   predicting.
#' * delta_base is the 95% (`coverage`) donor-weighted quantile of
#'   |v - g(u)| over all SKCS-core cells.
#' * delta_adapt(u) = delta_base x (1 + kappa x K(u)), with the Gaussian kernel
#'   K(u) = exp(-(u - IP_u)^2 / (2 bw^2)), which is 1 at IP_u. The band is
#'   therefore widest at the IP, and at IP_u +- bw the widening term has decayed
#'   to exp(-0.5) = 60.65% of its peak.
#' * A cell passes Stage 1 when |v - g(u)| <= delta_adapt(u). SKCS-core cells
#'   that pass are the reference set (Original + Refined, OriRef); non-SKCS
#'   cells that pass are the Stage 2 candidates; SKCS-core cells that fail are
#'   "Original only" and leave the Module.
#'
#' @param skcs A `cepip_skcs` object (with its SKCS-core fit).
#' @param ip A `cepip_ip` object of the same cells.
#' @param kappa Height of the widening at IP_u, as a fraction of delta_base.
#' @param bw Bandwidth (standard deviation) of the Gaussian kernel, in u units.
#' @param coverage Quantile of |v - g(u)| that sets delta_base.
#' @param spine_k `"prss"` or a fixed k for the spine.
#' @param min_spine_cells Minimum number of SKCS-core TREP cells for the spine.
#' @param min_donors Minimum number of donors for the donor random intercept.
#' @param y_spine Optional composite values (one per cell, same order) whose
#'   ranks give the v coordinates of the spine fit. By default the composite of
#'   `skcs` is used. The published SEA-AD analysis fitted the spine on the
#'   composite column of its data export, which is stored in single precision;
#'   passing that column reproduces the published spine exactly.
#' @param workers Number of parallel workers for the spine k search.
#' @param verbose Print progress messages.
#' @return An object of class `cepip_band`.
#' @examples
#' \donttest{
#' d <- cepip_toy
#' s <- cepip_skcs(d$x, d$y, donor = d$donor, id = d$id)
#' ip <- cepip_ip(d$x, d$y, donor = d$donor, id = d$id)
#' b <- cepip_band(s, ip)
#' b$delta_base
#' table(b$cells$is_candidate)
#' }
#' @export
cepip_band <- function(skcs, ip, kappa = 0.5, bw = 0.1, coverage = 0.95,
                       spine_k = "prss", min_spine_cells = 50L, min_donors = 3L,
                       y_spine = NULL, workers = 1L, verbose = FALSE) {
    check_class(skcs, "cepip_skcs", "skcs")
    check_class(ip, "cepip_ip", "ip")
    if (!isTRUE(ip$has_ip))
        stop("No IP was found (", ip$status, "), so the CMRF band cannot be built.",
             call. = FALSE)
    if (is.null(skcs$ep))
        stop("skcs has no SKCS-core fit; run cepip_skcs() with fit = TRUE.", call. = FALSE)

    cells <- skcs$cells
    n <- nrow(cells)
    has_donor <- !all(is.na(cells$donor))
    u <- rank(cells$x, ties.method = "average") / (n + 1)
    v <- rank(cells$y, ties.method = "average") / (n + 1)
    IP_u <- sum(cells$x <= ip$ip_x) / (n + 1)
    IP_v <- sum(cells$y <= ip$ip_y) / (n + 1)

    core <- cells$is_core
    is_trep <- rep(FALSE, n)
    is_trep[core] <- skcs$ep$is_TREP
    trep <- which(is_trep)
    if (length(trep) < min_spine_cells)
        stop(sprintf("Only %d SKCS-core TREP cells; at least %d are needed for the spine.",
                     length(trep), min_spine_cells), call. = FALSE)

    v_fit <- v
    if (!is.null(y_spine)) {
        if (length(y_spine) != n) stop("y_spine must have one value per cell.", call. = FALSE)
        v_fit <- rank(as.numeric(y_spine), ties.method = "average") / (n + 1)
    }
    spine <- fit_spine(u[trep], v_fit[trep], if (has_donor) cells$donor[trep] else NULL,
                       k = spine_k, min_donors = min_donors, workers = workers,
                       verbose = verbose)
    g_u <- spine_predict(spine$model, u)

    w_core <- if (has_donor) donor_weights(cells$donor[core]) else rep(1, sum(core))
    delta_base <- weighted_quantile(abs(v[core] - g_u[core]), w_core, coverage)

    K <- exp(-(u - IP_u)^2 / (2 * bw^2))
    delta_adapt <- delta_base * (1 + kappa * K)
    resid_abs <- abs(v - g_u)
    pass1 <- resid_abs <= delta_adapt

    # Four-group labels of the SKCS-core cells: IP side (axis of the IP step) and TREP.
    ip_side <- ip_labels(cells$x, cells$y, ip, ip$ip_axis)
    group_orig <- ifelse(core, paste(ip_side, ifelse(is_trep, "TREP", "non-TREP")),
                         "Non-selected")

    out_cells <- data.frame(
        id = cells$id, x = cells$x, y = cells$y, donor = cells$donor,
        Si_norm = cells$Si_norm, Ki = cells$Ki, is_core = core, is_TREP = is_trep,
        ip_side = ip_side, group_orig = group_orig,
        u = u, v = v, g_u = g_u, resid_abs = resid_abs, K = K,
        delta_adapt = delta_adapt, pass_stage1 = pass1,
        is_anchor = core & pass1, is_candidate = !core & pass1,
        is_original_only = core & !pass1, stringsAsFactors = FALSE)
    msg(verbose, "  Band: delta_base = %.6f | OriRef %d | candidates %d | Original only %d",
        delta_base, sum(out_cells$is_anchor), sum(out_cells$is_candidate),
        sum(out_cells$is_original_only))

    ug <- seq(0, 1, length.out = 1001L)
    gg <- spine_predict(spine$model, ug)
    structure(list(
        cells = out_cells, delta_base = delta_base, kappa = kappa, bw = bw,
        coverage = coverage, IP_u = IP_u, IP_v = IP_v, IP_x = ip$ip_x, IP_y = ip$ip_y,
        ip_axis = ip$ip_axis, spine = spine,
        spine_grid = data.frame(u = ug, g = gg),
        spine_monotone = all(diff(gg[seq(1L, 1001L, by = 2L)]) >= -1e-8),
        n_spine_cells = length(trep),
        qc = band_qc(delta_base, kappa, bw)),
        class = "cepip_band")
}

# Spine GAM v ~ s(u) (+ donor random intercept), k by PRSS as in fit_frame().
# A failed fit is retried at k - 1, as in the published pipeline.
fit_spine <- function(u, v, donor = NULL, k = "prss", k_max = 10L, gamma = 1.5,
                      select = TRUE, min_donors = 3L, workers = 1L, verbose = FALSE) {
    d <- data.frame(u_trep = as.numeric(u), v_resp = as.numeric(v))
    if (!is.null(donor)) {
        d$donor <- as_donor_factor(donor)
        d$w <- donor_weights(d$donor)
    } else d$w <- rep(1, nrow(d))
    use_re <- use_donor_re(d, min_donors)
    max_k <- max(3L, min(as.integer(k_max), length(unique(d$u_trep)) - 1L))
    ks <- if (identical(k, "prss")) seq.int(3L, max_k) else as.integer(k)
    args <- list(d = d, use_re = use_re, gamma = gamma, select = select, resp = "v_resp",
                 var = "u_trep", linear = FALSE, retry_lower = TRUE)
    res <- if (workers > 1L && length(ks) > 1L) {
        cl <- parallel::makePSOCKcluster(min(as.integer(workers), length(ks)))
        on.exit(parallel::stopCluster(cl), add = TRUE)
        do.call(parallel::parLapply, c(list(cl = cl, X = ks, fun = fit_k_worker), args))
    } else lapply(ks, function(kk) do.call(fit_k_worker, c(list(k = kk), args)))
    ok <- !vapply(res, is.null, logical(1))
    if (!any(ok)) stop("All spine fits failed.", call. = FALSE)
    tab <- data.frame(k = ks, PRSS = vapply(res, function(r) if (is.null(r)) NA_real_ else r$PRSS, 0),
                      RSS = vapply(res, function(r) if (is.null(r)) NA_real_ else r$RSS, 0),
                      lambda = vapply(res, function(r) if (is.null(r)) NA_real_ else r$lambda, 0))
    j <- which.min(tab$PRSS)
    tab$selected <- seq_len(nrow(tab)) == j
    for (jj in which(ok)) msg(verbose, "  spine k = %d | PRSS = %.6f", tab$k[jj], tab$PRSS[jj])
    list(model = res[[j]]$model, k = ks[j], prss = tab, donor_re = use_re,
         n_donors = n_donor_levels(d$donor), n_cells = nrow(d))
}

spine_predict <- function(model, u) {
    nd <- data.frame(u_trep = as.numeric(u))
    if (model_has_re(model)) {
        ref <- model$var.summary[["donor"]]
        ref <- if (is.factor(ref)) levels(ref)[1] else as.character(ref)[1]
        nd$donor <- factor(ref, levels = levels(model$model$donor))
        return(as.numeric(stats::predict(model, newdata = nd, exclude = DONOR_RE_LABEL,
                                         type = "response")))
    }
    as.numeric(stats::predict(model, newdata = nd, type = "response"))
}

# QC metrics of the adaptive band (they follow from delta_base, kappa and bw):
# QC1 the full width at IP_u, QC2 the share of the widening left at IP_u +- bw,
# QC3 the largest slope of the kernel, reached at IP_u +- bw.
band_qc <- function(delta_base, kappa, bw) {
    list(width_at_ip = 2 * delta_base * (1 + kappa),
         width_base = 2 * delta_base,
         widening_left_at_bw = exp(-0.5),
         width_at_bw = 2 * delta_base * (1 + kappa * exp(-0.5)),
         max_kernel_slope = 1 / (bw * sqrt(exp(1))),
         active_halfwidth = bw * sqrt(-2 * log(0.02)))
}
