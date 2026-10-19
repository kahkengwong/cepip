#' Copula monotone recovery filter (CMRF)
#'
#' Runs CMRF Stage 2, the repeated stratified k-fold cross-validation
#' (RSk-fCV) gate, on the Stage 1 band of [cepip_band()].
#'
#' * Curve A is an isotonic (PAVA) fit of v on u for the reference cells
#'   (OriRef) alone; curve B is the same fit for OriRef plus the active
#'   candidates. Each curve is cross-validated with `folds` folds stratified
#'   into deciles of u, repeated `repeats` times, which gives
#'   `folds` x `repeats` test RMSE values per curve.
#' * The gate compares B with A by a one-sided paired Wilcoxon test
#'   (B greater than A) and by dRMSE = median(B - A) / median(A) against the
#'   background variability BV = IQR(A) / median(A).
#' * While p < `alpha` and dRMSE > BV, the active candidate with the largest
#'   |v - PAVA(u)| (PAVA fitted once on OriRef plus all candidates) is removed
#'   and curve B is refitted. The candidates left when the gate closes are the
#'   Recovered cells.
#'
#' Classes: Original + Refined (OriRef; SKCS-core, inside the band),
#' Recovered (non-SKCS, inside the band, passed the gate), Original only
#' (SKCS-core, outside the band) and Excluded (all other cells, including
#' candidates removed by the gate). The Module is OriRef plus Recovered.
#'
#' The fold assignment of repeat r uses the seed `seed` + r x 1000 + stratum,
#' so results do not depend on `workers`.
#'
#' @param skcs,ip `cepip_skcs` and `cepip_ip` objects of the same cells. Not
#'   needed when `band` is given.
#' @param band Optional `cepip_band` object. When it is missing, the band is
#'   built with `kappa` and `bw`.
#' @param kappa,bw Band settings used when `band` is missing.
#' @param folds,repeats Number of folds and of repeats of RSk-fCV.
#' @param alpha Significance level of the gate.
#' @param seed Base seed of the fold assignment.
#' @param workers Number of parallel workers for the repeats.
#' @param tests Run the class tests (needs the Suggests packages clinfun,
#'   Hmisc and energy).
#' @param nperm_jt,nperm_hd,R_dcor Permutations of the Jonckheere-Terpstra,
#'   Hoeffding D and distance correlation tests.
#' @param verbose Print progress messages.
#' @return An object of class `cepip_cmrf`.
#' @examples
#' \donttest{
#' d <- cepip_toy
#' s <- cepip_skcs(d$x, d$y, donor = d$donor, id = d$id)
#' ip <- cepip_ip(d$x, d$y, donor = d$donor, id = d$id)
#' b <- cepip_band(s, ip)
#' # a short gate of 5 repeats (the default is 40)
#' m <- cepip_cmrf(band = b, repeats = 5, tests = FALSE)
#' m$counts
#' }
#' @export
cepip_cmrf <- function(skcs = NULL, ip = NULL, band = NULL, kappa = 0.5, bw = 0.1,
                       folds = 5L, repeats = 40L, alpha = 0.05, seed = 42L,
                       workers = 1L, tests = TRUE, nperm_jt = 10000L,
                       nperm_hd = 9999L, R_dcor = 999L, verbose = FALSE) {
    if (is.null(band)) band <- cepip_band(skcs, ip, kappa = kappa, bw = bw, verbose = verbose)
    check_class(band, "cepip_band", "band")
    rng <- save_rng()
    on.exit(restore_rng(rng), add = TRUE)
    cl <- band$cells
    a <- which(cl$is_anchor)
    cand <- which(cl$is_candidate)
    gate <- rskfcv_gate(cl$u[a], cl$v[a], cl$u[cand], cl$v[cand], folds = folds,
                        repeats = repeats, alpha = alpha, seed = seed,
                        workers = workers, verbose = verbose)
    gate$anchor_row <- a
    gate$cand_row <- cand

    pass2 <- rep(FALSE, nrow(cl))
    pass2[a] <- TRUE
    pass2[cand] <- gate$pass
    removed <- rep(FALSE, nrow(cl))
    removed[cand] <- !gate$pass
    pass <- cl$pass_stage1 & pass2
    cls <- ifelse(cl$is_core & pass, "Original + Refined",
           ifelse(cl$is_core, "Original only",
           ifelse(pass, "Recovered", "Excluded")))
    cl$pass_stage2 <- pass2
    cl$removed_by_stage2 <- removed
    cl$cell_class <- cls
    cl$is_recovered <- cls == "Recovered"
    cl$in_module <- cls %in% c("Original + Refined", "Recovered")
    # Recovered cells are labelled Pre-IP or Post-IP on the IP axis.
    rec_side <- ip_labels(cl$x, cl$y, list(has_ip = TRUE, ip_x = band$IP_x, ip_y = band$IP_y),
                          band$ip_axis)
    cl$group_refined <- cl$group_orig
    cl$group_refined[cl$is_recovered] <- paste("Recovered", rec_side[cl$is_recovered])
    cl$ip_side_all <- rec_side
    # SKCS coordinates in copula space (pitcher plot): ranks of S_i,norm and K_i.
    n <- nrow(cl)
    cl$u_SK <- rank(cl$Si_norm, ties.method = "average") / (n + 1)
    cl$v_SK <- rank(cl$Ki, ties.method = "average") / (n + 1)

    counts <- table(factor(cls, levels = c("Original + Refined", "Recovered",
                                           "Original only", "Excluded")))
    msg(verbose, "  CMRF: OriRef %d | Recovered %d | Original only %d | Excluded %d",
        counts[1], counts[2], counts[3], counts[4])
    out <- list(cells = cl, band = band[setdiff(names(band), "cells")], gate = gate,
                counts = counts, tests = NULL,
                gates_SK = c(u = mean(cl$Si_norm <= 0.05), v = mean(cl$Ki <= 0.5)),
                settings = list(folds = folds, repeats = repeats, alpha = alpha,
                                seed = seed))
    if (isTRUE(tests)) out$tests <- cmrf_tests(cl, nperm_jt, nperm_hd, R_dcor, seed)
    structure(out, class = "cepip_cmrf")
}

# One RSk-fCV repeat. Self-contained (base and stats only) so it can be sent to
# PSOCK workers without the package.
rskfcv_repeat <- function(r, u_set, v_set, folds, seed) {
    n <- length(u_set)
    rep_seed <- seed + r * 1000L
    n_strat <- max(1L, min(10L, floor(n / folds)))
    rk <- rank(u_set, ties.method = "first")
    strat <- pmin(ceiling(rk / n * n_strat), n_strat)
    fa <- integer(n)
    for (s in seq_len(n_strat)) {
        idx <- which(strat == s)
        set.seed(rep_seed + s)
        idx <- idx[sample.int(length(idx))]
        fa[idx] <- ((seq_along(idx) - 1L) %% folds) + 1L
    }
    rmse <- vapply(seq_len(folds), function(k) {
        te <- fa == k
        tr <- !te
        if (sum(tr) < 2L || sum(te) < 1L) return(NA_real_)
        o <- order(u_set[tr])
        us <- u_set[tr][o]
        yf <- stats::isoreg(us, v_set[tr][o])$yf
        i <- findInterval(u_set[te], us, rightmost.closed = TRUE)
        i[i < 1L] <- 1L
        i[i > length(us)] <- length(us)
        sqrt(mean((v_set[te] - yf[i])^2))
    }, numeric(1))
    list(rmse = rmse, folds = fa)
}
environment(rskfcv_repeat) <- baseenv()

# Isotonic (PAVA) prediction at u_new: step lookup into the fitted values.
iso_predict <- function(u_train, v_train, u_new) {
    o <- order(u_train)
    us <- u_train[o]
    yf <- stats::isoreg(us, v_train[o])$yf
    i <- findInterval(u_new, us, rightmost.closed = TRUE)
    i[i < 1L] <- 1L
    i[i > length(us)] <- length(us)
    yf[i]
}

run_rskfcv <- function(u, v, folds, repeats, seed, cl = NULL) {
    if (length(u) < folds + 1L)
        return(list(rmse = rep(NA_real_, folds * repeats), folds = NULL))
    res <- if (is.null(cl)) lapply(seq_len(repeats), rskfcv_repeat, u_set = u, v_set = v,
                                   folds = folds, seed = seed)
           else parallel::parLapply(cl, seq_len(repeats), rskfcv_repeat, u_set = u,
                                    v_set = v, folds = folds, seed = seed)
    list(rmse = unlist(lapply(res, `[[`, "rmse"), use.names = FALSE),
         folds = do.call(cbind, lapply(res, `[[`, "folds")))
}

rskfcv_gate <- function(u_a, v_a, u_c, v_c, folds = 5L, repeats = 40L, alpha = 0.05,
                        seed = 42L, workers = 1L, verbose = FALSE) {
    cl <- NULL
    if (workers > 1L) {
        cl <- parallel::makePSOCKcluster(workers)
        on.exit(parallel::stopCluster(cl), add = TRUE)
    }
    A <- run_rskfcv(u_a, v_a, folds, repeats, seed, cl)
    ra <- A$rmse
    ok_a <- ra[!is.na(ra)]
    med_a <- stats::median(ok_a)
    bv <- if (med_a > 0) stats::IQR(ok_a) / med_a else 0
    nc <- length(u_c)
    out <- list(rmse_A = ra, folds_A = A$folds, bv = bv, median_rmse_A = med_a,
                n_candidates = nc, folds = folds, repeats = repeats, alpha = alpha,
                seed = seed)
    if (nc == 0L) {
        return(c(out, list(pass = logical(0), removal_order = integer(0),
                           resid_pava = numeric(0), p = numeric(0), effect = numeric(0),
                           n_active = integer(0), rmse_B_first = ra, rmse_B_last = ra,
                           folds_B_first = A$folds, folds_B_last = A$folds,
                           p_final = NA_real_, effect_final = 0, n_removed = 0L)))
    }
    # Removal ranking from one PAVA fit on OriRef plus all candidates.
    resid_pava <- abs(v_c - iso_predict(c(u_a, u_c), c(v_a, v_c), u_c))
    active <- rep(TRUE, nc)
    removal <- integer(0); p_it <- numeric(0); eff_it <- numeric(0); n_it <- integer(0)
    first <- NULL; B <- NULL
    repeat {
        idx <- which(active)
        B <- run_rskfcv(c(u_a, u_c[idx]), c(v_a, v_c[idx]), folds, repeats, seed, cl)
        if (is.null(first)) first <- B
        p <- tryCatch(suppressWarnings(stats::wilcox.test(B$rmse, ra, paired = TRUE,
                                                          alternative = "greater",
                                                          exact = FALSE))$p.value,
                      error = function(e) 1)
        ok <- !is.na(B$rmse) & !is.na(ra)
        eff <- if (any(ok)) stats::median(B$rmse[ok] - ra[ok]) / med_a else NA_real_
        p_it <- c(p_it, p); eff_it <- c(eff_it, eff); n_it <- c(n_it, length(idx))
        msg(verbose, "  gate iteration %d | active %d | p = %.3e | dRMSE = %.5f | BV = %.5f",
            length(p_it), length(idx), p, eff, bv)
        if (!(isTRUE(p < alpha) && isTRUE(eff > bv))) break
        if (length(idx) <= 1L) {
            active[idx[1L]] <- FALSE
            removal <- c(removal, idx[1L])
            break
        }
        worst <- idx[which.max(resid_pava[idx])]
        active[worst] <- FALSE
        removal <- c(removal, worst)
    }
    c(out, list(pass = active, removal_order = removal, resid_pava = resid_pava,
                p = p_it, effect = eff_it, n_active = n_it,
                rmse_B_first = first$rmse, rmse_B_last = B$rmse,
                folds_B_first = first$folds, folds_B_last = B$folds,
                p_final = utils::tail(p_it, 1L), effect_final = utils::tail(eff_it, 1L),
                n_removed = length(removal)))
}

# Tests of the final classes: Jonckheere-Terpstra trend of |v - g(u)| across
# Excluded < Original only < Recovered < Original + Refined (decreasing), and
# Hoeffding D and distance correlation of (u, v) within the Recovered cells.
cmrf_tests <- function(cl, nperm_jt, nperm_hd, R_dcor, seed) {
    lev <- c("Excluded", "Original only", "Recovered", "Original + Refined")
    out <- data.frame(test = character(0), n = integer(0), statistic = numeric(0),
                      p_value = numeric(0), nperm = integer(0), stringsAsFactors = FALSE)
    if (requireNamespace("clinfun", quietly = TRUE)) {
        g <- factor(cl$cell_class, levels = lev, ordered = TRUE)
        ok <- !is.na(g) & !is.na(cl$resid_abs)
        set.seed(seed)
        jt <- tryCatch(clinfun::jonckheere.test(cl$resid_abs[ok], as.integer(g[ok]),
                                                alternative = "decreasing", nperm = nperm_jt),
                       error = function(e) NULL)
        if (!is.null(jt))
            out[nrow(out) + 1L, ] <- list("Jonckheere-Terpstra |v - g(u)| by class", sum(ok),
                                          as.numeric(jt$statistic), jt$p.value, nperm_jt)
    } else message("clinfun is not installed: the Jonckheere-Terpstra test is skipped.")
    rec <- cl[cl$cell_class == "Recovered", ]
    if (nrow(rec) > 4L) {
        if (requireNamespace("Hmisc", quietly = TRUE)) {
            d0 <- Hmisc::hoeffd(rec$u, rec$v)$D[1, 2]
            set.seed(seed)
            dp <- vapply(seq_len(nperm_hd), function(i)
                Hmisc::hoeffd(sample(rec$u), rec$v)$D[1, 2], numeric(1))
            out[nrow(out) + 1L, ] <- list("Hoeffding D (u, v) in Recovered", nrow(rec), d0,
                                          (sum(dp >= d0) + 1) / (nperm_hd + 1), nperm_hd)
        } else message("Hmisc is not installed: the Hoeffding D test is skipped.")
        if (requireNamespace("energy", quietly = TRUE)) {
            set.seed(seed)
            dc <- energy::dcor.test(rec$u, rec$v, R = R_dcor)
            out[nrow(out) + 1L, ] <- list("Distance correlation (u, v) in Recovered",
                                          nrow(rec), unname(dc$statistic)[[1]],
                                          unname(dc$p.value)[[1]], R_dcor)
        } else message("energy is not installed: the distance correlation test is skipped.")
    }
    out
}
