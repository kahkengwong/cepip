#' Fit the CEP-IP GAM
#'
#' Fits `y ~ x + s(x, bs = "tp", k)` by REML with `select = TRUE` and
#' `gamma = 1.5`, as in the published CEP-IP. When `donor` is given, a donor
#' random intercept `s(donor, bs = "re")` is added (if there are at least
#' `min_donors` donors) and each cell is weighted by 1 / (cells of its donor),
#' so that every donor carries the same total weight.
#'
#' The basis dimension k is chosen by the lowest penalized residual sum of
#' squares (PRSS) over k = 3 to min(`k_max`, unique x - 1), with ties going to
#' the smaller k. The published pipeline ran a second phase of up to five
#' iterations as a sanity check of this search; it cannot change the selected
#' k because the first phase already fits every k in the range and each fit is
#' deterministic, so it is not repeated here.
#'
#' The package does not transform the data. The published CEP-IP analyses
#' used log2(value + 1) of Seurat log-normalized data; the SEA-AD analysis used
#' ln(CP10K + 1) directly.
#'
#' @param x Expression of the gene of interest (GOI), one value per cell.
#' @param y Composite expression of the partner genes (see [cepip_composite()]).
#' @param donor Optional donor label per cell.
#' @param id Optional cell identifiers (for example barcodes).
#' @param k `"prss"` for the PRSS search, or a fixed integer k.
#' @param k_max Largest k in the PRSS search.
#' @param gamma,select Passed to [mgcv::gam()].
#' @param min_donors Minimum number of donors for the donor random intercept.
#' @param workers Number of parallel workers for the k search.
#' @param verbose Print progress messages.
#' @return An object of class `cepip_fit`: the model, the selected k, the PRSS
#'   table, the deviance explained (DE) and the fitting frame.
#' @examples
#' d <- cepip_toy
#' f <- cepip_fit(d$x, d$y, donor = d$donor, id = d$id)
#' f$k
#' f$de
#' @export
cepip_fit <- function(x, y, donor = NULL, id = NULL, k = "prss", k_max = 10L,
                      gamma = 1.5, select = TRUE, min_donors = 3L, workers = 1L,
                      verbose = FALSE) {
    d <- make_frame(x, y, donor, id)
    fit_frame(d, k = k, k_max = k_max, gamma = gamma, select = select,
              min_donors = min_donors, workers = workers, verbose = verbose)
}

gam_formula <- function(k, use_re, resp = "Expression", linear = TRUE, var = "GOI") {
    rhs <- sprintf("s(%s, bs = \"tp\", k = %d)", var, as.integer(k))
    if (linear) rhs <- paste(var, "+", rhs)
    if (use_re) rhs <- paste0(rhs, " + s(donor, bs = \"re\")")
    stats::as.formula(paste(resp, "~", rhs))
}


# One GAM fit and its PRSS. Self-contained (mgcv and base only) so it can run
# on PSOCK workers.
fit_k_worker <- function(k, d, use_re, gamma, select, resp = "Expression", var = "GOI",
                         linear = TRUE, retry_lower = FALSE) {
    fit_at <- function(kk) {
        rhs <- sprintf("s(%s, bs = \"tp\", k = %d)", var, as.integer(kk))
        if (linear) rhs <- paste(var, "+", rhs)
        if (use_re) rhs <- paste0(rhs, " + s(donor, bs = \"re\")")
        mgcv::gam(stats::as.formula(paste(resp, "~", rhs)), data = d, weights = d$w,
                  method = "REML", select = select, gamma = gamma)
    }
    m <- tryCatch(fit_at(k), error = function(e) {
        if (retry_lower) tryCatch(fit_at(max(3L, k - 1L)), error = function(e2) NULL) else NULL
    })
    if (is.null(m)) return(NULL)
    labs <- vapply(m$smooth, function(s) s$label, character(1))
    i <- which(labs == paste0("s(", var, ")"))[1]
    sm <- m$smooth[[i]]
    lambda <- unname(m$sp[sm$first.sp])
    b <- stats::coef(m)[sm$first.para:sm$last.para]
    rss <- sum(stats::residuals(m)^2)
    pen <- as.numeric(t(b) %*% sm$S[[1]] %*% b)
    list(model = m, k = k, PRSS = rss + lambda * pen, RSS = rss, lambda = lambda)
}
environment(fit_k_worker) <- baseenv()

# Fits the GAM on a prepared frame (GOI, Expression, optional donor, w).
fit_frame <- function(d, k = "prss", k_max = 10L, gamma = 1.5, select = TRUE,
                      min_donors = 3L, workers = 1L, verbose = FALSE) {
    unique_x <- length(unique(d$GOI))
    if (unique_x < 4L) stop("x has fewer than 4 unique values.", call. = FALSE)
    use_re <- use_donor_re(d, min_donors)
    k_hi <- min(as.integer(k_max), unique_x - 1L)
    ks <- if (identical(k, "prss")) seq.int(3L, max(3L, k_hi)) else as.integer(k)
    dd <- d[, intersect(c("GOI", "Expression", "donor", "w"), names(d))]
    res <- if (workers > 1L && length(ks) > 1L) {
        cl <- parallel::makePSOCKcluster(min(as.integer(workers), length(ks)))
        on.exit(parallel::stopCluster(cl), add = TRUE)
        parallel::parLapply(cl, ks, fit_k_worker, d = dd, use_re = use_re, gamma = gamma,
                            select = select)
    } else lapply(ks, fit_k_worker, d = dd, use_re = use_re, gamma = gamma, select = select)

    nd <- tryCatch(null_deviance(d, min_donors), error = function(e) NA_real_)
    de_of <- function(m) if (is.finite(nd) && nd > 0) (nd - stats::deviance(m)) / nd else NA_real_
    rows <- lapply(seq_along(ks), function(j) {
        r <- res[[j]]
        if (is.null(r)) return(data.frame(k = ks[j], PRSS = NA_real_, RSS = NA_real_,
                                          lambda = NA_real_, edf = NA_real_, DE = NA_real_))
        data.frame(k = ks[j], PRSS = r$PRSS, RSS = r$RSS, lambda = r$lambda,
                   edf = smooth_stat(r$model, "s(GOI)", "edf"), DE = de_of(r$model))
    })
    prss_table <- do.call(rbind, rows)
    if (all(is.na(prss_table$PRSS))) stop("All GAM fits failed.", call. = FALSE)
    # Lowest PRSS; which.min() returns the first minimum, so ties go to the smaller k.
    j <- which.min(prss_table$PRSS)
    best <- res[[j]]$model
    best_k <- ks[j]
    prss_table$selected <- seq_len(nrow(prss_table)) == j
    for (jj in seq_len(nrow(prss_table)))
        msg(verbose, "  k = %d | PRSS = %.6f", prss_table$k[jj], prss_table$PRSS[jj])

    structure(list(
        model = best, k = best_k, prss = prss_table,
        de = de_of(best), de_flat = as.numeric(summary(best)$dev.expl),
        donor_re = use_re, n_donors = n_donor_levels(d$donor),
        n_cells = nrow(d), data = d,
        settings = list(k = k, k_max = k_max, gamma = gamma, select = select,
                        min_donors = min_donors)),
        class = "cepip_fit")
}

#' Cell explanatory power and TREP cells
#'
#' Computes the explanatory power of each cell,
#' EP_i = 1 - (y_i - f(x_i))^2 / (y_i - mean(y))^2, where f is the
#' population-level curve (donor term excluded) and the null is the plain mean
#' of y. The top round(n x DE) cells by EP are the top-ranked explanatory power
#' (TREP) cells. DE here is the weighted population-level DE,
#' 1 - RSS_w / TSS_w, from the same curve; without donors it equals the usual
#' deviance explained.
#'
#' @param fit A `cepip_fit` object.
#' @return A data frame with one row per cell: id, x, y, donor, fitted value,
#'   EP, EP rank and TREP flag. The attribute `de_pop` holds the DE used for
#'   the TREP count.
#' @examples
#' d <- cepip_toy
#' f <- cepip_fit(d$x, d$y, donor = d$donor, id = d$id)
#' ep <- cepip_ep(f)
#' head(ep)
#' sum(ep$is_TREP)
#' @export
cepip_ep <- function(fit) {
    check_class(fit, "cepip_fit", "fit")
    d <- fit$data
    m <- fit$model
    n <- nrow(d)
    fitted <- as.numeric(stats::predict(m, newdata = d, exclude = DONOR_RE_LABEL,
                                        type = "response"))
    null_fit <- stats::fitted(mgcv::gam(Expression ~ 1, data = d))
    null_sq <- (d$Expression - null_fit)^2
    mod_sq <- (d$Expression - fitted)^2
    ep <- ifelse(null_sq > 0, 1 - mod_sq / null_sq, 0)

    de_pop <- de_population(d$Expression, fitted, d$w)
    if (!is.finite(de_pop)) de_pop <- fit$de_flat
    frac <- max(0, min(1, de_pop))
    is_trep <- rep(FALSE, n)
    # TREP cells are defined only when the model explains some deviance.
    if (is.finite(fit$de) && fit$de > 0) {
        n_trep <- max(1L, min(round(n * frac), n))
        is_trep[order(ep, decreasing = TRUE)[seq_len(n_trep)]] <- TRUE
    }
    out <- data.frame(id = d$id, x = d$GOI, y = d$Expression,
                      donor = if ("donor" %in% names(d)) as.character(d$donor) else NA_character_,
                      fitted = fitted, EP = ep,
                      EP_rank = match(seq_len(n), order(ep, decreasing = TRUE)),
                      is_TREP = is_trep, stringsAsFactors = FALSE)
    attr(out, "de_pop") <- de_pop
    out
}

#' Composite expression of partner genes
#'
#' The composite is the arithmetic mean of the partner genes per cell, on the
#' scale of the input values.
#'
#' @param expr Matrix or data frame of expression values (cells x genes).
#' @param partners Names (or column indices) of the partner genes.
#' @return A numeric vector, one value per cell.
#' @examples
#' d <- cepip_toy
#' y <- cepip_composite(d, paste0("P", 1:8))
#' head(y)
#' @export
cepip_composite <- function(expr, partners) {
    if (is.character(partners)) {
        miss <- setdiff(partners, colnames(expr))
        if (length(miss)) stop("Partner genes not found: ", paste(miss, collapse = ", "),
                               call. = FALSE)
    }
    sub <- expr[, partners, drop = FALSE]
    if (inherits(sub, "Matrix")) return(as.numeric(Matrix::rowMeans(sub)))
    rowMeans(as.matrix(sub))
}
