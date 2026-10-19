#' Spearman-Kendall cell-level scores (SKCS)
#'
#' Splits the global Spearman and Kendall correlations between x and y into
#' per-cell contributions.
#'
#' * S_i = (r_x,i - mean r_x)(r_y,i - mean r_y) / (n - 1), from average ranks;
#'   S_i,norm = S_i / ((n - 1) / 4).
#' * K_i = (C_i - D_i) / sqrt((n - 1 - t_x,i)(n - 1 - t_y,i)), where C_i and
#'   D_i are the concordant and discordant pairs of cell i and t_x,i and t_y,i
#'   its ties in x and in y (tau-b form).
#'
#' SKCS-core cells pass both gates, S_i,norm > `s_min` and K_i > `k_min`. The
#' CEP-IP GAM is then refitted on the SKCS-core cells (with the same settings
#' as [cepip_fit()]) and the TREP cells are defined within the core.
#'
#' The gates S_i,norm > 0.05 and K_i > 0.5 were introduced in the SEA-AD
#' analysis; the published CEP-IP did not use cell-level scores.
#'
#' @inheritParams cepip_fit
#' @param s_min,k_min SKCS-core gates.
#' @param fit Refit the GAM on the SKCS-core cells and define TREP cells there.
#' @param ... Passed to [cepip_fit()] for the SKCS-core fit (for example `k`).
#' @return An object of class `cepip_skcs`: a per-cell table (`cells`) with
#'   S_i, S_i,norm, K_i and the SKCS-core flag, the global rho and tau, and,
#'   when `fit = TRUE`, the SKCS-core fit and its EP and TREP table.
#' @examples
#' d <- cepip_toy
#' s <- cepip_skcs(d$x, d$y, donor = d$donor, id = d$id)
#' s$n_core
#' @export
cepip_skcs <- function(x, y, donor = NULL, id = NULL, s_min = 0.05, k_min = 0.5,
                       fit = TRUE, verbose = FALSE, ...) {
    n <- length(x)
    if (n < 3L) stop("At least 3 cells are needed.", call. = FALSE)
    x <- as.numeric(x); y <- as.numeric(y)
    sc <- skcs_scores(x, y)
    is_core <- (sc$Si_norm > s_min) & (sc$Ki > k_min)
    cells <- data.frame(id = if (is.null(id)) seq_len(n) else id, x = x, y = y,
                        donor = if (is.null(donor)) NA_character_ else as.character(donor),
                        Si = sc$Si, Si_norm = sc$Si_norm, Ki = sc$Ki,
                        is_core = is_core, stringsAsFactors = FALSE)
    msg(verbose, "  SKCS-core: %d of %d cells", sum(is_core), n)
    out <- list(cells = cells, n_core = sum(is_core),
                rho = stats::cor(x, y, method = "spearman"),
                tau = stats::cor(x, y, method = "kendall"),
                s_min = s_min, k_min = k_min, fit = NULL, ep = NULL)
    if (isTRUE(fit) && sum(is_core) >= 4L) {
        core <- which(is_core)
        f <- cepip_fit(x[core], y[core],
                       donor = if (is.null(donor)) NULL else donor[core],
                       id = cells$id[core], verbose = verbose, ...)
        out$fit <- f
        out$ep <- cepip_ep(f)
    }
    structure(out, class = "cepip_skcs")
}

# Per-cell Spearman and Kendall contributions. Kendall is computed in blocks
# of rows to limit memory; the result does not depend on the block size.
skcs_scores <- function(x, y, block = 500L) {
    n <- length(x)
    rx <- rank(x, ties.method = "average")
    ry <- rank(y, ties.method = "average")
    Si <- (rx - mean(rx)) * (ry - mean(ry)) / (n - 1)
    Si_norm <- Si / ((n - 1) / 4)
    cd <- numeric(n); tx <- numeric(n); ty <- numeric(n)
    for (s in seq.int(1L, n, by = block)) {
        idx <- s:min(s + block - 1L, n)
        sx <- sign(outer(x[idx], x, "-"))
        sy <- sign(outer(y[idx], y, "-"))
        p <- sx * sy
        self <- cbind(seq_along(idx), idx)
        p[self] <- 0
        cd[idx] <- rowSums(p > 0) - rowSums(p < 0)
        zx <- sx == 0; zx[self] <- FALSE
        zy <- sy == 0; zy[self] <- FALSE
        tx[idx] <- rowSums(zx)
        ty[idx] <- rowSums(zy)
    }
    den <- sqrt((n - 1 - tx) * (n - 1 - ty))
    den[den == 0] <- 1
    list(Si = Si, Si_norm = Si_norm, Ki = cd / den)
}
