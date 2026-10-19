#' Validate the CMRF Module by refitting the CEP-IP GAM (Validation1)
#'
#' Refits the CEP-IP GAM (same settings as [cepip_fit()], k by PRSS) on the
#' Module cells (OriRef plus Recovered), recomputes EP and the TREP cells, and
#' labels every Module cell Pre-IP or Post-IP on the IP axis. This gives the
#' Recovered cells their own four-group labels. The agreement between the
#' refitted and the original four-group labels of the OriRef cells is
#' reported.
#'
#' The refitted labels are kept in `group_ref`. Plots of the CMRF classes use
#' the original labels (`group_orig`), because the Original-only cells have no
#' refitted label.
#'
#' @param cmrf A `cepip_cmrf` object.
#' @param min_donors Minimum number of donors for the donor random intercept.
#' @param workers Number of parallel workers for the k search.
#' @param verbose Print progress messages.
#' @return An object of class `cepip_validate`: the refit, the per-cell table of
#'   the Module cells, the label agreement and the cross-tabulation.
#' @examples
#' \donttest{
#' d <- cepip_toy
#' # CMRF with a short gate of 5 repeats (the default is 40)
#' res <- cepip(d$x, d$y, donor = d$donor, id = d$id, cmrf = TRUE, repeats = 5,
#'              tests = FALSE, verbose = FALSE)
#' v <- cepip_validate(res$cmrf)
#' v$de
#' v$agreement
#' }
#' @export
cepip_validate <- function(cmrf, min_donors = 3L, workers = 1L, verbose = FALSE) {
    check_class(cmrf, "cepip_cmrf", "cmrf")
    cl <- cmrf$cells
    mod <- which(cl$in_module)
    if (length(mod) < 10L) stop("Fewer than 10 Module cells.", call. = FALSE)
    has_donor <- !all(is.na(cl$donor))
    d <- make_frame(cl$x[mod], cl$y[mod], if (has_donor) cl$donor[mod] else NULL, cl$id[mod])
    fit <- fit_frame(d, min_donors = min_donors, workers = workers, verbose = verbose)

    fitted <- as.numeric(predict_population(fit$model, d$GOI, d, se = FALSE))
    null_fit <- stats::fitted(mgcv::gam(Expression ~ 1, data = d))
    null_sq <- (d$Expression - null_fit)^2
    ep <- ifelse(null_sq > 0, 1 - (d$Expression - fitted)^2 / null_sq, 0)
    de_pop <- de_population(d$Expression, fitted, d$w)
    frac <- if (is.finite(de_pop)) max(0, min(1, de_pop)) else 0
    n_trep <- max(0L, min(as.integer(round(nrow(d) * frac)), nrow(d)))
    is_trep <- rep(FALSE, nrow(d))
    if (n_trep > 0L) is_trep[order(ep, decreasing = TRUE)[seq_len(n_trep)]] <- TRUE

    ip <- list(has_ip = TRUE, ip_x = cmrf$band$IP_x, ip_y = cmrf$band$IP_y)
    side <- ip_labels(d$GOI, d$Expression, ip, cmrf$band$ip_axis)
    group_ref <- paste(side, ifelse(is_trep, "TREP", "non-TREP"))
    cells <- data.frame(id = d$id, cell_class = cl$cell_class[mod],
                        group_orig = cl$group_orig[mod], group_ref = group_ref,
                        ip_side = side, EP = ep, is_TREP = is_trep,
                        stringsAsFactors = FALSE)
    ori <- cells$cell_class == "Original + Refined"
    agree <- cells$group_ref[ori] == cells$group_orig[ori]
    lev <- c("Pre-IP TREP", "Pre-IP non-TREP", "Post-IP TREP", "Post-IP non-TREP")
    xtab <- table(Original = factor(cells$group_orig[ori], levels = lev),
                  Refitted = factor(cells$group_ref[ori], levels = lev))
    msg(verbose, "  Validation1: k = %d | DE = %.4f | agreement %.1f%%", fit$k, fit$de,
        100 * mean(agree))
    structure(list(fit = fit, cells = cells, de = fit$de, de_pop = de_pop,
                   n_trep = sum(is_trep), agreement = mean(agree), n_oriref = sum(ori),
                   crosstab = xtab),
              class = "cepip_validate")
}
