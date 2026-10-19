#' Run CEP-IP, optionally with SKCS and CMRF
#'
#' Runs the whole method on one population of cells.
#'
#' * Always: the CEP-IP GAM with the PRSS k search ([cepip_fit()]), the EP and
#'   TREP cells ([cepip_ep()]) and the inflection point with its IPRS
#'   ([cepip_ip()]). With `skcs = FALSE` and `cmrf = FALSE` this is the
#'   published CEP-IP, except that cells are labelled Pre-IP or Post-IP on the
#'   response axis (`ip_axis = "y"`); set `ip_axis = "x"` for the published
#'   split on the GOI axis.
#' * `skcs = TRUE`: the Spearman-Kendall cell-level scores, the SKCS-core cells
#'   and the SKCS-core fit ([cepip_skcs()]).
#' * `cmrf = TRUE`: also SKCS, the CMRF band ([cepip_band()]), the RSk-fCV gate
#'   ([cepip_cmrf()]) and, by default, Validation1 ([cepip_validate()]). CMRF
#'   always uses the CEP-IP inflection point.
#'
#' When no IP is found, the run stops after the IP step and the object reports
#' why; no error is raised.
#'
#' @inheritParams cepip_fit
#' @inheritParams cepip_ip
#' @param skcs Add the SKCS step.
#' @param cmrf Add the CMRF steps (implies `skcs = TRUE`).
#' @param validate Run Validation1 after CMRF.
#' @param band Optional `cepip_band` object to use instead of building the band.
#' @param kappa,bw Band settings (see [cepip_band()]).
#' @param y_spine Optional composite values for the spine fit (see [cepip_band()]).
#' @param workers Number of parallel workers (GAM k searches and the RSk-fCV gate).
#' @param repeats Number of RSk-fCV repeats of the gate (see [cepip_cmrf()]).
#' @param tests Run the CMRF class tests.
#' @param verbose Print progress messages.
#' @return An object of class `cepip`.
#' @examples
#' d <- cepip_toy
#' res <- cepip(d$x, d$y, donor = d$donor, id = d$id, ip_k = 6)
#' res
#' \donttest{
#' # With SKCS, CMRF and Validation1 (a short gate of 5 repeats)
#' res2 <- cepip(d$x, d$y, donor = d$donor, id = d$id, cmrf = TRUE, repeats = 5,
#'               tests = FALSE)
#' res2
#' }
#' @export
cepip <- function(x, y, donor = NULL, id = NULL, skcs = FALSE, cmrf = FALSE,
                  validate = cmrf, ip_k = "prss", ip_axis = c("y", "x"), band = NULL,
                  kappa = 0.5, bw = 0.1, y_spine = NULL, workers = 1L, repeats = 40L,
                  tests = TRUE, verbose = TRUE) {
    ip_axis <- match.arg(ip_axis)
    if (cmrf) skcs <- TRUE
    if (is.null(id)) id <- seq_along(x)
    t0 <- proc.time()[["elapsed"]]
    msg(verbose, "CEP-IP: %d cells%s", length(x),
        if (is.null(donor)) "" else sprintf(", %d donors", n_donor_levels(donor)))
    fit <- cepip_fit(x, y, donor = donor, id = id, workers = workers)
    msg(verbose, "  GAM: k = %d | DE = %.4f", fit$k, fit$de)
    ep <- cepip_ep(fit)
    ip_fit <- if (identical(ip_k, "prss")) fit else NULL
    if (is.null(ip_fit))
        ip_fit <- cepip_fit(x, y, donor = donor, id = id, k = ip_k, workers = workers)
    ip <- cepip_ip(x, y, donor = donor, id = id, ip_k = ip_k, ip_axis = ip_axis,
                   fit = ip_fit)
    out <- list(fit = fit, ep = ep, ip = ip, skcs = NULL, band = NULL, cmrf = NULL,
                validate = NULL, status = ip$status,
                settings = list(skcs = skcs, cmrf = cmrf, ip_k = ip_k, ip_axis = ip_axis,
                                kappa = kappa, bw = bw, repeats = repeats))
    if (!isTRUE(ip$has_ip)) {
        msg(verbose, "  No IP: %s. The run stops here.", ip$status)
        return(structure(out, class = "cepip"))
    }
    msg(verbose, "  IP: x = %.4f [%.4f, %.4f] | y = %.4f | IPRS = %.3f (%s)", ip$ip_x,
        ip$ip_ci[1], ip$ip_ci[2], ip$ip_y, ip$iprs$IPRS, ip$iprs$tier)
    if (skcs) {
        out$skcs <- cepip_skcs(x, y, donor = donor, id = id, workers = workers)
        msg(verbose, "  SKCS-core: %d cells | TREP %d", out$skcs$n_core,
            sum(out$skcs$ep$is_TREP))
    }
    if (cmrf) {
        out$band <- if (is.null(band))
            cepip_band(out$skcs, ip, kappa = kappa, bw = bw, y_spine = y_spine,
                       workers = workers)
        else band
        msg(verbose, "  Band: delta_base = %.5f | candidates %d", out$band$delta_base,
            sum(out$band$cells$is_candidate))
        out$cmrf <- cepip_cmrf(band = out$band, repeats = repeats, workers = workers,
                               tests = tests)
        msg(verbose, "  CMRF: %s", paste(names(out$cmrf$counts), out$cmrf$counts,
                                         sep = " ", collapse = " | "))
        if (validate) {
            out$validate <- cepip_validate(out$cmrf, workers = workers)
            msg(verbose, "  Validation1: DE = %.4f | agreement %.1f%%", out$validate$de,
                100 * out$validate$agreement)
        }
    }
    msg(verbose, "Done in %.0f s.", proc.time()[["elapsed"]] - t0)
    structure(out, class = "cepip")
}

#' @export
print.cepip <- function(x, ...) {
    cat("CEP-IP analysis of", x$fit$n_cells, "cells")
    if (x$fit$n_donors > 0) cat(" from", x$fit$n_donors, "donors")
    cat("\n")
    cat(sprintf("  GAM: k = %d, DE = %.4f, TREP cells = %d\n", x$fit$k, x$fit$de,
                sum(x$ep$is_TREP)))
    if (isTRUE(x$ip$has_ip)) {
        cat(sprintf("  IP: x = %.4f (95%% CI %.4f to %.4f), y = %.4f, IPRS = %.3f (%s)\n",
                    x$ip$ip_x, x$ip$ip_ci[1], x$ip$ip_ci[2], x$ip$ip_y, x$ip$iprs$IPRS,
                    x$ip$iprs$tier))
    } else cat("  No IP:", x$ip$status, "\n")
    if (!is.null(x$skcs)) cat(sprintf("  SKCS-core cells: %d\n", x$skcs$n_core))
    if (!is.null(x$cmrf)) {
        cat("  CMRF classes:\n")
        for (k in names(x$cmrf$counts)) cat(sprintf("    %-20s %d\n", k, x$cmrf$counts[[k]]))
        cat(sprintf("  Module: %d cells\n", sum(x$cmrf$cells$in_module)))
    }
    if (!is.null(x$validate))
        cat(sprintf("  Validation1: DE = %.4f, label agreement %.1f%%\n", x$validate$de,
                    100 * x$validate$agreement))
    invisible(x)
}

#' @export
summary.cepip <- function(object, ...) {
    s <- data.frame(item = character(0), value = character(0), stringsAsFactors = FALSE)
    add <- function(i, v) s[nrow(s) + 1L, ] <<- list(i, as.character(v))
    add("cells", object$fit$n_cells)
    add("donors", object$fit$n_donors)
    add("GAM k", object$fit$k)
    add("DE", signif(object$fit$de, 6))
    add("TREP cells", sum(object$ep$is_TREP))
    add("IP found", object$ip$has_ip)
    if (isTRUE(object$ip$has_ip)) {
        add("IPx", signif(object$ip$ip_x, 8))
        add("IPx CI", paste(signif(object$ip$ip_ci, 8), collapse = " to "))
        add("IPy", signif(object$ip$ip_y, 8))
        add("IPRS", signif(object$ip$iprs$IPRS, 6))
        add("IPRS tier", object$ip$iprs$tier)
    }
    if (!is.null(object$skcs)) add("SKCS-core cells", object$skcs$n_core)
    if (!is.null(object$band)) add("delta_base", signif(object$band$delta_base, 8))
    if (!is.null(object$cmrf)) {
        for (k in names(object$cmrf$counts)) add(k, object$cmrf$counts[[k]])
        add("Module", sum(object$cmrf$cells$in_module))
        add("gate removals", object$cmrf$gate$n_removed)
    }
    if (!is.null(object$validate)) add("Validation1 agreement", round(object$validate$agreement, 4))
    s
}

#' @export
print.cepip_fit <- function(x, ...) {
    cat(sprintf("CEP-IP GAM: %d cells, k = %d, DE = %.4f%s\n", x$n_cells, x$k, x$de,
                if (x$donor_re) sprintf(", donor random intercept (%d donors)", x$n_donors) else ""))
    invisible(x)
}

#' @export
print.cepip_ip <- function(x, ...) {
    if (isTRUE(x$has_ip))
        cat(sprintf("IP: x = %.4f (95%% CI %.4f to %.4f), y = %.4f, IPRS = %.3f (%s)\n",
                    x$ip_x, x$ip_ci[1], x$ip_ci[2], x$ip_y, x$iprs$IPRS, x$iprs$tier))
    else cat("No IP:", x$status, "\n")
    invisible(x)
}

#' @export
print.cepip_cmrf <- function(x, ...) {
    cat("CMRF classes:\n")
    for (k in names(x$counts)) cat(sprintf("  %-20s %d\n", k, x$counts[[k]]))
    cat(sprintf("Gate: %d of %d candidates removed in %d iterations; final p = %.3g, dRMSE = %.5f, BV = %.5f\n",
                x$gate$n_removed, x$gate$n_candidates, length(x$gate$p), x$gate$p_final,
                x$gate$effect_final, x$gate$bv))
    invisible(x)
}
