#' Inflection point (IP) and IP reliability score (IPRS)
#'
#' Finds the CEP-IP inflection point.
#'
#' 1. The CEP-IP GAM f is fitted on all cells (`ip_k`), and the top
#'    round(n x DE) cells by explanatory power are kept (see [cepip_ep()]).
#' 2. The residuals y - f(x) of these cells are modelled by a residual GAM
#'    R(x) = s(x, bs = "tp", k = min(`resid_k`, unique x - 1)), fitted by REML
#'    (with the donor random intercept and inverse-donor weights when `donor`
#'    is given; the donor term is excluded when predicting).
#' 3. R(x) is evaluated on a grid of `grid` points. IPx is the first
#'    negative-to-positive zero crossing, located by linear interpolation
#'    between the two grid points.
#' 4. The 95% CI bounds are the roots of R(x) + 1.96 se(x) (lower) and
#'    R(x) - 1.96 se(x) (upper), found by Brent's method ([stats::uniroot()])
#'    from the band sign change nearest the IP. When a band has no sign change,
#'    the grid point where the band is closest to zero on that side is used.
#' 5. IPy = f(IPx), the population-level curve at IPx.
#' 6. IPRS = (C1 + C2 + C3 + C4 + C5) / 5: CI precision, a one-sided test of
#'    post-IP against pre-IP residuals, the slope of R(x) at the IP, the
#'    ribbon width of R(x), and the smoothness of R(x) combined with the share
#'    of post-IP residuals that are positive.
#'
#' Cells are labelled Pre-IP or Post-IP on the axis given by `ip_axis`. The
#' default `"y"` splits on the response, y < IPy. The published CEP-IP analyses
#' split on x (x < IPx), which `ip_axis = "x"` reproduces.
#'
#' @inheritParams cepip_fit
#' @param ip_k k of the GAM f: `"prss"` (default) or a fixed integer. The
#'   published SEA-AD analysis used 10.
#' @param resid_k Largest k of the residual GAM R(x).
#' @param grid Number of grid points for R(x).
#' @param z Multiplier of the standard error for the R(x) bands.
#' @param ip_axis `"y"` or `"x"`, the axis that splits Pre-IP from Post-IP cells.
#' @param min_cells Minimum number of top-EP cells for the residual GAM.
#' @param fit Optional `cepip_fit` of the same cells, to reuse instead of refitting.
#' @return An object of class `cepip_ip`. `has_ip` is `FALSE` when R(x) has no
#'   negative-to-positive crossing; the Pre-IP and Post-IP labels are then `NA`.
#' @examples
#' d <- cepip_toy
#' ip <- cepip_ip(d$x, d$y, donor = d$donor, id = d$id)
#' ip$ip_x
#' ip$iprs
#' @export
cepip_ip <- function(x, y, donor = NULL, id = NULL, ip_k = "prss", resid_k = 20L,
                     grid = 2000L, z = 1.96, ip_axis = c("y", "x"), min_cells = 10L,
                     gamma = 1.5, select = TRUE, min_donors = 3L, fit = NULL,
                     verbose = FALSE) {
    ip_axis <- match.arg(ip_axis)
    if (is.null(fit))
        fit <- cepip_fit(x, y, donor = donor, id = id, k = ip_k, gamma = gamma,
                         select = select, min_donors = min_donors, verbose = verbose)
    check_class(fit, "cepip_fit", "fit")
    d <- fit$data
    m <- fit$model
    ep <- cepip_ep(fit)
    n <- nrow(d)

    # Top-EP cells, sized by the weighted population DE of the same curve.
    de_pop <- attr(ep, "de_pop")
    n_top <- max(1L, min(round(n * max(0, min(1, de_pop))), n))
    is_top <- rep(FALSE, n)
    is_top[order(ep$EP, decreasing = TRUE)[seq_len(n_top)]] <- TRUE

    resid_all <- d$Expression -
        as.numeric(predict_population(m, d$GOI, d, se = FALSE))
    rd <- data.frame(GOI = d$GOI[is_top], residual = resid_all[is_top])
    if ("donor" %in% names(d)) rd$donor <- d$donor[is_top]
    rd <- reframe(rd)

    out <- list(fit = fit, ep = ep, is_top = is_top, de_pop = de_pop,
                residual_data = rd, has_ip = FALSE, ip_x = NA_real_,
                ip_ci = c(NA_real_, NA_real_), ci_method = c(NA_character_, NA_character_),
                ip_y = NA_real_, ip_y_ci = c(NA_real_, NA_real_),
                n_crossings = 0L, test = NULL, iprs = NULL, grid = NULL,
                resid_model = NULL, ip_axis = ip_axis, status = "ok",
                settings = list(ip_k = ip_k, resid_k = resid_k, grid = grid, z = z,
                                min_cells = min_cells, min_donors = min_donors))
    if (nrow(rd) < min_cells || length(unique(rd$GOI)) < 3L) {
        out$status <- "too few top-EP cells"
        out$iprs <- iprs_failed()
        out$labels <- rep(NA_character_, n)
        return(structure(out, class = "cepip_ip"))
    }

    # Residual GAM R(x); a k = 8 refit is the fallback if the first fit fails.
    n_unique <- length(unique(rd$GOI))
    resid_re <- use_donor_re(rd, min_donors)
    rm <- tryCatch(
        mgcv::gam(gam_formula(min(resid_k, n_unique - 1L), resid_re, resp = "residual",
                              linear = FALSE),
                  data = rd, weights = w, method = "REML"),
        error = function(e)
            mgcv::gam(gam_formula(min(8L, n_unique - 1L), resid_re, resp = "residual",
                                  linear = FALSE),
                      data = rd, weights = w, method = "REML"))
    out$resid_model <- rm

    x_seq <- seq(min(rd$GOI), max(rd$GOI), length.out = grid)
    pr <- predict_population(rm, x_seq, rd, se = TRUE)
    r_hat <- as.numeric(pr$fit)
    r_se <- as.numeric(pr$se.fit)
    out$grid <- data.frame(x = x_seq, r_hat = r_hat, se = r_se,
                           r_upper = r_hat + z * r_se, r_lower = r_hat - z * r_se)

    # First negative-to-positive crossing, by linear interpolation.
    sc_all <- which(diff(sign(r_hat)) != 0L)
    sc_np <- sc_all[r_hat[sc_all] < 0 & r_hat[sc_all + 1L] > 0]
    out$n_crossings <- length(sc_np)
    ip <- NA_real_
    if (length(sc_np) > 0L) {
        i <- sc_np[1L]
        ip <- x_seq[i] - r_hat[i] * (x_seq[i + 1L] - x_seq[i]) / (r_hat[i + 1L] - r_hat[i])
    }

    anchor <- if (!is.na(ip)) ip else stats::median(x_seq)
    lo <- tryCatch(ci_brent(rm, x_seq, r_hat, r_se, +z, anchor, rd), error = function(e) NA_real_)
    hi <- tryCatch(ci_brent(rm, x_seq, r_hat, r_se, -z, anchor, rd), error = function(e) NA_real_)
    meth <- c(if (!is.na(lo)) "brent" else NA_character_,
              if (!is.na(hi)) "brent" else NA_character_)
    if (is.na(lo) && !is.na(ip)) {
        fb <- ci_proximity(x_seq, r_hat, r_se, +z, ip, "lower")
        if (!is.na(fb)) { lo <- fb; meth[1] <- "proximity" }
    }
    if (is.na(hi) && !is.na(ip)) {
        fb <- ci_proximity(x_seq, r_hat, r_se, -z, ip, "upper")
        if (!is.na(fb)) { hi <- fb; meth[2] <- "proximity" }
    }
    if (!is.na(lo) && !is.na(hi) && lo > hi) {
        tmp <- lo; lo <- hi; hi <- tmp
        meth <- rev(meth)
    }

    # Post-IP against pre-IP residuals of the top-EP cells, split on x.
    test <- if (!is.na(ip)) ip_test(rd$residual[rd$GOI < ip], rd$residual[rd$GOI >= ip])
            else list(test = NA_character_, p = NA_real_, both_normal = NA)
    iprs <- compute_iprs(ip, lo, hi, test$p, x_seq, r_hat,
                         out$grid$r_upper, out$grid$r_lower, rd)

    ip_y <- fitted_at(m, ip)
    ip_y_ci <- fitted_at(m, c(lo, hi))
    if (all(is.finite(ip_y_ci))) ip_y_ci <- sort(ip_y_ci)

    out$has_ip <- !is.na(ip)
    out$ip_x <- ip
    out$ip_ci <- c(lo, hi)
    out$ci_method <- meth
    out$ip_y <- ip_y
    out$ip_y_ci <- ip_y_ci
    out$test <- test
    out$iprs <- iprs
    if (!out$has_ip) out$status <- "no IP: R(x) has no negative-to-positive crossing"
    out$labels <- ip_labels(d$GOI, d$Expression, out, ip_axis)
    structure(out, class = "cepip_ip")
}

# Pre-IP / Post-IP labels on the chosen axis.
ip_labels <- function(x, y, ip, axis = "y") {
    if (!isTRUE(ip$has_ip)) return(rep(NA_character_, length(x)))
    if (axis == "y") ifelse(y < ip$ip_y, "Pre-IP", "Post-IP")
    else ifelse(x < ip$ip_x, "Pre-IP", "Post-IP")
}

# Population-level fitted response at x (donor term excluded).
fitted_at <- function(model, x) {
    out <- rep(NA_real_, length(x))
    ok <- is.finite(x)
    if (!any(ok)) return(out)
    nd <- data.frame(GOI = x[ok])
    dl <- model$var.summary[["donor"]]
    if (!is.null(dl)) {
        lv <- if (is.factor(dl)) levels(dl) else as.character(dl)
        nd$donor <- factor(lv[1], levels = lv)
    }
    out[ok] <- as.numeric(stats::predict(model, newdata = nd, exclude = DONOR_RE_LABEL,
                                         type = "response"))
    out
}

ci_brent <- function(model, x_seq, r_hat, r_se, mult, near_x, fit_data) {
    band <- r_hat + mult * r_se
    chg <- which(diff(sign(band)) != 0L)
    if (length(chg) == 0L) return(NA_real_)
    if (length(chg) > 1L) chg <- chg[which.min(abs(x_seq[chg] - near_x))]
    br <- c(x_seq[chg], x_seq[chg + 1L])
    tryCatch(
        stats::uniroot(function(v) {
            p <- predict_population(model, v, fit_data, se = TRUE)
            as.numeric(p$fit) + mult * as.numeric(p$se.fit)
        }, interval = br, tol = .Machine$double.eps^0.5)$root,
        error = function(e) {
            y1 <- band[chg]; y2 <- band[chg + 1L]
            br[1] + (0 - y1) * (br[2] - br[1]) / (y2 - y1)
        })
}

ci_proximity <- function(x_seq, r_hat, r_se, mult, ip, side) {
    band <- r_hat + mult * r_se
    idx <- if (side == "lower") which(x_seq <= ip) else which(x_seq >= ip)
    if (length(idx) < 3L) return(NA_real_)
    b <- band[idx]
    if (any(diff(sign(b)) != 0L)) return(NA_real_)
    x_seq[idx[which.min(abs(b))]]
}

# One-sided test, post-IP greater than pre-IP: Welch t if both groups pass
# Shapiro-Wilk (p > 0.05), else Mann-Whitney U. Shapiro-Wilk is skipped above
# 5,000 cells, so the Mann-Whitney U test is used.
ip_test <- function(pre, post, min_pre = 3L, min_post = 5L, alpha = 0.05) {
    out <- list(test = NA_character_, p = NA_real_, both_normal = NA)
    if (length(pre) < min_pre || length(post) < min_post) return(out)
    if (length(pre) > 5000L || length(post) > 5000L) {
        both <- FALSE
    } else {
        sp <- suppressWarnings(stats::shapiro.test(pre))$p.value
        so <- suppressWarnings(stats::shapiro.test(post))$p.value
        both <- !is.na(sp) && !is.na(so) && sp > alpha && so > alpha
    }
    out$both_normal <- both
    if (both) {
        out$test <- "Welch t"
        out$p <- suppressWarnings(stats::t.test(post, pre, alternative = "greater",
                                                var.equal = FALSE))$p.value
    } else {
        out$test <- "Mann-Whitney U"
        out$p <- suppressWarnings(stats::wilcox.test(post, pre, alternative = "greater",
                                                     exact = FALSE))$p.value
    }
    out
}

iprs_tier <- function(v) {
    if (v >= 0.80) "Very Strong" else if (v >= 0.60) "Strong" else
        if (v >= 0.40) "Moderate" else if (v >= 0.20) "Weak" else "Failed"
}

iprs_failed <- function() {
    list(C1 = 0, C2 = 0, C3 = 0, C4 = 0, C5 = 0, IPRS = 0, tier = "Failed",
         flag = "No IP - do not proceed")
}

compute_iprs <- function(ip, lo, hi, p, x_seq, r_hat, r_upper, r_lower, rd) {
    rng <- max(x_seq) - min(x_seq)
    have_ip <- !is.na(ip)
    lo_ok <- have_ip && !is.na(lo)
    hi_ok <- have_ip && !is.na(hi)
    ci_width <- if (!is.na(lo) && !is.na(hi)) abs(hi - lo) else NA_real_

    # C1: CI precision relative to the x range.
    C1 <- if (lo_ok && hi_ok && rng > 0) max(0, 1 - ci_width / rng)
          else if ((lo_ok || hi_ok) && rng > 0)
              0.5 * max(0, 1 - (if (lo_ok) abs(ip - lo) else abs(hi - ip)) / rng)
          else 0
    # C2: p of the post-IP against pre-IP test; p = 0.05 gives 0.5.
    C2 <- if (!is.na(p)) tanh(atanh(0.5) / (-log10(0.05)) * (-log10(max(p, 1e-300)))) else 0
    # C3: absolute slope of R(x) at the grid point nearest the IP.
    C3 <- 0
    if (have_ip && length(x_seq) >= 3L) {
        i <- which.min(abs(x_seq - ip)); n <- length(x_seq)
        dv <- if (i > 1L && i < n) (r_hat[i + 1L] - r_hat[i - 1L]) / (x_seq[i + 1L] - x_seq[i - 1L])
              else if (i == 1L) (r_hat[2L] - r_hat[1L]) / (x_seq[2L] - x_seq[1L])
              else (r_hat[n] - r_hat[n - 1L]) / (x_seq[n] - x_seq[n - 1L])
        C3 <- tanh(abs(dv))
    }
    # C4: mean ribbon width of R(x) relative to its amplitude.
    amp <- diff(range(r_hat, na.rm = TRUE))
    C4 <- exp(-mean(r_upper - r_lower, na.rm = TRUE) / (amp + 1e-6))
    # C5: smoothness of R(x) (total variation over amplitude) times the squared
    # share of post-IP residuals that are not negative.
    tv <- sum(abs(diff(r_hat)), na.rm = TRUE)
    neg <- 0
    if (have_ip && nrow(rd) > 0L) {
        post <- rd$residual[rd$GOI >= ip]
        neg <- if (length(post)) sum(post < 0) / length(post) else 0
    }
    C5 <- if (amp < 1e-6) 0 else exp(-max(0, tv / (amp + 1e-6) - 1)) * (1 - neg)^2

    v <- (C1 + C2 + C3 + C4 + C5) / 5
    tier <- iprs_tier(v)
    flag <- if (any(c(C1, C2, C3, C4, C5) == 0)) "Visual review needed [C=0]" else
        if (tier %in% c("Very Strong", "Strong")) "None" else
            if (tier %in% c("Moderate", "Weak")) "Visual review needed" else "No IP - do not proceed"
    if (!have_ip) { tier <- "Failed"; flag <- "No IP - do not proceed" }
    list(C1 = C1, C2 = C2, C3 = C3, C4 = C4, C5 = C5, IPRS = v, tier = tier, flag = flag,
         ci_width = ci_width, x_range = rng, amplitude = amp, total_variation = tv)
}
