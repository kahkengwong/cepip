# Internal helpers shared by the fitting steps.
# Fitting frames use the column names GOI (x), Expression (y), donor and w,
# so that the smooth labels are s(GOI) and s(donor) in every model.

DONOR_RE_LABEL <- "s(donor)"

# Donor as a factor with unused levels dropped.
as_donor_factor <- function(x) droplevels(factor(as.character(x)))

n_donor_levels <- function(donor) {
    if (is.null(donor)) return(0L)
    nlevels(as_donor_factor(donor))
}

# w = 1 / (number of cells of this donor in this frame), so each donor sums to 1.
donor_weights <- function(donor) {
    f <- as_donor_factor(donor)
    as.numeric(1 / as.numeric(table(f)[f]))
}

# Builds a fitting frame. Without donors, weights are 1 and no random effect is fitted.
make_frame <- function(x, y, donor = NULL, id = NULL, resp = "Expression") {
    x <- as.numeric(x)
    y <- as.numeric(y)
    if (length(x) != length(y))
        stop("x and y must have the same length.", call. = FALSE)
    d <- data.frame(GOI = x)
    d[[resp]] <- y
    if (!is.null(donor)) {
        if (length(donor) != length(x))
            stop("donor must have the same length as x.", call. = FALSE)
        if (anyNA(donor)) stop("donor has missing values.", call. = FALSE)
        d$donor <- as_donor_factor(donor)
        d$w <- donor_weights(d$donor)
    } else {
        d$w <- rep(1, length(x))
    }
    d$id <- if (is.null(id)) seq_along(x) else id
    d
}

# Recomputes donor levels and weights for a subset of a frame.
reframe <- function(d) {
    if ("donor" %in% names(d)) {
        d$donor <- as_donor_factor(d$donor)
        d$w <- donor_weights(d$donor)
    } else {
        d$w <- rep(1, nrow(d))
    }
    d
}

use_donor_re <- function(d, min_donors) {
    "donor" %in% names(d) && n_donor_levels(d$donor) >= min_donors
}

model_has_re <- function(model) {
    !is.null(model$smooth) &&
        any(vapply(model$smooth, function(s) s$label, character(1)) == DONOR_RE_LABEL)
}

# Index of the s(GOI) smooth (the first non-donor smooth otherwise).
goi_smooth_idx <- function(model, label = "s(GOI)") {
    if (is.null(model$smooth) || length(model$smooth) == 0) return(NA_integer_)
    labs <- vapply(model$smooth, function(s) s$label, character(1))
    i <- which(labs == label)
    if (length(i) == 1L) return(i)
    i <- which(labs != DONOR_RE_LABEL)
    if (length(i) >= 1L) return(i[1])
    NA_integer_
}

smooth_lambda <- function(model, i) {
    if (is.na(i) || is.null(model$sp) || length(model$sp) == 0) return(NA_real_)
    sm <- model$smooth[[i]]
    if (!is.null(sm$first.sp) && sm$first.sp >= 1L && sm$first.sp <= length(model$sp))
        return(unname(model$sp[sm$first.sp]))
    hit <- which(startsWith(names(model$sp), sm$label))
    if (length(hit) > 0) return(unname(model$sp[hit[1]]))
    unname(model$sp[1])
}

# PRSS = RSS + lambda * b'Sb for the s(GOI) smooth. RSS is the sum of squared
# deviance residuals (weighted). Only the wiggliness penalty S[[1]] is used, as
# in the published CEP-IP code; the null-space penalty added by select = TRUE
# does not enter PRSS.
calc_prss <- function(model, label = "s(GOI)") {
    rss <- sum(stats::residuals(model)^2)
    i <- goi_smooth_idx(model, label)
    if (is.na(i)) return(list(PRSS = rss, RSS = rss, lambda = NA_real_, penalty = NA_real_))
    sm <- model$smooth[[i]]
    lambda <- smooth_lambda(model, i)
    beta <- stats::coef(model)[sm$first.para:sm$last.para]
    pen <- as.numeric(t(beta) %*% sm$S[[1]] %*% beta)
    list(PRSS = rss + lambda * pen, RSS = rss, lambda = lambda, penalty = pen)
}

smooth_stat <- function(model, label, column) {
    st <- summary(model)$s.table
    if (is.null(st) || nrow(st) == 0) return(NA_real_)
    i <- match(label, rownames(st))
    if (is.na(i)) i <- 1L
    if (!(column %in% colnames(st))) return(NA_real_)
    unname(st[i, column])
}

# Population-level prediction: the donor term is supplied at an observed level
# and then excluded, so the curve is the fixed part of the model.
predict_population <- function(model, x, fit_data, se = TRUE, var = "GOI") {
    nd <- data.frame(x)
    names(nd) <- var
    if (model_has_re(model)) {
        nd$donor <- as_donor_factor(fit_data$donor)[1]
        return(stats::predict(model, newdata = nd, se.fit = se, exclude = DONOR_RE_LABEL))
    }
    stats::predict(model, newdata = nd, se.fit = se)
}

# Null deviance of the donor-aware null model with the same weights.
null_deviance <- function(d, min_donors, resp = "Expression") {
    f <- if (use_donor_re(d, min_donors))
        stats::as.formula(paste(resp, "~ 1 + s(donor, bs = \"re\")"))
    else stats::as.formula(paste(resp, "~ 1"))
    if (use_donor_re(d, min_donors))
        stats::deviance(mgcv::gam(f, data = d, weights = w, method = "REML"))
    else stats::deviance(mgcv::gam(f, data = d, weights = w))
}

# Deviance explained against the donor-aware null. Without donors this equals
# summary(model)$dev.expl.
de_donor <- function(model, d, min_donors, resp = "Expression") {
    nd <- tryCatch(null_deviance(d, min_donors, resp), error = function(e) NA_real_)
    if (!is.finite(nd) || nd <= 0) return(NA_real_)
    as.numeric((nd - stats::deviance(model)) / nd)
}

# Weighted population-level DE, 1 - RSS_w / TSS_w, from the donor-free prediction.
de_population <- function(y, fitted, w) {
    ok <- is.finite(y) & is.finite(fitted) & is.finite(w) & w > 0
    if (sum(ok) < 2L) return(NA_real_)
    ybar <- sum(w[ok] * y[ok]) / sum(w[ok])
    den <- sum(w[ok] * (y[ok] - ybar)^2)
    if (!is.finite(den) || den <= 0) return(NA_real_)
    as.numeric(1 - sum(w[ok] * (y[ok] - fitted[ok])^2) / den)
}

# Weighted quantile by linear interpolation of the cumulative weights.
# With equal weights it is the type-7 quantile.
weighted_quantile <- function(x, w, prob) {
    ok <- is.finite(x) & is.finite(w) & w > 0
    x <- x[ok]; w <- w[ok]
    if (!length(x)) return(NA_real_)
    if (length(unique(w)) == 1L)
        return(unname(stats::quantile(x, prob, type = 7, names = FALSE)))
    o <- order(x)
    x <- x[o]; w <- w[o]
    cw <- cumsum(w) / sum(w)
    stats::approx(cw, x, xout = prob, method = "linear", rule = 2, ties = "ordered")$y
}

msg <- function(verbose, ...) if (isTRUE(verbose)) message(sprintf(...))

# The gate and the screen subsample set their own seeds. These two helpers keep
# the random number state of the user's session unchanged.
save_rng <- function() {
    g <- globalenv()
    if (exists(".Random.seed", envir = g, inherits = FALSE))
        get(".Random.seed", envir = g, inherits = FALSE) else NULL
}

restore_rng <- function(state) {
    g <- globalenv()
    if (is.null(state)) {
        if (exists(".Random.seed", envir = g, inherits = FALSE))
            rm(".Random.seed", envir = g)
    } else assign(".Random.seed", state, envir = g)
}

check_class <- function(obj, cls, arg) {
    if (!inherits(obj, cls))
        stop(sprintf("%s must be a %s object.", arg, cls), call. = FALSE)
}
