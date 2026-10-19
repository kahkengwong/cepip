# Plots of the CEP-IP and CMRF steps. Every function returns a ggplot object
# (one panel) or a gtable from gridExtra::arrangeGrob() (several panels; draw
# it with grid::grid.draw()). No function writes files.

utils::globalVariables(c(
    "x", "y", "u", "v", "g", "grp", "fit", "lo", "hi", "k", "PRSS", "rs", "tau",
    "gene", "is_dfg", "r_hat", "r_lower", "r_upper", "GOI", "residual",
    "Si_norm", "Ki", "u_SK", "v_SK", "f_lo", "f_up", "a_lo", "a_up",
    "x0", "x1", "y0", "y1", "fold", "role", "count", "stratum", "uu", "yy",
    "vv", "ff", "rmse", "curve", "rep_id", "fold_id", "d", "rank_i", "r",
    "removed", "p", "eff", "iter", "cls", "side", "n", "resid_abs", "Original",
    "Refitted", "lab", "donor", "pct_PostIP", "outcome_value", "avg_log2FC",
    "nlp", "dir", "outcome", "beta_std", "ci_low", "ci_high", "delta_R2", "sig",
    "txt", "label", "nm_cls", "ycen", "vj"))

# == Palette (as in the paper) ================================================
COL_FOUR <- c("Pre-IP TREP" = "#E69F00", "Pre-IP non-TREP" = "#1B3A8C",
              "Post-IP TREP" = "#009E73", "Post-IP non-TREP" = "#D94F4F")
COL_GRAY <- "gray65"
COL_BG <- "gray80"
COL_SPINE <- "#00695C"
COL_IP <- "#1A237E"
COL_IPV <- "gray10"
COL_FLAT <- "#0072B2"
COL_ADAPT <- "#D55E00"
COL_ADAPT_BORDER <- "#B84A00"
COL_BLUE_CAND <- "#2171B5"
COL_GAP_CAND <- "#2C2C2C"
COL_ORIONLY <- "#8B6914"
COL_RECOV <- "#7C3AED"
COL_REMOVED_X <- "#404040"
COL_QC <- c("#8B0000", "#4A148C", "#202020")
COL_A <- "#009E73"
COL_B <- "#1A237E"
COL_TEST <- "#D94F4F"
COL_TRAIN <- grDevices::rgb(0, 130, 176, maxColorValue = 255)
COL_TESTF <- grDevices::rgb(238, 207, 206, maxColorValue = 255)
COL_BAR_REM <- "#000000"
COL_BAR_ACC <- "#5B21B6"
COL_STICK <- "#8A9BAD"
COL_CLASS <- c("Original + Refined" = "#009E73", "Recovered" = "#7C3AED",
               "Original only" = "#A0522D", "Excluded" = "#8A9BAD")
COL_SINGLE <- "#5B21B6"
COL_STRATA <- grDevices::rgb(
    c(7, 40, 82, 120, 159, 195, 225, 245, 252, 246),
    c(5, 11, 14, 28, 42, 59, 86, 125, 168, 215),
    c(27, 83, 109, 109, 99, 79, 53, 21, 13, 70), maxColorValue = 255)
CLASS_LEVELS <- c("Original + Refined", "Recovered", "Original only", "Excluded")

# == Small helpers ============================================================
theme_cepip <- function(bs = 10) {
    ggplot2::theme_minimal(base_size = bs) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = bs + 1),
            plot.subtitle = ggplot2::element_text(size = bs - 1.5, colour = "gray30",
                                                  lineheight = 1.1),
            axis.title = ggplot2::element_text(size = bs - 0.5),
            axis.text = ggplot2::element_text(size = bs - 2),
            panel.grid.major = ggplot2::element_line(colour = "gray92"),
            panel.grid.minor = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(colour = "black", fill = NA,
                                                 linewidth = 0.35),
            plot.background = ggplot2::element_rect(fill = "white", colour = NA),
            panel.background = ggplot2::element_rect(fill = "white", colour = NA),
            legend.text = ggplot2::element_text(size = bs - 2),
            legend.title = ggplot2::element_text(size = bs - 1.5, face = "bold"),
            legend.key.size = grid::unit(0.38, "cm"))
}

# Legend inside the panel, for ggplot2 before and after 3.5.0.
legend_inside <- function(pos, just = pos) {
    box <- ggplot2::theme(
        legend.background = ggplot2::element_rect(fill = grDevices::adjustcolor("white", 0.85),
                                                  colour = NA),
        legend.justification = just)
    if (utils::packageVersion("ggplot2") >= "3.5.0")
        box + ggplot2::theme(legend.position = "inside", legend.position.inside = pos)
    else box + ggplot2::theme(legend.position = pos)
}

fmt_p <- function(p) {
    if (length(p) == 0L || is.na(p)) return("NA")
    if (p < 1e-3) formatC(p, format = "e", digits = 2) else sprintf("%.3f", p)
}

fmt_n <- function(n) formatC(n, format = "d", big.mark = ",")

with_counts <- function(levels, values) {
    tab <- table(factor(values, levels = levels))
    stats::setNames(sprintf("%s (n=%s)", levels, fmt_n(as.integer(tab))), levels)
}

empty_panel <- function(title, text) {
    ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = text, size = 3.5,
                          colour = "gray30") +
        ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
        ggplot2::labs(title = title, x = NULL, y = NULL) +
        theme_cepip() +
        ggplot2::theme(axis.text = ggplot2::element_blank(),
                       panel.grid = ggplot2::element_blank())
}

# Combines panels into one gtable with an overall title.
combine_panels <- function(grobs, layout, title, width, height) {
    top <- if (is.null(title)) NULL else
        grid::textGrob(title, x = grid::unit(0.01, "npc"), hjust = 0,
                       gp = grid::gpar(fontsize = 15, fontface = "bold"))
    g <- gridExtra::arrangeGrob(grobs = grobs, layout_matrix = layout, top = top)
    attr(g, "cepip_size") <- c(width = width, height = height)
    class(g) <- c("cepip_figure", class(g))
    g
}

#' Draw a multi-panel cepip figure
#'
#' The multi-panel plots ([plot_band()], [plot_rskfcv()], [plot_validate()])
#' return a gtable with the extra class `cepip_figure`, so that printing it (at
#' the console or in a report) draws it on a new page.
#'
#' @param x A `cepip_figure` object.
#' @param ... Not used.
#' @return `x`, invisibly.
#' @export
print.cepip_figure <- function(x, ...) {
    grid::grid.newpage()
    grid::grid.draw(x)
    invisible(x)
}

# Short y title (v, v_SK) written horizontally.
theme_vtitle <- function() ggplot2::theme(axis.title.y = ggplot2::element_text(angle = 0, vjust = 0.5))

# The same legend settings for colour, fill and shape, so the three merge.
merged_guides <- function(size = 2.2, order = 1) {
    ggplot2::guides(colour = ggplot2::guide_legend(order = order,
                                                   override.aes = list(size = size, alpha = 1)),
                    fill = ggplot2::guide_legend(order = order),
                    shape = ggplot2::guide_legend(order = order))
}

`%|||%` <- function(a, b) if (is.null(a)) b else a

# Band parameters, the per-cell table and the CMRF object from any of the
# objects that hold them.
resolve_band <- function(band, cmrf = NULL) {
    if (inherits(band, "cepip")) {
        if (is.null(cmrf)) cmrf <- band$cmrf
        band <- band$band
    }
    if (inherits(band, "cepip_cmrf")) {
        if (is.null(cmrf)) cmrf <- band
        band <- band$band
    }
    if (is.null(band)) stop("No CMRF band was found.", call. = FALSE)
    if (!is.null(cmrf)) check_class(cmrf, "cepip_cmrf", "cmrf")
    cells <- if (!is.null(cmrf)) cmrf$cells else band$cells
    if (is.null(cells)) stop("The band has no per-cell table.", call. = FALSE)
    list(b = band, cells = cells, cmrf = cmrf)
}

# Spine, flat and adaptive band edges on the spine grid, clipped to [0, 1].
band_grid <- function(b) {
    u <- b$spine_grid$u
    g <- b$spine_grid$g
    K <- exp(-(u - b$IP_u)^2 / (2 * b$bw^2))
    da <- b$delta_base * (1 + b$kappa * K)
    data.frame(u = u, g = g,
               f_lo = pmax(g - b$delta_base, 0), f_up = pmin(g + b$delta_base, 1),
               a_lo = pmax(g - da, 0), a_up = pmin(g + da, 1))
}

layer_flat <- function(gr) list(
    ggplot2::geom_ribbon(data = gr, ggplot2::aes(x = u, ymin = f_lo, ymax = f_up),
                         fill = COL_FLAT, alpha = 0.18, inherit.aes = FALSE),
    ggplot2::geom_line(data = gr, ggplot2::aes(x = u, y = f_up), colour = COL_FLAT,
                       linewidth = 0.45, linetype = "dashed", inherit.aes = FALSE),
    ggplot2::geom_line(data = gr, ggplot2::aes(x = u, y = f_lo), colour = COL_FLAT,
                       linewidth = 0.45, linetype = "dashed", inherit.aes = FALSE))

layer_adapt_fill <- function(gr, fill = COL_ADAPT, alpha = 0.18)
    ggplot2::geom_ribbon(data = gr, ggplot2::aes(x = u, ymin = a_lo, ymax = a_up),
                         fill = fill, alpha = alpha, inherit.aes = FALSE)

layer_adapt_border <- function(gr, lw = 0.6) list(
    ggplot2::geom_line(data = gr, ggplot2::aes(x = u, y = a_up), colour = COL_ADAPT_BORDER,
                       linewidth = lw, inherit.aes = FALSE),
    ggplot2::geom_line(data = gr, ggplot2::aes(x = u, y = a_lo), colour = COL_ADAPT_BORDER,
                       linewidth = lw, inherit.aes = FALSE))

layer_spine <- function(gr) ggplot2::geom_line(data = gr, ggplot2::aes(x = u, y = g),
                                               colour = COL_SPINE, linewidth = 1.1,
                                               inherit.aes = FALSE)

layer_ip <- function(b, ip_v = TRUE) {
    out <- list(ggplot2::geom_vline(xintercept = b$IP_u, colour = COL_IP,
                                    linewidth = 0.55, linetype = "dashed"))
    if (ip_v && is.finite(b$IP_v))
        out <- c(out, list(ggplot2::geom_hline(yintercept = b$IP_v, colour = COL_IPV,
                                               linewidth = 0.45, linetype = "dashed")))
    out
}

# Gray non-SKCS cells and four-group SKCS-core cells.
layer_cells <- function(cl, size = 1.1, alpha = 0.5) {
    gray <- cl[!cl$group_orig %in% names(COL_FOUR), , drop = FALSE]
    col <- cl[cl$group_orig %in% names(COL_FOUR), , drop = FALSE]
    list(ggplot2::geom_point(data = gray, ggplot2::aes(x = u, y = v), colour = COL_GRAY,
                             size = 0.5, alpha = 0.18, inherit.aes = FALSE),
         ggplot2::geom_point(data = col, ggplot2::aes(x = u, y = v, colour = group_orig),
                             size = size, alpha = alpha, inherit.aes = FALSE))
}
utils::globalVariables("group_orig")

scale_four <- function(guide = "none", labels = ggplot2::waiver(), name = NULL)
    ggplot2::scale_colour_manual(values = COL_FOUR, limits = names(COL_FOUR),
                                 labels = labels, name = name, guide = guide,
                                 drop = FALSE)

uv_panel <- function(p, title, subtitle = NULL, xlim = c(0, 1), ylim = c(-0.03, 1.03),
                     xbreaks = seq(0, 1, 0.25)) {
    p + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, expand = TRUE) +
        ggplot2::scale_x_continuous(breaks = xbreaks) +
        ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.25)) +
        ggplot2::labs(title = title, subtitle = subtitle, x = "u", y = "v") +
        theme_cepip() + ggplot2::theme(aspect.ratio = 1) + theme_vtitle()
}

# =============================================================================
# 1. Screen
# =============================================================================

#' Plot the partner-gene screen of one GOI
#'
#' Scatter of Spearman rs (x) against Kendall tau-b (y) for the partner genes
#' of one GOI that passed the Spearman step of [cepip_screen()]. The dual-filter
#' genes (DFGs) are highlighted and labelled; dashed lines mark the rs and tau
#' thresholds.
#'
#' @param screen A `cepip_screen` object.
#' @param goi The GOI to show. By default the top-ranked GOI.
#' @param title Optional plot title.
#' @return A ggplot object.
#' @examples
#' d <- cepip_toy
#' expr <- t(as.matrix(d[, c("x", paste0("P", 1:8))]))
#' colnames(expr) <- d$id
#' sc <- cepip_screen(expr, goi = "x", donor = d$donor)
#' plot_screen(sc)
#' @export
plot_screen <- function(screen, goi = NULL, title = NULL) {
    check_class(screen, "cepip_screen", "screen")
    if (is.null(goi)) {
        if (!NROW(screen$ranking)) stop("The screen found no DFGs; give goi.", call. = FALSE)
        goi <- screen$ranking$GOI[1]
    }
    sp <- screen$spearman[screen$spearman$GOI == goi, , drop = FALSE]
    if (!nrow(sp)) stop("No partner of ", goi, " passed the Spearman step.", call. = FALSE)
    dfg <- screen$partners$partner[screen$partners$GOI == goi]
    d <- data.frame(gene = sp$partner, rs = sp$rs, tau = sp$tau,
                    is_dfg = ifelse(sp$partner %in% dfg, "DFG", "Spearman only"))
    st <- screen$settings
    neg <- isTRUE(st$negative) && any(d$rs < 0)
    vx <- if (neg) c(-st$rs, st$rs) else st$rs
    hy <- if (neg) c(-st$tau, st$tau) else st$tau
    lab <- d[d$is_dfg == "DFG", , drop = FALSE]
    lab <- lab[order(-abs(lab$rs)), , drop = FALSE]
    ggplot2::ggplot(d, ggplot2::aes(x = rs, y = tau)) +
        ggplot2::geom_vline(xintercept = vx, colour = "gray40", linetype = "dashed",
                            linewidth = 0.4) +
        ggplot2::geom_hline(yintercept = hy, colour = "gray40", linetype = "dashed",
                            linewidth = 0.4) +
        ggplot2::geom_point(ggplot2::aes(colour = is_dfg), size = 2, alpha = 0.8) +
        ggplot2::geom_text(data = lab, ggplot2::aes(label = gene), size = 2.8,
                           vjust = -0.7, check_overlap = TRUE, colour = "gray15") +
        ggplot2::scale_colour_manual(values = c("DFG" = "#D94F4F", "Spearman only" = COL_GRAY),
                                     name = NULL) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.12))) +
        ggplot2::labs(title = title %|||% sprintf("Partner genes of %s", goi),
                      subtitle = sprintf("%d partners pass |rs| >= %.2f, %d are DFGs (tau-b threshold %.2f). %s cells.",
                                         nrow(d), st$rs, nrow(lab), st$tau, fmt_n(screen$n_cells)),
                      x = "Spearman rs", y = "Kendall tau-b") +
        theme_cepip() + ggplot2::theme(legend.position = "bottom")
}

# =============================================================================
# 2. PRSS
# =============================================================================

#' Plot the PRSS search over the basis dimension k
#'
#' @param fit A `cepip_fit` object, or a `cepip` or `cepip_ip` object (its fit
#'   is used).
#' @param title Optional plot title.
#' @return A ggplot object.
#' @examples
#' d <- cepip_toy
#' res <- cepip(d$x, d$y, donor = d$donor, id = d$id, verbose = FALSE)
#' plot_prss(res)
#' @export
plot_prss <- function(fit, title = NULL) {
    if (inherits(fit, "cepip") || inherits(fit, "cepip_ip")) fit <- fit$fit
    check_class(fit, "cepip_fit", "fit")
    d <- fit$prss[is.finite(fit$prss$PRSS), , drop = FALSE]
    if (!nrow(d)) stop("The PRSS table is empty.", call. = FALSE)
    sel <- d[d$k == fit$k, , drop = FALSE]
    ggplot2::ggplot(d, ggplot2::aes(x = k, y = PRSS)) +
        ggplot2::geom_line(colour = "gray40", linewidth = 0.6) +
        ggplot2::geom_point(size = 2, colour = "gray25") +
        ggplot2::geom_point(data = sel, size = 4, shape = 21, fill = "#D94F4F",
                            colour = "black", stroke = 0.6) +
        ggplot2::geom_text(data = sel, ggplot2::aes(label = sprintf("k = %d", k)),
                           vjust = -1.2, size = 3.2, colour = "#D94F4F", fontface = "bold") +
        ggplot2::scale_x_continuous(breaks = d$k) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.15))) +
        ggplot2::labs(title = title %|||% "PRSS by basis dimension k",
                      subtitle = sprintf("Selected k = %d (lowest PRSS). DE = %.4f, %s cells.",
                                         fit$k, fit$de, fmt_n(fit$n_cells)),
                      x = "k", y = "PRSS") +
        theme_cepip()
}

# =============================================================================
# 3. GAM fit
# =============================================================================

#' Plot the CEP-IP GAM with the inflection point
#'
#' GOI (x) against the composite (y). Cells are coloured by the four groups
#' (Pre-IP or Post-IP, TREP or non-TREP) when an IP was found; otherwise TREP
#' cells are coloured and non-TREP cells are gray. The population curve
#' (donor term excluded) is drawn with its 95% band, the IPx with its 95% CI
#' as a shaded strip, and IPy as a horizontal line.
#'
#' @param x A `cepip` object (uses `x$fit`, `x$ep` and `x$ip`), a `cepip_ip`
#'   object, or a `cepip_fit` object (no IP lines).
#' @param title Optional plot title.
#' @return A ggplot object.
#' @examples
#' d <- cepip_toy
#' res <- cepip(d$x, d$y, donor = d$donor, id = d$id, verbose = FALSE)
#' plot_fit(res)
#' @export
plot_fit <- function(x, title = NULL) {
    if (inherits(x, "cepip")) { fit <- x$fit; ep <- x$ep; ip <- x$ip }
    else if (inherits(x, "cepip_ip")) { fit <- x$fit; ep <- x$ep; ip <- x }
    else if (inherits(x, "cepip_fit")) { fit <- x; ep <- cepip_ep(x); ip <- NULL }
    else stop("x must be a cepip, cepip_ip or cepip_fit object.", call. = FALSE)
    has_ip <- !is.null(ip) && isTRUE(ip$has_ip)
    d <- data.frame(x = ep$x, y = ep$y, trep = ep$is_TREP)
    if (has_ip) {
        side <- ip_labels(ep$x, ep$y, ip, ip$ip_axis)
        d$grp <- paste(side, ifelse(d$trep, "TREP", "non-TREP"))
        cols <- COL_FOUR
    } else {
        d$grp <- ifelse(d$trep, "TREP", "non-TREP")
        cols <- c("TREP" = COL_SINGLE, "non-TREP" = COL_GRAY)
    }
    d$grp <- factor(d$grp, levels = names(cols))
    d <- d[order(d$trep), ]
    xs <- seq(min(d$x), max(d$x), length.out = 400L)
    pr <- predict_population(fit$model, xs, fit$data, se = TRUE)
    cv <- data.frame(x = xs, fit = as.numeric(pr$fit),
                     lo = as.numeric(pr$fit) - 1.96 * as.numeric(pr$se.fit),
                     hi = as.numeric(pr$fit) + 1.96 * as.numeric(pr$se.fit))
    p <- ggplot2::ggplot()
    sub <- sprintf("%s cells, k = %d, DE = %.4f, TREP cells = %s.", fmt_n(nrow(d)), fit$k,
                   fit$de, fmt_n(sum(d$trep)))
    if (has_ip) {
        if (all(is.finite(ip$ip_ci)))
            p <- p + ggplot2::annotate("rect", xmin = ip$ip_ci[1], xmax = ip$ip_ci[2],
                                       ymin = -Inf, ymax = Inf, fill = COL_IP, alpha = 0.10)
        sub <- paste0(sub, sprintf("\nIPx = %.3f (95%% CI %.3f to %.3f), IPy = %.3f, IPRS = %.3f (%s).",
                                   ip$ip_x, ip$ip_ci[1], ip$ip_ci[2], ip$ip_y,
                                   ip$iprs$IPRS, ip$iprs$tier))
    } else if (!is.null(ip)) sub <- paste0(sub, "\nNo IP: ", ip$status, ".")
    p <- p +
        ggplot2::geom_point(data = d, ggplot2::aes(x = x, y = y, colour = grp),
                            size = 0.9, alpha = 0.55) +
        ggplot2::geom_ribbon(data = cv, ggplot2::aes(x = x, ymin = lo, ymax = hi),
                             fill = "gray20", alpha = 0.2) +
        ggplot2::geom_line(data = cv, ggplot2::aes(x = x, y = fit), colour = "gray10",
                           linewidth = 1)
    if (has_ip)
        p <- p + ggplot2::geom_vline(xintercept = ip$ip_x, colour = COL_IP, linetype = "dashed",
                                     linewidth = 0.6) +
            ggplot2::geom_hline(yintercept = ip$ip_y, colour = COL_IPV, linetype = "dashed",
                                linewidth = 0.5)
    p + ggplot2::scale_colour_manual(values = cols, limits = names(cols), name = NULL,
                                     labels = with_counts(names(cols), d$grp), drop = FALSE) +
        ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2.2, alpha = 1))) +
        ggplot2::labs(title = title %|||% "CEP-IP GAM fit", subtitle = sub,
                      x = "GOI", y = "Composite") +
        theme_cepip() + ggplot2::theme(legend.position = "bottom")
}

# =============================================================================
# 4. Residual GAM R(x)
# =============================================================================

#' Plot the residual GAM R(x) and the IP
#'
#' R(x) with its 95% band, the residuals of the top-EP cells it was fitted to
#' (faint points), the zero line and the IPx with its CI. The subtitle gives
#' the IPRS, its tier and the components C1 to C5.
#'
#' @param ip A `cepip_ip` object, or a `cepip` object (its `ip` is used).
#' @param title Optional plot title.
#' @return A ggplot object.
#' @examples
#' d <- cepip_toy
#' res <- cepip(d$x, d$y, donor = d$donor, id = d$id, verbose = FALSE)
#' plot_residual(res)
#' @export
plot_residual <- function(ip, title = NULL) {
    if (inherits(ip, "cepip")) ip <- ip$ip
    check_class(ip, "cepip_ip", "ip")
    if (is.null(ip$grid)) stop("The IP object has no R(x) grid (", ip$status, ").", call. = FALSE)
    gr <- ip$grid
    rd <- ip$residual_data
    rng <- range(c(gr$r_lower, gr$r_upper, stats::quantile(rd$residual, c(0.01, 0.99))),
                 na.rm = TRUE)
    p <- ggplot2::ggplot()
    if (isTRUE(ip$has_ip) && all(is.finite(ip$ip_ci)))
        p <- p + ggplot2::annotate("rect", xmin = ip$ip_ci[1], xmax = ip$ip_ci[2],
                                   ymin = -Inf, ymax = Inf, fill = COL_IP, alpha = 0.10)
    p <- p +
        ggplot2::geom_point(data = rd, ggplot2::aes(x = GOI, y = residual), colour = "gray45",
                            size = 0.6, alpha = 0.2) +
        ggplot2::geom_hline(yintercept = 0, colour = "gray20", linewidth = 0.4) +
        ggplot2::geom_ribbon(data = gr, ggplot2::aes(x = x, ymin = r_lower, ymax = r_upper),
                             fill = COL_FLAT, alpha = 0.22) +
        ggplot2::geom_line(data = gr, ggplot2::aes(x = x, y = r_hat), colour = COL_FLAT,
                           linewidth = 1)
    s <- ip$iprs
    comp <- sprintf("C1 %.2f, C2 %.2f, C3 %.2f, C4 %.2f, C5 %.2f", s$C1, s$C2, s$C3, s$C4, s$C5)
    if (isTRUE(ip$has_ip)) {
        p <- p + ggplot2::geom_vline(xintercept = ip$ip_x, colour = COL_IP, linetype = "dashed",
                                     linewidth = 0.6)
        sub <- sprintf("IPx = %.3f (95%% CI %.3f to %.3f). IPRS = %.3f (%s).\n%s. %s top-EP cells.",
                       ip$ip_x, ip$ip_ci[1], ip$ip_ci[2], s$IPRS, s$tier, comp,
                       fmt_n(nrow(rd)))
    } else sub <- sprintf("No IP: %s.\n%s.", ip$status, comp)
    p + ggplot2::coord_cartesian(ylim = rng) +
        ggplot2::labs(title = title %|||% "Residual GAM R(x)", subtitle = sub,
                      x = "GOI", y = "Residual y - f(x)") +
        theme_cepip()
}

# =============================================================================
# 5. Pitcher plot (score space)
# =============================================================================

#' Pitcher plot of the SKCS scores
#'
#' K_i (y) against S_i,norm (x) for every cell, with the S_i,norm and K_i
#' gates. SKCS-core cells are coloured by the four groups when `ip` is given,
#' otherwise in one colour.
#'
#' @param skcs A `cepip_skcs` object, or a `cepip` object (its `skcs` is used,
#'   and its band or IP gives the groups).
#' @param ip Optional `cepip_band`, `cepip_cmrf` or `cepip_ip` object of the
#'   same cells, for the four-group colours.
#' @param title Optional plot title.
#' @return A ggplot object.
#' @examples
#' d <- cepip_toy
#' res <- cepip(d$x, d$y, donor = d$donor, id = d$id, skcs = TRUE, verbose = FALSE)
#' plot_pitcher(res)
#' @export
plot_pitcher <- function(skcs, ip = NULL, title = NULL) {
    if (inherits(skcs, "cepip")) {
        if (is.null(ip)) ip <- skcs$band %|||% skcs$ip
        skcs <- skcs$skcs
    }
    check_class(skcs, "cepip_skcs", "skcs")
    cl <- skcs$cells
    grp <- rep(NA_character_, nrow(cl))
    if (inherits(ip, "cepip_band") || inherits(ip, "cepip_cmrf")) {
        grp <- ip$cells$group_orig[match(cl$id, ip$cells$id)]
    } else if (inherits(ip, "cepip_ip") && isTRUE(ip$has_ip) && !is.null(skcs$ep)) {
        trep <- rep(FALSE, nrow(cl))
        trep[cl$is_core] <- skcs$ep$is_TREP
        side <- ip_labels(cl$x, cl$y, ip, ip$ip_axis)
        grp <- paste(side, ifelse(trep, "TREP", "non-TREP"))
    }
    four <- !all(is.na(grp[cl$is_core]))
    d <- data.frame(Si_norm = cl$Si_norm, Ki = cl$Ki, core = cl$is_core)
    if (four) {
        d$grp <- ifelse(cl$is_core, grp, "Non-SKCS-core")
        cols <- c(COL_FOUR, "Non-SKCS-core" = COL_GRAY)
    } else {
        d$grp <- ifelse(cl$is_core, "SKCS-core", "Non-SKCS-core")
        cols <- c("SKCS-core" = COL_SINGLE, "Non-SKCS-core" = COL_GRAY)
    }
    d$grp <- factor(d$grp, levels = names(cols))
    d <- d[order(d$core), ]
    ggplot2::ggplot(d, ggplot2::aes(x = Si_norm, y = Ki, colour = grp)) +
        ggplot2::geom_point(data = d[!d$core, ], size = 0.6, alpha = 0.35) +
        ggplot2::geom_point(data = d[d$core, ], size = 0.9, alpha = 0.55) +
        ggplot2::geom_vline(xintercept = skcs$s_min, colour = "#009E73", linewidth = 0.6) +
        ggplot2::geom_hline(yintercept = skcs$k_min, colour = "#009E73", linewidth = 0.6) +
        ggplot2::annotate("text", x = skcs$s_min, y = -Inf, hjust = -0.05, vjust = -0.6,
                          label = sprintf("Si > %.2f", skcs$s_min), size = 3, fontface = "bold") +
        ggplot2::annotate("text", x = -Inf, y = skcs$k_min, hjust = -0.05, vjust = -0.5,
                          label = sprintf("Ki > %.2f", skcs$k_min), size = 3, fontface = "bold") +
        ggplot2::scale_colour_manual(values = cols, limits = names(cols), drop = FALSE,
                                     labels = with_counts(names(cols), d$grp), name = NULL) +
        ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2.2, alpha = 1))) +
        ggplot2::labs(title = title %|||% "Pitcher plot",
                      subtitle = sprintf("SKCS-core: %s of %s cells. Spearman rho = %.3f, Kendall tau = %.3f.",
                                         fmt_n(sum(cl$is_core)), fmt_n(nrow(cl)), skcs$rho, skcs$tau),
                      x = "Spearman contribution Si (normalized)",
                      y = "Kendall contribution Ki") +
        theme_cepip() + legend_inside(c(0.99, 0.01), c(1, 0))
}

# =============================================================================
# 6. Copula space
# =============================================================================

#' Plot the cells and the CMRF band in copula space
#'
#' (u, v) of all cells, the SKCS-core cells coloured by their four groups, the
#' spine g(u), the flat delta_base band, the adaptive delta_adapt band and the
#' IP_u and IP_v lines.
#'
#' @param band A `cepip_band` object, or a `cepip_cmrf` or `cepip` object.
#' @param title Optional plot title.
#' @return A ggplot object.
#' @examples
#' \donttest{
#' d <- cepip_toy
#' # CMRF with a short gate of 5 repeats (the default is 40)
#' res <- cepip(d$x, d$y, donor = d$donor, id = d$id, cmrf = TRUE, repeats = 5,
#'              tests = FALSE, verbose = FALSE)
#' plot_copula(res)
#' }
#' @export
plot_copula <- function(band, title = NULL) {
    rb <- resolve_band(band)
    b <- rb$b; cl <- rb$cells
    gr <- band_grid(b)
    p <- ggplot2::ggplot() + layer_flat(gr)[1] + layer_adapt_fill(gr) +
        layer_cells(cl) + layer_spine(gr) + layer_flat(gr)[2:3] + layer_adapt_border(gr) +
        layer_ip(b) +
        scale_four(guide = "legend", labels = with_counts(names(COL_FOUR), cl$group_orig)) +
        ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2.2, alpha = 1)))
    uv_panel(p, title %|||% "Copula space",
             sprintf("delta_base = %.4f, kappa = %.2f, bw = %.2f. IP_u = %.3f, IP_v = %.3f.",
                     b$delta_base, b$kappa, b$bw, b$IP_u, b$IP_v)) +
        legend_inside(c(0.99, 0.01), c(1, 0))
}

# =============================================================================
# 7. Band formation (Figure 6)
# =============================================================================

#' Plot the formation of the CMRF band (Figure 6)
#'
#' Panels: (1) the spine fitted on the SKCS-core TREP cells, (2) all SKCS-core
#' cells in their four groups, (3) the flat delta_base band, (4) the adaptive
#' delta_adapt band over delta_base, a QC panel zoomed to IP_u +- 0.25 (QC1 the
#' width at IP_u, QC2 the widths at IP_u +- bw, where 60.65 percent of the
#' widening is left, and QC3 the tangents of the band edges at IP_u +- bw),
#' (5) the Stage 1 candidates, (6) the SKCS-core cells outside the band
#' (Original only) and, when `cmrf` is given, (7) the cells recovered by the
#' gate and those it removed and (8) the final Module. This plot is useful to
#' choose `kappa` and `bw` with [cepip_band()] before running [cepip_cmrf()].
#'
#' @param band A `cepip_band` object, or a `cepip_cmrf` or `cepip` object.
#' @param cmrf Optional `cepip_cmrf` object; adds panels 7 and 8.
#' @param title Optional overall title.
#' @return A gtable of class `cepip_figure` (printing it draws it; it can
#'   also be drawn with [grid::grid.draw()]). The attribute `cepip_size`
#'   gives a suitable width and height in inches.
#' @examples
#' \donttest{
#' d <- cepip_toy
#' # CMRF with a short gate of 5 repeats (the default is 40)
#' res <- cepip(d$x, d$y, donor = d$donor, id = d$id, cmrf = TRUE, repeats = 5,
#'              tests = FALSE, verbose = FALSE)
#' plot_band(res)
#' }
#' @export
plot_band <- function(band, cmrf = NULL, title = NULL) {
    rb <- resolve_band(band, cmrf)
    b <- rb$b; cl <- rb$cells; cmrf <- rb$cmrf
    gr <- band_grid(b)
    n_core <- sum(cl$is_core)
    cand_flat <- cl[cl$is_candidate & cl$resid_abs <= b$delta_base, , drop = FALSE]
    cand_gap <- cl[cl$is_candidate & cl$resid_abs > b$delta_base, , drop = FALSE]
    oo <- cl[cl$is_original_only, , drop = FALSE]
    pts_cand <- list(
        ggplot2::geom_point(data = cand_flat, ggplot2::aes(x = u, y = v), shape = 21, size = 2.2,
                            fill = COL_BLUE_CAND, colour = "white", stroke = 0.4, alpha = 0.88),
        ggplot2::geom_point(data = cand_gap, ggplot2::aes(x = u, y = v), shape = 21, size = 2.2,
                            fill = COL_GAP_CAND, colour = "white", stroke = 0.4, alpha = 0.88))
    pts_oo <- ggplot2::geom_point(data = oo, ggplot2::aes(x = u, y = v), shape = 25, size = 2.2,
                                  fill = COL_ORIONLY, colour = "white", stroke = 0.4, alpha = 0.85)

    # 1 Spine on the SKCS-core TREP cells
    tr <- cl[cl$is_core & cl$is_TREP, , drop = FALSE]
    p1 <- ggplot2::ggplot() +
        ggplot2::geom_point(data = cl, ggplot2::aes(x = u, y = v), colour = COL_BG,
                            size = 0.5, alpha = 0.35) +
        ggplot2::geom_point(data = tr, ggplot2::aes(x = u, y = v, colour = group_orig),
                            size = 1.1, alpha = 0.5) +
        layer_spine(gr) + layer_ip(b) + scale_four()
    p1 <- uv_panel(p1, "1 Spline calibrated by the SKCS-core TREP cells",
                   sprintf("Spine fitted on %s TREP cells (k = %d).", fmt_n(nrow(tr)), b$spine$k))

    # 2 All SKCS-core cells
    lev2 <- c(names(COL_FOUR), "Non-SKCS-core")
    g2 <- ifelse(cl$group_orig %in% names(COL_FOUR), cl$group_orig, "Non-SKCS-core")
    p2 <- ggplot2::ggplot() + layer_cells(cl) + layer_spine(gr) + layer_ip(b) +
        ggplot2::scale_colour_manual(values = c(COL_FOUR, "Non-SKCS-core" = COL_GRAY),
                                     limits = lev2, labels = with_counts(lev2, g2),
                                     name = NULL, drop = FALSE) +
        ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2, alpha = 1)))
    p2 <- uv_panel(p2, "2 Overlay of all SKCS-core cells",
                   sprintf("%s SKCS-core cells of %s.", fmt_n(n_core), fmt_n(nrow(cl)))) +
        legend_inside(c(0.99, 0.01), c(1, 0)) +
        ggplot2::theme(legend.text = ggplot2::element_text(size = 7),
                       legend.key.size = grid::unit(0.3, "cm"))

    # 3 Flat band
    p3 <- ggplot2::ggplot() + layer_flat(gr)[1] + layer_cells(cl) + layer_spine(gr) +
        layer_flat(gr)[2:3] + layer_ip(b) + scale_four()
    p3 <- uv_panel(p3, sprintf("3 delta_base band covering %.0f percent of SKCS-core cells",
                               100 * b$coverage),
                   sprintf("g(u) +- delta_base, delta_base = %.4f (donor-weighted).", b$delta_base))

    # QC zoom window
    zx <- c(max(0, b$IP_u - 0.25), min(1, b$IP_u + 0.25))
    in_z <- gr$u >= zx[1] & gr$u <= zx[2]
    zy <- c(max(-0.03, min(gr$a_lo[in_z]) - 0.08), min(1.03, max(gr$a_up[in_z]) + 0.05))

    # 4 Adaptive band over the flat band
    p4 <- ggplot2::ggplot() + layer_flat(gr)[1] + layer_adapt_fill(gr) + layer_cells(cl) +
        layer_spine(gr) + layer_adapt_border(gr) + layer_flat(gr)[2:3] + layer_ip(b) +
        ggplot2::annotate("rect", xmin = zx[1], xmax = zx[2], ymin = zy[1], ymax = zy[2],
                          fill = NA, colour = "black", linewidth = 0.4) +
        ggplot2::annotate("text", x = b$IP_u, y = zy[1] + 0.03, label = "IP_u", colour = COL_IP,
                          size = 2.8, fontface = "bold", hjust = -0.1) +
        scale_four()
    p4 <- uv_panel(p4, "4 delta_adapt band overlaid on delta_base",
                   sprintf("delta_adapt(u) = delta_base x (1 + %.2f K(u)), widest at IP_u. Box: QC zoom.",
                           b$kappa))

    # QC panel
    q <- b$qc
    at <- function(col, x) stats::approx(gr$u, gr[[col]], xout = x, rule = 2)$y
    uq <- c(b$IP_u - b$bw, b$IP_u, b$IP_u + b$bw)
    eps <- 0.005
    slope <- function(col, x) (at(col, x + eps) - at(col, x - eps)) / (2 * eps)
    th <- 0.04
    tang <- do.call(rbind, lapply(c("a_up", "a_lo"), function(col) {
        xs <- uq[c(1, 3)]
        y0 <- at(col, xs); s <- slope(col, xs)
        data.frame(x0 = xs - th, x1 = xs + th, y0 = y0 - s * th, y1 = y0 + s * th,
                   u = xs, v = y0)
    }))
    arw <- grid::arrow(angle = 20, length = grid::unit(0.07, "inches"), ends = "both",
                       type = "closed")
    off <- 0.008
    pq <- ggplot2::ggplot() + layer_flat(gr)[1] + layer_adapt_fill(gr) + layer_cells(cl, 1.3, 0.45) +
        layer_spine(gr) + layer_adapt_border(gr, 0.55) + layer_flat(gr)[2:3] + layer_ip(b) +
        ggplot2::geom_vline(xintercept = uq[c(1, 3)], colour = "gray50", linewidth = 0.35,
                            linetype = "dotdash") +
        ggplot2::annotate("segment", x = uq[2], xend = uq[2], y = at("a_lo", uq[2]),
                          yend = at("a_up", uq[2]), colour = COL_QC[1], linewidth = 0.7,
                          arrow = arw) +
        ggplot2::annotate("segment", x = uq[c(1, 3)] + c(-off, off), xend = uq[c(1, 3)] + c(-off, off),
                          y = at("a_lo", uq[c(1, 3)]), yend = at("a_up", uq[c(1, 3)]),
                          colour = COL_QC[2], linewidth = 0.55, arrow = arw) +
        ggplot2::geom_segment(data = tang, ggplot2::aes(x = x0, xend = x1, y = y0, yend = y1),
                              colour = COL_QC[3], linewidth = 0.5) +
        ggplot2::geom_point(data = tang, ggplot2::aes(x = u, y = v), colour = COL_QC[3], size = 2) +
        ggplot2::annotate("label", x = uq[2], y = mean(c(at("a_lo", uq[2]), at("a_up", uq[2]))),
                          label = "QC1", colour = COL_QC[1], size = 2.8, fontface = "bold",
                          fill = grDevices::adjustcolor("white", 0.85), label.size = 0) +
        ggplot2::annotate("label", x = uq[c(1, 3)] + c(-off, off),
                          y = (at("a_lo", uq[c(1, 3)]) + at("g", uq[c(1, 3)])) / 2,
                          label = "QC2", colour = COL_QC[2], size = 2.6, fontface = "bold",
                          fill = grDevices::adjustcolor("white", 0.85), label.size = 0) +
        ggplot2::annotate("text", x = tang$x0[1:2], y = tang$y0[1:2], label = "QC3",
                          colour = COL_QC[3], size = 2.6, fontface = "bold", hjust = 1.1) +
        ggplot2::annotate("text", x = uq[c(1, 3)], y = zy[1], label = c("-bw", "+bw"),
                          colour = "gray25", size = 2.8, vjust = 0, hjust = -0.15) +
        scale_four()
    pq <- uv_panel(pq, "QC metrics of the delta_adapt band",
                   sprintf(paste0("QC1 width at IP_u = %.4f (flat band %.4f).\n",
                                  "QC2 width at IP_u +- bw = %.4f (%.2f percent of the widening left).\n",
                                  "QC3 tangents at IP_u +- bw (largest kernel slope %.2f)."),
                           q$width_at_ip, q$width_base, q$width_at_bw,
                           100 * q$widening_left_at_bw, q$max_kernel_slope),
                   xlim = zx, ylim = zy, xbreaks = seq(0, 1, 0.1)) +
        ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 7.5))

    # 5 Stage 1 candidates
    p5 <- ggplot2::ggplot() + layer_flat(gr)[1] + layer_adapt_fill(gr) + layer_cells(cl) +
        layer_spine(gr) + pts_cand + layer_adapt_border(gr) + layer_flat(gr)[2:3] +
        layer_ip(b) + scale_four()
    p5 <- uv_panel(p5, "5 Stage 1 candidates captured by delta_base and delta_adapt",
                   sprintf("Blue: non-SKCS cells within delta_base (n=%s).\nDark: in delta_adapt but outside delta_base (n=%s).",
                           fmt_n(nrow(cand_flat)), fmt_n(nrow(cand_gap))))

    # 6 Original only
    p6 <- ggplot2::ggplot() + layer_flat(gr)[1] + layer_adapt_fill(gr) + layer_cells(cl) +
        layer_spine(gr) + pts_cand + pts_oo + layer_adapt_border(gr) + layer_flat(gr)[2:3] +
        layer_ip(b) + scale_four()
    p6 <- uv_panel(p6, "6 SKCS-core cells outside delta_adapt leave the Module",
                   sprintf("Tan triangles: Original only (n=%s).\nCandidates (n=%s) go to the RSk-fCV gate.",
                           fmt_n(nrow(oo)), fmt_n(sum(cl$is_candidate))))
    panels <- list(p1, p2, p3, p4, pq, p5, p6)
    lay <- rbind(1:3, 4:6, c(7, NA, NA))
    if (!is.null(cmrf)) {
        rec <- cl[cl$cell_class == "Recovered", , drop = FALSE]
        rem <- cl[cl$removed_by_stage2, , drop = FALSE]
        p7 <- ggplot2::ggplot() + layer_flat(gr)[1] + layer_adapt_fill(gr) + layer_cells(cl) +
            layer_spine(gr) + pts_oo +
            ggplot2::geom_point(data = rem, ggplot2::aes(x = u, y = v), shape = 4, size = 2,
                                colour = COL_REMOVED_X, stroke = 0.7, alpha = 0.75) +
            ggplot2::geom_point(data = rec, ggplot2::aes(x = u, y = v), shape = 23, size = 2.2,
                                fill = COL_RECOV, colour = "white", stroke = 0.4, alpha = 0.88) +
            layer_adapt_border(gr) + layer_flat(gr)[2:3] + layer_ip(b) + scale_four()
        p7 <- uv_panel(p7, "7 Cells recovered by the RSk-fCV gate",
                       sprintf("Purple diamonds: Recovered (n=%s).\nCrosses: removed by the gate (n=%s).",
                               fmt_n(nrow(rec)), fmt_n(nrow(rem))))
        lev8 <- c(names(COL_FOUR), "Recovered", "Non-Module cells")
        c8 <- ifelse(cl$cell_class == "Original + Refined", cl$group_orig,
              ifelse(cl$cell_class == "Recovered", "Recovered", "Non-Module cells"))
        d8 <- data.frame(u = cl$u, v = cl$v, grp = factor(c8, levels = lev8))
        nm <- d8[d8$grp == "Non-Module cells", , drop = FALSE]
        ori <- d8[d8$grp %in% names(COL_FOUR), , drop = FALSE]
        rc <- d8[d8$grp == "Recovered", , drop = FALSE]
        lab8 <- with_counts(lev8, c8)
        p8 <- ggplot2::ggplot() + layer_adapt_fill(gr, "gray70", 0.25) +
            ggplot2::geom_point(data = nm, ggplot2::aes(x = u, y = v, colour = grp, shape = grp,
                                                        fill = grp), size = 0.5, alpha = 0.25) +
            ggplot2::geom_point(data = ori, ggplot2::aes(x = u, y = v, colour = grp, shape = grp,
                                                         fill = grp), size = 1.1, alpha = 0.5) +
            layer_spine(gr) +
            ggplot2::geom_point(data = rc, ggplot2::aes(x = u, y = v, colour = grp, shape = grp,
                                                        fill = grp), size = 2.2, stroke = 0.4,
                                alpha = 0.88) +
            layer_adapt_border(gr, 0.75) + layer_ip(b) +
            ggplot2::scale_colour_manual(values = c(COL_FOUR, "Recovered" = "white",
                                                    "Non-Module cells" = COL_GRAY),
                                         limits = lev8, labels = lab8, name = "Final Module",
                                         drop = FALSE) +
            ggplot2::scale_fill_manual(values = c(COL_FOUR, "Recovered" = COL_RECOV,
                                                  "Non-Module cells" = COL_GRAY),
                                       limits = lev8, labels = lab8, name = "Final Module",
                                       drop = FALSE) +
            ggplot2::scale_shape_manual(values = c(16, 16, 16, 16, 23, 16), limits = lev8,
                                        labels = lab8, name = "Final Module", drop = FALSE) +
            merged_guides()
        p8 <- uv_panel(p8, "8 Final Module with the recovered cells",
                       sprintf("Module = OriRef (n=%s) plus Recovered (n=%s).",
                               fmt_n(sum(cl$cell_class == "Original + Refined")), fmt_n(nrow(rc)))) +
            legend_inside(c(0.99, 0.01), c(1, 0)) +
            ggplot2::theme(legend.text = ggplot2::element_text(size = 7),
                           legend.title = ggplot2::element_text(size = 7.5),
                           legend.key.size = grid::unit(0.3, "cm"))
        panels <- c(panels, list(p7, p8))
        lay <- rbind(1:3, 4:6, 7:9)
    }
    combine_panels(panels, lay,
                   title %|||% sprintf("CMRF band formation (n=%s cells)", fmt_n(nrow(cl))),
                   width = 15, height = if (is.null(cmrf)) 15.5 else 15.5)
}

# =============================================================================
# 8. RSk-fCV gate (Figure 7)
# =============================================================================

#' Plot the RSk-fCV gate of CMRF Stage 2 (Figure 7)
#'
#' Panels for the first or the last gate iteration: (A) the fold layout of
#' repeat 1, (B) the u strata of the pool, (C) curve A (isotonic fit on the
#' OriRef cells) and curve B (on the pool) with the candidates' residual
#' sticks, (D) the per-fold RMSE of repeat 1, (E) the RMSE_A and RMSE_B
#' heatmaps (repeat x fold), (F) the RMSE distributions with the Wilcoxon p and
#' the paired differences against BV, (G) the candidates ranked by
#' |v - PAVA(u)|, (I) the gate convergence over all iterations and (Q) the
#' class counts by IP side. When the class tests were run, (P) |v - g(u)| by
#' class with the Jonckheere-Terpstra p and (O) the Recovered cells with the
#' Hoeffding D and distance correlation p are added.
#'
#' The pool of the first iteration is the OriRef cells plus all candidates;
#' the pool of the last iteration is the OriRef cells plus the candidates
#' still active (the Recovered cells).
#'
#' @param cmrf A `cepip_cmrf` object, or a `cepip` object.
#' @param iteration `"first"` or `"last"`.
#' @param title Optional overall title.
#' @return A gtable of class `cepip_figure` (printing it draws it; it can
#'   also be drawn with [grid::grid.draw()]). The attribute `cepip_size`
#'   gives a suitable width and height in inches.
#' @examples
#' \donttest{
#' d <- cepip_toy
#' # CMRF with a short gate of 5 repeats (the default is 40)
#' res <- cepip(d$x, d$y, donor = d$donor, id = d$id, cmrf = TRUE, repeats = 5,
#'              tests = FALSE, verbose = FALSE)
#' plot_rskfcv(res, "first")
#' plot_rskfcv(res, "last")
#' }
#' @export
plot_rskfcv <- function(cmrf, iteration = c("first", "last"), title = NULL) {
    if (inherits(cmrf, "cepip")) cmrf <- cmrf$cmrf
    check_class(cmrf, "cepip_cmrf", "cmrf")
    iteration <- match.arg(iteration)
    gt <- cmrf$gate; cl <- cmrf$cells
    K <- gt$folds; R <- gt$repeats
    a <- gt$anchor_row; cand <- gt$cand_row
    nc <- length(cand)
    first <- iteration == "first"
    act <- if (first) rep(TRUE, nc) else gt$pass
    it <- if (nc == 0L) NA_integer_ else if (first) 1L else length(gt$p)
    it_lab <- if (first) "first iteration" else sprintf("iteration %d (last)", it)
    pool <- c(a, cand[act])
    u_p <- cl$u[pool]; v_p <- cl$v[pool]
    is_a <- seq_along(pool) <= length(a)
    fm <- if (first) gt$folds_B_first else gt$folds_B_last
    rB <- if (first) gt$rmse_B_first else gt$rmse_B_last
    rA <- gt$rmse_A
    p_it <- if (is.na(it)) NA_real_ else gt$p[it]
    e_it <- if (is.na(it)) NA_real_ else gt$effect[it]

    # A Fold layout, repeat 1
    nte <- if (!is.null(fm)) as.integer(table(factor(fm[, 1], levels = seq_len(K))))
           else rep(NA_integer_, K)
    dA <- data.frame(fold = factor(rep(sprintf("Fold %d", seq_len(K)), 2),
                                   levels = rev(sprintf("Fold %d", seq_len(K)))),
                     role = factor(rep(c("Train", "Test"), each = K), levels = c("Test", "Train")),
                     count = c(length(pool) - nte, nte))
    pA <- ggplot2::ggplot(dA, ggplot2::aes(x = count, y = fold, fill = role)) +
        ggplot2::geom_col(width = 0.7, colour = "white", linewidth = 0.5) +
        ggplot2::geom_text(ggplot2::aes(label = count), position = ggplot2::position_stack(vjust = 0.5),
                           size = 3, colour = ifelse(dA$role == "Train", "white", "gray20"),
                           fontface = "bold") +
        ggplot2::scale_fill_manual(values = c(Train = COL_TRAIN, Test = COL_TESTF), name = NULL,
                                   breaks = c("Train", "Test")) +
        ggplot2::labs(title = "A Fold layout, repeat 1",
                      subtitle = sprintf("Pool = OriRef (%s) + candidates (%s) = %s cells.",
                                         fmt_n(length(a)), fmt_n(sum(act)), fmt_n(length(pool))),
                      x = "Cells", y = NULL) +
        theme_cepip() + ggplot2::theme(legend.position = "bottom", panel.grid.major.y = ggplot2::element_blank())

    # B Strata of u
    n <- length(u_p)
    ns <- max(1L, min(10L, floor(n / K)))
    st <- pmin(ceiling(rank(u_p, ties.method = "first") / n * ns), ns)
    dB <- data.frame(u = u_p, v = v_p, stratum = factor(sprintf("S%d", st), levels = sprintf("S%d", seq_len(ns))))
    brk <- vapply(seq_len(ns - 1L), function(s) max(u_p[st == s]), numeric(1))
    pB <- ggplot2::ggplot(dB, ggplot2::aes(x = u, y = v, colour = stratum)) +
        ggplot2::geom_vline(xintercept = brk, colour = "gray60", linetype = "dotted", linewidth = 0.3) +
        ggplot2::geom_point(size = 0.7, alpha = 0.65) +
        ggplot2::scale_colour_manual(values = stats::setNames(COL_STRATA[seq_len(ns)], levels(dB$stratum)),
                                     name = "Stratum") +
        ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2, override.aes = list(size = 2))) +
        ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
        ggplot2::labs(title = "B Strata of u in the pool",
                      subtitle = sprintf("%d strata; each stratum is spread over the %d folds.", ns, K),
                      x = "u", y = "v") +
        theme_cepip() + theme_vtitle() + ggplot2::theme(aspect.ratio = 1, legend.position = "bottom")

    # C Curves A and B with the candidate residual sticks
    step_df <- function(uu, vv) { o <- order(uu); data.frame(uu = uu[o], yy = stats::isoreg(uu[o], vv[o])$yf) }
    sA <- step_df(u_p[is_a], v_p[is_a])
    sB <- step_df(u_p, v_p)
    dc <- data.frame(uu = u_p[!is_a], vv = v_p[!is_a])
    dc$ff <- if (nrow(dc)) iso_predict(u_p, v_p, dc$uu) else numeric(0)
    pC <- ggplot2::ggplot() +
        ggplot2::geom_point(data = data.frame(uu = u_p[is_a], vv = v_p[is_a]),
                            ggplot2::aes(x = uu, y = vv), colour = COL_GRAY, size = 0.7, alpha = 0.35) +
        ggplot2::geom_segment(data = dc, ggplot2::aes(x = uu, xend = uu, y = vv, yend = ff),
                              colour = COL_STICK, linewidth = 0.35, alpha = 0.8) +
        ggplot2::geom_point(data = dc, ggplot2::aes(x = uu, y = vv), colour = COL_RECOV,
                            size = 1, alpha = 0.55) +
        ggplot2::geom_step(data = sA, ggplot2::aes(x = uu, y = yy), colour = COL_A, linewidth = 0.9) +
        ggplot2::geom_step(data = sB, ggplot2::aes(x = uu, y = yy), colour = COL_B, linewidth = 0.9)
    if (first && length(gt$removal_order)) {
        w <- gt$removal_order[1]
        dw <- data.frame(uu = cl$u[cand[w]], vv = cl$v[cand[w]])
        pC <- pC + ggplot2::geom_point(data = dw, ggplot2::aes(x = uu, y = vv), shape = 4, size = 3.5,
                                       stroke = 1.4, colour = COL_TEST)
    }
    pC <- pC + ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
        ggplot2::labs(title = "C Curve A and curve B",
                      subtitle = paste0("Green: PAVA on OriRef (A). Navy: PAVA on the pool (B).\n",
                                        "Purple: candidates, sticks |v - PAVA(u)|",
                                        if (first && length(gt$removal_order)) ". Red cross: first removed." else "."),
                      x = "u", y = "v") +
        theme_cepip() + theme_vtitle() + ggplot2::theme(aspect.ratio = 1)

    # D RMSE per fold, repeat 1
    dD <- data.frame(fold = factor(rep(seq_len(K), 2)),
                     curve = factor(rep(c("Curve A", "Curve B"), each = K), levels = c("Curve A", "Curve B")),
                     rmse = c(rA[seq_len(K)], rB[seq_len(K)]))
    pD <- ggplot2::ggplot(dD, ggplot2::aes(x = fold, y = rmse, fill = curve)) +
        ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7), width = 0.6) +
        ggplot2::scale_fill_manual(values = c("Curve A" = COL_A, "Curve B" = COL_B), name = NULL) +
        ggplot2::labs(title = "D RMSE per fold, repeat 1",
                      subtitle = "One test RMSE per fold for each curve.", x = "Fold", y = "RMSE") +
        theme_cepip() + ggplot2::theme(legend.position = "bottom")

    # E Heatmaps
    heat <- function(vals, pal, ttl) {
        dd <- data.frame(fold_id = factor(rep(seq_len(K), R)),
                         rep_id = factor(rep(seq_len(R), each = K), levels = rev(seq_len(R))),
                         rmse = vals)
        ybr <- unique(c(1, seq(5, R, by = 5)))
        ggplot2::ggplot(dd, ggplot2::aes(x = fold_id, y = rep_id, fill = rmse)) +
            ggplot2::geom_tile(colour = "black", linewidth = 0.12) +
            ggplot2::scale_fill_gradientn(colours = rev(grDevices::hcl.colors(256, pal)), name = "RMSE") +
            ggplot2::scale_x_discrete(expand = c(0, 0)) +
            ggplot2::scale_y_discrete(expand = c(0, 0), breaks = as.character(ybr)) +
            ggplot2::labs(title = ttl, x = "Fold", y = "Repeat") +
            theme_cepip() + ggplot2::theme(panel.grid = ggplot2::element_blank(),
                                           panel.border = ggplot2::element_blank(),
                                           axis.ticks = ggplot2::element_blank(),
                                           legend.key.width = grid::unit(0.3, "cm"))
    }
    pE <- gridExtra::arrangeGrob(heat(rA, "Teal", "E RMSE_A"), heat(rB, "Reds", "RMSE_B"), ncol = 2)

    # F Distributions and paired differences
    dF <- data.frame(curve = factor(rep(c("Curve A", "Curve B"), each = length(rA)),
                                    levels = c("Curve A", "Curve B")), rmse = c(rA, rB))
    dF <- dF[is.finite(dF$rmse), ]
    pF1 <- ggplot2::ggplot(dF, ggplot2::aes(x = curve, y = rmse, fill = curve)) +
        ggplot2::geom_jitter(width = 0.2, size = 1, alpha = 0.35, shape = 21, stroke = 0.25,
                             colour = "gray35") +
        ggplot2::geom_boxplot(width = 0.55, colour = "black", linewidth = 0.5, outlier.shape = NA,
                              alpha = 0.35) +
        ggplot2::scale_fill_manual(values = c(COL_A, COL_RECOV), guide = "none") +
        ggplot2::labs(title = "F RMSE by curve",
                      subtitle = sprintf("Paired Wilcoxon (B > A)\np = %s", fmt_p(p_it)),
                      x = NULL, y = "Test RMSE") +
        theme_cepip()
    med_a <- gt$median_rmse_A
    dd <- data.frame(g = "B - A", d = (rB - rA) / med_a)
    dd <- dd[is.finite(dd$d), ]
    pF2 <- ggplot2::ggplot(dd, ggplot2::aes(x = g, y = d)) +
        ggplot2::geom_hline(yintercept = 0, colour = "gray30", linewidth = 0.4) +
        ggplot2::geom_hline(yintercept = gt$bv, colour = COL_A, linetype = "dashed", linewidth = 0.7) +
        ggplot2::geom_jitter(width = 0.2, size = 1, alpha = 0.35, shape = 21, stroke = 0.25,
                             colour = "gray35", fill = COL_TEST) +
        ggplot2::geom_boxplot(width = 0.55, colour = "black", linewidth = 0.5, outlier.shape = NA,
                              alpha = 0.35, fill = COL_TEST) +
        ggplot2::annotate("text", x = 1.45, y = gt$bv, label = "BV", colour = COL_A, vjust = -0.5,
                          size = 3, fontface = "bold") +
        ggplot2::labs(title = "Paired differences",
                      subtitle = sprintf("dRMSE = %.4f\nBV = %.4f", e_it, gt$bv),
                      x = NULL, y = "(RMSE_B - RMSE_A) / median(RMSE_A)") +
        theme_cepip()
    pF <- gridExtra::arrangeGrob(pF1, pF2, ncol = 2)

    # G Candidate ranking
    if (nc > 0L) {
        o <- order(gt$resid_pava, decreasing = TRUE)
        hl <- if (first) seq_len(nc) %in% gt$removal_order[seq_len(min(1L, length(gt$removal_order)))]
              else !gt$pass
        dG <- data.frame(rank_i = seq_len(nc), r = gt$resid_pava[o],
                         removed = factor(ifelse(hl[o], "Removed", "Kept"), levels = c("Removed", "Kept")))
        pG <- ggplot2::ggplot(dG, ggplot2::aes(x = rank_i, y = r, fill = removed)) +
            ggplot2::geom_col(width = 1) +
            ggplot2::scale_fill_manual(values = c(Removed = COL_BAR_REM, Kept = COL_BAR_ACC), name = NULL,
                                       drop = FALSE) +
            ggplot2::labs(title = "G Candidates ranked by PAVA residual",
                          subtitle = if (first) "Black: the first candidate removed."
                                     else sprintf("Black: removed (%s of %s). Purple: kept.",
                                                  fmt_n(gt$n_removed), fmt_n(nc)),
                          x = "Rank, descending |v - PAVA(u)|", y = "|v - PAVA(u)|") +
            theme_cepip() + ggplot2::theme(legend.position = "bottom")
        if (first && any(hl))
            pG <- pG + ggplot2::geom_point(data = dG[dG$removed == "Removed", , drop = FALSE],
                                           ggplot2::aes(x = rank_i, y = r), inherit.aes = FALSE,
                                           shape = 25, size = 3, fill = COL_BAR_REM, colour = "white",
                                           position = ggplot2::position_nudge(y = 0.04 * max(dG$r)))

        # I Gate convergence
        dI <- data.frame(p = pmax(gt$p, 1e-300), eff = gt$effect, iter = seq_along(gt$p))
        pI <- ggplot2::ggplot(dI, ggplot2::aes(x = p, y = eff)) +
            ggplot2::annotate("rect", xmin = min(dI$p, gt$alpha) / 10, xmax = gt$alpha, ymin = gt$bv, ymax = Inf,
                              fill = COL_TEST, alpha = 0.08) +
            ggplot2::geom_hline(yintercept = gt$bv, colour = COL_A, linetype = "dashed", linewidth = 0.6) +
            ggplot2::geom_vline(xintercept = gt$alpha, colour = COL_A, linetype = "dashed", linewidth = 0.6) +
            ggplot2::geom_path(colour = COL_B, linewidth = 0.35, alpha = 0.6) +
            ggplot2::geom_point(ggplot2::aes(colour = iter), size = 1.4) +
            ggplot2::geom_point(data = dI[it, , drop = FALSE], shape = 21, size = 4, stroke = 0.9,
                                colour = "black", fill = NA) +
            ggplot2::scale_colour_gradientn(colours = grDevices::hcl.colors(256, "Plasma"),
                                            name = "Iteration") +
            ggplot2::scale_x_log10() +
            ggplot2::labs(title = "I Gate convergence",
                          subtitle = sprintf("%d iterations, %s removed. Ring: %s.\nRemoval: p < %.2f and dRMSE > BV (shaded).",
                                             length(gt$p), fmt_n(gt$n_removed), it_lab, gt$alpha),
                          x = "Wilcoxon p (log scale)", y = "dRMSE") +
            theme_cepip() + ggplot2::theme(legend.position = "bottom",
                                           legend.key.width = grid::unit(0.9, "cm"),
                                           legend.key.height = grid::unit(0.25, "cm"))
    } else {
        pG <- empty_panel("G Candidates ranked by PAVA residual", "No candidates")
        pI <- empty_panel("I Gate convergence", "No candidates")
    }

    # Q Class counts by IP side
    dQ <- as.data.frame(table(cls = factor(cl$cell_class, levels = CLASS_LEVELS),
                              side = factor(cl$ip_side_all, levels = c("Pre-IP", "Post-IP"))),
                        responseName = "n")
    pQ <- ggplot2::ggplot(dQ, ggplot2::aes(x = cls, y = n, fill = cls)) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::geom_text(ggplot2::aes(label = fmt_n(n)), vjust = -0.4, size = 3) +
        ggplot2::facet_wrap(~side) +
        ggplot2::scale_fill_manual(values = COL_CLASS, guide = "none") +
        ggplot2::scale_x_discrete(labels = c("Original + Refined" = "OriRef", "Recovered" = "Recovered",
                                             "Original only" = "Original\nonly", "Excluded" = "Excluded")) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
        ggplot2::labs(title = "Q Cells per class and IP side",
                      subtitle = sprintf("Module = %s cells (OriRef plus Recovered).",
                                         fmt_n(sum(cl$in_module))),
                      x = NULL, y = "Cells") +
        theme_cepip() + ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))

    grobs <- list(pA, pB, pC, pD, pE, pF, pG, pI, pQ)
    lay <- rbind(c(1, 2, 3, 4), c(5, 5, 6, 6), c(7, 8, 9, 9))
    hgt <- 13.5
    tests <- cmrf$tests
    if (!is.null(tests) && nrow(tests)) {
        trow <- function(pat) tests[grepl(pat, tests$test), , drop = FALSE]
        jt <- trow("^Jonckheere"); hd <- trow("^Hoeffding"); dc_ <- trow("^Distance")
        ptxt <- function(rw) if (!nrow(rw)) "not run" else
            sprintf("%.4g, p = %s (%s permutations)", rw$statistic[1], fmt_p(rw$p_value[1]),
                    fmt_n(rw$nperm[1]))
        dP <- cl[, c("cell_class", "resid_abs")]
        dP$cls <- factor(dP$cell_class, levels = CLASS_LEVELS)
        nP <- table(dP$cls)
        pP <- ggplot2::ggplot(dP, ggplot2::aes(x = cls, y = resid_abs, fill = cls)) +
            ggplot2::geom_jitter(ggplot2::aes(colour = cls), width = 0.18, size = 0.7, alpha = 0.15,
                                 show.legend = FALSE) +
            ggplot2::geom_boxplot(width = 0.45, colour = "black", linewidth = 0.5, outlier.shape = NA,
                                  alpha = 0.35) +
            ggplot2::geom_hline(yintercept = cmrf$band$delta_base, colour = COL_FLAT,
                                linetype = "dashed", linewidth = 0.6) +
            ggplot2::annotate("text", x = -Inf, y = cmrf$band$delta_base, label = "delta_base",
                              colour = COL_FLAT, hjust = -0.1, vjust = -0.5, size = 2.8) +
            ggplot2::scale_fill_manual(values = COL_CLASS, guide = "none") +
            ggplot2::scale_colour_manual(values = COL_CLASS, guide = "none") +
            ggplot2::scale_x_discrete(labels = sprintf("%s\nn=%s",
                                                       c("OriRef", "Recovered", "Original only", "Excluded"),
                                                       fmt_n(as.integer(nP)))) +
            ggplot2::labs(title = "P Spine residual by class",
                          subtitle = sprintf("Jonckheere-Terpstra JT = %s", ptxt(jt)),
                          x = NULL, y = "|v - g(u)|") +
            theme_cepip()
        rec <- cl[cl$cell_class == "Recovered", , drop = FALSE]
        pO <- if (nrow(rec) > 0L) {
            ggplot2::ggplot(rec, ggplot2::aes(x = u, y = v, colour = resid_abs)) +
                ggplot2::geom_point(size = 1.3, alpha = 0.8) +
                ggplot2::geom_vline(xintercept = cmrf$band$IP_u, colour = COL_IP, linetype = "dashed",
                                    linewidth = 0.5) +
                ggplot2::geom_hline(yintercept = cmrf$band$IP_v, colour = COL_IPV, linetype = "dashed",
                                    linewidth = 0.45) +
                ggplot2::scale_colour_gradientn(colours = grDevices::hcl.colors(256, "Inferno")[13:243],
                                                limits = c(0, NA), name = "|v - g(u)|") +
                ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
                ggplot2::labs(title = "O Recovered cells",
                              subtitle = sprintf("Hoeffding D = %s\ndCor = %s", ptxt(hd), ptxt(dc_)),
                              x = "u", y = "v") +
                theme_cepip() + theme_vtitle() + ggplot2::theme(aspect.ratio = 1)
        } else empty_panel("O Recovered cells", "No Recovered cells")
        grobs <- c(grobs, list(pP, pO))
        lay <- rbind(lay, c(10, 10, 11, NA))
        hgt <- 18
    }
    combine_panels(grobs, lay,
                   title %|||% sprintf("RSk-fCV gate, %s (%d folds x %d repeats)", it_lab, K, R),
                   width = 17, height = hgt)
}

# =============================================================================
# 9. Module pitcher plot in copula space (Figure 8A)
# =============================================================================

#' Pitcher plot of the CMRF classes in copula space (Figure 8A)
#'
#' u_SK = F(S_i,norm) (x) against v_SK = F(K_i) (y). Module cells are the
#' OriRef cells in their four groups and the Recovered cells (purple
#' diamonds); non-Module cells are the Excluded cells (gray) and the Original
#' only cells (tan triangles). The S_i,norm and K_i gates are drawn in copula
#' coordinates.
#'
#' @param cmrf A `cepip_cmrf` object, or a `cepip` object.
#' @param title Optional plot title.
#' @return A ggplot object.
#' @examples
#' \donttest{
#' d <- cepip_toy
#' # CMRF with a short gate of 5 repeats (the default is 40)
#' res <- cepip(d$x, d$y, donor = d$donor, id = d$id, cmrf = TRUE, repeats = 5,
#'              tests = FALSE, verbose = FALSE)
#' plot_module(res)
#' }
#' @export
plot_module <- function(cmrf, title = NULL) {
    if (inherits(cmrf, "cepip")) cmrf <- cmrf$cmrf
    check_class(cmrf, "cepip_cmrf", "cmrf")
    cl <- cmrf$cells
    lev_m <- c(names(COL_FOUR), "Recovered")
    lev_n <- c("Excluded", "Original only")
    ex <- cl[cl$cell_class == "Excluded", , drop = FALSE]
    oo <- cl[cl$cell_class == "Original only", , drop = FALSE]
    nm <- rbind(data.frame(u_SK = ex$u_SK, v_SK = ex$v_SK, nm_cls = "Excluded"),
                data.frame(u_SK = oo$u_SK, v_SK = oo$v_SK, nm_cls = "Original only"))
    nm$nm_cls <- factor(nm$nm_cls, levels = lev_n)
    ori <- cl[cl$cell_class == "Original + Refined", , drop = FALSE]
    ori <- data.frame(u_SK = ori$u_SK, v_SK = ori$v_SK, grp = factor(ori$group_orig, levels = lev_m))
    rc <- cl[cl$cell_class == "Recovered", , drop = FALSE]
    rc <- data.frame(u_SK = rc$u_SK, v_SK = rc$v_SK, grp = factor("Recovered", levels = lev_m))
    lab_m <- with_counts(lev_m, c(as.character(ori$grp), as.character(rc$grp)))
    lab_n <- with_counts(lev_n, as.character(nm$nm_cls))
    gs <- cmrf$gates_SK
    ggplot2::ggplot() +
        ggplot2::geom_point(data = nm, ggplot2::aes(x = u_SK, y = v_SK, alpha = nm_cls),
                            shape = ifelse(nm$nm_cls == "Excluded", 16, 25),
                            colour = ifelse(nm$nm_cls == "Excluded", COL_BG, "white"),
                            fill = COL_ORIONLY,
                            size = ifelse(nm$nm_cls == "Excluded", 0.6, 2.0),
                            stroke = ifelse(nm$nm_cls == "Excluded", 0.5, 0.15)) +
        ggplot2::geom_point(data = ori, ggplot2::aes(x = u_SK, y = v_SK, colour = grp, fill = grp,
                                                     shape = grp), size = 1.2, alpha = 0.55) +
        ggplot2::geom_point(data = rc, ggplot2::aes(x = u_SK, y = v_SK, colour = grp, fill = grp,
                                                    shape = grp), size = 2.0, alpha = 0.93, stroke = 0.15) +
        ggplot2::geom_vline(xintercept = gs[["u"]], colour = "#009E73", linewidth = 0.6) +
        ggplot2::geom_hline(yintercept = gs[["v"]], colour = "#009E73", linewidth = 0.6) +
        ggplot2::annotate("text", x = gs[["u"]], y = 0, label = "Si > 0.05", hjust = -0.05,
                          vjust = 0, size = 3, fontface = "bold") +
        ggplot2::annotate("text", x = 0, y = gs[["v"]], label = "Ki > 0.5", hjust = 0,
                          vjust = -0.4, size = 3, fontface = "bold") +
        ggplot2::scale_colour_manual(values = c(COL_FOUR, Recovered = "white"), limits = lev_m,
                                     labels = lab_m, name = "Module cells", drop = FALSE) +
        ggplot2::scale_fill_manual(values = c(COL_FOUR, Recovered = COL_RECOV), limits = lev_m,
                                   labels = lab_m, name = "Module cells", drop = FALSE) +
        ggplot2::scale_shape_manual(values = c(16, 16, 16, 16, 23), limits = lev_m,
                                    labels = lab_m, name = "Module cells", drop = FALSE) +
        ggplot2::scale_alpha_manual(values = c(Excluded = 0.3, "Original only" = 0.88), limits = lev_n,
                                    labels = lab_n, name = "Non-Module cells", drop = FALSE) +
        merged_guides(2.5) +
        ggplot2::guides(
            alpha = ggplot2::guide_legend(order = 2, override.aes = list(
                shape = c(16, 25), colour = c("gray70", "white"), fill = c(NA, COL_ORIONLY),
                size = c(2.2, 2.5), stroke = c(0.5, 0.3), alpha = 1))) +
        ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
        ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.25)) +
        ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.25)) +
        ggplot2::labs(title = title %|||% "Pitcher plot in copula space",
                      subtitle = sprintf("u_SK = F(Si), v_SK = F(Ki), %s cells. Module = %s cells.",
                                         fmt_n(nrow(cl)), fmt_n(sum(cl$in_module))),
                      x = "u_SK", y = "v_SK") +
        theme_cepip() + theme_vtitle() + ggplot2::theme(aspect.ratio = 1, legend.position = "right")
}

# =============================================================================
# 10. Validation1
# =============================================================================

#' Plot the Validation1 refit of the CMRF Module
#'
#' Left: the Module cells (x against y) in their refitted four groups
#' (Recovered cells as diamonds) with the refitted population curve and its
#' 95% band. Right: the original against the refitted four-group labels of the
#' OriRef cells, with the agreement in the title.
#'
#' @param validate A `cepip_validate` object, or a `cepip` object (its
#'   `validate` and `cmrf` are used).
#' @param cmrf Optional `cepip_cmrf` object, for the IPx and IPy lines.
#' @param title Optional overall title.
#' @return A gtable of class `cepip_figure` (printing it draws it; it can
#'   also be drawn with [grid::grid.draw()]). The attribute `cepip_size`
#'   gives a suitable width and height in inches.
#' @examples
#' \donttest{
#' d <- cepip_toy
#' # CMRF with a short gate of 5 repeats (the default is 40)
#' res <- cepip(d$x, d$y, donor = d$donor, id = d$id, cmrf = TRUE, repeats = 5,
#'              tests = FALSE, verbose = FALSE)
#' plot_validate(res)
#' }
#' @export
plot_validate <- function(validate, cmrf = NULL, title = NULL) {
    if (inherits(validate, "cepip")) {
        if (is.null(cmrf)) cmrf <- validate$cmrf
        validate <- validate$validate
    }
    check_class(validate, "cepip_validate", "validate")
    fit <- validate$fit; vc <- validate$cells
    d <- data.frame(x = fit$data$GOI, y = fit$data$Expression,
                    grp = factor(vc$group_ref, levels = names(COL_FOUR)),
                    rec = vc$cell_class == "Recovered")
    xs <- seq(min(d$x), max(d$x), length.out = 400L)
    pr <- predict_population(fit$model, xs, fit$data, se = TRUE)
    cv <- data.frame(x = xs, fit = as.numeric(pr$fit),
                     lo = as.numeric(pr$fit) - 1.96 * as.numeric(pr$se.fit),
                     hi = as.numeric(pr$fit) + 1.96 * as.numeric(pr$se.fit))
    pl <- ggplot2::ggplot() +
        ggplot2::geom_point(data = d[!d$rec, ], ggplot2::aes(x = x, y = y, colour = grp),
                            size = 0.9, alpha = 0.5) +
        ggplot2::geom_point(data = d[d$rec, ], ggplot2::aes(x = x, y = y, fill = grp), shape = 23,
                            colour = "white", stroke = 0.2, size = 1.8, alpha = 0.85) +
        ggplot2::geom_ribbon(data = cv, ggplot2::aes(x = x, ymin = lo, ymax = hi), fill = "gray20",
                             alpha = 0.2) +
        ggplot2::geom_line(data = cv, ggplot2::aes(x = x, y = fit), colour = "gray10", linewidth = 1)
    if (!is.null(cmrf)) {
        b <- if (inherits(cmrf, "cepip_cmrf")) cmrf$band else cmrf
        pl <- pl + ggplot2::geom_vline(xintercept = b$IP_x, colour = COL_IP, linetype = "dashed",
                                       linewidth = 0.6) +
            ggplot2::geom_hline(yintercept = b$IP_y, colour = COL_IPV, linetype = "dashed",
                                linewidth = 0.5)
    }
    pl <- pl + ggplot2::scale_colour_manual(values = COL_FOUR, limits = names(COL_FOUR), drop = FALSE,
                                            labels = with_counts(names(COL_FOUR), d$grp), name = NULL) +
        ggplot2::scale_fill_manual(values = COL_FOUR, limits = names(COL_FOUR), guide = "none") +
        ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2.2, alpha = 1))) +
        ggplot2::labs(title = "Refit on the Module cells",
                      subtitle = sprintf("%s Module cells, k = %d, DE = %.4f. Diamonds: Recovered (n=%s).",
                                         fmt_n(nrow(d)), fit$k, validate$de, fmt_n(sum(d$rec))),
                      x = "GOI", y = "Composite") +
        theme_cepip() + legend_inside(c(0.01, 0.99), c(0, 1))
    xt <- as.data.frame(validate$crosstab, responseName = "n")
    xt$Original <- factor(xt$Original, levels = rev(names(COL_FOUR)))
    xt$Refitted <- factor(xt$Refitted, levels = names(COL_FOUR))
    pr_ <- ggplot2::ggplot(xt, ggplot2::aes(x = Refitted, y = Original, fill = n)) +
        ggplot2::geom_tile(colour = "white", linewidth = 0.8) +
        ggplot2::geom_text(ggplot2::aes(label = fmt_n(n),
                                        colour = n > 0.6 * max(n))) +
        ggplot2::scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "gray10"), guide = "none") +
        ggplot2::scale_fill_gradientn(colours = rev(grDevices::hcl.colors(256, "Blues 3")), name = "Cells") +
        ggplot2::scale_x_discrete(labels = function(z) sub(" ", "\n", z)) +
        ggplot2::labs(title = sprintf("Label agreement %.1f percent", 100 * validate$agreement),
                      subtitle = sprintf("Four-group labels of the %s OriRef cells.", fmt_n(validate$n_oriref)),
                      x = "Refitted", y = "Original") +
        theme_cepip() + ggplot2::theme(panel.grid = ggplot2::element_blank(), aspect.ratio = 1)
    combine_panels(list(pl, pr_), matrix(1:2, nrow = 1), title %|||% "Validation1",
                   width = 13, height = 6)
}

# =============================================================================
# 11. Donors
# =============================================================================

#' Plot the per-donor Module Post-IP percentage
#'
#' With only the table of [cepip_donor()], a sorted bar chart of `pct_PostIP`
#' per donor. With `data` and `outcome`, a scatter of `pct_PostIP` against the
#' outcome with a linear fit.
#'
#' @param donor_table Output of [cepip_donor()], or a `cepip`, `cepip_cmrf` or
#'   `cepip_ip` object (then [cepip_donor()] is run on it).
#' @param outcome Optional column name of `data`.
#' @param data Optional donor table with a `donor` column and the outcome.
#' @param title Optional plot title.
#' @return A ggplot object.
#' @examples
#' d <- cepip_toy
#' res <- cepip(d$x, d$y, donor = d$donor, id = d$id, verbose = FALSE)
#' plot_donor(res)
#' @export
plot_donor <- function(donor_table, outcome = NULL, data = NULL, title = NULL) {
    if (inherits(donor_table, c("cepip", "cepip_cmrf", "cepip_ip")))
        donor_table <- cepip_donor(donor_table)
    if (!all(c("donor", "pct_PostIP") %in% names(donor_table)))
        stop("donor_table needs the columns donor and pct_PostIP.", call. = FALSE)
    if (is.null(outcome) || is.null(data)) {
        d <- donor_table[, c("donor", "pct_PostIP")]
        d$donor <- factor(d$donor, levels = d$donor[order(d$pct_PostIP)])
        return(ggplot2::ggplot(d, ggplot2::aes(x = pct_PostIP, y = donor)) +
            ggplot2::geom_col(fill = "#009E73", width = 0.75) +
            ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
            ggplot2::labs(title = title %|||% "Module Post-IP cells per donor",
                          subtitle = sprintf("%d donors. Median %.1f percent.", nrow(d),
                                             stats::median(d$pct_PostIP)),
                          x = "Module Post-IP (percent of the donor's cells)", y = NULL) +
            theme_cepip() + ggplot2::theme(
                panel.grid.major.y = ggplot2::element_blank(),
                axis.text.y = ggplot2::element_text(size = if (nrow(d) > 40) 5 else 7)))
    }
    if (!all(c("donor", outcome) %in% names(data)))
        stop("data needs the columns donor and ", outcome, ".", call. = FALSE)
    d <- merge(donor_table[, c("donor", "pct_PostIP")],
               data.frame(donor = data$donor, outcome_value = as.numeric(data[[outcome]])),
               by = "donor")
    d <- d[is.finite(d$pct_PostIP) & is.finite(d$outcome_value), ]
    if (nrow(d) < 3L) stop("Fewer than 3 donors with both values.", call. = FALSE)
    ct <- stats::cor.test(d$pct_PostIP, d$outcome_value)
    ggplot2::ggplot(d, ggplot2::aes(x = pct_PostIP, y = outcome_value)) +
        ggplot2::geom_smooth(method = "lm", formula = y ~ x, colour = COL_IP, fill = COL_IP,
                             alpha = 0.12, linewidth = 0.8) +
        ggplot2::geom_point(size = 2, shape = 21, fill = "#009E73", colour = "white", stroke = 0.3) +
        ggplot2::labs(title = title %|||% sprintf("Module Post-IP and %s", outcome),
                      subtitle = sprintf("%d donors. Pearson r = %.3f, p = %s.", nrow(d),
                                         unname(ct$estimate), fmt_p(ct$p.value)),
                      x = "Module Post-IP (percent of the donor's cells)", y = outcome) +
        theme_cepip()
}

# =============================================================================
# 12. DEG volcano
# =============================================================================

#' Volcano plot of one differential expression contrast
#'
#' avg_log2FC (x) against -log10 of the BH-adjusted p (y) for one contrast of
#' [cepip_deg()]. Significant genes are coloured by direction and the top 10
#' are labelled; dashed lines mark the `lfc` and `fdr` thresholds.
#'
#' @param deg A `cepip_deg` object.
#' @param contrast Name of the contrast. By default the first one that was
#'   tested.
#' @param title Optional plot title.
#' @return A ggplot object.
#' @examples
#' d <- cepip_toy
#' res <- cepip(d$x, d$y, donor = d$donor, id = d$id, verbose = FALSE)
#' expr <- t(as.matrix(d[, c("x", paste0("P", 1:8))]))
#' colnames(expr) <- d$id
#' deg <- cepip_deg(res, expr, contrast = "PostIP_vs_PreIP")
#' plot_deg(deg)
#' @export
plot_deg <- function(deg, contrast = NULL, title = NULL) {
    check_class(deg, "cepip_deg", "deg")
    ok <- names(deg)[!vapply(deg, is.null, logical(1))]
    if (!length(ok)) stop("No contrast of deg was tested.", call. = FALSE)
    if (is.null(contrast)) contrast <- ok[1]
    if (!contrast %in% ok) stop("Contrast ", contrast, " is missing or was not tested.", call. = FALSE)
    st <- attr(deg, "settings")
    lfc <- st$lfc %|||% 0.25; fdr <- st$fdr %|||% 0.05
    r <- deg[[contrast]]
    r <- r[!is.na(r$p_val_adj) & is.finite(r$avg_log2FC), , drop = FALSE]
    pmin_ <- suppressWarnings(min(r$p_val_adj[r$p_val_adj > 0]))
    if (!is.finite(pmin_)) pmin_ <- 1e-300
    parts <- strsplit(contrast, "_vs_", fixed = TRUE)[[1]]
    up1 <- sprintf("Up in %s", parts[1]); up2 <- sprintf("Up in %s", parts[2])
    d <- data.frame(gene = r$gene, avg_log2FC = r$avg_log2FC,
                    nlp = -log10(pmax(r$p_val_adj, pmin_ / 10)),
                    dir = ifelse(!r$significant, "Not significant",
                                 ifelse(r$avg_log2FC > 0, up1, up2)))
    d$dir <- factor(d$dir, levels = c(up1, up2, "Not significant"))
    top <- d[d$dir != "Not significant", , drop = FALSE]
    top <- top[order(-top$nlp, -abs(top$avg_log2FC)), , drop = FALSE][seq_len(min(10L, nrow(top))), ]
    # Labels of genes with similar values are staggered vertically.
    top <- top[order(top$avg_log2FC), , drop = FALSE]
    top$vj <- rep(c(-0.7, 1.7, -2.1), length.out = nrow(top))
    cols <- stats::setNames(c("#D94F4F", "#1B3A8C", COL_GRAY), levels(d$dir))
    ggplot2::ggplot(d, ggplot2::aes(x = avg_log2FC, y = nlp, colour = dir)) +
        ggplot2::geom_vline(xintercept = c(-lfc, lfc), colour = "gray40", linetype = "dashed",
                            linewidth = 0.4) +
        ggplot2::geom_hline(yintercept = -log10(fdr), colour = "gray40", linetype = "dashed",
                            linewidth = 0.4) +
        ggplot2::geom_point(data = d[d$dir == "Not significant", ], size = 0.9, alpha = 0.5) +
        ggplot2::geom_point(data = d[d$dir != "Not significant", ], size = 1.4, alpha = 0.85) +
        ggplot2::geom_text(data = top, ggplot2::aes(label = gene), size = 2.8, vjust = top$vj,
                           check_overlap = TRUE, show.legend = FALSE) +
        ggplot2::scale_colour_manual(values = cols, labels = with_counts(levels(d$dir), d$dir),
                                     name = NULL, drop = FALSE) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.1))) +
        ggplot2::labs(title = title %|||% sprintf("Differential expression, %s", gsub("_", " ", contrast)),
                      subtitle = sprintf("%s genes tested (%s unit). Significant: BH p < %.2f and |log2FC| >= %.2f.",
                                         fmt_n(nrow(d)), st$unit %|||% "donor", fdr, lfc),
                      x = "Average log2 fold change", y = "-log10 adjusted p") +
        theme_cepip() + ggplot2::theme(legend.position = "bottom")
}

# =============================================================================
# 13. Association forest plot
# =============================================================================

#' Forest plot of donor-level associations
#'
#' The standardized beta of the predictor with its 95% CI for each outcome of
#' [cepip_associate()], in the order given. Point size shows delta R2; filled
#' points have FDR < 0.05. delta R2 and the FDR are written beside each row.
#'
#' @param assoc Output of [cepip_associate()].
#' @param title Optional plot title.
#' @return A ggplot object.
#' @examples
#' d <- cepip_toy
#' res <- cepip(d$x, d$y, donor = d$donor, id = d$id, verbose = FALSE)
#' dt <- cepip_donor(res)
#' set.seed(1)
#' dt$score <- 0.05 * dt$pct_PostIP + rnorm(nrow(dt))
#' a <- cepip_associate(dt, "pct_PostIP", outcomes = "score", min_donors = 10)
#' plot_associate(a)
#' @export
plot_associate <- function(assoc, title = NULL) {
    need <- c("outcome", "beta_std", "ci_low", "ci_high", "delta_R2", "FDR")
    if (!all(need %in% names(assoc)))
        stop("assoc needs the columns ", paste(need, collapse = ", "), ".", call. = FALSE)
    d <- assoc
    d$outcome <- factor(d$outcome, levels = rev(unique(d$outcome)))
    d <- d[is.finite(d$beta_std), , drop = FALSE]
    if (!nrow(d)) stop("No outcome has an estimate.", call. = FALSE)
    d$sig <- factor(ifelse(!is.na(d$FDR) & d$FDR < 0.05, "FDR < 0.05", "FDR >= 0.05"),
                    levels = c("FDR < 0.05", "FDR >= 0.05"))
    d$txt <- sprintf("dR2 %.3f, FDR %s", d$delta_R2, vapply(d$FDR, fmt_p, character(1)))
    xr <- range(c(d$ci_low, d$ci_high, 0), na.rm = TRUE)
    xt <- xr[2] + 0.05 * diff(xr)
    ggplot2::ggplot(d, ggplot2::aes(y = outcome)) +
        ggplot2::geom_vline(xintercept = 0, colour = "gray40", linewidth = 0.4) +
        ggplot2::geom_segment(ggplot2::aes(x = ci_low, xend = ci_high, yend = outcome),
                              colour = "gray20", linewidth = 0.6) +
        ggplot2::geom_point(ggplot2::aes(x = beta_std, size = delta_R2, fill = sig), shape = 21,
                            colour = "black", stroke = 0.5) +
        ggplot2::geom_text(ggplot2::aes(x = xt, label = txt), hjust = 0, size = 2.8, colour = "gray25") +
        ggplot2::scale_fill_manual(values = c("FDR < 0.05" = "#D94F4F", "FDR >= 0.05" = "white"),
                                   name = NULL, drop = FALSE) +
        ggplot2::scale_size_continuous(range = c(1.8, 5), name = "delta R2") +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.45))) +
        ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 3))) +
        ggplot2::labs(title = title %|||% sprintf("Association with %s", d$predictor[1] %|||% "the predictor"),
                      subtitle = "Standardized beta with 95% CI; covariate-adjusted OLS on donors.",
                      x = "Standardized beta", y = NULL) +
        theme_cepip() + ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}
