# readme_figure.R
# Draws the six-panel example figure of the cepip README from the SEA-AD example data:
# PVM_3 NeoC cells, MT-CO3 against the mean of its eight partner genes, with the
# settings of the published run (reproduce_seaad.R).
#
# The colours, point sizes, and line styles follow the SEA-AD scripts that drew the
# paper figures (named at the start of the plotting section), with stronger grey
# cells in all panels and larger points in panels A and F.
#
#   A  SKCS-core cells (four groups) and non-SKCS-core cells (grey), with the SKCS-core GAM
#   B  Residual GAM R(x): binned mean residuals (bars) and the fitted R(x), with the IP
#   C  CMRF band in copula space with the candidate cells of the RSk-fCV gate
#   D  Candidate cells recovered or excluded by the RSk-fCV gate
#   E  Final Module, with the non-Module cells in grey
#   F  Pitcher plot in copula space
#
# Output folder (first argument, default readme_fig):
#   cepip_example_figure.png (all six panels), cepip_example_AB.png, cepip_example_CD.png,
#   cepip_example_EF.png (two panels each, for the worked example of the README),
#   panel_A.png ... panel_F.png, readme_figure_log.txt
#   The log holds IPx and IPRS for checking against the published values
#   (PVM_3 NeoC: IPx 2.4953, IPRS 0.791) and the cell counts of each panel.
#
# Usage inside the container, from /data/Package:
#   R_LIBS=/data/Package/Rlib Rscript readme_figure.R readme_fig 8 cepip_example
#   arguments: output folder, workers, folder with the example files (optional;
#   without it the files are downloaded with cepip_example())

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[1] else "readme_fig"
workers <- if (length(args) >= 2) as.integer(args[2]) else 4L
src <- if (length(args) >= 3) args[3] else NA
dir.create(out, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(out, "readme_figure_log.txt")
zz <- file(log_file, open = "wt")
sink(zz, split = TRUE); sink(zz, type = "message")

suppressPackageStartupMessages({ library(cepip); library(ggplot2) })
cat("cepip ", as.character(packageVersion("cepip")), " | R ", R.version$major, ".", R.version$minor, "\n", sep = "")

if (!is.na(src) && file.exists(file.path(src, "pvm3_cells.csv"))) {
    pvm3 <- utils::read.csv(file.path(src, "pvm3_cells.csv"), check.names = FALSE)
    cat("Data: ", file.path(src, "pvm3_cells.csv"), "\n", sep = "")
} else {
    pvm3 <- cepip_example("pvm3")
    if (is.null(pvm3)) stop("The example data could not be loaded.")
    cat("Data: cepip_example(\"pvm3\")\n")
}

mt8 <- c("MT-ATP6", "MT-CO2", "MT-ND3", "MT-ND4", "MT-CYB", "MT-ND1", "MT-ND2", "MT-CO1")
d <- pvm3[pvm3$region_set == "NeoC", ]
y <- cepip_composite(d, mt8)
cat("PVM_3 NeoC:", nrow(d), "cells,", length(unique(d$donor)), "donors\n")

set.seed(1)
res <- cepip(d[["MT-CO3"]], y, donor = d$donor, id = d$barcode, cmrf = TRUE, ip_k = 10,
             y_spine = d$composite_part8, workers = workers)
print(res)

# ---------- styles copied from the SEA-AD scripts ----------
# A: Part10b_CP10K_SKCS_TREP_code.r (p_curve_4c); B: Part10b_CP10K_IPRS_code.r (create_rx);
# C to E: Part10f_CP10K_Band_Formation_code.r (BF09, BF10, BF11); F: Part10e_CP10K_Validation1_code.r (p5_sk).
FOUR_COLS <- c("Pre-IP TREP" = "#E69F00", "Pre-IP non-TREP" = "#1B3A8C",
               "Post-IP TREP" = "#009E73", "Post-IP non-TREP" = "#D94F4F")
GRAY_COL <- "gray65"; COL_SPINE <- "#00695C"; COL_ADAPT <- "#D55E00"; COL_ADAPT_BORDER <- "#B84A00"
COL_FLAT <- "#0072B2"; BLUE_FILL <- "#2171B5"; GAP_FILL <- "#2C2C2C"; TAN_FILL <- "#8B6914"
RECOV_FILL <- "#7C3AED"; S2R_COL <- "#404040"
DRAW_ORD <- c("Pre-IP non-TREP" = 1, "Post-IP non-TREP" = 2, "Pre-IP TREP" = 3, "Post-IP TREP" = 4)

BS <- 14   # base font size of each panel in the combined figure
XLAB <- expression(italic("MT-CO3"))
YLAB <- "Mean of the eight partner genes"
fmt <- function(n) formatC(n, format = "d", big.mark = ",")
cnt <- function(levels, values) setNames(sprintf("%s (n=%s)", levels,
                                         fmt(as.integer(table(factor(values, levels = levels))))), levels)
theme_fig <- function()
    theme_minimal(base_size = BS) +
    theme(plot.title = element_text(face = "bold", size = BS + 1),
          plot.subtitle = element_text(size = BS - 2, colour = "gray30"),
          axis.title = element_text(face = "bold", size = BS - 1),
          axis.text = element_text(size = BS - 2),
          panel.grid.major = element_line(colour = "gray90"),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.35),
          plot.background = element_rect(fill = "white", colour = NA),
          panel.background = element_rect(fill = "white", colour = NA),
          aspect.ratio = 1,
          legend.title = element_blank(), legend.text = element_text(size = BS - 3),
          legend.key.size = unit(0.55, "cm"),
          legend.background = element_rect(fill = alpha("white", 0.85), colour = NA))
inside <- function(pos, just) theme(legend.position = pos, legend.justification = just)
uv <- function(xlab = "u", ylab = "v")
    list(coord_cartesian(xlim = c(0, 1), ylim = c(-0.04, 1.04)),
         scale_x_continuous(breaks = seq(0, 1, 0.25)),
         scale_y_continuous(breaks = seq(0, 1, 0.25), expand = expansion(mult = 0.01)),
         labs(x = xlab, y = ylab))

cl <- res$cmrf$cells
b <- res$cmrf$band
gr <- getFromNamespace("band_grid", "cepip")(b)
core <- cl$is_core
four <- cl$group_orig %in% names(FOUR_COLS)
cl$ord <- ifelse(four, DRAW_ORD[cl$group_orig], 0)

# ---------- A: SKCS-core cells in expression space (p_curve_4c, no IP line) ----------
sk <- res$skcs$fit
xs <- seq(min(sk$data$GOI), max(sk$data$GOI), length.out = 1000)
pr <- getFromNamespace("predict_population", "cepip")(sk$model, xs, sk$data, se = TRUE)
cvA <- data.frame(x = xs, fit = as.numeric(pr$fit), se = as.numeric(pr$se.fit))
levA <- c(names(FOUR_COLS), "Non-SKCS-core")
gA <- ifelse(core & four, cl$group_orig, "Non-SKCS-core")
dA <- data.frame(x = cl$x, y = cl$y, grp = factor(gA, levels = levA), ord = cl$ord)
dAc <- dA[dA$grp != "Non-SKCS-core", ]; dAc <- dAc[order(dAc$ord), ]
pA <- ggplot() +
    geom_point(data = dA[dA$grp == "Non-SKCS-core", ], aes(x, y, colour = grp), alpha = 0.45, size = 0.8) +
    geom_point(data = dAc, aes(x, y, colour = grp), alpha = 0.55, size = 1.2,
               position = position_jitter(width = 0.01, height = 0, seed = 123)) +
    geom_line(data = cvA, aes(x, fit), colour = "#4B0082", linewidth = 1.3) +
    geom_ribbon(data = cvA, aes(x = x, ymin = fit - 1.96 * se, ymax = fit + 1.96 * se), fill = "#4B0082", alpha = 0.12) +
    scale_colour_manual(values = c(FOUR_COLS, "Non-SKCS-core" = "#A9A9A9"), limits = levA, labels = cnt(levA, gA), drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(size = 2.4, alpha = 1))) +
    labs(title = "A   SKCS-core cells", x = XLAB, y = YLAB,
         subtitle = sprintf("%s of %s cells pass the Spearman and Kendall gates.", fmt(sum(core)), fmt(nrow(cl)))) +
    theme_fig() + inside(c(0.01, 0.99), c(0, 1))

# ---------- B: residual GAM R(x) (create_rx) ----------
rd <- res$ip$residual_data          # top-EP (purple) cells
brk <- seq(min(rd$GOI), max(rd$GOI), length.out = 101)
bid <- cut(rd$GOI, breaks = brk, include.lowest = TRUE, labels = FALSE)
ctr <- (brk[-length(brk)] + brk[-1]) / 2
bars <- data.frame(goi_center = ctr, mean_residual = vapply(seq_along(ctr), function(i) {
    idx <- which(bid == i); if (length(idx)) mean(rd$residual[idx]) else NA_real_ }, numeric(1)))
bars <- bars[!is.na(bars$mean_residual), ]
bar_w <- if (nrow(bars) > 1) min(diff(bars$goi_center)) * 0.8 else (brk[2] - brk[1]) * 0.8
g2 <- res$ip$grid
lc <- "#2E8B57"
pB <- ggplot() +
    annotate("rect", xmin = res$ip$ip_ci[1], xmax = res$ip$ip_ci[2], ymin = -Inf, ymax = Inf, fill = lc, alpha = 0.08) +
    geom_bar(data = bars, aes(x = goi_center, y = mean_residual), stat = "identity",
             fill = "#4B0082", colour = "#3D0065", alpha = 0.40, width = bar_w, linewidth = 0.4) +
    geom_ribbon(data = g2, aes(x = x, ymin = r_lower, ymax = r_upper), fill = "#FFCC99", alpha = 0.35) +
    geom_line(data = g2, aes(x = x, y = r_hat), colour = "#FFCC99", linewidth = 1.5) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "gray40", linewidth = 0.5) +
    geom_vline(xintercept = res$ip$ip_ci, linetype = "dotdash", colour = lc, linewidth = 0.7) +
    geom_vline(xintercept = res$ip$ip_x, linetype = "dashed", colour = lc, linewidth = 1.0) +
    labs(title = "B   Residual GAM R(x) and the inflection point", x = XLAB,
         y = "Expected signed residual, R(x)",
         subtitle = bquote(IP[x] * " = " * .(sprintf("%.3f", res$ip$ip_x)) * " (95% CI " *
                           .(sprintf("%.3f", res$ip$ip_ci[1])) * " to " * .(sprintf("%.3f", res$ip$ip_ci[2])) *
                           "), IPRS = " * .(sprintf("%.3f", res$ip$iprs$IPRS)) * ".")) +
    theme_fig()

# Shared copula-space layers (BF09 to BF11)
ribbons <- list(
    geom_ribbon(data = gr, aes(x = u, ymin = f_lo, ymax = f_up), fill = COL_FLAT, alpha = 0.18, inherit.aes = FALSE),
    geom_ribbon(data = gr, aes(x = u, ymin = a_lo, ymax = a_up), fill = COL_ADAPT, alpha = 0.18, inherit.aes = FALSE))
edges <- list(
    geom_line(data = gr, aes(x = u, y = a_up), colour = COL_ADAPT, linewidth = 0.55, inherit.aes = FALSE),
    geom_line(data = gr, aes(x = u, y = a_lo), colour = COL_ADAPT, linewidth = 0.55, inherit.aes = FALSE),
    geom_line(data = gr, aes(x = u, y = f_up), colour = COL_FLAT, linewidth = 0.45, linetype = "dashed", inherit.aes = FALSE),
    geom_line(data = gr, aes(x = u, y = f_lo), colour = COL_FLAT, linewidth = 0.45, linetype = "dashed", inherit.aes = FALSE))
spine <- geom_line(data = gr, aes(x = u, y = g), colour = COL_SPINE, linewidth = 1.1, inherit.aes = FALSE)
base_pts <- function(d) {
    gry <- d[!(d$is_core & d$group_orig %in% names(FOUR_COLS)), ]
    col <- d[d$is_core & d$group_orig %in% names(FOUR_COLS), ]; col <- col[order(col$ord), ]
    list(geom_point(data = gry, aes(u, v), colour = GRAY_COL, size = 0.6, alpha = 0.35),
         geom_point(data = col, aes(u, v), colour = FOUR_COLS[col$group_orig], size = 1.25, alpha = 0.55))
}

# ---------- C: candidates of the RSk-fCV gate (BF09) ----------
flat_cand <- !core & cl$pass_stage1 & cl$resid_abs <= b$delta_base
gap_cand <- !core & cl$pass_stage1 & cl$resid_abs > b$delta_base
oo <- cl$cell_class == "Original only"
levC <- c("flat", "gap", "oo")
dC <- rbind(data.frame(u = cl$u[flat_cand], v = cl$v[flat_cand], k = "flat"),
            data.frame(u = cl$u[gap_cand], v = cl$v[gap_cand], k = "gap"),
            data.frame(u = cl$u[oo], v = cl$v[oo], k = "oo"))
dC$k <- factor(dC$k, levels = levC)
labC <- c(flat = bquote("Stage 1 non-SKCS in " * delta[base] * " (n=" * .(fmt(sum(flat_cand))) * ")"),
          gap = bquote("Stage 1 non-SKCS in " * delta[adapt] * " only (n=" * .(fmt(sum(gap_cand))) * ")"),
          oo = bquote("Original only (n=" * .(fmt(sum(oo))) * ")"))
pC <- ggplot() + ribbons + base_pts(cl[!gap_cand, ]) + spine +
    geom_point(data = dC[dC$k != "oo", ], aes(u, v, fill = k, shape = k), size = 3.0, colour = "white", stroke = 0.5, alpha = 0.88) +
    geom_point(data = dC[dC$k == "oo", ], aes(u, v, fill = k, shape = k), size = 3.0, colour = "white", stroke = 0.4, alpha = 0.85) +
    edges +
    scale_fill_manual(values = c(flat = BLUE_FILL, gap = GAP_FILL, oo = TAN_FILL), limits = levC, labels = labC) +
    scale_shape_manual(values = c(flat = 21, gap = 21, oo = 25), limits = levC, labels = labC) +
    guides(fill = guide_legend(override.aes = list(size = 3.4, alpha = 1))) +
    uv() + labs(title = "C   Candidate cells for the RSk-fCV gate",
                subtitle = "Non-SKCS-core cells inside the CMRF bands.") +
    theme_fig() + inside(c(0.99, 0.01), c(1, 0))

# ---------- D: cells recovered or excluded by the gate (BF10) ----------
rec <- cl$cell_class == "Recovered"
rem <- !core & cl$pass_stage1 & cl$cell_class == "Excluded"
levD <- c("Recovered", "Excluded", "Original only")
dD <- rbind(data.frame(u = cl$u[rem], v = cl$v[rem], k = "Excluded"),
            data.frame(u = cl$u[oo], v = cl$v[oo], k = "Original only"),
            data.frame(u = cl$u[rec], v = cl$v[rec], k = "Recovered"))
dD$k <- factor(dD$k, levels = levD)
labD <- cnt(levD, dD$k)
pD <- ggplot() + ribbons + base_pts(cl) + spine +
    geom_point(data = dD[dD$k == "Excluded", ], aes(u, v), shape = 4, size = 2.8, colour = S2R_COL, stroke = 0.7, alpha = 0.75) +
    geom_point(data = dD[dD$k == "Original only", ], aes(u, v), shape = 25, size = 3.0, fill = TAN_FILL, colour = "white", stroke = 0.4, alpha = 0.85) +
    geom_point(data = dD[dD$k == "Recovered", ], aes(u, v), shape = 23, size = 3.0, fill = RECOV_FILL, colour = "white", stroke = 0.5, alpha = 0.88) +
    edges +
    geom_point(data = data.frame(u = -1, v = -1, k = factor(levD, levels = levD)), aes(u, v, shape = k, fill = k, colour = k)) +
    scale_shape_manual(values = c(23, 4, 25), limits = levD, labels = labD) +
    scale_fill_manual(values = c(RECOV_FILL, NA, TAN_FILL), limits = levD, labels = labD) +
    scale_colour_manual(values = c("white", S2R_COL, "white"), limits = levD, labels = labD) +
    guides(fill = "none", shape = "none",
           colour = guide_legend(override.aes = list(shape = c(23, 4, 25), fill = c(RECOV_FILL, NA, TAN_FILL),
                                                     colour = c("white", S2R_COL, "white"),
                                                     size = c(3.4, 3, 3.4), stroke = c(0.5, 0.9, 0.4), alpha = 1))) +
    uv() + labs(title = "D   Cells recovered or excluded by the gate",
                subtitle = "Candidates that passed (Recovered) or failed (Excluded) the RSk-fCV gate.") +
    theme_fig() + inside(c(0.99, 0.01), c(1, 0))

# ---------- E: final Module (BF11) ----------
lvE <- c(names(FOUR_COLS), "Recovered", "Non-Module cells")
gE <- ifelse(rec, "Recovered", ifelse(cl$cell_class == "Original + Refined" & four, cl$group_orig, "Non-Module cells"))
dE <- data.frame(u = cl$u, v = cl$v, grp = factor(gE, levels = lvE), ord = cl$ord)
dEc <- dE[dE$grp %in% names(FOUR_COLS), ]; dEc <- dEc[order(dEc$ord), ]
labE <- cnt(lvE, gE)
pE <- ggplot() +
    geom_ribbon(data = gr, aes(x = u, ymin = a_lo, ymax = a_up), fill = "gray70", alpha = 0.25, inherit.aes = FALSE) +
    geom_point(data = dE[dE$grp == "Non-Module cells", ], aes(u, v), colour = GRAY_COL, size = 0.6, alpha = 0.35) +
    geom_point(data = dEc, aes(u, v), colour = FOUR_COLS[as.character(dEc$grp)], size = 1.25, alpha = 0.55) +
    spine +
    geom_point(data = dE[dE$grp == "Recovered", ], aes(u, v), shape = 23, size = 3.0, fill = RECOV_FILL, colour = "white", stroke = 0.5, alpha = 0.88) +
    geom_line(data = gr, aes(x = u, y = a_up), colour = COL_ADAPT_BORDER, linewidth = 0.75, inherit.aes = FALSE) +
    geom_line(data = gr, aes(x = u, y = a_lo), colour = COL_ADAPT_BORDER, linewidth = 0.75, inherit.aes = FALSE) +
    geom_point(data = data.frame(u = -1, v = -1, k = factor(lvE, levels = lvE)), aes(u, v, shape = k, fill = k, colour = k)) +
    scale_shape_manual(values = c(16, 16, 16, 16, 23, 16), limits = lvE, labels = labE) +
    scale_colour_manual(values = c(FOUR_COLS, "white", GRAY_COL), limits = lvE, labels = labE) +
    scale_fill_manual(values = c(FOUR_COLS, RECOV_FILL, GRAY_COL), limits = lvE, labels = labE) +
    guides(fill = "none", shape = "none",
           colour = guide_legend(override.aes = list(shape = c(16, 16, 16, 16, 23, 16),
                                                     fill = c(FOUR_COLS, RECOV_FILL, GRAY_COL),
                                                     colour = c(FOUR_COLS, "white", GRAY_COL),
                                                     size = c(rep(2.6, 4), 3.4, 2.6), stroke = 0.5, alpha = 1))) +
    uv() + labs(title = "E   Final Module",
                subtitle = sprintf("%s Module cells; non-Module cells in grey.", fmt(sum(cl$in_module)))) +
    theme_fig() + inside(c(0.99, 0.01), c(1, 0))

# ---------- F: pitcher plot in copula space (p5_sk) ----------
gs <- res$cmrf$gates_SK
ori <- cl[cl$cell_class == "Original + Refined" & four, ]; ori <- ori[order(ori$ord), ]
lvF <- c(names(FOUR_COLS), "Recovered", "Excluded", "Original only")
gF <- ifelse(cl$cell_class == "Original + Refined" & four, cl$group_orig, cl$cell_class)
labF <- cnt(lvF, gF)
pF <- ggplot() +
    geom_point(data = cl[cl$cell_class == "Excluded", ], aes(u_SK, v_SK), colour = "gray60", alpha = 0.6, size = 1.0) +
    geom_point(data = cl[oo, ], aes(u_SK, v_SK), colour = "white", fill = TAN_FILL, shape = 25, size = 3.0, alpha = 0.95, stroke = 0.2) +
    geom_point(data = ori, aes(u_SK, v_SK), colour = FOUR_COLS[ori$group_orig], size = 1.9, alpha = 0.75, shape = 16) +
    geom_point(data = cl[rec, ], aes(u_SK, v_SK), colour = "white", fill = RECOV_FILL, shape = 23, size = 3.0, alpha = 1, stroke = 0.2) +
    geom_vline(xintercept = gs[["u"]], colour = "#009E73", linewidth = 0.6) +
    geom_hline(yintercept = gs[["v"]], colour = "#009E73", linewidth = 0.6) +
    geom_point(data = data.frame(u = -1, v = -1, k = factor(lvF, levels = lvF)), aes(u, v, shape = k, fill = k, colour = k)) +
    scale_shape_manual(values = c(16, 16, 16, 16, 23, 16, 25), limits = lvF, labels = labF) +
    scale_colour_manual(values = c(FOUR_COLS, "white", "gray70", "white"), limits = lvF, labels = labF) +
    scale_fill_manual(values = c(FOUR_COLS, RECOV_FILL, "gray70", TAN_FILL), limits = lvF, labels = labF) +
    guides(fill = "none", shape = "none",
           colour = guide_legend(override.aes = list(shape = c(16, 16, 16, 16, 23, 16, 25),
                                                     fill = c(FOUR_COLS, RECOV_FILL, "gray60", TAN_FILL),
                                                     colour = c(FOUR_COLS, "white", "gray60", "white"),
                                                     size = c(rep(2.6, 4), 3.2, 2.6, 3.2), stroke = 0.3, alpha = 1))) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    scale_x_continuous(breaks = seq(0, 1, 0.25)) + scale_y_continuous(breaks = seq(0, 1, 0.25)) +
    labs(title = "F   Pitcher plot in copula space",
         subtitle = "Module cells on the Spearman and Kendall cell-level scores.",
         x = expression(u[SK]), y = expression(v[SK])) +
    theme_fig() + inside(c(0.01, 0.99), c(0, 1))

cat("\nPanel counts:\n")
cat("  A: SKCS-core", sum(core), "| non-SKCS-core", sum(!core), "\n")
cat("  C: candidates in delta_base", sum(flat_cand), "| delta_adapt only", sum(gap_cand), "| original only", sum(oo), "\n")
cat("  D: recovered", sum(rec), "| excluded by the gate", sum(rem), "\n")
cat("  E: Module", sum(cl$in_module), "| non-Module", sum(!cl$in_module), "\n")

panels <- list(A = pA, B = pB, C = pC, D = pD, E = pE, F = pF)
for (n in names(panels))
    ggsave(file.path(out, paste0("panel_", n, ".png")), panels[[n]], width = 7.5, height = 7.5, dpi = 200, bg = "white")
g <- gridExtra::arrangeGrob(grobs = panels, ncol = 2)
ggsave(file.path(out, "cepip_example_figure.png"), g, width = 15, height = 22.5, dpi = 130, bg = "white")
pairs <- list(AB = c("A", "B"), CD = c("C", "D"), EF = c("E", "F"))
for (n in names(pairs))
    ggsave(file.path(out, paste0("cepip_example_", n, ".png")),
           gridExtra::arrangeGrob(grobs = panels[pairs[[n]]], ncol = 2),
           width = 15, height = 7.5, dpi = 130, bg = "white")
cat("\nSaved:", file.path(out, "cepip_example_figure.png"), "and",
    paste0("cepip_example_", names(pairs), ".png", collapse = ", "), "\n")

sink(type = "message"); sink(); close(zz)
