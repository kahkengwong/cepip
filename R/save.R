#' Save the tables and plots of one arm
#'
#' Writes all tables and plots of one arm into one folder, with no
#' subfolders. Several region sets (or other subsets) of the same arm go into
#' the same folder, each with its name in the file and sheet names.
#'
#' * `00_settings.txt`: settings, package versions and the summary of each
#'   result. `00_log.txt`: the text passed in `log` (for example the console
#'   output of the run, captured with [utils::capture.output()] or [sink()]).
#' * `<arm>_cepip_tables.xlsx`: one sheet per table, numbered in pipeline order.
#'   When the cell-level tables are large (more than `max_cell_rows` rows in
#'   total), they go into a second workbook, `<arm>_cepip_cells.xlsx`.
#' * Plots as .png and .pdf, numbered in pipeline order.
#'
#' @param results A `cepip` object, or a named list of `cepip` objects (one per
#'   region set, for example `list(NeoC = a, NonNeoC = b)`).
#' @param dir Output folder (created if needed).
#' @param arm Arm name used in the workbook names. Defaults to the folder name.
#' @param extras Optional named list (same names as `results`) of lists with
#'   any of `screen` ([cepip_screen()]), `donor` ([cepip_donor()]), `deg`
#'   ([cepip_deg()]), `associate` ([cepip_associate()]) and `outcomes` (a donor
#'   table for [plot_donor()]) for that region set.
#' @param plots Write the plots.
#' @param formats Plot file formats.
#' @param log Optional character vector written to `00_log.txt`.
#' @param max_cell_rows Row limit for keeping cell-level sheets in the main workbook.
#' @return The paths of the files written, invisibly.
#' @examples
#' d <- cepip_toy
#' res <- cepip(d$x, d$y, donor = d$donor, id = d$id, verbose = FALSE)
#' out <- file.path(tempdir(), "toy")
#' files <- cepip_save(res, out, plots = FALSE)
#' basename(files)
#' unlink(out, recursive = TRUE)
#' @export
cepip_save <- function(results, dir, arm = basename(dir), extras = NULL, plots = TRUE,
                       formats = c("png", "pdf"), log = NULL, max_cell_rows = 20000L) {
    if (inherits(results, "cepip")) results <- list(all = results)
    if (is.null(names(results)) || any(!nzchar(names(results))))
        stop("results must be a named list.", call. = FALSE)
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    files <- character(0)

    tabs <- list(); cells <- list()
    for (s in names(results)) {
        r <- results[[s]]; ex <- extras[[s]]
        if (!is.null(ex$screen)) {
            tabs[[sheet("01_screen", s)]] <- ex$screen$ranking
            tabs[[sheet("01_partners", s)]] <- ex$screen$partners
        }
        tabs[[sheet("02_prss", s)]] <- r$fit$prss
        cells[[sheet("03_fit", s)]] <- r$ep
        tabs[[sheet("04_ip", s)]] <- ip_table(r$ip)
        if (!is.null(r$skcs)) cells[[sheet("05_skcs", s)]] <- skcs_table(r$skcs)
        if (!is.null(r$band)) tabs[[sheet("07_band", s)]] <- band_table(r$band)
        if (!is.null(r$cmrf)) {
            tabs[[sheet("08_rskfcv_trace", s)]] <- gate_table(r$cmrf)
            cells[[sheet("10_module_cells", s)]] <- r$cmrf$cells
            tabs[[sheet("10_class_tests", s)]] <- class_table(r$cmrf)
        }
        if (!is.null(r$validate)) {
            tabs[[sheet("11_validate", s)]] <- validate_table(r$validate)
            cells[[sheet("11_validate_cells", s)]] <- r$validate$cells
        }
        if (!is.null(ex$donor)) tabs[[sheet("12_donor", s)]] <- ex$donor
        if (!is.null(ex$deg)) for (ct in names(ex$deg)) if (!is.null(ex$deg[[ct]]))
            tabs[[sheet(paste0("13_deg_", ct), s)]] <- ex$deg[[ct]]
        if (!is.null(ex$associate)) tabs[[sheet("14_associate", s)]] <- ex$associate
    }
    n_cell_rows <- sum(vapply(cells, nrow, integer(1)))
    if (n_cell_rows <= max_cell_rows) {
        tabs <- c(tabs, cells)
        cells <- list()
    }
    ord <- function(l) l[order(names(l))]
    files <- c(files, write_book(ord(tabs), file.path(dir, paste0(arm, "_cepip_tables.xlsx"))))
    if (length(cells))
        files <- c(files, write_book(ord(cells), file.path(dir, paste0(arm, "_cepip_cells.xlsx"))))

    set_txt <- c(sprintf("cepip %s | R %s | mgcv %s | written %s",
                         as.character(utils::packageVersion("cepip")),
                         paste(R.version$major, R.version$minor, sep = "."),
                         as.character(utils::packageVersion("mgcv")), format(Sys.time())))
    for (s in names(results)) {
        set_txt <- c(set_txt, "", paste0("== ", s),
                     utils::capture.output(print(results[[s]])),
                     paste("settings:", paste(names(results[[s]]$settings),
                                              unlist(lapply(results[[s]]$settings, format)),
                                              sep = " = ", collapse = "; ")))
    }
    writeLines(set_txt, file.path(dir, "00_settings.txt"))
    files <- c(files, file.path(dir, "00_settings.txt"))
    if (!is.null(log)) {
        writeLines(log, file.path(dir, "00_log.txt"))
        files <- c(files, file.path(dir, "00_log.txt"))
    }

    if (plots) {
        for (s in names(results)) {
            r <- results[[s]]; ex <- extras[[s]]
            pl <- list()
            if (!is.null(ex$screen)) pl[["01_screen"]] <- function() plot_screen(ex$screen)
            pl[["02_prss"]] <- function() plot_prss(r$fit)
            pl[["03_fit"]] <- function() plot_fit(r)
            if (isTRUE(r$ip$has_ip)) pl[["04_ip_residual"]] <- function() plot_residual(r$ip)
            if (!is.null(r$skcs)) pl[["05_skcs_pitcher"]] <- function() plot_pitcher(r$skcs, ip = r$ip)
            if (!is.null(r$band)) {
                pl[["06_copula"]] <- function() plot_copula(r$band)
                pl[["07_band"]] <- function() plot_band(r$band, r$cmrf)
            }
            if (!is.null(r$cmrf)) {
                pl[["08_rskfcv_first"]] <- function() plot_rskfcv(r$cmrf, "first")
                pl[["09_rskfcv_last"]] <- function() plot_rskfcv(r$cmrf, "last")
                pl[["10_module"]] <- function() plot_module(r$cmrf)
            }
            if (!is.null(r$validate)) pl[["11_validate"]] <- function() plot_validate(r$validate, r$cmrf)
            if (!is.null(ex$donor)) pl[["12_donor"]] <- function() plot_donor(ex$donor)
            if (!is.null(ex$deg)) for (ct in names(ex$deg)) if (!is.null(ex$deg[[ct]])) {
                local({ cc <- ct; pl[[paste0("13_deg_", cc)]] <<- function() plot_deg(ex$deg, cc) })
            }
            if (!is.null(ex$associate)) pl[["14_associate"]] <- function() plot_associate(ex$associate)
            for (nm in names(pl)) {
                p <- tryCatch(pl[[nm]](), error = function(e) {
                    warning(sprintf("Plot %s for %s failed: %s", nm, s, conditionMessage(e)),
                            call. = FALSE)
                    NULL
                })
                if (is.null(p)) next
                size <- attr(p, "cepip_size") %||% plot_size(nm)
                base <- file.path(dir, if (s == "all") nm else paste0(nm, "_", s))
                files <- c(files, save_plot(p, base, formats, size[1], size[2]))
            }
        }
    }
    invisible(files)
}

sheet <- function(step, set) {
    nm <- if (set == "all") step else paste0(step, "_", set)
    substr(nm, 1L, 31L)
}

write_book <- function(tabs, path) {
    wb <- openxlsx::createWorkbook()
    bold <- openxlsx::createStyle(textDecoration = "bold")
    for (nm in names(tabs)) {
        d <- as.data.frame(tabs[[nm]])
        openxlsx::addWorksheet(wb, nm)
        openxlsx::writeData(wb, nm, d, headerStyle = bold)
        openxlsx::freezePane(wb, nm, firstRow = TRUE)
    }
    openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
    path
}

plot_size <- function(nm) {
    if (grepl("^07_band", nm)) return(c(16, 12))
    if (grepl("^0[89]_rskfcv", nm)) return(c(18, 14))
    if (grepl("^11_validate", nm)) return(c(12, 5))
    c(7, 5.5)
}

save_plot <- function(p, base, formats, w, h) {
    out <- character(0)
    for (fm in formats) {
        f <- paste0(base, ".", fm)
        if (fm == "png") grDevices::png(f, width = w, height = h, units = "in", res = 200)
        else if (fm == "pdf") grDevices::pdf(f, width = w, height = h)
        else next
        if (inherits(p, "ggplot")) print(p) else { grid::grid.newpage(); grid::grid.draw(p) }
        grDevices::dev.off()
        out <- c(out, f)
    }
    out
}

ip_table <- function(ip) {
    i <- ip$iprs
    data.frame(has_IP = ip$has_ip, status = ip$status, IP_x = ip$ip_x,
               IP_x_CI_lower = ip$ip_ci[1], IP_x_CI_upper = ip$ip_ci[2],
               CI_lower_method = ip$ci_method[1], CI_upper_method = ip$ci_method[2],
               IP_y = ip$ip_y, IP_y_CI_lower = ip$ip_y_ci[1], IP_y_CI_upper = ip$ip_y_ci[2],
               n_crossings = ip$n_crossings, n_top_EP_cells = sum(ip$is_top),
               DE_pop = ip$de_pop, test = if (is.null(ip$test)) NA else ip$test$test,
               test_p = if (is.null(ip$test)) NA else ip$test$p,
               C1 = i$C1, C2 = i$C2, C3 = i$C3, C4 = i$C4, C5 = i$C5, IPRS = i$IPRS,
               IPRS_tier = i$tier, IPRS_flag = i$flag, IP_axis = ip$ip_axis,
               stringsAsFactors = FALSE)
}

skcs_table <- function(sk) {
    out <- sk$cells
    if (!is.null(sk$ep)) {
        m <- match(out$id, sk$ep$id)
        out$EP_core <- sk$ep$EP[m]
        out$is_TREP <- !is.na(m) & sk$ep$is_TREP[m] %in% TRUE
    }
    out
}

band_table <- function(b) {
    q <- b$qc
    data.frame(item = c("delta_base", "kappa", "bw", "coverage", "IP_x", "IP_y", "IP_u", "IP_v",
                        "spine_k", "spine_cells", "spine_monotone", "QC1_width_at_IP_u",
                        "width_without_widening", "QC2_widening_left_at_bw", "width_at_IP_u_pm_bw",
                        "QC3_max_kernel_slope", "n_OriRef", "n_candidates", "n_Original_only"),
               value = c(b$delta_base, b$kappa, b$bw, b$coverage, b$IP_x, b$IP_y, b$IP_u, b$IP_v,
                         b$spine$k, b$n_spine_cells, as.numeric(b$spine_monotone), q$width_at_ip,
                         q$width_base, q$widening_left_at_bw, q$width_at_bw, q$max_kernel_slope,
                         sum(b$cells$is_anchor), sum(b$cells$is_candidate),
                         sum(b$cells$is_original_only)),
               stringsAsFactors = FALSE)
}

gate_table <- function(cm) {
    g <- cm$gate
    n <- length(g$p)
    removed <- c(cm$cells$id[g$cand_row[g$removal_order]], rep(NA, max(0L, n - g$n_removed)))
    data.frame(iteration = seq_len(n), n_active = g$n_active, wilcoxon_p = g$p,
               dRMSE = g$effect, BV = g$bv, removed_after = removed[seq_len(n)],
               stringsAsFactors = FALSE)
}

class_table <- function(cm) {
    out <- data.frame(item = paste0("n_", names(cm$counts)), statistic = as.numeric(cm$counts),
                      p_value = NA_real_, stringsAsFactors = FALSE)
    out <- rbind(out, data.frame(item = c("n_Module", "gate_removals", "gate_p_final",
                                          "gate_dRMSE_final", "gate_BV"),
                                 statistic = c(sum(cm$cells$in_module), cm$gate$n_removed,
                                               cm$gate$p_final, cm$gate$effect_final, cm$gate$bv),
                                 p_value = NA_real_))
    if (!is.null(cm$tests) && nrow(cm$tests))
        out <- rbind(out, data.frame(item = cm$tests$test, statistic = cm$tests$statistic,
                                     p_value = cm$tests$p_value))
    out
}

validate_table <- function(v) {
    data.frame(item = c("k", "DE", "DE_pop", "n_Module", "n_TREP", "n_OriRef", "agreement"),
               value = c(v$fit$k, v$de, v$de_pop, nrow(v$cells), v$n_trep, v$n_oriref,
                         v$agreement),
               stringsAsFactors = FALSE)
}
