#' Spearman-Kendall dual-filter screen for partner genes
#'
#' Screens candidate genes of interest (GOIs) for partner genes (dual-filter
#' genes, DFGs) across the cells of one population, as in the published CEP-IP
#' and the SEA-AD analysis.
#'
#' 1. Donors with fewer than `min_donor_cells` cells are dropped.
#' 2. A GOI is screened if it is expressed (value > 0) in at least `min_pct` of
#'    the cells and its mean over the expressing cells is at least `min_mean`.
#'    In the SEA-AD data (ln(CP10K + 1)) these defaults kept 2,500 to 4,000
#'    genes per population. Adjust them for other data and normalizations.
#' 3. Partner genes must be expressed in at least `min_partner_cells` cells.
#' 4. Spearman rs is computed for every GOI-partner pair (average ranks, zeros
#'    included). Kendall tau-b is computed for the pairs with |rs| >= `rs`.
#' 5. A partner is a DFG when rs >= `rs` and tau-b >= `tau` (positive), or, with
#'    `negative = TRUE`, rs <= -`rs` and tau-b <= -`tau` (negative).
#'
#' GOIs are ranked by the number of DFGs, then by the largest rs. The choice of
#' the GOI and of the DFG panel for [cepip_composite()] is left to the user.
#'
#' Run time grows with cells x GOIs x partners. For large populations,
#' `max_cells` screens a random subsample (stratified by `strata` when given).
#' Kendall tau-b uses `pcaPP::cor.fk()` when pcaPP is installed (fast, same
#' result), else [stats::cor()].
#'
#' @param expr Expression matrix, genes x cells (dense or sparse), with gene
#'   names as row names.
#' @param goi Candidate GOIs. `NULL` screens every gene.
#' @param donor Optional donor label per cell.
#' @param strata Optional stratum per cell (for example region) for `max_cells`.
#' @param min_pct,min_mean GOI expression gate.
#' @param min_partner_cells Minimum number of expressing cells of a partner.
#' @param min_donor_cells Minimum number of cells per donor.
#' @param rs,tau Spearman and Kendall thresholds.
#' @param negative Keep negative DFGs.
#' @param max_cells Optional number of cells to screen.
#' @param seed Seed of the subsample.
#' @param goi_block,partner_block Block sizes for the Spearman step (memory).
#' @param verbose Print progress messages.
#' @return An object of class `cepip_screen` with `ranking` (one row per
#'   screened GOI with at least one DFG), `partners` (one row per GOI-DFG
#'   pair), `gate` (the GOI gate of every candidate) and `cells` (the cells
#'   screened).
#' @examples
#' d <- cepip_toy
#' expr <- t(as.matrix(d[, c("x", paste0("P", 1:8))]))
#' colnames(expr) <- d$id
#' sc <- cepip_screen(expr, goi = "x", donor = d$donor)
#' sc$ranking
#' head(sc$partners)
#' @export
cepip_screen <- function(expr, goi = NULL, donor = NULL, strata = NULL,
                         min_pct = 0.25, min_mean = 0.30, min_partner_cells = 10L,
                         min_donor_cells = 5L, rs = 0.6, tau = 0.5, negative = TRUE,
                         max_cells = NULL, seed = 42L, goi_block = 500L,
                         partner_block = 2000L, verbose = FALSE) {
    if (is.null(rownames(expr))) stop("expr needs gene names as row names.", call. = FALSE)
    rng <- save_rng()
    on.exit(restore_rng(rng), add = TRUE)
    n0 <- ncol(expr)
    keep <- rep(TRUE, n0)
    if (!is.null(donor)) {
        tab <- table(donor)
        keep <- as.character(donor) %in% names(tab)[tab >= min_donor_cells]
    }
    cells <- which(keep)
    if (!is.null(max_cells) && length(cells) > max_cells)
        cells <- subsample_cells(cells, if (is.null(strata)) NULL else strata[cells],
                                 max_cells, seed)
    X <- expr[, cells, drop = FALSE]
    n <- ncol(X)
    msg(verbose, "  Screen: %d cells", n)

    genes <- rownames(X)
    cand <- if (is.null(goi)) genes else intersect(goi, genes)
    nnz <- if (inherits(X, "Matrix")) Matrix::rowSums(X > 0) else rowSums(X > 0)
    sums <- if (inherits(X, "Matrix")) Matrix::rowSums(X) else rowSums(X)
    gate <- data.frame(GOI = cand, pct_expr = as.numeric(nnz[cand]) / n,
                       mean_expressing = ifelse(nnz[cand] > 0,
                                                as.numeric(sums[cand]) / as.numeric(nnz[cand]), 0),
                       stringsAsFactors = FALSE)
    gate$passes <- gate$pct_expr >= min_pct & gate$mean_expressing >= min_mean
    gois <- gate$GOI[gate$passes]
    partners <- genes[nnz >= min_partner_cells]
    msg(verbose, "  %d GOIs pass the gate | %d partner genes", length(gois), length(partners))

    dense_cols <- function(g) {
        m <- X[g, , drop = FALSE]
        t(as.matrix(m))
    }
    rank_centre <- function(M) {
        R <- apply(M, 2L, rank, ties.method = "average")
        if (is.null(dim(R))) R <- matrix(R, ncol = ncol(M))
        R <- sweep(R, 2L, colMeans(R))
        nr <- sqrt(colSums(R^2))
        nr[nr == 0] <- 1
        sweep(R, 2L, nr, "/")
    }
    hits <- list()
    for (ps in seq.int(1L, length(partners), by = partner_block)) {
        pg <- partners[ps:min(ps + partner_block - 1L, length(partners))]
        P <- rank_centre(dense_cols(pg))
        for (gs in seq.int(1L, length(gois), by = goi_block)) {
            gg <- gois[gs:min(gs + goi_block - 1L, length(gois))]
            G <- rank_centre(dense_cols(gg))
            C <- crossprod(G, P)
            dimnames(C) <- list(gg, pg)
            w <- which(abs(C) >= rs, arr.ind = TRUE)
            if (nrow(w)) {
                h <- data.frame(GOI = gg[w[, 1]], partner = pg[w[, 2]], rs = C[w],
                                stringsAsFactors = FALSE)
                hits[[length(hits) + 1L]] <- h[h$GOI != h$partner, ]
            }
        }
        msg(verbose, "  partners %d to %d done", ps, ps + length(pg) - 1L)
    }
    sp <- if (length(hits)) do.call(rbind, hits) else
        data.frame(GOI = character(0), partner = character(0), rs = numeric(0))
    if (!negative) sp <- sp[sp$rs >= rs, ]

    ktau <- if (requireNamespace("pcaPP", quietly = TRUE)) pcaPP::cor.fk
            else function(a, b) stats::cor(a, b, method = "kendall")
    sp$tau <- if (nrow(sp)) vapply(seq_len(nrow(sp)), function(i) {
        a <- as.numeric(X[sp$GOI[i], ]); b <- as.numeric(X[sp$partner[i], ])
        if (stats::sd(b) == 0) 0 else ktau(a, b)
    }, numeric(1)) else numeric(0)
    sp$direction <- ifelse(sp$rs > 0, "positive", "negative")
    dfg <- sp[(sp$rs >= rs & sp$tau >= tau) | (sp$rs <= -rs & sp$tau <= -tau), ]

    rk <- if (nrow(dfg)) {
        s <- split(dfg, dfg$GOI)
        r <- do.call(rbind, lapply(names(s), function(g) {
            spg <- sp[sp$GOI == g, ]
            data.frame(GOI = g, n_DFG = nrow(s[[g]]),
                       n_positive = sum(s[[g]]$direction == "positive"),
                       n_negative = sum(s[[g]]$direction == "negative"),
                       max_rs = max(spg$rs), max_tau = max(spg$tau),
                       stringsAsFactors = FALSE)
        }))
        r <- merge(r, gate[, c("GOI", "pct_expr", "mean_expressing")], by = "GOI")
        r[order(-r$n_DFG, -r$max_rs), ]
    } else data.frame()
    rownames(rk) <- NULL
    dfg <- dfg[order(dfg$GOI, -dfg$rs), ]
    rownames(dfg) <- NULL
    structure(list(ranking = rk, partners = dfg, spearman = sp, gate = gate,
                   cells = colnames(X) %||% cells, n_cells = n,
                   settings = list(min_pct = min_pct, min_mean = min_mean,
                                   min_partner_cells = min_partner_cells,
                                   min_donor_cells = min_donor_cells, rs = rs, tau = tau,
                                   negative = negative, max_cells = max_cells, seed = seed)),
              class = "cepip_screen")
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# Subsample of max_cells cells. With strata, each stratum gets an equal share;
# strata smaller than the share keep all their cells and the remainder is
# split over the larger strata.
subsample_cells <- function(cells, strata, max_cells, seed) {
    set.seed(seed)
    if (is.null(strata)) return(sort(sample(cells, max_cells)))
    lv <- sort(unique(as.character(strata)))
    sz <- table(factor(as.character(strata), levels = lv))
    target <- max_cells %/% length(lv)
    small <- names(sz)[sz <= target]
    big <- setdiff(lv, small)
    rest <- max_cells - sum(sz[small])
    share <- stats::setNames(rep(rest %/% max(1L, length(big)), length(big)), big)
    extra <- rest - sum(share)
    if (extra > 0L) share[seq_len(extra)] <- share[seq_len(extra)] + 1L
    out <- unlist(lapply(lv, function(s) {
        idx <- cells[as.character(strata) == s]
        if (s %in% small) idx else sample(idx, min(length(idx), share[[s]]))
    }))
    sort(out)
}
