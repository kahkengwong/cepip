toy <- cepip_toy

test_that("PRSS search picks the lowest PRSS, ties to the smaller k", {
    f <- cepip_fit(toy$x, toy$y, donor = toy$donor)
    expect_s3_class(f, "cepip_fit")
    expect_equal(f$k, f$prss$k[which.min(f$prss$PRSS)])
    expect_true(f$donor_re)
    expect_equal(sum(f$data$w), length(unique(toy$donor)))
    f3 <- cepip_fit(toy$x, toy$y, donor = toy$donor, k = 5)
    expect_equal(f3$k, 5L)
})

test_that("without donors, DE equals the usual deviance explained", {
    f <- cepip_fit(toy$x, toy$y)
    expect_false(f$donor_re)
    expect_equal(f$de, f$de_flat, tolerance = 1e-10)
    ep <- cepip_ep(f)
    expect_equal(sum(ep$is_TREP), round(nrow(toy) * attr(ep, "de_pop")))
})

test_that("SKCS scores match their definitions", {
    set.seed(3)
    x <- c(rnorm(40), rep(0.5, 5)); y <- x + rnorm(45, 0, 0.5)
    sc <- cepip:::skcs_scores(x, y, block = 7)
    n <- length(x)
    # The Spearman contributions add up to the rank covariance.
    expect_equal(sum(sc$Si) / (stats::sd(rank(x)) * stats::sd(rank(y))),
                 stats::cor(x, y, method = "spearman"))
    # Each Kendall contribution lies in [-1, 1].
    expect_true(all(abs(sc$Ki) <= 1 + 1e-12))
    # Kendall contributions computed in blocks equal the full computation.
    sc2 <- cepip:::skcs_scores(x, y, block = 1000)
    expect_equal(sc$Ki, sc2$Ki)
})

test_that("the IP step finds an IP on the toy data and labels the cells", {
    ip <- cepip_ip(toy$x, toy$y, donor = toy$donor, ip_k = 6)
    expect_true(ip$has_ip)
    expect_true(ip$ip_ci[1] <= ip$ip_x && ip$ip_x <= ip$ip_ci[2])
    expect_equal(ip$iprs$IPRS, mean(unlist(ip$iprs[c("C1", "C2", "C3", "C4", "C5")])))
    expect_setequal(unique(ip$labels), c("Pre-IP", "Post-IP"))
    expect_equal(ip$labels == "Pre-IP", toy$y < ip$ip_y)
    ipx <- cepip_ip(toy$x, toy$y, donor = toy$donor, ip_k = 6, ip_axis = "x")
    expect_equal(ipx$labels == "Pre-IP", toy$x < ipx$ip_x)
})

test_that("a dataset without an IP is reported without an error", {
    set.seed(5)
    x <- runif(300, 0, 3); y <- 0.5 * x + rnorm(300, 0, 0.05)
    res <- cepip(x, y, cmrf = TRUE, verbose = FALSE)
    if (!isTRUE(res$ip$has_ip)) {
        expect_null(res$cmrf)
        expect_match(res$status, "no IP|too few")
    } else succeed()
})

test_that("the full method runs on the toy data", {
    res <- cepip(toy$x, toy$y, donor = toy$donor, id = toy$id, cmrf = TRUE, tests = FALSE,
                 verbose = FALSE)
    expect_s3_class(res, "cepip")
    cl <- res$cmrf$cells
    expect_equal(sum(res$cmrf$counts), nrow(toy))
    expect_true(all(cl$cell_class[cl$is_core & cl$pass_stage1] == "Original + Refined"))
    expect_true(all(cl$cell_class[cl$is_core & !cl$pass_stage1] == "Original only"))
    expect_true(all(!cl$is_core[cl$cell_class == "Recovered"]))
    g <- res$cmrf$gate
    expect_equal(length(g$rmse_A), 5 * 40)
    expect_equal(g$n_removed, sum(!g$pass))
    expect_true(g$p_final >= 0.05 || g$effect_final <= g$bv)
    expect_s3_class(res$validate, "cepip_validate")
    grp <- cepip_groups(res)
    expect_setequal(unique(grp$group), c("PostIP", "PreIP", "NMPost", "NMPre"))
    dp <- cepip_donor(res)
    expect_equal(nrow(dp), 12L)
    expect_equal(dp$pct_PostIP + dp$pct_PreIP + dp$pct_NMPost + dp$pct_NMPre, rep(100, 12))
})

test_that("the band follows its definition", {
    res <- cepip(toy$x, toy$y, donor = toy$donor, skcs = TRUE, verbose = FALSE)
    b <- cepip_band(res$skcs, res$ip, kappa = 0.3, bw = 0.15)
    cl <- b$cells
    expect_equal(cl$delta_adapt, b$delta_base * (1 + 0.3 * exp(-(cl$u - b$IP_u)^2 / (2 * 0.15^2))))
    expect_equal(cl$pass_stage1, cl$resid_abs <= cl$delta_adapt)
    expect_equal(max(cl$K), 1, tolerance = 0.05)
    expect_equal(b$qc$widening_left_at_bw, exp(-0.5))
})

test_that("the gate does not depend on the number of workers", {
    skip_on_cran()
    res <- cepip(toy$x, toy$y, donor = toy$donor, skcs = TRUE, verbose = FALSE)
    b <- cepip_band(res$skcs, res$ip)
    a1 <- cepip_cmrf(band = b, repeats = 10, tests = FALSE)
    a2 <- cepip_cmrf(band = b, repeats = 10, tests = FALSE, workers = 2)
    expect_equal(a1$cells$cell_class, a2$cells$cell_class)
    expect_equal(a1$gate$p, a2$gate$p)
})

test_that("association models reproduce lm() on z-scores", {
    set.seed(7)
    d <- data.frame(p = runif(40), a = rnorm(40), s = rbinom(40, 1, 0.5))
    d$o <- 2 * d$p + d$a + rnorm(40)
    r <- cepip_associate(d, "p", c("o", "a"), c("a", "s"))
    z <- as.data.frame(scale(d))
    m <- stats::lm(o ~ p + a + s, data = z)
    expect_equal(r$beta_std[1], unname(coef(m)["p"]))
    expect_equal(r$covariates[2], "s")
    expect_equal(r$FDR, stats::p.adjust(r$p_value, "BH"))
    r0 <- cepip_associate(d, "p", "o")
    m0 <- stats::lm(o ~ p, data = z)
    expect_equal(r0$beta_std, unname(coef(m0)["p"]))
    expect_equal(r0$delta_R2, summary(m0)$r.squared)
})

test_that("donor-unit DEGs follow the paired design", {
    res <- cepip(toy$x, toy$y, donor = toy$donor, id = toy$id, cmrf = TRUE, tests = FALSE,
                 verbose = FALSE)
    set.seed(11)
    g <- cepip_groups(res)
    m <- matrix(rpois(30 * nrow(toy), 1), nrow = 30,
                dimnames = list(paste0("G", 1:30), toy$id))
    # Three planted genes, so that BH over 30 genes can reach FDR < 0.05 with 12 donors.
    m[1:3, g$group == "PostIP"] <- m[1:3, g$group == "PostIP"] + 6
    m <- log1p(m)
    deg <- cepip_deg(res, m, contrast = "PostIP_vs_PreIP")
    r <- deg$PostIP_vs_PreIP
    expect_true(r$significant[r$gene == "G1"])
    expect_true(r$avg_log2FC[r$gene == "G1"] > 0)
})

test_that("the screen finds planted partner genes", {
    set.seed(2)
    n <- 300
    z <- rexp(n)
    m <- rbind(GOI = z, A = z + rnorm(n, 0, 0.1), B = z + rnorm(n, 0, 0.1),
               C = rexp(n), D = -z + rnorm(n, 0, 0.1) + 5)
    sc <- cepip_screen(m, goi = "GOI", min_mean = 0)
    expect_setequal(sc$partners$partner, c("A", "B", "D"))
    expect_equal(sc$ranking$n_DFG[sc$ranking$GOI == "GOI"], 3)
    sc2 <- cepip_screen(m, goi = "GOI", min_mean = 0, negative = FALSE)
    expect_setequal(sc2$partners$partner, c("A", "B"))
})
