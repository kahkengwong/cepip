# Reproduces the published HKG8 control arm (about a minute). Runs only when
# CEPIP_EXAMPLE_DATA points to a folder with hkg8_cells.csv and donors.csv.
# The full reproduction of both arms is inst/scripts/reproduce_seaad.R.
test_that("the HKG8 control arm reproduces the published values", {
    dd <- Sys.getenv("CEPIP_EXAMPLE_DATA")
    skip_if(!nzchar(dd) || !file.exists(file.path(dd, "hkg8_cells.csv")), "example data not set")
    hk <- utils::read.csv(file.path(dd, "hkg8_cells.csv"), check.names = FALSE)
    ref <- utils::read.csv(system.file("reference", "seaad_reference.csv", package = "cepip"))
    ref <- ref[ref$arm == "HKG8", ]
    val <- function(set, item) ref$value[ref$region_set == set & ref$item == item]
    HK8 <- c("EEF1A1", "RPL13A", "B2M", "YWHAZ", "FTL", "PSAP", "RAB7A", "HNRNPA2B1")
    for (set in c("NeoC", "NonNeoC")) {
        d <- hk[hk$region_set == set, ]
        r <- cepip(d$ACTB, cepip_composite(d, HK8), donor = d$donor, id = d$barcode,
                   cmrf = TRUE, ip_k = 10, y_spine = d$composite_part8, tests = FALSE,
                   verbose = FALSE)
        expect_equal(r$ip$ip_x, val(set, "IP_Value"), tolerance = 1e-8)
        expect_equal(r$ip$iprs$IPRS, val(set, "IPRS"), tolerance = 1e-8)
        expect_equal(r$skcs$n_core, val(set, "n_core"))
        expect_equal(r$band$delta_base, val(set, "delta_base"), tolerance = 1e-8)
        for (cc in names(r$cmrf$counts))
            expect_equal(r$cmrf$counts[[cc]], val(set, paste0("n_", cc)))
        expect_equal(r$cmrf$gate$n_removed, val(set, "n_removed"))
        expect_equal(100 * r$validate$agreement, val(set, "agreement_pct"), tolerance = 1e-3)
    }
})
