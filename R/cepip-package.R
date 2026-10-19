#' cepip: CEP-IP with CMRF for single-cell co-expression
#'
#' The main function is [cepip()]. The steps can also be run one by one:
#' [cepip_screen()], [cepip_composite()], [cepip_fit()], [cepip_ep()],
#' [cepip_ip()], [cepip_skcs()], [cepip_band()], [cepip_cmrf()] and
#' [cepip_validate()], then [cepip_donor()], [cepip_deg()] and
#' [cepip_associate()]. [cepip_save()] writes the tables and plots of one
#' analysis to a folder.
#'
#' @keywords internal
#' @importFrom stats median
#' @importFrom mgcv gam
"_PACKAGE"

#' Simulated example data
#'
#' 1,200 simulated cells from 12 donors with a gene of interest `x`, eight
#' partner genes `P1` to `P8` and their composite `y` (the row mean of the
#' partners). The data are simulated (see `data-raw/make_toy.R`) and run
#' through the whole method in under a minute.
#'
#' @format A data frame with 1,200 rows and columns `id`, `donor`, `region`,
#'   `x`, `P1` to `P8` and `y`.
"cepip_toy"

utils::globalVariables(c("w"))
