#' Download the SEA-AD example data
#'
#' Downloads the example data used in the paper from Hugging Face once and
#' keeps them in a local cache. The data come from the Seattle Alzheimer's
#' Disease Brain Cell Atlas (SEA-AD) and are subject to the Allen Institute
#' Terms of Use (noncommercial research use, with citation); see
#' <https://alleninstitute.org/legal/terms-of-use> and
#' <https://alleninstitute.org/legal/citation-policy>.
#'
#' Files:
#' * `"pvm3"`: the 24,463 Micro-PVM_3-SEAAD cells with MT-CO3 > 0 (NeoC 13,331,
#'   NonNeoC 11,132): barcode, region set, region, donor, MT-CO3, the eight
#'   mitochondrial partner genes, the composite of the data export
#'   (`composite_part8`), ACTB and the eight other housekeeping genes, all as
#'   ln(CP10K + 1).
#' * `"hkg8"`: the 2,511 cells of the housekeeping-gene control arm of the
#'   paper (ACTB > 0 cells of the M4 arm): ACTB, the eight housekeeping genes
#'   and their export composite.
#' * `"donors"`: donor metadata (covariates and outcomes) and the total number
#'   of Micro-PVM_3-SEAAD cells per donor and region set.
#' * `"deg_NeoC"`, `"deg_NonNeoC"`: sparse ln(CP10K + 1) matrices (genes x cells)
#'   of the PVM_3 cells for [cepip_deg()].
#'
#' The files are kept in `cache` (by default the user cache folder of R, see
#' [tools::R_user_dir()]), so they are downloaded only once.
#' [cepip_example_clear()] deletes them. When the download fails (for example
#' without an internet connection), the function gives a message and returns
#' `NULL`.
#'
#' @param which Which file to load.
#' @param cache Cache folder.
#' @param repo Hugging Face dataset repository.
#' @param quiet Suppress the download progress.
#' @return A data frame (or a sparse matrix for the DEG files), or `NULL` when
#'   the download fails.
#' @examples
#' \donttest{
#' hk <- cepip_example("hkg8", cache = file.path(tempdir(), "cepip_cache"))
#' if (!is.null(hk)) dim(hk)
#' }
#' @export
cepip_example <- function(which = c("pvm3", "hkg8", "donors", "deg_NeoC", "deg_NonNeoC"),
                          cache = tools::R_user_dir("cepip", "cache"),
                          repo = "kahkengwong/cepip_example", quiet = FALSE) {
    which <- match.arg(which)
    file <- example_files[[which]]
    path <- file.path(cache, file)
    if (!file.exists(path)) {
        dir.create(cache, showWarnings = FALSE, recursive = TRUE)
        url <- sprintf("https://huggingface.co/datasets/%s/resolve/main/%s", repo, file)
        old <- options(timeout = max(3600, getOption("timeout")))
        on.exit(options(old), add = TRUE)
        tmp <- paste0(path, ".part")
        ok <- tryCatch(utils::download.file(url, tmp, mode = "wb", quiet = quiet) == 0L,
                       error = function(e) FALSE, warning = function(w) FALSE)
        if (!ok || !file.exists(tmp)) {
            unlink(tmp)
            message("The example data could not be downloaded from ", url,
                    ". Check the internet connection and try again.")
            return(invisible(NULL))
        }
        file.rename(tmp, path)
    }
    message("SEA-AD data: Allen Institute for Brain Science. Use is subject to the Allen ",
            "Institute Terms of Use and Citation Policy. Please cite SEA-AD (Gabitto et al., ",
            "Nature Neuroscience 2024) and the Allen Brain Map (https://portal.brain-map.org).")
    if (grepl("\\.rds$", file)) return(readRDS(path))
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

example_files <- list(pvm3 = "pvm3_cells.csv", hkg8 = "hkg8_cells.csv",
                      donors = "donors.csv", deg_NeoC = "pvm3_deg_NeoC.rds",
                      deg_NonNeoC = "pvm3_deg_NonNeoC.rds")

#' Delete the cached example data
#'
#' Deletes the files that [cepip_example()] saved in the cache folder, and the
#' folder itself when it is then empty.
#'
#' @param cache Cache folder (the same as in [cepip_example()]).
#' @return The paths of the deleted files, invisibly.
#' @examples
#' cepip_example_clear(cache = file.path(tempdir(), "cepip_cache"))
#' @export
cepip_example_clear <- function(cache = tools::R_user_dir("cepip", "cache")) {
    paths <- file.path(cache, c(unlist(example_files, use.names = FALSE),
                                paste0(unlist(example_files, use.names = FALSE), ".part")))
    gone <- paths[file.exists(paths)]
    unlink(gone)
    if (dir.exists(cache) && !length(list.files(cache, all.files = TRUE, no.. = TRUE)))
        unlink(cache, recursive = TRUE)
    message(length(gone), " cached file(s) deleted from ", cache, ".")
    invisible(gone)
}
