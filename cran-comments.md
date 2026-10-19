## Submission

This is a new submission.

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE: New submission.

## Test environments

* Fedora, Podman container, R 4.6.1
* win-builder, R-devel

## Notes for the reviewers

* The example data of the paper (about 85 MB) are not in the package. `cepip_example()` downloads them from a Hugging Face dataset into `tools::R_user_dir("cepip", "cache")`, and `cepip_example_clear()` deletes them. When the download fails, `cepip_example()` gives a message and returns `NULL`. Its example is in `\donttest{}` and writes to `tempdir()`.
* The help examples use the simulated dataset `cepip_toy` that is bundled with the package. Examples that run the cross-validation gate take more than 5 seconds and are in `\donttest{}`.
* Words such as CEP-IP, SKCS, CMRF, REML and PRSS are names and abbreviations of the method, not misspellings.
