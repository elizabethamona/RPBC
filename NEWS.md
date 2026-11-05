# RPBC 0.99.1

- Added an Installation section to the package vignette with Bioconductor installation instructions using `BiocManager`.
- Replaced `protected_rows <- 1:2` with `protected_rows <- c(1, 2)` in `R/rpbc-engine.R`; this does not change the statistical behavior of RPBC.

# RPBC 0.99.0

- Initial Bioconductor submission.
- Implemented Ordinary RPBC for preserving the patient-balanced population-average regional effect while harmonizing patient-associated nuisance variation.
- Implemented Adaptive RPBC with group-size-dependent stabilization of patient-region scale adjustments.
- Implemented Covariate-adjusted RPBC for protecting prespecified patient-level covariate-associated variation.
- Added support for `SingleCellExperiment` and `SpatialExperiment` objects.
- Added the `haley_gbm_example` dataset, package vignette, documentation, and test suite.
