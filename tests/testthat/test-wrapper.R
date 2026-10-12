###############################################################################
# RPBC WRAPPER TESTS
#
# Tests:
#   - ordinary RPBC
#   - input transformations: none, asinh, log1p
#   - covariate-adjusted RPBC
#   - Adaptive RPBC
#   - SingleCellExperiment
#   - SpatialExperiment
#   - wrapper-level error handling
###############################################################################


.make_rpbc_wrapper_object <- function(
    spatial = FALSE
) {

  set.seed(20261007)

  n_patients <- 4L
  cells_per_region <- 3L
  n_features <- 8L

  patient_levels <- paste0(
    "P",
    seq_len(n_patients)
  )

  patient <- factor(
    rep(
      patient_levels,
      each = 2L * cells_per_region
    ),
    levels = patient_levels
  )

  region <- factor(
    rep(
      rep(
        c("Core", "Edge"),
        each = cells_per_region
      ),
      times = n_patients
    ),
    levels = c(
      "Core",
      "Edge"
    )
  )

  n_cells <- length(patient)

  age_patient <- c(
    45,
    52,
    61,
    48
  )

  sex_patient <- factor(
    c(
      "F",
      "M",
      "F",
      "M"
    ),
    levels = c(
      "F",
      "M"
    )
  )

  patient_index <- match(
    as.character(patient),
    patient_levels
  )

  dat <- matrix(
    rnorm(
      n_features * n_cells,
      sd = 0.8
    ),
    nrow = n_features,
    ncol = n_cells
  )

  rownames(dat) <- paste0(
    "Feature",
    seq_len(n_features)
  )

  colnames(dat) <- paste0(
    "Cell",
    seq_len(n_cells)
  )

  region_code <- ifelse(
    region == "Edge",
    1,
    -1
  )

  age_centered <-
    age_patient -
    mean(age_patient)

  for (g in seq_len(n_features)) {

    dat[g, ] <- dat[g, ] +
      5 +
      0.5 * region_code +
      0.04 *
      age_centered[patient_index]
  }

  metadata <- data.frame(
    patient = patient,
    region = region,
    age = age_patient[patient_index],
    sex = sex_patient[patient_index],
    row.names = colnames(dat)
  )

  if (!spatial) {

    return(
      SingleCellExperiment::SingleCellExperiment(
        assays = list(
          counts = dat
        ),
        colData = metadata
      )
    )
  }

  coordinates <- cbind(
    x = seq_len(n_cells),
    y = rev(
      seq_len(n_cells)
    )
  )

  rownames(coordinates) <-
    colnames(dat)

  SpatialExperiment::SpatialExperiment(
    assays = list(
      counts = dat
    ),
    colData = metadata,
    spatialCoords = coordinates
  )
}


###############################################################################
# Ordinary RPBC
###############################################################################

test_that("rpbc wrapper runs ordinary RPBC", {

  sce <- .make_rpbc_wrapper_object()

  original <- SummarizedExperiment::assay(
    sce,
    "counts"
  )

  result <- rpbc(
    object = sce,
    patient = "patient",
    region = "region"
  )

  expect_s4_class(
    result,
    "SingleCellExperiment"
  )

  expect_true(
    "rpbc" %in%
      SummarizedExperiment::assayNames(
        result
      )
  )

  harmonized <- SummarizedExperiment::assay(
    result,
    "rpbc"
  )

  expect_equal(
    dim(harmonized),
    dim(original)
  )

  expect_identical(
    dimnames(harmonized),
    dimnames(original)
  )

  expect_true(
    all(
      is.finite(
        harmonized
      )
    )
  )

  expect_equal(
    SummarizedExperiment::assay(
      result,
      "counts"
    ),
    original
  )
})

###############################################################################
# Input transformations
###############################################################################

test_that("rpbc wrapper uses the selected assay unchanged with transform = none", {

  sce <- .make_rpbc_wrapper_object()

  counts_input <- SummarizedExperiment::assay(
    sce,
    "counts"
  )

  wrapped <- rpbc(
    object = sce,
    patient = "patient",
    region = "region",
    assay_name = "counts",
    transform = "none",
    output_assay = "rpbc_none",
    return_fit = TRUE
  )

  direct <- .rpbc_fit_matrix(
    dat = counts_input,
    patient = SummarizedExperiment::colData(sce)$patient,
    region = SummarizedExperiment::colData(sce)$region
  )

  expect_equal(
    wrapped$fit$data,
    direct$data,
    tolerance = 1e-12
  )
})


test_that("rpbc wrapper applies asinh transformation correctly", {

  sce <- .make_rpbc_wrapper_object()

  counts_input <- SummarizedExperiment::assay(
    sce,
    "counts"
  )

  wrapped <- rpbc(
    object = sce,
    patient = "patient",
    region = "region",
    assay_name = "counts",
    transform = "asinh",
    cofactor = 5,
    output_assay = "rpbc_asinh",
    return_fit = TRUE
  )

  direct <- .rpbc_fit_matrix(
    dat = asinh(
      counts_input / 5
    ),
    patient = SummarizedExperiment::colData(sce)$patient,
    region = SummarizedExperiment::colData(sce)$region
  )

  expect_equal(
    wrapped$fit$data,
    direct$data,
    tolerance = 1e-12
  )
})


test_that("rpbc wrapper applies log1p transformation correctly", {

  sce <- .make_rpbc_wrapper_object()

  counts_input <- SummarizedExperiment::assay(
    sce,
    "counts"
  )

  wrapped <- rpbc(
    object = sce,
    patient = "patient",
    region = "region",
    assay_name = "counts",
    transform = "log1p",
    output_assay = "rpbc_log1p",
    return_fit = TRUE
  )

  direct <- .rpbc_fit_matrix(
    dat = log1p(
      counts_input
    ),
    patient = SummarizedExperiment::colData(sce)$patient,
    region = SummarizedExperiment::colData(sce)$region
  )

  expect_equal(
    wrapped$fit$data,
    direct$data,
    tolerance = 1e-12
  )
})

test_that("rpbc wrapper rejects an invalid asinh cofactor", {

  sce <- .make_rpbc_wrapper_object()

  expect_error(
    rpbc(
      object = sce,
      patient = "patient",
      region = "region",
      transform = "asinh",
      cofactor = 0
    ),
    "positive finite"
  )
})


test_that("rpbc wrapper rejects negative input for log1p transformation", {

  sce <- .make_rpbc_wrapper_object()

  dat_bad <- SummarizedExperiment::assay(
    sce,
    "counts"
  )

  dat_bad[1, 1] <- -1

  SummarizedExperiment::assay(
    sce,
    "counts"
  ) <- dat_bad

  expect_error(
    rpbc(
      object = sce,
      patient = "patient",
      region = "region",
      transform = "log1p"
    ),
    "requires non-negative input"
  )
})

###############################################################################
# Covariate-adjusted RPBC
###############################################################################

test_that("rpbc wrapper runs covariate-adjusted RPBC", {

  sce <- .make_rpbc_wrapper_object()

  result <- rpbc(
    object = sce,
    patient = "patient",
    region = "region",
    covariates = c(
      "age",
      "sex"
    ),
    output_assay = "rpbc_covariate",
    return_fit = TRUE
  )

  expect_named(
    result,
    c(
      "object",
      "fit"
    )
  )

  expect_s4_class(
    result$object,
    "SingleCellExperiment"
  )

  expect_true(
    result$fit$covariate_adjusted
  )

  expect_equal(
    result$fit$q_covariates,
    2L
  )

  expect_true(
    "rpbc_covariate" %in%
      SummarizedExperiment::assayNames(
        result$object
      )
  )

  expect_true(
    all(
      is.finite(
        SummarizedExperiment::assay(
          result$object,
          "rpbc_covariate"
        )
      )
    )
  )
})


###############################################################################
# Adaptive RPBC
###############################################################################

test_that("rpbc wrapper runs Adaptive RPBC", {

  sce <- .make_rpbc_wrapper_object()

  result <- rpbc(
    object = sce,
    patient = "patient",
    region = "region",
    adaptive = TRUE,
    output_assay = "rpbc_adaptive",
    return_fit = TRUE
  )

  expect_named(
    result,
    c(
      "object",
      "fit"
    )
  )

  expect_true(
    "rpbc_adaptive" %in%
      SummarizedExperiment::assayNames(
        result$object
      )
  )

  expect_equal(
    result$fit$delta_adaptive_k,
    1
  )

  expect_equal(
    result$fit$delta_bounds,
    c(
      0.5,
      2.0
    )
  )

  expect_true(
    all(
      result$fit$delta_star >= 0.5
    )
  )

  expect_true(
    all(
      result$fit$delta_star <= 2.0
    )
  )

  expect_true(
    all(
      is.finite(
        SummarizedExperiment::assay(
          result$object,
          "rpbc_adaptive"
        )
      )
    )
  )
})


###############################################################################
# SpatialExperiment
###############################################################################

test_that("rpbc wrapper accepts SpatialExperiment objects", {

  spe <- .make_rpbc_wrapper_object(
    spatial = TRUE
  )

  result <- rpbc(
    object = spe,
    patient = "patient",
    region = "region"
  )

  expect_s4_class(
    result,
    "SpatialExperiment"
  )

  expect_true(
    "rpbc" %in%
      SummarizedExperiment::assayNames(
        result
      )
  )

  expect_equal(
    dim(
      SummarizedExperiment::assay(
        result,
        "rpbc"
      )
    ),
    dim(
      SummarizedExperiment::assay(
        spe,
        "counts"
      )
    )
  )
})


###############################################################################
# Adaptive and covariate-adjusted modes cannot be combined
###############################################################################

test_that("rpbc wrapper rejects Adaptive plus covariates", {

  sce <- .make_rpbc_wrapper_object()

  expect_error(
    rpbc(
      object = sce,
      patient = "patient",
      region = "region",
      covariates = "age",
      adaptive = TRUE
    ),
    "separate extensions"
  )
})


###############################################################################
# Missing assay
###############################################################################

test_that("rpbc wrapper rejects a missing assay", {

  sce <- .make_rpbc_wrapper_object()

  expect_error(
    rpbc(
      object = sce,
      patient = "patient",
      region = "region",
      assay_name = "does_not_exist"
    ),
    "was not found"
  )
})


###############################################################################
# Missing patient metadata
###############################################################################

test_that("rpbc wrapper rejects a missing patient column", {

  sce <- .make_rpbc_wrapper_object()

  expect_error(
    rpbc(
      object = sce,
      patient = "wrong_patient",
      region = "region"
    ),
    "was not found"
  )
})


###############################################################################
# Missing region metadata
###############################################################################

test_that("rpbc wrapper rejects a missing region column", {

  sce <- .make_rpbc_wrapper_object()

  expect_error(
    rpbc(
      object = sce,
      patient = "patient",
      region = "wrong_region"
    ),
    "was not found"
  )
})


###############################################################################
# Existing output assay
###############################################################################

test_that("rpbc wrapper protects an existing assay from overwrite", {

  sce <- .make_rpbc_wrapper_object()

  SummarizedExperiment::assay(
    sce,
    "rpbc"
  ) <- SummarizedExperiment::assay(
    sce,
    "counts"
  )

  expect_error(
    rpbc(
      object = sce,
      patient = "patient",
      region = "region"
    ),
    "already exists"
  )
})
