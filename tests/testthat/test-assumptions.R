###############################################################################
# RPBC ASSUMPTION AND ERROR-HANDLING TESTS
#
# Covers:
#   - ordinary RPBC input assumptions
#   - paired Core/Edge design
#   - patient-region sample-size requirement
#   - covariate input assumptions
#   - covariate design identifiability
#   - explicit covariates = NULL equivalence
###############################################################################


###############################################################################
# Reusable valid test data
###############################################################################

.make_rpbc_assumption_data <- function() {

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

  dat <- matrix(
    rnorm(
      n_features * length(patient)
    ),
    nrow = n_features,
    ncol = length(patient)
  )

  rownames(dat) <- paste0(
    "Feature",
    seq_len(n_features)
  )

  colnames(dat) <- paste0(
    "Cell",
    seq_len(ncol(dat))
  )

  covariates <- data.frame(
    age = c(
      45,
      52,
      61,
      48
    ),
    sex = factor(
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
    ),
    row.names = patient_levels
  )

  list(
    dat = dat,
    patient = patient,
    region = region,
    covariates = covariates,
    patient_levels = patient_levels
  )
}


###############################################################################
# Ordinary RPBC assumptions
###############################################################################

test_that("RPBC rejects non-numeric feature data", {

  x <- .make_rpbc_assumption_data()

  dat_bad <- matrix(
    "not_numeric",
    nrow = nrow(x$dat),
    ncol = ncol(x$dat)
  )

  expect_error(
    .rpbc_fit_matrix(
      dat = dat_bad,
      patient = x$patient,
      region = x$region
    ),
    "numeric matrix"
  )
})


test_that("RPBC rejects mismatched patient and region lengths", {

  x <- .make_rpbc_assumption_data()

  expect_error(
    .rpbc_fit_matrix(
      dat = x$dat,
      patient = x$patient[-1],
      region = x$region
    ),
    "one value for every column"
  )

  expect_error(
    .rpbc_fit_matrix(
      dat = x$dat,
      patient = x$patient,
      region = x$region[-1]
    ),
    "one value for every column"
  )
})


test_that("RPBC rejects non-finite feature values", {

  x <- .make_rpbc_assumption_data()

  dat_bad <- x$dat
  dat_bad[1, 1] <- Inf

  expect_error(
    .rpbc_fit_matrix(
      dat = dat_bad,
      patient = x$patient,
      region = x$region
    ),
    "NA, NaN, or Inf"
  )
})


test_that("RPBC rejects missing patient or region labels", {

  x <- .make_rpbc_assumption_data()

  patient_bad <- as.character(
    x$patient
  )

  patient_bad[1] <- NA

  expect_error(
    .rpbc_fit_matrix(
      dat = x$dat,
      patient = patient_bad,
      region = x$region
    ),
    "cannot contain missing values"
  )

  region_bad <- as.character(
    x$region
  )

  region_bad[1] <- NA

  expect_error(
    .rpbc_fit_matrix(
      dat = x$dat,
      patient = x$patient,
      region = region_bad
    ),
    "cannot contain missing values"
  )
})


test_that("RPBC requires exactly the specified two region levels", {

  x <- .make_rpbc_assumption_data()

  region_bad <- as.character(
    x$region
  )

  region_bad[
    region_bad == "Edge"
  ] <- "Margin"

  expect_error(
    .rpbc_fit_matrix(
      dat = x$dat,
      patient = x$patient,
      region = region_bad,
      core_level = "Core",
      edge_level = "Edge"
    ),
    "must contain exactly the two specified levels"
  )
})


test_that("RPBC requires paired Core and Edge observations for every patient", {

  x <- .make_rpbc_assumption_data()

  region_bad <- as.character(
    x$region
  )

  region_bad[
    as.character(x$patient) == "P1"
  ] <- "Core"

  expect_error(
    .rpbc_fit_matrix(
      dat = x$dat,
      patient = x$patient,
      region = region_bad,
      core_level = "Core",
      edge_level = "Edge"
    ),
    "Every patient must have observations in both"
  )
})


test_that("RPBC requires at least two cells in every patient-region group", {

  x <- .make_rpbc_assumption_data()

  keep <- rep(
    TRUE,
    length(x$patient)
  )

  p1_edge <- which(
    as.character(x$patient) == "P1" &
      as.character(x$region) == "Edge"
  )

  keep[
    p1_edge[-1]
  ] <- FALSE

  expect_error(
    .rpbc_fit_matrix(
      dat = x$dat[
        ,
        keep,
        drop = FALSE
      ],
      patient = x$patient[keep],
      region = x$region[keep],
      core_level = "Core",
      edge_level = "Edge"
    ),
    "at least two cells"
  )
})


###############################################################################
# Ordinary RPBC equivalence with covariates = NULL
###############################################################################

test_that("explicit covariates = NULL reproduces ordinary RPBC", {

  x <- .make_rpbc_assumption_data()

  fit_default <- .rpbc_fit_matrix(
    dat = x$dat,
    patient = x$patient,
    region = x$region
  )

  fit_null <- .rpbc_fit_matrix(
    dat = x$dat,
    patient = x$patient,
    region = x$region,
    covariates = NULL
  )

  expect_false(
    fit_null$covariate_adjusted
  )

  expect_equal(
    fit_null$q_covariates,
    0L
  )

  expect_equal(
    fit_null$residual_patient_df,
    3L
  )

  expect_equal(
    fit_default$data,
    fit_null$data,
    tolerance = 1e-12
  )

  expect_equal(
    fit_default$gamma_star,
    fit_null$gamma_star,
    tolerance = 1e-12
  )

  expect_equal(
    fit_default$kappa_star,
    fit_null$kappa_star,
    tolerance = 1e-12
  )

  expect_equal(
    fit_default$delta_star,
    fit_null$delta_star,
    tolerance = 1e-12
  )
})


###############################################################################
# Covariate assumptions
###############################################################################

test_that("covariate-adjusted RPBC rejects missing covariate values", {

  x <- .make_rpbc_assumption_data()

  cov_bad <- x$covariates
  cov_bad$age[1] <- NA_real_

  expect_error(
    .rpbc_fit_matrix(
      dat = x$dat,
      patient = x$patient,
      region = x$region,
      covariates = cov_bad
    ),
    "contains missing values"
  )
})


test_that("covariate-adjusted RPBC rejects non-finite numeric covariates", {

  x <- .make_rpbc_assumption_data()

  cov_bad <- x$covariates
  cov_bad$age[1] <- Inf

  expect_error(
    .rpbc_fit_matrix(
      dat = x$dat,
      patient = x$patient,
      region = x$region,
      covariates = cov_bad
    ),
    "non-finite numeric values"
  )
})


test_that("covariate rows must correspond to patients or cells", {

  x <- .make_rpbc_assumption_data()

  cov_bad <- data.frame(
    age = c(
      40,
      50,
      60
    )
  )

  expect_error(
    .rpbc_fit_matrix(
      dat = x$dat,
      patient = x$patient,
      region = x$region,
      covariates = cov_bad
    ),
    "one row per cell"
  )
})


test_that("patient-level covariate row names must match patient labels", {

  x <- .make_rpbc_assumption_data()

  cov_bad <- x$covariates

  rownames(cov_bad) <- c(
    "A",
    "B",
    "C",
    "D"
  )

  expect_error(
    .rpbc_fit_matrix(
      dat = x$dat,
      patient = x$patient,
      region = x$region,
      covariates = cov_bad
    ),
    "row names"
  )
})


test_that("cell-level covariates must be constant within patient", {

  x <- .make_rpbc_assumption_data()

  cov_cell <- x$covariates[
    match(
      as.character(x$patient),
      rownames(x$covariates)
    ),
    ,
    drop = FALSE
  ]

  rownames(cov_cell) <- NULL

  cov_cell$age[2] <-
    cov_cell$age[2] + 1

  expect_error(
    .rpbc_fit_matrix(
      dat = x$dat,
      patient = x$patient,
      region = x$region,
      covariates = cov_cell
    ),
    "not constant within patient"
  )
})


test_that("constant covariates are rejected because X0 is not full rank", {

  x <- .make_rpbc_assumption_data()

  cov_bad <- data.frame(
    constant = rep(
      1,
      4
    ),
    row.names = x$patient_levels
  )

  expect_error(
    .rpbc_fit_matrix(
      dat = x$dat,
      patient = x$patient,
      region = x$region,
      covariates = cov_bad
    ),
    "not full rank"
  )
})


test_that("linearly dependent covariates are rejected", {

  x <- .make_rpbc_assumption_data()

  cov_bad <- data.frame(
    x1 = c(
      1,
      2,
      3,
      4
    ),
    x2 = c(
      2,
      4,
      6,
      8
    ),
    row.names = x$patient_levels
  )

  expect_error(
    .rpbc_fit_matrix(
      dat = x$dat,
      patient = x$patient,
      region = x$region,
      covariates = cov_bad
    ),
    "not full rank"
  )
})


###############################################################################
# Boundary case: covariates exhaust residual patient degrees of freedom
###############################################################################

test_that("RPBC handles zero residual patient degrees of freedom", {

  x <- .make_rpbc_assumption_data()

  cov_full <- data.frame(
    z1 = c(
      1,
      0,
      0,
      0
    ),
    z2 = c(
      0,
      1,
      0,
      0
    ),
    z3 = c(
      0,
      0,
      1,
      0
    ),
    row.names = x$patient_levels
  )

  expect_warning(
    cov_info <- .rpbc_prepare_covariates(
      covariates = cov_full,
      patient = x$patient
    ),
    "residual patient effect"
  )

  expect_true(
    cov_info$adjusted
  )

  expect_equal(
    cov_info$q,
    3L
  )

  expect_equal(
    cov_info$residual_df,
    0L
  )

  expect_equal(
    dim(cov_info$X),
    c(
      4L,
      3L
    )
  )

  expect_equal(
    qr(cov_info$X0)$rank,
    ncol(cov_info$X0)
  )
})
