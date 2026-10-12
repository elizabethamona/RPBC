test_that("covariate-adjusted RPBC satisfies design and constraint properties", {

  set.seed(2026)

  n_patients <- 6L
  cells_per_region <- 4L
  n_features <- 20L

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
    levels = c("Core", "Edge")
  )

  n_cells <- length(patient)

  covariates <- data.frame(
    age = c(
      45,
      52,
      61,
      48,
      57,
      66
    ),
    sex = factor(
      c(
        "F",
        "M",
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

  age_centered <-
    covariates$age -
    mean(covariates$age)

  sex_numeric <-
    as.numeric(
      covariates$sex == "M"
    )

  patient_index <- match(
    as.character(patient),
    patient_levels
  )

  region_code <- ifelse(
    region == "Edge",
    1,
    -1
  )

  dat <- matrix(
    rnorm(
      n_features * n_cells,
      mean = 0,
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

  for (g in seq_len(n_features)) {

    dat[g, ] <- dat[g, ] +
      5 +
      0.5 * region_code +
      0.04 *
      age_centered[patient_index] +
      0.35 *
      sex_numeric[patient_index]
  }

  fit <- .rpbc_fit_matrix(
    dat = dat,
    patient = patient,
    region = region,
    covariates = covariates,
    core_level = "Core",
    edge_level = "Edge"
  )


  # -------------------------------------------------------------------------
  # Covariate-adjusted mode
  # -------------------------------------------------------------------------

  expect_true(
    fit$covariate_adjusted
  )

  expect_equal(
    fit$q_covariates,
    2L
  )

  expect_equal(
    fit$residual_patient_df,
    3L
  )


  # -------------------------------------------------------------------------
  # Covariate design
  # -------------------------------------------------------------------------

  expect_equal(
    dim(
      fit$covariate_design_patient
    ),
    c(
      6L,
      2L
    )
  )

  expect_equal(
    dim(
      fit$augmented_design_patient
    ),
    c(
      6L,
      3L
    )
  )

  expect_equal(
    dim(
      fit$beta_covariate_hat
    ),
    c(
      2L,
      n_features
    )
  )

  expect_lt(
    max(
      abs(
        colMeans(
          fit$covariate_design_patient
        )
      )
    ),
    1e-10
  )


  # -------------------------------------------------------------------------
  # Null-space construction
  # -------------------------------------------------------------------------

  expect_lt(
    fit$nullspace_error,
    1e-8
  )

  expect_lt(
    max(
      abs(
        t(
          fit$augmented_design_patient
        ) %*%
          fit$residual_patient_basis
      )
    ),
    1e-8
  )


  # -------------------------------------------------------------------------
  # Initial identifying constraints
  # -------------------------------------------------------------------------

  expect_lt(
    fit$gamma_hat_constraint_error,
    1e-7
  )

  expect_lt(
    fit$kappa_hat_constraint_error,
    1e-7
  )

  expect_lt(
    max(
      abs(
        t(
          fit$augmented_design_patient
        ) %*%
          fit$gamma_hat
      )
    ),
    1e-7
  )

  expect_lt(
    max(
      abs(
        colSums(
          fit$kappa_hat
        )
      )
    ),
    1e-7
  )


  # -------------------------------------------------------------------------
  # Posterior identifying constraints
  # -------------------------------------------------------------------------

  expect_lt(
    fit$gamma_star_constraint_error,
    1e-7
  )

  expect_lt(
    fit$kappa_star_constraint_error,
    1e-7
  )

  expect_lt(
    max(
      abs(
        t(
          fit$augmented_design_patient
        ) %*%
          fit$gamma_star
      )
    ),
    1e-7
  )

  expect_lt(
    max(
      abs(
        colSums(
          fit$kappa_star
        )
      )
    ),
    1e-7
  )


  # -------------------------------------------------------------------------
  # Model rank
  # -------------------------------------------------------------------------

  expect_equal(
    fit$design_rank,
    fit$design_columns
  )

  expect_equal(
    fit$constraint_rank,
    fit$constraint_rows
  )


  # -------------------------------------------------------------------------
  # Harmonized output
  # -------------------------------------------------------------------------

  expect_equal(
    dim(fit$data),
    dim(dat)
  )

  expect_identical(
    dimnames(fit$data),
    dimnames(dat)
  )

  expect_true(
    all(
      is.finite(
        fit$data
      )
    )
  )
})
