test_that("ordinary RPBC returns a valid harmonized matrix", {

  set.seed(2026)

  patient <- rep(
    paste0("P", 1:4),
    each = 10
  )

  region <- rep(
    rep(c("Core", "Edge"), each = 5),
    times = 4
  )

  dat <- matrix(
    rnorm(6 * 40, mean = 0, sd = 0.5),
    nrow = 6,
    ncol = 40
  )

  rownames(dat) <- paste0("Feature", 1:6)
  colnames(dat) <- paste0("Cell", 1:40)

  dat <- sweep(
    dat,
    MARGIN = 1,
    STATS = seq(1, 6),
    FUN = "+"
  )

  dat <- sweep(
    dat,
    MARGIN = 2,
    STATS = ifelse(region == "Edge", 0.8, -0.8),
    FUN = "+"
  )

  patient_shift <- c(
    P1 = -0.7,
    P2 = 0.4,
    P3 = 1.0,
    P4 = -0.7
  )

  dat <- sweep(
    dat,
    MARGIN = 2,
    STATS = patient_shift[patient],
    FUN = "+"
  )

  fit <- .rpbc_fit_matrix(
    dat = dat,
    patient = patient,
    region = region
  )

  expect_equal(
    dim(fit$data),
    dim(dat)
  )

  expect_true(
    all(is.finite(fit$data))
  )

  expect_true(
    all(fit$converged)
  )
})

test_that("ordinary RPBC enforces identifiability constraints", {

  set.seed(2026)

  patient <- rep(
    paste0("P", 1:4),
    each = 10
  )

  region <- rep(
    rep(c("Core", "Edge"), each = 5),
    times = 4
  )

  dat <- matrix(
    rnorm(6 * 40, mean = 0, sd = 0.5),
    nrow = 6,
    ncol = 40
  )

  rownames(dat) <- paste0("Feature", 1:6)
  colnames(dat) <- paste0("Cell", 1:40)

  dat <- sweep(
    dat,
    MARGIN = 1,
    STATS = seq(1, 6),
    FUN = "+"
  )

  dat <- sweep(
    dat,
    MARGIN = 2,
    STATS = ifelse(region == "Edge", 0.8, -0.8),
    FUN = "+"
  )

  patient_shift <- c(
    P1 = -0.7,
    P2 = 0.4,
    P3 = 1.0,
    P4 = -0.7
  )

  dat <- sweep(
    dat,
    MARGIN = 2,
    STATS = patient_shift[patient],
    FUN = "+"
  )

  fit <- .rpbc_fit_matrix(
    dat = dat,
    patient = patient,
    region = region
  )

  expect_equal(
    unname(colSums(fit$gamma_star)),
    rep(0, nrow(dat)),
    tolerance = 1e-10
  )

  expect_equal(
    unname(colSums(fit$kappa_star)),
    rep(0, nrow(dat)),
    tolerance = 1e-10
  )
})

test_that("ordinary RPBC reduces patient-wide nuisance location variation", {

  set.seed(2026)

  patient <- rep(
    paste0("P", 1:4),
    each = 10
  )

  region <- rep(
    rep(c("Core", "Edge"), each = 5),
    times = 4
  )

  dat <- matrix(
    rnorm(6 * 40, mean = 0, sd = 0.5),
    nrow = 6,
    ncol = 40
  )

  rownames(dat) <- paste0("Feature", 1:6)
  colnames(dat) <- paste0("Cell", 1:40)

  dat <- sweep(
    dat,
    MARGIN = 1,
    STATS = seq(1, 6),
    FUN = "+"
  )

  dat <- sweep(
    dat,
    MARGIN = 2,
    STATS = ifelse(region == "Edge", 0.8, -0.8),
    FUN = "+"
  )

  patient_shift <- c(
    P1 = -0.7,
    P2 = 0.4,
    P3 = 1.0,
    P4 = -0.7
  )

  dat <- sweep(
    dat,
    MARGIN = 2,
    STATS = patient_shift[patient],
    FUN = "+"
  )

  fit <- .rpbc_fit_matrix(
    dat = dat,
    patient = patient,
    region = region
  )

  patient_means_before <- sapply(
    unique(patient),
    function(p) {
      rowMeans(
        dat[, patient == p, drop = FALSE]
      )
    }
  )

  patient_means_after <- sapply(
    unique(patient),
    function(p) {
      rowMeans(
        fit$data[, patient == p, drop = FALSE]
      )
    }
  )

  sd_before <- apply(
    patient_means_before,
    1,
    sd
  )

  sd_after <- apply(
    patient_means_after,
    1,
    sd
  )

  expect_true(
    all(sd_after < sd_before)
  )
})

test_that("ordinary RPBC reduces patient-region scale heterogeneity", {

  set.seed(2027)

  patient <- rep(
    paste0("P", 1:4),
    each = 20
  )

  region <- rep(
    rep(c("Core", "Edge"), each = 10),
    times = 4
  )

  group_sd <- c(
    "P1_Core" = 0.30,
    "P1_Edge" = 1.20,
    "P2_Core" = 0.60,
    "P2_Edge" = 1.50,
    "P3_Core" = 1.30,
    "P3_Edge" = 0.40,
    "P4_Core" = 0.80,
    "P4_Edge" = 1.10
  )

  patient_shift <- c(
    P1 = -0.8,
    P2 = 0.3,
    P3 = 1.0,
    P4 = -0.5
  )

  dat <- matrix(
    0,
    nrow = 6,
    ncol = length(patient)
  )

  rownames(dat) <- paste0("Feature", 1:6)
  colnames(dat) <- paste0(
    "Cell",
    seq_len(ncol(dat))
  )

  for (j in seq_len(ncol(dat))) {

    p <- patient[j]
    r <- region[j]

    sd_j <- group_sd[
      paste0(p, "_", r)
    ]

    region_effect <- ifelse(
      r == "Edge",
      0.8,
      -0.8
    )

    dat[, j] <- seq(1, 6) +
      region_effect +
      patient_shift[p] +
      rnorm(6, mean = 0, sd = sd_j)
  }

  fit <- .rpbc_fit_matrix(
    dat = dat,
    patient = patient,
    region = region
  )

  expect_true(
    all(fit$converged)
  )

  patient_region <- interaction(
    patient,
    region,
    sep = "_",
    drop = TRUE
  )

  group_variance <- function(x) {

    sapply(
      levels(patient_region),
      function(gr) {
        matrixStats::rowVars(
          x[, patient_region == gr, drop = FALSE]
        )
      }
    )
  }

  variance_before <- group_variance(dat)
  variance_after <- group_variance(fit$data)

  group_sd_before <- sqrt(
    colMeans(variance_before)
  )

  group_sd_after <- sqrt(
    colMeans(variance_after)
  )

  expect_lt(
    sd(group_sd_after),
    sd(group_sd_before)
  )
})

test_that("ordinary RPBC preserves the population-average regional contrast", {

  set.seed(2027)

  patient <- rep(
    paste0("P", 1:4),
    each = 20
  )

  region <- rep(
    rep(c("Core", "Edge"), each = 10),
    times = 4
  )

  group_sd <- c(
    "P1_Core" = 0.30,
    "P1_Edge" = 1.20,
    "P2_Core" = 0.60,
    "P2_Edge" = 1.50,
    "P3_Core" = 1.30,
    "P3_Edge" = 0.40,
    "P4_Core" = 0.80,
    "P4_Edge" = 1.10
  )

  patient_shift <- c(
    P1 = -0.8,
    P2 = 0.3,
    P3 = 1.0,
    P4 = -0.5
  )

  dat <- matrix(
    0,
    nrow = 6,
    ncol = length(patient)
  )

  rownames(dat) <- paste0("Feature", 1:6)
  colnames(dat) <- paste0(
    "Cell",
    seq_len(ncol(dat))
  )

  for (j in seq_len(ncol(dat))) {

    p <- patient[j]
    r <- region[j]

    sd_j <- group_sd[
      paste0(p, "_", r)
    ]

    region_effect <- ifelse(
      r == "Edge",
      0.8,
      -0.8
    )

    dat[, j] <- seq(1, 6) +
      region_effect +
      patient_shift[p] +
      rnorm(6, mean = 0, sd = sd_j)
  }

  fit <- .rpbc_fit_matrix(
    dat = dat,
    patient = patient,
    region = region
  )

  patient_balanced_contrast <- function(x) {

    per_patient <- sapply(
      unique(patient),
      function(p) {

        rowMeans(
          x[, patient == p & region == "Edge", drop = FALSE]
        ) -
          rowMeans(
            x[, patient == p & region == "Core", drop = FALSE]
          )
      }
    )

    rowMeans(per_patient)
  }

  contrast_before <- patient_balanced_contrast(dat)
  contrast_after <- patient_balanced_contrast(fit$data)

  expect_lt(
    max(abs(contrast_after - contrast_before)),
    0.10
  )
})

test_that("ordinary RPBC rejects invalid inputs", {

  set.seed(2028)

  patient <- rep(
    paste0("P", 1:4),
    each = 10
  )

  region <- rep(
    rep(c("Core", "Edge"), each = 5),
    times = 4
  )

  dat <- matrix(
    rnorm(6 * 40),
    nrow = 6,
    ncol = 40
  )

  # -------------------------------------------------------------------------
  # Missing value in the data
  # -------------------------------------------------------------------------

  dat_na <- dat
  dat_na[1, 1] <- NA

  expect_error(
    .rpbc_fit_matrix(
      dat = dat_na,
      patient = patient,
      region = region
    ),
    "contains NA, NaN, or Inf"
  )

  # -------------------------------------------------------------------------
  # Only one patient
  # -------------------------------------------------------------------------

  one_patient <- patient == "P1"

  expect_error(
    .rpbc_fit_matrix(
      dat = dat[, one_patient, drop = FALSE],
      patient = patient[one_patient],
      region = region[one_patient]
    ),
    "at least two patients"
  )

  # -------------------------------------------------------------------------
  # Incorrect region labels
  # -------------------------------------------------------------------------

  region_bad <- ifelse(
    region == "Edge",
    "Margin",
    "Core"
  )

  expect_error(
    .rpbc_fit_matrix(
      dat = dat,
      patient = patient,
      region = region_bad
    ),
    "must contain exactly the two specified levels"
  )

  # -------------------------------------------------------------------------
  # Unpaired patient
  # -------------------------------------------------------------------------

  keep_paired <- !(
    patient == "P4" &
      region == "Edge"
  )

  expect_error(
    .rpbc_fit_matrix(
      dat = dat[, keep_paired, drop = FALSE],
      patient = patient[keep_paired],
      region = region[keep_paired]
    ),
    "Every patient must have observations in both"
  )

  # -------------------------------------------------------------------------
  # Fewer than two cells in one patient-region group
  # -------------------------------------------------------------------------

  p1_core <- which(
    patient == "P1" &
      region == "Core"
  )

  keep_cells <- rep(TRUE, length(patient))
  keep_cells[p1_core[-1]] <- FALSE

  expect_error(
    .rpbc_fit_matrix(
      dat = dat[, keep_cells, drop = FALSE],
      patient = patient[keep_cells],
      region = region[keep_cells]
    ),
    "must contain at least two cells"
  )

  # -------------------------------------------------------------------------
  # Fewer than two features
  # -------------------------------------------------------------------------

  expect_error(
    .rpbc_fit_matrix(
      dat = dat[1, , drop = FALSE],
      patient = patient,
      region = region
    ),
    "at least two features"
  )
})
