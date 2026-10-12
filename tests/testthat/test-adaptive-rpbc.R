# ============================================================================
# Automated tests for Adaptive RPBC
# ============================================================================


# ----------------------------------------------------------------------------
# Helper: reproducible Adaptive RPBC toy data
# ----------------------------------------------------------------------------

make_adaptive_test_data <- function(
    group_n = rep(5L, 8L),
    seed = 20261007
) {

  set.seed(seed)

  group_info <- data.frame(
    patient = c(
      "P1", "P1",
      "P2", "P2",
      "P3", "P3",
      "P4", "P4"
    ),
    region = c(
      "Core", "Edge",
      "Core", "Edge",
      "Core", "Edge",
      "Core", "Edge"
    ),
    n = group_n
  )

  patient <- rep(
    group_info$patient,
    times = group_info$n
  )

  region <- rep(
    group_info$region,
    times = group_info$n
  )

  n_features <- 6L

  alpha <- c(
    5.0, 6.0, 7.0,
    8.0, 9.0, 10.0
  )

  beta_region <- c(
    0.60, 0.40, 0.75,
    0.50, 0.65, 0.35
  )

  gamma_true <- rbind(
    P1 = c( 1.00,  0.80,  1.20,  0.70,  1.10,  0.90),
    P2 = c(-0.70, -0.60, -0.90, -0.50, -0.80, -0.60),
    P3 = c( 0.50,  0.40,  0.60,  0.30,  0.50,  0.40),
    P4 = c(-0.80, -0.60, -0.90, -0.50, -0.80, -0.70)
  )

  kappa_true <- rbind(
    P1 = c( 0.30,  0.20,  0.25,  0.30,  0.20,  0.25),
    P2 = c(-0.20, -0.15, -0.20, -0.20, -0.15, -0.20),
    P3 = c( 0.15,  0.10,  0.15,  0.10,  0.15,  0.10),
    P4 = c(-0.25, -0.15, -0.20, -0.20, -0.20, -0.15)
  )

  group_sd <- c(
    P1_Core = 0.40,
    P1_Edge = 1.40,
    P2_Core = 0.65,
    P2_Edge = 1.20,
    P3_Core = 1.50,
    P3_Edge = 0.50,
    P4_Core = 0.75,
    P4_Edge = 1.60
  )

  n_cells <- length(patient)

  dat <- matrix(
    NA_real_,
    nrow = n_features,
    ncol = n_cells
  )

  rownames(dat) <- paste0(
    "Feature",
    seq_len(n_features)
  )

  cell_counter <- 1L

  for (h in seq_len(nrow(group_info))) {

    p <- group_info$patient[h]
    r <- group_info$region[h]
    nh <- group_info$n[h]

    R_code <- ifelse(
      r == "Edge",
      1,
      -1
    )

    sd_h <- group_sd[
      paste0(p, "_", r)
    ]

    idx <- cell_counter:(
      cell_counter + nh - 1L
    )

    for (g in seq_len(n_features)) {

      mu <- alpha[g] +
        beta_region[g] * R_code +
        gamma_true[p, g] +
        kappa_true[p, g] * R_code

      dat[g, idx] <- rnorm(
        nh,
        mean = mu,
        sd = sd_h
      )
    }

    cell_counter <- cell_counter + nh
  }

  list(
    dat = dat,
    patient = patient,
    region = region
  )
}


# ----------------------------------------------------------------------------
# Helper: patient-balanced Edge-Core contrast
# ----------------------------------------------------------------------------

patient_balanced_contrast <- function(
    dat,
    patient,
    region
) {

  patients <- unique(patient)

  vapply(
    seq_len(nrow(dat)),
    function(g) {

      patient_contrasts <- vapply(
        patients,
        function(p) {

          mean(
            dat[
              g,
              patient == p &
                region == "Edge"
            ]
          ) -
            mean(
              dat[
                g,
                patient == p &
                  region == "Core"
              ]
            )
        },
        numeric(1)
      )

      mean(patient_contrasts)
    },
    numeric(1)
  )
}


# ----------------------------------------------------------------------------
# Helper: patient-region feature SDs
# ----------------------------------------------------------------------------

patient_region_sd <- function(
    dat,
    patient,
    region
) {

  patient_region <- interaction(
    factor(patient),
    factor(region),
    drop = TRUE,
    sep = "__"
  )

  groups <- levels(patient_region)

  out <- sapply(
    groups,
    function(gr) {

      idx <- patient_region == gr

      apply(
        dat[, idx, drop = FALSE],
        1,
        sd
      )
    }
  )

  t(out)
}


# ============================================================================
# Test 1
# ============================================================================

test_that(
  "Adaptive RPBC returns finite output and converges",
  {

    x <- make_adaptive_test_data()

    fit <- .rpbc_fit_matrix_adaptive(
      dat = x$dat,
      patient = x$patient,
      region = x$region
    )

    expect_equal(
      dim(fit$data),
      dim(x$dat)
    )

    expect_true(
      all(is.finite(fit$data))
    )

    expect_true(
      all(fit$converged)
    )
  }
)


# ============================================================================
# Test 2
# ============================================================================

test_that(
  "Adaptive RPBC enforces posterior identifiability constraints",
  {

    x <- make_adaptive_test_data()

    fit <- .rpbc_fit_matrix_adaptive(
      dat = x$dat,
      patient = x$patient,
      region = x$region
    )

    expect_lt(
      max(
        abs(
          colSums(fit$gamma_star)
        )
      ),
      1e-10
    )

    expect_lt(
      max(
        abs(
          colSums(fit$kappa_star)
        )
      ),
      1e-10
    )
  }
)


# ============================================================================
# Test 3
# Equation (15)
# ============================================================================

test_that(
  "Adaptive RPBC computes group-size weights correctly",
  {

    x <- make_adaptive_test_data()

    fit <- .rpbc_fit_matrix_adaptive(
      dat = x$dat,
      patient = x$patient,
      region = x$region
    )

    expected_weight <-
      fit$n_group /
      (
        fit$n_group +
          fit$delta_adaptive_k
      )

    expect_equal(
      as.numeric(fit$n_group),
      rep(5, 8)
    )

    expect_equal(
      unname(
        fit$delta_shrinkage_weight
      ),
      matrix(
        rep(
          expected_weight,
          ncol(fit$delta_shrinkage_weight)
        ),
        nrow = length(expected_weight)
      ),
      tolerance = 1e-12
    )
  }
)


# ============================================================================
# Test 4
# Equations (16)-(17)
# ============================================================================

test_that(
  "Adaptive RPBC stabilization matches Equations 16 and 17",
  {

    x <- make_adaptive_test_data()

    fit <- .rpbc_fit_matrix_adaptive(
      dat = x$dat,
      patient = x$patient,
      region = x$region
    )

    expected_delta <- exp(
      fit$delta_shrinkage_weight *
        log(fit$delta_star_raw)
    )

    expected_delta <- pmin(
      pmax(
        expected_delta,
        fit$delta_bounds[1]
      ),
      fit$delta_bounds[2]
    )

    expect_equal(
      fit$delta_star,
      expected_delta,
      tolerance = 1e-12
    )

    expect_true(
      all(
        fit$delta_star >=
          fit$delta_bounds[1]
      )
    )

    expect_true(
      all(
        fit$delta_star <=
          fit$delta_bounds[2]
      )
    )
  }
)


# ============================================================================
# Test 5
# Equation (18)
# ============================================================================

test_that(
  "Adaptive RPBC output reproduces Equation 18",
  {

    x <- make_adaptive_test_data()

    fit <- .rpbc_fit_matrix_adaptive(
      dat = x$dat,
      patient = x$patient,
      region = x$region
    )

    R_check <- ifelse(
      x$region == "Edge",
      1,
      -1
    )

    preserved_mean <- outer(
      fit$alpha_hat,
      rep(
        1,
        length(x$region)
      )
    ) +
      outer(
        fit$beta_region_hat,
        R_check
      )

    s_check <- sweep(
      x$dat - preserved_mean,
      MARGIN = 1,
      STATS = sqrt(fit$var_pooled),
      FUN = "/"
    )

    patient_region <- interaction(
      factor(x$patient),
      factor(x$region),
      drop = TRUE,
      sep = "__"
    )

    group_indices <- split(
      seq_along(x$patient),
      patient_region
    )

    expected_standardized <- s_check

    for (h in seq_along(fit$group_names)) {

      group_name <- fit$group_names[h]

      idx <- group_indices[[group_name]]

      p_index <- match(
        fit$group_patient[h],
        rownames(fit$gamma_star)
      )

      R_h <- ifelse(
        fit$group_region[h] == "Edge",
        1,
        -1
      )

      nuisance_mean <-
        fit$gamma_star[p_index, ] +
        fit$kappa_star[p_index, ] *
        R_h

      expected_group <- sweep(
        s_check[
          ,
          idx,
          drop = FALSE
        ],
        MARGIN = 1,
        STATS = nuisance_mean,
        FUN = "-"
      )

      expected_group <- sweep(
        expected_group,
        MARGIN = 1,
        STATS = sqrt(
          fit$delta_star[h, ]
        ),
        FUN = "/"
      )

      expected_standardized[
        ,
        idx
      ] <- expected_group
    }

    expected_data <- sweep(
      expected_standardized,
      MARGIN = 1,
      STATS = sqrt(fit$var_pooled),
      FUN = "*"
    ) +
      preserved_mean

    expect_equal(
      fit$data,
      expected_data,
      tolerance = 1e-12
    )
  }
)


# ============================================================================
# Test 6
# Preservation target
# ============================================================================

test_that(
  "Adaptive RPBC closely preserves the patient-balanced regional contrast",
  {

    x <- make_adaptive_test_data()

    fit <- .rpbc_fit_matrix_adaptive(
      dat = x$dat,
      patient = x$patient,
      region = x$region
    )

    contrast_before <- patient_balanced_contrast(
      x$dat,
      x$patient,
      x$region
    )

    contrast_after <- patient_balanced_contrast(
      fit$data,
      x$patient,
      x$region
    )

    # Fixture-specific tolerance:
    # this is not a universal theoretical tolerance.
    expect_lt(
      max(
        abs(
          contrast_after -
            contrast_before
        )
      ),
      0.10
    )
  }
)


# ============================================================================
# Test 7
# Patient-region scale harmonization
# ============================================================================

test_that(
  "Adaptive RPBC reduces patient-region scale heterogeneity",
  {

    x <- make_adaptive_test_data()

    fit <- .rpbc_fit_matrix_adaptive(
      dat = x$dat,
      patient = x$patient,
      region = x$region
    )

    sd_before <- patient_region_sd(
      x$dat,
      x$patient,
      x$region
    )

    sd_after <- patient_region_sd(
      fit$data,
      x$patient,
      x$region
    )

    heterogeneity_before <- apply(
      log(sd_before),
      2,
      sd
    )

    heterogeneity_after <- apply(
      log(sd_after),
      2,
      sd
    )

    expect_true(
      all(
        heterogeneity_after <
          heterogeneity_before
      )
    )

    expect_lt(
      mean(heterogeneity_after),
      mean(heterogeneity_before)
    )
  }
)


# ============================================================================
# Test 8
# Unequal group sizes
# ============================================================================

test_that(
  "Adaptive stabilization weakens as patient-region group size increases",
  {

    x <- make_adaptive_test_data(
      group_n = c(
        3, 5,
        8, 12,
        15, 20,
        25, 30
      ),
      seed = 20261008
    )

    fit <- .rpbc_fit_matrix_adaptive(
      dat = x$dat,
      patient = x$patient,
      region = x$region
    )

    expected_weight <-
      fit$n_group /
      (
        fit$n_group +
          fit$delta_adaptive_k
      )

    observed_weight <-
      fit$delta_shrinkage_weight[
        ,
        1
      ]

    expect_equal(
      as.numeric(observed_weight),
      as.numeric(expected_weight),
      tolerance = 1e-12
    )

    order_n <- order(
      fit$n_group
    )

    expect_true(
      all(
        diff(
          observed_weight[
            order_n
          ]
        ) > 0
      )
    )
  }
)


# ============================================================================
# Test 9
# Adaptive-specific argument validation
# ============================================================================

test_that(
  "Adaptive RPBC rejects invalid stabilization settings",
  {

    x <- make_adaptive_test_data()

    expect_error(
      .rpbc_fit_matrix_adaptive(
        dat = x$dat,
        patient = x$patient,
        region = x$region,
        delta_adaptive_k = 0
      ),
      "positive finite"
    )

    expect_error(
      .rpbc_fit_matrix_adaptive(
        dat = x$dat,
        patient = x$patient,
        region = x$region,
        delta_bounds = c(2, 0.5)
      ),
      "lower < upper"
    )
  }
)
