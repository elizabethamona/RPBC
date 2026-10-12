# Internal RPBC statistical engine
#
# Supports ordinary and covariate-adjusted RPBC.
#
# This function operates on an expression matrix and paired patient/region
# vectors. It contains the statistical RPBC algorithm and is intentionally
# independent of SingleCellExperiment and SpatialExperiment objects.

.rpbc_fit_matrix <- function(
    dat,
    patient,
    region,
    covariates = NULL,
    core_level = "Core",
    edge_level = "Edge",
    conv = 1e-4,
    max_iter = 100L,
    eps = 1e-8
) {
  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------

  dat <- as.matrix(dat)

  if (!is.numeric(dat)) {
    stop("`dat` must be a numeric matrix.")
  }

  if (nrow(dat) < 2L || ncol(dat) < 1L) {
    stop("`dat` must contain at least two features and one cell.")
  }

  if (ncol(dat) != length(patient) || ncol(dat) != length(region)) {
    stop("`patient` and `region` must have one value for every column of `dat`.")
  }

  if (any(!is.finite(dat))) {
    stop("`dat` contains NA, NaN, or Inf values.")
  }

  if (anyNA(patient) || anyNA(region)) {
    stop("`patient` and `region` cannot contain missing values.")
  }

  if (!is.character(core_level) || length(core_level) != 1L ||
      !is.character(edge_level) || length(edge_level) != 1L ||
      identical(core_level, edge_level)) {
    stop("`core_level` and `edge_level` must be two different region labels.")
  }

  patient <- droplevels(factor(patient))
  region <- droplevels(factor(region))

  if (nlevels(patient) < 2L) {
    stop("RPBC requires at least two patients.")
  }

  if (!setequal(levels(region), c(core_level, edge_level))) {
    stop(
      "`region` must contain exactly the two specified levels: `",
      core_level, "` and `", edge_level, "`."
    )
  }

  if (!is.numeric(conv) || length(conv) != 1L ||
      !is.finite(conv) || conv <= 0) {
    stop("`conv` must be a single positive finite value.")
  }

  if (!is.numeric(max_iter) || length(max_iter) != 1L ||
      !is.finite(max_iter) || max_iter < 1 ||
      max_iter != as.integer(max_iter)) {
    stop("`max_iter` must be a positive whole number.")
  }

  max_iter <- as.integer(max_iter)

  if (!is.numeric(eps) || length(eps) != 1L ||
      !is.finite(eps) || eps <= 0) {
    stop("`eps` must be a single positive finite value.")
  }

  # RPBC requires a paired design: every patient must contribute observations
  # from both prespecified regions.
  patient_region_table <- table(patient, region)

  if (any(patient_region_table == 0L)) {
    stop(
      "Every patient must have observations in both `",
      core_level, "` and `", edge_level, "`."
    )
  }

  n_patients <- nlevels(patient)
  n_features <- nrow(dat)
  n_cells <- ncol(dat)

  # ---------------------------------------------------------------------------
  # Region coding
  #
  # Reference region = -1
  # Contrast region  = +1
  # ---------------------------------------------------------------------------

  R <- ifelse(region == edge_level, 1, -1)

  # ---------------------------------------------------------------------------
  # Patient effect coding
  #
  # Sum-to-zero coding:
  # P patients -> P - 1 independent contrast columns
  # ---------------------------------------------------------------------------

  patient_contrast <- stats::contr.sum(n_patients)
  rownames(patient_contrast) <- levels(patient)

  P <- patient_contrast[
    as.character(patient),
    ,
    drop = FALSE
  ]
  # ---------------------------------------------------------------------------
  # Prepare patient-level covariate design
  #
  # covariates = NULL:
  #   ordinary RPBC
  #
  # covariates supplied:
  #   covariate-adjusted RPBC
  # ---------------------------------------------------------------------------

  covariate_info <- .rpbc_prepare_covariates(
    covariates = covariates,
    patient = patient
  )

  covariate_adjusted <- covariate_info$adjusted
  X_patient <- covariate_info$X
  X0_patient <- covariate_info$X0
  q_covariates <- covariate_info$q
  residual_patient_df <- covariate_info$residual_df

  # ---------------------------------------------------------------------------
  # Patient x region interaction
  #
  # Represents patient-specific deviations from the population-average
  # regional contrast.
  # ---------------------------------------------------------------------------

  PR <- sweep(
    P,
    MARGIN = 1,
    STATS = R,
    FUN = "*"
  )

  # ---------------------------------------------------------------------------
  # Full location model
  #
  # Ordinary RPBC:
  #   intercept + region + patient + patient x region
  #
  # Covariate-adjusted RPBC:
  #   intercept + region + covariates
  #   + residual patient effect + patient x region
  # ---------------------------------------------------------------------------

  if (!covariate_adjusted) {

    # Preserve the original ordinary-RPBC parameterization exactly.
    design_full <- cbind(
      intercept = 1,
      region = R,
      P,
      PR
    )

    protected_rows <- 1:2
    covariate_rows <- integer(0)

    residual_patient_rows <-
      2L + seq_len(n_patients - 1L)

    interaction_rows <-
      2L +
      (n_patients - 1L) +
      seq_len(n_patients - 1L)

    residual_patient_basis <-
      patient_contrast

  } else {

    # -------------------------------------------------------------------------
    # Null-space basis for residual patient effects
    #
    # Columns of residual_patient_basis span Null(X0^T).
    # -------------------------------------------------------------------------

    qr_X0 <- qr(
      X0_patient
    )

    Q_X0 <- qr.Q(
      qr_X0,
      complete = TRUE
    )

    if (residual_patient_df > 0L) {

      residual_patient_basis <- Q_X0[
        ,
        (ncol(X0_patient) + 1L):n_patients,
        drop = FALSE
      ]

      colnames(residual_patient_basis) <- paste0(
        "residual_patient_basis_",
        seq_len(residual_patient_df)
      )

    } else {

      residual_patient_basis <- matrix(
        numeric(0),
        nrow = n_patients,
        ncol = 0,
        dimnames = list(
          levels(patient),
          NULL
        )
      )
    }

    rownames(residual_patient_basis) <-
      levels(patient)

    X_cell <- X_patient[
      as.character(patient),
      ,
      drop = FALSE
    ]

    residual_patient_cell <-
      residual_patient_basis[
        as.character(patient),
        ,
        drop = FALSE
      ]

    X_cell_design <- X_cell

    colnames(X_cell_design) <- paste0(
      "covariate__",
      colnames(X_patient)
    )

    residual_patient_cell_design <-
      residual_patient_cell

    if (ncol(residual_patient_cell_design) > 0L) {

      colnames(residual_patient_cell_design) <- paste0(
        "patient_residual__",
        seq_len(
          ncol(residual_patient_cell_design)
        )
      )
    }

    PR_design <- PR

    colnames(PR_design) <- paste0(
      "patient_region__",
      seq_len(ncol(PR_design))
    )

    design_full <- cbind(
      intercept = 1,
      region = R,
      X_cell_design,
      residual_patient_cell_design,
      PR_design
    )

    covariate_rows <-
      2L + seq_len(q_covariates)

    protected_rows <- c(
      1L,
      2L,
      covariate_rows
    )

    residual_patient_rows <-
      if (residual_patient_df > 0L) {

        2L +
          q_covariates +
          seq_len(residual_patient_df)

      } else {

        integer(0)
      }

    interaction_rows <-
      2L +
      q_covariates +
      residual_patient_df +
      seq_len(n_patients - 1L)
  }


  # ---------------------------------------------------------------------------
  # Full-rank model-design check
  # ---------------------------------------------------------------------------

  design_qr <- qr(
    design_full
  )

  if (design_qr$rank != ncol(design_full)) {
    stop(
      "RPBC design matrix is not full rank."
    )
  }


  # ---------------------------------------------------------------------------
  # Covariate-design orthogonality check
  #
  # Verify X0^T Q_perp = 0.
  # ---------------------------------------------------------------------------

  if (
    covariate_adjusted &&
    residual_patient_df > 0L
  ) {

    nullspace_error <- max(
      abs(
        t(X0_patient) %*%
          residual_patient_basis
      )
    )

    if (
      !is.finite(nullspace_error) ||
      nullspace_error > 1e-8
    ) {
      stop(
        "Failed to construct the residual patient-effect null-space basis."
      )
    }

  } else {

    nullspace_error <- 0
  }

  # ---------------------------------------------------------------------------
  # Initial least-squares estimates
  #
  # Fit the complete location model simultaneously for all features.
  # ---------------------------------------------------------------------------

  B_hat <- qr.solve(
    design_full,
    t(dat)
  )

  fitted_full <- design_full %*% B_hat
  residual_full <- t(dat) - fitted_full

  # ---------------------------------------------------------------------------
  # Pooled feature variance
  #
  # var_pooled[g] estimates the residual variance for feature g.
  # ---------------------------------------------------------------------------

  var_pooled <- colMeans(residual_full^2)

  if (any(!is.finite(var_pooled)) || any(var_pooled <= 0)) {
    stop("At least one feature has non-positive pooled variance.")
  }

  # ---------------------------------------------------------------------------
  # Protected location component
  #
  # Ordinary RPBC:
  #   alpha_g + beta_Rg R_r
  #
  # Covariate-adjusted RPBC:
  #   alpha_g + beta_Rg R_r + x_p^T beta_Xg
  # ---------------------------------------------------------------------------

  preserved_mean <- design_full[
    ,
    protected_rows,
    drop = FALSE
  ] %*%
    B_hat[
      protected_rows,
      ,
      drop = FALSE
    ]

  preserved_mean <- t(
    preserved_mean
  )


  # ---------------------------------------------------------------------------
  # Standardize data after removing the protected location component
  # ---------------------------------------------------------------------------

  s_data <- sweep(
    dat - preserved_mean,
    MARGIN = 1,
    STATS = sqrt(var_pooled),
    FUN = "/"
  )


  # ---------------------------------------------------------------------------
  # Initial residual patient effects
  # ---------------------------------------------------------------------------

  if (!covariate_adjusted) {

    # Preserve the original ordinary-RPBC reconstruction exactly.
    gamma_contrast <- B_hat[
      residual_patient_rows,
      ,
      drop = FALSE
    ]

    gamma_contrast <- sweep(
      gamma_contrast,
      MARGIN = 2,
      STATS = sqrt(var_pooled),
      FUN = "/"
    )

    gamma_hat <- rbind(
      gamma_contrast,
      -colSums(gamma_contrast)
    )

  } else if (residual_patient_df > 0L) {

    gamma_basis_coeff <- B_hat[
      residual_patient_rows,
      ,
      drop = FALSE
    ]

    gamma_basis_coeff <- sweep(
      gamma_basis_coeff,
      MARGIN = 2,
      STATS = sqrt(var_pooled),
      FUN = "/"
    )

    gamma_hat <- residual_patient_basis %*%
      gamma_basis_coeff

  } else {

    gamma_hat <- matrix(
      0,
      nrow = n_patients,
      ncol = n_features
    )
  }


  # ---------------------------------------------------------------------------
  # Initial patient x region interaction effects
  #
  # The original RPBC sum-to-zero interaction coding is retained
  # in both ordinary and covariate-adjusted RPBC.
  # ---------------------------------------------------------------------------

  kappa_contrast <- B_hat[
    interaction_rows,
    ,
    drop = FALSE
  ]

  kappa_contrast <- sweep(
    kappa_contrast,
    MARGIN = 2,
    STATS = sqrt(var_pooled),
    FUN = "/"
  )

  kappa_hat <- rbind(
    kappa_contrast,
    -colSums(kappa_contrast)
  )

  rownames(gamma_hat) <-
    levels(patient)

  rownames(kappa_hat) <-
    levels(patient)


  # ---------------------------------------------------------------------------
  # Covariate coefficients on original feature scale
  # ---------------------------------------------------------------------------

  if (covariate_adjusted) {

    beta_covariate_hat <- B_hat[
      covariate_rows,
      ,
      drop = FALSE
    ]

    rownames(beta_covariate_hat) <-
      colnames(X_patient)

  } else {

    beta_covariate_hat <- matrix(
      numeric(0),
      nrow = 0,
      ncol = n_features,
      dimnames = list(
        NULL,
        rownames(dat)
      )
    )
  }


  # ---------------------------------------------------------------------------
  # Verify initial identifying constraints
  #
  # Residual patient effect:
  #   X0^T gamma_g = 0
  #
  # Patient x region interaction:
  #   sum_p kappa_pg = 0
  # ---------------------------------------------------------------------------

  gamma_hat_constraint_error <- max(
    abs(
      t(X0_patient) %*%
        gamma_hat
    )
  )

  kappa_hat_constraint_error <- max(
    abs(
      colSums(kappa_hat)
    )
  )

  if (
    !is.finite(gamma_hat_constraint_error) ||
    gamma_hat_constraint_error > 1e-7
  ) {
    stop(
      "Initial residual patient-effect estimates violate X0^T u_g = 0."
    )
  }

  if (
    !is.finite(kappa_hat_constraint_error) ||
    kappa_hat_constraint_error > 1e-7
  ) {
    stop(
      "Initial patient-by-region estimates violate the sum-to-zero constraint."
    )
  }

  # ---------------------------------------------------------------------------
  # Patient x region groups
  #
  # Used for delta_pr,g^2, the relative patient-region variance multiplier.
  # ---------------------------------------------------------------------------

  patient_region <- interaction(
    patient,
    region,
    drop = TRUE,
    sep = "__"
  )

  group_indices <- split(
    seq_len(n_cells),
    patient_region
  )

  group_names <- names(group_indices)
  n_groups <- length(group_indices)

  # ---------------------------------------------------------------------------
  # Patient and region represented by each group
  # ---------------------------------------------------------------------------

  group_patient <- vapply(
    group_indices,
    function(idx) as.character(patient[idx[1]]),
    character(1)
  )

  group_region <- vapply(
    group_indices,
    function(idx) as.character(region[idx[1]]),
    character(1)
  )

  group_patient_index <- match(
    group_patient,
    levels(patient)
  )

  group_R <- ifelse(
    group_region == edge_level,
    1,
    -1
  )

  n_group <- lengths(group_indices)

  if (any(n_group < 2L)) {
    stop("Every patient-region group must contain at least two cells.")
  }

  # ---------------------------------------------------------------------------
  # Patient-region sufficient statistics
  # ---------------------------------------------------------------------------

  sum_z <- matrix(
    0,
    nrow = n_groups,
    ncol = n_features
  )

  sum_z2 <- matrix(
    0,
    nrow = n_groups,
    ncol = n_features
  )

  delta_hat <- matrix(
    0,
    nrow = n_groups,
    ncol = n_features
  )

  rownames(sum_z) <- group_names
  rownames(sum_z2) <- group_names
  rownames(delta_hat) <- group_names

  for (h in seq_len(n_groups)) {
    idx <- group_indices[[h]]
    z_h <- s_data[, idx, drop = FALSE]

    sum_z[h, ] <- rowSums(z_h)
    sum_z2[h, ] <- rowSums(z_h^2)

    # delta_hat stores the standardized patient-region variance multiplier,
    # i.e. delta_pr,g^2 in the Methods notation.
    delta_hat[h, ] <- matrixStats::rowVars(z_h)
  }

  if (any(!is.finite(delta_hat)) || any(delta_hat <= 0)) {
    stop(
      "At least one feature has zero or non-finite variance ",
      "inside a patient-region group."
    )
  }

  # ---------------------------------------------------------------------------
  # Empirical-Bayes priors for patient effects
  #
  # gamma_pg ~ N(gamma_bar_p, tau_gamma_p^2)
  # ---------------------------------------------------------------------------

  gamma_bar <- rowMeans(gamma_hat)

  tau2_gamma <- pmax(
    matrixStats::rowVars(gamma_hat),
    eps
  )

  # ---------------------------------------------------------------------------
  # Empirical-Bayes priors for patient x region effects
  #
  # kappa_pg ~ N(kappa_bar_p, tau_kappa_p^2)
  # ---------------------------------------------------------------------------

  kappa_bar <- rowMeans(kappa_hat)

  tau2_kappa <- pmax(
    matrixStats::rowVars(kappa_hat),
    eps
  )

  # ---------------------------------------------------------------------------
  # Inverse-gamma hyperparameter estimators
  #
  # Applied to delta_hat = delta_pr,g^2 across features within each
  # patient-region group.
  # ---------------------------------------------------------------------------

  aprior_local <- function(x) {
    m <- mean(x)
    s2 <- max(stats::var(x), eps)

    (2 * s2 + m^2) / s2
  }

  bprior_local <- function(x) {
    m <- mean(x)
    s2 <- max(stats::var(x), eps)

    (m * s2 + m^3) / s2
  }

  a_prior <- apply(
    delta_hat,
    1,
    aprior_local
  )

  b_prior <- apply(
    delta_hat,
    1,
    bprior_local
  )

  # ---------------------------------------------------------------------------
  # Initialize empirical-Bayes posterior estimates
  # ---------------------------------------------------------------------------

  gamma_star <- gamma_hat
  kappa_star <- kappa_hat
  delta_star <- delta_hat

  iterations <- integer(n_features)
  converged <- logical(n_features)

  # ---------------------------------------------------------------------------
  # Identifiability constraints for empirical-Bayes updates
  #
  # Ordinary RPBC:
  #   sum_p gamma_pg = 0
  #   sum_p kappa_pg = 0
  #
  # Covariate-adjusted RPBC:
  #   X0^T gamma_g = 0
  #   sum_p kappa_pg = 0
  # ---------------------------------------------------------------------------

  if (!covariate_adjusted) {

    # Preserve the original ordinary-RPBC constraints exactly.
    constraint <- matrix(
      0,
      nrow = 2,
      ncol = 2 * n_patients
    )

    constraint[
      1,
      seq_len(n_patients)
    ] <- 1

    constraint[
      2,
      n_patients + seq_len(n_patients)
    ] <- 1

  } else {

    n_gamma_constraints <-
      ncol(X0_patient)

    constraint <- matrix(
      0,
      nrow = n_gamma_constraints + 1L,
      ncol = 2 * n_patients
    )

    constraint[
      seq_len(n_gamma_constraints),
      seq_len(n_patients)
    ] <- t(X0_patient)

    constraint[
      n_gamma_constraints + 1L,
      n_patients + seq_len(n_patients)
    ] <- 1
  }


  # ---------------------------------------------------------------------------
  # Constraint-matrix rank check
  # ---------------------------------------------------------------------------

  constraint_rank <- qr(
    constraint
  )$rank

  if (
    constraint_rank != nrow(constraint)
  ) {
    stop(
      "The RPBC identifying constraint matrix is not full row rank."
    )
  }

  # ---------------------------------------------------------------------------
  # Prior mean and variance
  # ---------------------------------------------------------------------------

  prior_mean <- c(
    gamma_bar,
    kappa_bar
  )

  prior_var <- c(
    tau2_gamma,
    tau2_kappa
  )

  message("Finding RPBC parametric empirical-Bayes adjustments...")

  # ---------------------------------------------------------------------------
  # Empirical-Bayes estimation feature by feature
  # ---------------------------------------------------------------------------

  for (g in seq_len(n_features)) {

    theta_old <- c(
      gamma_hat[, g],
      kappa_hat[, g]
    )

    delta_old <- pmax(
      delta_hat[, g],
      eps
    )
    change <- Inf

    for (iter in seq_len(max_iter)) {

      # -----------------------------------------------------------------------
      # Posterior precision
      # -----------------------------------------------------------------------

      precision <- diag(1 / prior_var)
      rhs <- prior_mean / prior_var

      # -----------------------------------------------------------------------
      # Likelihood contribution from each patient x region group
      # -----------------------------------------------------------------------

      for (h in seq_len(n_groups)) {
        p_index <- group_patient_index[h]

        a <- numeric(2 * n_patients)
        a[p_index] <- 1
        a[n_patients + p_index] <- group_R[h]

        precision <- precision +
          (n_group[h] / delta_old[h]) * tcrossprod(a)

        rhs <- rhs +
          (sum_z[h, g] / delta_old[h]) * a
      }

      # -----------------------------------------------------------------------
      # Posterior mean and covariance
      # -----------------------------------------------------------------------

      posterior_var <- solve(precision)
      posterior_mean <- as.vector(
        posterior_var %*% rhs
      )

      # -----------------------------------------------------------------------
      # Enforce identifying constraints
      # -----------------------------------------------------------------------

      VC <- posterior_var %*% t(constraint)

      correction <- VC %*%
        solve(
          constraint %*% VC,
          constraint %*% posterior_mean
        )

      theta_new <- posterior_mean -
        as.vector(correction)

      # -----------------------------------------------------------------------
      # Update patient-region scale parameters
      # -----------------------------------------------------------------------

      delta_new <- numeric(n_groups)

      for (h in seq_len(n_groups)) {
        p_index <- group_patient_index[h]

        mu_h <- theta_new[p_index] +
          theta_new[n_patients + p_index] * group_R[h]

        SSE <- sum_z2[h, g] -
          2 * mu_h * sum_z[h, g] +
          n_group[h] * mu_h^2

        SSE <- max(SSE, 0)

        delta_new[h] <- (
          0.5 * SSE + b_prior[h]
        ) / (
          n_group[h] / 2 + a_prior[h] - 1
        )
      }

      delta_new <- pmax(
        delta_new,
        eps
      )

      # -----------------------------------------------------------------------
      # Convergence
      # -----------------------------------------------------------------------

      change <- max(
        c(
          abs(theta_new - theta_old) /
            (abs(theta_old) + eps),
          abs(delta_new - delta_old) /
            (abs(delta_old) + eps)
        )
      )

      theta_old <- theta_new
      delta_old <- delta_new

      if (change < conv) {
        break
      }
    }

    # -------------------------------------------------------------------------
    # Store posterior estimates
    # -------------------------------------------------------------------------

    gamma_star[, g] <- theta_old[
      seq_len(n_patients)
    ]

    kappa_star[, g] <- theta_old[
      n_patients + seq_len(n_patients)
    ]

    delta_star[, g] <- delta_old

    iterations[g] <- iter
    converged[g] <- is.finite(change) &&
      change < conv
  }

  # ---------------------------------------------------------------------------
  # Verify final identifying constraints
  #
  # Residual patient effect:
  #   X0^T gamma_star = 0
  #
  # Patient x region interaction:
  #   sum_p kappa_star_pg = 0
  # ---------------------------------------------------------------------------

  gamma_star_constraint_error <- max(
    abs(
      t(X0_patient) %*%
        gamma_star
    )
  )

  kappa_star_constraint_error <- max(
    abs(
      colSums(kappa_star)
    )
  )

  if (
    !is.finite(gamma_star_constraint_error) ||
    gamma_star_constraint_error > 1e-7
  ) {
    stop(
      "Posterior residual patient effects violate X0^T u_g = 0."
    )
  }

  if (
    !is.finite(kappa_star_constraint_error) ||
    kappa_star_constraint_error > 1e-7
  ) {
    stop(
      "Posterior patient-by-region effects violate the sum-to-zero constraint."
    )
  }

  # ---------------------------------------------------------------------------
  # Convergence warning
  # ---------------------------------------------------------------------------

  if (any(!converged)) {
    warning(
      sum(!converged),
      " feature(s) did not converge within `max_iter`."
    )
  }

  # ---------------------------------------------------------------------------
  # Harmonization of standardized data
  # ---------------------------------------------------------------------------

  harmonized_standardized <- s_data

  for (h in seq_len(n_groups)) {

    idx <- group_indices[[h]]
    p_index <- group_patient_index[h]

    # -------------------------------------------------------------------------
    # Patient + patient x region nuisance location effect
    # -------------------------------------------------------------------------

    nuisance_mean <- gamma_star[p_index, ] +
      kappa_star[p_index, ] * group_R[h]

    # -------------------------------------------------------------------------
    # Remove nuisance location effect
    # -------------------------------------------------------------------------

    harmonized_group <- sweep(
      s_data[, idx, drop = FALSE],
      MARGIN = 1,
      STATS = nuisance_mean,
      FUN = "-"
    )

    # -------------------------------------------------------------------------
    # Remove patient-region scale effect
    #
    # delta_star contains the variance multiplier delta_pr,g^2,
    # so harmonization divides by sqrt(delta_star).
    # -------------------------------------------------------------------------

    harmonized_group <- sweep(
      harmonized_group,
      MARGIN = 1,
      STATS = sqrt(delta_star[h, ]),
      FUN = "/"
    )

    harmonized_standardized[, idx] <- harmonized_group
  }

  # ---------------------------------------------------------------------------
  # Reconstruction to original feature scale
  #
  # Restore the pooled feature scale and add back the complete
  # protected location component:
  #
  # ordinary RPBC:
  #   alpha_g + beta_Rg R_r
  #
  # covariate-adjusted RPBC:
  #   alpha_g + beta_Rg R_r + x_p^T beta_Xg
  # ---------------------------------------------------------------------------

  harmonized_data <- sweep(
    harmonized_standardized,
    MARGIN = 1,
    STATS = sqrt(var_pooled),
    FUN = "*"
  ) + preserved_mean

  dimnames(harmonized_data) <- dimnames(dat)

  # ---------------------------------------------------------------------------
  # Return RPBC results
  # ---------------------------------------------------------------------------

  list(

    # Harmonized feature matrix
    data = harmonized_data,

    # Protected location model
    alpha_hat = B_hat[1, ],
    beta_region_hat = B_hat[2, ],
    beta_covariate_hat = beta_covariate_hat,
    protected_mean = preserved_mean,

    # Covariate-adjustment information
    covariate_adjusted = covariate_adjusted,
    covariate_patient_data = covariate_info$patient_data,
    covariate_design_patient = X_patient,
    augmented_design_patient = X0_patient,
    covariate_original_names = covariate_info$original_names,
    covariate_design_names = covariate_info$design_names,
    q_covariates = q_covariates,
    residual_patient_df = residual_patient_df,
    residual_patient_basis = residual_patient_basis,

    # Feature-specific baseline residual variance
    var_pooled = var_pooled,

    # Initial standardized nuisance effects
    gamma_hat = gamma_hat,
    kappa_hat = kappa_hat,

    # Empirical-Bayes hyperparameters
    gamma_bar = gamma_bar,
    tau2_gamma = tau2_gamma,
    kappa_bar = kappa_bar,
    tau2_kappa = tau2_kappa,

    # Initial patient-region variance multipliers
    # and inverse-gamma hyperparameters
    delta_hat = delta_hat,
    a_prior = a_prior,
    b_prior = b_prior,

    # Posterior empirical-Bayes estimates
    gamma_star = gamma_star,
    kappa_star = kappa_star,
    delta_star = delta_star,

    # Patient-region group information
    group_names = group_names,
    group_patient = group_patient,
    group_region = group_region,
    n_group = n_group,

    # Model and constraint diagnostics
    design_rank = design_qr$rank,
    design_columns = ncol(design_full),
    constraint_rank = constraint_rank,
    constraint_rows = nrow(constraint),
    nullspace_error = nullspace_error,
    gamma_hat_constraint_error = gamma_hat_constraint_error,
    kappa_hat_constraint_error = kappa_hat_constraint_error,
    gamma_star_constraint_error = gamma_star_constraint_error,
    kappa_star_constraint_error = kappa_star_constraint_error,

    # Convergence diagnostics
    iterations = iterations,
    converged = converged,
    patient_region_table = patient_region_table
  )

}


