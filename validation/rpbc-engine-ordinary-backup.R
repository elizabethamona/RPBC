# Internal ordinary RPBC statistical engine
#
# This function operates on an expression matrix and paired patient/region
# vectors. It contains the statistical RPBC algorithm and is intentionally
# independent of SingleCellExperiment and SpatialExperiment objects.

.rpbc_fit_matrix <- function(
    dat,
    patient,
    region,
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

  patient_contrast <- contr.sum(n_patients)
  rownames(patient_contrast) <- levels(patient)

  P <- patient_contrast[
    as.character(patient),
    ,
    drop = FALSE
  ]

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
  # alpha_g
  # + population-average region effect
  # + patient effect
  # + patient x region interaction
  # ---------------------------------------------------------------------------

  design_full <- cbind(
    intercept = 1,
    region = R,
    P,
    PR
  )

  if (qr(design_full)$rank != ncol(design_full)) {
    stop("RPBC design matrix is not full rank.")
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
  # Biological component to preserve
  #
  # Preserve the feature-specific intercept and population-average
  # regional effect.
  # ---------------------------------------------------------------------------

  preserved_mean <- design_full[
    ,
    1:2,
    drop = FALSE
  ] %*% B_hat[
    1:2,
    ,
    drop = FALSE
  ]

  preserved_mean <- t(preserved_mean)

  # ---------------------------------------------------------------------------
  # Standardize data
  #
  # Remove the protected biological mean and scale each feature by its
  # pooled residual standard deviation.
  # ---------------------------------------------------------------------------

  s_data <- sweep(
    dat - preserved_mean,
    MARGIN = 1,
    STATS = sqrt(var_pooled),
    FUN = "/"
  )

  # ---------------------------------------------------------------------------
  # Initial patient and patient x region coefficient positions
  # ---------------------------------------------------------------------------

  q <- n_patients - 1L

  patient_rows <- 2L + seq_len(q)

  interaction_rows <- 2L + q + seq_len(q)

  # ---------------------------------------------------------------------------
  # Initial patient effects
  #
  # Convert the raw-scale patient coefficients to standardized gamma_pg.
  # ---------------------------------------------------------------------------

  gamma_contrast <- B_hat[
    patient_rows,
    ,
    drop = FALSE
  ]

  gamma_contrast <- sweep(
    gamma_contrast,
    MARGIN = 2,
    STATS = sqrt(var_pooled),
    FUN = "/"
  )

  # ---------------------------------------------------------------------------
  # Initial patient x region interaction effects
  #
  # Convert the raw-scale interaction coefficients to standardized kappa_pg.
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

  # ---------------------------------------------------------------------------
  # Reconstruct effects for all patients
  #
  # Enforce:
  # sum_p gamma_pg = 0
  # sum_p kappa_pg = 0
  # ---------------------------------------------------------------------------

  gamma_hat <- rbind(
    gamma_contrast,
    -colSums(gamma_contrast)
  )

  kappa_hat <- rbind(
    kappa_contrast,
    -colSums(kappa_contrast)
  )

  rownames(gamma_hat) <- levels(patient)
  rownames(kappa_hat) <- levels(patient)

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
    s2 <- max(var(x), eps)

    (2 * s2 + m^2) / s2
  }

  bprior_local <- function(x) {
    m <- mean(x)
    s2 <- max(var(x), eps)

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
  # Identifiability constraints
  #
  # Enforce:
  # sum_p gamma_pg = 0
  # sum_p kappa_pg = 0
  # ---------------------------------------------------------------------------

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
      # Enforce sum-to-zero constraints
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
  # Restore the pooled feature scale and add back the preserved
  # population-average regional location component.
  # ---------------------------------------------------------------------------

  harmonized_data <- sweep(
    harmonized_standardized,
    MARGIN = 1,
    STATS = sqrt(var_pooled),
    FUN = "*"
  ) + preserved_mean

  dimnames(harmonized_data) <- dimnames(dat)

  # ---------------------------------------------------------------------------
  # Return ordinary RPBC results
  # ---------------------------------------------------------------------------

  list(
    # Harmonized feature matrix
    data = harmonized_data,

    # Initial raw-scale location-model estimates
    alpha_hat = B_hat[1, ],
    beta_region_hat = B_hat[2, ],

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

    # Model and convergence diagnostics
    design_rank = qr(design_full)$rank,
    design_columns = ncol(design_full),
    iterations = iterations,
    converged = converged,
    patient_region_table = patient_region_table
  )

}


