################################################################################
# COVARIATE PREPARATION FOR RPBC
#
# Internal package helper
#
# Prepares the patient-level covariate design used by RPBC.
# With covariates = NULL, the helper returns the ordinary-RPBC configuration.
# Covariates may be supplied once per patient or once per cell,
# provided cell-level values are constant within patient.
################################################################################


###############################################################################
# PART 1: PATIENT-LEVEL COVARIATE PREPARATION
#
# Covariates may be supplied:
#   1. once per patient; or
#   2. once per cell, provided they are constant within patient.
#
# Categorical covariates are expanded with model.matrix().
# Design columns are centered across patients.
#
# X0 = [1, X]
#
# The residual patient effect is constrained by
#
#   X0^T u_g = 0.
###############################################################################

.rpbc_prepare_covariates <- function(
    covariates,
    patient
) {

  patient <- droplevels(
    factor(patient)
  )

  patient_levels <- levels(patient)
  n_patients <- length(patient_levels)
  n_cells <- length(patient)

  if (is.null(covariates)) {

    return(
      list(
        adjusted = FALSE,
        patient_data = NULL,
        original_names = character(0),
        X = matrix(
          numeric(0),
          nrow = n_patients,
          ncol = 0,
          dimnames = list(
            patient_levels,
            NULL
          )
        ),
        X0 = matrix(
          1,
          nrow = n_patients,
          ncol = 1,
          dimnames = list(
            patient_levels,
            "intercept"
          )
        ),
        design_names = character(0),
        q = 0L,
        residual_df = n_patients - 1L
      )
    )
  }


  ###########################################################################
  # Convert a single vector/factor to a one-column data.frame
  ###########################################################################

  if (
    is.atomic(covariates) &&
    is.null(dim(covariates))
  ) {

    covariates <- data.frame(
      covariate = covariates,
      check.names = FALSE
    )
  }


  ###########################################################################
  # Basic input checks
  ###########################################################################

  if (
    !is.data.frame(covariates) &&
    !is.matrix(covariates)
  ) {
    stop(
      "covariates must be a vector/factor, matrix, or data.frame."
    )
  }

  covariates <- as.data.frame(
    covariates,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (ncol(covariates) == 0L) {
    return(
      .rpbc_prepare_covariates(
        covariates = NULL,
        patient = patient
      )
    )
  }

  if (any(!nzchar(names(covariates)))) {
    names(covariates)[!nzchar(names(covariates))] <-
      paste0(
        "covariate",
        which(!nzchar(names(covariates)))
      )
  }

  names(covariates) <-
    make.unique(
      names(covariates)
    )

  unsupported <- vapply(
    covariates,
    function(x) {
      is.list(x) &&
        !is.factor(x)
    },
    logical(1)
  )

  if (any(unsupported)) {
    stop(
      "covariates contains unsupported list-column(s): ",
      paste(
        names(covariates)[unsupported],
        collapse = ", "
      ),
      "."
    )
  }

  if (anyNA(covariates)) {
    stop(
      "covariates contains missing values."
    )
  }

  numeric_bad <- vapply(
    covariates,
    function(x) {
      is.numeric(x) &&
        any(!is.finite(x))
    },
    logical(1)
  )

  if (any(numeric_bad)) {
    stop(
      "covariates contains non-finite numeric values."
    )
  }


  ###########################################################################
  # Convert cell-level input to one row per patient
  ###########################################################################

  if (nrow(covariates) == n_cells) {

    first_index <- match(
      patient_levels,
      as.character(patient)
    )

    for (j in seq_len(ncol(covariates))) {

      for (p in patient_levels) {

        idx <- which(
          as.character(patient) == p
        )

        values_p <- covariates[[j]][idx]

        if (
          length(unique(values_p)) != 1L
        ) {
          stop(
            "Covariate '",
            names(covariates)[j],
            "' is not constant within patient '",
            p,
            "'. RPBC covariates must be patient-level."
          )
        }
      }
    }

    patient_covariates <- covariates[
      first_index,
      ,
      drop = FALSE
    ]

    rownames(patient_covariates) <-
      patient_levels


    ###########################################################################
    # Patient-level input
    ###########################################################################

  } else if (
    nrow(covariates) == n_patients
  ) {

    supplied_rownames <-
      rownames(covariates)

    default_rownames <- identical(
      supplied_rownames,
      as.character(
        seq_len(n_patients)
      )
    )

    if (
      !is.null(supplied_rownames) &&
      length(unique(supplied_rownames)) == n_patients &&
      setequal(
        supplied_rownames,
        patient_levels
      )
    ) {

      patient_covariates <- covariates[
        match(
          patient_levels,
          supplied_rownames
        ),
        ,
        drop = FALSE
      ]

    } else if (
      is.null(supplied_rownames) ||
      default_rownames
    ) {

      patient_covariates <- covariates

    } else {

      stop(
        "Patient-level covariate row names must match patient labels."
      )
    }

    rownames(patient_covariates) <-
      patient_levels

  } else {

    stop(
      "covariates must have either one row per cell (",
      n_cells,
      ") or one row per patient (",
      n_patients,
      ")."
    )
  }


  ###########################################################################
  # Prepare categorical variables
  ###########################################################################

  for (j in seq_len(ncol(patient_covariates))) {

    if (
      is.character(
        patient_covariates[[j]]
      )
    ) {

      patient_covariates[[j]] <-
        factor(
          patient_covariates[[j]]
        )

    } else if (
      is.factor(
        patient_covariates[[j]]
      )
    ) {

      patient_covariates[[j]] <-
        droplevels(
          patient_covariates[[j]]
        )
    }
  }

  original_names <-
    names(patient_covariates)

  names(patient_covariates) <-
    make.names(
      original_names,
      unique = TRUE
    )


  ###########################################################################
  # Construct covariate design matrix
  ###########################################################################

  mm <- tryCatch(
    stats::model.matrix(
      ~ .,
      data = patient_covariates
    ),
    error = function(e) {
      stop(
        "Unable to construct the patient-level covariate design: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  X <- mm[
    ,
    colnames(mm) != "(Intercept)",
    drop = FALSE
  ]

  if (ncol(X) == 0L) {
    stop(
      "covariates did not produce any estimable design columns."
    )
  }


  ###########################################################################
  # Center design columns equally across patients
  ###########################################################################

  X <- sweep(
    X,
    MARGIN = 2,
    STATS = colMeans(X),
    FUN = "-"
  )

  rownames(X) <-
    patient_levels

  colnames(X) <-
    make.unique(
      colnames(X)
    )

  if (any(!is.finite(X))) {
    stop(
      "The coded/centered covariate design contains non-finite values."
    )
  }


  ###########################################################################
  # Augmented design
  #
  # X0 = [1, X]
  ###########################################################################

  X0 <- cbind(
    intercept = 1,
    X
  )

  if (
    qr(X0)$rank !=
    ncol(X0)
  ) {
    stop(
      "The augmented patient-level design X0 = [1, X] is not full rank."
    )
  }


  ###########################################################################
  # Remaining residual-patient degrees of freedom
  ###########################################################################

  q_covariates <- ncol(X)

  residual_patient_df <-
    n_patients -
    ncol(X0)

  if (residual_patient_df < 0L) {
    stop(
      "The covariate design contains too many independent columns ",
      "for the number of patients."
    )
  }

  if (residual_patient_df == 0L) {
    warning(
      "The covariate design exhausts the patient-effect degrees of freedom. ",
      "The residual patient effect will therefore be zero."
    )
  }


  ###########################################################################
  # Return
  ###########################################################################

  list(
    adjusted = TRUE,
    patient_data = patient_covariates,
    original_names = original_names,
    X = X,
    X0 = X0,
    design_names = colnames(X),
    q = q_covariates,
    residual_df = residual_patient_df
  )
}


