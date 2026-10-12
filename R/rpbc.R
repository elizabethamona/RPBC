#' Region-Preserving Batch Correction
#'
#' Apply Region-Preserving Batch Correction (RPBC) to paired
#' region-resolved single-cell or spatial data.
#'
#' RPBC harmonizes patient and patient-by-region nuisance variation while
#' preserving the prespecified population-average regional effect.
#'
#' The function supports:
#'
#' - ordinary RPBC;
#' - covariate-adjusted RPBC; and
#' - Adaptive RPBC.
#'
#' Ordinary and covariate-adjusted RPBC use the same internal statistical
#' engine. Adaptive RPBC uses a separate internal engine because its
#' patient-region scale update includes group-size-dependent stabilization.
#'
#' @param object A `SingleCellExperiment` or `SpatialExperiment` object.
#'
#' @param patient Character string naming the column in `colData(object)`
#'   containing patient identifiers.
#'
#' @param region Character string naming the column in `colData(object)`
#'   containing the paired region labels.
#'
#' @param assay_name Character string naming the assay to use as input.
#'   Default is `"counts"`. For count-based sequencing data, an appropriately
#'   normalized and log-transformed assay such as `"logcounts"` (e.g.,
#'   log2(CPM + 1)) may be supplied with `transform = "none"`.
#'
#' @param transform Character string specifying the transformation applied
#'   to the selected assay before RPBC harmonization. One of `"none"`,
#'   `"asinh"` or `"log1p"`. Default is `"none"`.
#'
#' @param cofactor Positive numeric value used when
#'   `transform = "asinh"`. The transformation is
#'   `asinh(x / cofactor)`. Default is `5`.
#'
#' @param covariates Optional character vector naming patient-level
#'   covariates in `colData(object)`. Covariates supplied at the cell level
#'   must be constant within patient. If `NULL`, ordinary RPBC is used.
#'
#' @param adaptive Logical. If `TRUE`, use Adaptive RPBC. Adaptive RPBC
#'   cannot currently be combined with covariate adjustment.
#'
#' @param output_assay Character string giving the name of the assay in which
#'   the harmonized expression matrix will be stored. Default is `"rpbc"`.
#'
#' @param core_level Character string identifying the Core region.
#'
#' @param edge_level Character string identifying the Edge region.
#'
#' @param conv Positive convergence tolerance for empirical-Bayes estimation.
#'
#' @param max_iter Positive integer giving the maximum number of
#'   empirical-Bayes iterations.
#'
#' @param eps Small positive numerical constant used for numerical stability.
#'
#' @param delta_adaptive_k Positive Adaptive RPBC stabilization constant.
#'   The validated default is `1`.
#'
#' @param delta_bounds Numeric vector of length two giving the lower and
#'   upper bounds for Adaptive RPBC patient-region variance multipliers.
#'   The validated default is `c(0.5, 2.0)`.
#'
#' @param return_fit Logical. If `FALSE`, return the modified Bioconductor
#'   object. If `TRUE`, return a list containing both the modified object
#'   and the full internal RPBC fit.
#'
#' @return By default, the input object with the harmonized expression matrix
#'   stored in `output_assay`. If `return_fit = TRUE`, a list containing
#'   `object` and `fit`.
#'
#' @examples
#' data("haley_gbm_example", package = "RPBC")
#'
#' # Ordinary RPBC
#' example_ordinary <- rpbc(
#'   object = haley_gbm_example,
#'   patient = "patient_id",
#'   region = "region_type",
#'   assay_name = "counts",
#'   transform = "asinh",
#'   cofactor = 5,
#'   adaptive = FALSE,
#'   output_assay = "rpbc"
#' )
#'
#' # Adaptive RPBC
#' example_adaptive <- rpbc(
#'   object = haley_gbm_example,
#'   patient = "patient_id",
#'   region = "region_type",
#'   assay_name = "counts",
#'   transform = "asinh",
#'   cofactor = 5,
#'   adaptive = TRUE,
#'   output_assay = "rpbc_adaptive"
#' )
#'
#' # Covariate-adjusted RPBC
#' example_covariate <- rpbc(
#'   object = haley_gbm_example,
#'   patient = "patient_id",
#'   region = "region_type",
#'   assay_name = "counts",
#'   transform = "asinh",
#'   cofactor = 5,
#'   covariates = c("age", "sex"),
#'   adaptive = FALSE,
#'   output_assay = "rpbc_covariate"
#' )
#'
#' @export
rpbc <- function(
    object,
    patient,
    region,
    assay_name = "counts",
    transform = c(
      "none",
      "asinh",
      "log1p"
    ),
    cofactor = 5,
    covariates = NULL,
    adaptive = FALSE,
    output_assay = "rpbc",
    core_level = "Core",
    edge_level = "Edge",
    conv = 1e-4,
    max_iter = 100L,
    eps = 1e-8,
    delta_adaptive_k = 1,
    delta_bounds = c(0.5, 2.0),
    return_fit = FALSE
) {

  # -------------------------------------------------------------------------
  # Validate object class
  # -------------------------------------------------------------------------

  if (
    !inherits(object, "SingleCellExperiment") &&
    !inherits(object, "SpatialExperiment")
  ) {
    stop(
      "`object` must be a SingleCellExperiment or SpatialExperiment object."
    )
  }


  # -------------------------------------------------------------------------
  # Validate column-name arguments
  # -------------------------------------------------------------------------

  if (
    !is.character(patient) ||
    length(patient) != 1L ||
    is.na(patient) ||
    !nzchar(patient)
  ) {
    stop(
      "`patient` must be a single non-empty column name."
    )
  }

  if (
    !is.character(region) ||
    length(region) != 1L ||
    is.na(region) ||
    !nzchar(region)
  ) {
    stop(
      "`region` must be a single non-empty column name."
    )
  }

  if (
    !is.character(assay_name) ||
    length(assay_name) != 1L ||
    is.na(assay_name) ||
    !nzchar(assay_name)
  ) {
    stop(
      "`assay_name` must be a single non-empty assay name."
    )
  }

  if (
    !is.character(output_assay) ||
    length(output_assay) != 1L ||
    is.na(output_assay) ||
    !nzchar(output_assay)
  ) {
    stop(
      "`output_assay` must be a single non-empty assay name."
    )
  }

  # -------------------------------------------------------------------------
  # Validate transformation
  # -------------------------------------------------------------------------

  transform <- match.arg(
    transform
  )

  if (
    !is.numeric(cofactor) ||
    length(cofactor) != 1L ||
    !is.finite(cofactor) ||
    cofactor <= 0
  ) {
    stop(
      "`cofactor` must be a single positive finite numeric value."
    )
  }

  # -------------------------------------------------------------------------
  # Validate logical arguments
  # -------------------------------------------------------------------------

  if (
    !is.logical(adaptive) ||
    length(adaptive) != 1L ||
    is.na(adaptive)
  ) {
    stop(
      "`adaptive` must be TRUE or FALSE."
    )
  }

  if (
    !is.logical(return_fit) ||
    length(return_fit) != 1L ||
    is.na(return_fit)
  ) {
    stop(
      "`return_fit` must be TRUE or FALSE."
    )
  }


  # -------------------------------------------------------------------------
  # Check assay
  # -------------------------------------------------------------------------

  available_assays <-
    SummarizedExperiment::assayNames(
      object
    )

  if (
    !assay_name %in%
    available_assays
  ) {
    stop(
      "Assay '",
      assay_name,
      "' was not found in `object`. Available assays: ",
      paste(
        available_assays,
        collapse = ", "
      ),
      "."
    )
  }

  if (
    output_assay %in%
    available_assays
  ) {
    stop(
      "An assay named '",
      output_assay,
      "' already exists. Choose a different `output_assay` name."
    )
  }


  # -------------------------------------------------------------------------
  # Extract cell metadata
  # -------------------------------------------------------------------------

  cell_data <- as.data.frame(
    SummarizedExperiment::colData(
      object
    )
  )

  if (
    !patient %in%
    names(cell_data)
  ) {
    stop(
      "Patient column '",
      patient,
      "' was not found in `colData(object)`."
    )
  }

  if (
    !region %in%
    names(cell_data)
  ) {
    stop(
      "Region column '",
      region,
      "' was not found in `colData(object)`."
    )
  }

  patient_vector <-
    cell_data[[patient]]

  region_vector <-
    cell_data[[region]]


  # -------------------------------------------------------------------------
  # Prepare optional covariates
  # -------------------------------------------------------------------------

  covariate_data <- NULL

  if (!is.null(covariates)) {

    if (
      !is.character(covariates) ||
      length(covariates) < 1L ||
      anyNA(covariates) ||
      any(!nzchar(covariates))
    ) {
      stop(
        "`covariates` must be NULL or a character vector of column names."
      )
    }

    covariates <- unique(
      covariates
    )

    missing_covariates <- setdiff(
      covariates,
      names(cell_data)
    )

    if (
      length(missing_covariates) > 0L
    ) {
      stop(
        "The following covariate column(s) were not found in ",
        "`colData(object)`: ",
        paste(
          missing_covariates,
          collapse = ", "
        ),
        "."
      )
    }

    if (
      any(
        covariates %in%
        c(
          patient,
          region
        )
      )
    ) {
      stop(
        "`covariates` cannot include the patient or region column."
      )
    }

    covariate_data <- cell_data[
      ,
      covariates,
      drop = FALSE
    ]
  }


  # -------------------------------------------------------------------------
  # Adaptive and covariate-adjusted RPBC are separate extensions
  # -------------------------------------------------------------------------

  if (
    adaptive &&
    !is.null(covariate_data)
  ) {
    stop(
      "Adaptive RPBC and covariate-adjusted RPBC are currently separate ",
      "extensions. Use either `adaptive = TRUE` or supply `covariates`, ",
      "but not both."
    )
  }


  # -------------------------------------------------------------------------
  # Extract feature-by-cell assay matrix
  # -------------------------------------------------------------------------

  dat <- as.matrix(
    SummarizedExperiment::assay(
      object,
      assay_name
    )
  )

  # -------------------------------------------------------------------------
  # Validate selected assay
  # -------------------------------------------------------------------------

  if (!is.numeric(dat)) {
    stop(
      "The selected assay must contain numeric values."
    )
  }

  if (any(!is.finite(dat))) {
    stop(
      "The selected assay contains NA, NaN, or Inf values."
    )
  }


  # -------------------------------------------------------------------------
  # Apply optional transformation before RPBC harmonization
  #
  # none  : use the selected assay exactly as supplied
  # asinh : asinh(x / cofactor)
  # log1p : log(x + 1)
  # -------------------------------------------------------------------------

  if (transform == "none") {

    # No transformation is applied.

  } else if (transform == "asinh") {

    dat <- asinh(
      dat / cofactor
    )

  } else if (transform == "log1p") {

    if (any(dat < 0)) {
      stop(
        "`transform = \"log1p\"` requires non-negative input values."
      )
    }

    dat <- log1p(
      dat
    )
  }
  # -------------------------------------------------------------------------
  # Run appropriate RPBC engine
  # -------------------------------------------------------------------------

  if (adaptive) {

    fit <- .rpbc_fit_matrix_adaptive(
      dat = dat,
      patient = patient_vector,
      region = region_vector,
      core_level = core_level,
      edge_level = edge_level,
      conv = conv,
      max_iter = max_iter,
      eps = eps,
      delta_adaptive_k = delta_adaptive_k,
      delta_bounds = delta_bounds
    )

  } else {

    fit <- .rpbc_fit_matrix(
      dat = dat,
      patient = patient_vector,
      region = region_vector,
      covariates = covariate_data,
      core_level = core_level,
      edge_level = edge_level,
      conv = conv,
      max_iter = max_iter,
      eps = eps
    )
  }


  # -------------------------------------------------------------------------
  # Store harmonized matrix as a new assay
  # -------------------------------------------------------------------------

  SummarizedExperiment::assay(
    object,
    output_assay
  ) <- fit$data


  # -------------------------------------------------------------------------
  # Return
  # -------------------------------------------------------------------------

  if (return_fit) {

    return(
      list(
        object = object,
        fit = fit
      )
    )
  }

  object
}
