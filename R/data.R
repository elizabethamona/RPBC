###############################################################################
# Package datasets
###############################################################################

#' Haley GBM IMC Example Dataset
#'
#' A representative example dataset for demonstrating Region-Preserving Batch
#' Correction (RPBC) using glioblastoma imaging mass cytometry (IMC) data.
#' The dataset retains all seven patients and the paired Core/Edge design, with
#' 30 cells sampled from each patient-region combination.
#'
#' The dataset contains 10 representative biological markers measured across
#' 420 cells. Raw IMC intensities are stored in the `"counts"` assay, allowing
#' preprocessing transformations to be applied directly through [rpbc()].
#'
#' @format A `SpatialExperiment` object with 10 markers (rows) and
#'   420 cells (columns). The `"counts"` assay contains raw IMC intensities.
#'
#' Cell-level metadata include:
#' \describe{
#'   \item{cell_id}{Unique cell identifier.}
#'   \item{ROI}{Region-of-interest identifier.}
#'   \item{patient_id}{Patient identifier.}
#'   \item{region_type}{Paired tumor region, coded as `"Core"` or `"Edge"`.}
#'   \item{age}{Patient age.}
#'   \item{sex}{Patient sex.}
#'   \item{mgmt_status}{MGMT status.}
#' }
#'
#' @details
#' This example dataset was constructed by sampling 30 cells from each
#' patient-region combination while retaining all seven patients and both
#' Core and Edge regions. It is intended for package examples, testing, and
#' vignettes and does not represent the complete Haley GBM IMC dataset.
#'
#' The 10 included markers are CD45, CD3, Cd8a, CD68, CD163, GFAP, Sox2,
#' Olig2, Ki67, and Vimentin.
#'
#' @source
#' Haley MJ, et al. (2024). Hypoxia coordinates the spatial landscape of
#' myeloid cells within glioblastoma to affect survival.
#' \emph{Science Advances}, 10(20), eadj3301.
#' \url{https://doi.org/10.1126/sciadv.adj3301}
#'
#' Haley M, Couper K. (2024). Imaging mass cytometry data from IDH wildtype
#' glioblastomas. Dryad.
#' \url{https://doi.org/10.5061/dryad.2z34tmprw}
#'
#' @usage data(haley_gbm_example)
#'
#' @keywords datasets
#'
"haley_gbm_example"
