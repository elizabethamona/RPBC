# RPBC

## Region-Preserving Batch Correction for Paired Region-Resolved Single-Cell Data

**RPBC** is an R/Bioconductor package implementing **Region-Preserving Batch Correction**, an empirical-Bayes location-scale framework for paired region-resolved single-cell and spatial molecular data.

RPBC is designed for studies in which patient-associated nuisance variation should be reduced while a prespecified **patient-balanced population-average regional contrast** is retained. Rather than treating all between-patient and patient-region variation as interchangeable batch effects, RPBC explicitly separates the regional effect targeted for preservation from patient-wide shifts, patient-specific regional deviations, and patient-region scale heterogeneity.

The package currently supports three related workflows:

- **Ordinary RPBC** — harmonizes patient and patient-by-region nuisance variation while preserving the population-average regional effect.
- **Adaptive RPBC** — adds group-size-dependent stabilization of patient-region scale adjustments, with stronger stabilization for smaller patient-region groups.
- **Covariate-adjusted RPBC** — protects prespecified patient-level covariate-associated variation while harmonizing residual patient-associated variation.

RPBC is designed to interoperate directly with Bioconductor `SingleCellExperiment` and `SpatialExperiment` objects.

---

## Why RPBC?

In paired region-resolved studies, biological structure of interest can coexist with substantial patient-associated heterogeneity. A harmonization procedure that removes patient-related structure without explicitly defining what should be preserved may also attenuate meaningful regional differences.

RPBC addresses this problem by defining the **population-average regional effect as part of the preserved mean structure**. The remaining patient-associated location and scale components are estimated within an empirical-Bayes framework and harmonized around that target.

For the default paired Core/Edge setting, the objective is therefore not simply to make patients more similar. It is to reduce nuisance variation while retaining the population-level Core–Edge contrast specified by the design.

---

## Method at a glance

RPBC uses a crossed patient-by-region representation to distinguish:

1. the overall mean;
2. the population-average regional effect;
3. patient-wide location shifts;
4. patient-specific regional deviations; and
5. patient-region scale heterogeneity.

Empirical-Bayes estimation pools information across features to stabilize nuisance-component estimates. Harmonization and reconstruction then adjust the estimated nuisance components while retaining the prespecified protected mean structure.

The three RPBC workflows share the same regional preservation target but address different analysis settings:

| Workflow | Intended use |
|---|---|
| **Ordinary RPBC** | Paired region-resolved data when the population-average regional effect is the primary protected structure |
| **Adaptive RPBC** | Settings in which patient-region group sizes may differ and additional stabilization of scale estimates is desired |
| **Covariate-adjusted RPBC** | Settings in which known patient-level covariates should remain part of the protected mean structure |

Adaptive RPBC and covariate-adjusted RPBC are currently implemented as separate extensions and should not be requested simultaneously.

---

## Installation

RPBC is currently under review for inclusion in **Bioconductor**.

Until a Bioconductor release is available, the development version can be installed from GitHub:

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}

remotes::install_github("elizabethamona/RPBC")
```

RPBC requires R >= 4.6.0 and uses Bioconductor infrastructure including `SingleCellExperiment`, `SpatialExperiment`, and `SummarizedExperiment`.

After installation:

```r
library(RPBC)
```

---

## Basic usage

The primary user-facing function is `rpbc()`.

The function accepts a `SingleCellExperiment` or `SpatialExperiment` object and stores the harmonized expression matrix as a **new assay**, leaving the selected input assay available in the object.

The package includes `haley_gbm_example`, a representative paired Core/Edge imaging mass cytometry dataset that can be used to explore the interface.

```r
library(RPBC)

data("haley_gbm_example", package = "RPBC")
```

### Ordinary RPBC

```r
fit <- rpbc(
    object = haley_gbm_example,
    patient = "patient_id",
    region = "region_type",
    assay_name = "counts",
    transform = "asinh",
    cofactor = 5,
    adaptive = FALSE,
    output_assay = "rpbc"
)
```

The harmonized values are stored in the new `"rpbc"` assay:

```r
SummarizedExperiment::assayNames(fit)

rpbc_matrix <- SummarizedExperiment::assay(
    fit,
    "rpbc"
)
```

---

## Adaptive RPBC

Adaptive RPBC uses the same regional preservation target while introducing group-size-dependent stabilization of patient-region scale adjustment.

```r
fit_adaptive <- rpbc(
    object = haley_gbm_example,
    patient = "patient_id",
    region = "region_type",
    assay_name = "counts",
    transform = "asinh",
    cofactor = 5,
    adaptive = TRUE,
    output_assay = "rpbc_adaptive"
)
```

The default settings are:

```r
delta_adaptive_k = 1
delta_bounds = c(0.5, 2)
```

---

## Covariate-adjusted RPBC

Covariate-adjusted RPBC can be used when prespecified patient-level variables should remain part of the protected mean structure.

For example:

```r
fit_covariate <- rpbc(
    object = haley_gbm_example,
    patient = "patient_id",
    region = "region_type",
    assay_name = "counts",
    transform = "asinh",
    cofactor = 5,
    covariates = c("age", "sex"),
    adaptive = FALSE,
    output_assay = "rpbc_covariate"
)
```

Covariates supplied through `covariates` must be constant within patient.

---

## Input and transformation

RPBC separates **assay selection** from **transformation**.

Available transformations are:

- `"none"` — use the selected assay as supplied;
- `"asinh"` — apply `asinh(x / cofactor)`;
- `"log1p"` — apply `log(x + 1)` to non-negative values.

For imaging mass cytometry data, an `asinh` transformation may be appropriate. For sequencing data that have already been appropriately normalized and transformed, such as a `"logcounts"` assay, users can provide that assay with:

```r
transform = "none"
```

The `"log1p"` option performs only the mathematical transformation `log(x + 1)` and is **not a replacement for sequencing-data normalization**.

---

## Design requirements

RPBC is intended for **paired region-resolved designs**. Each patient must contribute observations from both specified region levels.

By default:

```r
core_level = "Core"
edge_level = "Edge"
```

Alternative labels can be supplied through `core_level` and `edge_level`.

For covariate-adjusted RPBC, patient-level covariates supplied at the cell level must remain constant within patient.

---

## Output and downstream analysis

By default, `rpbc()` returns the input Bioconductor object with the harmonized matrix stored in the assay specified by `output_assay`.

The harmonized matrix retains the same feature-by-cell dimensions and dimnames as the selected input assay.

RPBC does not prescribe a particular downstream analysis workflow. Harmonized values can be used for exploratory and descriptive analyses such as:

- dimensionality reduction;
- clustering;
- visualization; and
- cell-state or phenotype characterization.

For formal inferential analyses, users should employ statistical methods appropriate to the study design and scientific question.

Users who require fitted model quantities and convergence information can set:

```r
return_fit = TRUE
```

which returns both the modified object and the internal RPBC fit.

---

## Documentation

A complete worked example is provided in the package vignette:

```r
vignette("RPBC")
```

Function-level documentation is available through:

```r
?rpbc
```

---

## Example data

The bundled `haley_gbm_example` object is a representative subset of published imaging mass cytometry data from IDH-wildtype glioblastoma and retains the paired Core/Edge structure used to demonstrate RPBC.

The source data are described in:

Haley MJ, et al. (2024). *Hypoxia coordinates the spatial landscape of myeloid cells within glioblastoma to affect survival*. **Science Advances**, 10(20), eadj3301.  
https://doi.org/10.1126/sciadv.adj3301

Haley M, Couper K. (2024). *Imaging mass cytometry data from IDH wildtype glioblastomas*. Dryad.  
https://doi.org/10.5061/dryad.2z34tmprw

---

## Associated manuscript

The methodological manuscript describing RPBC, its extensions, theoretical foundations, simulation evaluation, and real-data application is currently in preparation.

The complete manuscript citation and manuscript author list will be added here when finalized.

---

## Software author and maintainer

**Elizabeth B. Amona, PhD**  
Biostatistics Core, Brown Cancer Center  
University of Louisville  
ORCID: https://orcid.org/0000-0002-7575-788X

---

## License

RPBC is released under the **MIT License**.

Copyright © 2026 Elizabeth B. Amona.
