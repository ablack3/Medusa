# Copyright 2026 Observational Health Data Sciences and Informatics
#
# This file is part of Medusa
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# --- MR Theme for ggplot2 ---

#' Medusa ggplot2 Theme
#'
#' @title Consistent plotting theme for Medusa package
#'
#' @description
#' A clean ggplot2 theme applied to all plots produced by the Medusa package.
#' Based on \code{theme_minimal} with consistent font sizes, axis formatting,
#' and color scheme.
#'
#' @param baseSize Base font size in points. Default is 12.
#'
#' @return A ggplot2 theme object.
#'
#' @examples
#' library(ggplot2)
#' ggplot(data.frame(x = 1:10, y = rnorm(10)), aes(x, y)) +
#'   geom_point() +
#'   mrTheme()
#'
#' @export
mrTheme <- function(baseSize = 12) {
  ggplot2::theme_minimal(base_size = baseSize) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = baseSize + 2),
      plot.subtitle = ggplot2::element_text(color = "grey40", size = baseSize),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "grey30"),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        color = "grey80",
        fill = NA,
        linewidth = 0.5
      ),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

# --- Color Palette ---

#' @keywords internal
MR_COLORS <- list(
  primary = "#1B4F72",
  secondary = "#85C1E9",
  site = "#AED6F1",
  combined = "#1B4F72",
  ci = "#D4E6F1",
  warning = "#E74C3C",
  success = "#27AE60",
  neutral = "#95A5A6",
  methods = c(
    IVW = "#1B4F72",
    MREgger = "#E74C3C",
    WeightedMedian = "#27AE60",
    Steiger = "#F39C12",
    LeaveOneOut = "#8E44AD"
  )
)


# --- Input Validation Helpers ---

#' Validate Instrument Table Structure
#'
#' @description Checks that an instrument table has the required columns and
#'   valid data types.
#'
#' @param instrumentTable A data frame to validate.
#'
#' @return Invisible TRUE if valid; throws an error otherwise.
#'
#' @keywords internal
validateInstrumentTable <- function(instrumentTable) {
  checkmate::assertDataFrame(instrumentTable, min.rows = 1,
                             .var.name = "instrumentTable")
  requiredCols <- c("snp_id", "effect_allele", "other_allele",
                     "beta_ZX", "se_ZX", "pval_ZX", "eaf")
  missingCols <- setdiff(requiredCols, names(instrumentTable))
  if (length(missingCols) > 0) {
    stop(sprintf("instrumentTable is missing required columns: %s",
                 paste(missingCols, collapse = ", ")))
  }
  checkmate::assertNumeric(instrumentTable$beta_ZX, any.missing = FALSE,
                           .var.name = "instrumentTable$beta_ZX")
  checkmate::assertNumeric(instrumentTable$se_ZX, lower = 0, any.missing = FALSE,
                           .var.name = "instrumentTable$se_ZX")
  checkmate::assertNumeric(instrumentTable$eaf, lower = 0, upper = 1,
                           .var.name = "instrumentTable$eaf")
  invisible(TRUE)
}


#' Validate Beta Grid
#'
#' @description Checks that a beta grid is a sorted numeric vector with
#'   sufficient resolution.
#'
#' @param betaGrid Numeric vector of grid points.
#'
#' @return Invisible TRUE if valid; throws an error otherwise.
#'
#' @keywords internal
validateBetaGrid <- function(betaGrid) {
  checkmate::assertNumeric(betaGrid, min.len = 10, any.missing = FALSE,
                           sorted = TRUE, .var.name = "betaGrid")
  invisible(TRUE)
}


#' Validate Site Profile Object
#'
#' @description Checks that a site profile object has the required structure.
#'
#' @param siteProfile A list to validate.
#'
#' @return Invisible TRUE if valid; throws an error otherwise.
#'
#' @keywords internal
validateSiteProfile <- function(siteProfile) {
  checkmate::assertList(siteProfile, .var.name = "siteProfile")
  requiredElements <- c("siteId", "betaGrid", "logLikProfile",
                         "nCases", "nControls", "snpIds")
  missingElements <- setdiff(requiredElements, names(siteProfile))
  if (length(missingElements) > 0) {
    stop(sprintf("Site profile is missing required elements: %s",
                 paste(missingElements, collapse = ", ")))
  }
  checkmate::assertNumeric(siteProfile$logLikProfile, any.missing = FALSE,
                           .var.name = "logLikProfile")
  if (length(siteProfile$betaGrid) != length(siteProfile$logLikProfile)) {
    stop("betaGrid and logLikProfile must have the same length.")
  }
  if (!is.null(siteProfile$scoreDefinition)) {
    checkmate::assertList(siteProfile$scoreDefinition,
                          .var.name = "siteProfile$scoreDefinition")
    checkmate::assertSubset(c("snpIds", "scoreWeights", "betaZX", "seZX"),
                            names(siteProfile$scoreDefinition))
    checkmate::assertCharacter(siteProfile$scoreDefinition$snpIds,
                               any.missing = FALSE)
    checkmate::assertNumeric(siteProfile$scoreDefinition$scoreWeights,
                             any.missing = FALSE)
    if (length(siteProfile$scoreDefinition$snpIds) !=
        length(siteProfile$scoreDefinition$scoreWeights)) {
      stop("scoreDefinition$snpIds and scoreDefinition$scoreWeights must have the same length.")
    }
  }
  invisible(TRUE)
}


#' Compute Approximate F-statistic from GWAS Summary Statistics
#'
#' @description Computes the approximation F = (beta_ZX / se_ZX)^2 for each SNP.
#'   This is the standard approximation used when individual-level data is not
#'   available.
#'
#' @param betaZX Numeric vector of SNP-exposure effect estimates.
#' @param seZX Numeric vector of standard errors for SNP-exposure effects.
#'
#' @return Numeric vector of approximate F-statistics.
#'
#' @examples
#' computeApproxFStatistic(c(0.5, 0.3), c(0.05, 0.1))
#'
#' @export
computeApproxFStatistic <- function(betaZX, seZX) {
  checkmate::assertNumeric(betaZX, any.missing = FALSE)
  checkmate::assertNumeric(seZX, lower = 0, any.missing = FALSE,
                           len = length(betaZX))
  (betaZX / seZX)^2
}


#' Build the canonical cohort-data column name for an instrument SNP
#'
#' @param snpId Character vector of SNP identifiers.
#'
#' @return Character vector of column names in the form \code{snp_<sanitized_id>}.
#'
#' @keywords internal
makeSnpColumnName <- function(snpId) {
  checkmate::assertCharacter(snpId, any.missing = FALSE)
  paste0("snp_", gsub("[^a-zA-Z0-9]", "_", snpId))
}


#' Check for Strand-Ambiguous SNPs
#'
#' @description Identifies SNPs with ambiguous strand (A/T or G/C allele pairs).
#'
#' @param effectAllele Character vector of effect alleles.
#' @param otherAllele Character vector of other alleles.
#'
#' @return Logical vector indicating which SNPs are strand-ambiguous.
#'
#' @examples
#' isStrandAmbiguous(c("A", "G", "A"), c("T", "C", "C"))
#'
#' @export
isStrandAmbiguous <- function(effectAllele, otherAllele) {
  checkmate::assertCharacter(effectAllele, any.missing = FALSE)
  checkmate::assertCharacter(otherAllele, any.missing = FALSE,
                             len = length(effectAllele))
  ea <- toupper(effectAllele)
  oa <- toupper(otherAllele)
  pair <- paste0(ea, oa)
  pair %in% c("AT", "TA", "GC", "CG")
}


#' Load and Render SQL from Package
#'
#' @description Convenience wrapper around SqlRender functions that loads SQL
#'   from the package inst/sql directory, renders parameters, and translates
#'   to the target dialect.
#'
#' @param sqlFileName Name of the SQL file in inst/sql/sql_server/.
#' @param dbms Target database management system (e.g., "postgresql", "sql server").
#' @param ... Named parameters to substitute into the SQL template.
#'
#' @return Translated SQL string.
#'
#' @keywords internal
loadRenderTranslateSql <- function(sqlFileName, dbms, ...) {
  if (!requireNamespace("SqlRender", quietly = TRUE)) {
    stop("Package 'SqlRender' is required for SQL operations. ",
         "Install it with: remotes::install_github('ohdsi/SqlRender')",
         call. = FALSE)
  }

  sqlPath <- system.file("sql", "sql_server", sqlFileName, package = "Medusa")
  if (identical(sqlPath, "")) {
    sqlPath <- file.path("inst", "sql", "sql_server", sqlFileName)
  }
  if (!file.exists(sqlPath)) {
    stop(sprintf("SQL file not found: %s", sqlFileName), call. = FALSE)
  }

  sqlTemplate <- paste(readLines(sqlPath, warn = FALSE), collapse = "\n")
  renderedSql <- SqlRender::render(sqlTemplate, ...)
  translatedSql <- SqlRender::translate(renderedSql, targetDialect = dbms)

  translatedSql
}
