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

#' Launch Interactive Results Explorer
#'
#' @title Shiny dashboard for MR diagnostics and results
#'
#' @description
#' Opens an interactive bslib dashboard displaying instrument diagnostics
#' (F-statistics, PheWAS, allele frequencies, missingness, heterogeneity)
#' and MR results (primary estimate, sensitivity analyses, scatter/funnel/
#' leave-one-out plots). Styled to match the Medusa package website.
#'
#' @param mrEstimate Output of \code{\link{computeMREstimate}}.
#' @param sensitivityResults Output of \code{\link{runSensitivityAnalyses}}.
#'   Can be NULL.
#' @param diagnosticResults Output of \code{\link{runInstrumentDiagnostics}}.
#'   Can be NULL.
#' @param combinedProfile Output of \code{\link{poolLikelihoodProfiles}}.
#' @param siteProfileList Optional named list of site profile objects from
#'   \code{\link{fitOutcomeModel}}.
#' @param instrumentTable Output of \code{\link{getMRInstruments}}.
#' @param perSnpEstimates Optional data frame of per-SNP summary statistics
#'   (output of \code{\link{fitOutcomeModel}} with \code{analysisType = "perSNP"}).
#'   Needed for scatter and funnel plots.
#' @param negativeControlResults Output of
#'   \code{\link{runNegativeControlAnalysis}}. Can be NULL.
#' @param launch.browser Logical. Whether to open a browser window. Default
#'   is TRUE.
#' @param port Integer or NULL. Port for the Shiny server. NULL uses a random
#'   available port.
#'
#' @return Invisible NULL. Launches a Shiny app.
#'
#' @examples
#' \dontrun{
#' launchResultsExplorer(
#'   mrEstimate = estimate,
#'   sensitivityResults = sensitivity,
#'   diagnosticResults = diagnostics,
#'   combinedProfile = combined,
#'   instrumentTable = instruments,
#'   perSnpEstimates = perSnpEst
#' )
#' }
#'
#' @seealso \code{\link{computeMREstimate}}, \code{\link{runSensitivityAnalyses}},
#'   \code{\link{runInstrumentDiagnostics}}, \code{\link{generateMRReport}}
#'
#' @export
launchResultsExplorer <- function(mrEstimate,
                                   sensitivityResults = NULL,
                                   diagnosticResults = NULL,
                                   combinedProfile = NULL,
                                   siteProfileList = NULL,
                                   instrumentTable = NULL,
                                   perSnpEstimates = NULL,
                                   negativeControlResults = NULL,
                                   launch.browser = TRUE,
                                   port = NULL) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required for launchResultsExplorer(). ",
         "Install it with: install.packages('shiny')",
         call. = FALSE)
  }
  if (!requireNamespace("bslib", quietly = TRUE)) {
    stop("Package 'bslib' is required for launchResultsExplorer(). ",
         "Install it with: install.packages('bslib')",
         call. = FALSE)
  }
  if (!requireNamespace("bsicons", quietly = TRUE)) {
    stop("Package 'bsicons' is required for launchResultsExplorer(). ",
         "Install it with: install.packages('bsicons')",
         call. = FALSE)
  }

  # Use combined profile from mrEstimate if not supplied

  if (is.null(combinedProfile) && !is.null(mrEstimate$combinedProfile)) {
    combinedProfile <- mrEstimate$combinedProfile
  }

  medusaResults <- list(
    mrEstimate         = mrEstimate,
    sensitivityResults = sensitivityResults,
    diagnosticResults  = diagnosticResults,
    combinedProfile    = combinedProfile,
    siteProfileList    = siteProfileList,
    instrumentTable          = instrumentTable,
    perSnpEstimates          = perSnpEstimates,
    negativeControlResults   = negativeControlResults
  )

  appDir <- system.file("shiny", "MedusaExplorer", package = "Medusa")
  if (appDir == "") {
    # Fallback for development
    appDir <- file.path("inst", "shiny", "MedusaExplorer")
  }
  if (!dir.exists(appDir)) {
    stop("Cannot find MedusaExplorer Shiny app directory.", call. = FALSE)
  }

  shiny::shinyOptions(medusaResults = medusaResults)
  shiny::runApp(
    appDir,
    launch.browser = launch.browser,
    port = port,
    display.mode = "normal"
  )
}
