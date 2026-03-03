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

#' Validate a Single-Site Medusa Estimate Against TwoSampleMR
#'
#' @title Single-site cross-check against standard summary-data MR
#'
#' @description
#' Runs Medusa's single-site profile-likelihood estimate and compares it with
#' \pkg{TwoSampleMR} estimates computed from the same site's per-SNP summary
#' statistics. This is a calibration check: agreement is expected when the
#' profile is close to quadratic and the allele-score approximation aligns with
#' standard summarized-data MR assumptions.
#'
#' @param siteProfile A single \code{medusaSiteProfile} object produced by
#'   \code{\link{fitOutcomeModel}} with \code{analysisType = "perSNP"}.
#' @param instrumentTable Data frame. Output of \code{\link{getMRInstruments}}.
#'   Used both for the Medusa Wald ratio and, if needed, to backfill allele
#'   columns missing from older \code{perSnpEstimates} objects.
#' @param methods Character vector of summary-data methods to compare. Default
#'   is \code{c("IVW", "MREgger", "WeightedMedian")}.
#' @param outcomeSampleSize Optional integer. Passed to
#'   \code{\link{runSensitivityAnalyses}}.
#' @param exposureSampleSize Optional integer. Passed to
#'   \code{\link{runSensitivityAnalyses}}.
#' @param outcomeType Character. Outcome summary-statistic scale represented by
#'   \code{siteProfile$perSnpEstimates}. Default is \code{"binary"}.
#' @param ciLevel Numeric. Confidence interval level for the Medusa estimate.
#'   Default is 0.95.
#'
#' @return A data frame with one row for the Medusa profile-likelihood estimate
#'   and one row for each requested \pkg{TwoSampleMR} method. The
#'   \code{delta_vs_medusa} column reports each summary-data estimate minus the
#'   Medusa estimate.
#'
#' @examples
#' simData <- simulateMRData(n = 3000, nSnps = 5, trueEffect = 0.3)
#' siteProfile <- fitOutcomeModel(
#'   cohortData = simData$data,
#'   instrumentTable = simData$instrumentTable,
#'   analysisType = "perSNP",
#'   betaGrid = seq(-2, 2, by = 0.05)
#' )
#' comparison <- validateAgainstTwoSampleMR(siteProfile, simData$instrumentTable)
#' comparison
#'
#' @seealso \code{\link{fitOutcomeModel}}, \code{\link{computeMREstimate}},
#'   \code{\link{runSensitivityAnalyses}}
#'
#' @export
validateAgainstTwoSampleMR <- function(siteProfile,
                                       instrumentTable,
                                       methods = c("IVW", "MREgger", "WeightedMedian"),
                                       outcomeSampleSize = NULL,
                                       exposureSampleSize = NULL,
                                       outcomeType = "binary",
                                       ciLevel = 0.95) {
  validateSiteProfile(siteProfile)
  validateInstrumentTable(instrumentTable)
  checkmate::assertSubset(methods,
                          c("IVW", "MREgger", "WeightedMedian",
                            "Steiger", "LeaveOneOut"))
  checkmate::assertChoice(outcomeType, c("binary", "continuous"))
  checkmate::assertNumber(ciLevel, lower = 0.5, upper = 0.999)

  if (is.null(siteProfile$perSnpEstimates)) {
    stop(
      "siteProfile must contain perSnpEstimates. Re-run fitOutcomeModel() with analysisType = 'perSNP'."
    )
  }

  perSnpEstimates <- ensureHarmonisationColumns(
    perSnpEstimates = siteProfile$perSnpEstimates,
    instrumentTable = instrumentTable
  )

  combinedProfile <- poolLikelihoodProfiles(list(siteProfile))
  medusaEstimate <- computeMREstimate(
    combinedProfile = combinedProfile,
    instrumentTable = instrumentTable,
    ciLevel = ciLevel
  )
  sensitivityResults <- runSensitivityAnalyses(
    perSnpEstimates = perSnpEstimates,
    methods = methods,
    outcomeSampleSize = outcomeSampleSize,
    exposureSampleSize = exposureSampleSize,
    outcomeType = outcomeType,
    engine = "TwoSampleMR"
  )

  medusaRow <- data.frame(
    source = "Medusa",
    method = "Profile likelihood",
    beta_MR = medusaEstimate$betaMR,
    se_MR = medusaEstimate$seMR,
    ci_lower = medusaEstimate$ciLower,
    ci_upper = medusaEstimate$ciUpper,
    pval = medusaEstimate$pValue,
    delta_vs_medusa = 0,
    stringsAsFactors = FALSE
  )

  summaryRows <- sensitivityResults$summary
  if (!is.null(summaryRows) && nrow(summaryRows) > 0) {
    summaryRows <- data.frame(
      source = "TwoSampleMR",
      method = summaryRows$method,
      beta_MR = summaryRows$beta_MR,
      se_MR = summaryRows$se_MR,
      ci_lower = summaryRows$ci_lower,
      ci_upper = summaryRows$ci_upper,
      pval = summaryRows$pval,
      delta_vs_medusa = summaryRows$beta_MR - medusaEstimate$betaMR,
      stringsAsFactors = FALSE
    )
  }

  comparison <- rbind(
    medusaRow,
    if (!is.null(summaryRows) && nrow(summaryRows) > 0) summaryRows
  )
  rownames(comparison) <- NULL
  comparison
}


#' @keywords internal
ensureHarmonisationColumns <- function(perSnpEstimates,
                                       instrumentTable) {
  requiredCols <- c("effect_allele", "other_allele", "eaf")
  if (all(requiredCols %in% names(perSnpEstimates))) {
    return(perSnpEstimates)
  }

  matchIdx <- match(perSnpEstimates$snp_id, instrumentTable$snp_id)
  if (anyNA(matchIdx)) {
    stop(
      "perSnpEstimates are missing harmonisation columns and could not be matched to instrumentTable by snp_id."
    )
  }

  enriched <- perSnpEstimates
  for (col in requiredCols) {
    if (!col %in% names(enriched)) {
      enriched[[col]] <- instrumentTable[[col]][matchIdx]
    }
  }
  if (!"pval_ZX" %in% names(enriched) && "pval_ZX" %in% names(instrumentTable)) {
    enriched$pval_ZX <- instrumentTable$pval_ZX[matchIdx]
  }

  enriched
}
