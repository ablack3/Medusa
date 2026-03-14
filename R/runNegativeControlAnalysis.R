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

#' Run Negative Control Outcome Analysis for Empirical Calibration
#'
#' @title Empirical method validation using negative control outcomes
#'
#' @description
#' Implements OHDSI-style empirical calibration for MR estimates. For each
#' negative control outcome (an outcome with no expected causal relationship
#' to the exposure), the same allele score model is fitted to obtain a
#' null MR estimate. The distribution of these null estimates is then used
#' to assess systematic bias and, optionally, to calibrate the primary
#' estimate's p-value and confidence interval via the \pkg{EmpiricalCalibration}
#' package.
#'
#' @param cohortData Data frame. Output of \code{\link{buildMRCohort}}, expected
#'   to contain negative control outcome columns (named \code{nc_outcome_<id>}
#'   or as specified in \code{negativeControlColumns}).
#' @param instrumentTable Data frame. Output of \code{\link{getMRInstruments}}.
#' @param covariateData Optional covariate data for adjustment. Same format as
#'   accepted by \code{\link{fitOutcomeModel}}.
#' @param negativeControlColumns Optional character vector of column names in
#'   \code{cohortData} that contain negative control outcomes. If NULL,
#'   auto-detects columns matching \code{nc_outcome_*}.
#' @param primaryEstimate Optional list. Output of
#'   \code{\link{computeMREstimate}}. If provided, the primary estimate is
#'   calibrated using the negative control distribution.
#' @param modelBackend Character. Model fitting backend: \code{"glm"} or
#'   \code{"cyclops"}. Default is \code{"glm"}.
#'
#' @return A list with class \code{"medusaNegativeControls"} containing:
#'   \describe{
#'     \item{ncEstimates}{Data frame with one row per negative control outcome:
#'       \code{outcome_id}, \code{beta_ZY}, \code{se_ZY}, \code{beta_MR},
#'       \code{se_MR}, \code{pval}, \code{log_rr}, \code{se_log_rr}.}
#'     \item{calibration}{Output of
#'       \code{EmpiricalCalibration::fitSystematicErrorModel()} if the package is
#'       available, otherwise NULL.}
#'     \item{calibratedPrimary}{Named list with \code{calibratedP},
#'       \code{calibratedCiLower}, \code{calibratedCiUpper} for the primary
#'       estimate, or NULL if \code{primaryEstimate} was not provided or
#'       calibration is unavailable.}
#'     \item{biasDetected}{Logical. TRUE if the systematic error model indicates
#'       the null distribution mean is significantly different from zero.}
#'   }
#'
#' @details
#' Negative control outcomes should be conditions that are biologically
#' implausible as consequences of the exposure. A well-calibrated MR analysis
#' should produce null estimates for these outcomes. Systematic deviation from
#' the null suggests residual pleiotropy, population stratification, or other
#' bias.
#'
#' When \pkg{EmpiricalCalibration} is installed, the function fits a systematic
#' error model to the negative control estimates and uses it to produce
#' calibrated p-values and confidence intervals for the primary estimate.
#' This approach is described in Schuemie et al. (2014, 2018) and is standard
#' practice in OHDSI observational studies.
#'
#' @references
#' Schuemie, M. J., et al. (2014). Interpreting observational studies: why
#' empirical calibration is needed to correct p-values. \emph{Statistics in
#' Medicine}, 33(2), 209-218.
#'
#' Schuemie, M. J., et al. (2018). Empirical confidence interval calibration
#' for population-level effect estimation studies in observational healthcare
#' data. \emph{PNAS}, 115(11), 2571-2577.
#'
#' @examples
#' simData <- simulateMRData(n = 2000, nSnps = 5, seed = 123)
#' cohortData <- simulateNegativeControlOutcomes(simData$data, nControls = 10)
#' ncResults <- runNegativeControlAnalysis(
#'   cohortData = cohortData,
#'   instrumentTable = simData$instrumentTable
#' )
#' ncResults$ncEstimates
#'
#' @seealso \code{\link{runInstrumentDiagnostics}},
#'   \code{\link{computeMREstimate}}, \code{\link{buildMRCohort}}
#'
#' @export
runNegativeControlAnalysis <- function(cohortData,
                                       instrumentTable,
                                       covariateData = NULL,
                                       negativeControlColumns = NULL,
                                       primaryEstimate = NULL,
                                       modelBackend = "glm") {
  checkmate::assertDataFrame(cohortData, min.rows = 1)
  validateInstrumentTable(instrumentTable)
  checkmate::assertChoice(modelBackend, c("glm", "cyclops"))

  # Identify NC outcome columns
  if (is.null(negativeControlColumns)) {
    negativeControlColumns <- grep("^nc_outcome_", names(cohortData), value = TRUE)
  }
  checkmate::assertCharacter(negativeControlColumns, min.len = 1)
  checkmate::assertSubset(negativeControlColumns, names(cohortData))

  nNC <- length(negativeControlColumns)
  message(sprintf("Running negative control analysis for %d outcomes...", nNC))

  # Build allele score once
  alignment <- alignInstrumentColumns(cohortData, instrumentTable)
  alignedInstruments <- alignment$instrumentTable
  snpCols <- alignment$snpColumns
  scoreWeights <- computeAlleleScoreWeights(alignedInstruments)

  snpMatrix <- as.matrix(cohortData[, snpCols, drop = FALSE])
  if (anyNA(snpMatrix)) {
    expectedDosage <- 2 * alignedInstruments$eaf
    for (j in seq_len(ncol(snpMatrix))) {
      missingMask <- is.na(snpMatrix[, j])
      if (any(missingMask)) {
        snpMatrix[missingMask, j] <- expectedDosage[[j]]
      }
    }
  }
  alleleScore <- as.numeric(snpMatrix %*% scoreWeights)

  # Score-level exposure effect for Wald ratio
  betaZX <- sum(scoreWeights * alignedInstruments$beta_ZX)
  seZX <- sqrt(sum((scoreWeights^2) * (alignedInstruments$se_ZX^2)))

  if (!is.finite(betaZX) || abs(betaZX) < .Machine$double.eps^0.5) {
    stop("Allele score exposure effect is too close to zero for Wald ratio.")
  }

  # Fit model for each NC outcome
  ncResults <- vector("list", nNC)

  for (i in seq_len(nNC)) {
    ncCol <- negativeControlColumns[i]
    outcomeId <- sub("^nc_outcome_", "", ncCol)

    # Build model data with NC outcome
    modelData <- data.frame(
      outcome = cohortData[[ncCol]],
      alleleScore = alleleScore
    )
    modelParts <- appendCovariatesToModelData(modelData, cohortData, covariateData)
    modelData <- modelParts$modelData
    covCols <- modelParts$covariateColumns

    completeIdx <- complete.cases(modelData)
    modelData <- modelData[completeIdx, , drop = FALSE]

    nCases <- sum(modelData$outcome == 1, na.rm = TRUE)

    if (nCases < 5 || nCases == nrow(modelData)) {
      message(sprintf("  NC outcome %s: skipped (cases = %d).", outcomeId, nCases))
      ncResults[[i]] <- data.frame(
        outcome_id = outcomeId,
        beta_ZY = NA_real_, se_ZY = NA_real_,
        beta_MR = NA_real_, se_MR = NA_real_,
        pval = NA_real_,
        log_rr = NA_real_, se_log_rr = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }

    fit <- tryCatch({
      fitBinaryOutcomeCoefficient(
        modelData = modelData,
        exposureColumn = "alleleScore",
        covariateColumns = covCols,
        modelBackend = modelBackend,
        regularizationVariance = Inf,
        instrumentRegularization = FALSE
      )
    }, error = function(e) {
      message(sprintf("  NC outcome %s: model failed (%s).", outcomeId,
                      conditionMessage(e)))
      NULL
    })

    if (is.null(fit)) {
      ncResults[[i]] <- data.frame(
        outcome_id = outcomeId,
        beta_ZY = NA_real_, se_ZY = NA_real_,
        beta_MR = NA_real_, se_MR = NA_real_,
        pval = NA_real_,
        log_rr = NA_real_, se_log_rr = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }

    betaZY <- fit$betaHat
    seZY <- fit$seHat

    # Wald ratio + delta method SE
    betaMR <- betaZY / betaZX
    seMR <- sqrt((seZY / betaZX)^2 + (betaZY * seZX / betaZX^2)^2)
    zStat <- betaMR / seMR
    pval <- 2 * stats::pnorm(-abs(zStat))

    ncResults[[i]] <- data.frame(
      outcome_id = outcomeId,
      beta_ZY = betaZY, se_ZY = seZY,
      beta_MR = betaMR, se_MR = seMR,
      pval = pval,
      log_rr = betaMR, se_log_rr = seMR,
      stringsAsFactors = FALSE
    )
  }

  ncEstimates <- do.call(rbind, ncResults)
  ncEstimates <- ncEstimates[!is.na(ncEstimates$beta_MR), , drop = FALSE]

  nSuccessful <- nrow(ncEstimates)
  message(sprintf("  %d of %d negative control models fitted successfully.",
                  nSuccessful, nNC))

  # Empirical calibration
  calibration <- NULL
  calibratedPrimary <- NULL
  biasDetected <- FALSE

  if (nSuccessful >= 5) {
    biasDetected <- isBiasDetected(ncEstimates)

    if (requireNamespace("EmpiricalCalibration", quietly = TRUE)) {
      calibration <- tryCatch({
        EmpiricalCalibration::fitSystematicErrorModel(
          logRr = ncEstimates$log_rr,
          seLogRr = ncEstimates$se_log_rr
        )
      }, error = function(e) {
        warning(sprintf("EmpiricalCalibration::fitSystematicErrorModel failed: %s",
                        conditionMessage(e)))
        NULL
      })

      if (!is.null(calibration) && !is.null(primaryEstimate)) {
        calibratedPrimary <- tryCatch({
          calP <- EmpiricalCalibration::calibrateP(
            logRr = primaryEstimate$betaMR,
            seLogRr = primaryEstimate$seMR,
            model = calibration
          )
          calCI <- EmpiricalCalibration::calibrateConfidenceInterval(
            logRr = primaryEstimate$betaMR,
            seLogRr = primaryEstimate$seMR,
            model = calibration
          )
          list(
            calibratedP = as.numeric(calP),
            calibratedCiLower = calCI$logLb95Rr,
            calibratedCiUpper = calCI$logUb95Rr
          )
        }, error = function(e) {
          warning(sprintf("Calibration of primary estimate failed: %s",
                          conditionMessage(e)))
          NULL
        })
      }
    } else {
      message("  Package 'EmpiricalCalibration' not installed. ",
              "Skipping systematic error model. ",
              "Install with: remotes::install_github('ohdsi/EmpiricalCalibration')")
    }
  } else if (nSuccessful > 0) {
    message(sprintf(
      "  Only %d successful NC estimates (minimum 5 needed for calibration).",
      nSuccessful
    ))
  }

  result <- list(
    ncEstimates = ncEstimates,
    calibration = calibration,
    calibratedPrimary = calibratedPrimary,
    biasDetected = biasDetected
  )
  class(result) <- "medusaNegativeControls"

  message("Negative control analysis complete.")
  result
}


#' Detect systematic bias from negative control estimates
#'
#' Tests whether the mean of negative control MR estimates is significantly
#' different from zero using a one-sample t-test.
#'
#' @param ncEstimates Data frame of negative control estimates with columns
#'   \code{beta_MR} and \code{se_MR}.
#'
#' @return Logical. TRUE if bias is detected (p < 0.05).
#'
#' @keywords internal
isBiasDetected <- function(ncEstimates) {
  validEstimates <- ncEstimates$beta_MR[is.finite(ncEstimates$beta_MR)]
  if (length(validEstimates) < 3) return(FALSE)

  tTest <- tryCatch(
    stats::t.test(validEstimates, mu = 0),
    error = function(e) NULL
  )
  if (is.null(tTest)) return(FALSE)

  tTest$p.value < 0.05
}


#' Simulate Negative Control Outcome Columns
#'
#' Adds negative control outcome columns to existing cohort data. These
#' outcomes are generated independently of genotype (true null effects)
#' with varying prevalence.
#'
#' @param cohortData Data frame with at least a \code{person_id} column.
#' @param nControls Integer. Number of negative control outcomes to simulate.
#' @param seed Random seed for reproducibility.
#'
#' @return The input \code{cohortData} with appended \code{nc_outcome_<i>}
#'   columns (binary 0/1).
#'
#' @examples
#' simData <- simulateMRData(n = 1000, nSnps = 5)
#' withNC <- simulateNegativeControlOutcomes(simData$data, nControls = 10)
#' head(withNC)
#'
#' @export
simulateNegativeControlOutcomes <- function(cohortData, nControls = 20,
                                             seed = 99) {
  checkmate::assertDataFrame(cohortData, min.rows = 1)
  checkmate::assertCount(nControls, positive = TRUE)

  set.seed(seed)
  n <- nrow(cohortData)

  for (i in seq_len(nControls)) {
    prevalence <- stats::runif(1, 0.03, 0.15)
    cohortData[[paste0("nc_outcome_", i)]] <- stats::rbinom(n, 1, prevalence)
  }

  cohortData
}
