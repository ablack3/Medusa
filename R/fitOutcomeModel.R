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

#' Fit Outcome Model and Evaluate Profile Log-Likelihood
#'
#' @title Outcome model with exact grid likelihood evaluation
#'
#' @description
#' Fits a logistic regression of the binary outcome on SNP genotype(s) plus
#' covariates, then evaluates the profile log-likelihood across a pre-specified
#' grid of beta_ZY values. This is the core methodological function for
#' federated MR: each site runs this locally and shares only the resulting
#' log-likelihood profile vector.
#'
#' Two analysis modes are supported: \code{alleleScore} fits a single model
#' using a weighted allele score as the genetic exposure variable, while
#' \code{perSNP} fits separate models for each SNP (needed for multi-SNP
#' sensitivity analyses).
#'
#' @param cohortData Data frame. Output of \code{\link{buildMRCohort}}.
#' @param covariateData Covariate data. Output of \code{\link{buildMRCovariates}}
#'   or a data frame with person_id and covariate columns.
#' @param instrumentTable Data frame. Output of \code{\link{getMRInstruments}}.
#' @param betaGrid Numeric vector. Grid of beta_ZY values at which to evaluate
#'   the profile log-likelihood. Default is \code{seq(-3, 3, by = 0.01)} (601
#'   grid points).
#' @param regularizationVariance Numeric. Prior variance for Cyclops
#'   regularization of covariate coefficients. Default is 0.1.
#' @param instrumentRegularization Logical. If FALSE (default), the SNP/allele
#'   score coefficient is NOT regularized (no shrinkage). Covariates are still
#'   regularized.
#' @param outcomeType Character. Type of outcome: "binary" for logistic
#'   regression. Default is "binary".
#' @param analysisType Character. "alleleScore" for single weighted score model
#'   or "perSNP" for separate per-SNP models. Default is "alleleScore".
#' @param siteId Character. Identifier for this site. Included in the returned
#'   profile object. Default is "site_1".
#'
#' @return A list with class "medusaSiteProfile" containing:
#'   \describe{
#'     \item{siteId}{Character site identifier.}
#'     \item{betaGrid}{Numeric vector of grid points (same as input).}
#'     \item{logLikProfile}{Numeric vector of profile log-likelihood values,
#'       same length as betaGrid. In \code{perSNP} mode this remains the valid
#'       allele-score profile used for pooling.}
#'     \item{nCases}{Number of outcome cases.}
#'     \item{nControls}{Number of controls.}
#'     \item{snpIds}{Character vector of SNP IDs used.}
#'     \item{diagnosticFlags}{List of diagnostic flags from the model fit.}
#'     \item{betaHat}{Point estimate of beta_ZY (MLE from profile).}
#'     \item{seHat}{Approximate standard error from profile curvature.}
#'     \item{perSnpEstimates}{When \code{analysisType = "perSNP"}, a data frame
#'       of per-SNP beta_ZY / se_ZY estimates for summarized-data sensitivity
#'       analyses.}
#'   }
#'   This object contains no individual-level data and is safe to share.
#'
#' @details
#' The profile log-likelihood evaluation works as follows:
#' \enumerate{
#'   \item Fit the unconstrained model to obtain the MLE and its Wald SE.
#'   \item At each grid point, fix beta_ZY to that value and re-fit the
#'     nuisance parameters using the fixed term as an offset.
#'   \item Record the maximized constrained log-likelihood at each grid point.
#'     This is the exact profile likelihood for the coefficient in a logistic
#'     generalized linear model.
#' }
#'
#' When using the allele score, weights are
#' \eqn{w_j = \gamma_j / \mathrm{SE}(\gamma_j)^2}, normalized by
#' \eqn{\sum_j |w_j|}. The same weights are reused downstream so that the
#' MR denominator matches the exact score fitted at each site.
#'
#' @references
#' Suchard, M. A., et al. (2013). Massive parallelization of serial inference
#' algorithms for a complex generalized linear model. \emph{ACM Transactions
#' on Modeling and Computer Simulation}, 23(1), 1-17.
#'
#' @examples
#' simData <- simulateMRData(n = 2000, nSnps = 5, trueEffect = 0.3)
#' profile <- fitOutcomeModel(
#'   cohortData = simData$data,
#'   covariateData = NULL,
#'   instrumentTable = simData$instrumentTable,
#'   betaGrid = seq(-2, 2, by = 0.05)
#' )
#' plot(profile$betaGrid, profile$logLikProfile, type = "l",
#'      xlab = "beta_ZY", ylab = "Profile log-likelihood")
#'
#' @seealso \code{\link{poolLikelihoodProfiles}}, \code{\link{computeMREstimate}},
#'   \code{\link{buildMRCohort}}
#'
#' @export
fitOutcomeModel <- function(cohortData,
                            covariateData = NULL,
                            instrumentTable,
                            betaGrid = seq(-3, 3, by = 0.01),
                            regularizationVariance = 0.1,
                            instrumentRegularization = FALSE,
                            outcomeType = "binary",
                            analysisType = "alleleScore",
                            siteId = "site_1") {
  # Input validation
  checkmate::assertDataFrame(cohortData, min.rows = 1)
  validateInstrumentTable(instrumentTable)
  validateBetaGrid(betaGrid)
  checkmate::assertNumber(regularizationVariance, lower = 0)
  checkmate::assertLogical(instrumentRegularization, len = 1)
  checkmate::assertChoice(outcomeType, c("binary"))
  checkmate::assertChoice(analysisType, c("alleleScore", "perSNP"))
  checkmate::assertString(siteId)

  checkmate::assertSubset("outcome", names(cohortData))
  snpCols <- grep("^snp_", names(cohortData), value = TRUE)
  if (length(snpCols) == 0) {
    stop("No SNP columns (snp_*) found in cohortData.")
  }

  nCases <- sum(cohortData$outcome == 1, na.rm = TRUE)
  nControls <- sum(cohortData$outcome == 0, na.rm = TRUE)

  if (nCases == 0) {
    stop("No outcome cases found in cohortData.")
  }
  if (nCases < 50) {
    warning(sprintf(
      "Only %d cases in outcome cohort. Results may be unstable. Consider expanding cohort definition or adding sites.",
      nCases
    ))
  }

  message(sprintf("Fitting outcome model at site '%s' (%d cases, %d controls)...",
                  siteId, nCases, nControls))

  if (analysisType == "alleleScore") {
    result <- fitAlleleScoreModel(cohortData, covariateData, instrumentTable,
                                   betaGrid, regularizationVariance,
                                   instrumentRegularization)
  } else {
    result <- fitPerSNPModels(cohortData, covariateData, instrumentTable,
                               betaGrid, regularizationVariance,
                               instrumentRegularization)
  }

  # Check for grid boundary MLE
  gridBoundaryMLE <- FALSE
  if (is.numeric(result$logLikProfile) && length(result$logLikProfile) > 0) {
    peakIdx <- which.max(result$logLikProfile)
    if (peakIdx == 1 || peakIdx == length(result$logLikProfile)) {
      gridBoundaryMLE <- TRUE
      warning(sprintf(
        "MLE is at grid boundary (%.2f). Expand betaGrid range. Current range: [%.2f, %.2f].",
        betaGrid[peakIdx], min(betaGrid), max(betaGrid)
      ))
    }
  }

  # Check for flat profile
  if (is.numeric(result$logLikProfile)) {
    llRange <- diff(range(result$logLikProfile, na.rm = TRUE))
    if (llRange < 0.5) {
      warning("Profile likelihood is flat -- instrument may be too weak to estimate beta_ZY. Check F-statistic.")
    }
  }

  profile <- list(
    siteId = siteId,
    betaGrid = betaGrid,
    logLikProfile = result$logLikProfile,
    nCases = nCases,
    nControls = nControls,
    snpIds = instrumentTable$snp_id,
    diagnosticFlags = list(
      weakInstruments = any(instrumentTable$fStatistic < 10, na.rm = TRUE),
      lowCaseCount = nCases < 50,
      gridBoundaryMLE = gridBoundaryMLE
    ),
    betaHat = result$betaHat,
    seHat = result$seHat
  )
  if (!is.null(result$perSnpEstimates)) {
    profile$perSnpEstimates <- result$perSnpEstimates
  }
  class(profile) <- "medusaSiteProfile"

  message(sprintf("Site '%s': beta_ZY_hat = %.4f (SE = %.4f).",
                  siteId, result$betaHat, result$seHat))

  profile
}


#' @keywords internal
fitAlleleScoreModel <- function(cohortData, covariateData, instrumentTable,
                                 betaGrid, regularizationVariance,
                                 instrumentRegularization) {
  snpCols <- grep("^snp_", names(cohortData), value = TRUE)
  nSnps <- min(length(snpCols), nrow(instrumentTable))

  # Burgess et al. (2013): the allele score is a weighted sum of genotypes,
  # using external SNP-exposure effects to define fixed score weights.
  weights <- computeAlleleScoreWeights(
    instrumentTable[seq_len(nSnps), , drop = FALSE]
  )

  snpMatrix <- as.matrix(cohortData[, snpCols[seq_len(nSnps)], drop = FALSE])
  # Replace NA with 0 for score computation (common imputation for missing genotypes)
  snpMatrixImputed <- snpMatrix
  snpMatrixImputed[is.na(snpMatrixImputed)] <- 0

  alleleScore <- as.numeric(snpMatrixImputed %*% weights)

  # Build model data frame
  modelData <- data.frame(
    outcome = cohortData$outcome,
    alleleScore = alleleScore
  )
  modelParts <- appendCovariatesToModelData(modelData, cohortData, covariateData)
  modelData <- modelParts$modelData
  covCols <- modelParts$covariateColumns

  # Use one common complete-case sample for the unconstrained model and every
  # constrained fit. Otherwise the profile compares likelihoods on different
  # patient sets.
  completeIdx <- complete.cases(modelData)
  modelData <- modelData[completeIdx, , drop = FALSE]

  # Fit unconstrained logistic regression
  formula <- as.formula(paste("outcome ~ alleleScore",
                               if (length(covCols) > 0) paste("+", paste(covCols, collapse = " + ")) else ""))

  tryCatch({
    fit <- glm(formula, data = modelData, family = binomial())
    betaHat <- coef(fit)["alleleScore"]
    seHat <- summary(fit)$coefficients["alleleScore", "Std. Error"]

    # Exact profile likelihood for a fixed coefficient:
    #   l_p(beta) = max_alpha l(alpha, beta)
    # We implement this by re-fitting the nuisance parameters with
    # beta * alleleScore supplied as an offset.
    logLikProfile <- evaluateBinaryProfile(
      modelData = modelData,
      exposureColumn = "alleleScore",
      covariateColumns = covCols,
      betaGrid = betaGrid
    )

    list(
      logLikProfile = logLikProfile,
      betaHat = unname(betaHat),
      seHat = unname(seHat)
    )
  }, error = function(e) {
    warning(sprintf("Model fitting failed: %s. Returning flat profile.", conditionMessage(e)))
    list(
      logLikProfile = rep(0, length(betaGrid)),
      betaHat = 0,
      seHat = Inf
    )
  })
}


#' @keywords internal
fitPerSNPModels <- function(cohortData, covariateData, instrumentTable,
                             betaGrid, regularizationVariance,
                             instrumentRegularization) {
  snpCols <- grep("^snp_", names(cohortData), value = TRUE)
  nSnps <- min(length(snpCols), nrow(instrumentTable))

  # A site can provide per-SNP summary estimates for IVW / MR-Egger / weighted
  # median, but those single-SNP regressions do not form a joint likelihood that
  # can be summed across SNPs. Keep the pooled site profile tied to the valid
  # allele-score model and expose the per-SNP fits as auxiliary outputs only.
  alleleScoreFit <- fitAlleleScoreModel(
    cohortData = cohortData,
    covariateData = covariateData,
    instrumentTable = instrumentTable,
    betaGrid = betaGrid,
    regularizationVariance = regularizationVariance,
    instrumentRegularization = instrumentRegularization
  )

  profileMatrix <- matrix(NA_real_, nrow = length(betaGrid), ncol = nSnps)
  betaHats <- numeric(nSnps)
  seHats <- numeric(nSnps)

  for (j in seq_len(nSnps)) {
    snpCol <- snpCols[j]
    modelData <- data.frame(
      outcome = cohortData$outcome,
      snp = cohortData[[snpCol]]
    )
    modelParts <- appendCovariatesToModelData(modelData, cohortData, covariateData)
    modelData <- modelParts$modelData
    covCols <- modelParts$covariateColumns

    completeIdx <- complete.cases(modelData)
    modelData <- modelData[completeIdx, , drop = FALSE]

    if (nrow(modelData) < 20 || length(unique(modelData$outcome)) < 2) {
      profileMatrix[, j] <- rep(-Inf, length(betaGrid))
      betaHats[j] <- NA_real_
      seHats[j] <- NA_real_
      next
    }

    formula <- as.formula(paste("outcome ~ snp",
                                 if (length(covCols) > 0) paste("+", paste(covCols, collapse = " + ")) else ""))

    tryCatch({
      fit <- glm(formula, data = modelData, family = binomial())
      betaHat <- coef(fit)["snp"]
      seHat <- summary(fit)$coefficients["snp", "Std. Error"]

      profileMatrix[, j] <- evaluateBinaryProfile(
        modelData = modelData,
        exposureColumn = "snp",
        covariateColumns = covCols,
        betaGrid = betaGrid
      )
      betaHats[j] <- unname(betaHat)
      seHats[j] <- unname(seHat)
    }, error = function(e) {
      profileMatrix[, j] <<- rep(-Inf, length(betaGrid))
      betaHats[j] <<- NA_real_
      seHats[j] <<- NA_real_
    })
  }

  colnames(profileMatrix) <- instrumentTable$snp_id[seq_len(nSnps)]
  perSnpEstimates <- data.frame(
    snp_id = instrumentTable$snp_id[seq_len(nSnps)],
    beta_ZY = betaHats,
    se_ZY = seHats,
    beta_ZX = instrumentTable$beta_ZX[seq_len(nSnps)],
    se_ZX = instrumentTable$se_ZX[seq_len(nSnps)],
    stringsAsFactors = FALSE
  )

  list(
    logLikProfile = alleleScoreFit$logLikProfile,
    perSnpProfiles = profileMatrix,
    betaHat = alleleScoreFit$betaHat,
    seHat = alleleScoreFit$seHat,
    perSnpBetaHats = betaHats,
    perSnpSEHats = seHats,
    perSnpEstimates = perSnpEstimates[stats::complete.cases(perSnpEstimates), , drop = FALSE]
  )
}


#' @keywords internal
appendCovariatesToModelData <- function(modelData, cohortData, covariateData) {
  covCols <- character(0)

  if (!is.null(covariateData) && is.data.frame(covariateData)) {
    covCols <- setdiff(names(covariateData), c("person_id", "personId", "rowId"))

    if ("personId" %in% names(covariateData) && "personId" %in% names(cohortData)) {
      mergedCov <- dplyr::left_join(
        data.frame(personId = cohortData$personId),
        covariateData,
        by = "personId"
      )
    } else {
      mergedCov <- covariateData
    }

    for (col in covCols) {
      if (col %in% names(mergedCov)) {
        modelData[[col]] <- mergedCov[[col]]
      }
    }
  }

  for (confCol in c("confounder_1", "confounder_2")) {
    if (confCol %in% names(cohortData)) {
      modelData[[confCol]] <- cohortData[[confCol]]
      covCols <- c(covCols, confCol)
    }
  }

  list(
    modelData = modelData,
    covariateColumns = unique(covCols)
  )
}


#' @keywords internal
evaluateBinaryProfile <- function(modelData,
                                  exposureColumn,
                                  covariateColumns,
                                  betaGrid) {
  checkmate::assertDataFrame(modelData, min.rows = 1)
  checkmate::assertString(exposureColumn)
  if (!is.character(covariateColumns)) {
    stop("covariateColumns must be a character vector.")
  }

  profileFormula <- if (length(covariateColumns) > 0) {
    as.formula(paste("outcome ~", paste(covariateColumns, collapse = " + ")))
  } else {
    stats::as.formula("outcome ~ 1")
  }

  response <- modelData$outcome
  designMatrix <- stats::model.matrix(profileFormula, data = modelData)
  exposure <- modelData[[exposureColumn]]
  logLikProfile <- vapply(betaGrid, function(betaFixed) {
    fit <- tryCatch(
      suppressWarnings(
        stats::glm.fit(
          x = designMatrix,
          y = response,
          family = stats::binomial(),
          offset = betaFixed * exposure
        )
      ),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      return(-Inf)
    }

    fittedProb <- pmin(pmax(fit$fitted.values, 1e-12), 1 - 1e-12)
    fitLogLik <- sum(stats::dbinom(response, size = 1, prob = fittedProb, log = TRUE))
    if (!is.finite(fitLogLik)) {
      return(-Inf)
    }

    fitLogLik
  }, numeric(1))

  if (all(!is.finite(logLikProfile))) {
    stop("Profile likelihood failed at all grid points.")
  }

  logLikProfile
}


#' @keywords internal
computeAlleleScoreWeights <- function(instrumentTable) {
  rawWeights <- instrumentTable$beta_ZX / (instrumentTable$se_ZX^2)
  scalingConstant <- sum(abs(rawWeights))

  if (!is.finite(scalingConstant) || scalingConstant <= 0) {
    stop("Instrument weights are not finite. Check beta_ZX and se_ZX.")
  }

  rawWeights / scalingConstant
}


#' Estimate Standard Error from Profile Log-Likelihood Curvature
#'
#' @description Uses the second derivative of the profile log-likelihood at the
#'   MLE to estimate the standard error.
#'
#' @param betaGrid Numeric vector of grid points.
#' @param logLikProfile Numeric vector of log-likelihood values.
#'
#' @return Numeric standard error estimate.
#'
#' @keywords internal
estimateSEFromProfile <- function(betaGrid, logLikProfile) {
  peakIdx <- which.max(logLikProfile)
  if (peakIdx <= 1 || peakIdx >= length(betaGrid)) {
    return(Inf)
  }

  # Finite difference approximation of second derivative at peak
  gridStep <- betaGrid[peakIdx + 1] - betaGrid[peakIdx]
  d2 <- (logLikProfile[peakIdx + 1] - 2 * logLikProfile[peakIdx] +
           logLikProfile[peakIdx - 1]) / (gridStep^2)

  if (d2 >= 0) {
    return(Inf)
  }

  sqrt(-1 / d2)
}
