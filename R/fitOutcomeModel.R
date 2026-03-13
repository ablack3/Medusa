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
#' @param regularizationVariance Numeric. When
#'   \code{modelBackend = "cyclops"}, this is the variance of the normal prior
#'   applied to nuisance coefficients. Smaller values imply stronger shrinkage;
#'   use \code{Inf} for an unpenalized Cyclops fit. Ignored by the
#'   \code{"glm"} backend. Default is 0.1.
#' @param instrumentRegularization Logical. When
#'   \code{modelBackend = "cyclops"}, whether the exposure coefficient is
#'   included in the Cyclops prior. Default is FALSE so the profiled exposure
#'   coefficient remains unpenalized. Ignored by the \code{"glm"} backend.
#' @param outcomeType Character. Type of outcome: "binary" for logistic
#'   regression. Default is "binary".
#' @param analysisType Character. "alleleScore" for single weighted score model
#'   or "perSNP" for separate per-SNP models. Default is "alleleScore".
#' @param siteId Character. Identifier for this site. Included in the returned
#'   profile object. Default is "site_1".
#' @param modelBackend Character. Outcome-model fitting backend:
#'   \code{"glm"} uses base R logistic regression, while \code{"cyclops"} uses
#'   \pkg{Cyclops} for scalable logistic regression with optional Gaussian
#'   shrinkage on nuisance covariates. Default is \code{"glm"}.
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
#' When \code{modelBackend = "cyclops"} and \code{regularizationVariance} is
#' finite, the nuisance parameters are estimated under a Gaussian prior at both
#' the unconstrained optimum and every grid point. In that configuration the
#' returned objective is a penalized profile, not an unpenalized MLE profile.
#'
#' When using the allele score, weights are
#' \eqn{w_j = \gamma_j / \mathrm{SE}(\gamma_j)^2}, normalized by
#' \eqn{\sum_j |w_j|}. The same weights are reused downstream so that the
#' MR denominator matches the exact score fitted at each site.
#'
#' Missing SNP dosages in the allele-score model are imputed to the expected
#' dosage \eqn{2 \times \mathrm{EAF}_j} from the instrument table, rather than
#' being treated as homozygous reference. This avoids a systematic downward
#' bias in the score when genotype missingness is present.
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
                            siteId = "site_1",
                            modelBackend = "glm") {
  # Input validation
  checkmate::assertDataFrame(cohortData, min.rows = 1)
  validateInstrumentTable(instrumentTable)
  validateBetaGrid(betaGrid)
  checkmate::assertNumber(regularizationVariance, lower = 0, finite = FALSE)
  checkmate::assertLogical(instrumentRegularization, len = 1)
  checkmate::assertChoice(outcomeType, c("binary"))
  checkmate::assertChoice(analysisType, c("alleleScore", "perSNP"))
  checkmate::assertString(siteId)
  checkmate::assertChoice(modelBackend, c("glm", "cyclops"))

  checkmate::assertSubset("outcome", names(cohortData))
  if (length(grep("^snp_", names(cohortData), value = TRUE)) == 0) {
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
                                   instrumentRegularization, modelBackend)
  } else {
    result <- fitPerSNPModels(cohortData, covariateData, instrumentTable,
                               betaGrid, regularizationVariance,
                               instrumentRegularization, modelBackend)
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
    snpIds = result$alignedInstruments$snp_id,
    diagnosticFlags = list(
      weakInstruments = "fStatistic" %in% names(result$alignedInstruments) &&
        any(result$alignedInstruments$fStatistic < 10, na.rm = TRUE),
      lowCaseCount = nCases < 50,
      gridBoundaryMLE = gridBoundaryMLE
    ),
    betaHat = result$betaHat,
    seHat = result$seHat,
    scoreDefinition = result$scoreDefinition
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
                                 instrumentRegularization,
                                 modelBackend) {
  alignment <- alignInstrumentColumns(cohortData, instrumentTable)
  alignedInstruments <- alignment$instrumentTable
  snpCols <- alignment$snpColumns
  # Burgess et al. (2013): the allele score is a weighted sum of genotypes,
  # using external SNP-exposure effects to define fixed score weights.
  weights <- computeAlleleScoreWeights(alignedInstruments)

  snpMatrix <- as.matrix(cohortData[, snpCols, drop = FALSE])
  # Impute missing additive dosages to their expected allele count (2 * EAF),
  # rather than coding them as homozygous reference (0).
  snpMatrixImputed <- snpMatrix
  if (anyNA(snpMatrixImputed)) {
    expectedDosage <- 2 * alignedInstruments$eaf
    for (j in seq_len(ncol(snpMatrixImputed))) {
      missingMask <- is.na(snpMatrixImputed[, j])
      if (any(missingMask)) {
        snpMatrixImputed[missingMask, j] <- expectedDosage[[j]]
      }
    }
  }

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
  tryCatch({
    coefficientFit <- fitBinaryOutcomeCoefficient(
      modelData = modelData,
      exposureColumn = "alleleScore",
      covariateColumns = covCols,
      modelBackend = modelBackend,
      regularizationVariance = regularizationVariance,
      instrumentRegularization = instrumentRegularization
    )
    betaHat <- coefficientFit$betaHat
    seHat <- coefficientFit$seHat

    # Exact profile likelihood for a fixed coefficient:
    #   l_p(beta) = max_alpha l(alpha, beta)
    # We implement this by re-fitting the nuisance parameters with
    # beta * alleleScore supplied as an offset.
    logLikProfile <- evaluateBinaryProfile(
      modelData = modelData,
      exposureColumn = "alleleScore",
      covariateColumns = covCols,
      betaGrid = betaGrid,
      modelBackend = modelBackend,
      regularizationVariance = regularizationVariance,
      instrumentRegularization = instrumentRegularization
    )

    list(
      logLikProfile = logLikProfile,
      betaHat = unname(betaHat),
      seHat = unname(seHat),
      alignedInstruments = alignedInstruments,
      scoreDefinition = list(
        snpIds = alignedInstruments$snp_id,
        scoreWeights = weights,
        betaZX = sum(weights * alignedInstruments$beta_ZX),
        seZX = sqrt(sum((weights^2) * (alignedInstruments$se_ZX^2)))
      )
    )
  }, error = function(e) {
    warning(sprintf("Model fitting failed: %s. Returning flat profile.", conditionMessage(e)))
    list(
      logLikProfile = rep(0, length(betaGrid)),
      betaHat = 0,
      seHat = Inf,
      alignedInstruments = alignedInstruments,
      scoreDefinition = list(
        snpIds = alignedInstruments$snp_id,
        scoreWeights = weights,
        betaZX = sum(weights * alignedInstruments$beta_ZX),
        seZX = sqrt(sum((weights^2) * (alignedInstruments$se_ZX^2)))
      )
    )
  })
}


#' @keywords internal
fitPerSNPModels <- function(cohortData, covariateData, instrumentTable,
                             betaGrid, regularizationVariance,
                             instrumentRegularization,
                             modelBackend) {
  alignment <- alignInstrumentColumns(cohortData, instrumentTable)
  alignedInstruments <- alignment$instrumentTable
  snpCols <- alignment$snpColumns
  nSnps <- nrow(alignedInstruments)

  # A site can provide per-SNP summary estimates for IVW / MR-Egger / weighted
  # median, but those single-SNP regressions do not form a joint likelihood that
  # can be summed across SNPs. Keep the pooled site profile tied to the valid
  # allele-score model and expose the per-SNP fits as auxiliary outputs only.
  alleleScoreFit <- fitAlleleScoreModel(
    cohortData = cohortData,
    covariateData = covariateData,
    instrumentTable = alignedInstruments,
    betaGrid = betaGrid,
    regularizationVariance = regularizationVariance,
    instrumentRegularization = instrumentRegularization,
    modelBackend = modelBackend
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
      reason <- if (nrow(modelData) < 20) {
        sprintf("too few complete observations (%d < 20)", nrow(modelData))
      } else {
        "outcome has no variation (all cases or all controls)"
      }
      message(sprintf("  Skipping SNP %s: %s.",
                      alignedInstruments$snp_id[j], reason))
      profileMatrix[, j] <- rep(-Inf, length(betaGrid))
      betaHats[j] <- NA_real_
      seHats[j] <- NA_real_
      next
    }

    tryCatch({
      coefficientFit <- fitBinaryOutcomeCoefficient(
        modelData = modelData,
        exposureColumn = "snp",
        covariateColumns = covCols,
        modelBackend = modelBackend,
        regularizationVariance = regularizationVariance,
        instrumentRegularization = instrumentRegularization
      )
      betaHat <- coefficientFit$betaHat
      seHat <- coefficientFit$seHat

      profileMatrix[, j] <- evaluateBinaryProfile(
        modelData = modelData,
        exposureColumn = "snp",
        covariateColumns = covCols,
        betaGrid = betaGrid,
        modelBackend = modelBackend,
        regularizationVariance = regularizationVariance,
        instrumentRegularization = instrumentRegularization
      )
      betaHats[j] <- unname(betaHat)
      seHats[j] <- unname(seHat)
    }, error = function(e) {
      profileMatrix[, j] <<- rep(-Inf, length(betaGrid))
      betaHats[j] <<- NA_real_
      seHats[j] <<- NA_real_
    })
  }

  colnames(profileMatrix) <- alignedInstruments$snp_id
  perSnpEstimates <- data.frame(
    snp_id = alignedInstruments$snp_id,
    effect_allele = alignedInstruments$effect_allele,
    other_allele = alignedInstruments$other_allele,
    eaf = alignedInstruments$eaf,
    beta_ZY = betaHats,
    se_ZY = seHats,
    beta_ZX = alignedInstruments$beta_ZX,
    se_ZX = alignedInstruments$se_ZX,
    pval_ZX = if ("pval_ZX" %in% names(alignedInstruments)) {
      alignedInstruments$pval_ZX
    } else {
      2 * stats::pnorm(-abs(alignedInstruments$beta_ZX / alignedInstruments$se_ZX))
    },
    pval_ZY = 2 * stats::pnorm(-abs(betaHats / seHats)),
    stringsAsFactors = FALSE
  )

  list(
    logLikProfile = alleleScoreFit$logLikProfile,
    perSnpProfiles = profileMatrix,
    betaHat = alleleScoreFit$betaHat,
    seHat = alleleScoreFit$seHat,
    perSnpBetaHats = betaHats,
    perSnpSEHats = seHats,
    perSnpEstimates = perSnpEstimates[stats::complete.cases(perSnpEstimates), , drop = FALSE],
    alignedInstruments = alignedInstruments,
    scoreDefinition = alleleScoreFit$scoreDefinition
  )
}


#' @keywords internal
alignInstrumentColumns <- function(cohortData, instrumentTable) {
  validateInstrumentTable(instrumentTable)
  expectedColumns <- makeSnpColumnName(instrumentTable$snp_id)
  present <- expectedColumns %in% names(cohortData)

  if (!any(present)) {
    stop(
      paste0(
        "No cohortData SNP columns match the instrumentTable SNP IDs. Expected columns like ",
        paste(utils::head(expectedColumns, 3), collapse = ", "),
        ". Ensure the cohort was built with the same instrument table."
      )
    )
  }

  if (any(!present)) {
    warning(
      sprintf(
        "Dropping %d instrument(s) missing from cohortData: %s",
        sum(!present),
        paste(instrumentTable$snp_id[!present], collapse = ", ")
      )
    )
  }

  list(
    instrumentTable = instrumentTable[present, , drop = FALSE],
    snpColumns = expectedColumns[present]
  )
}


#' @keywords internal
appendCovariatesToModelData <- function(modelData, cohortData, covariateData) {
  covCols <- character(0)

  # Unwrap medusaCovariateData to a plain data frame
  if (inherits(covariateData, "medusaCovariateData")) {
    covariateData <- extractCovariateDataFrame(covariateData, cohortData)
  }

  if (!is.null(covariateData) && is.data.frame(covariateData)) {
    covCols <- setdiff(names(covariateData), c("person_id", "personId", "rowId"))
    joinKey <- intersect(
      c("person_id", "personId", "rowId"),
      intersect(names(covariateData), names(cohortData))
    )

    if (length(joinKey) > 0) {
      joinKey <- joinKey[[1]]
      mergeIndex <- data.frame(joinValue = cohortData[[joinKey]])
      names(mergeIndex) <- joinKey
      mergedCov <- dplyr::left_join(
        mergeIndex,
        covariateData,
        by = joinKey
      )
    } else {
      if (nrow(covariateData) != nrow(cohortData)) {
        stop(
          "covariateData must either share a join key (person_id, personId, or rowId) with cohortData or have the same number of rows for explicit row-wise alignment."
        )
      }
      warning(
        "covariateData has no join key in common with cohortData; aligning covariates by row order."
      )
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
                                  betaGrid,
                                  modelBackend = "glm",
                                  regularizationVariance = 0.1,
                                  instrumentRegularization = FALSE) {
  checkmate::assertDataFrame(modelData, min.rows = 1)
  checkmate::assertString(exposureColumn)
  if (!is.character(covariateColumns)) {
    stop("covariateColumns must be a character vector.")
  }

  modelBackend <- match.arg(modelBackend, c("glm", "cyclops"))
  exposure <- modelData[[exposureColumn]]
  logLikProfile <- vapply(betaGrid, function(betaFixed) {
    evaluateBinaryProfilePoint(
      modelData = modelData,
      exposureColumn = exposureColumn,
      covariateColumns = covariateColumns,
      offsetVector = betaFixed * exposure,
      modelBackend = modelBackend,
      regularizationVariance = regularizationVariance,
      instrumentRegularization = instrumentRegularization
    )
  }, numeric(1))

  if (all(!is.finite(logLikProfile))) {
    stop("Profile likelihood failed at all grid points.")
  }

  logLikProfile
}


#' @keywords internal
fitBinaryOutcomeCoefficient <- function(modelData,
                                        exposureColumn,
                                        covariateColumns,
                                        modelBackend,
                                        regularizationVariance,
                                        instrumentRegularization) {
  formula <- as.formula(paste(
    "outcome ~",
    exposureColumn,
    if (length(covariateColumns) > 0) {
      paste("+", paste(covariateColumns, collapse = " + "))
    } else {
      ""
    }
  ))
  modelBackend <- match.arg(modelBackend, c("glm", "cyclops"))

  if (identical(modelBackend, "glm")) {
    fit <- glm(formula, data = modelData, family = binomial())
    betaHat <- coef(fit)[[exposureColumn]]
    seHat <- summary(fit)$coefficients[exposureColumn, "Std. Error"]
    return(list(betaHat = unname(betaHat), seHat = unname(seHat)))
  }

    fit <- fitCyclopsLogistic(
      formula = formula,
      data = modelData,
      offsetVector = NULL,
      regularizationVariance = regularizationVariance,
      excludeIndices = if (isTRUE(instrumentRegularization)) {
        1L
    } else {
      c(1L, 2L)
    }
  )
  coefficients <- stats::coef(fit)
  covariance <- stats::vcov(fit)
  betaHat <- coefficients[[exposureColumn]]
  seHat <- sqrt(diag(covariance))[[exposureColumn]]

  list(betaHat = unname(betaHat), seHat = unname(seHat))
}


#' @keywords internal
evaluateBinaryProfilePoint <- function(modelData,
                                       exposureColumn,
                                       covariateColumns,
                                       offsetVector,
                                       modelBackend,
                                       regularizationVariance,
                                       instrumentRegularization) {
  response <- modelData$outcome
  modelBackend <- match.arg(modelBackend, c("glm", "cyclops"))

  profileFormula <- if (length(covariateColumns) > 0) {
    stats::as.formula(paste("outcome ~", paste(covariateColumns, collapse = " + ")))
  } else {
    stats::as.formula("outcome ~ 1")
  }

  if (identical(modelBackend, "glm")) {
    designMatrix <- stats::model.matrix(profileFormula, data = modelData)
    fit <- tryCatch(
      suppressWarnings(
        stats::glm.fit(
          x = designMatrix,
          y = response,
          family = stats::binomial(),
          offset = offsetVector
        )
      ),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      return(-Inf)
    }
    if (!isTRUE(fit$converged) || isTRUE(fit$boundary)) {
      return(-Inf)
    }

    fittedProb <- pmin(pmax(fit$fitted.values, 1e-12), 1 - 1e-12)
    fitLogLik <- sum(stats::dbinom(response, size = 1, prob = fittedProb, log = TRUE))
    if (!is.finite(fitLogLik)) {
      return(-Inf)
    }

    return(fitLogLik)
  }

  fit <- tryCatch(
    fitCyclopsLogistic(
      formula = profileFormula,
      data = modelData,
      offsetVector = offsetVector,
      regularizationVariance = regularizationVariance,
      excludeIndices = 1L
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(-Inf)
  }

  # Cyclops predict() returns linear predictors (log-odds) for logistic
  # regression. Always apply the logistic transformation to obtain
  # probabilities, rather than guessing from the value range.
  fittedProb <- tryCatch({
    pred <- stats::predict(fit)
    1 / (1 + exp(-pred))
  }, error = function(e) NULL)
  if (is.null(fittedProb)) {
    return(-Inf)
  }
  fittedProb <- pmin(pmax(fittedProb, 1e-12), 1 - 1e-12)
  fitLogLik <- sum(stats::dbinom(response, size = 1, prob = fittedProb, log = TRUE))
  if (!is.finite(fitLogLik)) {
    return(-Inf)
  }

  fitLogLik
}


#' @keywords internal
fitCyclopsLogistic <- function(formula,
                               data,
                               offsetVector = NULL,
                               regularizationVariance,
                               excludeIndices = 1L) {
  cyclopsFormula <- formula
  cyclopsDataFrame <- data
  if (!is.null(offsetVector)) {
    cyclopsDataFrame$.__medusa_offset__ <- offsetVector
    cyclopsFormula <- stats::update.formula(
      formula,
      . ~ . + offset(.__medusa_offset__)
    )
  }

  cyclopsData <- Cyclops::createCyclopsData(
    formula = cyclopsFormula,
    modelType = "lr",
    data = cyclopsDataFrame
  )

  prior <- createCyclopsPrior(
    regularizationVariance = regularizationVariance,
    excludeIndices = excludeIndices
  )
  control <- Cyclops::createControl(noiseLevel = "silent")

  Cyclops::fitCyclopsModel(
    cyclopsData = cyclopsData,
    prior = prior,
    control = control,
    warnings = FALSE
  )
}


#' @keywords internal
createCyclopsPrior <- function(regularizationVariance,
                               excludeIndices = 1L) {
  if (!is.finite(regularizationVariance)) {
    return(Cyclops::createPrior("none"))
  }
  checkmate::assertNumber(regularizationVariance, lower = 0)
  Cyclops::createPrior(
    priorType = "normal",
    variance = regularizationVariance,
    exclude = excludeIndices
  )
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

  # Use a local quadratic fit over up to 5 points around the peak, which is

  # robust to non-uniform grid spacing (e.g., after interpolation).
  loIdx <- max(1, peakIdx - 2)
  hiIdx <- min(length(betaGrid), peakIdx + 2)
  localBeta <- betaGrid[loIdx:hiIdx]
  localLL <- logLikProfile[loIdx:hiIdx]

  if (length(localBeta) < 3) {
    return(Inf)
  }

  # Fit y = a * x^2 + b * x + c; second derivative is 2a
  quadFit <- tryCatch(
    lm(localLL ~ poly(localBeta, 2, raw = TRUE)),
    error = function(e) NULL
  )
  if (is.null(quadFit)) {
    return(Inf)
  }

  d2 <- 2 * unname(coef(quadFit)[3])

  if (is.na(d2) || d2 >= 0) {
    return(Inf)
  }

  sqrt(-1 / d2)
}


#' @keywords internal
extractCovariateDataFrame <- function(medusaCovData, cohortData) {
  # Convert a medusaCovariateData list (from buildMRCovariates) into a plain

  # data frame that appendCovariatesToModelData can merge.
  covDf <- NULL

  if (!is.null(medusaCovData$covariateData)) {
    covObj <- medusaCovData$covariateData
    # FeatureExtraction returns Andromeda-backed objects; collect to data frame
    if ("covariates" %in% names(covObj)) {
      covDf <- tryCatch(
        as.data.frame(covObj$covariates),
        error = function(e) NULL
      )
    }
  }

  # Append ancestry PCs if present
  if (!is.null(medusaCovData$ancestryPCs) &&
      is.data.frame(medusaCovData$ancestryPCs)) {
    pcDf <- medusaCovData$ancestryPCs
    if (is.null(covDf)) {
      covDf <- pcDf
    } else {
      joinCol <- intersect(names(covDf), names(pcDf))
      joinCol <- joinCol[joinCol %in% c("personId", "person_id", "rowId")]
      if (length(joinCol) > 0) {
        covDf <- merge(covDf, pcDf, by = joinCol[1], all.x = TRUE)
      }
    }
  }

  covDf
}
