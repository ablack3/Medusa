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

#' Run LD-aware cis-MR on correlated variants
#'
#' @title Correlated cis-variant Mendelian randomization
#'
#' @description
#' Performs cis-MR using correlated SNPs. The primary method is generalized IVW
#' with an optional generalized MR-Egger sensitivity analysis for invalid-IV
#' tolerance. The function expects harmonized exposure and outcome summary
#' statistics at the same variants.
#'
#' @param exposureSummary Data frame with \code{snp_id}, \code{beta_ZX}, and
#'   \code{se_ZX}.
#' @param outcomeSummary Data frame with \code{snp_id}, \code{beta_ZY}, and
#'   \code{se_ZY}.
#' @param ldMatrix Optional SNP correlation matrix. When NULL, SNPs are treated
#'   as independent.
#' @param methods Character vector containing any of \code{"gls_ivw"} and
#'   \code{"gls_egger"}.
#' @param ciLevel Numeric confidence level. Default is 0.95.
#'
#' @return A list with class \code{"medusaCisMr"}.
#'
#' @export
runCisMr <- function(exposureSummary,
                     outcomeSummary,
                     ldMatrix = NULL,
                     methods = c("gls_ivw", "gls_egger"),
                     ciLevel = 0.95) {
  exposureSummary <- normalizeCisExposureSummary(exposureSummary)
  outcomeSummary <- normalizeCisOutcomeSummary(outcomeSummary)
  methods <- unique(methods)
  checkmate::assertSubset(methods, c("gls_ivw", "gls_egger"))
  checkmate::assertNumber(ciLevel, lower = 0.5, upper = 0.999)

  merged <- merge(exposureSummary, outcomeSummary, by = "snp_id")
  if (nrow(merged) < 2) {
    stop("runCisMr() requires at least two overlapping SNPs.")
  }

  sigma <- buildOutcomeCovariance(merged$se_ZY, ldMatrix, merged$snp_id)
  weightMatrix <- solve(sigma)
  alpha <- 1 - ciLevel
  zCrit <- stats::qnorm(1 - alpha / 2)

  summaryRows <- list()

  if ("gls_ivw" %in% methods) {
    x <- matrix(merged$beta_ZX, ncol = 1)
    y <- matrix(merged$beta_ZY, ncol = 1)
    xtwx <- as.numeric(t(x) %*% weightMatrix %*% x)
    betaHat <- as.numeric((t(x) %*% weightMatrix %*% y) / xtwx)
    seHat <- sqrt(1 / xtwx)
    residual <- y - x * betaHat
    qStatistic <- as.numeric(t(residual) %*% weightMatrix %*% residual)
    pValue <- 2 * stats::pnorm(-abs(betaHat / seHat))

    summaryRows[["gls_ivw"]] <- data.frame(
      method = "GLS IVW",
      beta = betaHat,
      se = seHat,
      ci_lower = betaHat - zCrit * seHat,
      ci_upper = betaHat + zCrit * seHat,
      pval = pValue,
      intercept = NA_real_,
      intercept_se = NA_real_,
      q_statistic = qStatistic,
      stringsAsFactors = FALSE
    )
  }

  if ("gls_egger" %in% methods) {
    design <- cbind(1, merged$beta_ZX)
    xtwx <- t(design) %*% weightMatrix %*% design
    betaVec <- solve(xtwx, t(design) %*% weightMatrix %*% matrix(merged$beta_ZY, ncol = 1))
    covariance <- solve(xtwx)
    slope <- as.numeric(betaVec[2, 1])
    intercept <- as.numeric(betaVec[1, 1])
    seSlope <- sqrt(covariance[2, 2])
    seIntercept <- sqrt(covariance[1, 1])
    pValue <- 2 * stats::pnorm(-abs(slope / seSlope))

    summaryRows[["gls_egger"]] <- data.frame(
      method = "GLS Egger",
      beta = slope,
      se = seSlope,
      ci_lower = slope - zCrit * seSlope,
      ci_upper = slope + zCrit * seSlope,
      pval = pValue,
      intercept = intercept,
      intercept_se = seIntercept,
      q_statistic = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  result <- list(
    summary = do.call(rbind, summaryRows),
    mergedData = merged,
    ldMatrix = if (is.null(ldMatrix)) diag(nrow(merged)) else ldMatrix,
    fStatistics = computeApproxFStatistic(merged$beta_ZX, merged$se_ZX),
    nVariants = nrow(merged)
  )
  class(result) <- "medusaCisMr"
  result
}


#' Run single-signal coloc using Wakefield approximate Bayes factors
#'
#' @title Colocalization wrapper for target validation
#'
#' @description
#' Implements a single-signal approximate-Bayes-factor colocalization analysis
#' similar in spirit to \code{coloc.abf}. This provides a native Medusa
#' colocalization path without requiring external infrastructure.
#'
#' @param exposureSummary Data frame with \code{snp_id}, \code{beta}, and
#'   \code{se} or \code{beta_ZX} and \code{se_ZX}.
#' @param outcomeSummary Data frame with \code{snp_id}, \code{beta}, and
#'   \code{se} or \code{beta_ZY} and \code{se_ZY}.
#' @param locusWindow Optional character label for the locus window.
#' @param priors Named numeric vector with \code{p1}, \code{p2}, and
#'   \code{p12}. Defaults are conservative coloc-style priors.
#' @param priorSd Numeric prior standard deviation for causal effects.
#' @param ldReference Optional LD matrix or metadata object. Currently recorded
#'   for provenance only.
#'
#' @return A list with class \code{"medusaColoc"}.
#'
#' @export
runColocalization <- function(exposureSummary,
                              outcomeSummary,
                              locusWindow = NULL,
                              priors = c(p1 = 1e-4, p2 = 1e-4, p12 = 1e-5),
                              priorSd = 0.15,
                              ldReference = NULL) {
  exposureSummary <- normalizeColocSummary(exposureSummary, trait = "exposure")
  outcomeSummary <- normalizeColocSummary(outcomeSummary, trait = "outcome")
  checkmate::assertNumeric(priors, names = "named", len = 3, any.missing = FALSE)
  checkmate::assertNumber(priorSd, lower = 0)

  merged <- merge(exposureSummary, outcomeSummary, by = "snp_id")
  if (nrow(merged) < 2) {
    stop("runColocalization() requires at least two overlapping SNPs.")
  }

  abfExposure <- computeWakefieldAbf(merged$beta_exposure, merged$se_exposure, priorSd)
  abfOutcome <- computeWakefieldAbf(merged$beta_outcome, merged$se_outcome, priorSd)
  abfShared <- abfExposure * abfOutcome

  h0 <- 1
  h1 <- priors[["p1"]] * sum(abfExposure)
  h2 <- priors[["p2"]] * sum(abfOutcome)
  h4 <- priors[["p12"]] * sum(abfShared)
  h3 <- priors[["p1"]] * priors[["p2"]] *
    (sum(abfExposure) * sum(abfOutcome) - sum(abfShared))

  posterior <- c(H0 = h0, H1 = h1, H2 = h2, H3 = h3, H4 = h4)
  posterior <- posterior / sum(posterior)

  snpPosteriorH4 <- abfShared / sum(abfShared)
  topIdx <- which.max(snpPosteriorH4)

  result <- list(
    posterior = posterior,
    ppH4 = posterior[["H4"]],
    topSnp = merged$snp_id[topIdx],
    topSnpPosterior = snpPosteriorH4[topIdx],
    snpPosteriorH4 = data.frame(
      snp_id = merged$snp_id,
      posteriorH4 = snpPosteriorH4,
      stringsAsFactors = FALSE
    ),
    locusWindow = locusWindow,
    model = "single_signal_abf",
    ldReference = ldReference
  )
  class(result) <- "medusaColoc"
  result
}


#' Run overlap diagnostics for two-sample MR
#'
#' @title Sample-overlap and winner's-curse diagnostics
#'
#' @description
#' Flags overlap risk between exposure and outcome samples using simple,
#' transparent heuristics that can be embedded in Medusa reports. An optional
#' approximate correction factor is returned when an observed MR estimate is
#' supplied.
#'
#' @param exposureSampleSize Integer. Exposure GWAS sample size.
#' @param outcomeSampleSize Integer. Outcome sample size.
#' @param overlapCount Optional integer count of overlapping participants.
#' @param overlapProportion Optional numeric proportion of overlapping
#'   participants. Overrides \code{overlapCount} when supplied.
#' @param meanFStatistic Optional numeric mean F-statistic across instruments.
#' @param observedEstimate Optional numeric observed MR estimate for approximate
#'   de-biasing.
#'
#' @return A list with class \code{"medusaOverlapDiagnostics"}.
#'
#' @export
runOverlapDiagnostics <- function(exposureSampleSize,
                                  outcomeSampleSize,
                                  overlapCount = NULL,
                                  overlapProportion = NULL,
                                  meanFStatistic = NULL,
                                  observedEstimate = NULL) {
  checkmate::assertCount(exposureSampleSize)
  checkmate::assertCount(outcomeSampleSize)
  if (!is.null(overlapCount)) {
    checkmate::assertCount(overlapCount, positive = FALSE)
  }
  if (!is.null(overlapProportion)) {
    checkmate::assertNumber(overlapProportion, lower = 0, upper = 1)
  }
  if (!is.null(meanFStatistic)) {
    checkmate::assertNumber(meanFStatistic, lower = 0)
  }

  if (is.null(overlapProportion)) {
    overlapProportion <- if (is.null(overlapCount)) {
      0
    } else {
      overlapCount / min(exposureSampleSize, outcomeSampleSize)
    }
  }

  biasIndex <- if (!is.null(meanFStatistic) && meanFStatistic > 0) {
    overlapProportion / meanFStatistic
  } else {
    NA_real_
  }

  riskLevel <- if (overlapProportion == 0) {
    "low"
  } else if (is.na(biasIndex) || biasIndex >= 0.05 || overlapProportion >= 0.3) {
    "high"
  } else if (biasIndex >= 0.01 || overlapProportion >= 0.1) {
    "moderate"
  } else {
    "low"
  }

  correctionFactor <- if (!is.na(biasIndex)) max(1 - biasIndex, 0.5) else NA_real_
  correctedEstimate <- if (!is.null(observedEstimate) && is.finite(correctionFactor)) {
    observedEstimate / correctionFactor
  } else {
    NA_real_
  }

  structure(
    list(
      overlapProportion = overlapProportion,
      meanFStatistic = meanFStatistic,
      biasIndex = biasIndex,
      riskLevel = riskLevel,
      correctionFactor = correctionFactor,
      correctedEstimate = correctedEstimate
    ),
    class = "medusaOverlapDiagnostics"
  )
}


#' Run ancestry-composition diagnostics
#'
#' @title Exposure/outcome/LD population-matching diagnostics
#'
#' @description
#' Compares ancestry composition across exposure summary statistics, outcome
#' cohorts, and LD references. The check fails closed when the dominant
#' population mismatch exceeds the specified threshold.
#'
#' @param exposurePopulation Character or named numeric vector.
#' @param outcomePopulation Character or named numeric vector.
#' @param ldReferencePopulation Optional character or named numeric vector.
#' @param mismatchThreshold Numeric threshold on total-variation distance.
#'
#' @return A list with class \code{"medusaAncestryDiagnostics"}.
#'
#' @export
runAncestryDiagnostics <- function(exposurePopulation,
                                   outcomePopulation,
                                   ldReferencePopulation = NULL,
                                   mismatchThreshold = 0.25) {
  checkmate::assertNumber(mismatchThreshold, lower = 0, upper = 1)

  exposureComp <- normalizePopulationComposition(exposurePopulation)
  outcomeComp <- normalizePopulationComposition(outcomePopulation)
  ldComp <- if (is.null(ldReferencePopulation)) {
    NULL
  } else {
    normalizePopulationComposition(ldReferencePopulation)
  }

  exposureOutcomeDistance <- populationDistance(exposureComp, outcomeComp)
  outcomeLdDistance <- if (is.null(ldComp)) NA_real_ else populationDistance(outcomeComp, ldComp)
  exposureLdDistance <- if (is.null(ldComp)) NA_real_ else populationDistance(exposureComp, ldComp)

  failClosed <- isTRUE(exposureOutcomeDistance > mismatchThreshold) ||
    isTRUE(outcomeLdDistance > mismatchThreshold) ||
    isTRUE(exposureLdDistance > mismatchThreshold)

  structure(
    list(
      exposureOutcomeDistance = exposureOutcomeDistance,
      outcomeLdDistance = outcomeLdDistance,
      exposureLdDistance = exposureLdDistance,
      mismatchThreshold = mismatchThreshold,
      failClosed = failClosed,
      riskLevel = if (failClosed) "high" else if (exposureOutcomeDistance > mismatchThreshold / 2) "moderate" else "low"
    ),
    class = "medusaAncestryDiagnostics"
  )
}


#' Classify a target decision from Medusa evidence
#'
#' @title Convert study evidence into advance/restrict/stop recommendations
#'
#' @param studySpec A \code{medusaStudySpec}.
#' @param cisMrResults Optional \code{medusaCisMr} result.
#' @param colocResults Optional \code{medusaColoc} result.
#' @param diagnosticResults Optional diagnostics list.
#' @param ancestryDiagnostics Optional ancestry diagnostics result.
#'
#' @return Named list with recommendation and reasons.
#'
#' @export
classifyTargetDecision <- function(studySpec,
                                   cisMrResults = NULL,
                                   colocResults = NULL,
                                   diagnosticResults = NULL,
                                   ancestryDiagnostics = NULL) {
  validateMedusaStudySpec(studySpec)
  thresholds <- studySpec$decisionThresholds
  reasons <- character(0)
  recommendation <- "restrict"

  meanF <- if (!is.null(cisMrResults) && !is.null(cisMrResults$fStatistics)) {
    mean(cisMrResults$fStatistics, na.rm = TRUE)
  } else {
    NA_real_
  }

  if (!is.na(meanF) && meanF < thresholds$minFStatistic) {
    reasons <- c(reasons, sprintf("Mean F-statistic %.2f is below %.2f.", meanF, thresholds$minFStatistic))
    recommendation <- "stop"
  }

  ppH4 <- if (!is.null(colocResults)) colocResults$ppH4 else NA_real_
  if (!is.na(ppH4) && ppH4 < thresholds$minColocPp4) {
    reasons <- c(reasons, sprintf("Colocalization PP4 %.2f is below %.2f.", ppH4, thresholds$minColocPp4))
    recommendation <- "stop"
  }

  concordantDatasets <- cisMrResults$concordantDatasets %||% 1L
  if (concordantDatasets < thresholds$minConcordantDatasets) {
    reasons <- c(
      reasons,
      sprintf(
        "Only %d concordant dataset(s) available; %d required.",
        concordantDatasets,
        thresholds$minConcordantDatasets
      )
    )
    if (!identical(recommendation, "stop")) {
      recommendation <- "restrict"
    }
  }

  if (isTRUE(thresholds$failOnAncestryMismatch) &&
      !is.null(ancestryDiagnostics) &&
      isTRUE(ancestryDiagnostics$failClosed)) {
    reasons <- c(reasons, "Ancestry mismatch diagnostics failed closed.")
    recommendation <- "stop"
  }

  negativeControlFailure <- !is.null(diagnosticResults) &&
    !is.null(diagnosticResults$diagnosticFlags) &&
    isTRUE(diagnosticResults$diagnosticFlags[["negativeControlFailure"]])
  if (isTRUE(thresholds$failOnNegativeControls) && negativeControlFailure) {
    reasons <- c(reasons, "Negative control analysis suggests residual bias.")
    recommendation <- "stop"
  }

  liabilityFlag <- !is.null(diagnosticResults) &&
    !is.null(diagnosticResults$diagnosticFlags) &&
    isTRUE(diagnosticResults$diagnosticFlags[["phewasSignificant"]])
  if (liabilityFlag && identical(recommendation, "advance")) {
    recommendation <- "restrict"
  } else if (liabilityFlag && !identical(recommendation, "stop")) {
    reasons <- c(reasons, "Phenome-wide liabilities require indication restriction.")
  }

  if (length(reasons) == 0) {
    recommendation <- "advance"
    reasons <- "All pre-specified evidence gates passed."
  }

  list(
    recommendation = recommendation,
    reasons = reasons,
    meanFStatistic = meanF,
    ppH4 = ppH4,
    concordantDatasets = concordantDatasets
  )
}


#' Create the built-in oncology outcome recipe library
#'
#' @title Binary endpoint recipes for the Medusa oncology program
#'
#' @return Named list of endpoint definitions.
#'
#' @export
getOncologyOutcomeLibrary <- function() {
  list(
    incidentCancer = list(
      label = "Incident cancer",
      requiresTreatmentData = FALSE,
      washoutDays = 365L,
      recommendedDomains = c("condition_occurrence", "observation_period")
    ),
    metastaticAtPresentation = list(
      label = "Metastatic-at-presentation proxy",
      requiresTreatmentData = FALSE,
      washoutDays = 365L,
      recommendedDomains = c("condition_occurrence", "measurement", "observation")
    ),
    systemicTreatmentWithin12Months = list(
      label = "Systemic treatment within 12 months",
      requiresTreatmentData = TRUE,
      lookforwardDays = 365L,
      recommendedDomains = c("drug_exposure", "procedure_occurrence")
    ),
    secondLineWithin12Months = list(
      label = "Second-line treatment within 12 months",
      requiresTreatmentData = TRUE,
      lookforwardDays = 365L,
      recommendedDomains = c("drug_exposure", "procedure_occurrence")
    )
  )
}


#' Create an oncology outcome definition
#'
#' @param cancer Character cancer label.
#' @param endpointType Character endpoint recipe name from
#'   \code{\link{getOncologyOutcomeLibrary}}.
#' @param cohortId Optional integer cohort identifier.
#'
#' @return Named list combining the cancer label and endpoint recipe.
#'
#' @export
createOncologyOutcomeDefinition <- function(cancer,
                                            endpointType = "incidentCancer",
                                            cohortId = NULL) {
  checkmate::assertString(cancer)
  libraryDef <- getOncologyOutcomeLibrary()
  checkmate::assertChoice(endpointType, names(libraryDef))

  definition <- c(
    list(
      cancer = cancer,
      endpointType = endpointType,
      cohortId = cohortId
    ),
    libraryDef[[endpointType]]
  )
  definition
}


#' Generate a target-validation report
#'
#' @title HTML report for oncology target triage
#'
#' @param studySpec Study specification object.
#' @param cisMrResults Optional \code{runCisMr()} output.
#' @param colocResults Optional \code{runColocalization()} output.
#' @param diagnosticResults Optional diagnostics result.
#' @param ancestryDiagnostics Optional ancestry diagnostics result.
#' @param overlapDiagnostics Optional overlap diagnostics result.
#' @param trialContext Optional character vector of trial notes.
#' @param outputPath Character output path for the HTML report.
#'
#' @return Output path, invisibly.
#'
#' @export
generateTargetValidationReport <- function(studySpec,
                                           cisMrResults = NULL,
                                           colocResults = NULL,
                                           diagnosticResults = NULL,
                                           ancestryDiagnostics = NULL,
                                           overlapDiagnostics = NULL,
                                           trialContext = NULL,
                                           outputPath = "./Medusa_target_validation.html") {
  validateMedusaStudySpec(studySpec)
  checkmate::assertString(outputPath)
  decision <- classifyTargetDecision(
    studySpec = studySpec,
    cisMrResults = cisMrResults,
    colocResults = colocResults,
    diagnosticResults = diagnosticResults,
    ancestryDiagnostics = ancestryDiagnostics
  )

  renderMedusaHtmlReport(
    templateName = "TargetValidationReport.Rmd",
    outputPath = outputPath,
    params = list(
      studySpec = studySpec,
      cisMrResults = cisMrResults,
      colocResults = colocResults,
      diagnosticResults = diagnosticResults,
      ancestryDiagnostics = ancestryDiagnostics,
      overlapDiagnostics = overlapDiagnostics,
      trialContext = trialContext,
      decision = decision
    ),
    contextLabel = "Target validation report"
  )
}


#' @keywords internal
normalizeCisExposureSummary <- function(exposureSummary) {
  checkmate::assertDataFrame(exposureSummary, min.rows = 1)
  snpCol <- detectExistingColumn(names(exposureSummary), c("snp_id", "snpId"))
  betaCol <- detectExistingColumn(names(exposureSummary), c("beta_ZX", "beta"))
  seCol <- detectExistingColumn(names(exposureSummary), c("se_ZX", "se"))
  if (any(vapply(list(snpCol, betaCol, seCol), is.null, logical(1)))) {
    stop("Exposure summary must contain SNP, beta_ZX/beta, and se_ZX/se columns.")
  }
  data.frame(
    snp_id = exposureSummary[[snpCol]],
    beta_ZX = exposureSummary[[betaCol]],
    se_ZX = exposureSummary[[seCol]],
    stringsAsFactors = FALSE
  )
}


#' @keywords internal
normalizeCisOutcomeSummary <- function(outcomeSummary) {
  checkmate::assertDataFrame(outcomeSummary, min.rows = 1)
  snpCol <- detectExistingColumn(names(outcomeSummary), c("snp_id", "snpId"))
  betaCol <- detectExistingColumn(names(outcomeSummary), c("beta_ZY", "beta"))
  seCol <- detectExistingColumn(names(outcomeSummary), c("se_ZY", "se"))
  if (any(vapply(list(snpCol, betaCol, seCol), is.null, logical(1)))) {
    stop("Outcome summary must contain SNP, beta_ZY/beta, and se_ZY/se columns.")
  }
  data.frame(
    snp_id = outcomeSummary[[snpCol]],
    beta_ZY = outcomeSummary[[betaCol]],
    se_ZY = outcomeSummary[[seCol]],
    stringsAsFactors = FALSE
  )
}


#' @keywords internal
normalizeColocSummary <- function(summaryData, trait) {
  checkmate::assertDataFrame(summaryData, min.rows = 1)
  snpCol <- detectExistingColumn(names(summaryData), c("snp_id", "snpId"))
  betaCol <- detectExistingColumn(names(summaryData), c("beta", if (trait == "exposure") "beta_ZX" else "beta_ZY"))
  seCol <- detectExistingColumn(names(summaryData), c("se", if (trait == "exposure") "se_ZX" else "se_ZY"))
  if (any(vapply(list(snpCol, betaCol, seCol), is.null, logical(1)))) {
    stop(sprintf("Colocalization %s summary is missing SNP/beta/se columns.", trait))
  }
  result <- data.frame(
    snp_id = summaryData[[snpCol]],
    beta = summaryData[[betaCol]],
    se = summaryData[[seCol]],
    stringsAsFactors = FALSE
  )
  names(result)[2:3] <- paste0(c("beta_", "se_"), trait)
  result
}


#' @keywords internal
detectExistingColumn <- function(columnNames, candidates) {
  match <- intersect(candidates, columnNames)
  if (length(match) == 0) {
    return(NULL)
  }
  match[[1]]
}


#' @keywords internal
buildOutcomeCovariance <- function(seZY, ldMatrix, snpIds) {
  if (is.null(ldMatrix)) {
    return(diag(seZY^2))
  }

  if (is.data.frame(ldMatrix)) {
    ldMatrix <- as.matrix(ldMatrix)
  }
  checkmate::assertMatrix(ldMatrix, mode = "numeric", nrows = length(seZY), ncols = length(seZY))
  rownames(ldMatrix) <- rownames(ldMatrix) %||% snpIds
  colnames(ldMatrix) <- colnames(ldMatrix) %||% snpIds
  diag(seZY) %*% ldMatrix %*% diag(seZY)
}


#' @keywords internal
computeWakefieldAbf <- function(beta, se, priorSd) {
  variance <- se^2
  priorVariance <- priorSd^2
  z <- beta / se
  shrinkage <- priorVariance / (variance + priorVariance)
  exp(0.5 * (log(1 - shrinkage) + shrinkage * z^2))
}


#' @keywords internal
normalizePopulationComposition <- function(population) {
  if (is.character(population) && length(population) == 1) {
    composition <- 1
    names(composition) <- population
    return(composition)
  }
  checkmate::assertNumeric(population, lower = 0, upper = 1, names = "named")
  total <- sum(population)
  if (total <= 0) {
    stop("Population composition must sum to a positive value.")
  }
  population / total
}


#' @keywords internal
populationDistance <- function(x, y) {
  populations <- union(names(x), names(y))
  xFull <- stats::setNames(rep(0, length(populations)), populations)
  yFull <- xFull
  xFull[names(x)] <- x
  yFull[names(y)] <- y
  sum(abs(xFull - yFull)) / 2
}
