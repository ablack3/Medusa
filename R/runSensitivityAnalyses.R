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

#' Run MR Sensitivity Analyses
#'
#' @title Standard Mendelian Randomization sensitivity analyses
#'
#' @description
#' Implements multiple MR estimation methods that are robust to different
#' patterns of instrument invalidity. Requires per-SNP beta_ZY estimates
#' (available when \code{analysisType = "perSNP"} in \code{\link{fitOutcomeModel}}).
#'
#' Methods include: Inverse Variance Weighted (IVW), MR-Egger regression,
#' Weighted Median, Steiger directional filtering, and Leave-One-Out analysis.
#' Concordance across methods strengthens causal evidence.
#'
#' @param perSnpEstimates Data frame with per-SNP summary statistics. Required
#'   columns: \code{snp_id}, \code{beta_ZY}, \code{se_ZY}, \code{beta_ZX},
#'   \code{se_ZX}, \code{effect_allele}, \code{other_allele}, and \code{eaf}.
#'   These allele fields make the object harmonisation-ready for delegation to
#'   \pkg{TwoSampleMR}. Optional columns \code{pval_ZX} and \code{pval_ZY} are
#'   used when present and otherwise derived from the reported estimates.
#' @param methods Character vector of methods to run. Default is
#'   \code{c("IVW", "MREgger", "WeightedMedian", "Steiger", "LeaveOneOut")}.
#' @param outcomeSampleSize Optional integer. Total sample size for outcome
#'   GWAS/analysis. Needed for Steiger filtering.
#' @param exposureSampleSize Optional integer. Total sample size for exposure
#'   GWAS. Needed for Steiger filtering.
#' @param outcomeType Character. Type of outcome summary statistics represented
#'   in \code{perSnpEstimates}. Medusa's internal Steiger implementation is
#'   available only for \code{"continuous"} outcomes, while
#'   \code{engine = "TwoSampleMR"} delegates to \pkg{TwoSampleMR}'s Steiger
#'   implementation on the harmonised data. Default is \code{"binary"}.
#' @param engine Character. \code{"auto"} prefers \pkg{TwoSampleMR} when
#'   installed and otherwise falls back to Medusa's internal implementations.
#'   Use \code{"TwoSampleMR"} to require delegation, or \code{"internal"} to
#'   force Medusa's built-in implementations. Default is \code{"auto"}.
#'
#' @return A named list with class "medusaSensitivity" containing:
#'   \describe{
#'     \item{ivw}{Data frame with method, beta_MR, se_MR, ci_lower, ci_upper, pval.}
#'     \item{mrEgger}{Data frame with beta_MR, se_MR, pval, plus intercept,
#'       intercept_se, intercept_pval.}
#'     \item{weightedMedian}{Data frame with beta_MR, se_MR, ci_lower, ci_upper, pval.}
#'     \item{steiger}{Data frame with IVW results after Steiger filtering, plus
#'       n_removed (number of SNPs removed).}
#'     \item{leaveOneOut}{Data frame with snp_removed, beta_MR, se_MR, pval
#'       for each SNP dropped.}
#'     \item{summary}{Data frame comparing all methods side by side.}
#'   }
#'
#' @details
#' \strong{IVW}: Weighted regression of beta_ZY on beta_ZX through the origin,
#' with weights 1/se_ZY^2. This is the primary MR estimate assuming all
#' instruments are valid.
#'
#' \strong{MR-Egger}: Weighted regression of beta_ZY on beta_ZX with an
#' intercept, after orienting all SNPs so that beta_ZX is positive. A non-zero
#' intercept indicates directional pleiotropy. The slope is a consistent causal
#' estimate even under directional pleiotropy, provided the InSIDE assumption
#' holds.
#'
#' \strong{Weighted Median}: Weighted median of per-SNP Wald ratio estimates
#' using first-order inverse-variance weights for the ratio estimates
#' (\eqn{w_j = \beta_{Xj}^2 / \mathrm{SE}(\beta_{Yj})^2}). Produces a consistent
#' estimate if fewer than 50% of the total weight comes from invalid instruments.
#'
#' \strong{Steiger}: Tests whether each SNP explains more variance in the
#' exposure than the outcome (the expected causal direction). With
#' \code{engine = "internal"}, Medusa applies the continuous-trait
#' correlation-based approximation described by Hemani et al. (2017). With
#' \code{engine = "TwoSampleMR"}, Medusa delegates to
#' \code{TwoSampleMR::steiger_filtering()} on the harmonised summary data.
#'
#' \strong{Leave-One-Out}: Drops each SNP in turn and recomputes the IVW
#' estimate. Identifies influential outlier instruments.
#'
#' @references
#' Bowden, J., Davey Smith, G., & Burgess, S. (2015). Mendelian randomization
#' with invalid instruments: effect estimation and bias detection through Egger
#' regression. \emph{International Journal of Epidemiology}, 44(2), 512-525.
#' doi:10.1093/ije/dyv080.
#' Abstract: https://pubmed.ncbi.nlm.nih.gov/26050253/
#' Full text (if available to you): https://academic.oup.com/ije/article-pdf/44/2/512/18631079/dyv080.pdf
#'
#' Bowden, J., et al. (2016). Consistent estimation in Mendelian randomization
#' with some invalid instruments using a weighted median estimator.
#' \emph{Genetic Epidemiology}, 40(4), 304-314.
#' doi:10.1002/gepi.21965.
#' Open access: https://pmc.ncbi.nlm.nih.gov/articles/PMC4849733/
#'
#' Hemani, G., Tilling, K., & Davey Smith, G. (2017). Orienting the causal
#' relationship between imprecisely measured traits using GWAS summary data.
#' \emph{PLoS Genetics}, 13(11), e1007081.
#' doi:10.1371/journal.pgen.1007081.
#' Open access: https://pmc.ncbi.nlm.nih.gov/articles/PMC5711033/
#'
#' @examples
#' # Simulate per-SNP estimates
#' set.seed(42)
#' nSnps <- 10
#' betaZX <- rnorm(nSnps, 0.3, 0.05)
#' betaZY <- 0.5 * betaZX + rnorm(nSnps, 0, 0.02)
#' perSnp <- data.frame(
#'   snp_id = paste0("rs", 1:nSnps),
#'   effect_allele = c("A", "C", "G", "T", "A", "C", "G", "T", "A", "C"),
#'   other_allele = c("C", "G", "T", "A", "C", "G", "T", "A", "C", "G"),
#'   eaf = seq(0.1, 0.55, length.out = nSnps),
#'   beta_ZY = betaZY,
#'   se_ZY = rep(0.02, nSnps),
#'   beta_ZX = betaZX,
#'   se_ZX = rep(0.05, nSnps)
#' )
#' results <- runSensitivityAnalyses(
#'   perSnp,
#'   outcomeSampleSize = 10000,
#'   exposureSampleSize = 10000,
#'   outcomeType = "continuous"
#' )
#' results$summary
#'
#' @seealso \code{\link{computeMREstimate}}, \code{\link{fitOutcomeModel}},
#'   \code{\link{generateMRReport}}
#'
#' @export
runSensitivityAnalyses <- function(perSnpEstimates,
                                   methods = c("IVW", "MREgger",
                                               "WeightedMedian", "Steiger",
                                               "LeaveOneOut"),
                                   outcomeSampleSize = NULL,
                                   exposureSampleSize = NULL,
                                   outcomeType = "binary",
                                   engine = "auto") {
  # Input validation
  validatePerSnpSummaryData(perSnpEstimates)
  checkmate::assertSubset(methods,
                          c("IVW", "MREgger", "WeightedMedian",
                            "Steiger", "LeaveOneOut"))
  checkmate::assertChoice(outcomeType, c("binary", "continuous"))
  checkmate::assertChoice(engine, c("auto", "TwoSampleMR", "internal"))

  nSnps <- nrow(perSnpEstimates)
  engine <- resolveSensitivityEngine(engine)
  message(sprintf("Running sensitivity analyses with %d SNPs...", nSnps))
  message(sprintf("  Engine: %s", engine))

  results <- list()
  harmonisedData <- NULL
  if (identical(engine, "TwoSampleMR")) {
    harmonisedData <- harmonisePerSnpEstimates(
      perSnpEstimates = perSnpEstimates,
      outcomeSampleSize = outcomeSampleSize,
      exposureSampleSize = exposureSampleSize,
      outcomeType = outcomeType
    )
  }

  # --- IVW ---
  if ("IVW" %in% methods) {
    message("  IVW...")
    results$ivw <- if (identical(engine, "TwoSampleMR")) {
      computeIVWTwoSampleMR(harmonisedData)
    } else {
      computeIVW(perSnpEstimates)
    }
  }

  # --- MR-Egger ---
  if ("MREgger" %in% methods) {
    if (nSnps < 3) {
      message("  MR-Egger: skipped (requires >= 3 SNPs).")
      results$mrEgger <- NULL
    } else {
      message("  MR-Egger...")
      results$mrEgger <- if (identical(engine, "TwoSampleMR")) {
        computeMREggerTwoSampleMR(harmonisedData)
      } else {
        computeMREgger(perSnpEstimates)
      }
    }
  }

  # --- Weighted Median ---
  if ("WeightedMedian" %in% methods) {
    if (nSnps < 3) {
      message("  Weighted Median: skipped (requires >= 3 SNPs).")
      results$weightedMedian <- NULL
    } else {
      message("  Weighted Median...")
      results$weightedMedian <- if (identical(engine, "TwoSampleMR")) {
        computeWeightedMedianTwoSampleMR(harmonisedData)
      } else {
        computeWeightedMedian(perSnpEstimates)
      }
    }
  }

  # --- Steiger ---
  if ("Steiger" %in% methods) {
    message("  Steiger filtering...")
    results$steiger <- if (identical(engine, "TwoSampleMR")) {
      computeSteigerTwoSampleMR(
        harmonisedData = harmonisedData,
        outcomeType = outcomeType
      )
    } else {
      computeSteiger(perSnpEstimates, outcomeSampleSize,
                     exposureSampleSize, outcomeType)
    }
  }

  # --- Leave-One-Out ---
  if ("LeaveOneOut" %in% methods) {
    if (nSnps < 3) {
      message("  Leave-One-Out: skipped (requires >= 3 SNPs).")
      results$leaveOneOut <- NULL
    } else {
      message("  Leave-One-Out...")
      results$leaveOneOut <- if (identical(engine, "TwoSampleMR")) {
        computeLeaveOneOutTwoSampleMR(harmonisedData)
      } else {
        computeLeaveOneOut(perSnpEstimates)
      }
    }
  }

  # --- Summary table ---
  results$summary <- buildSensitivitySummary(results)

  class(results) <- "medusaSensitivity"

  message("Sensitivity analyses complete.")
  results
}


#' @keywords internal
validatePerSnpSummaryData <- function(perSnpEstimates) {
  checkmate::assertDataFrame(perSnpEstimates, min.rows = 1)
  requiredCols <- c(
    "snp_id", "beta_ZY", "se_ZY", "beta_ZX", "se_ZX",
    "effect_allele", "other_allele", "eaf"
  )
  checkmate::assertSubset(requiredCols, names(perSnpEstimates))
  checkmate::assertCharacter(perSnpEstimates$snp_id, any.missing = FALSE)
  checkmate::assertCharacter(perSnpEstimates$effect_allele, any.missing = FALSE)
  checkmate::assertCharacter(perSnpEstimates$other_allele, any.missing = FALSE)
  checkmate::assertNumeric(perSnpEstimates$beta_ZY, any.missing = FALSE)
  checkmate::assertNumeric(perSnpEstimates$se_ZY, lower = 0, any.missing = FALSE)
  checkmate::assertNumeric(perSnpEstimates$beta_ZX, any.missing = FALSE)
  checkmate::assertNumeric(perSnpEstimates$se_ZX, lower = 0, any.missing = FALSE)
  checkmate::assertNumeric(perSnpEstimates$eaf, lower = 0, upper = 1,
                           any.missing = FALSE)
  invisible(TRUE)
}


#' @keywords internal
resolveSensitivityEngine <- function(engine) {
  engine <- match.arg(engine, c("auto", "TwoSampleMR", "internal"))
  if (identical(engine, "auto")) {
    if (requireNamespace("TwoSampleMR", quietly = TRUE)) {
      return("TwoSampleMR")
    }
    warning(
      "Package 'TwoSampleMR' is not installed; falling back to Medusa's internal sensitivity estimators."
    )
    return("internal")
  }
  if (identical(engine, "TwoSampleMR") &&
      !requireNamespace("TwoSampleMR", quietly = TRUE)) {
    stop(
      "Package 'TwoSampleMR' is required when engine = 'TwoSampleMR'. ",
      "Install it with: remotes::install_github('mrcieu/TwoSampleMR')",
      call. = FALSE
    )
  }
  engine
}


#' @keywords internal
harmonisePerSnpEstimates <- function(perSnpEstimates,
                                     outcomeSampleSize = NULL,
                                     exposureSampleSize = NULL,
                                     outcomeType = "binary") {
  exposureData <- buildTwoSampleMRData(
    perSnpEstimates = perSnpEstimates,
    type = "exposure",
    sampleSize = exposureSampleSize,
    outcomeType = outcomeType
  )
  outcomeData <- buildTwoSampleMRData(
    perSnpEstimates = perSnpEstimates,
    type = "outcome",
    sampleSize = outcomeSampleSize,
    outcomeType = outcomeType
  )

  exposureFormatted <- TwoSampleMR::format_data(
    exposureData,
    type = "exposure",
    phenotype_col = "phenotype",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    eaf_col = "eaf",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval",
    samplesize_col = "samplesize",
    ncase_col = "ncase",
    ncontrol_col = "ncontrol",
    id_col = "id"
  )
  outcomeFormatted <- TwoSampleMR::format_data(
    outcomeData,
    type = "outcome",
    phenotype_col = "phenotype",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    eaf_col = "eaf",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval",
    samplesize_col = "samplesize",
    ncase_col = "ncase",
    ncontrol_col = "ncontrol",
    id_col = "id"
  )

  harmonised <- TwoSampleMR::harmonise_data(
    exposure_dat = exposureFormatted,
    outcome_dat = outcomeFormatted,
    action = 2
  )
  keepMask <- !is.na(harmonised$mr_keep) & harmonised$mr_keep
  harmonised <- harmonised[keepMask, , drop = FALSE]
  harmonised <- addTwoSampleMRUnits(
    harmonised = harmonised,
    perSnpEstimates = perSnpEstimates,
    outcomeType = outcomeType
  )

  if (nrow(harmonised) == 0) {
    stop(
      "No SNPs remained after harmonisation. Check allele coding and strand ambiguity."
    )
  }

  harmonised
}


#' @keywords internal
buildTwoSampleMRData <- function(perSnpEstimates,
                                 type = c("exposure", "outcome"),
                                 sampleSize = NULL,
                                 outcomeType = "binary") {
  type <- match.arg(type)
  if (identical(type, "exposure")) {
    beta <- perSnpEstimates$beta_ZX
    se <- perSnpEstimates$se_ZX
    pval <- if ("pval_ZX" %in% names(perSnpEstimates)) {
      perSnpEstimates$pval_ZX
    } else {
      2 * stats::pnorm(-abs(beta / se))
    }
    phenotype <- if ("exposure_label" %in% names(perSnpEstimates)) {
      perSnpEstimates$exposure_label
    } else {
      rep("exposure", nrow(perSnpEstimates))
    }
    identifier <- if ("id_exposure" %in% names(perSnpEstimates)) {
      perSnpEstimates$id_exposure
    } else {
      rep("medusa_exposure", nrow(perSnpEstimates))
    }
    ncase <- rep(NA_real_, nrow(perSnpEstimates))
    ncontrol <- rep(NA_real_, nrow(perSnpEstimates))
  } else {
    beta <- perSnpEstimates$beta_ZY
    se <- perSnpEstimates$se_ZY
    pval <- if ("pval_ZY" %in% names(perSnpEstimates)) {
      perSnpEstimates$pval_ZY
    } else {
      2 * stats::pnorm(-abs(beta / se))
    }
    phenotype <- if ("outcome_label" %in% names(perSnpEstimates)) {
      perSnpEstimates$outcome_label
    } else {
      rep("outcome", nrow(perSnpEstimates))
    }
    identifier <- if ("id_outcome" %in% names(perSnpEstimates)) {
      perSnpEstimates$id_outcome
    } else {
      rep("medusa_outcome", nrow(perSnpEstimates))
    }
    ncase <- if ("ncase_outcome" %in% names(perSnpEstimates)) {
      perSnpEstimates$ncase_outcome
    } else {
      rep(NA_real_, nrow(perSnpEstimates))
    }
    ncontrol <- if ("ncontrol_outcome" %in% names(perSnpEstimates)) {
      perSnpEstimates$ncontrol_outcome
    } else {
      rep(NA_real_, nrow(perSnpEstimates))
    }
  }

  data.frame(
    phenotype = as.character(phenotype),
    SNP = perSnpEstimates$snp_id,
    beta = beta,
    se = se,
    effect_allele = perSnpEstimates$effect_allele,
    other_allele = perSnpEstimates$other_allele,
    eaf = perSnpEstimates$eaf,
    pval = pval,
    samplesize = rep_len(
      if (is.null(sampleSize)) NA_real_ else as.numeric(sampleSize),
      nrow(perSnpEstimates)
    ),
    ncase = ncase,
    ncontrol = ncontrol,
    id = as.character(identifier),
    stringsAsFactors = FALSE
  )
}


#' @keywords internal
addTwoSampleMRUnits <- function(harmonised,
                                perSnpEstimates,
                                outcomeType = "binary") {
  exposureUnits <- if ("units_exposure" %in% names(perSnpEstimates)) {
    perSnpEstimates$units_exposure
  } else {
    rep("unit increase", nrow(perSnpEstimates))
  }
  outcomeUnits <- if ("units_outcome" %in% names(perSnpEstimates)) {
    perSnpEstimates$units_outcome
  } else if (identical(outcomeType, "binary")) {
    rep("log odds", nrow(perSnpEstimates))
  } else {
    rep("unit increase", nrow(perSnpEstimates))
  }

  matchIdx <- match(harmonised$SNP, perSnpEstimates$snp_id)
  harmonised$units.exposure <- exposureUnits[matchIdx]
  harmonised$units.outcome <- outcomeUnits[matchIdx]
  harmonised
}


#' @keywords internal
computeIVWTwoSampleMR <- function(harmonisedData) {
  methodId <- if (nrow(harmonisedData) == 1) "mr_wald_ratio" else "mr_ivw"
  mrResult <- TwoSampleMR::mr(harmonisedData, method_list = methodId)
  row <- mrResult[1, , drop = FALSE]

  heterogeneity <- NULL
  if (nrow(harmonisedData) > 1) {
    heterogeneity <- tryCatch(
      TwoSampleMR::mr_heterogeneity(harmonisedData, method_list = "mr_ivw"),
      error = function(e) NULL
    )
  }

  data.frame(
    method = "IVW",
    beta_MR = row$b,
    se_MR = row$se,
    ci_lower = row$b - 1.96 * row$se,
    ci_upper = row$b + 1.96 * row$se,
    pval = row$pval,
    cochran_Q = if (!is.null(heterogeneity) && nrow(heterogeneity) > 0) {
      heterogeneity$Q[1]
    } else {
      NA_real_
    },
    cochran_Q_pval = if (!is.null(heterogeneity) && nrow(heterogeneity) > 0) {
      heterogeneity$Q_pval[1]
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
}


#' @keywords internal
computeMREggerTwoSampleMR <- function(harmonisedData) {
  mrResult <- TwoSampleMR::mr(harmonisedData, method_list = "mr_egger_regression")
  pleiotropy <- TwoSampleMR::mr_pleiotropy_test(harmonisedData)
  row <- mrResult[1, , drop = FALSE]
  interceptRow <- pleiotropy[1, , drop = FALSE]

  data.frame(
    method = "MR-Egger",
    beta_MR = row$b,
    se_MR = row$se,
    ci_lower = row$b - 1.96 * row$se,
    ci_upper = row$b + 1.96 * row$se,
    pval = row$pval,
    intercept = interceptRow$egger_intercept,
    intercept_se = interceptRow$se,
    intercept_pval = interceptRow$pval,
    stringsAsFactors = FALSE
  )
}


#' @keywords internal
computeWeightedMedianTwoSampleMR <- function(harmonisedData) {
  mrResult <- TwoSampleMR::mr(harmonisedData, method_list = "mr_weighted_median")
  row <- mrResult[1, , drop = FALSE]

  data.frame(
    method = "Weighted Median",
    beta_MR = row$b,
    se_MR = row$se,
    ci_lower = row$b - 1.96 * row$se,
    ci_upper = row$b + 1.96 * row$se,
    pval = row$pval,
    stringsAsFactors = FALSE
  )
}


#' @keywords internal
computeSteigerTwoSampleMR <- function(harmonisedData,
                                      outcomeType = "binary") {
  steigerData <- tryCatch(
    TwoSampleMR::steiger_filtering(harmonisedData),
    error = function(e) {
      warning(
        sprintf(
          "Steiger filtering could not be evaluated by TwoSampleMR: %s",
          conditionMessage(e)
        )
      )
      NULL
    }
  )

  if (is.null(steigerData)) {
    return(buildUnavailableSteigerResult())
  }

  if (!"steiger_dir" %in% names(steigerData)) {
    warning("TwoSampleMR::steiger_filtering() did not return steiger_dir. Returning NA result.")
    return(buildUnavailableSteigerResult())
  }

  steigerPass <- !is.na(steigerData$steiger_dir) & steigerData$steiger_dir
  nRemoved <- sum(!steigerPass, na.rm = TRUE)

  if (!any(steigerPass, na.rm = TRUE)) {
    warning("All SNPs failed Steiger filter. Steiger-filtered estimate not available. Investigate instrument validity.")
    result <- buildUnavailableSteigerResult()
    result$n_removed <- nRemoved
    result$n_remaining <- 0L
    if ("steiger_pval" %in% names(steigerData)) {
      result$min_steiger_pval <- min(steigerData$steiger_pval, na.rm = TRUE)
    }
    return(result)
  }

  filtered <- steigerData[steigerPass, , drop = FALSE]
  ivwResult <- computeIVWTwoSampleMR(filtered)
  ivwResult$method <- "Steiger-filtered IVW"
  ivwResult$n_removed <- nRemoved
  ivwResult$n_remaining <- nrow(filtered)
  ivwResult$min_steiger_pval <- if ("steiger_pval" %in% names(steigerData)) {
    min(steigerData$steiger_pval, na.rm = TRUE)
  } else {
    NA_real_
  }

  ivwResult
}


#' @keywords internal
buildUnavailableSteigerResult <- function() {
  data.frame(
    method = "Steiger-filtered IVW",
    beta_MR = NA_real_,
    se_MR = NA_real_,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    pval = NA_real_,
    n_removed = NA_integer_,
    n_remaining = NA_integer_,
    stringsAsFactors = FALSE
  )
}


#' @keywords internal
computeLeaveOneOutTwoSampleMR <- function(harmonisedData) {
  loo <- TwoSampleMR::mr_leaveoneout(harmonisedData, method = TwoSampleMR::mr_ivw)
  loo <- loo[loo$SNP != "All", , drop = FALSE]

  data.frame(
    snp_removed = loo$SNP,
    beta_MR = loo$b,
    se_MR = loo$se,
    pval = loo$p,
    stringsAsFactors = FALSE
  )
}


#' @keywords internal
computeIVW <- function(perSnpEstimates) {
  weights <- 1 / (perSnpEstimates$se_ZY^2)
  betaMR <- sum(weights * perSnpEstimates$beta_ZY * perSnpEstimates$beta_ZX) /
    sum(weights * perSnpEstimates$beta_ZX^2)
  seMR <- sqrt(1 / sum(weights * perSnpEstimates$beta_ZX^2))
  zStat <- betaMR / seMR
  pval <- 2 * pnorm(-abs(zStat))

  # Cochran's Q for heterogeneity
  ratioSummary <- computeRatioSummary(perSnpEstimates)
  ratioEstimates <- ratioSummary$ratioEstimate
  ratioSE <- ratioSummary$ratioSe
  ratioWeights <- 1 / ratioSE^2
  cochranQ <- sum(ratioWeights * (ratioEstimates - betaMR)^2)
  cochranQPval <- pchisq(cochranQ, df = nrow(perSnpEstimates) - 1,
                           lower.tail = FALSE)

  data.frame(
    method = "IVW",
    beta_MR = betaMR,
    se_MR = seMR,
    ci_lower = betaMR - 1.96 * seMR,
    ci_upper = betaMR + 1.96 * seMR,
    pval = pval,
    cochran_Q = cochranQ,
    cochran_Q_pval = cochranQPval,
    stringsAsFactors = FALSE
  )
}


#' @keywords internal
computeMREgger <- function(perSnpEstimates) {
  # Bowden et al. (2015) recommend orienting all instruments so the
  # SNP-exposure association is positive; otherwise the Egger intercept depends
  # on arbitrary allele coding.
  perSnpEstimates <- orientToPositiveExposure(perSnpEstimates)
  weights <- 1 / (perSnpEstimates$se_ZY^2)

  # Weighted regression with intercept: beta_ZY = alpha + beta * beta_ZX
  # Note: lm() with weights produces a multiplicative random-effects SE that
  # incorporates the residual variance (sigma-hat-squared). This matches the
  # default in TwoSampleMR. The fixed-effect version would assume sigma^2 = 1.
  fit <- lm(beta_ZY ~ beta_ZX,
            data = perSnpEstimates,
            weights = weights)

  coefs <- suppressWarnings(summary(fit)$coefficients)
  betaMR <- coefs["beta_ZX", "Estimate"]
  seMR <- coefs["beta_ZX", "Std. Error"]
  pval <- coefs["beta_ZX", "Pr(>|t|)"]
  intercept <- coefs["(Intercept)", "Estimate"]
  interceptSE <- coefs["(Intercept)", "Std. Error"]
  interceptPval <- coefs["(Intercept)", "Pr(>|t|)"]

  data.frame(
    method = "MR-Egger",
    beta_MR = betaMR,
    se_MR = seMR,
    ci_lower = betaMR - 1.96 * seMR,
    ci_upper = betaMR + 1.96 * seMR,
    pval = pval,
    intercept = intercept,
    intercept_se = interceptSE,
    intercept_pval = interceptPval,
    stringsAsFactors = FALSE
  )
}


#' @keywords internal
computeWeightedMedian <- function(perSnpEstimates, nBoot = 1000) {
  ratioSummary <- computeRatioSummary(perSnpEstimates)
  ratioEstimates <- ratioSummary$ratioEstimate
  ratioSE <- ratioSummary$ratioSe
  ratioWeights <- (1 / ratioSE^2) / sum(1 / ratioSE^2)

  # Bowden et al. (2016): estimate the 50th percentile of the weighted ratio
  # distribution, using linear interpolation on the weighted empirical CDF.
  betaMR <- weightedMedian(ratioEstimates, ratioWeights)

  # Parametric bootstrap on the ratio scale, using the same first-order ratio
  # standard errors that define the weighted-median weights in Bowden et al.
  hadSeed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (hadSeed) {
    oldSeed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (hadSeed) {
      assign(".Random.seed", oldSeed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(42)
  bootEstimates <- numeric(nBoot)
  for (b in seq_len(nBoot)) {
    bootRatio <- stats::rnorm(
      n = length(ratioEstimates),
      mean = ratioEstimates,
      sd = ratioSE
    )
    bootEstimates[b] <- weightedMedian(bootRatio, ratioWeights)
  }
  seMR <- sd(bootEstimates)
  zStat <- betaMR / seMR
  pval <- 2 * pnorm(-abs(zStat))

  data.frame(
    method = "Weighted Median",
    beta_MR = betaMR,
    se_MR = seMR,
    ci_lower = betaMR - 1.96 * seMR,
    ci_upper = betaMR + 1.96 * seMR,
    pval = pval,
    stringsAsFactors = FALSE
  )
}


#' @keywords internal
weightedMedian <- function(values, weights) {
  ord <- order(values)
  values <- values[ord]
  weights <- weights[ord] / sum(weights)

  # Bowden et al. define the weighted median as the weighted empirical 50th
  # percentile. Using midpoint positions enables linear interpolation between
  # adjacent ratio estimates when the median falls inside a weight interval.
  midpoint <- cumsum(weights) - 0.5 * weights
  stats::approx(
    x = midpoint,
    y = values,
    xout = 0.5,
    method = "linear",
    ties = "ordered",
    rule = 2
  )$y
}


#' @keywords internal
computeSteiger <- function(perSnpEstimates, outcomeSampleSize = NULL,
                            exposureSampleSize = NULL,
                            outcomeType = "binary") {
  if (identical(outcomeType, "binary")) {
    warning(
      "Steiger filtering is not implemented for binary outcomes because logistic-regression coefficients cannot be converted to correlations without additional scale assumptions. Returning NA result."
    )
    return(data.frame(
      method = "Steiger-filtered IVW",
      beta_MR = NA_real_,
      se_MR = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      pval = NA_real_,
      n_removed = NA_integer_,
      n_remaining = NA_integer_,
      stringsAsFactors = FALSE
    ))
  }

  if (is.null(outcomeSampleSize) || is.null(exposureSampleSize)) {
    warning(
      "Steiger filtering requires both outcomeSampleSize and exposureSampleSize. Returning NA result."
    )
    return(data.frame(
      method = "Steiger-filtered IVW",
      beta_MR = NA_real_,
      se_MR = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      pval = NA_real_,
      n_removed = NA_integer_,
      n_remaining = NA_integer_,
      stringsAsFactors = FALSE
    ))
  }

  checkmate::assertCount(outcomeSampleSize, positive = TRUE)
  checkmate::assertCount(exposureSampleSize, positive = TRUE)
  nSnps <- nrow(perSnpEstimates)

  # Hemani et al. (2017): compare the absolute SNP-exposure and SNP-outcome
  # correlations, estimated from the summary-statistic t values and sample
  # sizes. In the two-sample setting, this reduces to a Steiger comparison of
  # two independent correlations, implemented here with the standard Fisher-z
  # test.
  rExposure <- summaryStatToCorrelation(
    beta = perSnpEstimates$beta_ZX,
    se = perSnpEstimates$se_ZX,
    sampleSize = exposureSampleSize
  )
  rOutcome <- summaryStatToCorrelation(
    beta = perSnpEstimates$beta_ZY,
    se = perSnpEstimates$se_ZY,
    sampleSize = outcomeSampleSize
  )

  fisherDenominator <- sqrt(
    (1 / pmax(exposureSampleSize - 3, 1)) +
      (1 / pmax(outcomeSampleSize - 3, 1))
  )
  steigerZ <- (atanh(abs(rExposure)) - atanh(abs(rOutcome))) / fisherDenominator
  steigerPval <- 2 * pnorm(-abs(steigerZ))
  # Hemani et al. (2017) recommend requiring both the point estimate
  # comparison AND statistical significance to establish causal direction.
  steigerPass <- abs(rExposure) > abs(rOutcome) & steigerPval < 0.05

  nRemoved <- sum(!steigerPass)

  if (nRemoved == nSnps) {
    warning("All SNPs failed Steiger filter. Steiger-filtered estimate not available. Investigate instrument validity.")
    return(data.frame(
      method = "Steiger-filtered IVW",
      beta_MR = NA_real_,
      se_MR = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      pval = NA_real_,
      n_removed = nRemoved,
      n_remaining = 0L,
      min_steiger_pval = min(steigerPval, na.rm = TRUE),
      stringsAsFactors = FALSE
    ))
  }

  if (nRemoved > 0) {
    message(sprintf("    Steiger filter removed %d of %d SNPs.", nRemoved, nSnps))
  }

  filtered <- perSnpEstimates[steigerPass, , drop = FALSE]

  if (nrow(filtered) < 2) {
    ivwResult <- data.frame(
      method = "Steiger-filtered IVW",
      beta_MR = filtered$beta_ZY[1] / filtered$beta_ZX[1],
      se_MR = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      pval = NA_real_,
      stringsAsFactors = FALSE
    )
  } else {
    ivwResult <- computeIVW(filtered)
    ivwResult$method <- "Steiger-filtered IVW"
  }

  ivwResult$n_removed <- nRemoved
  ivwResult$n_remaining <- nrow(filtered)
  ivwResult$min_steiger_pval <- min(steigerPval, na.rm = TRUE)

  ivwResult
}


#' @keywords internal
computeRatioSummary <- function(perSnpEstimates) {
  betaZX <- perSnpEstimates$beta_ZX
  nearZero <- abs(betaZX) < .Machine$double.eps^0.5
  if (any(nearZero)) {
    stop("Cannot compute ratio estimates when beta_ZX is zero or numerically indistinguishable from zero.")
  }

  ratioEstimate <- perSnpEstimates$beta_ZY / betaZX

  # Burgess et al. (2013) and Bowden et al. (2016) use the first-order delta
  # approximation for summarized-data MR: se(ratio_j) ≈ se(betaY_j) / |betaX_j|.
  # This yields the same unstandardized weights as IVW.
  ratioSe <- perSnpEstimates$se_ZY / abs(betaZX)

  list(
    ratioEstimate = ratioEstimate,
    ratioSe = ratioSe
  )
}


#' @keywords internal
orientToPositiveExposure <- function(perSnpEstimates) {
  oriented <- perSnpEstimates
  flip <- oriented$beta_ZX < 0
  oriented$beta_ZX[flip] <- -oriented$beta_ZX[flip]
  oriented$beta_ZY[flip] <- -oriented$beta_ZY[flip]
  oriented
}


#' @keywords internal
summaryStatToCorrelation <- function(beta, se, sampleSize) {
  tStatistic <- beta / se
  df <- pmax(sampleSize - 2, 1)
  correlation <- sign(beta) * sqrt(tStatistic^2 / (tStatistic^2 + df))
  pmin(pmax(correlation, -1 + 1e-12), 1 - 1e-12)
}


#' @keywords internal
computeLeaveOneOut <- function(perSnpEstimates) {
  nSnps <- nrow(perSnpEstimates)
  results <- data.frame(
    snp_removed = character(nSnps),
    beta_MR = numeric(nSnps),
    se_MR = numeric(nSnps),
    pval = numeric(nSnps),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nSnps)) {
    subset <- perSnpEstimates[-i, , drop = FALSE]
    ivw <- computeIVW(subset)
    results$snp_removed[i] <- perSnpEstimates$snp_id[i]
    results$beta_MR[i] <- ivw$beta_MR
    results$se_MR[i] <- ivw$se_MR
    results$pval[i] <- ivw$pval
  }

  results
}


#' @keywords internal
buildSensitivitySummary <- function(results) {
  summaryRows <- list()

  if (!is.null(results$ivw)) {
    summaryRows <- c(summaryRows, list(data.frame(
      method = "IVW",
      beta_MR = results$ivw$beta_MR,
      se_MR = results$ivw$se_MR,
      ci_lower = results$ivw$ci_lower,
      ci_upper = results$ivw$ci_upper,
      pval = results$ivw$pval,
      stringsAsFactors = FALSE
    )))
  }

  if (!is.null(results$mrEgger)) {
    summaryRows <- c(summaryRows, list(data.frame(
      method = "MR-Egger",
      beta_MR = results$mrEgger$beta_MR,
      se_MR = results$mrEgger$se_MR,
      ci_lower = results$mrEgger$ci_lower,
      ci_upper = results$mrEgger$ci_upper,
      pval = results$mrEgger$pval,
      stringsAsFactors = FALSE
    )))
  }

  if (!is.null(results$weightedMedian)) {
    summaryRows <- c(summaryRows, list(data.frame(
      method = "Weighted Median",
      beta_MR = results$weightedMedian$beta_MR,
      se_MR = results$weightedMedian$se_MR,
      ci_lower = results$weightedMedian$ci_lower,
      ci_upper = results$weightedMedian$ci_upper,
      pval = results$weightedMedian$pval,
      stringsAsFactors = FALSE
    )))
  }

  if (!is.null(results$steiger) && !is.na(results$steiger$beta_MR)) {
    summaryRows <- c(summaryRows, list(data.frame(
      method = "Steiger-filtered IVW",
      beta_MR = results$steiger$beta_MR,
      se_MR = results$steiger$se_MR,
      ci_lower = results$steiger$ci_lower,
      ci_upper = results$steiger$ci_upper,
      pval = results$steiger$pval,
      stringsAsFactors = FALSE
    )))
  }

  if (length(summaryRows) > 0) {
    do.call(rbind, summaryRows)
  } else {
    data.frame(
      method = character(0),
      beta_MR = numeric(0),
      se_MR = numeric(0),
      ci_lower = numeric(0),
      ci_upper = numeric(0),
      pval = numeric(0),
      stringsAsFactors = FALSE
    )
  }
}
