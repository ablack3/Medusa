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
#'   columns: snp_id, beta_ZY, se_ZY, beta_ZX, se_ZX.
#' @param methods Character vector of methods to run. Default is
#'   \code{c("IVW", "MREgger", "WeightedMedian", "Steiger", "LeaveOneOut")}.
#' @param outcomeSampleSize Optional integer. Total sample size for outcome
#'   GWAS/analysis. Needed for Steiger filtering.
#' @param exposureSampleSize Optional integer. Total sample size for exposure
#'   GWAS. Needed for Steiger filtering.
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
#' intercept. A non-zero intercept indicates directional pleiotropy. The slope
#' is a consistent causal estimate even under directional pleiotropy, provided
#' the InSIDE assumption holds.
#'
#' \strong{Weighted Median}: Bootstrap-based weighted median of per-SNP Wald
#' ratio estimates. Produces a consistent estimate if fewer than 50% of the
#' weight comes from invalid instruments.
#'
#' \strong{Steiger}: Tests whether each SNP explains more variance in the
#' exposure than the outcome (the expected causal direction). SNPs failing
#' this test are removed and IVW is re-run.
#'
#' \strong{Leave-One-Out}: Drops each SNP in turn and recomputes the IVW
#' estimate. Identifies influential outlier instruments.
#'
#' @references
#' Bowden, J., Davey Smith, G., & Burgess, S. (2015). Mendelian randomization
#' with invalid instruments: effect estimation and bias detection through Egger
#' regression. \emph{International Journal of Epidemiology}, 44(2), 512-525.
#'
#' Bowden, J., et al. (2016). Consistent estimation in Mendelian randomization
#' with some invalid instruments using a weighted median estimator.
#' \emph{Genetic Epidemiology}, 40(4), 304-314.
#'
#' Hemani, G., Tilling, K., & Davey Smith, G. (2017). Orienting the causal
#' relationship between imprecisely measured traits using GWAS summary data.
#' \emph{PLoS Genetics}, 13(11), e1007081.
#'
#' @examples
#' # Simulate per-SNP estimates
#' set.seed(42)
#' nSnps <- 10
#' betaZX <- rnorm(nSnps, 0.3, 0.05)
#' betaZY <- 0.5 * betaZX + rnorm(nSnps, 0, 0.02)
#' perSnp <- data.frame(
#'   snp_id = paste0("rs", 1:nSnps),
#'   beta_ZY = betaZY,
#'   se_ZY = rep(0.02, nSnps),
#'   beta_ZX = betaZX,
#'   se_ZX = rep(0.05, nSnps)
#' )
#' results <- runSensitivityAnalyses(perSnp)
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
                                   exposureSampleSize = NULL) {
  # Input validation
  checkmate::assertDataFrame(perSnpEstimates, min.rows = 1)
  requiredCols <- c("snp_id", "beta_ZY", "se_ZY", "beta_ZX", "se_ZX")
  checkmate::assertSubset(requiredCols, names(perSnpEstimates))
  checkmate::assertSubset(methods,
                          c("IVW", "MREgger", "WeightedMedian",
                            "Steiger", "LeaveOneOut"))

  nSnps <- nrow(perSnpEstimates)
  message(sprintf("Running sensitivity analyses with %d SNPs...", nSnps))

  results <- list()

  # --- IVW ---
  if ("IVW" %in% methods) {
    message("  IVW...")
    results$ivw <- computeIVW(perSnpEstimates)
  }

  # --- MR-Egger ---
  if ("MREgger" %in% methods) {
    if (nSnps < 3) {
      message("  MR-Egger: skipped (requires >= 3 SNPs).")
      results$mrEgger <- NULL
    } else {
      message("  MR-Egger...")
      results$mrEgger <- computeMREgger(perSnpEstimates)
    }
  }

  # --- Weighted Median ---
  if ("WeightedMedian" %in% methods) {
    if (nSnps < 3) {
      message("  Weighted Median: skipped (requires >= 3 SNPs).")
      results$weightedMedian <- NULL
    } else {
      message("  Weighted Median...")
      results$weightedMedian <- computeWeightedMedian(perSnpEstimates)
    }
  }

  # --- Steiger ---
  if ("Steiger" %in% methods) {
    message("  Steiger filtering...")
    results$steiger <- computeSteiger(perSnpEstimates, outcomeSampleSize,
                                      exposureSampleSize)
  }

  # --- Leave-One-Out ---
  if ("LeaveOneOut" %in% methods) {
    if (nSnps < 3) {
      message("  Leave-One-Out: skipped (requires >= 3 SNPs).")
      results$leaveOneOut <- NULL
    } else {
      message("  Leave-One-Out...")
      results$leaveOneOut <- computeLeaveOneOut(perSnpEstimates)
    }
  }

  # --- Summary table ---
  results$summary <- buildSensitivitySummary(results)

  class(results) <- "medusaSensitivity"

  message("Sensitivity analyses complete.")
  results
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
  ratioEstimates <- perSnpEstimates$beta_ZY / perSnpEstimates$beta_ZX
  ratioSE <- perSnpEstimates$se_ZY / abs(perSnpEstimates$beta_ZX)
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
  weights <- 1 / (perSnpEstimates$se_ZY^2)

  # Weighted regression with intercept: beta_ZY = alpha + beta * beta_ZX
  fit <- lm(beta_ZY ~ beta_ZX,
            data = perSnpEstimates,
            weights = weights)

  coefs <- summary(fit)$coefficients
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
  # Per-SNP ratio estimates
  ratioEstimates <- perSnpEstimates$beta_ZY / perSnpEstimates$beta_ZX
  ratioWeights <- (perSnpEstimates$beta_ZX^2) / (perSnpEstimates$se_ZY^2)
  ratioWeights <- ratioWeights / sum(ratioWeights)

  # Weighted median
  betaMR <- weightedMedian(ratioEstimates, ratioWeights)

  # Bootstrap SE
  set.seed(42)
  bootEstimates <- numeric(nBoot)
  for (b in seq_len(nBoot)) {
    bootBetaZY <- perSnpEstimates$beta_ZY + rnorm(nrow(perSnpEstimates)) *
      perSnpEstimates$se_ZY
    bootBetaZX <- perSnpEstimates$beta_ZX + rnorm(nrow(perSnpEstimates)) *
      perSnpEstimates$se_ZX
    bootRatio <- bootBetaZY / bootBetaZX
    bootWeights <- (bootBetaZX^2) / (perSnpEstimates$se_ZY^2)
    bootWeights <- bootWeights / sum(bootWeights)
    bootEstimates[b] <- weightedMedian(bootRatio, bootWeights)
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
  weights <- weights[ord]
  cumWeights <- cumsum(weights)
  idx <- which(cumWeights >= 0.5)[1]
  values[idx]
}


#' @keywords internal
computeSteiger <- function(perSnpEstimates, outcomeSampleSize = NULL,
                            exposureSampleSize = NULL) {
  nSnps <- nrow(perSnpEstimates)

  # R-squared for exposure and outcome
  r2Exposure <- (perSnpEstimates$beta_ZX^2) /
    (perSnpEstimates$beta_ZX^2 + perSnpEstimates$se_ZX^2)
  r2Outcome <- (perSnpEstimates$beta_ZY^2) /
    (perSnpEstimates$beta_ZY^2 + perSnpEstimates$se_ZY^2)

  # Steiger test: SNP should explain more variance in exposure than outcome
  steigerPass <- r2Exposure > r2Outcome

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

  ivwResult
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
