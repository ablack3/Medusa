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

#' Compute Mendelian Randomization Causal Estimate
#'
#' @title Wald ratio MR estimate from pooled profile likelihood
#'
#' @description
#' Takes the pooled profile log-likelihood from \code{\link{poolLikelihoodProfiles}}
#' and the instrument table to compute the MR causal estimate via the Wald ratio.
#' The SNP-outcome association (beta_ZY) is estimated from the profile
#' likelihood peak, and the causal estimate is beta_MR = beta_ZY / beta_ZX.
#' For multi-SNP analyses, beta_ZX is the association of the exact weighted
#' allele score used at each site, not a simple mean of the component SNP
#' effects. Uncertainty is propagated using the delta method.
#'
#' Confidence intervals are computed both from the profile likelihood (likelihood
#' ratio method) and via the delta method (normal approximation). The likelihood-
#' based CI is preferred as it does not assume normality.
#'
#' @param combinedProfile Output of \code{\link{poolLikelihoodProfiles}}.
#' @param instrumentTable Data frame. Output of \code{\link{getMRInstruments}}.
#'   Used for beta_ZX and se_ZX values.
#' @param ciLevel Numeric. Confidence interval level. Default is 0.95.
#'
#' @return A list with class "medusaMREstimate" containing:
#'   \describe{
#'     \item{betaZY}{Point estimate of SNP-outcome association (profile MLE).}
#'     \item{seZY}{Standard error of betaZY from profile curvature.}
#'     \item{betaMR}{Causal effect estimate: betaZY / betaZX.}
#'     \item{seMR}{Standard error of betaMR via delta method.}
#'     \item{ciLower}{Lower confidence limit for betaMR.}
#'     \item{ciUpper}{Upper confidence limit for betaMR.}
#'     \item{pValue}{Two-sided p-value for betaMR.}
#'     \item{oddsRatio}{Exp(betaMR) — causal odds ratio for binary outcomes.}
#'     \item{orCiLower}{Lower CI for odds ratio.}
#'     \item{orCiUpper}{Upper CI for odds ratio.}
#'     \item{ciLevel}{The confidence level used.}
#'     \item{betaZX}{The exposure effect of the fitted allele score.}
#'     \item{seZX}{The SE of betaZX used.}
#'     \item{nInstruments}{Number of instruments used.}
#'     \item{combinedProfile}{The full combined profile for plotting.}
#'   }
#'
#' @details
#' \strong{Wald ratio}: For a single instrument, beta_MR = beta_ZY / beta_ZX.
#' When multiple instruments are pooled into a single allele score, the
#' denominator must be the exposure effect of that same score:
#' \deqn{\beta_{ZX, score} = \sum_j w_j \gamma_j}
#' where \eqn{w_j} are the fixed score weights used at the sites and
#' \eqn{\gamma_j} are the SNP-exposure coefficients.
#'
#' \strong{Delta method SE}: se_MR = sqrt( (se_ZY/beta_ZX)^2 +
#' (beta_ZY * se_ZX / beta_ZX^2)^2 ).
#'
#' \strong{Likelihood-based CI}: The region of betaZY values where the
#' log-likelihood does not drop below peak - qchisq(ciLevel, df=1)/2.
#' This CI is then transformed to the MR scale via division by beta_ZX.
#'
#' @references
#' Burgess, S., Butterworth, A., & Thompson, S. G. (2013). Mendelian
#' randomization analysis with multiple genetic variants using summarized
#' data. \emph{Genetic Epidemiology}, 37(7), 658-665.
#' doi:10.1002/gepi.21758.
#' Open access: https://pmc.ncbi.nlm.nih.gov/articles/PMC4377079/
#'
#' @examples
#' profiles <- simulateSiteProfiles(nSites = 3, trueBeta = 0.5)
#' combined <- poolLikelihoodProfiles(profiles)
#' instruments <- simulateInstrumentTable(nSnps = 5)
#' estimate <- computeMREstimate(combined, instruments)
#' print(estimate$betaMR)
#'
#' @seealso \code{\link{poolLikelihoodProfiles}}, \code{\link{runSensitivityAnalyses}},
#'   \code{\link{generateMRReport}}
#'
#' @export
computeMREstimate <- function(combinedProfile,
                              instrumentTable,
                              ciLevel = 0.95) {
  # Input validation
  checkmate::assertList(combinedProfile)
  checkmate::assertSubset(c("betaGrid", "logLikProfile"), names(combinedProfile))
  validateInstrumentTable(instrumentTable)
  checkmate::assertNumber(ciLevel, lower = 0.5, upper = 0.999)

  betaGrid <- combinedProfile$betaGrid
  logLikProfile <- combinedProfile$logLikProfile

  # Check for multiple local maxima
  # Detect by finding sign changes in the first derivative
  diffs <- diff(logLikProfile)
  signChanges <- sum(diffs[-length(diffs)] > 0 & diffs[-1] < 0)
  if (signChanges > 1) {
    warning("Profile log-likelihood has multiple local maxima. Results may be unreliable.")
  }

  # Step 1: Find MLE of beta_ZY
  peakIdx <- which.max(logLikProfile)
  betaZYHat <- betaGrid[peakIdx]

  # Step 2: Estimate SE from profile curvature
  seZY <- estimateSEFromProfile(betaGrid, logLikProfile)

  # Step 3: Likelihood-based CI for beta_ZY
  chiSqThreshold <- qchisq(ciLevel, df = 1) / 2
  ciMask <- logLikProfile >= (max(logLikProfile) - chiSqThreshold)
  ciIndices <- which(ciMask)

  if (length(ciIndices) > 0) {
    betaZYCILower <- betaGrid[min(ciIndices)]
    betaZYCIUpper <- betaGrid[max(ciIndices)]
  } else {
    betaZYCILower <- betaZYHat - qnorm(1 - (1 - ciLevel) / 2) * seZY
    betaZYCIUpper <- betaZYHat + qnorm(1 - (1 - ciLevel) / 2) * seZY
  }

  # Step 4: Compute the exposure association of the exact allele score fitted
  # at each site. Following Burgess et al. (2013), if S = sum_j w_j G_j, then
  # gamma_S = sum_j w_j gamma_j.
  scoreWeights <- computeAlleleScoreWeights(instrumentTable)
  betaZX <- sum(scoreWeights * instrumentTable$beta_ZX)

  # Under the standard LD-pruned approximation for uncorrelated variants,
  # Var(gamma_S) = sum_j w_j^2 Var(gamma_j).
  seZX <- sqrt(sum((scoreWeights^2) * (instrumentTable$se_ZX^2)))

  if (!is.finite(betaZX) || abs(betaZX) < .Machine$double.eps^0.5) {
    stop("The fitted allele score has an exposure association too close to zero for a stable Wald ratio.")
  }

  # Step 5: Wald ratio
  betaMR <- betaZYHat / betaZX

  # Step 6: Delta method SE
  seMR <- sqrt(
    (seZY / betaZX)^2 +
      (betaZYHat * seZX / betaZX^2)^2
  )

  # Step 7: Transform CI to MR scale
  ciLower <- betaZYCILower / betaZX
  ciUpper <- betaZYCIUpper / betaZX
  # Ensure correct ordering if betaZX is negative
  if (betaZX < 0) {
    tmp <- ciLower
    ciLower <- ciUpper
    ciUpper <- tmp
  }

  # Step 8: P-value
  zStat <- betaMR / seMR
  pValue <- 2 * pnorm(-abs(zStat))

  # Odds ratio
  oddsRatio <- exp(betaMR)
  orCiLower <- exp(ciLower)
  orCiUpper <- exp(ciUpper)

  message(sprintf("MR estimate: beta = %.4f (%.0f%% CI: %.4f, %.4f), p = %.2e",
                  betaMR, ciLevel * 100, ciLower, ciUpper, pValue))
  message(sprintf("Odds ratio: %.3f (%.0f%% CI: %.3f, %.3f)",
                  oddsRatio, ciLevel * 100, orCiLower, orCiUpper))

  result <- list(
    betaZY = betaZYHat,
    seZY = seZY,
    betaMR = betaMR,
    seMR = seMR,
    ciLower = ciLower,
    ciUpper = ciUpper,
    pValue = pValue,
    oddsRatio = oddsRatio,
    orCiLower = orCiLower,
    orCiUpper = orCiUpper,
    ciLevel = ciLevel,
    betaZX = betaZX,
    seZX = seZX,
    nInstruments = nrow(instrumentTable),
    combinedProfile = combinedProfile
  )
  class(result) <- "medusaMREstimate"

  result
}
