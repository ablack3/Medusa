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

#' Pool Profile Log-Likelihoods Across Sites
#'
#' @title Federated one-shot likelihood pooling
#'
#' @description
#' Combines profile log-likelihood vectors from multiple federated sites by
#' pointwise summation. This is the core of the one-shot federated analysis:
#' summing log-likelihoods is equivalent to multiplying likelihoods, yielding
#' the joint likelihood under independence across sites.
#'
#' Before summing, the function validates that all sites used identical beta
#' grid vectors. If grids differ, linear interpolation to a common grid over
#' the overlapping beta range is applied with a warning. Pooling never
#' extrapolates beyond a site's evaluated range.
#'
#' @param siteProfileList Named list of site profile objects, each the output
#'   of \code{\link{fitOutcomeModel}}.
#' @param validateGridAlignment Logical. If TRUE (default), checks that all
#'   sites used the same betaGrid and warns if interpolation is needed.
#'
#' @return A list with class "medusaCombinedProfile" containing:
#'   \describe{
#'     \item{betaGrid}{The common beta grid used for pooling.}
#'     \item{logLikProfile}{Numeric vector of summed log-likelihood values.}
#'     \item{siteContributions}{Data frame with siteId, nCases, nControls,
#'       and diagnostic flags per site.}
#'     \item{nSites}{Number of sites pooled.}
#'     \item{totalCases}{Total number of cases across sites.}
#'     \item{totalControls}{Total number of controls across sites.}
#'   }
#'
#' @details
#' The pooling operation is mathematically straightforward: for each grid point
#' \eqn{b_k}, the combined log-likelihood is
#' \deqn{\ell_{combined}(b_k) = \sum_{s=1}^{S} \ell_s(b_k)}
#' where \eqn{\ell_s} is the profile log-likelihood from site \eqn{s}.
#'
#' This pooling is exact under the assumption that observations are independent
#' across sites and every site uses the same beta grid. If interpolation is
#' needed, pooling remains a close numerical approximation on the common grid,
#' but only over the beta range evaluated by every site. No iterative
#' communication protocol is needed: each site computes its profile once and
#' shares only the numeric vector.
#'
#' @references
#' Luo, Y., et al. (2022). dPQL: a lossless distributed algorithm for
#' generalized linear mixed model with application to privacy-preserving
#' hospital profiling. \emph{Journal of the American Medical Informatics
#' Association}, 29(8), 1366-1373.
#'
#' @examples
#' profiles <- simulateSiteProfiles(nSites = 3, trueBeta = 0.5)
#' combined <- poolLikelihoodProfiles(profiles)
#' plot(combined$betaGrid, combined$logLikProfile, type = "l",
#'      xlab = "beta_ZY", ylab = "Pooled profile log-likelihood")
#'
#' @seealso \code{\link{fitOutcomeModel}}, \code{\link{computeMREstimate}}
#'
#' @export
poolLikelihoodProfiles <- function(siteProfileList,
                                   validateGridAlignment = TRUE) {
  # Input validation
  checkmate::assertList(siteProfileList, min.len = 1)

  for (i in seq_along(siteProfileList)) {
    validateSiteProfile(siteProfileList[[i]])
  }

  nSites <- length(siteProfileList)
  message(sprintf("Pooling profile likelihoods from %d site(s)...", nSites))

  # Get reference grid from first site
  referenceGrid <- siteProfileList[[1]]$betaGrid
  referenceScoreDefinition <- siteProfileList[[1]]$scoreDefinition

  # Validate grid alignment
  needsInterpolation <- FALSE
  if (nSites > 1) {
    for (i in 2:nSites) {
      siteGrid <- siteProfileList[[i]]$betaGrid
      if (validateGridAlignment &&
          (length(siteGrid) != length(referenceGrid) ||
           !all(abs(siteGrid - referenceGrid) < 1e-10))) {
        needsInterpolation <- TRUE
        warning(sprintf(
          paste0("Site '%s' used different betaGrid. Interpolating to common grid. ",
                 "Consider specifying identical betaGrid at all sites."),
          siteProfileList[[i]]$siteId
        ))
      }

      if (!identicalScoreDefinition(siteProfileList[[i]]$scoreDefinition,
                                    referenceScoreDefinition)) {
        stop(
          sprintf(
            "Site '%s' used a different allele-score definition. All pooled sites must use the same SNP set and score weights.",
            siteProfileList[[i]]$siteId
          )
        )
      }
    }
  }

  # Determine common grid
  if (needsInterpolation) {
    gridLower <- max(vapply(siteProfileList, function(profile) {
      min(profile$betaGrid)
    }, numeric(1)))
    gridUpper <- min(vapply(siteProfileList, function(profile) {
      max(profile$betaGrid)
    }, numeric(1)))
    if (!is.finite(gridLower) || !is.finite(gridUpper) || gridLower >= gridUpper) {
      stop(
        "Site betaGrid ranges do not overlap. Pooling requires a shared beta range."
      )
    }

    gridStep <- min(vapply(siteProfileList, function(profile) {
      stats::median(diff(profile$betaGrid))
    }, numeric(1)))
    commonGrid <- seq(from = gridLower, to = gridUpper, by = gridStep)
    if (tail(commonGrid, 1) < gridUpper - (gridStep / 1000)) {
      commonGrid <- c(commonGrid, gridUpper)
    }
  } else {
    commonGrid <- referenceGrid
  }

  # Sum log-likelihood profiles
  combinedLogLik <- rep(0, length(commonGrid))
  siteContributions <- data.frame(
    siteId = character(nSites),
    nCases = integer(nSites),
    nControls = integer(nSites),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(siteProfileList)) {
    profile <- siteProfileList[[i]]
    siteContributions$siteId[i] <- profile$siteId
    siteContributions$nCases[i] <- profile$nCases
    siteContributions$nControls[i] <- profile$nControls

    if (needsInterpolation &&
        (length(profile$betaGrid) != length(commonGrid) ||
         !all(abs(profile$betaGrid - commonGrid) < 1e-10))) {
      # Interpolate to common grid
      interpolated <- approx(
        x = profile$betaGrid,
        y = profile$logLikProfile,
        xout = commonGrid,
        method = "linear",
        rule = 1,
        ties = "ordered"
      )
      if (anyNA(interpolated$y)) {
        stop(
          sprintf(
            "Site '%s' does not cover the common beta range after interpolation.",
            profile$siteId
          )
        )
      }
      combinedLogLik <- combinedLogLik + interpolated$y
    } else {
      combinedLogLik <- combinedLogLik + profile$logLikProfile
    }
  }

  # Normalize so max is 0 (for numerical stability)
  combinedLogLik <- combinedLogLik - max(combinedLogLik)

  totalCases <- sum(siteContributions$nCases)
  totalControls <- sum(siteContributions$nControls)

  message(sprintf("Pooling complete: %d sites, %d total cases, %d total controls.",
                  nSites, totalCases, totalControls))

  result <- list(
    betaGrid = commonGrid,
    logLikProfile = combinedLogLik,
    siteContributions = siteContributions,
    nSites = nSites,
    totalCases = totalCases,
    totalControls = totalControls,
    scoreDefinition = referenceScoreDefinition
  )
  class(result) <- "medusaCombinedProfile"

  result
}


#' @keywords internal
identicalScoreDefinition <- function(lhs, rhs) {
  if (is.null(lhs) || is.null(rhs)) {
    return(is.null(lhs) && is.null(rhs))
  }

  identical(lhs$snpIds, rhs$snpIds) &&
    isTRUE(all.equal(lhs$scoreWeights, rhs$scoreWeights, tolerance = 1e-10)) &&
    isTRUE(all.equal(lhs$betaZX, rhs$betaZX, tolerance = 1e-10)) &&
    isTRUE(all.equal(lhs$seZX, rhs$seZX, tolerance = 1e-10))
}
