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

#' Export a Site Profile to CSV
#'
#' @description
#' Writes a site profile (output of \code{\link{fitOutcomeModel}}) to
#' human-readable CSV files: one for the profile log-likelihood grid, one for
#' site metadata, a pair of allele-score definition files (weights plus the
#' aggregate score association), and optionally one for per-SNP estimates.
#' These CSV files are the artifacts shared between sites and the coordinator
#' in a federated Medusa analysis.
#'
#' CSV is used instead of binary formats so that every value leaving a site is
#' human-readable and auditable.
#'
#' @param profile A site profile object (output of \code{\link{fitOutcomeModel}}).
#' @param outputDir Character. Directory to write files to. Default is current
#'   working directory.
#' @param prefix Character. Filename prefix. Default is "medusa".
#'
#' @return A named character vector with the paths to the written files
#'   (invisibly).
#'
#' @examples
#' simData <- simulateMRData(n = 500, nSnps = 3, trueEffect = 0.3)
#' profile <- fitOutcomeModel(
#'   cohortData = simData$data,
#'   covariateData = NULL,
#'   instrumentTable = simData$instrumentTable,
#'   betaGrid = seq(-2, 2, by = 0.1),
#'   siteId = "example_site"
#' )
#' \dontrun{
#' exportSiteProfile(profile, outputDir = tempdir())
#' }
#'
#' @seealso \code{\link{importSiteProfile}}, \code{\link{fitOutcomeModel}}
#'
#' @export
exportSiteProfile <- function(profile,
                              outputDir = ".",
                              prefix = "medusa") {
  checkmate::assertList(profile)
  checkmate::assertString(outputDir)
  checkmate::assertString(prefix)

  requiredFields <- c("siteId", "betaGrid", "logLikProfile",
                       "nCases", "nControls", "snpIds")
  missing <- setdiff(requiredFields, names(profile))
  if (length(missing) > 0) {
    stop(sprintf("Profile is missing required fields: %s",
                 paste(missing, collapse = ", ")))
  }

  if (!dir.exists(outputDir)) {
    dir.create(outputDir, recursive = TRUE)
  }

  siteId <- profile$siteId
  paths <- character(0)

  # --- Profile CSV: beta grid and log-likelihood values ---
  profilePath <- file.path(outputDir,
                            sprintf("%s_profile_%s.csv", prefix, siteId))
  profileDf <- data.frame(
    beta = profile$betaGrid,
    log_likelihood = profile$logLikProfile
  )
  utils::write.csv(profileDf, profilePath, row.names = FALSE)
  paths["profile"] <- profilePath

  # --- Metadata CSV: site summary and diagnostic flags ---
  metadataPath <- file.path(outputDir,
                             sprintf("%s_metadata_%s.csv", prefix, siteId))

  flags <- profile$diagnosticFlags %||% list()

  metaDf <- data.frame(
    site_id = siteId,
    n_cases = profile$nCases,
    n_controls = profile$nControls,
    beta_hat = profile$betaHat %||% NA_real_,
    se_hat = profile$seHat %||% NA_real_,
    snp_ids = paste(profile$snpIds, collapse = ";"),
    weak_instruments = flags$weakInstruments %||% NA,
    low_case_count = flags$lowCaseCount %||% NA,
    grid_boundary_mle = flags$gridBoundaryMLE %||% NA,
    stringsAsFactors = FALSE
  )
  utils::write.csv(metaDf, metadataPath, row.names = FALSE)
  paths["metadata"] <- metadataPath

  # --- Score definition CSV: allele-score weights and exposure summary ---
  if (!is.null(profile$scoreDefinition)) {
    scorePath <- file.path(outputDir,
                            sprintf("%s_score_%s.csv", prefix, siteId))
    scoreDef <- profile$scoreDefinition
    scoreDf <- data.frame(
      snp_id = scoreDef$snpIds,
      score_weight = scoreDef$scoreWeights,
      stringsAsFactors = FALSE
    )
    # Store betaZX and seZX as attributes in a one-row summary at the top
    utils::write.csv(scoreDf, scorePath, row.names = FALSE)

    # Write a companion summary file with the aggregate score statistics
    scoreSummaryPath <- file.path(outputDir,
                                   sprintf("%s_score_summary_%s.csv", prefix, siteId))
    scoreSummaryDf <- data.frame(
      beta_zx_score = scoreDef$betaZX,
      se_zx_score = scoreDef$seZX,
      stringsAsFactors = FALSE
    )
    utils::write.csv(scoreSummaryDf, scoreSummaryPath, row.names = FALSE)
    paths["score"] <- scorePath
    paths["score_summary"] <- scoreSummaryPath
  }

  # --- Per-SNP estimates CSV (optional) ---
  if (!is.null(profile$perSnpEstimates) &&
      is.data.frame(profile$perSnpEstimates) &&
      nrow(profile$perSnpEstimates) > 0) {
    perSnpPath <- file.path(outputDir,
                             sprintf("%s_per_snp_%s.csv", prefix, siteId))
    utils::write.csv(profile$perSnpEstimates, perSnpPath, row.names = FALSE)
    paths["per_snp"] <- perSnpPath
  }

  message(sprintf("Exported site profile to:\n  %s",
                  paste(paths, collapse = "\n  ")))
  invisible(paths)
}


#' Import a Site Profile from CSV
#'
#' @description
#' Reads a site profile from CSV files previously written by
#' \code{\link{exportSiteProfile}}. Reconstructs the profile list object
#' that can be passed to \code{\link{poolLikelihoodProfiles}}, restoring the
#' optional allele-score definition and per-SNP summary sidecar files when
#' they are present.
#'
#' @param profilePath Character. Path to the profile CSV file (the file
#'   containing beta and log_likelihood columns).
#' @param metadataPath Character. Path to the metadata CSV file. If NULL
#'   (default), the function infers the path by replacing "_profile_" with
#'   "_metadata_" in \code{profilePath}.
#'
#' @return A list with the same structure as \code{\link{fitOutcomeModel}}
#'   output, suitable for passing to \code{\link{poolLikelihoodProfiles}}.
#'   If the companion \code{*_per_snp_*.csv} file exists, the returned object
#'   also includes \code{perSnpEstimates}.
#'
#' @examples
#' \dontrun{
#' profile <- importSiteProfile("medusa_profile_site_A.csv")
#' }
#'
#' @seealso \code{\link{exportSiteProfile}}, \code{\link{poolLikelihoodProfiles}}
#'
#' @export
importSiteProfile <- function(profilePath,
                              metadataPath = NULL) {
  checkmate::assertString(profilePath)
  if (!file.exists(profilePath)) {
    stop(sprintf("Profile file not found: %s", profilePath))
  }

  # Infer companion file paths from the profile path
  if (is.null(metadataPath)) {
    metadataPath <- sub("_profile_", "_metadata_", profilePath)
  }
  scorePath <- sub("_profile_", "_score_", profilePath)
  scoreSummaryPath <- sub("_profile_", "_score_summary_", profilePath)
  perSnpPath <- sub("_profile_", "_per_snp_", profilePath)

  # Read profile grid
  profileDf <- utils::read.csv(profilePath, stringsAsFactors = FALSE)
  requiredCols <- c("beta", "log_likelihood")
  missingCols <- setdiff(requiredCols, names(profileDf))
  if (length(missingCols) > 0) {
    stop(sprintf("Profile CSV is missing required columns: %s",
                 paste(missingCols, collapse = ", ")))
  }

  # Read metadata
  siteId <- "unknown"
  nCases <- NA_integer_
  nControls <- NA_integer_
  betaHat <- NA_real_
  seHat <- NA_real_
  snpIds <- character(0)
  diagnosticFlags <- list(
    weakInstruments = FALSE,
    lowCaseCount = FALSE,
    gridBoundaryMLE = FALSE
  )

  if (file.exists(metadataPath)) {
    metaDf <- utils::read.csv(metadataPath, stringsAsFactors = FALSE)

    if ("site_id" %in% names(metaDf)) siteId <- metaDf$site_id[1]
    if ("n_cases" %in% names(metaDf)) nCases <- as.integer(metaDf$n_cases[1])
    if ("n_controls" %in% names(metaDf)) nControls <- as.integer(metaDf$n_controls[1])
    if ("beta_hat" %in% names(metaDf)) betaHat <- metaDf$beta_hat[1]
    if ("se_hat" %in% names(metaDf)) seHat <- metaDf$se_hat[1]
    if ("snp_ids" %in% names(metaDf)) {
      snpIds <- strsplit(as.character(metaDf$snp_ids[1]), ";")[[1]]
    }
    if ("weak_instruments" %in% names(metaDf)) {
      diagnosticFlags$weakInstruments <- as.logical(metaDf$weak_instruments[1])
    }
    if ("low_case_count" %in% names(metaDf)) {
      diagnosticFlags$lowCaseCount <- as.logical(metaDf$low_case_count[1])
    }
    if ("grid_boundary_mle" %in% names(metaDf)) {
      diagnosticFlags$gridBoundaryMLE <- as.logical(metaDf$grid_boundary_mle[1])
    }
  } else {
    warning(sprintf("Metadata file not found: %s. Using defaults.", metadataPath))
  }

  # Read score definition if present
  scoreDefinition <- NULL
  if (file.exists(scorePath)) {
    scoreDf <- utils::read.csv(scorePath, stringsAsFactors = FALSE)
    scoreDefinition <- list(
      snpIds = scoreDf$snp_id,
      scoreWeights = scoreDf$score_weight,
      betaZX = NA_real_,
      seZX = NA_real_
    )
    if (file.exists(scoreSummaryPath)) {
      summDf <- utils::read.csv(scoreSummaryPath, stringsAsFactors = FALSE)
      if ("beta_zx_score" %in% names(summDf)) {
        scoreDefinition$betaZX <- summDf$beta_zx_score[1]
      }
      if ("se_zx_score" %in% names(summDf)) {
        scoreDefinition$seZX <- summDf$se_zx_score[1]
      }
    }
  }

  profile <- list(
    siteId = siteId,
    betaGrid = profileDf$beta,
    logLikProfile = profileDf$log_likelihood,
    nCases = nCases,
    nControls = nControls,
    snpIds = snpIds,
    diagnosticFlags = diagnosticFlags,
    betaHat = betaHat,
    seHat = seHat,
    scoreDefinition = scoreDefinition
  )
  if (file.exists(perSnpPath)) {
    profile$perSnpEstimates <- utils::read.csv(
      perSnpPath,
      stringsAsFactors = FALSE
    )
  }
  class(profile) <- "medusaSiteProfile"

  message(sprintf("Imported site profile '%s' (%d grid points).",
                  siteId, length(profile$betaGrid)))
  profile
}
