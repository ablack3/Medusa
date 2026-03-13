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

#' Assemble Covariate Matrix for Mendelian Randomization
#'
#' @title Covariate assembly using FeatureExtraction
#'
#' @description
#' Leverages the OHDSI FeatureExtraction package to assemble a rich covariate
#' object from OMOP CDM data. This object serves two purposes: (1) providing
#' covariates for the outcome model to control for confounders and population
#' stratification, and (2) exposing the full phenome for the instrument PheWAS
#' diagnostic that checks for pleiotropic associations.
#'
#' If ancestry principal components (PCs) are available, they are merged into
#' the covariate matrix. Ancestry PCs are critical for controlling population
#' stratification in genetic association studies.
#'
#' @param connectionDetails A \code{DatabaseConnector::connectionDetails} object.
#' @param cdmDatabaseSchema Character. Schema containing OMOP CDM tables.
#' @param cohortDatabaseSchema Character. Schema containing the cohort table.
#' @param cohortTable Character. Name of the cohort table.
#' @param outcomeCohortId Integer. Cohort definition ID for the outcome.
#' @param covariateSettings A FeatureExtraction covariate settings object. If
#'   NULL (default), uses \code{\link{createDefaultMRCovariateSettings}}.
#' @param ancestryPCsTable Character or NULL. Name of a table containing
#'   person_id and ancestry principal components (PC1 through PC_K). If NULL,
#'   ancestry PCs are not included.
#' @param ancestryPCsSchema Character or NULL. Schema containing the ancestry
#'   PCs table.
#' @param numAncestryPCs Integer. Number of ancestry PCs to include (1 through
#'   this value). Default is 10.
#'
#' @return A list with class "medusaCovariateData" containing:
#'   \describe{
#'     \item{covariateData}{The \pkg{FeatureExtraction} covariate object
#'       returned by \code{getDbCovariateData()}. Its main tables are
#'       \code{covariateData$covariates} and \code{covariateData$covariateRef}.}
#'     \item{ancestryPCs}{Data frame of ancestry PCs if provided, NULL otherwise.}
#'     \item{settings}{The covariate settings object used.}
#'   }
#'
#' @details
#' Default covariate settings include: conditions in 365-day lookback (binary),
#' drug exposures in 365-day lookback (binary), most recent measurement values,
#' and demographics (age group, sex, index year). These are assembled using
#' standard FeatureExtraction covariate setting objects.
#'
#' @examples
#' \dontrun{
#' covData <- buildMRCovariates(
#'   connectionDetails = connDetails,
#'   cdmDatabaseSchema = "cdm",
#'   cohortDatabaseSchema = "results",
#'   cohortTable = "cohort",
#'   outcomeCohortId = 1234
#' )
#' }
#'
#' @seealso \code{\link{createDefaultMRCovariateSettings}},
#'   \code{\link{buildMRCohort}}, \code{\link{runInstrumentDiagnostics}}
#'
#' @export
buildMRCovariates <- function(connectionDetails,
                              cdmDatabaseSchema,
                              cohortDatabaseSchema,
                              cohortTable,
                              outcomeCohortId,
                              covariateSettings = NULL,
                              ancestryPCsTable = NULL,
                              ancestryPCsSchema = NULL,
                              numAncestryPCs = 10) {
  # Input validation
  checkmate::assertClass(connectionDetails, "connectionDetails")
  checkmate::assertString(cdmDatabaseSchema)
  checkmate::assertString(cohortDatabaseSchema)
  checkmate::assertString(cohortTable)
  checkmate::assertCount(outcomeCohortId, positive = TRUE)
  checkmate::assertCount(numAncestryPCs, positive = TRUE)

  if (xor(is.null(ancestryPCsTable), is.null(ancestryPCsSchema))) {
    stop("Both 'ancestryPCsTable' and 'ancestryPCsSchema' must be provided together, or both must be NULL.")
  }

  if (!requireNamespace("DatabaseConnector", quietly = TRUE)) {
    stop("Package 'DatabaseConnector' is required for buildMRCovariates(). ",
         "Install it with: remotes::install_github('ohdsi/DatabaseConnector')",
         call. = FALSE)
  }
  if (!requireNamespace("FeatureExtraction", quietly = TRUE)) {
    stop("Package 'FeatureExtraction' is required for buildMRCovariates(). ",
         "Install it with: remotes::install_github('ohdsi/FeatureExtraction')",
         call. = FALSE)
  }

  if (is.null(covariateSettings)) {
    covariateSettings <- createDefaultMRCovariateSettings()
  }

  connection <- DatabaseConnector::connect(connectionDetails)
  on.exit(DatabaseConnector::disconnect(connection), add = TRUE)

  # Extract covariates via FeatureExtraction
  message("Extracting covariates via FeatureExtraction...")
  covariateData <- FeatureExtraction::getDbCovariateData(
    connection = connection,
    cdmDatabaseSchema = cdmDatabaseSchema,
    cohortDatabaseSchema = cohortDatabaseSchema,
    cohortTable = cohortTable,
    cohortId = outcomeCohortId,
    covariateSettings = covariateSettings,
    aggregated = FALSE
  )

  # Andromeda-backed tables need collect/pull to compute length(unique())
  nCovariates <- tryCatch(
    nrow(as.data.frame(covariateData$covariateRef)),
    error = function(e) NA_integer_
  )
  nPersons <- tryCatch(
    length(unique(as.data.frame(covariateData$covariates)$rowId)),
    error = function(e) NA_integer_
  )
  message(sprintf("Extracted %s covariates for %s persons.",
                  if (is.na(nCovariates)) "unknown" else as.character(nCovariates),
                  if (is.na(nPersons)) "unknown" else as.character(nPersons)))

  # Extract ancestry PCs if available
  ancestryPCs <- NULL
  if (!is.null(ancestryPCsTable) && !is.null(ancestryPCsSchema)) {
    message(sprintf("Extracting %d ancestry principal components...", numAncestryPCs))
    pcCols <- paste(paste0("pc", seq_len(numAncestryPCs)), collapse = ", ")
    pcSql <- sprintf(
      "SELECT person_id, %s FROM %s.%s",
      pcCols, ancestryPCsSchema, ancestryPCsTable
    )
    pcSql <- SqlRender::translateSql(pcSql, targetDialect = connectionDetails$dbms)$sql
    ancestryPCs <- DatabaseConnector::querySql(connection, pcSql,
                                                snakeCaseToCamelCase = TRUE)
    message(sprintf("Retrieved ancestry PCs for %d persons.", nrow(ancestryPCs)))
  }

  result <- list(
    covariateData = covariateData,
    ancestryPCs = ancestryPCs,
    settings = covariateSettings
  )
  class(result) <- "medusaCovariateData"

  result
}


#' Create Default Covariate Settings for MR Analysis
#'
#' @title Default FeatureExtraction settings for Medusa
#'
#' @description
#' Creates a FeatureExtraction covariate settings object tailored for Mendelian
#' Randomization analysis. Includes demographics, conditions, drug exposures,
#' and measurements in standard lookback windows.
#'
#' @return A FeatureExtraction \code{covariateSettings} object.
#'
#' @examples
#' settings <- createDefaultMRCovariateSettings()
#'
#' @seealso \code{\link{buildMRCovariates}}
#'
#' @export
createDefaultMRCovariateSettings <- function() {
  if (!requireNamespace("FeatureExtraction", quietly = TRUE)) {
    stop("Package 'FeatureExtraction' is required for createDefaultMRCovariateSettings(). ",
         "Install it with: remotes::install_github('ohdsi/FeatureExtraction')",
         call. = FALSE)
  }
  FeatureExtraction::createCovariateSettings(
    useDemographicsGender = TRUE,
    useDemographicsAge = TRUE,
    useDemographicsAgeGroup = TRUE,
    useDemographicsIndexYear = TRUE,
    useConditionGroupEraLongTerm = TRUE,
    useDrugGroupEraLongTerm = TRUE,
    useMeasurementValueLongTerm = TRUE,
    longTermStartDays = -365,
    endDays = 0
  )
}
