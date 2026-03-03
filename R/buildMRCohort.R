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

#' Build Study Cohort for Mendelian Randomization at a Single Site
#'
#' @title Site-level cohort extraction from OMOP CDM
#'
#' @description
#' Runs locally at each OMOP CDM site. Executes parameterized SQL to extract
#' outcome cohort data and genotype data, performs allele harmonization to ensure
#' genotype coding matches the instrument table, and returns a local R data frame
#' suitable for downstream analysis. The returned data never leaves the site.
#'
#' The function queries: PERSON (age, sex), CONDITION_OCCURRENCE (outcome),
#' OBSERVATION_PERIOD (eligibility), and the genomic linkage table (genotypes).
#' All SQL is rendered and translated via SqlRender for cross-dialect compatibility.
#'
#' @param connectionDetails A \code{DatabaseConnector::connectionDetails} object
#'   specifying the database connection.
#' @param cdmDatabaseSchema Character. Schema containing the OMOP CDM tables.
#' @param cohortDatabaseSchema Character. Schema containing the cohort table.
#' @param cohortTable Character. Name of the cohort table.
#' @param outcomeCohortId Integer. Cohort definition ID for the outcome of interest
#'   (e.g., incident colorectal cancer).
#' @param instrumentTable Data frame. Output of \code{\link{getMRInstruments}} or
#'   \code{\link{createInstrumentTable}} containing the instrument SNPs.
#' @param genomicLinkageSchema Character. Schema containing the genomic linkage table.
#' @param genomicLinkageTable Character. Name of the table linking person_id to
#'   SNP genotypes. Must contain columns for person identifier, snp_id, and
#'   genotype (coded as 0/1/2 count of effect alleles).
#' @param indexDateOffset Integer. Days offset from cohort start date for
#'   defining the index date. Default is 0.
#' @param washoutPeriod Integer. Minimum days of prior observation required
#'   before index date. Default is 365.
#' @param excludePriorOutcome Logical. If TRUE, persons with the outcome before
#'   their index date are excluded. Default is TRUE.
#' @param genomicPersonIdColumn Character. Name of the person identifier column
#'   in the genomic linkage table if it differs from "person_id". Default is
#'   "person_id".
#'
#' @return A data frame with one row per person and columns:
#'   \describe{
#'     \item{person_id}{Integer person identifier.}
#'     \item{outcome}{Integer 0/1 outcome status.}
#'     \item{age_at_index}{Age in years at index date.}
#'     \item{gender_concept_id}{OMOP concept ID for gender.}
#'     \item{index_date}{Date of cohort entry.}
#'     \item{snp_<sanitized rsID>}{Integer genotype values (0, 1, or 2) for
#'       each instrument SNP, coded as count of the effect allele after
#'       harmonization. Missing genotypes are NA, not 0.}
#'   }
#'
#' @details
#' Genotype coding: genotypes are coded as 0, 1, 2 representing the count of
#' effect alleles. The function performs allele harmonization by comparing the
#' effect allele in the instrument table to the allele coding in the genotype
#' data. If alleles are swapped, genotypes are flipped (2 - genotype).
#'
#' @references
#' Hripcsak, G., et al. (2015). Observational Health Data Sciences and Informatics
#' (OHDSI): Opportunities for Observational Researchers. \emph{Studies in Health
#' Technology and Informatics}, 216, 574-578.
#'
#' @examples
#' \dontrun{
#' connectionDetails <- DatabaseConnector::createConnectionDetails(
#'   dbms = "postgresql",
#'   server = "localhost/ohdsi",
#'   user = "user",
#'   password = "password"
#' )
#' instruments <- getMRInstruments("ieu-a-1119")
#' cohort <- buildMRCohort(
#'   connectionDetails = connectionDetails,
#'   cdmDatabaseSchema = "cdm",
#'   cohortDatabaseSchema = "results",
#'   cohortTable = "cohort",
#'   outcomeCohortId = 1234,
#'   instrumentTable = instruments,
#'   genomicLinkageSchema = "genomics",
#'   genomicLinkageTable = "genotype_data"
#' )
#' }
#'
#' @seealso \code{\link{getMRInstruments}}, \code{\link{harmonizeAlleles}},
#'   \code{\link{buildMRCovariates}}
#'
#' @export
buildMRCohort <- function(connectionDetails,
                          cdmDatabaseSchema,
                          cohortDatabaseSchema,
                          cohortTable,
                          outcomeCohortId,
                          instrumentTable,
                          genomicLinkageSchema,
                          genomicLinkageTable,
                          indexDateOffset = 0,
                          washoutPeriod = 365,
                          excludePriorOutcome = TRUE,
                          genomicPersonIdColumn = "person_id") {
  # Input validation
  checkmate::assertClass(connectionDetails, "connectionDetails")
  checkmate::assertString(cdmDatabaseSchema)
  checkmate::assertString(cohortDatabaseSchema)
  checkmate::assertString(cohortTable)
  checkmate::assertCount(outcomeCohortId, positive = TRUE)
  validateInstrumentTable(instrumentTable)
  checkmate::assertString(genomicLinkageSchema)
  checkmate::assertString(genomicLinkageTable)
  checkmate::assertInt(indexDateOffset)
  checkmate::assertCount(washoutPeriod)
  checkmate::assertLogical(excludePriorOutcome, len = 1)
  checkmate::assertString(genomicPersonIdColumn)

  if (!requireNamespace("DatabaseConnector", quietly = TRUE)) {
    stop("Package 'DatabaseConnector' is required for buildMRCohort(). ",
         "Install it with: remotes::install_github('ohdsi/DatabaseConnector')",
         call. = FALSE)
  }
  if (!requireNamespace("SqlRender", quietly = TRUE)) {
    stop("Package 'SqlRender' is required for buildMRCohort(). ",
         "Install it with: remotes::install_github('ohdsi/SqlRender')",
         call. = FALSE)
  }

  connection <- DatabaseConnector::connect(connectionDetails)
  on.exit(DatabaseConnector::disconnect(connection), add = TRUE)

  # Step 1: Extract outcome cohort
  message("Extracting outcome cohort...")
  sql <- loadRenderTranslateSql(
    sqlFileName = "extractOutcomeCohort.sql",
    dbms = connectionDetails$dbms,
    cdm_database_schema = cdmDatabaseSchema,
    cohort_database_schema = cohortDatabaseSchema,
    cohort_table = cohortTable,
    outcome_cohort_id = outcomeCohortId,
    washout_days = washoutPeriod,
    exclude_prior_outcome = as.integer(excludePriorOutcome)
  )
  DatabaseConnector::executeSql(connection, sql, progressBar = FALSE,
                                reportOverallTime = FALSE)

  cohortData <- DatabaseConnector::querySql(
    connection,
    "SELECT * FROM #mr_cohort;",
    snakeCaseToCamelCase = TRUE
  )

  if (nrow(cohortData) == 0) {
    stop("No persons found in outcome cohort. Check cohort definition and database schema.")
  }

  nCases <- sum(cohortData$outcome == 1)
  nControls <- sum(cohortData$outcome == 0)
  message(sprintf("Cohort extracted: %d cases, %d controls.", nCases, nControls))

  if (nCases < 50) {
    warning(sprintf(
      "Only %d cases in outcome cohort. Results may be unstable. Consider expanding cohort definition or adding sites.",
      nCases
    ))
  }

  # Step 2: Extract genotypes
  message("Extracting genotype data...")
  snpIdList <- paste(paste0("'", instrumentTable$snp_id, "'"), collapse = ", ")
  sql <- loadRenderTranslateSql(
    sqlFileName = "extractGenotypes.sql",
    dbms = connectionDetails$dbms,
    genomic_linkage_schema = genomicLinkageSchema,
    genomic_linkage_table = genomicLinkageTable,
    genomic_person_id_column = genomicPersonIdColumn,
    cohort_database_schema = cohortDatabaseSchema,
    cohort_table = cohortTable,
    outcome_cohort_id = outcomeCohortId,
    snp_ids = snpIdList
  )
  genotypeData <- DatabaseConnector::querySql(
    connection, sql,
    snakeCaseToCamelCase = TRUE
  )

  if (nrow(genotypeData) == 0) {
    stop(sprintf(
      paste0("No persons in cohort have genotype data in %s.%s. ",
             "Verify genomic linkage table schema and person_id join."),
      genomicLinkageSchema, genomicLinkageTable
    ))
  }

  # Clean up temp tables
  DatabaseConnector::executeSql(
    connection,
    "TRUNCATE TABLE #mr_cohort; DROP TABLE #mr_cohort;",
    progressBar = FALSE, reportOverallTime = FALSE
  )

  # Step 3: Reshape genotypes to wide format and merge
  message("Reshaping genotype data...")
  genotypeWide <- reshapeGenotypes(genotypeData, instrumentTable)

  # Step 4: Merge cohort with genotypes
  result <- dplyr::left_join(
    cohortData,
    genotypeWide,
    by = "personId"
  )

  # Report genotype missingness
  snpCols <- grep("^snp_", names(result), value = TRUE)
  for (col in snpCols) {
    nMissing <- sum(is.na(result[[col]]))
    pctMissing <- 100 * nMissing / nrow(result)
    if (pctMissing > 10) {
      warning(sprintf("SNP %s has %.1f%% missing genotypes.", col, pctMissing))
    }
  }

  nGenotyped <- sum(complete.cases(result[, snpCols, drop = FALSE]))
  message(sprintf("Cohort build complete: %d persons, %d fully genotyped.",
                  nrow(result), nGenotyped))

  result
}


#' Reshape Long-Format Genotype Data to Wide Format
#'
#' @description Converts genotype data from long format (person_id, snp_id,
#'   genotype) to wide format with one column per SNP. Missing genotypes
#'   are coded as NA.
#'
#' @param genotypeData Data frame in long format with columns personId, snpId,
#'   genotype.
#' @param instrumentTable Instrument table for SNP ordering.
#'
#' @return Data frame in wide format with personId and one column per SNP.
#'
#' @keywords internal
reshapeGenotypes <- function(genotypeData, instrumentTable) {
  # Validate genotype values
  validGenotypes <- c(0L, 1L, 2L, NA_integer_)
  genotypeData$genotype <- as.integer(genotypeData$genotype)
  invalidGeno <- !genotypeData$genotype %in% validGenotypes & !is.na(genotypeData$genotype)
  if (any(invalidGeno)) {
    warning(sprintf("Found %d invalid genotype values (not 0, 1, or 2). Setting to NA.",
                    sum(invalidGeno)))
    genotypeData$genotype[invalidGeno] <- NA_integer_
  }

  # Pivot to wide format
  personIds <- unique(genotypeData$personId)
  snpIds <- instrumentTable$snp_id

  wideData <- data.frame(personId = personIds)
  for (snp in snpIds) {
    colName <- makeSnpColumnName(snp)
    snpSubset <- genotypeData[genotypeData$snpId == snp, ]
    snpMerge <- data.frame(
      personId = snpSubset$personId,
      geno = snpSubset$genotype,
      stringsAsFactors = FALSE
    )
    wideData <- dplyr::left_join(wideData, snpMerge, by = "personId")
    names(wideData)[names(wideData) == "geno"] <- colName
  }

  wideData
}
