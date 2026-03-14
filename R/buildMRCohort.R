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
#' Genotype data is extracted from the \strong{VARIANT_OCCURRENCE} table defined
#' by the OMOP CDM Genomic Extension. The minimal required columns from that
#' table are:
#' \itemize{
#'   \item \code{person_id} — links variants to persons
#'   \item \code{rs_id} — dbSNP rs identifier for the variant
#'   \item \code{genotype} — genotype call (VCF-style "0/0", "0/1", "1/1" or
#'     plain integer "0", "1", "2")
#' }
#' Additionally, \code{reference_allele} and \code{alternate_allele} are used
#' for allele harmonization when available.
#'
#' The function also queries: PERSON (age, sex), CONDITION_OCCURRENCE (outcome),
#' and OBSERVATION_PERIOD (eligibility). All SQL is rendered and translated via
#' SqlRender for cross-dialect compatibility.
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
#' @param genomicDatabaseSchema Character. Schema containing the VARIANT_OCCURRENCE
#'   table from the OMOP Genomic Extension. Defaults to \code{cdmDatabaseSchema}.
#' @param indexDateOffset Integer. Days offset from cohort start date for
#'   defining the index date. Default is 0.
#' @param washoutPeriod Integer. Minimum days of prior observation required
#'   before index date. Default is 365.
#' @param excludePriorOutcome Logical. If TRUE, persons with the outcome before
#'   their index date are excluded. Default is TRUE.
#' @param negativeControlCohortIds Optional integer vector of cohort definition
#'   IDs for negative control outcomes. When provided, the function extracts
#'   negative control outcome flags and adds \code{nc_outcome_<id>} columns
#'   (binary 0/1) to the returned data frame. These columns are used by
#'   \code{\link{runNegativeControlAnalysis}} for empirical calibration.
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
#' data (alternate allele from VARIANT_OCCURRENCE). If alleles are swapped,
#' genotypes are flipped (2 - genotype) and instrument beta_ZX is negated.
#'
#' @references
#' OHDSI Genomic CDM: \url{https://github.com/OHDSI/Genomic-CDM}
#'
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
#'   genomicDatabaseSchema = "genomics"
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
                          genomicDatabaseSchema = cdmDatabaseSchema,
                          indexDateOffset = 0,
                          washoutPeriod = 365,
                          excludePriorOutcome = TRUE,
                          negativeControlCohortIds = NULL) {
  # Input validation
  checkmate::assertClass(connectionDetails, "connectionDetails")
  checkmate::assertString(cdmDatabaseSchema)
  checkmate::assertString(cohortDatabaseSchema)
  checkmate::assertString(cohortTable)
  checkmate::assertCount(outcomeCohortId, positive = TRUE)
  validateInstrumentTable(instrumentTable)
  checkmate::assertString(genomicDatabaseSchema)
  checkmate::assertInt(indexDateOffset)
  checkmate::assertCount(washoutPeriod)
  checkmate::assertLogical(excludePriorOutcome, len = 1)

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

  # Step 2: Extract genotypes from VARIANT_OCCURRENCE (OMOP Genomic Extension)
  message("Extracting genotype data from VARIANT_OCCURRENCE...")
  invalidSnpIds <- instrumentTable$snp_id[!grepl("^rs[0-9]+$", instrumentTable$snp_id)]
  if (length(invalidSnpIds) > 0) {
    stop(sprintf(
      "Invalid SNP IDs detected (expected rs<number> format): %s",
      paste(utils::head(invalidSnpIds, 5), collapse = ", ")
    ))
  }
  snpIdList <- paste(paste0("'", instrumentTable$snp_id, "'"), collapse = ", ")
  sql <- loadRenderTranslateSql(
    sqlFileName = "extractGenotypes.sql",
    dbms = connectionDetails$dbms,
    genomic_schema = genomicDatabaseSchema,
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
      paste0("No persons in cohort have genotype data in %s.VARIANT_OCCURRENCE. ",
             "Verify the OMOP Genomic Extension tables exist and contain data ",
             "for the requested rs IDs."),
      genomicDatabaseSchema
    ))
  }

  # Step 2b: Extract negative control outcomes (optional)
  ncData <- NULL
  if (!is.null(negativeControlCohortIds)) {
    checkmate::assertIntegerish(negativeControlCohortIds, lower = 1,
                                 any.missing = FALSE, min.len = 1)
    message(sprintf("Extracting %d negative control outcomes...",
                    length(negativeControlCohortIds)))
    ncSql <- loadRenderTranslateSql(
      sqlFileName = "extractNegativeControls.sql",
      dbms = connectionDetails$dbms,
      cohort_database_schema = cohortDatabaseSchema,
      cohort_table = cohortTable,
      negative_control_cohort_ids = paste(negativeControlCohortIds, collapse = ", "),
      use_index_date = 0
    )
    ncData <- DatabaseConnector::querySql(
      connection, ncSql, snakeCaseToCamelCase = TRUE
    )
  }

  # Clean up temp tables
  DatabaseConnector::executeSql(
    connection,
    "TRUNCATE TABLE #mr_cohort; DROP TABLE #mr_cohort;",
    progressBar = FALSE, reportOverallTime = FALSE
  )

  # Step 3: Convert genotype strings to integer dosage
  message("Converting genotype values...")
  genotypeData$genotype <- convertGenotypeString(genotypeData$genotypeRaw)

  # Step 4: Build allele table from genotype data for harmonization
  alleleInfo <- unique(genotypeData[, c("snpId", "referenceAllele", "alternateAllele")])
  # In VARIANT_OCCURRENCE, genotype counts the alternate allele
  genotypeAlleles <- data.frame(
    snp_id = alleleInfo$snpId,
    allele_coded = alleleInfo$alternateAllele,
    allele_noncoded = alleleInfo$referenceAllele,
    stringsAsFactors = FALSE
  )

  # Step 5: Harmonize alleles
  message("Harmonizing alleles...")
  cohortAlleleFreqs <- computeCohortAlleleFrequencies(genotypeData)
  harmonized <- harmonizeAlleles(
    instrumentTable,
    genotypeAlleles,
    cohortAlleleFrequencies = cohortAlleleFreqs
  )
  instrumentTable <- harmonized$instrumentTable

  if (nrow(instrumentTable) == 0) {
    stop("No instruments remain after allele harmonization. ",
         "Check that VARIANT_OCCURRENCE alleles match the instrument table.")
  }

  # Step 6: Reshape genotypes to wide format and merge
  message("Reshaping genotype data...")
  genotypeWide <- reshapeGenotypes(genotypeData, instrumentTable)

  # Step 7: Merge cohort with genotypes
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

  # Step 8: Merge negative control outcomes (optional)
  if (!is.null(ncData) && nrow(ncData) > 0) {
    ncOutcomeIds <- unique(ncData$outcomeCohortId)
    for (ncId in ncOutcomeIds) {
      ncSubset <- ncData[ncData$outcomeCohortId == ncId, , drop = FALSE]
      colName <- paste0("nc_outcome_", ncId)
      result[[colName]] <- 0L
      matchIdx <- match(result$personId, ncSubset$personId)
      hasOutcome <- !is.na(matchIdx) &
        ncSubset$hasOutcome[matchIdx] == 1
      hasOutcome[is.na(hasOutcome)] <- FALSE
      result[[colName]][hasOutcome] <- 1L
    }
    message(sprintf("  Added %d negative control outcome columns.",
                    length(ncOutcomeIds)))
  }

  message(sprintf("Cohort build complete: %d persons, %d fully genotyped.",
                  nrow(result), nGenotyped))

  result
}


#' Compute Cohort Allele Frequencies from Genotype Data
#'
#' @description Computes the coded (alternate) allele frequency for each SNP
#'   from the extracted genotype data. Used for palindromic SNP resolution
#'   during allele harmonization.
#'
#' @param genotypeData Data frame with columns snpId and genotype (integer).
#'
#' @return Named numeric vector of allele frequencies, keyed by SNP ID.
#'
#' @keywords internal
computeCohortAlleleFrequencies <- function(genotypeData) {
  snpIds <- unique(genotypeData$snpId)
  freqs <- vapply(snpIds, function(snp) {
    geno <- genotypeData$genotype[genotypeData$snpId == snp]
    validGeno <- geno[!is.na(geno)]
    if (length(validGeno) == 0) return(NA_real_)
    mean(validGeno) / 2
  }, numeric(1))
  names(freqs) <- snpIds
  freqs
}


#' Reshape Long-Format Genotype Data to Wide Format
#'
#' @description Converts genotype data from long format (personId, snpId,
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
