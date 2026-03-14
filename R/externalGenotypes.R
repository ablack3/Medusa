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

#' Format an external dosage table for use with Medusa
#'
#' @title Standardize PLINK/VCF/Hail dosage exports to Medusa SNP columns
#'
#' @description
#' Converts externally extracted dosage data into the internal format used by
#' Medusa: one row per person with SNP columns named \code{snp_<rsid>}. This is
#' intended for datasets such as All of Us, UK Biobank, or other genomics
#' resources where dosages are available outside the OMOP Genomic CDM
#' \code{VARIANT_OCCURRENCE} table.
#'
#' The input can be either:
#' \itemize{
#'   \item wide format: a person identifier column plus one column per SNP, or
#'   \item long format: person identifier, SNP identifier, and dosage columns.
#' }
#'
#' When allele columns are provided for long-format input, the instrument table
#' is harmonized so downstream MR uses the same allele coding as the external
#' dosage extract.
#'
#' @param dosageTable Data frame containing external dosage data.
#' @param instrumentTable Instrument table produced by
#'   \code{\link{getMRInstruments}} or \code{\link{createInstrumentTable}}.
#' @param inputFormat Character. One of \code{"auto"}, \code{"wide"}, or
#'   \code{"long"}. Default is \code{"auto"}.
#' @param personIdColumn Optional character name of the person identifier
#'   column. When NULL, Medusa searches for \code{personId} or
#'   \code{person_id}.
#' @param snpIdColumn Optional character name of the SNP identifier column for
#'   long-format input. When NULL, Medusa searches for \code{snpId} or
#'   \code{snp_id}.
#' @param dosageColumn Optional character name of the dosage column for
#'   long-format input. When NULL, Medusa searches for \code{dosage},
#'   \code{genotype}, or \code{genotypeRaw}.
#' @param codedAlleleColumn Optional character column naming the coded allele in
#'   long-format input.
#' @param noncodedAlleleColumn Optional character column naming the non-coded
#'   allele in long-format input.
#'
#' @return A list with class \code{"medusaExternalDosage"} containing:
#'   \describe{
#'     \item{dosageTable}{Wide data frame with \code{personId} and SNP columns.}
#'     \item{instrumentTable}{Instrument table aligned to the dosage coding.}
#'     \item{metadata}{List describing the detected format and harmonization
#'       steps.}
#'   }
#'
#' @examples
#' instruments <- createInstrumentTable(
#'   snpId = c("rs1", "rs2"),
#'   effectAllele = c("A", "G"),
#'   otherAllele = c("C", "T"),
#'   betaZX = c(0.2, 0.1),
#'   seZX = c(0.05, 0.04),
#'   pvalZX = c(1e-8, 2e-7),
#'   eaf = c(0.3, 0.4)
#' )
#' dosage <- data.frame(
#'   person_id = c(1, 2),
#'   rs1 = c(0, 1),
#'   rs2 = c(2, 1),
#'   stringsAsFactors = FALSE
#' )
#' formatted <- formatExternalDosageTable(dosage, instruments)
#'
#' @export
formatExternalDosageTable <- function(dosageTable,
                                      instrumentTable,
                                      inputFormat = c("auto", "wide", "long"),
                                      personIdColumn = NULL,
                                      snpIdColumn = NULL,
                                      dosageColumn = NULL,
                                      codedAlleleColumn = NULL,
                                      noncodedAlleleColumn = NULL) {
  checkmate::assertDataFrame(dosageTable, min.rows = 1)
  validateInstrumentTable(instrumentTable)
  inputFormat <- match.arg(inputFormat)

  personIdColumn <- personIdColumn %||% detectExternalColumn(
    names(dosageTable),
    c("personId", "person_id")
  )

  if (is.null(personIdColumn)) {
    stop("Could not identify a person ID column in dosageTable.")
  }

  resolvedFormat <- detectExternalDosageFormat(
    dosageTable = dosageTable,
    instrumentTable = instrumentTable,
    inputFormat = inputFormat,
    personIdColumn = personIdColumn,
    snpIdColumn = snpIdColumn,
    dosageColumn = dosageColumn
  )

  formatted <- switch(
    resolvedFormat,
    wide = formatWideDosageTable(
      dosageTable = dosageTable,
      instrumentTable = instrumentTable,
      personIdColumn = personIdColumn
    ),
    long = formatLongDosageTable(
      dosageTable = dosageTable,
      instrumentTable = instrumentTable,
      personIdColumn = personIdColumn,
      snpIdColumn = snpIdColumn,
      dosageColumn = dosageColumn,
      codedAlleleColumn = codedAlleleColumn,
      noncodedAlleleColumn = noncodedAlleleColumn
    )
  )

  structure(
    list(
      dosageTable = formatted$dosageTable,
      instrumentTable = formatted$instrumentTable,
      metadata = c(
        list(
          inputFormat = resolvedFormat,
          personIdColumn = personIdColumn
        ),
        formatted$metadata
      )
    ),
    class = "medusaExternalDosage"
  )
}


#' @keywords internal
detectExternalColumn <- function(columnNames, candidates) {
  match <- intersect(candidates, columnNames)
  if (length(match) == 0) {
    return(NULL)
  }
  match[[1]]
}


#' @keywords internal
detectExternalDosageFormat <- function(dosageTable,
                                       instrumentTable,
                                       inputFormat,
                                       personIdColumn,
                                       snpIdColumn,
                                       dosageColumn) {
  if (!identical(inputFormat, "auto")) {
    return(inputFormat)
  }

  expectedWide <- unique(c(makeSnpColumnName(instrumentTable$snp_id), instrumentTable$snp_id))
  nonIdColumns <- setdiff(names(dosageTable), personIdColumn)
  if (length(intersect(nonIdColumns, expectedWide)) > 0) {
    return("wide")
  }

  snpIdColumn <- snpIdColumn %||% detectExternalColumn(
    names(dosageTable),
    c("snpId", "snp_id")
  )
  dosageColumn <- dosageColumn %||% detectExternalColumn(
    names(dosageTable),
    c("dosage", "genotype", "genotypeRaw")
  )
  if (!is.null(snpIdColumn) && !is.null(dosageColumn)) {
    return("long")
  }

  stop(
    paste0(
      "Could not infer dosageTable format. Supply inputFormat = 'wide' or 'long', ",
      "or provide standard column names."
    )
  )
}


#' @keywords internal
formatWideDosageTable <- function(dosageTable,
                                  instrumentTable,
                                  personIdColumn) {
  result <- dosageTable
  names(result)[names(result) == personIdColumn] <- "personId"

  expectedCols <- makeSnpColumnName(instrumentTable$snp_id)
  rawCols <- setNames(instrumentTable$snp_id, expectedCols)

  for (i in seq_along(expectedCols)) {
    if (!expectedCols[[i]] %in% names(result) && rawCols[[i]] %in% names(result)) {
      names(result)[names(result) == rawCols[[i]]] <- expectedCols[[i]]
    }
  }

  availableCols <- intersect(expectedCols, names(result))
  if (length(availableCols) == 0) {
    stop("No SNP columns in wide dosageTable match instrumentTable SNP IDs.")
  }

  standardized <- result[, c("personId", availableCols), drop = FALSE]
  standardized <- coerceWideDosageNumeric(standardized, availableCols)

  list(
    dosageTable = standardized,
    instrumentTable = instrumentTable[instrumentTable$snp_id %in%
                                        sub("^snp_", "", availableCols), ,
                                      drop = FALSE],
    metadata = list(
      harmonized = FALSE,
      availableSnps = sub("^snp_", "", availableCols)
    )
  )
}


#' @keywords internal
formatLongDosageTable <- function(dosageTable,
                                  instrumentTable,
                                  personIdColumn,
                                  snpIdColumn,
                                  dosageColumn,
                                  codedAlleleColumn,
                                  noncodedAlleleColumn) {
  snpIdColumn <- snpIdColumn %||% detectExternalColumn(
    names(dosageTable),
    c("snpId", "snp_id")
  )
  dosageColumn <- dosageColumn %||% detectExternalColumn(
    names(dosageTable),
    c("dosage", "genotype", "genotypeRaw")
  )

  if (is.null(snpIdColumn) || is.null(dosageColumn)) {
    stop("Long-format dosageTable must include SNP ID and dosage columns.")
  }

  genotypeData <- data.frame(
    personId = dosageTable[[personIdColumn]],
    snpId = dosageTable[[snpIdColumn]],
    genotype = dosageTable[[dosageColumn]],
    stringsAsFactors = FALSE
  )
  genotypeData <- genotypeData[genotypeData$snpId %in% instrumentTable$snp_id, , drop = FALSE]
  if (nrow(genotypeData) == 0) {
    stop("Long-format dosageTable does not contain any instrument SNPs.")
  }

  genotypeData$genotype <- normalizeExternalDosage(genotypeData$genotype)

  harmonizedInstruments <- instrumentTable
  harmonized <- FALSE
  if (!is.null(codedAlleleColumn) && !is.null(noncodedAlleleColumn)) {
    alleleInfo <- unique(data.frame(
      snp_id = dosageTable[[snpIdColumn]],
      allele_coded = dosageTable[[codedAlleleColumn]],
      allele_noncoded = dosageTable[[noncodedAlleleColumn]],
      stringsAsFactors = FALSE
    ))
    alleleInfo <- alleleInfo[alleleInfo$snp_id %in% instrumentTable$snp_id, , drop = FALSE]
    harmonizedResult <- harmonizeAlleles(harmonizedInstruments, alleleInfo)
    harmonizedInstruments <- harmonizedResult$instrumentTable
    genotypeData <- genotypeData[genotypeData$snpId %in% harmonizedInstruments$snp_id, ,
                                 drop = FALSE]
    harmonized <- TRUE
  }

  wideData <- reshapeGenotypes(genotypeData, harmonizedInstruments)

  list(
    dosageTable = wideData,
    instrumentTable = harmonizedInstruments,
    metadata = list(
      harmonized = harmonized,
      availableSnps = harmonizedInstruments$snp_id
    )
  )
}


#' @keywords internal
normalizeExternalDosage <- function(values) {
  if (is.character(values)) {
    values <- convertGenotypeString(values)
  }
  values <- suppressWarnings(as.numeric(values))
  invalid <- !is.na(values) & (values < 0 | values > 2)
  if (any(invalid)) {
    warning(
      sprintf(
        "Found %d external dosage values outside [0, 2]; coercing them to NA.",
        sum(invalid)
      )
    )
    values[invalid] <- NA_real_
  }
  as.integer(round(values))
}


#' @keywords internal
coerceWideDosageNumeric <- function(dosageTable, snpColumns) {
  for (col in snpColumns) {
    dosageTable[[col]] <- normalizeExternalDosage(dosageTable[[col]])
  }
  dosageTable
}
