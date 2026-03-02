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

#' Retrieve and Clump Genetic Instruments for Mendelian Randomization
#'
#' @title Instrument assembly from OpenGWAS
#'
#' @description
#' Queries the IEU OpenGWAS database via the \code{ieugwasr} package to retrieve
#' GWAS summary statistics for a specified exposure trait, applies LD clumping
#' to obtain independent instruments, computes approximate F-statistics, flags
#' strand-ambiguous SNPs, and returns a clean instrument table ready for
#' distribution to federated analysis sites.
#'
#' This function runs at the coordinator node. The returned instrument table
#' is serialized to disk and distributed to all sites unchanged.
#'
#' @param exposureTraitId Character string. IEU OpenGWAS trait ID
#'   (e.g., "ieu-a-1119" for IL-6 receptor levels).
#' @param pThreshold Numeric. Genome-wide significance p-value threshold for
#'   selecting SNPs. Default is 5e-8.
#' @param r2Threshold Numeric. LD clumping r-squared threshold. SNP pairs with
#'   r-squared above this value will be pruned, keeping the more significant SNP.
#'   Default is 0.001.
#' @param kb Numeric. LD clumping window in kilobases. Default is 10000.
#' @param ancestryPopulation Character. Reference panel ancestry population for
#'   LD clumping. One of "EUR", "EAS", "AFR", "SAS", "AMR". Default is "EUR".
#' @param additionalSnps Optional character vector of SNP rsIDs to force-include
#'   in the instrument set (added after clumping, not subject to LD pruning).
#'   Default is NULL.
#'
#' @return A data frame with one row per independent instrument SNP and columns:
#'   \describe{
#'     \item{snp_id}{rsID of the SNP (e.g., "rs2228145").}
#'     \item{effect_allele}{The allele associated with increased exposure level.}
#'     \item{other_allele}{The non-effect allele.}
#'     \item{beta_ZX}{Effect size of the SNP on the exposure (log scale).}
#'     \item{se_ZX}{Standard error of beta_ZX.}
#'     \item{pval_ZX}{P-value for the SNP-exposure association.}
#'     \item{eaf}{Effect allele frequency in the GWAS reference population.}
#'     \item{gene_region}{Nearest gene or genomic region annotation.}
#'     \item{fStatistic}{Approximate F-statistic: (beta_ZX / se_ZX)^2.}
#'     \item{strandAmbiguous}{Logical. TRUE if the SNP is strand-ambiguous
#'       (A/T or G/C allele pair).}
#'   }
#'   The data frame also carries the following attributes:
#'   \describe{
#'     \item{retrievalTimestamp}{POSIXct timestamp of when instruments were retrieved.}
#'     \item{exposureTraitId}{The trait ID used for retrieval.}
#'     \item{parameters}{List of all parameter values used.}
#'   }
#'
#' @details
#' The function queries the IEU OpenGWAS API for all SNP associations with the
#' specified trait below \code{pThreshold}, then applies LD clumping using the
#' specified reference panel to retain only independent instruments. The
#' approximate F-statistic is computed as (beta / SE)^2 for each SNP. SNPs
#' with F < 10 are flagged as potentially weak instruments via a warning message
#' but are not removed automatically.
#'
#' If the OpenGWAS API is unavailable, the function throws an informative error
#' suggesting the user provide a cached instrument table instead.
#'
#' @references
#' Hemani, G., et al. (2018). The MR-Base platform supports systematic causal
#' inference across the human phenome. \emph{eLife}, 7, e34408.
#'
#' @examples
#' # Using simulated data (no API call)
#' instruments <- simulateInstrumentTable(nSnps = 10)
#' head(instruments)
#'
#' \dontrun{
#' # Real API call (requires internet)
#' instruments <- getMRInstruments(
#'   exposureTraitId = "ieu-a-1119",
#'   pThreshold = 5e-8,
#'   r2Threshold = 0.001,
#'   kb = 10000,
#'   ancestryPopulation = "EUR"
#' )
#' }
#'
#' @seealso \code{\link{harmonizeAlleles}}, \code{\link{buildMRCohort}},
#'   \code{\link{computeApproxFStatistic}}
#'
#' @export
getMRInstruments <- function(exposureTraitId,
                             pThreshold = 5e-8,
                             r2Threshold = 0.001,
                             kb = 10000,
                             ancestryPopulation = "EUR",
                             additionalSnps = NULL) {
  # Input validation
  checkmate::assertString(exposureTraitId)
  checkmate::assertNumber(pThreshold, lower = 0, upper = 1)
  checkmate::assertNumber(r2Threshold, lower = 0, upper = 1)
  checkmate::assertNumber(kb, lower = 1)
  checkmate::assertChoice(ancestryPopulation, c("EUR", "EAS", "AFR", "SAS", "AMR"))
  if (!is.null(additionalSnps)) {
    checkmate::assertCharacter(additionalSnps, min.len = 1)
  }

  if (!requireNamespace("ieugwasr", quietly = TRUE)) {
    stop("Package 'ieugwasr' is required for getMRInstruments(). ",
         "Install it with: remotes::install_github('mrcieu/ieugwasr')",
         call. = FALSE)
  }

  message(sprintf("Querying OpenGWAS for trait '%s' with p < %g...",
                  exposureTraitId, pThreshold))

  # Query GWAS associations
  associations <- tryCatch(
    {
      ieugwasr::associations(
        variants = NULL,
        id = exposureTraitId,
        proxies = 0
      )
    },
    error = function(e) {
      stop(sprintf(
        paste0("Failed to query OpenGWAS for trait '%s'. ",
               "Error: %s\n",
               "Check your internet connection or provide a cached instrument ",
               "table via the instrumentTable parameter in downstream functions."),
        exposureTraitId, conditionMessage(e)
      ))
    }
  )

  if (is.null(associations) || nrow(associations) == 0) {
    stop(sprintf(
      paste0("No genome-wide significant SNPs found for trait '%s'. ",
             "Consider relaxing pThreshold or checking trait ID."),
      exposureTraitId
    ))
  }

  # Filter by p-value threshold
  associations <- associations[associations$p < pThreshold, , drop = FALSE]

  if (nrow(associations) == 0) {
    stop(sprintf(
      paste0("No genome-wide significant SNPs found for trait '%s' at p < %g. ",
             "Consider relaxing pThreshold or checking trait ID."),
      exposureTraitId, pThreshold
    ))
  }

  message(sprintf("Found %d SNPs at p < %g. Applying LD clumping...",
                  nrow(associations), pThreshold))

  # LD clumping
  clumped <- tryCatch(
    {
      ieugwasr::ld_clump(
        dat = data.frame(
          rsid = associations$rsid,
          pval = associations$p,
          id = associations$id,
          stringsAsFactors = FALSE
        ),
        clump_r2 = r2Threshold,
        clump_kb = kb,
        pop = ancestryPopulation
      )
    },
    error = function(e) {
      stop(sprintf(
        paste0("LD clumping failed for trait '%s'. Error: %s\n",
               "Check your internet connection or the OpenGWAS LD reference ",
               "panel availability."),
        exposureTraitId, conditionMessage(e)
      ))
    }
  )

  if (is.null(clumped) || nrow(clumped) == 0) {
    stop(sprintf(
      paste0("No SNPs remained after LD clumping for trait '%s'. ",
             "This is unexpected -- check input data."),
      exposureTraitId
    ))
  }

  # Filter associations to clumped SNPs
  associations <- associations[associations$rsid %in% clumped$rsid, , drop = FALSE]

  message(sprintf("%d independent instruments after LD clumping (r2 < %g, kb = %d).",
                  nrow(associations), r2Threshold, kb))

  # Build clean instrument table
  instrumentTable <- data.frame(
    snp_id = associations$rsid,
    effect_allele = toupper(associations$ea),
    other_allele = toupper(associations$nea),
    beta_ZX = associations$beta,
    se_ZX = associations$se,
    pval_ZX = associations$p,
    eaf = associations$eaf,
    gene_region = if ("gene" %in% names(associations)) associations$gene else NA_character_,
    stringsAsFactors = FALSE
  )

  # Add forced-include SNPs
  if (!is.null(additionalSnps)) {
    newSnps <- setdiff(additionalSnps, instrumentTable$snp_id)
    if (length(newSnps) > 0) {
      message(sprintf("Force-including %d additional SNPs: %s",
                      length(newSnps), paste(newSnps, collapse = ", ")))
      additionalAssoc <- tryCatch(
        {
          ieugwasr::associations(
            variants = newSnps,
            id = exposureTraitId,
            proxies = 0
          )
        },
        error = function(e) {
          warning(sprintf("Failed to retrieve additional SNPs: %s",
                          conditionMessage(e)))
          NULL
        }
      )
      if (!is.null(additionalAssoc) && nrow(additionalAssoc) > 0) {
        additionalRows <- data.frame(
          snp_id = additionalAssoc$rsid,
          effect_allele = toupper(additionalAssoc$ea),
          other_allele = toupper(additionalAssoc$nea),
          beta_ZX = additionalAssoc$beta,
          se_ZX = additionalAssoc$se,
          pval_ZX = additionalAssoc$p,
          eaf = additionalAssoc$eaf,
          gene_region = if ("gene" %in% names(additionalAssoc)) {
            additionalAssoc$gene
          } else {
            NA_character_
          },
          stringsAsFactors = FALSE
        )
        instrumentTable <- rbind(instrumentTable, additionalRows)
      }
    }
  }

  # Compute F-statistics
  instrumentTable$fStatistic <- computeApproxFStatistic(
    instrumentTable$beta_ZX,
    instrumentTable$se_ZX
  )

  # Flag strand-ambiguous SNPs
  instrumentTable$strandAmbiguous <- isStrandAmbiguous(
    instrumentTable$effect_allele,
    instrumentTable$other_allele
  )

  # Warnings
  nAmbiguous <- sum(instrumentTable$strandAmbiguous)
  if (nAmbiguous > 0) {
    warning(sprintf(
      paste0("%d strand-ambiguous SNP(s) detected (A/T or G/C allele pairs): %s. ",
             "These may cause harmonization issues at sites with different strand ",
             "conventions. Consider removing them or using allele frequency to ",
             "infer strand."),
      nAmbiguous,
      paste(instrumentTable$snp_id[instrumentTable$strandAmbiguous], collapse = ", ")
    ))
  }

  nWeak <- sum(instrumentTable$fStatistic < 10)
  if (nWeak > 0) {
    warning(sprintf(
      paste0("%d instrument(s) have F-statistic < 10 (potentially weak): %s. ",
             "Weak instruments may bias MR estimates toward the null."),
      nWeak,
      paste(instrumentTable$snp_id[instrumentTable$fStatistic < 10], collapse = ", ")
    ))
  }

  if (nrow(instrumentTable) < 3) {
    warning(sprintf(
      paste0("Only %d instruments available. Sensitivity analyses requiring ",
             ">= 3 SNPs will be skipped."),
      nrow(instrumentTable)
    ))
  }

  # Add metadata attributes
  attr(instrumentTable, "retrievalTimestamp") <- Sys.time()
  attr(instrumentTable, "exposureTraitId") <- exposureTraitId
  attr(instrumentTable, "parameters") <- list(
    pThreshold = pThreshold,
    r2Threshold = r2Threshold,
    kb = kb,
    ancestryPopulation = ancestryPopulation,
    additionalSnps = additionalSnps
  )

  message(sprintf("Instrument assembly complete: %d SNPs (%d strand-ambiguous, %d weak).",
                  nrow(instrumentTable), nAmbiguous, nWeak))

  instrumentTable
}


#' Create Instrument Table from Local Data
#'
#' @title Build instrument table from user-provided data
#'
#' @description
#' Creates a properly formatted instrument table from user-provided GWAS summary
#' statistics, bypassing the OpenGWAS API. Useful when working offline or with
#' custom GWAS results not available in OpenGWAS.
#'
#' @param snpId Character vector of SNP rsIDs.
#' @param effectAllele Character vector of effect alleles.
#' @param otherAllele Character vector of other alleles.
#' @param betaZX Numeric vector of SNP-exposure effect estimates.
#' @param seZX Numeric vector of standard errors.
#' @param pvalZX Numeric vector of p-values.
#' @param eaf Numeric vector of effect allele frequencies.
#' @param geneRegion Optional character vector of gene/region annotations.
#'
#' @return A data frame with the same structure as \code{\link{getMRInstruments}}
#'   output.
#'
#' @examples
#' instruments <- createInstrumentTable(
#'   snpId = c("rs1234", "rs5678"),
#'   effectAllele = c("A", "G"),
#'   otherAllele = c("G", "C"),
#'   betaZX = c(0.5, 0.3),
#'   seZX = c(0.05, 0.08),
#'   pvalZX = c(1e-10, 1e-5),
#'   eaf = c(0.3, 0.45)
#' )
#'
#' @seealso \code{\link{getMRInstruments}}
#'
#' @export
createInstrumentTable <- function(snpId,
                                  effectAllele,
                                  otherAllele,
                                  betaZX,
                                  seZX,
                                  pvalZX,
                                  eaf,
                                  geneRegion = NULL) {
  n <- length(snpId)
  checkmate::assertCharacter(snpId, min.len = 1, any.missing = FALSE)
  checkmate::assertCharacter(effectAllele, len = n, any.missing = FALSE)
  checkmate::assertCharacter(otherAllele, len = n, any.missing = FALSE)
  checkmate::assertNumeric(betaZX, len = n, any.missing = FALSE)
  checkmate::assertNumeric(seZX, len = n, lower = 0, any.missing = FALSE)
  checkmate::assertNumeric(pvalZX, len = n, lower = 0, upper = 1, any.missing = FALSE)
  checkmate::assertNumeric(eaf, len = n, lower = 0, upper = 1, any.missing = FALSE)
  if (!is.null(geneRegion)) {
    checkmate::assertCharacter(geneRegion, len = n)
  }

  instrumentTable <- data.frame(
    snp_id = snpId,
    effect_allele = toupper(effectAllele),
    other_allele = toupper(otherAllele),
    beta_ZX = betaZX,
    se_ZX = seZX,
    pval_ZX = pvalZX,
    eaf = eaf,
    gene_region = geneRegion %||% rep(NA_character_, n),
    stringsAsFactors = FALSE
  )

  instrumentTable$fStatistic <- computeApproxFStatistic(betaZX, seZX)
  instrumentTable$strandAmbiguous <- isStrandAmbiguous(effectAllele, otherAllele)

  attr(instrumentTable, "retrievalTimestamp") <- Sys.time()
  attr(instrumentTable, "exposureTraitId") <- "user_provided"
  attr(instrumentTable, "parameters") <- list(source = "createInstrumentTable")

  instrumentTable
}
