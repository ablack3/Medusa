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

#' Harmonize Alleles Between Instrument Table and Genotype Data
#'
#' @title Allele harmonization for Mendelian Randomization
#'
#' @description
#' Ensures that the effect allele coding in the instrument table (from GWAS)
#' matches the allele coding in the genotype data at each site. Handles allele
#' flipping (when alleles are swapped) and detects strand-ambiguous SNPs (A/T
#' and G/C pairs) that cannot be reliably harmonized.
#'
#' When alleles need to be flipped, beta_ZX is multiplied by -1 and EAF is
#' replaced by 1 - EAF. Palindromic SNPs with EAF near 0.5 (within
#' \code{eafPalindromicThreshold} of 0.5) are dropped because strand cannot
#' be reliably inferred from allele frequency.
#'
#' @param instrumentTable Data frame with columns: snp_id, effect_allele,
#'   other_allele, beta_ZX, se_ZX, eaf. Typically the output of
#'   \code{\link{getMRInstruments}}.
#' @param genotypeAlleles Data frame with columns: snp_id, allele_coded,
#'   allele_noncoded. Describes the allele coding used in the genotype data
#'   at a particular site.
#' @param eafPalindromicThreshold Numeric threshold for removing palindromic
#'   SNPs. Palindromic SNPs with EAF between (0.5 - threshold) and
#'   (0.5 + threshold) are removed. Default is 0.08 (i.e., EAF between
#'   0.42 and 0.58).
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{instrumentTable}{The harmonized instrument table with updated
#'       effect_allele, other_allele, beta_ZX, and eaf where allele flips
#'       were applied. A logical column \code{flipped} is added.}
#'     \item{removedSnps}{Data frame of SNPs that were removed, with a column
#'       \code{reason} indicating why (e.g., "palindromic_ambiguous_eaf",
#'       "allele_mismatch").}
#'   }
#'
#' @details
#' The harmonization algorithm proceeds as follows for each SNP:
#' \enumerate{
#'   \item If the instrument effect_allele matches the genotype allele_coded:
#'     no action needed.
#'   \item If the instrument effect_allele matches the genotype allele_noncoded:
#'     flip the coding (beta_ZX *= -1, eaf = 1 - eaf).
#'   \item If neither direct match is found, try complement alleles
#'     (A<->T, G<->C):
#'     \itemize{
#'       \item If the SNP is palindromic (A/T or G/C) and EAF is near 0.5,
#'         remove the SNP.
#'       \item If palindromic but EAF clearly differs from 0.5, use EAF to
#'         infer strand orientation.
#'     }
#'   \item If no match is found even with complements: remove the SNP with
#'     reason "allele_mismatch".
#' }
#'
#' @references
#' Hartwig, F. P., Davies, N. M., Hemani, G., & Davey Smith, G. (2016).
#' Two-sample Mendelian randomization: avoiding the downsides of a powerful,
#' widely applicable but potentially fallible technique. \emph{International
#' Journal of Epidemiology}, 45(6), 1717-1726.
#'
#' @examples
#' instruments <- data.frame(
#'   snp_id = c("rs1", "rs2", "rs3"),
#'   effect_allele = c("A", "G", "A"),
#'   other_allele = c("G", "T", "C"),
#'   beta_ZX = c(0.5, -0.3, 0.2),
#'   se_ZX = c(0.05, 0.08, 0.06),
#'   pval_ZX = c(1e-10, 1e-5, 1e-8),
#'   eaf = c(0.3, 0.45, 0.6),
#'   stringsAsFactors = FALSE
#' )
#' genotypeAlleles <- data.frame(
#'   snp_id = c("rs1", "rs2", "rs3"),
#'   allele_coded = c("G", "G", "A"),
#'   allele_noncoded = c("A", "T", "C"),
#'   stringsAsFactors = FALSE
#' )
#' result <- harmonizeAlleles(instruments, genotypeAlleles)
#'
#' @seealso \code{\link{isStrandAmbiguous}}, \code{\link{getMRInstruments}},
#'   \code{\link{buildMRCohort}}
#'
#' @export
harmonizeAlleles <- function(instrumentTable,
                             genotypeAlleles,
                             eafPalindromicThreshold = 0.08) {
  checkmate::assertDataFrame(instrumentTable, min.rows = 1)
  checkmate::assertSubset(c("snp_id", "effect_allele", "other_allele",
                            "beta_ZX", "se_ZX", "eaf"),
                          names(instrumentTable))
  checkmate::assertDataFrame(genotypeAlleles, min.rows = 1)
  checkmate::assertSubset(c("snp_id", "allele_coded", "allele_noncoded"),
                          names(genotypeAlleles))
  checkmate::assertNumber(eafPalindromicThreshold, lower = 0, upper = 0.5)

  removed <- data.frame(
    snp_id = character(0),
    reason = character(0),
    stringsAsFactors = FALSE
  )
  instrumentTable$flipped <- FALSE

  for (i in seq_len(nrow(instrumentTable))) {
    snpId <- instrumentTable$snp_id[i]
    ea <- toupper(instrumentTable$effect_allele[i])
    oa <- toupper(instrumentTable$other_allele[i])
    eaf <- instrumentTable$eaf[i]

    genoRow <- genotypeAlleles[genotypeAlleles$snp_id == snpId, , drop = FALSE]
    if (nrow(genoRow) == 0) {
      removed <- rbind(removed, data.frame(
        snp_id = snpId, reason = "not_in_genotype_data",
        stringsAsFactors = FALSE
      ))
      instrumentTable$flipped[i] <- NA
      next
    }

    genoCodedAllele <- toupper(genoRow$allele_coded[1])
    genoNoncodedAllele <- toupper(genoRow$allele_noncoded[1])

    # Check for palindromic SNPs FIRST (before allele matching)
    isPalindromic <- isStrandAmbiguous(ea, oa)
    if (isPalindromic && abs(eaf - 0.5) < eafPalindromicThreshold) {
      removed <- rbind(removed, data.frame(
        snp_id = snpId, reason = "palindromic_ambiguous_eaf",
        stringsAsFactors = FALSE
      ))
      instrumentTable$flipped[i] <- NA
      message(sprintf(
        "Removing palindromic SNP %s (EAF = %.3f, too close to 0.5).",
        snpId, eaf
      ))
      next
    }

    # Case 1: Direct match — effect allele matches coded allele
    if (ea == genoCodedAllele && oa == genoNoncodedAllele) {
      next
    }

    # Case 2: Alleles swapped — need to flip
    if (ea == genoNoncodedAllele && oa == genoCodedAllele) {
      instrumentTable$beta_ZX[i] <- -instrumentTable$beta_ZX[i]
      instrumentTable$eaf[i] <- 1 - instrumentTable$eaf[i]
      instrumentTable$effect_allele[i] <- genoCodedAllele
      instrumentTable$other_allele[i] <- genoNoncodedAllele
      instrumentTable$flipped[i] <- TRUE
      next
    }

    # Case 3: Try complement alleles
    eaComp <- complementAllele(ea)
    oaComp <- complementAllele(oa)

    # Try complement match
    if (eaComp == genoCodedAllele && oaComp == genoNoncodedAllele) {
      next
    }
    if (eaComp == genoNoncodedAllele && oaComp == genoCodedAllele) {
      instrumentTable$beta_ZX[i] <- -instrumentTable$beta_ZX[i]
      instrumentTable$eaf[i] <- 1 - instrumentTable$eaf[i]
      instrumentTable$effect_allele[i] <- genoCodedAllele
      instrumentTable$other_allele[i] <- genoNoncodedAllele
      instrumentTable$flipped[i] <- TRUE
      next
    }

    # Case 4: No match possible
    removed <- rbind(removed, data.frame(
      snp_id = snpId, reason = "allele_mismatch",
      stringsAsFactors = FALSE
    ))
    instrumentTable$flipped[i] <- NA
  }

  # Remove SNPs that couldn't be harmonized
  keepIdx <- !is.na(instrumentTable$flipped)
  instrumentTable <- instrumentTable[keepIdx, , drop = FALSE]
  rownames(instrumentTable) <- NULL

  if (nrow(removed) > 0) {
    message(sprintf("Removed %d SNPs during harmonization: %s",
                    nrow(removed),
                    paste(removed$snp_id, collapse = ", ")))
  }

  list(
    instrumentTable = instrumentTable,
    removedSnps = removed
  )
}


#' Get Complement of a DNA Allele
#'
#' @description Returns the Watson-Crick complement: A<->T, G<->C.
#'
#' @param allele Character string representing an allele (A, T, G, or C).
#'
#' @return Character string with the complement allele.
#'
#' @examples
#' complementAllele("A")  # returns "T"
#' complementAllele("G")  # returns "C"
#'
#' @export
complementAllele <- function(allele) {
  allele <- toupper(allele)
  mapping <- c(A = "T", T = "A", G = "C", C = "G")
  result <- mapping[allele]
  if (any(is.na(result))) {
    badAlleles <- allele[is.na(result)]
    stop(sprintf("Invalid allele(s): %s. Must be A, T, G, or C.",
                 paste(unique(badAlleles), collapse = ", ")))
  }
  unname(result)
}
