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

#' Create a Medusa study specification
#'
#' @title Study-spec object for target-validation workflows
#'
#' @description
#' Creates a lightweight, JSON-serializable study specification that can drive
#' manifest-based Medusa analyses. The object captures the target, exposure
#' source, outcome definition, dataset role, and target-triage thresholds used
#' throughout the oncology program.
#'
#' @param studyId Character. Unique study identifier.
#' @param targetGene Character vector naming the primary gene or pathway target.
#' @param exposureSource Character. Human-readable source for the primary
#'   exposure summary statistics.
#' @param instrumentMode Character. Instrument mode for the study. Defaults to
#'   \code{"cis"}.
#' @param outcomeDefinition Named list describing the phenotype and endpoint.
#' @param datasetRole Character. One of \code{"discovery"},
#'   \code{"replication"}, or \code{"killStudy"}.
#' @param primaryPopulation Character. Primary ancestry/population label.
#' @param negativeControls Character vector of requested negative-control
#'   outcomes.
#' @param decisionThresholds Named list of evidence gates used for final target
#'   classification.
#' @param metadata Optional named list with study-level notes.
#'
#' @return A list with class \code{"medusaStudySpec"}.
#'
#' @export
createMedusaStudySpec <- function(studyId,
                                  targetGene,
                                  exposureSource,
                                  instrumentMode = "cis",
                                  outcomeDefinition,
                                  datasetRole = c("discovery", "replication", "killStudy"),
                                  primaryPopulation = "EUR",
                                  negativeControls = character(0),
                                  decisionThresholds = NULL,
                                  metadata = list()) {
  datasetRole <- match.arg(datasetRole)
  checkmate::assertString(studyId)
  checkmate::assertCharacter(targetGene, min.len = 1, any.missing = FALSE)
  checkmate::assertString(exposureSource)
  checkmate::assertString(instrumentMode)
  checkmate::assertList(outcomeDefinition, names = "named")
  checkmate::assertString(primaryPopulation)
  checkmate::assertCharacter(negativeControls, any.missing = FALSE)
  checkmate::assertList(metadata)

  decisionThresholds <- decisionThresholds %||% createDefaultDecisionThresholds()
  validateDecisionThresholds(decisionThresholds)

  spec <- list(
    studyId = studyId,
    targetGene = unname(targetGene),
    exposureSource = exposureSource,
    instrumentMode = instrumentMode,
    outcomeDefinition = outcomeDefinition,
    datasetRole = datasetRole,
    primaryPopulation = primaryPopulation,
    negativeControls = unname(negativeControls),
    decisionThresholds = decisionThresholds,
    metadata = metadata
  )
  class(spec) <- "medusaStudySpec"
  spec
}


#' Validate a Medusa study specification
#'
#' @param studySpec Study specification created by
#'   \code{\link{createMedusaStudySpec}}.
#'
#' @return Invisible TRUE.
#'
#' @export
validateMedusaStudySpec <- function(studySpec) {
  checkmate::assertList(studySpec, names = "named")
  required <- c(
    "studyId", "targetGene", "exposureSource", "instrumentMode",
    "outcomeDefinition", "datasetRole", "primaryPopulation",
    "negativeControls", "decisionThresholds", "metadata"
  )
  missing <- setdiff(required, names(studySpec))
  if (length(missing) > 0) {
    stop(sprintf(
      "studySpec is missing required fields: %s",
      paste(missing, collapse = ", ")
    ))
  }

  checkmate::assertString(studySpec$studyId)
  checkmate::assertCharacter(studySpec$targetGene, min.len = 1, any.missing = FALSE)
  checkmate::assertString(studySpec$exposureSource)
  checkmate::assertString(studySpec$instrumentMode)
  checkmate::assertChoice(studySpec$datasetRole, c("discovery", "replication", "killStudy"))
  checkmate::assertString(studySpec$primaryPopulation)
  checkmate::assertCharacter(studySpec$negativeControls, any.missing = FALSE)
  checkmate::assertList(studySpec$outcomeDefinition, names = "named")
  validateDecisionThresholds(studySpec$decisionThresholds)
  invisible(TRUE)
}


#' Create default evidence thresholds for target decisions
#'
#' @return Named list of decision thresholds.
#'
#' @export
createDefaultDecisionThresholds <- function() {
  list(
    minFStatistic = 10,
    minColocPp4 = 0.8,
    minConcordantDatasets = 2L,
    failOnAncestryMismatch = TRUE,
    failOnNegativeControls = TRUE
  )
}


#' Build the default 8-study oncology manifest
#'
#' @title Focused oncology target-validation portfolio
#'
#' @description
#' Returns the focused eight-study manifest described in the Medusa oncology
#' roadmap. The output is a list of \code{medusaStudySpec} objects that can be
#' passed to \code{\link{runStudyManifest}}.
#'
#' @return Named list of \code{medusaStudySpec} objects.
#'
#' @export
createOncologyStudyManifest <- function() {
  specs <- list(
    createMedusaStudySpec(
      studyId = "wave1_il6r_crc",
      targetGene = "IL6R",
      exposureSource = "cis-pQTL/eQTL",
      outcomeDefinition = list(
        label = "Incident colorectal cancer",
        cancer = "CRC",
        endpointType = "incidentCancer"
      ),
      datasetRole = "discovery",
      negativeControls = c("appendicitis", "cataract"),
      metadata = list(
        wave = 1L,
        datasets = c("UK Biobank", "FinnGen", "All of Us"),
        decisionUse = "Prioritize or deprioritize IL-6 blockade for CRC prevention/interception."
      )
    ),
    createMedusaStudySpec(
      studyId = "wave1_il6r_nsclc",
      targetGene = "IL6R",
      exposureSource = "cis-pQTL/eQTL",
      outcomeDefinition = list(
        label = "Incident NSCLC",
        cancer = "NSCLC",
        endpointType = "incidentCancer",
        secondaryEndpointType = "systemicTreatmentWithin12Months"
      ),
      datasetRole = "discovery",
      negativeControls = c("kidney_stone", "osteoarthritis"),
      metadata = list(
        wave = 1L,
        datasets = c("UK Biobank", "FinnGen", "All of Us"),
        decisionUse = "Disentangle etiologic from treatment-resistance IL-6 biology in lung cancer."
      )
    ),
    createMedusaStudySpec(
      studyId = "wave1_nt5e_pdac",
      targetGene = "NT5E",
      exposureSource = "cis-pQTL/eQTL",
      outcomeDefinition = list(
        label = "Incident pancreatic adenocarcinoma",
        cancer = "PDAC",
        endpointType = "incidentCancer"
      ),
      datasetRole = "discovery",
      negativeControls = c("diverticulitis", "migraine"),
      metadata = list(
        wave = 1L,
        datasets = c("UK Biobank", "FinnGen", "All of Us"),
        decisionUse = "Assess indication-specific support for CD73 inhibition in pancreas."
      )
    ),
    createMedusaStudySpec(
      studyId = "wave1_entpd1_panbili",
      targetGene = "ENTPD1",
      exposureSource = "cis-pQTL/eQTL",
      outcomeDefinition = list(
        label = "Pancreaticobiliary cancer composite",
        cancer = "Pancreaticobiliary",
        endpointType = "incidentCancer",
        sensitivityEndpoint = "cholangiocarcinoma"
      ),
      datasetRole = "discovery",
      negativeControls = c("eczema", "gallstones"),
      metadata = list(
        wave = 1L,
        datasets = c("FinnGen", "UK Biobank", "All of Us"),
        decisionUse = "Determine whether CD39 merits indication-specific development."
      )
    ),
    createMedusaStudySpec(
      studyId = "wave1_cxcl8_hcc",
      targetGene = c("CXCL8", "CXCR2"),
      exposureSource = "cis-eQTL/pQTL",
      outcomeDefinition = list(
        label = "Incident hepatocellular carcinoma",
        cancer = "HCC",
        endpointType = "incidentCancer"
      ),
      datasetRole = "discovery",
      negativeControls = c("benign_prostatic_hyperplasia", "herpes_zoster"),
      metadata = list(
        wave = 1L,
        datasets = c("UK Biobank", "FinnGen", "All of Us"),
        decisionUse = "Test whether IL-8 axis has etiologic support in HCC."
      )
    ),
    createMedusaStudySpec(
      studyId = "wave1_ccr2_csf1r_pdac",
      targetGene = c("CCR2", "CSF1R"),
      exposureSource = "cis-eQTL/pQTL",
      outcomeDefinition = list(
        label = "Incident pancreatic adenocarcinoma",
        cancer = "PDAC",
        endpointType = "incidentCancer"
      ),
      datasetRole = "discovery",
      negativeControls = c("gout", "retinal_detachment"),
      metadata = list(
        wave = 1L,
        datasets = c("UK Biobank", "FinnGen"),
        decisionUse = "Support or downgrade macrophage-targeting combination strategies in PDAC."
      )
    ),
    createMedusaStudySpec(
      studyId = "wave2_tgfb_crc_hb",
      targetGene = c("TGFB1", "TGFBR1", "TGFBR2"),
      exposureSource = "cis-eQTL/pQTL",
      outcomeDefinition = list(
        label = "CRC and hepatobiliary cancers",
        cancer = "CRC_Hepatobiliary",
        endpointType = "incidentCancer"
      ),
      datasetRole = "replication",
      negativeControls = c("allergic_rhinitis", "nephrolithiasis"),
      metadata = list(
        wave = 2L,
        datasets = c("UK Biobank", "FinnGen", "All of Us"),
        decisionUse = "Re-evaluate TGF-beta combination relevance after robust cis-MR and coloc."
      )
    ),
    createMedusaStudySpec(
      studyId = "wave2_tigit_nsclc",
      targetGene = c("TIGIT", "PVR", "CD226"),
      exposureSource = "cis-eQTL/pQTL",
      outcomeDefinition = list(
        label = "Incident NSCLC",
        cancer = "NSCLC",
        endpointType = "incidentCancer",
        stratifications = c("smoking", "sex")
      ),
      datasetRole = "killStudy",
      negativeControls = c("inguinal_hernia", "varicose_veins"),
      metadata = list(
        wave = 2L,
        datasets = c("UK Biobank", "FinnGen", "All of Us"),
        decisionUse = "Produce an explicit stop-or-restrict recommendation for TIGIT-class development."
      )
    )
  )

  names(specs) <- vapply(specs, `[[`, character(1), "studyId")
  specs
}


#' Run a study manifest with a user-supplied executor
#'
#' @title Manifest-based target-validation orchestration
#'
#' @description
#' Iterates over a list of study specifications, invokes a user-supplied
#' execution function for each study, classifies the resulting evidence, and
#' optionally renders target-validation reports.
#'
#' @param studyManifest List of \code{medusaStudySpec} objects.
#' @param executor Function taking one \code{studySpec} and returning a named
#'   list of result objects.
#' @param outputDir Character. Directory for report outputs.
#' @param generateReports Logical. Whether to render target-validation reports.
#'
#' @return A list with class \code{"medusaManifestResults"}.
#'
#' @export
runStudyManifest <- function(studyManifest,
                             executor,
                             outputDir = tempdir(),
                             generateReports = TRUE) {
  checkmate::assertList(studyManifest, min.len = 1)
  checkmate::assertFunction(executor)
  checkmate::assertString(outputDir)
  checkmate::assertLogical(generateReports, len = 1)

  if (!dir.exists(outputDir)) {
    dir.create(outputDir, recursive = TRUE)
  }

  results <- vector("list", length(studyManifest))
  names(results) <- names(studyManifest)

  for (i in seq_along(studyManifest)) {
    studySpec <- studyManifest[[i]]
    validateMedusaStudySpec(studySpec)

    studyResult <- executor(studySpec)
    checkmate::assertList(studyResult, names = "named")

    decision <- classifyTargetDecision(
      studySpec = studySpec,
      cisMrResults = studyResult$cisMrResults,
      colocResults = studyResult$colocResults,
      diagnosticResults = studyResult$diagnosticResults,
      ancestryDiagnostics = studyResult$ancestryDiagnostics
    )

    reportPath <- NULL
    if (isTRUE(generateReports)) {
      reportPath <- generateTargetValidationReport(
        studySpec = studySpec,
        cisMrResults = studyResult$cisMrResults,
        colocResults = studyResult$colocResults,
        diagnosticResults = studyResult$diagnosticResults,
        ancestryDiagnostics = studyResult$ancestryDiagnostics,
        overlapDiagnostics = studyResult$overlapDiagnostics,
        trialContext = studyResult$trialContext,
        outputPath = file.path(outputDir, paste0(studySpec$studyId, ".html"))
      )
    }

    results[[i]] <- c(
      list(
        studySpec = studySpec,
        decision = decision,
        reportPath = reportPath
      ),
      studyResult
    )
  }

  structure(
    list(
      results = results,
      summary = do.call(
        rbind,
        lapply(results, function(x) {
          data.frame(
            studyId = x$studySpec$studyId,
            recommendation = x$decision$recommendation,
            rationale = paste(x$decision$reasons, collapse = "; "),
            stringsAsFactors = FALSE
          )
        })
      )
    ),
    class = "medusaManifestResults"
  )
}


#' @keywords internal
validateDecisionThresholds <- function(decisionThresholds) {
  checkmate::assertList(decisionThresholds, names = "named")
  required <- c(
    "minFStatistic", "minColocPp4", "minConcordantDatasets",
    "failOnAncestryMismatch", "failOnNegativeControls"
  )
  missing <- setdiff(required, names(decisionThresholds))
  if (length(missing) > 0) {
    stop(sprintf(
      "decisionThresholds is missing required fields: %s",
      paste(missing, collapse = ", ")
    ))
  }
  invisible(TRUE)
}
