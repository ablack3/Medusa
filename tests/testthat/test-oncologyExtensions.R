test_that("formatExternalDosageTable standardizes wide dosage data", {
  instruments <- data.frame(
    snp_id = c("rs1", "rs2"),
    effect_allele = c("A", "G"),
    other_allele = c("C", "T"),
    beta_ZX = c(0.2, 0.3),
    se_ZX = c(0.05, 0.06),
    pval_ZX = c(1e-8, 1e-9),
    eaf = c(0.3, 0.4),
    stringsAsFactors = FALSE
  )
  dosageTable <- data.frame(
    person_id = c(1, 2),
    rs1 = c(0, 1),
    rs2 = c(2, 1),
    stringsAsFactors = FALSE
  )

  result <- formatExternalDosageTable(dosageTable, instruments)

  expect_s3_class(result, "medusaExternalDosage")
  expect_true(all(c("personId", "snp_rs1", "snp_rs2") %in% names(result$dosageTable)))
  expect_equal(result$dosageTable$snp_rs1, c(0L, 1L))
})

test_that("buildMRCohort accepts external dosage tables", {
  skip_if_not_installed("mockery")

  buildFn <- buildMRCohort
  connectionDetails <- structure(list(dbms = "postgresql"), class = "connectionDetails")
  fakeConnection <- structure(list(), class = "fakeConnection")
  instrumentTable <- data.frame(
    snp_id = c("rs1", "rs2"),
    effect_allele = c("A", "G"),
    other_allele = c("C", "T"),
    beta_ZX = c(0.2, 0.3),
    se_ZX = c(0.05, 0.06),
    pval_ZX = c(1e-8, 1e-9),
    eaf = c(0.3, 0.4),
    stringsAsFactors = FALSE
  )
  cohortFrame <- data.frame(
    personId = c(1L, 2L, 3L),
    outcome = c(1L, 0L, 1L),
    stringsAsFactors = FALSE
  )
  externalDosage <- data.frame(
    person_id = c(1L, 2L, 3L),
    rs1 = c(0, 1, 2),
    rs2 = c(2, 1, 0),
    stringsAsFactors = FALSE
  )

  mockery::stub(buildFn, "DatabaseConnector::connect", function(...) fakeConnection)
  mockery::stub(buildFn, "DatabaseConnector::disconnect", function(...) invisible(NULL))
  mockery::stub(buildFn, "loadRenderTranslateSql", function(...) "SELECT 1")
  mockery::stub(buildFn, "DatabaseConnector::executeSql", function(...) invisible(NULL))
  mockery::stub(buildFn, "DatabaseConnector::querySql", function(...) cohortFrame)

  result <- suppressMessages(
    suppressWarnings(
      buildFn(
        connectionDetails = connectionDetails,
        cdmDatabaseSchema = "cdm",
        cohortDatabaseSchema = "results",
        cohortTable = "cohort",
        outcomeCohortId = 1L,
        instrumentTable = instrumentTable,
        genotypeSource = "externalDosageTable",
        externalDosageTable = externalDosage
      )
    )
  )

  expect_equal(result$snp_rs1, c(0L, 1L, 2L))
  expect_identical(attr(result, "genotypeSource"), "externalDosageTable")
  expect_true(is.data.frame(attr(result, "harmonizedInstrumentTable")))
})

test_that("runCisMr returns LD-aware IVW and Egger summaries", {
  exposure <- data.frame(
    snp_id = c("rs1", "rs2", "rs3"),
    beta_ZX = c(0.3, 0.25, 0.2),
    se_ZX = c(0.05, 0.05, 0.04),
    stringsAsFactors = FALSE
  )
  outcome <- data.frame(
    snp_id = c("rs1", "rs2", "rs3"),
    beta_ZY = c(0.12, 0.10, 0.08),
    se_ZY = c(0.03, 0.03, 0.025),
    stringsAsFactors = FALSE
  )
  ld <- matrix(
    c(1, 0.2, 0.1,
      0.2, 1, 0.15,
      0.1, 0.15, 1),
    nrow = 3,
    byrow = TRUE
  )

  result <- runCisMr(exposure, outcome, ldMatrix = ld)

  expect_s3_class(result, "medusaCisMr")
  expect_true(all(c("GLS IVW", "GLS Egger") %in% result$summary$method))
  expect_equal(result$nVariants, 3)
})

test_that("runColocalization returns posterior probabilities", {
  exposure <- data.frame(
    snp_id = c("rs1", "rs2", "rs3"),
    beta = c(0.3, 0.18, 0.05),
    se = c(0.05, 0.05, 0.05),
    stringsAsFactors = FALSE
  )
  outcome <- data.frame(
    snp_id = c("rs1", "rs2", "rs3"),
    beta = c(0.28, 0.17, 0.03),
    se = c(0.06, 0.05, 0.05),
    stringsAsFactors = FALSE
  )

  result <- runColocalization(exposure, outcome)

  expect_s3_class(result, "medusaColoc")
  expect_true(all(c("H0", "H1", "H2", "H3", "H4") %in% names(result$posterior)))
  expect_true(result$ppH4 >= 0)
})

test_that("ancestry and overlap diagnostics classify high-risk settings", {
  ancestry <- runAncestryDiagnostics(
    exposurePopulation = c(EUR = 1),
    outcomePopulation = c(AFR = 1),
    ldReferencePopulation = c(EUR = 1)
  )
  overlap <- runOverlapDiagnostics(
    exposureSampleSize = 10000L,
    outcomeSampleSize = 8000L,
    overlapProportion = 0.4,
    meanFStatistic = 5
  )

  expect_true(ancestry$failClosed)
  expect_identical(overlap$riskLevel, "high")
})

test_that("study specs and decision classification support manifest execution", {
  spec <- createMedusaStudySpec(
    studyId = "demo",
    targetGene = "IL6R",
    exposureSource = "cis-pQTL",
    outcomeDefinition = createOncologyOutcomeDefinition("CRC"),
    datasetRole = "discovery"
  )
  decision <- classifyTargetDecision(
    studySpec = spec,
    cisMrResults = structure(
      list(
        fStatistics = c(15, 18),
        concordantDatasets = 2L
      ),
      class = "medusaCisMr"
    ),
    colocResults = structure(list(ppH4 = 0.9), class = "medusaColoc"),
    diagnosticResults = list(
      diagnosticFlags = c(
        negativeControlFailure = FALSE,
        phewasSignificant = FALSE
      )
    ),
    ancestryDiagnostics = structure(list(failClosed = FALSE), class = "medusaAncestryDiagnostics")
  )

  expect_identical(decision$recommendation, "advance")
  validateMedusaStudySpec(spec)
})

test_that("runStudyManifest wires executor output into summary rows", {
  manifest <- list(
    demo = createMedusaStudySpec(
      studyId = "demo",
      targetGene = "IL6R",
      exposureSource = "cis-pQTL",
      outcomeDefinition = createOncologyOutcomeDefinition("CRC"),
      datasetRole = "discovery"
    )
  )

  result <- runStudyManifest(
    studyManifest = manifest,
    executor = function(studySpec) {
      list(
        cisMrResults = structure(
          list(
            summary = data.frame(method = "GLS IVW", beta = 0.2, se = 0.05),
            fStatistics = c(15, 16),
            concordantDatasets = 2L
          ),
          class = "medusaCisMr"
        ),
        colocResults = structure(
          list(ppH4 = 0.9, posterior = c(H0 = 0.01, H1 = 0.02, H2 = 0.02, H3 = 0.05, H4 = 0.9)),
          class = "medusaColoc"
        ),
        diagnosticResults = list(
          diagnosticFlags = c(
            negativeControlFailure = FALSE,
            phewasSignificant = FALSE
          )
        ),
        ancestryDiagnostics = structure(list(failClosed = FALSE), class = "medusaAncestryDiagnostics"),
        overlapDiagnostics = NULL,
        trialContext = paste("Study", studySpec$studyId)
      )
    },
    generateReports = FALSE
  )

  expect_s3_class(result, "medusaManifestResults")
  expect_identical(result$summary$recommendation, "advance")
})

test_that("generateTargetValidationReport renders with a stubbed renderer", {
  skip_if_not_installed("mockery")

  reportFn <- generateTargetValidationReport
  outputPath <- file.path(tempdir(), "medusa-target-report.html")
  if (file.exists(outputPath)) {
    unlink(outputPath)
  }

  mockery::stub(
    reportFn,
    "renderMedusaHtmlReport",
    function(templateName, outputPath, ...) {
      writeLines("<html><body>ok</body></html>", outputPath)
      invisible(outputPath)
    }
  )

  result <- suppressMessages(
    reportFn(
      studySpec = createMedusaStudySpec(
        studyId = "demo",
        targetGene = "IL6R",
        exposureSource = "cis-pQTL",
        outcomeDefinition = createOncologyOutcomeDefinition("CRC"),
        datasetRole = "discovery"
      ),
      cisMrResults = structure(list(fStatistics = c(12, 14), concordantDatasets = 2L), class = "medusaCisMr"),
      colocResults = structure(list(ppH4 = 0.85), class = "medusaColoc"),
      diagnosticResults = list(
        diagnosticFlags = c(
          negativeControlFailure = FALSE,
          phewasSignificant = FALSE
        )
      ),
      ancestryDiagnostics = structure(list(failClosed = FALSE), class = "medusaAncestryDiagnostics"),
      outputPath = outputPath
    )
  )

  expect_equal(result, outputPath)
  expect_true(file.exists(outputPath))
})
