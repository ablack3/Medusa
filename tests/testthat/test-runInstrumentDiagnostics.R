test_that("runInstrumentDiagnostics returns diagnostic components for a clean cohort", {
  simData <- simulateMRData(n = 300, nSnps = 3, seed = 505)
  cohortData <- simData$data
  cohortData$personId <- cohortData$person_id
  covariateData <- data.frame(
    person_id = cohortData$person_id,
    binary_cov = rep(c(0, 1), length.out = nrow(cohortData)),
    continuous_cov = rnorm(nrow(cohortData)),
    stringsAsFactors = FALSE
  )

  diagnostics <- suppressMessages(
    runInstrumentDiagnostics(
      cohortData = cohortData,
      covariateData = covariateData,
      instrumentTable = simData$instrumentTable,
      exposureProxyConceptIds = 1L,
      negativeControlOutcomeIds = c(10L, 11L)
    )
  )

  expect_s3_class(diagnostics, "medusaDiagnostics")
  expect_true(all(c(
    "fStatistics", "phewasResults", "negativeControlResults",
    "afComparison", "missingnessReport", "diagnosticFlags"
  ) %in% names(diagnostics)))
  expect_equal(nrow(diagnostics$fStatistics), 3)
  expect_equal(nrow(diagnostics$afComparison), 3)
  expect_equal(nrow(diagnostics$missingnessReport), 3)
})

test_that("runInstrumentDiagnostics errors when cohortData has no SNP columns", {
  simData <- simulateMRData(n = 50, nSnps = 1, seed = 509)
  cohortData <- simData$data[, c("person_id", "outcome"), drop = FALSE]

  expect_error(
    runInstrumentDiagnostics(
      cohortData = cohortData,
      covariateData = NULL,
      instrumentTable = simData$instrumentTable
    ),
    "No SNP columns"
  )
})

test_that("computeFStatistics falls back to GWAS approximations when cohort fits are unavailable", {
  simData <- simulateMRData(n = 100, nSnps = 2, seed = 515)
  cohortData <- simData$data
  snpCols <- grep("^snp_", names(cohortData), value = TRUE)
  cohortData[[snpCols[2]]] <- NA_integer_

  fStats <- computeFStatistics(
    cohortData = cohortData,
    instrumentTable = simData$instrumentTable,
    exposureProxyConceptIds = 1L
  )

  expect_equal(fStats$source[1], "cohort_data")
  expect_equal(fStats$source[2], "gwas_approximation")
})

test_that("runInstrumentDiagnostics flags allele frequency and missingness problems", {
  simData <- simulateMRData(n = 120, nSnps = 2, seed = 606)
  cohortData <- simData$data
  cohortData$personId <- cohortData$person_id
  snpCols <- grep("^snp_", names(cohortData), value = TRUE)
  cohortData[[snpCols[1]]] <- 2L
  cohortData[[snpCols[2]]][1:20] <- NA_integer_
  simData$instrumentTable$eaf[1] <- 0.1
  covariateData <- data.frame(
    person_id = cohortData$person_id,
    dummy_cov = rep(c(0, 1), length.out = nrow(cohortData)),
    stringsAsFactors = FALSE
  )

  expect_warning(
    diagnostics <- suppressMessages(
      runInstrumentDiagnostics(
        cohortData = cohortData,
        covariateData = covariateData,
        instrumentTable = simData$instrumentTable
      )
    ),
    "Instrument diagnostics flagged potential issues"
  )

  expect_true(isTRUE(diagnostics$diagnosticFlags["alleleFreqDiscrepancy"]))
  expect_true(isTRUE(diagnostics$diagnosticFlags["highMissingness"]))
})

test_that("runInstrumentPheWAS supports medusaCovariateData inputs and empty cases", {
  simData <- simulateMRData(n = 150, nSnps = 2, seed = 707)
  cohortData <- simData$data
  cohortData$personId <- cohortData$person_id

  medusaCovariates <- list(
    covariateData = list(
      covariates = data.frame(
        person_id = cohortData$person_id,
        covariate_binary = rep(c(0, 1), length.out = nrow(cohortData)),
        stringsAsFactors = FALSE
      ),
      covariateRef = data.frame(
        covariateId = 1L,
        covariateName = "covariate_binary",
        domainId = "Condition",
        stringsAsFactors = FALSE
      )
    )
  )
  class(medusaCovariates) <- "medusaCovariateData"

  phewas <- runInstrumentPheWAS(
    cohortData = cohortData,
    covariateData = medusaCovariates,
    instrumentTable = simData$instrumentTable
  )
  emptyPhewas <- runInstrumentPheWAS(
    cohortData = cohortData,
    covariateData = NULL,
    instrumentTable = simData$instrumentTable
  )

  expect_s3_class(phewas, "data.frame")
  expect_true(all(c("snp_id", "covariate_name", "pval") %in% names(phewas)))
  expect_equal(nrow(emptyPhewas), 0)
})

test_that("runInstrumentPheWAS covers personId joins, no-join fallback, and skipped regressions", {
  simData <- simulateMRData(n = 60, nSnps = 1, seed = 709)
  cohortData <- simData$data
  cohortData$personId <- cohortData$person_id

  joinedPhewas <- runInstrumentPheWAS(
    cohortData = cohortData,
    covariateData = data.frame(
      personId = cohortData$personId,
      constant_binary = 1L,
      stringsAsFactors = FALSE
    ),
    instrumentTable = simData$instrumentTable
  )
  skippedPhewas <- runInstrumentPheWAS(
    cohortData = cohortData,
    covariateData = data.frame(
      no_join_key = seq_len(nrow(cohortData)),
      sparse_covariate = c(rep(NA_real_, nrow(cohortData) - 5L), rep(1, 5L)),
      stringsAsFactors = FALSE
    ),
    instrumentTable = simData$instrumentTable
  )

  expect_equal(nrow(joinedPhewas), 0)
  expect_equal(nrow(skippedPhewas), 0)
})

test_that("negative controls, allele frequencies, and missingness helpers return expected frames", {
  simData <- simulateMRData(n = 100, nSnps = 2, seed = 808)
  cohortData <- simData$data

  negativeControls <- testNegativeControls(
    cohortData = cohortData,
    instrumentTable = simData$instrumentTable,
    negativeControlOutcomeIds = c(10L, 11L)
  )
  afComparison <- compareAlleleFrequencies(cohortData, simData$instrumentTable)
  missingness <- summarizeGenotypeMissingness(cohortData, simData$instrumentTable)

  expect_equal(nrow(negativeControls), 0)
  expect_equal(nrow(afComparison), 2)
  expect_equal(nrow(missingness), 2)
  expect_true(all(missingness$n_total == nrow(cohortData)))
})
