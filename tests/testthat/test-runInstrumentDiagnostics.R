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

test_that("computeFStatistics matches SNPs by SNP ID rather than cohort column order", {
  simData <- simulateMRData(n = 200, nSnps = 2, seed = 516)
  snpCols <- grep("^snp_", names(simData$data), value = TRUE)
  reordered <- simData$data[, c(setdiff(names(simData$data), snpCols), rev(snpCols))]

  original <- computeFStatistics(
    cohortData = simData$data,
    instrumentTable = simData$instrumentTable,
    exposureProxyConceptIds = 1L
  )
  reorderedStats <- computeFStatistics(
    cohortData = reordered,
    instrumentTable = simData$instrumentTable,
    exposureProxyConceptIds = 1L
  )

  expect_equal(original$fStatistic, reorderedStats$fStatistic, tolerance = 1e-10)
  expect_equal(original$source, reorderedStats$source)
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

test_that("runInstrumentPheWAS handles sparse FeatureExtraction format", {
  simData <- simulateMRData(n = 150, nSnps = 2, seed = 710)
  cohortData <- simData$data
  cohortData$personId <- cohortData$person_id

  # Create sparse long format matching real FeatureExtraction output
  # Only persons with the condition have a row (covariateValue = 1)
  affectedIds <- cohortData$personId[cohortData$personId %% 3 == 0]

  medusaCovariates <- list(
    covariateData = list(
      covariates = data.frame(
        rowId = affectedIds,
        covariateId = rep(2001L, length(affectedIds)),
        covariateValue = rep(1, length(affectedIds)),
        stringsAsFactors = FALSE
      ),
      covariateRef = data.frame(
        covariateId = 2001L,
        covariateName = "condition_group_era:diabetes",
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

  expect_s3_class(phewas, "data.frame")
  expect_true(all(c("snp_id", "covariate_name", "pval") %in% names(phewas)))
  # One model with 2 SNPs → 2 result rows
  expect_equal(nrow(phewas), 2)
  expect_true(all(phewas$covariate_name == "condition_group_era:diabetes"))
})

test_that("runInstrumentPheWAS handles multiple sparse covariates", {
  simData <- simulateMRData(n = 200, nSnps = 2, seed = 711)
  cohortData <- simData$data
  cohortData$personId <- cohortData$person_id

  # Two covariates in sparse format
  ids1 <- cohortData$personId[cohortData$personId %% 2 == 0]
  ids2 <- cohortData$personId[cohortData$personId %% 4 == 0]

  medusaCovariates <- list(
    covariateData = list(
      covariates = data.frame(
        rowId = c(ids1, ids2),
        covariateId = c(rep(3001L, length(ids1)), rep(3002L, length(ids2))),
        covariateValue = rep(1, length(ids1) + length(ids2)),
        stringsAsFactors = FALSE
      ),
      covariateRef = data.frame(
        covariateId = c(3001L, 3002L),
        covariateName = c("condition:hypertension", "drug:aspirin"),
        domainId = c("Condition", "Drug"),
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

  expect_s3_class(phewas, "data.frame")
  # 2 covariates × 2 SNPs = 4 result rows
  expect_equal(nrow(phewas), 4)
  expect_true(all(c("condition:hypertension", "drug:aspirin") %in% phewas$covariate_name))
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

test_that("testNegativeControls honors requested negative control IDs", {
  simData <- simulateMRData(n = 300, nSnps = 3, seed = 809)
  cohortData <- simulateNegativeControlOutcomes(simData$data, nControls = 4, seed = 810)

  selected <- suppressMessages(
    testNegativeControls(
      cohortData = cohortData,
      instrumentTable = simData$instrumentTable,
      negativeControlOutcomeIds = c(2L, 4L)
    )
  )

  expect_equal(sort(selected$outcome_id), c("2", "4"))
  expect_false("1" %in% selected$outcome_id)
  expect_false("3" %in% selected$outcome_id)
})

test_that("runInstrumentDiagnostics uses aggregate negative-control bias flag", {
  skip_if_not_installed("mockery")

  diagnosticsFn <- runInstrumentDiagnostics
  simData <- simulateMRData(n = 80, nSnps = 2, seed = 811)
  cohortData <- simData$data
  cohortData$personId <- cohortData$person_id
  covariateData <- data.frame(
    person_id = cohortData$person_id,
    binary_cov = rep(c(0, 1), length.out = nrow(cohortData)),
    stringsAsFactors = FALSE
  )

  mockery::stub(
    diagnosticsFn,
    "computeFStatistics",
    function(...) data.frame(
      snp_id = simData$instrumentTable$snp_id,
      fStatistic = c(20, 25),
      source = "gwas_approximation",
      weakFlag = c(FALSE, FALSE),
      stringsAsFactors = FALSE
    )
  )
  mockery::stub(
    diagnosticsFn,
    "runInstrumentPheWAS",
    function(...) data.frame(
      snp_id = character(0),
      covariate_name = character(0),
      beta = numeric(0),
      se = numeric(0),
      pval = numeric(0),
      significant = logical(0),
      stringsAsFactors = FALSE
    )
  )
  mockery::stub(
    diagnosticsFn,
    "testNegativeControls",
    function(...) {
      out <- data.frame(
        outcome_id = "201",
        beta_ZY = 0.1,
        se_ZY = 0.05,
        beta_MR = 0.2,
        se_MR = 0.1,
        pval = 0.01,
        stringsAsFactors = FALSE
      )
      attr(out, "biasDetected") <- FALSE
      out
    }
  )
  mockery::stub(
    diagnosticsFn,
    "compareAlleleFrequencies",
    function(...) data.frame(
      snp_id = simData$instrumentTable$snp_id,
      eaf_gwas = simData$instrumentTable$eaf,
      eaf_cohort = simData$instrumentTable$eaf,
      eaf_diff = c(0, 0),
      discrepancyFlag = c(FALSE, FALSE),
      stringsAsFactors = FALSE
    )
  )
  mockery::stub(
    diagnosticsFn,
    "summarizeGenotypeMissingness",
    function(...) data.frame(
      snp_id = simData$instrumentTable$snp_id,
      n_total = nrow(cohortData),
      n_missing = c(0L, 0L),
      pct_missing = c(0, 0),
      highMissingFlag = c(FALSE, FALSE),
      stringsAsFactors = FALSE
    )
  )

  diagnostics <- suppressMessages(
    diagnosticsFn(
      cohortData = cohortData,
      covariateData = covariateData,
      instrumentTable = simData$instrumentTable,
      negativeControlOutcomeIds = 201L
    )
  )

  expect_false(isTRUE(diagnostics$diagnosticFlags["negativeControlFailure"]))
  expect_equal(diagnostics$negativeControlResults$pval, 0.01)
})
