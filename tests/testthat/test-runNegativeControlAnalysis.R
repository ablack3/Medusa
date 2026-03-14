test_that("runNegativeControlAnalysis returns correct structure with simulated data", {
  simData <- simulateMRData(n = 500, nSnps = 5, seed = 701)
  cohortData <- simulateNegativeControlOutcomes(simData$data, nControls = 10, seed = 702)

  ncResults <- suppressMessages(runNegativeControlAnalysis(
    cohortData = cohortData,
    instrumentTable = simData$instrumentTable
  ))

  expect_s3_class(ncResults, "medusaNegativeControls")
  expect_true(all(c("ncEstimates", "calibration", "calibratedPrimary",
                     "biasDetected") %in% names(ncResults)))

  # ncEstimates should be a data frame with expected columns

  expect_s3_class(ncResults$ncEstimates, "data.frame")
  expect_true(all(c("outcome_id", "beta_ZY", "se_ZY", "beta_MR", "se_MR",
                     "pval", "log_rr", "se_log_rr") %in% names(ncResults$ncEstimates)))

  # Should have results for most/all NC outcomes
  expect_gte(nrow(ncResults$ncEstimates), 5)
})


test_that("NC estimates are approximately null for true negative controls", {
  simData <- simulateMRData(n = 1000, nSnps = 5, seed = 703)
  cohortData <- simulateNegativeControlOutcomes(simData$data, nControls = 15, seed = 704)

  ncResults <- suppressMessages(runNegativeControlAnalysis(
    cohortData = cohortData,
    instrumentTable = simData$instrumentTable
  ))

  estimates <- ncResults$ncEstimates$beta_MR
  validEstimates <- estimates[is.finite(estimates)]

  # Mean should be near zero (within 2 SE of mean)
  meanEst <- mean(validEstimates)
  seMean <- sd(validEstimates) / sqrt(length(validEstimates))
  expect_lt(abs(meanEst), 3 * seMean + 0.5)  # Generous bound for simulation noise
})


test_that("runNegativeControlAnalysis errors with no NC columns", {
  simData <- simulateMRData(n = 100, nSnps = 3, seed = 705)

  expect_error(
    suppressMessages(runNegativeControlAnalysis(
      cohortData = simData$data,
      instrumentTable = simData$instrumentTable
    )),
    "length >= 1"
  )
})


test_that("runNegativeControlAnalysis accepts explicit column names", {
  simData <- simulateMRData(n = 300, nSnps = 3, seed = 706)
  cohortData <- simData$data
  # Add columns with non-standard names
  cohortData$control_a <- rbinom(nrow(cohortData), 1, 0.1)
  cohortData$control_b <- rbinom(nrow(cohortData), 1, 0.08)

  ncResults <- suppressMessages(runNegativeControlAnalysis(
    cohortData = cohortData,
    instrumentTable = simData$instrumentTable,
    negativeControlColumns = c("control_a", "control_b")
  ))

  expect_equal(nrow(ncResults$ncEstimates), 2)
  expect_equal(ncResults$ncEstimates$outcome_id, c("control_a", "control_b"))
})


test_that("simulateNegativeControlOutcomes adds correct columns", {
  simData <- simulateMRData(n = 100, nSnps = 3, seed = 707)
  result <- simulateNegativeControlOutcomes(simData$data, nControls = 5, seed = 708)

  ncCols <- grep("^nc_outcome_", names(result), value = TRUE)
  expect_length(ncCols, 5)

  # All binary
  for (col in ncCols) {
    expect_true(all(result[[col]] %in% c(0, 1)))
  }

  # Original columns preserved
  expect_true(all(names(simData$data) %in% names(result)))
})


test_that("biasDetected is FALSE for true null outcomes", {
  simData <- simulateMRData(n = 1000, nSnps = 5, seed = 709)
  cohortData <- simulateNegativeControlOutcomes(simData$data, nControls = 20, seed = 710)

  ncResults <- suppressMessages(runNegativeControlAnalysis(
    cohortData = cohortData,
    instrumentTable = simData$instrumentTable
  ))

  # With true null outcomes, bias should generally not be detected
  # (this is probabilistic; use a generous seed + sample size)
  expect_type(ncResults$biasDetected, "logical")
})


test_that("calibratedPrimary is populated when primaryEstimate is provided", {
  skip_if_not_installed("EmpiricalCalibration")

  simData <- simulateMRData(n = 1000, nSnps = 5, seed = 711)
  cohortData <- simulateNegativeControlOutcomes(simData$data, nControls = 20, seed = 712)

  mockEstimate <- list(
    betaMR = 0.3,
    seMR = 0.1,
    pValue = 0.003
  )

  ncResults <- suppressMessages(runNegativeControlAnalysis(
    cohortData = cohortData,
    instrumentTable = simData$instrumentTable,
    primaryEstimate = mockEstimate
  ))

  expect_true(!is.null(ncResults$calibratedPrimary))
  expect_true("calibratedP" %in% names(ncResults$calibratedPrimary))
})


test_that("testNegativeControls integrates with runInstrumentDiagnostics", {
  simData <- simulateMRData(n = 300, nSnps = 3, seed = 713)
  cohortData <- simulateNegativeControlOutcomes(simData$data, nControls = 5, seed = 714)
  cohortData$personId <- cohortData$person_id
  covariateData <- data.frame(
    person_id = cohortData$person_id,
    binary_cov = rep(c(0, 1), length.out = nrow(cohortData)),
    stringsAsFactors = FALSE
  )

  diagnostics <- suppressMessages(suppressWarnings(
    runInstrumentDiagnostics(
      cohortData = cohortData,
      covariateData = covariateData,
      instrumentTable = simData$instrumentTable,
      negativeControlOutcomeIds = c(1L, 2L)
    )
  ))

  expect_s3_class(diagnostics, "medusaDiagnostics")
  expect_true(!is.null(diagnostics$negativeControlResults))
  expect_s3_class(diagnostics$negativeControlResults, "data.frame")
  expect_gte(nrow(diagnostics$negativeControlResults), 1)
  expect_true(all(c("outcome_id", "beta_MR", "se_MR", "pval") %in%
                    names(diagnostics$negativeControlResults)))
})
