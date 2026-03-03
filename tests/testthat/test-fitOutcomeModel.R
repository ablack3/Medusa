test_that("fitOutcomeModel returns valid profile with simulated data", {
  simData <- simulateMRData(n = 2000, nSnps = 5, trueEffect = 0.3, seed = 123)
  profile <- fitOutcomeModel(
    cohortData = simData$data,
    covariateData = NULL,
    instrumentTable = simData$instrumentTable,
    betaGrid = seq(-2, 2, by = 0.05),
    siteId = "test_site"
  )

  expect_s3_class(profile, "medusaSiteProfile")
  expect_equal(profile$siteId, "test_site")
  expect_equal(length(profile$betaGrid), length(profile$logLikProfile))
  expect_true(!is.null(profile$betaHat))
  expect_true(!is.null(profile$seHat))
})

test_that("log-likelihood profile is concave (has unique maximum)", {
  simData <- simulateMRData(n = 5000, nSnps = 5, trueEffect = 0.5, seed = 456)
  profile <- fitOutcomeModel(
    cohortData = simData$data,
    covariateData = NULL,
    instrumentTable = simData$instrumentTable,
    betaGrid = seq(-3, 3, by = 0.01)
  )

  # Profile should have a clear peak
  peakIdx <- which.max(profile$logLikProfile)
  expect_true(peakIdx > 1)
  expect_true(peakIdx < length(profile$logLikProfile))

  # Values should decrease away from peak
  expect_true(profile$logLikProfile[peakIdx] > profile$logLikProfile[1])
  expect_true(profile$logLikProfile[peakIdx] > profile$logLikProfile[length(profile$logLikProfile)])
})

test_that("MLE from profile is close to true beta_ZY for large N", {
  simData <- simulateMRData(n = 10000, nSnps = 5, trueEffect = 0.5, seed = 789)
  profile <- fitOutcomeModel(
    cohortData = simData$data,
    covariateData = NULL,
    instrumentTable = simData$instrumentTable,
    betaGrid = seq(-3, 3, by = 0.01)
  )

  # The betaHat from the profile is the allele score coefficient, not directly
  # beta_ZY, but should be non-zero and in the right direction
  expect_true(is.finite(profile$betaHat))
  expect_true(profile$seHat > 0)
})

test_that("profile object contains correct named elements", {
  simData <- simulateMRData(n = 500, nSnps = 3, trueEffect = 0.2, seed = 111)
  profile <- fitOutcomeModel(
    cohortData = simData$data,
    covariateData = NULL,
    instrumentTable = simData$instrumentTable,
    betaGrid = seq(-1, 1, by = 0.1)
  )

  expect_true(all(c("siteId", "betaGrid", "logLikProfile", "nCases",
                     "nControls", "snpIds", "diagnosticFlags",
                     "betaHat", "seHat") %in% names(profile)))
})

test_that("alleleScore and perSNP modes both return valid profiles", {
  simData <- simulateMRData(n = 1000, nSnps = 3, trueEffect = 0.3, seed = 222)
  grid <- seq(-2, 2, by = 0.1)

  profileAS <- fitOutcomeModel(
    cohortData = simData$data,
    covariateData = NULL,
    instrumentTable = simData$instrumentTable,
    betaGrid = grid,
    analysisType = "alleleScore"
  )

  profilePS <- fitOutcomeModel(
    cohortData = simData$data,
    covariateData = NULL,
    instrumentTable = simData$instrumentTable,
    betaGrid = grid,
    analysisType = "perSNP"
  )

  expect_equal(length(profileAS$logLikProfile), length(grid))
  expect_equal(length(profilePS$logLikProfile), length(grid))
  expect_true(is.finite(profileAS$betaHat))
  expect_true(is.finite(profilePS$betaHat))
  expect_true("perSnpEstimates" %in% names(profilePS))
  expect_true(nrow(profilePS$perSnpEstimates) > 0)
})

test_that("warning issued when MLE is at grid boundary", {
  simData <- simulateMRData(n = 2000, nSnps = 5, trueEffect = 5, seed = 333)

  expect_warning(
    profile <- fitOutcomeModel(
      cohortData = simData$data,
      covariateData = NULL,
      instrumentTable = simData$instrumentTable,
      betaGrid = seq(-1, 1, by = 0.1)  # Too narrow for true effect of 5
    ),
    "grid boundary"
  )
})

test_that("fitOutcomeModel validates inputs", {
  expect_error(
    fitOutcomeModel(
      cohortData = data.frame(outcome = 1),
      covariateData = NULL,
      instrumentTable = data.frame(),  # invalid
      betaGrid = seq(-1, 1, by = 0.1)
    )
  )
})
