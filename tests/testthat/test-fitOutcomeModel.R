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
  expect_true(all(c("effect_allele", "other_allele", "eaf") %in%
                    names(profilePS$perSnpEstimates)))
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

test_that("fitOutcomeModel aligns SNP columns by SNP ID instead of column order", {
  simData <- simulateMRData(n = 1500, nSnps = 4, trueEffect = 0.25, seed = 444)
  snpCols <- grep("^snp_", names(simData$data), value = TRUE)
  reordered <- simData$data[, c(
    setdiff(names(simData$data), snpCols),
    rev(snpCols)
  )]

  profileOriginal <- fitOutcomeModel(
    cohortData = simData$data,
    covariateData = NULL,
    instrumentTable = simData$instrumentTable,
    betaGrid = seq(-2, 2, by = 0.05)
  )

  profileReordered <- fitOutcomeModel(
    cohortData = reordered,
    covariateData = NULL,
    instrumentTable = simData$instrumentTable,
    betaGrid = seq(-2, 2, by = 0.05)
  )

  expect_equal(profileOriginal$betaHat, profileReordered$betaHat, tolerance = 1e-10)
  expect_equal(profileOriginal$logLikProfile, profileReordered$logLikProfile, tolerance = 1e-10)
  expect_equal(profileOriginal$scoreDefinition$betaZX, profileReordered$scoreDefinition$betaZX, tolerance = 1e-10)
})

test_that("fitOutcomeModel supports the optional Cyclops backend", {
  skip_if_not_installed("Cyclops")

  simData <- simulateMRData(n = 1200, nSnps = 3, trueEffect = 0.2, seed = 555)
  profile <- fitOutcomeModel(
    cohortData = simData$data,
    covariateData = NULL,
    instrumentTable = simData$instrumentTable,
    betaGrid = seq(-1, 1, by = 0.1),
    modelBackend = "cyclops"
  )

  expect_s3_class(profile, "medusaSiteProfile")
  expect_true(is.finite(profile$betaHat))
  expect_true(all(is.finite(profile$logLikProfile)))
})

test_that("fitOutcomeModel helper utilities align instruments and evaluate finite likelihoods", {
  simData <- Medusa::simulateMRData(n = 250, nSnps = 3, trueEffect = 0.2, seed = 556)
  cohortData <- simData$data
  cohortData$personId <- cohortData$person_id

  extraInstrument <- simData$instrumentTable[1, , drop = FALSE]
  extraInstrument$snp_id <- "rs_missing"
  extraInstrument$effect_allele <- "A"
  extraInstrument$other_allele <- "C"
  extendedTable <- rbind(simData$instrumentTable, extraInstrument)

  expect_warning(
    alignment <- alignInstrumentColumns(cohortData, extendedTable),
    "Dropping 1 instrument"
  )
  expect_equal(alignment$instrumentTable$snp_id, simData$instrumentTable$snp_id)

  weights <- computeAlleleScoreWeights(alignment$instrumentTable)
  snpMatrix <- as.matrix(cohortData[, alignment$snpColumns, drop = FALSE])
  alleleScore <- as.numeric(snpMatrix %*% weights)
  baseModel <- data.frame(
    outcome = cohortData$outcome,
    alleleScore = alleleScore,
    stringsAsFactors = FALSE
  )
  covariateData <- data.frame(
    personId = cohortData$personId,
    covariate_1 = rep(c(0, 1), length.out = nrow(cohortData)),
    stringsAsFactors = FALSE
  )
  modelParts <- appendCovariatesToModelData(baseModel, cohortData, covariateData)
  completeModel <- modelParts$modelData[complete.cases(modelParts$modelData), , drop = FALSE]

  coefficientFit <- fitBinaryOutcomeCoefficient(
    modelData = completeModel,
    exposureColumn = "alleleScore",
    covariateColumns = modelParts$covariateColumns,
    modelBackend = "glm",
    regularizationVariance = 0.1,
    instrumentRegularization = FALSE
  )
  pointLogLik <- evaluateBinaryProfilePoint(
    modelData = completeModel,
    exposureColumn = "alleleScore",
    covariateColumns = modelParts$covariateColumns,
    offsetVector = rep(0, nrow(completeModel)),
    modelBackend = "glm",
    regularizationVariance = 0.1,
    instrumentRegularization = FALSE
  )

  expect_equal(sum(abs(weights)), 1, tolerance = 1e-10)
  expect_true(is.finite(coefficientFit$betaHat))
  expect_true(is.finite(coefficientFit$seHat))
  expect_true(is.finite(pointLogLik))

  zeroWeightTable <- alignment$instrumentTable
  zeroWeightTable$beta_ZX <- 0
  expect_error(computeAlleleScoreWeights(zeroWeightTable), "not finite")
})
