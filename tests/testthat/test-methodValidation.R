buildCleanPerSnpSummary <- function(nSnps = 10,
                                    effect = 0.25,
                                    intercept = 0.02,
                                    noiseSd = 0.004,
                                    seed = 9901) {
  set.seed(seed)
  betaZX <- runif(nSnps, min = 0.1, max = 0.5)
  allelePairs <- rep(
    list(c("A", "C"), c("C", "G"), c("G", "T"), c("T", "A")),
    length.out = nSnps
  )

  data.frame(
    snp_id = paste0("rs", seq_len(nSnps)),
    effect_allele = vapply(allelePairs, `[`, character(1), 1),
    other_allele = vapply(allelePairs, `[`, character(1), 2),
    eaf = seq(0.1, 0.7, length.out = nSnps),
    beta_ZX = betaZX,
    se_ZX = rep(0.03, nSnps),
    pval_ZX = 2 * stats::pnorm(-abs(betaZX / 0.03)),
    beta_ZY = intercept + effect * betaZX + stats::rnorm(nSnps, sd = noiseSd),
    se_ZY = rep(0.02, nSnps),
    pval_ZY = NA_real_,
    stringsAsFactors = FALSE
  )
}


buildAlleleScoreModelData <- function(cohortData,
                                      covariateData,
                                      instrumentTable) {
  alignment <- alignInstrumentColumns(cohortData, instrumentTable)
  weights <- computeAlleleScoreWeights(alignment$instrumentTable)
  snpMatrix <- as.matrix(cohortData[, alignment$snpColumns, drop = FALSE])
  snpMatrix[is.na(snpMatrix)] <- 0
  alleleScore <- as.numeric(snpMatrix %*% weights)

  modelParts <- appendCovariatesToModelData(
    modelData = data.frame(outcome = cohortData$outcome, alleleScore = alleleScore),
    cohortData = cohortData,
    covariateData = covariateData
  )
  modelData <- modelParts$modelData
  modelData <- modelData[complete.cases(modelData), , drop = FALSE]

  list(
    modelData = modelData,
    covariateColumns = modelParts$covariateColumns
  )
}


computeBruteForceGLMProfile <- function(modelData,
                                        covariateColumns,
                                        betaGrid) {
  profileFormula <- if (length(covariateColumns) > 0) {
    stats::as.formula(paste("~", paste(covariateColumns, collapse = " + ")))
  } else {
    stats::as.formula("~ 1")
  }
  designMatrix <- stats::model.matrix(profileFormula, data = modelData)
  response <- modelData$outcome
  exposure <- modelData$alleleScore

  unconstrainedFormula <- as.formula(paste(
    "outcome ~ alleleScore",
    if (length(covariateColumns) > 0) {
      paste("+", paste(covariateColumns, collapse = " + "))
    } else {
      ""
    }
  ))
  unconstrainedFit <- stats::glm(
    unconstrainedFormula,
    data = modelData,
    family = stats::binomial()
  )

  logLikProfile <- vapply(betaGrid, function(betaFixed) {
    fit <- stats::glm.fit(
      x = designMatrix,
      y = response,
      family = stats::binomial(),
      offset = betaFixed * exposure
    )
    fittedProb <- pmin(pmax(fit$fitted.values, 1e-12), 1 - 1e-12)
    sum(stats::dbinom(response, size = 1, prob = fittedProb, log = TRUE))
  }, numeric(1))

  peakIdx <- which.max(logLikProfile)
  ciMask <- logLikProfile >= max(logLikProfile) - stats::qchisq(0.95, df = 1) / 2
  ciBounds <- findProfileInterval(ciMask, peakIdx)

  list(
    logLikProfile = logLikProfile,
    betaHat = unname(stats::coef(unconstrainedFit)[["alleleScore"]]),
    seHat = unname(summary(unconstrainedFit)$coefficients["alleleScore", "Std. Error"]),
    ciLower = betaGrid[ciBounds$lower],
    ciUpper = betaGrid[ciBounds$upper]
  )
}


makeQuadraticSiteProfiles <- function(siteBetas,
                                      infos,
                                      betaGrid = seq(-2, 2, by = 0.02),
                                      sitePrefix = "site") {
  stopifnot(length(siteBetas) == length(infos))

  profiles <- lapply(seq_along(siteBetas), function(idx) {
    logLik <- -0.5 * infos[[idx]] * (betaGrid - siteBetas[[idx]])^2
    logLik <- logLik - max(logLik)

    list(
      siteId = sprintf("%s_%02d", sitePrefix, idx),
      betaGrid = betaGrid,
      logLikProfile = logLik,
      nCases = 100L + idx,
      nControls = 900L - idx,
      snpIds = "rs1"
    )
  })
  names(profiles) <- vapply(profiles, `[[`, character(1), "siteId")
  profiles
}


test_that("internal summary-data estimators agree with TwoSampleMR on clean harmonised data", {
  skip_if_not_installed("TwoSampleMR")

  perSnp <- buildCleanPerSnpSummary()

  internal <- suppressWarnings(
    suppressMessages(
      runSensitivityAnalyses(
        perSnp,
        methods = c("IVW", "MREgger", "WeightedMedian"),
        outcomeSampleSize = 10000,
        exposureSampleSize = 10000,
        outcomeType = "continuous",
        engine = "internal"
      )
    )
  )
  delegated <- suppressWarnings(
    suppressMessages(
      runSensitivityAnalyses(
        perSnp,
        methods = c("IVW", "MREgger", "WeightedMedian"),
        outcomeSampleSize = 10000,
        exposureSampleSize = 10000,
        outcomeType = "continuous",
        engine = "TwoSampleMR"
      )
    )
  )

  expect_lt(abs(internal$ivw$beta_MR - delegated$ivw$beta_MR), 0.01)
  expect_lt(abs(internal$weightedMedian$beta_MR - delegated$weightedMedian$beta_MR), 0.02)
  expect_lt(abs(internal$mrEgger$beta_MR - delegated$mrEgger$beta_MR), 0.02)
  expect_lt(abs(internal$mrEgger$intercept - delegated$mrEgger$intercept), 0.005)
})

test_that("fitOutcomeModel matches a brute-force constrained glm profile on a fixed dataset", {
  simData <- simulateMRData(n = 600, nSnps = 4, trueEffect = 0.2, seed = 9100)
  covariateData <- simData$data[, c("person_id", "confounder_1", "confounder_2"), drop = FALSE]
  betaGrid <- seq(-1.2, 1.2, by = 0.1)

  profile <- suppressWarnings(
    suppressMessages(
      fitOutcomeModel(
        cohortData = simData$data,
        covariateData = covariateData,
        instrumentTable = simData$instrumentTable,
        betaGrid = betaGrid,
        modelBackend = "glm"
      )
    )
  )

  modelParts <- buildAlleleScoreModelData(
    cohortData = simData$data,
    covariateData = covariateData,
    instrumentTable = simData$instrumentTable
  )
  brute <- computeBruteForceGLMProfile(
    modelData = modelParts$modelData,
    covariateColumns = modelParts$covariateColumns,
    betaGrid = betaGrid
  )

  profileMask <- profile$logLikProfile >= max(profile$logLikProfile) - stats::qchisq(0.95, 1) / 2
  profileBounds <- findProfileInterval(profileMask, which.max(profile$logLikProfile))

  expect_equal(profile$logLikProfile, brute$logLikProfile, tolerance = 1e-10)
  expect_equal(profile$betaHat, brute$betaHat, tolerance = 1e-10)
  expect_equal(profile$seHat, brute$seHat, tolerance = 1e-10)
  expect_equal(betaGrid[profileBounds$lower], brute$ciLower, tolerance = 1e-12)
  expect_equal(betaGrid[profileBounds$upper], brute$ciUpper, tolerance = 1e-12)
})

test_that("glm and unpenalized Cyclops backends agree on low-dimensional fits", {
  skip_if_not_installed("Cyclops")

  simData <- simulateMRData(n = 800, nSnps = 4, trueEffect = 0.2, seed = 9200)
  covariateData <- simData$data[, c("person_id", "confounder_1", "confounder_2"), drop = FALSE]
  betaGrid <- seq(-1.2, 1.2, by = 0.1)

  glmProfile <- suppressWarnings(
    suppressMessages(
      fitOutcomeModel(
        cohortData = simData$data,
        covariateData = covariateData,
        instrumentTable = simData$instrumentTable,
        betaGrid = betaGrid,
        modelBackend = "glm"
      )
    )
  )
  cyclopsProfile <- suppressWarnings(
    suppressMessages(
      fitOutcomeModel(
        cohortData = simData$data,
        covariateData = covariateData,
        instrumentTable = simData$instrumentTable,
        betaGrid = betaGrid,
        modelBackend = "cyclops",
        regularizationVariance = Inf
      )
    )
  )
  penalizedProfile <- suppressWarnings(
    suppressMessages(
      fitOutcomeModel(
        cohortData = simData$data,
        covariateData = covariateData,
        instrumentTable = simData$instrumentTable,
        betaGrid = betaGrid,
        modelBackend = "cyclops",
        regularizationVariance = 0.05,
        instrumentRegularization = TRUE
      )
    )
  )

  expect_equal(glmProfile$betaHat, cyclopsProfile$betaHat, tolerance = 1e-3)
  expect_equal(glmProfile$seHat, cyclopsProfile$seHat, tolerance = 1e-3)
  expect_equal(glmProfile$logLikProfile, cyclopsProfile$logLikProfile, tolerance = 1e-6)
  expect_gt(abs(penalizedProfile$betaHat - glmProfile$betaHat), 1e-4)
})

test_that("federated pooling recovers truth more closely with more sites and reduces uncertainty", {
  betaGrid <- seq(-2, 2, by = 0.01)
  instrumentTable <- data.frame(
    snp_id = "rs1",
    effect_allele = "A",
    other_allele = "C",
    beta_ZX = 1,
    se_ZX = 0.01,
    pval_ZX = 1e-30,
    eaf = 0.3,
    stringsAsFactors = FALSE
  )

  profiles3 <- makeQuadraticSiteProfiles(
    siteBetas = c(-0.10, 0.30, 0.30),
    infos = c(20, 18, 22),
    betaGrid = betaGrid,
    sitePrefix = "three"
  )
  profiles10 <- makeQuadraticSiteProfiles(
    siteBetas = c(-0.10, 0.30, 0.30, 0.18, 0.21, 0.19, 0.22, 0.20, 0.17, 0.23),
    infos = c(20, 18, 22, 21, 19, 20, 22, 21, 19, 20),
    betaGrid = betaGrid,
    sitePrefix = "ten"
  )

  pooled3 <- suppressMessages(poolLikelihoodProfiles(profiles3))
  pooled10 <- suppressMessages(poolLikelihoodProfiles(profiles10))
  estimate3 <- suppressMessages(computeMREstimate(pooled3, instrumentTable))
  estimate10 <- suppressMessages(computeMREstimate(pooled10, instrumentTable))

  singleEstimates10 <- lapply(profiles10, function(profile) {
    singleCombined <- suppressMessages(poolLikelihoodProfiles(list(profile)))
    suppressMessages(computeMREstimate(singleCombined, instrumentTable))
  })
  medianSingleSe <- stats::median(vapply(singleEstimates10, `[[`, numeric(1), "seMR"))

  expect_lt(abs(estimate10$betaMR - 0.2), abs(estimate3$betaMR - 0.2))
  expect_lt(estimate10$seMR, medianSingleSe)

  reversedEstimate10 <- suppressMessages(
    computeMREstimate(
      suppressMessages(poolLikelihoodProfiles(rev(profiles10))),
      instrumentTable
    )
  )
  expect_equal(estimate10$betaMR, reversedEstimate10$betaMR, tolerance = 1e-12)
  expect_equal(estimate10$seMR, reversedEstimate10$seMR, tolerance = 1e-12)
})

test_that("federated pooling is invariant to equivalent site partitioning", {
  betaGrid <- seq(-2, 2, by = 0.02)
  instrumentTable <- data.frame(
    snp_id = "rs1",
    effect_allele = "A",
    other_allele = "C",
    beta_ZX = 1,
    se_ZX = 0.01,
    pval_ZX = 1e-30,
    eaf = 0.3,
    stringsAsFactors = FALSE
  )

  coarseProfiles <- makeQuadraticSiteProfiles(
    siteBetas = c(0.1, 0.3),
    infos = c(40, 36),
    betaGrid = betaGrid,
    sitePrefix = "coarse"
  )
  splitProfiles <- makeQuadraticSiteProfiles(
    siteBetas = c(0.1, 0.1, 0.3, 0.3),
    infos = c(20, 20, 18, 18),
    betaGrid = betaGrid,
    sitePrefix = "split"
  )

  coarseEstimate <- suppressMessages(
    computeMREstimate(
      suppressMessages(poolLikelihoodProfiles(coarseProfiles)),
      instrumentTable
    )
  )
  splitEstimate <- suppressMessages(
    computeMREstimate(
      suppressMessages(poolLikelihoodProfiles(splitProfiles)),
      instrumentTable
    )
  )

  expect_equal(coarseEstimate$betaMR, splitEstimate$betaMR, tolerance = 1e-12)
  expect_equal(coarseEstimate$seMR, splitEstimate$seMR, tolerance = 1e-12)
  expect_equal(
    coarseEstimate$combinedProfile$logLikProfile,
    splitEstimate$combinedProfile$logLikProfile,
    tolerance = 1e-12
  )
})
