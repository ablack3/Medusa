test_that("Wald ratio formula: beta_MR = beta_ZY / beta_ZX", {
  betaGrid <- seq(-3, 3, by = 0.01)
  trueBetaZY <- 0.5
  se <- 0.1
  logLikProfile <- -0.5 * ((betaGrid - trueBetaZY) / se)^2
  logLikProfile <- logLikProfile - max(logLikProfile)

  combined <- list(
    betaGrid = betaGrid,
    logLikProfile = logLikProfile,
    siteContributions = data.frame(siteId = "A", nCases = 100, nControls = 900),
    nSites = 1,
    totalCases = 100,
    totalControls = 900
  )
  class(combined) <- "medusaCombinedProfile"

  instruments <- data.frame(
    snp_id = "rs1",
    effect_allele = "A",
    other_allele = "C",
    beta_ZX = 0.25,
    se_ZX = 0.02,
    pval_ZX = 1e-10,
    eaf = 0.3,
    stringsAsFactors = FALSE
  )

  result <- computeMREstimate(combined, instruments)

  # beta_MR should equal beta_ZY / beta_ZX = 0.5 / 0.25 = 2.0
  expect_equal(result$betaMR, trueBetaZY / 0.25, tolerance = 0.05)
  expect_equal(result$betaZY, trueBetaZY, tolerance = 0.02)
})

test_that("delta method standard error formula", {
  betaGrid <- seq(-3, 3, by = 0.01)
  betaZY <- 0.5
  seZY <- 0.1
  logLikProfile <- -0.5 * ((betaGrid - betaZY) / seZY)^2
  logLikProfile <- logLikProfile - max(logLikProfile)

  combined <- list(
    betaGrid = betaGrid,
    logLikProfile = logLikProfile,
    siteContributions = data.frame(siteId = "A", nCases = 100, nControls = 900),
    nSites = 1,
    totalCases = 100,
    totalControls = 900
  )

  betaZX <- 0.25
  seZX <- 0.02
  instruments <- data.frame(
    snp_id = "rs1",
    effect_allele = "A",
    other_allele = "C",
    beta_ZX = betaZX,
    se_ZX = seZX,
    pval_ZX = 1e-10,
    eaf = 0.3,
    stringsAsFactors = FALSE
  )

  result <- computeMREstimate(combined, instruments)

  # Expected delta method SE
  expectedSE <- sqrt((seZY / betaZX)^2 + (betaZY * seZX / betaZX^2)^2)
  expect_equal(result$seMR, expectedSE, tolerance = 0.1 * expectedSE)
})

test_that("95% CI from likelihood profile contains true parameter in simulation", {
  set.seed(42)
  nReps <- 50
  trueEffect <- 0.5
  covered <- logical(nReps)

  for (rep in seq_len(nReps)) {
    betaGrid <- seq(-3, 3, by = 0.01)
    # Simulate a betaZY estimate with some noise
    betaZYSim <- trueEffect + rnorm(1, 0, 0.15)
    seSim <- 0.15
    logLikProfile <- -0.5 * ((betaGrid - betaZYSim) / seSim)^2
    logLikProfile <- logLikProfile - max(logLikProfile)

    combined <- list(
      betaGrid = betaGrid,
      logLikProfile = logLikProfile,
      siteContributions = data.frame(siteId = "A", nCases = 200, nControls = 800),
      nSites = 1,
      totalCases = 200,
      totalControls = 800
    )

    instruments <- data.frame(
      snp_id = "rs1",
      effect_allele = "A",
      other_allele = "C",
      beta_ZX = 1.0,  # Use 1.0 so betaMR = betaZY
      se_ZX = 0.01,
      pval_ZX = 1e-50,
      eaf = 0.3,
      stringsAsFactors = FALSE
    )

    result <- suppressMessages(computeMREstimate(combined, instruments))
    covered[rep] <- result$ciLower <= trueEffect && result$ciUpper >= trueEffect
  }

  # Should cover at least 85% of the time (with some slack for simulation noise)
  coverage <- mean(covered)
  expect_true(coverage >= 0.80,
              info = sprintf("Coverage was %.1f%%, expected >= 80%%", coverage * 100))
})

test_that("CI level parameter correctly changes threshold", {
  betaGrid <- seq(-3, 3, by = 0.01)
  logLikProfile <- -0.5 * (betaGrid / 0.2)^2  # peak at 0
  logLikProfile <- logLikProfile - max(logLikProfile)

  combined <- list(
    betaGrid = betaGrid,
    logLikProfile = logLikProfile,
    siteContributions = data.frame(siteId = "A", nCases = 100, nControls = 900),
    nSites = 1,
    totalCases = 100,
    totalControls = 900
  )

  instruments <- data.frame(
    snp_id = "rs1",
    effect_allele = "A",
    other_allele = "C",
    beta_ZX = 1.0,
    se_ZX = 0.01,
    pval_ZX = 1e-50,
    eaf = 0.3,
    stringsAsFactors = FALSE
  )

  result90 <- suppressMessages(computeMREstimate(combined, instruments, ciLevel = 0.90))
  result99 <- suppressMessages(computeMREstimate(combined, instruments, ciLevel = 0.99))

  # 99% CI should be wider than 90% CI
  width90 <- result90$ciUpper - result90$ciLower
  width99 <- result99$ciUpper - result99$ciLower
  expect_true(width99 > width90)
})

test_that("odds ratio is computed correctly", {
  betaGrid <- seq(-3, 3, by = 0.01)
  logLikProfile <- -0.5 * ((betaGrid - 0.5) / 0.1)^2
  logLikProfile <- logLikProfile - max(logLikProfile)

  combined <- list(
    betaGrid = betaGrid,
    logLikProfile = logLikProfile,
    siteContributions = data.frame(siteId = "A", nCases = 100, nControls = 900),
    nSites = 1,
    totalCases = 100,
    totalControls = 900
  )

  instruments <- data.frame(
    snp_id = "rs1",
    effect_allele = "A",
    other_allele = "C",
    beta_ZX = 1.0,
    se_ZX = 0.01,
    pval_ZX = 1e-50,
    eaf = 0.3,
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(computeMREstimate(combined, instruments))

  expect_equal(result$oddsRatio, exp(result$betaMR))
  expect_equal(result$orCiLower, exp(result$ciLower))
  expect_equal(result$orCiUpper, exp(result$ciUpper))
})

test_that("warning issued for profile with multiple local maxima", {
  betaGrid <- seq(-3, 3, by = 0.01)
  # Create bimodal profile
  logLikProfile <- pmax(
    -0.5 * ((betaGrid - (-1)) / 0.2)^2,
    -0.5 * ((betaGrid - 1) / 0.2)^2
  )
  logLikProfile <- logLikProfile - max(logLikProfile)

  combined <- list(
    betaGrid = betaGrid,
    logLikProfile = logLikProfile,
    siteContributions = data.frame(siteId = "A", nCases = 100, nControls = 900),
    nSites = 1,
    totalCases = 100,
    totalControls = 900
  )

  instruments <- data.frame(
    snp_id = "rs1",
    effect_allele = "A",
    other_allele = "C",
    beta_ZX = 1.0,
    se_ZX = 0.01,
    pval_ZX = 1e-50,
    eaf = 0.3,
    stringsAsFactors = FALSE
  )

  expect_warning(
    computeMREstimate(combined, instruments),
    "multiple local maxima"
  )
})

test_that("multi-SNP estimate uses the fitted allele-score denominator", {
  betaGrid <- seq(-3, 3, by = 0.01)
  trueBetaZY <- 0.6
  se <- 0.1
  logLikProfile <- -0.5 * ((betaGrid - trueBetaZY) / se)^2
  logLikProfile <- logLikProfile - max(logLikProfile)

  combined <- list(
    betaGrid = betaGrid,
    logLikProfile = logLikProfile,
    siteContributions = data.frame(siteId = "A", nCases = 100, nControls = 900),
    nSites = 1,
    totalCases = 100,
    totalControls = 900
  )

  instruments <- data.frame(
    snp_id = c("rs1", "rs2"),
    effect_allele = c("A", "G"),
    other_allele = c("C", "T"),
    beta_ZX = c(0.2, 0.4),
    se_ZX = c(0.1, 0.1),
    pval_ZX = c(1e-20, 1e-30),
    eaf = c(0.3, 0.4),
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(computeMREstimate(combined, instruments))
  expectedWeights <- c(0.3333333, 0.6666667)
  expectedBetaZX <- sum(expectedWeights * instruments$beta_ZX)

  expect_equal(result$betaZX, expectedBetaZX, tolerance = 1e-6)
  expect_equal(result$betaMR, trueBetaZY / expectedBetaZX, tolerance = 0.05)
})

test_that("computeMREstimate falls back to Wald CI when no connected profile interval is found", {
  betaGrid <- seq(-2, 2, by = 0.01)
  betaZY <- 0.4
  seZY <- 0.15
  logLikProfile <- -0.5 * ((betaGrid - betaZY) / seZY)^2
  logLikProfile <- logLikProfile - max(logLikProfile)

  combined <- list(
    betaGrid = betaGrid,
    logLikProfile = logLikProfile,
    siteContributions = data.frame(siteId = "A", nCases = 100, nControls = 900),
    nSites = 1,
    totalCases = 100,
    totalControls = 900
  )

  instruments <- data.frame(
    snp_id = "rs1",
    effect_allele = "A",
    other_allele = "C",
    beta_ZX = 1.0,
    se_ZX = 0.01,
    pval_ZX = 1e-20,
    eaf = 0.2,
    stringsAsFactors = FALSE
  )

  computeMREstimateStub <- computeMREstimate
  mockery::stub(computeMREstimateStub, "findProfileInterval", function(ciMask, peakIdx) NULL)

  result <- suppressMessages(computeMREstimateStub(combined, instruments))
  expectedHalfWidth <- qnorm(0.975) * result$seZY

  expect_equal(result$ciLower, result$betaMR - expectedHalfWidth, tolerance = 0.02)
  expect_equal(result$ciUpper, result$betaMR + expectedHalfWidth, tolerance = 0.02)
})

test_that("computeMREstimate errors when the stored score denominator is effectively zero", {
  betaGrid <- seq(-1, 1, by = 0.01)
  logLikProfile <- -0.5 * (betaGrid / 0.1)^2
  logLikProfile <- logLikProfile - max(logLikProfile)

  combined <- list(
    betaGrid = betaGrid,
    logLikProfile = logLikProfile,
    siteContributions = data.frame(siteId = "A", nCases = 100, nControls = 900),
    nSites = 1,
    totalCases = 100,
    totalControls = 900,
    scoreDefinition = list(
      snp_id = "rs1",
      weights = 1,
      betaZX = 0,
      seZX = 0.01
    )
  )

  instruments <- data.frame(
    snp_id = "rs1",
    effect_allele = "A",
    other_allele = "C",
    beta_ZX = 0.2,
    se_ZX = 0.02,
    pval_ZX = 1e-20,
    eaf = 0.2,
    stringsAsFactors = FALSE
  )

  expect_error(
    suppressMessages(computeMREstimate(combined, instruments)),
    "too close to zero"
  )
})

test_that("computeMREstimate reorders confidence limits when the score denominator is negative", {
  betaGrid <- seq(-2, 2, by = 0.01)
  betaZY <- 0.3
  seZY <- 0.1
  logLikProfile <- -0.5 * ((betaGrid - betaZY) / seZY)^2
  logLikProfile <- logLikProfile - max(logLikProfile)

  combined <- list(
    betaGrid = betaGrid,
    logLikProfile = logLikProfile,
    siteContributions = data.frame(siteId = "A", nCases = 100, nControls = 900),
    nSites = 1,
    totalCases = 100,
    totalControls = 900,
    scoreDefinition = list(
      snp_id = "rs1",
      weights = 1,
      betaZX = -0.5,
      seZX = 0.02
    )
  )

  instruments <- data.frame(
    snp_id = "rs1",
    effect_allele = "A",
    other_allele = "C",
    beta_ZX = 0.2,
    se_ZX = 0.02,
    pval_ZX = 1e-20,
    eaf = 0.2,
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(computeMREstimate(combined, instruments))

  expect_lt(result$ciLower, result$ciUpper)
  expect_lt(result$betaZX, 0)
})

test_that("findProfileInterval returns NULL when the peak is outside the CI mask", {
  expect_null(findProfileInterval(c(FALSE, TRUE, TRUE), peakIdx = 1))
})

test_that("computeMREstimate is equivariant to rescaling the outcome profile", {
  baseGrid <- seq(-2, 2, by = 0.01)
  baseBetaZY <- 0.4
  baseSeZY <- 0.1
  scaledGrid <- seq(-4, 4, by = 0.02)
  scaleFactor <- 2

  baseProfile <- -0.5 * ((baseGrid - baseBetaZY) / baseSeZY)^2
  baseProfile <- baseProfile - max(baseProfile)
  scaledProfile <- -0.5 * ((scaledGrid - scaleFactor * baseBetaZY) /
                             (scaleFactor * baseSeZY))^2
  scaledProfile <- scaledProfile - max(scaledProfile)

  combinedBase <- list(
    betaGrid = baseGrid,
    logLikProfile = baseProfile,
    siteContributions = data.frame(siteId = "A", nCases = 100, nControls = 900),
    nSites = 1,
    totalCases = 100,
    totalControls = 900
  )
  combinedScaled <- list(
    betaGrid = scaledGrid,
    logLikProfile = scaledProfile,
    siteContributions = data.frame(siteId = "A", nCases = 100, nControls = 900),
    nSites = 1,
    totalCases = 100,
    totalControls = 900
  )

  instruments <- data.frame(
    snp_id = "rs1",
    effect_allele = "A",
    other_allele = "C",
    beta_ZX = 0.5,
    se_ZX = 0.02,
    pval_ZX = 1e-20,
    eaf = 0.3,
    stringsAsFactors = FALSE
  )

  baseResult <- suppressMessages(computeMREstimate(combinedBase, instruments))
  scaledResult <- suppressMessages(computeMREstimate(combinedScaled, instruments))

  expect_equal(scaledResult$betaMR, scaleFactor * baseResult$betaMR, tolerance = 0.05)
  expect_equal(scaledResult$seMR, scaleFactor * baseResult$seMR, tolerance = 0.05)
  expect_equal(
    scaledResult$ciUpper - scaledResult$ciLower,
    scaleFactor * (baseResult$ciUpper - baseResult$ciLower),
    tolerance = 0.1
  )
})

test_that("computeMREstimate returns tighter intervals for sharper profiles", {
  betaGrid <- seq(-2, 2, by = 0.01)
  betaZY <- 0.3
  broadProfile <- -0.5 * ((betaGrid - betaZY) / 0.2)^2
  broadProfile <- broadProfile - max(broadProfile)
  sharpProfile <- -0.5 * ((betaGrid - betaZY) / 0.08)^2
  sharpProfile <- sharpProfile - max(sharpProfile)

  combinedBroad <- list(
    betaGrid = betaGrid,
    logLikProfile = broadProfile,
    siteContributions = data.frame(siteId = "A", nCases = 100, nControls = 900),
    nSites = 1,
    totalCases = 100,
    totalControls = 900
  )
  combinedSharp <- list(
    betaGrid = betaGrid,
    logLikProfile = sharpProfile,
    siteContributions = data.frame(siteId = "A", nCases = 100, nControls = 900),
    nSites = 1,
    totalCases = 100,
    totalControls = 900
  )

  instruments <- data.frame(
    snp_id = "rs1",
    effect_allele = "A",
    other_allele = "C",
    beta_ZX = 0.5,
    se_ZX = 0.02,
    pval_ZX = 1e-20,
    eaf = 0.3,
    stringsAsFactors = FALSE
  )

  broadResult <- suppressMessages(computeMREstimate(combinedBroad, instruments))
  sharpResult <- suppressMessages(computeMREstimate(combinedSharp, instruments))

  expect_lt(sharpResult$ciUpper - sharpResult$ciLower,
            broadResult$ciUpper - broadResult$ciLower)
  expect_lt(sharpResult$seMR, broadResult$seMR)
})
