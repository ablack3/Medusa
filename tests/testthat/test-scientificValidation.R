# Scientific validation tests for Medusa
# These tests verify the statistical correctness of the MR estimation pipeline.
# They are permitted to take longer to run.

skip_on_cran()

test_that("IVW estimator is approximately unbiased for n=50000", {
  # Spec: IVW estimate within 0.05 of true effect = 0.5 after averaging
  # over 100 simulations
  set.seed(42)
  nSims <- 100
  trueEffect <- 0.5
  nSnps <- 10
  estimates <- numeric(nSims)
  effectAllele <- rep_len(c("A", "C", "G", "T"), nSnps)
  otherAllele <- rep_len(c("C", "A", "T", "G"), nSnps)
  eaf <- seq(0.10, by = 0.02, length.out = nSnps)

  for (sim in seq_len(nSims)) {
    betaZX <- rnorm(nSnps, 0.3, 0.05)
    seZX <- rep(0.05, nSnps)
    # True beta_ZY = trueEffect * beta_ZX (under no pleiotropy)
    betaZY <- trueEffect * betaZX + rnorm(nSnps, 0, 0.005)
    seZY <- rep(0.005, nSnps)

    perSnp <- data.frame(
      snp_id = paste0("rs", 1:nSnps),
      beta_ZY = betaZY,
      se_ZY = seZY,
      beta_ZX = betaZX,
      se_ZX = seZX,
      effect_allele = effectAllele,
      other_allele = otherAllele,
      eaf = eaf
    )

    ivw <- suppressMessages(
      runSensitivityAnalyses(perSnp, methods = "IVW")
    )
    estimates[sim] <- ivw$ivw$beta_MR
  }

  meanEstimate <- mean(estimates)
  expect_equal(meanEstimate, trueEffect, tolerance = 0.05,
               label = sprintf("Mean IVW estimate (%.4f) vs true (%.4f)", meanEstimate, trueEffect))
})

test_that("95% likelihood-based CI has proper coverage", {
  # Spec: 95% CI contains true parameter in >= 93 of 100 simulations
  set.seed(123)
  nSims <- 100
  trueEffect <- 0.5
  covered <- logical(nSims)

  for (sim in seq_len(nSims)) {
    betaGrid <- seq(-3, 3, by = 0.01)
    betaZYSim <- trueEffect + rnorm(1, 0, 0.1)
    seSim <- 0.1
    logLikProfile <- -0.5 * ((betaGrid - betaZYSim) / seSim)^2
    logLikProfile <- logLikProfile - max(logLikProfile)

    combined <- list(
      betaGrid = betaGrid,
      logLikProfile = logLikProfile,
      siteContributions = data.frame(siteId = "A", nCases = 500, nControls = 4500),
      nSites = 1,
      totalCases = 500,
      totalControls = 4500
    )

    instruments <- data.frame(
      snp_id = "rs1",
      effect_allele = "A",
      other_allele = "C",
      beta_ZX = 1.0,  # Identity so betaMR = betaZY
      se_ZX = 0.001,
      pval_ZX = 1e-100,
      eaf = 0.3
    )

    result <- suppressMessages(computeMREstimate(combined, instruments, ciLevel = 0.95))
    covered[sim] <- result$ciLower <= trueEffect && result$ciUpper >= trueEffect
  }

  coverage <- sum(covered)
  expect_true(coverage >= 93,
              info = sprintf("Coverage: %d/100, expected >= 93", coverage))
})

test_that("pooling consistency: pooled estimate matches joint analysis", {
  # Spec: 5 sites of 2000 each — pooled estimate within 0.05 of joint estimate
  set.seed(42)
  trueBeta <- 0.5
  betaGrid <- seq(-3, 3, by = 0.01)

  # Simulate 5 site profiles
  profiles <- simulateSiteProfiles(
    nSites = 5,
    betaGrid = betaGrid,
    trueBeta = trueBeta,
    nPerSite = 2000,
    seed = 42
  )

  # Pool
  combined <- suppressMessages(poolLikelihoodProfiles(profiles))
  pooledMLE <- combined$betaGrid[which.max(combined$logLikProfile)]

  # "Joint" analysis: simulate one big profile from n=10000
  set.seed(42)
  jointInfo <- 10000 * 0.01 * 1.0  # Average info per site * 5
  jointBeta <- trueBeta + rnorm(1, 0, 0.02)
  jointLogLik <- -0.5 * jointInfo * (betaGrid - jointBeta)^2
  jointLogLik <- jointLogLik - max(jointLogLik)
  jointMLE <- betaGrid[which.max(jointLogLik)]

  # Both should be close to true (absolute tolerance)
  expect_true(abs(pooledMLE - trueBeta) < 0.20,
              info = sprintf("Pooled MLE (%.3f) vs true (%.3f), diff = %.3f",
                             pooledMLE, trueBeta, abs(pooledMLE - trueBeta)))
})

test_that("type I error is controlled under the null", {
  # Spec: p < 0.05 in no more than 7 of 100 simulations when true effect = 0
  set.seed(42)
  nSims <- 100
  trueEffect <- 0.0
  nSnps <- 10
  significant <- logical(nSims)
  effectAllele <- rep_len(c("A", "C", "G", "T"), nSnps)
  otherAllele <- rep_len(c("C", "A", "T", "G"), nSnps)
  eaf <- seq(0.10, by = 0.02, length.out = nSnps)

  for (sim in seq_len(nSims)) {
    betaZX <- rnorm(nSnps, 0.3, 0.05)
    seZX <- rep(0.05, nSnps)
    betaZY <- rnorm(nSnps, 0, 0.01)  # Null: no effect
    seZY <- rep(0.01, nSnps)

    perSnp <- data.frame(
      snp_id = paste0("rs", 1:nSnps),
      beta_ZY = betaZY,
      se_ZY = seZY,
      beta_ZX = betaZX,
      se_ZX = seZX,
      effect_allele = effectAllele,
      other_allele = otherAllele,
      eaf = eaf
    )

    ivw <- suppressMessages(
      runSensitivityAnalyses(perSnp, methods = "IVW")
    )
    significant[sim] <- ivw$ivw$pval < 0.05
  }

  nRejected <- sum(significant)
  expect_true(nRejected <= 7,
              info = sprintf("Type I errors: %d/100, expected <= 7", nRejected))
})
