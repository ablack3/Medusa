test_that("IVW estimate matches manual weighted regression calculation", {
  set.seed(42)
  nSnps <- 10
  betaZX <- rnorm(nSnps, 0.3, 0.05)
  trueEffect <- 0.5
  betaZY <- trueEffect * betaZX + rnorm(nSnps, 0, 0.01)
  seZY <- rep(0.02, nSnps)

  perSnp <- data.frame(
    snp_id = paste0("rs", 1:nSnps),
    beta_ZY = betaZY,
    se_ZY = seZY,
    beta_ZX = betaZX,
    se_ZX = rep(0.05, nSnps),
    stringsAsFactors = FALSE
  )

  results <- suppressMessages(runSensitivityAnalyses(perSnp, methods = "IVW"))

  # Manual IVW calculation
  weights <- 1 / seZY^2
  manualBeta <- sum(weights * betaZY * betaZX) / sum(weights * betaZX^2)
  manualSE <- sqrt(1 / sum(weights * betaZX^2))

  expect_equal(results$ivw$beta_MR, manualBeta, tolerance = 1e-10)
  expect_equal(results$ivw$se_MR, manualSE, tolerance = 1e-10)
})

test_that("MR-Egger returns intercept estimate and its p-value", {
  set.seed(42)
  nSnps <- 10
  betaZX <- rnorm(nSnps, 0.3, 0.05)
  betaZY <- 0.5 * betaZX + 0.02 + rnorm(nSnps, 0, 0.01)  # Non-zero intercept

  perSnp <- data.frame(
    snp_id = paste0("rs", 1:nSnps),
    beta_ZY = betaZY,
    se_ZY = rep(0.02, nSnps),
    beta_ZX = betaZX,
    se_ZX = rep(0.05, nSnps),
    stringsAsFactors = FALSE
  )

  results <- suppressMessages(runSensitivityAnalyses(perSnp, methods = "MREgger"))

  expect_true("intercept" %in% names(results$mrEgger))
  expect_true("intercept_se" %in% names(results$mrEgger))
  expect_true("intercept_pval" %in% names(results$mrEgger))
  expect_true(is.numeric(results$mrEgger$intercept))
})

test_that("leave-one-out returns K estimates for K SNPs", {
  nSnps <- 8
  perSnp <- data.frame(
    snp_id = paste0("rs", 1:nSnps),
    beta_ZY = rnorm(nSnps, 0.15, 0.02),
    se_ZY = rep(0.02, nSnps),
    beta_ZX = rnorm(nSnps, 0.3, 0.05),
    se_ZX = rep(0.05, nSnps),
    stringsAsFactors = FALSE
  )

  results <- suppressMessages(
    runSensitivityAnalyses(perSnp, methods = "LeaveOneOut")
  )

  expect_equal(nrow(results$leaveOneOut), nSnps)
  expect_true(all(perSnp$snp_id %in% results$leaveOneOut$snp_removed))
})

test_that("Steiger filtering removes SNPs with wrong variance direction", {
  perSnp <- data.frame(
    snp_id = c("rs1", "rs2", "rs3", "rs4"),
    beta_ZY = c(0.01, 0.5, 0.02, 0.6),  # rs2, rs4 have large outcome effect
    se_ZY = rep(0.02, 4),
    beta_ZX = c(0.3, 0.01, 0.25, 0.01),  # rs2, rs4 have small exposure effect
    se_ZX = rep(0.05, 4),
    stringsAsFactors = FALSE
  )

  results <- suppressMessages(
    runSensitivityAnalyses(perSnp, methods = "Steiger")
  )

  # SNPs where r2_outcome > r2_exposure should be removed
  expect_true(results$steiger$n_removed > 0)
  expect_true(results$steiger$n_remaining < 4)
})

test_that("all methods handle single-SNP case gracefully", {
  perSnp <- data.frame(
    snp_id = "rs1",
    beta_ZY = 0.15,
    se_ZY = 0.02,
    beta_ZX = 0.3,
    se_ZX = 0.05,
    stringsAsFactors = FALSE
  )

  results <- suppressMessages(
    runSensitivityAnalyses(perSnp)
  )

  # IVW should work with single SNP
  expect_true(is.finite(results$ivw$beta_MR))
  # MR-Egger and Weighted Median need >= 3 SNPs
  expect_null(results$mrEgger)
  expect_null(results$weightedMedian)
  expect_null(results$leaveOneOut)
})

test_that("summary table contains all computed methods", {
  set.seed(42)
  nSnps <- 10
  perSnp <- data.frame(
    snp_id = paste0("rs", 1:nSnps),
    beta_ZY = rnorm(nSnps, 0.15, 0.02),
    se_ZY = rep(0.02, nSnps),
    beta_ZX = rnorm(nSnps, 0.3, 0.05),
    se_ZX = rep(0.05, nSnps),
    stringsAsFactors = FALSE
  )

  results <- suppressMessages(runSensitivityAnalyses(perSnp))

  expect_true("summary" %in% names(results))
  expect_true("IVW" %in% results$summary$method)
  expect_true("MR-Egger" %in% results$summary$method)
  expect_true("Weighted Median" %in% results$summary$method)
})

test_that("IVW p-value is computed correctly", {
  perSnp <- data.frame(
    snp_id = c("rs1", "rs2", "rs3"),
    beta_ZY = c(0, 0, 0),  # Null effect
    se_ZY = c(0.1, 0.1, 0.1),
    beta_ZX = c(0.3, 0.3, 0.3),
    se_ZX = c(0.05, 0.05, 0.05),
    stringsAsFactors = FALSE
  )

  results <- suppressMessages(runSensitivityAnalyses(perSnp, methods = "IVW"))

  # With true zero effect, beta_MR should be ~0 and pval should be non-significant
  expect_equal(results$ivw$beta_MR, 0, tolerance = 1e-10)
  expect_true(results$ivw$pval > 0.5)  # Should be 1.0 for exact zero
})

test_that("all SNPs failing Steiger produces warning", {
  perSnp <- data.frame(
    snp_id = c("rs1", "rs2", "rs3"),
    beta_ZY = c(0.5, 0.4, 0.6),  # Large outcome effect
    se_ZY = rep(0.02, 3),
    beta_ZX = c(0.01, 0.01, 0.01),  # Tiny exposure effect
    se_ZX = rep(0.05, 3),
    stringsAsFactors = FALSE
  )

  expect_warning(
    results <- suppressMessages(
      runSensitivityAnalyses(perSnp, methods = "Steiger")
    ),
    "All SNPs failed Steiger filter"
  )
  expect_true(is.na(results$steiger$beta_MR))
})
