with_harmonisation_columns <- function(df) {
  n <- nrow(df)
  allelePairs <- rep(
    list(c("A", "C"), c("C", "G"), c("G", "T"), c("T", "A")),
    length.out = n
  )
  df$effect_allele <- vapply(allelePairs, `[`, character(1), 1)
  df$other_allele <- vapply(allelePairs, `[`, character(1), 2)
  df$eaf <- seq(0.1, 0.8, length.out = n)
  df
}

test_that("IVW estimate matches manual weighted regression calculation", {
  set.seed(42)
  nSnps <- 10
  betaZX <- rnorm(nSnps, 0.3, 0.05)
  trueEffect <- 0.5
  betaZY <- trueEffect * betaZX + rnorm(nSnps, 0, 0.01)
  seZY <- rep(0.02, nSnps)

  perSnp <- with_harmonisation_columns(data.frame(
    snp_id = paste0("rs", 1:nSnps),
    beta_ZY = betaZY,
    se_ZY = seZY,
    beta_ZX = betaZX,
    se_ZX = rep(0.05, nSnps),
    stringsAsFactors = FALSE
  ))

  results <- suppressMessages(
    runSensitivityAnalyses(perSnp, methods = "IVW", engine = "internal")
  )

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
  betaZY <- 0.5 * betaZX + 0.02 + rnorm(nSnps, 0, 0.01)

  perSnp <- with_harmonisation_columns(data.frame(
    snp_id = paste0("rs", 1:nSnps),
    beta_ZY = betaZY,
    se_ZY = rep(0.02, nSnps),
    beta_ZX = betaZX,
    se_ZX = rep(0.05, nSnps),
    stringsAsFactors = FALSE
  ))

  results <- suppressMessages(
    runSensitivityAnalyses(perSnp, methods = "MREgger", engine = "internal")
  )

  expect_true("intercept" %in% names(results$mrEgger))
  expect_true("intercept_se" %in% names(results$mrEgger))
  expect_true("intercept_pval" %in% names(results$mrEgger))
  expect_true(is.numeric(results$mrEgger$intercept))
})

test_that("leave-one-out returns K estimates for K SNPs", {
  nSnps <- 8
  perSnp <- with_harmonisation_columns(data.frame(
    snp_id = paste0("rs", 1:nSnps),
    beta_ZY = rnorm(nSnps, 0.15, 0.02),
    se_ZY = rep(0.02, nSnps),
    beta_ZX = rnorm(nSnps, 0.3, 0.05),
    se_ZX = rep(0.05, nSnps),
    stringsAsFactors = FALSE
  ))

  results <- suppressMessages(
    runSensitivityAnalyses(perSnp, methods = "LeaveOneOut", engine = "internal")
  )

  expect_equal(nrow(results$leaveOneOut), nSnps)
  expect_true(all(perSnp$snp_id %in% results$leaveOneOut$snp_removed))
})

test_that("Steiger filtering removes SNPs with wrong variance direction", {
  perSnp <- with_harmonisation_columns(data.frame(
    snp_id = c("rs1", "rs2", "rs3", "rs4"),
    beta_ZY = c(0.01, 0.5, 0.02, 0.6),
    se_ZY = rep(0.02, 4),
    beta_ZX = c(0.3, 0.01, 0.25, 0.01),
    se_ZX = rep(0.05, 4),
    stringsAsFactors = FALSE
  ))

  results <- suppressMessages(
    runSensitivityAnalyses(
      perSnp,
      methods = "Steiger",
      outcomeSampleSize = 10000,
      exposureSampleSize = 10000,
      outcomeType = "continuous",
      engine = "internal"
    )
  )

  expect_true(results$steiger$n_removed > 0)
  expect_true(results$steiger$n_remaining < 4)
})

test_that("all methods handle single-SNP case gracefully", {
  perSnp <- with_harmonisation_columns(data.frame(
    snp_id = "rs1",
    beta_ZY = 0.15,
    se_ZY = 0.02,
    beta_ZX = 0.3,
    se_ZX = 0.05,
    stringsAsFactors = FALSE
  ))

  results <- suppressMessages(
    runSensitivityAnalyses(
      perSnp,
      methods = c("IVW", "MREgger", "WeightedMedian", "LeaveOneOut"),
      engine = "internal"
    )
  )

  expect_true(is.finite(results$ivw$beta_MR))
  expect_null(results$mrEgger)
  expect_null(results$weightedMedian)
  expect_null(results$leaveOneOut)
})

test_that("summary table contains all computed methods", {
  set.seed(42)
  nSnps <- 10
  perSnp <- with_harmonisation_columns(data.frame(
    snp_id = paste0("rs", 1:nSnps),
    beta_ZY = rnorm(nSnps, 0.15, 0.02),
    se_ZY = rep(0.02, nSnps),
    beta_ZX = rnorm(nSnps, 0.3, 0.05),
    se_ZX = rep(0.05, nSnps),
    stringsAsFactors = FALSE
  ))

  results <- suppressMessages(
    runSensitivityAnalyses(
      perSnp,
      methods = c("IVW", "MREgger", "WeightedMedian"),
      engine = "internal"
    )
  )

  expect_true("summary" %in% names(results))
  expect_true("IVW" %in% results$summary$method)
  expect_true("MR-Egger" %in% results$summary$method)
  expect_true("Weighted Median" %in% results$summary$method)
})

test_that("IVW p-value is computed correctly", {
  perSnp <- with_harmonisation_columns(data.frame(
    snp_id = c("rs1", "rs2", "rs3"),
    beta_ZY = c(0, 0, 0),
    se_ZY = c(0.1, 0.1, 0.1),
    beta_ZX = c(0.3, 0.3, 0.3),
    se_ZX = c(0.05, 0.05, 0.05),
    stringsAsFactors = FALSE
  ))

  results <- suppressMessages(
    runSensitivityAnalyses(perSnp, methods = "IVW", engine = "internal")
  )

  expect_equal(results$ivw$beta_MR, 0, tolerance = 1e-10)
  expect_true(results$ivw$pval > 0.5)
})

test_that("all SNPs failing Steiger produces warning", {
  perSnp <- with_harmonisation_columns(data.frame(
    snp_id = c("rs1", "rs2", "rs3"),
    beta_ZY = c(0.5, 0.4, 0.6),
    se_ZY = rep(0.02, 3),
    beta_ZX = c(0.01, 0.01, 0.01),
    se_ZX = rep(0.05, 3),
    stringsAsFactors = FALSE
  ))

  expect_warning(
    results <- suppressMessages(
      runSensitivityAnalyses(
        perSnp,
        methods = "Steiger",
        outcomeSampleSize = 10000,
        exposureSampleSize = 10000,
        outcomeType = "continuous",
        engine = "internal"
      )
    ),
    "All SNPs failed Steiger filter"
  )
  expect_true(is.na(results$steiger$beta_MR))
})

test_that("MR-Egger is invariant to SNP sign coding after orientation", {
  perSnp <- with_harmonisation_columns(data.frame(
    snp_id = c("rs1", "rs2", "rs3", "rs4"),
    beta_ZY = c(0.14, 0.16, 0.18, 0.20),
    se_ZY = rep(0.02, 4),
    beta_ZX = c(0.28, 0.32, 0.36, 0.40),
    se_ZX = rep(0.03, 4),
    stringsAsFactors = FALSE
  ))

  flipped <- perSnp
  flipped$beta_ZX[c(2, 4)] <- -flipped$beta_ZX[c(2, 4)]
  flipped$beta_ZY[c(2, 4)] <- -flipped$beta_ZY[c(2, 4)]

  baseResult <- suppressMessages(
    runSensitivityAnalyses(perSnp, methods = "MREgger", engine = "internal")
  )
  flippedResult <- suppressMessages(
    runSensitivityAnalyses(flipped, methods = "MREgger", engine = "internal")
  )

  expect_equal(baseResult$mrEgger$beta_MR, flippedResult$mrEgger$beta_MR, tolerance = 1e-10)
  expect_equal(baseResult$mrEgger$intercept, flippedResult$mrEgger$intercept, tolerance = 1e-10)
})

test_that("Steiger requires sample sizes for a formal directional test", {
  perSnp <- with_harmonisation_columns(data.frame(
    snp_id = c("rs1", "rs2", "rs3"),
    beta_ZY = c(0.05, 0.06, 0.07),
    se_ZY = rep(0.02, 3),
    beta_ZX = c(0.20, 0.22, 0.24),
    se_ZX = rep(0.03, 3),
    stringsAsFactors = FALSE
  ))

  expect_warning(
    results <- suppressMessages(
      runSensitivityAnalyses(perSnp, methods = "Steiger", engine = "internal")
    ),
    "not implemented for binary outcomes"
  )

  expect_true(is.na(results$steiger$beta_MR))
})

test_that("Steiger requires sample sizes once continuous mode is requested", {
  perSnp <- with_harmonisation_columns(data.frame(
    snp_id = c("rs1", "rs2", "rs3"),
    beta_ZY = c(0.05, 0.06, 0.07),
    se_ZY = rep(0.02, 3),
    beta_ZX = c(0.20, 0.22, 0.24),
    se_ZX = rep(0.03, 3),
    stringsAsFactors = FALSE
  ))

  expect_warning(
    results <- suppressMessages(
      runSensitivityAnalyses(
        perSnp,
        methods = "Steiger",
        outcomeType = "continuous",
        engine = "internal"
      )
    ),
    "requires both outcomeSampleSize and exposureSampleSize"
  )

  expect_true(is.na(results$steiger$beta_MR))
})

test_that("TwoSampleMR engine returns harmonised summary estimates", {
  skip_if_not_installed("TwoSampleMR")

  perSnp <- with_harmonisation_columns(data.frame(
    snp_id = paste0("rs", 1:4),
    beta_ZY = c(0.05, 0.06, 0.05, 0.055),
    se_ZY = rep(0.02, 4),
    beta_ZX = c(0.10, 0.12, 0.09, 0.11),
    se_ZX = rep(0.02, 4),
    stringsAsFactors = FALSE
  ))

  results <- suppressMessages(
    runSensitivityAnalyses(
      perSnp,
      methods = c("IVW", "MREgger", "WeightedMedian"),
      outcomeSampleSize = 1000,
      exposureSampleSize = 1000,
      outcomeType = "continuous",
      engine = "TwoSampleMR"
    )
  )

  expect_true(is.finite(results$ivw$beta_MR))
  expect_true(is.finite(results$mrEgger$intercept))
  expect_true(is.finite(results$weightedMedian$beta_MR))
})

test_that("TwoSampleMR helper builders populate summary-data metadata", {
  perSnp <- with_harmonisation_columns(data.frame(
    snp_id = c("rs1", "rs2"),
    beta_ZY = c(0.05, 0.06),
    se_ZY = c(0.02, 0.02),
    beta_ZX = c(0.10, 0.12),
    se_ZX = c(0.02, 0.02),
    units_exposure = c("SD", "SD"),
    units_outcome = c("log odds", "log odds"),
    ncase_outcome = c(100, 100),
    ncontrol_outcome = c(900, 900),
    stringsAsFactors = FALSE
  ))

  exposureDat <- buildTwoSampleMRData(
    perSnpEstimates = perSnp,
    type = "exposure",
    sampleSize = 1000,
    outcomeType = "binary"
  )
  outcomeDat <- buildTwoSampleMRData(
    perSnpEstimates = perSnp,
    type = "outcome",
    sampleSize = 1000,
    outcomeType = "binary"
  )
  harmonised <- data.frame(SNP = perSnp$snp_id, stringsAsFactors = FALSE)
  harmonised <- addTwoSampleMRUnits(harmonised, perSnp, outcomeType = "binary")

  expect_true(all(c("samplesize", "ncase", "ncontrol") %in% names(outcomeDat)))
  expect_false("units" %in% names(exposureDat))
  expect_false("units" %in% names(outcomeDat))
  expect_equal(unique(harmonised$units.exposure), "SD")
  expect_equal(unique(harmonised$units.outcome), "log odds")
  expect_true(identical(resolveSensitivityEngine("auto"), "TwoSampleMR"))
})

test_that("TwoSampleMR helper builders preserve explicit labels, ids, and default binary units", {
  perSnp <- with_harmonisation_columns(data.frame(
    snp_id = c("rs1", "rs2"),
    beta_ZY = c(0.05, 0.06),
    se_ZY = c(0.02, 0.02),
    beta_ZX = c(0.10, 0.12),
    se_ZX = c(0.02, 0.02),
    exposure_label = c("BMI", "BMI"),
    id_exposure = c("exp-1", "exp-1"),
    outcome_label = c("CAD", "CAD"),
    id_outcome = c("out-1", "out-1"),
    stringsAsFactors = FALSE
  ))

  exposureDat <- buildTwoSampleMRData(perSnp, type = "exposure", outcomeType = "binary")
  outcomeDat <- buildTwoSampleMRData(perSnp, type = "outcome", outcomeType = "binary")
  harmonised <- addTwoSampleMRUnits(
    harmonised = data.frame(SNP = perSnp$snp_id, stringsAsFactors = FALSE),
    perSnpEstimates = perSnp,
    outcomeType = "binary"
  )

  expect_equal(unique(exposureDat$phenotype), "BMI")
  expect_equal(unique(exposureDat$id), "exp-1")
  expect_equal(unique(outcomeDat$phenotype), "CAD")
  expect_equal(unique(outcomeDat$id), "out-1")
  expect_equal(unique(harmonised$units.outcome), "log odds")
})

test_that("runSensitivityAnalyses TwoSampleMR helpers cover Steiger and leave-one-out branches", {
  skip_if_not_installed("mockery")

  harmonised <- data.frame(
    SNP = c("rs1", "rs2"),
    b = c(0.1, 0.2),
    se = c(0.05, 0.06),
    p = c(0.01, 0.02),
    stringsAsFactors = FALSE
  )

  steigerFn <- computeSteigerTwoSampleMR
  mockery::stub(steigerFn, "TwoSampleMR::steiger_filtering", function(...) {
    data.frame(
      SNP = c("rs1", "rs2"),
      steiger_dir = c(TRUE, FALSE),
      steiger_pval = c(0.01, 0.2),
      stringsAsFactors = FALSE
    )
  })
  mockery::stub(steigerFn, "computeIVWTwoSampleMR", function(...) {
    data.frame(
      method = "IVW",
      beta_MR = 0.5,
      se_MR = 0.1,
      ci_lower = 0.3,
      ci_upper = 0.7,
      pval = 0.001,
      cochran_Q = NA_real_,
      cochran_Q_pval = NA_real_,
      stringsAsFactors = FALSE
    )
  })
  steigerResult <- suppressWarnings(steigerFn(harmonised))
  expect_equal(steigerResult$n_removed, 1)
  expect_equal(steigerResult$n_remaining, 1)
  expect_equal(steigerResult$method, "Steiger-filtered IVW")

  steigerErrorFn <- computeSteigerTwoSampleMR
  mockery::stub(steigerErrorFn, "TwoSampleMR::steiger_filtering", function(...) stop("boom"))
  expect_warning(
    unavailable <- steigerErrorFn(harmonised),
    "could not be evaluated"
  )
  expect_true(is.na(unavailable$beta_MR))

  steigerMissingFn <- computeSteigerTwoSampleMR
  mockery::stub(steigerMissingFn, "TwoSampleMR::steiger_filtering", function(...) {
    data.frame(SNP = "rs1", stringsAsFactors = FALSE)
  })
  expect_warning(
    missingDir <- steigerMissingFn(harmonised),
    "did not return steiger_dir"
  )
  expect_true(is.na(missingDir$beta_MR))

  steigerAllFailFn <- computeSteigerTwoSampleMR
  mockery::stub(steigerAllFailFn, "TwoSampleMR::steiger_filtering", function(...) {
    data.frame(
      SNP = c("rs1", "rs2"),
      steiger_dir = c(FALSE, FALSE),
      steiger_pval = c(0.4, 0.5),
      stringsAsFactors = FALSE
    )
  })
  expect_warning(
    allFail <- steigerAllFailFn(harmonised),
    "All SNPs failed Steiger filter"
  )
  expect_equal(allFail$n_remaining, 0)

  looFn <- computeLeaveOneOutTwoSampleMR
  mockery::stub(looFn, "TwoSampleMR::mr_leaveoneout", function(...) {
    data.frame(
      SNP = c("rs1", "rs2", "All"),
      b = c(0.1, 0.2, 0.15),
      se = c(0.05, 0.06, 0.04),
      p = c(0.05, 0.02, 0.01),
      stringsAsFactors = FALSE
    )
  })
  loo <- looFn(harmonised)
  expect_equal(nrow(loo), 2)
  expect_false("All" %in% loo$snp_removed)
})

test_that("runSensitivityAnalyses delegates Steiger and leave-one-out through TwoSampleMR hooks", {
  skip_if_not_installed("mockery")

  perSnp <- with_harmonisation_columns(data.frame(
    snp_id = c("rs1", "rs2", "rs3"),
    beta_ZY = c(0.04, 0.05, 0.06),
    se_ZY = c(0.02, 0.02, 0.02),
    beta_ZX = c(0.10, 0.11, 0.12),
    se_ZX = c(0.02, 0.02, 0.02),
    stringsAsFactors = FALSE
  ))

  runFn <- runSensitivityAnalyses
  mockery::stub(runFn, "resolveSensitivityEngine", function(...) "TwoSampleMR")
  mockery::stub(
    runFn,
    "harmonisePerSnpEstimates",
    function(...) data.frame(SNP = perSnp$snp_id, stringsAsFactors = FALSE)
  )
  mockery::stub(
    runFn,
    "computeSteigerTwoSampleMR",
    function(...) {
      list(
        method = "Steiger-filtered IVW",
        beta_MR = 0.4,
        se_MR = 0.1,
        ci_lower = 0.2,
        ci_upper = 0.6,
        pval = 0.01,
        n_removed = 1L,
        n_remaining = 2L
      )
    }
  )
  mockery::stub(
    runFn,
    "computeLeaveOneOutTwoSampleMR",
    function(...) {
      data.frame(
        snp_removed = c("rs1", "rs2", "rs3"),
        beta_MR = c(0.41, 0.40, 0.39),
        se_MR = c(0.1, 0.1, 0.1),
        pval = c(0.02, 0.03, 0.04),
        stringsAsFactors = FALSE
      )
    }
  )

  results <- suppressMessages(
    runFn(
      perSnpEstimates = perSnp,
      methods = c("Steiger", "LeaveOneOut"),
      outcomeSampleSize = 1000,
      exposureSampleSize = 1000,
      outcomeType = "continuous",
      engine = "TwoSampleMR"
    )
  )

  expect_equal(results$steiger$method, "Steiger-filtered IVW")
  expect_equal(nrow(results$leaveOneOut), 3)
})

test_that("runSensitivityAnalyses helper branches cover fallback and error paths", {
  skip_if_not_installed("mockery")

  resolveFn <- resolveSensitivityEngine
  mockery::stub(resolveFn, "requireNamespace", function(package, ...) FALSE)
  expect_warning(expect_equal(resolveFn("auto"), "internal"), "falling back")
  expect_error(resolveFn("TwoSampleMR"), "Package 'TwoSampleMR' is required")

  emptyHarmoniseFn <- harmonisePerSnpEstimates
  basePerSnp <- with_harmonisation_columns(data.frame(
    snp_id = c("rs1", "rs2"),
    beta_ZY = c(0.1, 0.2),
    se_ZY = c(0.05, 0.05),
    beta_ZX = c(0.2, 0.3),
    se_ZX = c(0.05, 0.05),
    stringsAsFactors = FALSE
  ))
  mockery::stub(emptyHarmoniseFn, "TwoSampleMR::format_data", function(dat, ...) dat)
  mockery::stub(emptyHarmoniseFn, "TwoSampleMR::harmonise_data", function(...) {
    data.frame(SNP = c("rs1", "rs2"), mr_keep = c(FALSE, FALSE), stringsAsFactors = FALSE)
  })
  expect_error(
    emptyHarmoniseFn(basePerSnp, outcomeType = "continuous"),
    "No SNPs remained after harmonisation"
  )

  steigerSingle <- with_harmonisation_columns(data.frame(
    snp_id = c("rs1", "rs2"),
    beta_ZY = c(0.1, 0.01),
    se_ZY = c(0.01, 0.01),
    beta_ZX = c(0.2, 0.001),
    se_ZX = c(0.01, 0.01),
    stringsAsFactors = FALSE
  ))
  singleSteiger <- suppressMessages(
    computeSteiger(
      steigerSingle,
      outcomeSampleSize = 10000,
      exposureSampleSize = 10000,
      outcomeType = "continuous"
    )
  )
  expect_true(is.na(singleSteiger$se_MR))

  expect_error(
    computeRatioSummary(
      with_harmonisation_columns(data.frame(
        snp_id = "rs1",
        beta_ZY = 0.1,
        se_ZY = 0.05,
        beta_ZX = 0,
        se_ZX = 0.05,
        stringsAsFactors = FALSE
      ))
    ),
    "Cannot compute ratio estimates"
  )
})
