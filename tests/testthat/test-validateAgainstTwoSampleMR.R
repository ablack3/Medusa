test_that("validateAgainstTwoSampleMR returns a comparison table", {
  # Single-site cross-package validation is useful before release, but not needed by default.
  skip_if_not_full_validation()
  skip_if_not_installed("TwoSampleMR")

  simData <- simulateMRData(n = 2500, nSnps = 4, trueEffect = 0.25, seed = 777)
  siteProfile <- suppressMessages(
    fitOutcomeModel(
      cohortData = simData$data,
      covariateData = NULL,
      instrumentTable = simData$instrumentTable,
      betaGrid = seq(-2, 2, by = 0.05),
      analysisType = "perSNP"
    )
  )

  comparison <- suppressMessages(
    validateAgainstTwoSampleMR(
      siteProfile = siteProfile,
      instrumentTable = simData$instrumentTable,
      methods = c("IVW", "WeightedMedian"),
      outcomeSampleSize = nrow(simData$data),
      exposureSampleSize = nrow(simData$data),
      outcomeType = "continuous"
    )
  )

  expect_s3_class(comparison, "data.frame")
  expect_true("Profile likelihood" %in% comparison$method)
  expect_true("IVW" %in% comparison$method)
  expect_true("delta_vs_medusa" %in% names(comparison))
})

test_that("validateAgainstTwoSampleMR requires per-SNP estimates on the site profile", {
  simData <- simulateMRData(n = 200, nSnps = 2, seed = 778)
  siteProfile <- list(
    siteId = "site_1",
    betaGrid = seq(-1, 1, by = 0.1),
    logLikProfile = rep(0, 21),
    nCases = 20,
    nControls = 180,
    snpIds = simData$instrumentTable$snp_id
  )

  expect_error(
    validateAgainstTwoSampleMR(
      siteProfile = siteProfile,
      instrumentTable = simData$instrumentTable
    ),
    "must contain perSnpEstimates"
  )
})

test_that("ensureHarmonisationColumns backfills allele columns and errors on unmatched SNPs", {
  instrumentTable <- simulateInstrumentTable(nSnps = 2, seed = 779)
  perSnp <- data.frame(
    snp_id = instrumentTable$snp_id,
    beta_ZY = c(0.1, 0.2),
    se_ZY = c(0.05, 0.05),
    beta_ZX = instrumentTable$beta_ZX,
    se_ZX = instrumentTable$se_ZX,
    stringsAsFactors = FALSE
  )

  enriched <- ensureHarmonisationColumns(perSnp, instrumentTable)
  badPerSnp <- perSnp
  badPerSnp$snp_id[1] <- "rs_missing"

  expect_true(all(c("effect_allele", "other_allele", "eaf", "pval_ZX") %in% names(enriched)))
  expect_error(
    ensureHarmonisationColumns(badPerSnp, instrumentTable),
    "could not be matched"
  )
})
