test_that("simulateMRData returns coherent simulated cohort and instrument data", {
  simData <- Medusa::simulateMRData(n = 250, nSnps = 4, trueEffect = 0.4, seed = 101)

  expect_named(simData, c("data", "instrumentTable", "trueEffect"))
  expect_equal(simData$trueEffect, 0.4)
  expect_equal(nrow(simData$data), 250)
  expect_equal(nrow(simData$instrumentTable), 4)
  expect_true(all(grepl("^snp_", grep("^snp_", names(simData$data), value = TRUE))))
  expect_true(all(c("fStatistic", "strandAmbiguous") %in% names(simData$instrumentTable)))
  expect_true(all(simData$instrumentTable$effect_allele != simData$instrumentTable$other_allele))
})

test_that("simulateInstrumentTable produces valid non-ambiguous instruments", {
  instruments <- Medusa::simulateInstrumentTable(nSnps = 6, seed = 202)

  expect_equal(nrow(instruments), 6)
  expect_true(all(instruments$eaf > 0.05))
  expect_true(all(instruments$eaf < 0.95))
  expect_false(any(isStrandAmbiguous(
    instruments$effect_allele,
    instruments$other_allele
  )))
})

test_that("simulateSiteProfiles and simulateCovariateData return expected structures", {
  profiles <- Medusa::simulateSiteProfiles(nSites = 3, trueBeta = 0.2, seed = 303)
  covariates <- Medusa::simulateCovariateData(n = 50, nCovariates = 6, seed = 404)

  expect_length(profiles, 3)
  expect_true(all(vapply(profiles, function(x) max(x$logLikProfile), numeric(1)) == 0))
  expect_true(all(vapply(profiles, function(x) x$nCases + x$nControls, numeric(1)) > 0))

  expect_equal(nrow(covariates), 50)
  expect_equal(ncol(covariates), 7)
  expect_true("person_id" %in% names(covariates))
})
