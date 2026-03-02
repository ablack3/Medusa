test_that("pointwise sum equals manual sum of vectors", {
  profiles <- simulateSiteProfiles(nSites = 3, trueBeta = 0.5, seed = 42)

  combined <- poolLikelihoodProfiles(profiles)

  # Manual sum
  manualSum <- profiles[[1]]$logLikProfile +
    profiles[[2]]$logLikProfile +
    profiles[[3]]$logLikProfile
  # Normalize both
  manualSum <- manualSum - max(manualSum)

  expect_equal(combined$logLikProfile, manualSum, tolerance = 1e-10)
})

test_that("pooled MLE recovers true beta within tolerance for large N", {
  profiles <- simulateSiteProfiles(
    nSites = 5,
    trueBeta = 0.5,
    nPerSite = 5000,
    seed = 123
  )

  combined <- poolLikelihoodProfiles(profiles)
  peakIdx <- which.max(combined$logLikProfile)
  pooledMLE <- combined$betaGrid[peakIdx]

  # With 5 sites of 5000 each, should be very close to 0.5
  expect_equal(pooledMLE, 0.5, tolerance = 0.15)
})

test_that("misaligned grids trigger warning and interpolation", {
  grid1 <- seq(-3, 3, by = 0.01)
  grid2 <- seq(-3, 3, by = 0.02)  # Different resolution

  profiles <- list(
    site_A = list(
      siteId = "site_A",
      betaGrid = grid1,
      logLikProfile = -0.5 * 10 * (grid1 - 0.5)^2,
      nCases = 100,
      nControls = 900,
      snpIds = "rs1"
    ),
    site_B = list(
      siteId = "site_B",
      betaGrid = grid2,
      logLikProfile = -0.5 * 10 * (grid2 - 0.5)^2,
      nCases = 80,
      nControls = 920,
      snpIds = "rs1"
    )
  )

  expect_warning(
    combined <- poolLikelihoodProfiles(profiles),
    "different betaGrid"
  )

  # Should still produce a valid result
  expect_true(length(combined$logLikProfile) > 0)
  peakIdx <- which.max(combined$logLikProfile)
  expect_equal(combined$betaGrid[peakIdx], 0.5, tolerance = 0.05)
})

test_that("single-site pooling returns same result as un-pooled profile", {
  profiles <- simulateSiteProfiles(nSites = 1, trueBeta = 0.3, seed = 42)

  combined <- poolLikelihoodProfiles(profiles)

  # Should be the same profile (after normalization)
  origProfile <- profiles[[1]]$logLikProfile - max(profiles[[1]]$logLikProfile)
  expect_equal(combined$logLikProfile, origProfile, tolerance = 1e-10)
})

test_that("site contribution table is correct", {
  profiles <- simulateSiteProfiles(nSites = 3, trueBeta = 0.5, seed = 42)

  combined <- poolLikelihoodProfiles(profiles)

  expect_equal(combined$nSites, 3)
  expect_equal(nrow(combined$siteContributions), 3)
  expect_equal(combined$totalCases,
               sum(combined$siteContributions$nCases))
  expect_equal(combined$totalControls,
               sum(combined$siteContributions$nControls))
})

test_that("poolLikelihoodProfiles validates input", {
  expect_error(
    poolLikelihoodProfiles(list()),
    "length >= 1"
  )

  expect_error(
    poolLikelihoodProfiles(list(list(siteId = "A"))),
    "missing required"
  )
})
