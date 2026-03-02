test_that("exportSiteProfile writes two CSV files", {
  profiles <- simulateSiteProfiles(nSites = 1, trueBeta = 0.5, seed = 99)
  profile <- profiles[[1]]

  tmpDir <- tempfile("medusa_export_")
  dir.create(tmpDir)
  on.exit(unlink(tmpDir, recursive = TRUE))

  paths <- exportSiteProfile(profile, outputDir = tmpDir)

  expect_true(file.exists(paths[["profile"]]))
  expect_true(file.exists(paths[["metadata"]]))
  expect_true(grepl("\\.csv$", paths[["profile"]]))
  expect_true(grepl("\\.csv$", paths[["metadata"]]))
})

test_that("exported profile CSV has expected columns", {
  profiles <- simulateSiteProfiles(nSites = 1, trueBeta = 0.3, seed = 101)
  profile <- profiles[[1]]

  tmpDir <- tempfile("medusa_export_")
  dir.create(tmpDir)
  on.exit(unlink(tmpDir, recursive = TRUE))

  paths <- exportSiteProfile(profile, outputDir = tmpDir)

  profileDf <- read.csv(paths[["profile"]])
  expect_true("beta" %in% names(profileDf))
  expect_true("log_likelihood" %in% names(profileDf))
  expect_equal(nrow(profileDf), length(profile$betaGrid))
})

test_that("exported metadata CSV has expected fields", {
  profiles <- simulateSiteProfiles(nSites = 1, trueBeta = 0.3, seed = 102)
  profile <- profiles[[1]]

  tmpDir <- tempfile("medusa_export_")
  dir.create(tmpDir)
  on.exit(unlink(tmpDir, recursive = TRUE))

  paths <- exportSiteProfile(profile, outputDir = tmpDir)

  metaDf <- read.csv(paths[["metadata"]], stringsAsFactors = FALSE)
  expect_true("site_id" %in% names(metaDf))
  expect_true("n_cases" %in% names(metaDf))
  expect_true("n_controls" %in% names(metaDf))
  expect_true("snp_ids" %in% names(metaDf))
  expect_equal(metaDf$site_id[1], profile$siteId)
  expect_equal(metaDf$n_cases[1], profile$nCases)
})

test_that("importSiteProfile round-trips correctly", {
  profiles <- simulateSiteProfiles(nSites = 1, trueBeta = 0.4, seed = 103)
  original <- profiles[[1]]

  tmpDir <- tempfile("medusa_export_")
  dir.create(tmpDir)
  on.exit(unlink(tmpDir, recursive = TRUE))

  paths <- exportSiteProfile(original, outputDir = tmpDir)
  imported <- importSiteProfile(paths[["profile"]])

  expect_equal(imported$siteId, original$siteId)
  expect_equal(imported$betaGrid, original$betaGrid)
  expect_equal(imported$logLikProfile, original$logLikProfile, tolerance = 1e-10)
  expect_equal(imported$nCases, original$nCases)
  expect_equal(imported$nControls, original$nControls)
  expect_equal(imported$snpIds, original$snpIds)
})

test_that("importSiteProfile works with poolLikelihoodProfiles", {
  profiles <- simulateSiteProfiles(nSites = 3, trueBeta = 0.5, seed = 104)

  tmpDir <- tempfile("medusa_export_")
  dir.create(tmpDir)
  on.exit(unlink(tmpDir, recursive = TRUE))

  # Export all profiles
  for (nm in names(profiles)) {
    exportSiteProfile(profiles[[nm]], outputDir = tmpDir)
  }

  # Import all profiles
  importedProfiles <- list()
  for (nm in names(profiles)) {
    siteId <- profiles[[nm]]$siteId
    path <- file.path(tmpDir, sprintf("medusa_profile_%s.csv", siteId))
    importedProfiles[[nm]] <- importSiteProfile(path)
  }

  # Pool the imported profiles
  combined <- poolLikelihoodProfiles(importedProfiles)
  expect_equal(combined$nSites, 3)
  expect_true(length(combined$logLikProfile) > 0)

  # Should match pooling from originals
  combinedOriginal <- poolLikelihoodProfiles(profiles)
  expect_equal(combined$logLikProfile, combinedOriginal$logLikProfile,
               tolerance = 1e-10)
})

test_that("importSiteProfile warns when metadata file is missing", {
  profiles <- simulateSiteProfiles(nSites = 1, trueBeta = 0.3, seed = 105)
  profile <- profiles[[1]]

  tmpDir <- tempfile("medusa_export_")
  dir.create(tmpDir)
  on.exit(unlink(tmpDir, recursive = TRUE))

  paths <- exportSiteProfile(profile, outputDir = tmpDir)

  # Remove metadata file
  file.remove(paths[["metadata"]])

  expect_warning(
    imported <- importSiteProfile(paths[["profile"]]),
    "Metadata file not found"
  )
  expect_equal(imported$siteId, "unknown")
  expect_equal(imported$betaGrid, profile$betaGrid)
})

test_that("importSiteProfile errors on missing profile file", {
  expect_error(
    importSiteProfile("nonexistent_file.csv"),
    "Profile file not found"
  )
})

test_that("exportSiteProfile creates output directory if needed", {
  profiles <- simulateSiteProfiles(nSites = 1, trueBeta = 0.3, seed = 106)
  profile <- profiles[[1]]

  tmpDir <- file.path(tempdir(), "medusa_nested", "subdir")
  on.exit(unlink(file.path(tempdir(), "medusa_nested"), recursive = TRUE))

  paths <- exportSiteProfile(profile, outputDir = tmpDir)
  expect_true(file.exists(paths[["profile"]]))
})

test_that("exportSiteProfile respects custom prefix", {
  profiles <- simulateSiteProfiles(nSites = 1, trueBeta = 0.3, seed = 107)
  profile <- profiles[[1]]

  tmpDir <- tempfile("medusa_export_")
  dir.create(tmpDir)
  on.exit(unlink(tmpDir, recursive = TRUE))

  paths <- exportSiteProfile(profile, outputDir = tmpDir, prefix = "mystudy")
  expect_true(grepl("^mystudy_profile_", basename(paths[["profile"]])))
  expect_true(grepl("^mystudy_metadata_", basename(paths[["metadata"]])))
})
