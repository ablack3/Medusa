test_that("generateMRReport renders a report with a stubbed renderer", {
  skip_if_not_installed("mockery")

  reportFn <- generateMRReport
  outputPath <- file.path(tempdir(), "medusa-test-report.html")
  if (file.exists(outputPath)) {
    unlink(outputPath)
  }

  mockery::stub(
    reportFn,
    "rmarkdown::render",
    function(input, output_file, output_dir, ...) {
      writeLines("<html><body>ok</body></html>", file.path(output_dir, output_file))
      invisible(file.path(output_dir, output_file))
    }
  )

  result <- suppressMessages(
    reportFn(
      mrEstimate = list(betaMR = 0.2),
      sensitivityResults = NULL,
      diagnosticResults = NULL,
      combinedProfile = list(betaGrid = c(-1, 0, 1), logLikProfile = c(-1, 0, -1)),
      outputPath = outputPath
    )
  )

  expect_equal(result, outputPath)
  expect_true(file.exists(outputPath))
})

test_that("generateMRReport wraps renderer failures with a report-specific error", {
  skip_if_not_installed("mockery")

  reportFn <- generateMRReport
  mockery::stub(
    reportFn,
    "rmarkdown::render",
    function(...) stop("render failed")
  )

  expect_error(
    suppressMessages(
      reportFn(
        mrEstimate = list(betaMR = 0.2),
        sensitivityResults = NULL,
        diagnosticResults = NULL,
        combinedProfile = list(betaGrid = c(-1, 0, 1), logLikProfile = c(-1, 0, -1)),
        outputPath = file.path(tempdir(), "medusa-test-report-fail.html")
      )
    ),
    "Report generation failed"
  )
})

test_that("plot helpers return ggplot objects for populated and empty inputs", {
  profiles <- simulateSiteProfiles(nSites = 2, trueBeta = 0.3, seed = 909)
  combined <- poolLikelihoodProfiles(profiles)
  sensitivityResults <- list(
    summary = data.frame(
      method = c("IVW", "Weighted Median"),
      beta_MR = c(0.2, 0.15),
      ci_lower = c(0.1, 0.05),
      ci_upper = c(0.3, 0.25),
      stringsAsFactors = FALSE
    )
  )
  emptySensitivity <- list(summary = data.frame())

  likelihoodPlot <- plotLikelihoodProfile(
    combinedProfile = combined,
    siteProfileList = profiles,
    mrEstimate = list(betaZY = 0.3, ciLevel = 0.95)
  )
  forestPlot <- plotSensitivityForest(sensitivityResults)
  emptyPlot <- plotSensitivityForest(emptySensitivity)

  expect_s3_class(likelihoodPlot, "ggplot")
  expect_s3_class(forestPlot, "ggplot")
  expect_s3_class(emptyPlot, "ggplot")
})
