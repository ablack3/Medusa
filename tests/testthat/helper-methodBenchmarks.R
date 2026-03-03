# Monte Carlo benchmark helpers for methodological calibration tests.

runMonteCarloMRBenchmark <- function(nRep = 8,
                                     n = 1000,
                                     nSnps = 6,
                                     trueEffect = 0.3,
                                     betaGrid = seq(-1.5, 1.5, by = 0.15),
                                     baseSeed = 5000,
                                     alpha = 0.05) {
  checkmate::assertCount(nRep, positive = TRUE)
  checkmate::assertCount(n, positive = TRUE)
  checkmate::assertCount(nSnps, positive = TRUE)
  checkmate::assertNumber(trueEffect)
  validateBetaGrid(betaGrid)
  checkmate::assertNumber(alpha, lower = 0.001, upper = 0.25)

  replicateRows <- lapply(seq_len(nRep), function(idx) {
    simData <- simulateMRData(
      n = n,
      nSnps = nSnps,
      trueEffect = trueEffect,
      seed = baseSeed + idx
    )
    covariateData <- simData$data[, c("person_id", "confounder_1", "confounder_2"), drop = FALSE]

    siteProfile <- suppressWarnings(
      suppressMessages(
        fitOutcomeModel(
          cohortData = simData$data,
          covariateData = covariateData,
          instrumentTable = simData$instrumentTable,
          betaGrid = betaGrid,
          analysisType = "perSNP",
          siteId = sprintf("benchmark_%02d", idx)
        )
      )
    )
    combinedProfile <- suppressWarnings(
      suppressMessages(poolLikelihoodProfiles(list(siteProfile)))
    )
    medusaEstimate <- suppressWarnings(
      suppressMessages(
        computeMREstimate(
          combinedProfile = combinedProfile,
          instrumentTable = simData$instrumentTable
        )
      )
    )
    sensitivityResults <- suppressWarnings(
      suppressMessages(
        runSensitivityAnalyses(
          perSnpEstimates = siteProfile$perSnpEstimates,
          methods = c("IVW", "WeightedMedian"),
          engine = "internal"
        )
      )
    )

    data.frame(
      replicate = idx,
      method = c("Medusa", "IVW", "Weighted Median"),
      beta_MR = c(
        medusaEstimate$betaMR,
        sensitivityResults$ivw$beta_MR,
        sensitivityResults$weightedMedian$beta_MR
      ),
      se_MR = c(
        medusaEstimate$seMR,
        sensitivityResults$ivw$se_MR,
        sensitivityResults$weightedMedian$se_MR
      ),
      ci_lower = c(
        medusaEstimate$ciLower,
        sensitivityResults$ivw$ci_lower,
        sensitivityResults$weightedMedian$ci_lower
      ),
      ci_upper = c(
        medusaEstimate$ciUpper,
        sensitivityResults$ivw$ci_upper,
        sensitivityResults$weightedMedian$ci_upper
      ),
      pval = c(
        medusaEstimate$pValue,
        sensitivityResults$ivw$pval,
        sensitivityResults$weightedMedian$pval
      ),
      stringsAsFactors = FALSE
    )
  })

  replicateResults <- do.call(rbind, replicateRows)
  rownames(replicateResults) <- NULL
  replicateResults$covered <- with(
    replicateResults,
    ci_lower <= trueEffect & ci_upper >= trueEffect
  )
  replicateResults$reject_null <- replicateResults$pval < alpha
  replicateResults$abs_error <- abs(replicateResults$beta_MR - trueEffect)
  replicateResults$sq_error <- (replicateResults$beta_MR - trueEffect)^2

  splitByMethod <- split(replicateResults, replicateResults$method)
  summaryRows <- lapply(names(splitByMethod), function(methodName) {
    methodRows <- splitByMethod[[methodName]]
    data.frame(
      method = methodName,
      n_rep = nrow(methodRows),
      mean_estimate = mean(methodRows$beta_MR),
      bias = mean(methodRows$beta_MR - trueEffect),
      empirical_sd = stats::sd(methodRows$beta_MR),
      mean_reported_se = mean(methodRows$se_MR),
      rmse = sqrt(mean(methodRows$sq_error)),
      coverage = mean(methodRows$covered),
      type1_error = mean(methodRows$reject_null),
      mean_pval = mean(methodRows$pval),
      stringsAsFactors = FALSE
    )
  })

  summaryTable <- do.call(rbind, summaryRows)
  rownames(summaryTable) <- NULL

  list(
    truth = trueEffect,
    alpha = alpha,
    replicates = replicateResults,
    summary = summaryTable
  )
}


getBenchmarkMetric <- function(benchmark, method, metric) {
  checkmate::assertList(benchmark)
  checkmate::assertDataFrame(benchmark$summary, min.rows = 1)
  checkmate::assertChoice(method, benchmark$summary$method)
  checkmate::assertChoice(metric, names(benchmark$summary))

  benchmark$summary[benchmark$summary$method == method, metric, drop = TRUE][[1]]
}
