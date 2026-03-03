test_that("Monte Carlo MR benchmark returns calibration summaries", {
  # This benchmark is intended for explicit pre-release validation only.
  skip_if_not_full_validation()

  benchmark <- runMonteCarloMRBenchmark(
    nRep = 4,
    n = 700,
    nSnps = 5,
    trueEffect = 0.25,
    betaGrid = seq(-1.25, 1.25, by = 0.2),
    baseSeed = 6100
  )

  expect_named(benchmark, c("truth", "alpha", "replicates", "summary"))
  expect_s3_class(benchmark$replicates, "data.frame")
  expect_s3_class(benchmark$summary, "data.frame")
  expect_equal(nrow(benchmark$replicates), 12)
  expect_equal(sort(benchmark$summary$method), sort(c("IVW", "Medusa", "Weighted Median")))
  expect_true(all(c(
    "bias", "rmse", "coverage", "rejection_rate", "type1_error",
    "power", "mean_reported_se", "empirical_sd"
  ) %in% names(benchmark$summary)))
  expect_true(all(is.finite(benchmark$summary$rmse)))
})

test_that("Monte Carlo benchmark shows lower RMSE at larger sample sizes", {
  # This benchmark is intended for explicit pre-release validation only.
  skip_if_not_full_validation()

  smallBenchmark <- runMonteCarloMRBenchmark(
    nRep = 8,
    n = 400,
    nSnps = 6,
    trueEffect = 0.3,
    betaGrid = seq(-1.5, 1.5, by = 0.15),
    baseSeed = 7000
  )
  largeBenchmark <- runMonteCarloMRBenchmark(
    nRep = 8,
    n = 1200,
    nSnps = 6,
    trueEffect = 0.3,
    betaGrid = seq(-1.5, 1.5, by = 0.15),
    baseSeed = 7000
  )

  expect_lt(
    getBenchmarkMetric(largeBenchmark, "IVW", "rmse"),
    getBenchmarkMetric(smallBenchmark, "IVW", "rmse")
  )
  expect_lt(
    getBenchmarkMetric(largeBenchmark, "Weighted Median", "rmse"),
    getBenchmarkMetric(smallBenchmark, "Weighted Median", "rmse")
  )
})

test_that("Monte Carlo benchmark gives broadly acceptable coverage and null calibration", {
  # This benchmark is intended for explicit pre-release validation only.
  skip_if_not_full_validation()

  signalBenchmark <- runMonteCarloMRBenchmark(
    nRep = 8,
    n = 900,
    nSnps = 6,
    trueEffect = 0.15,
    betaGrid = seq(-1.5, 1.5, by = 0.15),
    baseSeed = 7100
  )
  nullBenchmark <- runMonteCarloMRBenchmark(
    nRep = 8,
    n = 900,
    nSnps = 6,
    trueEffect = 0,
    betaGrid = seq(-1.5, 1.5, by = 0.15),
    baseSeed = 7200
  )

  expect_gte(getBenchmarkMetric(signalBenchmark, "IVW", "coverage"), 0.75)
  expect_gte(getBenchmarkMetric(signalBenchmark, "Weighted Median", "coverage"), 0.75)
  expect_lte(abs(getBenchmarkMetric(signalBenchmark, "IVW", "bias")), 0.10)
  expect_lte(abs(getBenchmarkMetric(signalBenchmark, "Weighted Median", "bias")), 0.10)
  expect_true(is.na(getBenchmarkMetric(signalBenchmark, "IVW", "type1_error")))
  expect_gte(getBenchmarkMetric(signalBenchmark, "IVW", "power"), 0.05)
  expect_gte(getBenchmarkMetric(signalBenchmark, "Weighted Median", "power"), 0.05)

  expect_lte(getBenchmarkMetric(nullBenchmark, "Medusa", "type1_error"), 0.20)
  expect_lte(getBenchmarkMetric(nullBenchmark, "IVW", "type1_error"), 0.20)
  expect_lte(getBenchmarkMetric(nullBenchmark, "Weighted Median", "type1_error"), 0.20)
  expect_true(is.na(getBenchmarkMetric(nullBenchmark, "Medusa", "power")))
  expect_lte(abs(getBenchmarkMetric(nullBenchmark, "IVW", "mean_estimate")), 0.10)
  expect_lte(abs(getBenchmarkMetric(nullBenchmark, "Weighted Median", "mean_estimate")), 0.10)
  expect_gt(
    getBenchmarkMetric(signalBenchmark, "Medusa", "mean_estimate"),
    getBenchmarkMetric(nullBenchmark, "Medusa", "mean_estimate")
  )
})
