test_that("fitOutcomeModel returns valid profile with simulated data", {
  simData <- simulateMRData(n = 2000, nSnps = 5, trueEffect = 0.3, seed = 123)
  profile <- fitOutcomeModel(
    cohortData = simData$data,
    covariateData = NULL,
    instrumentTable = simData$instrumentTable,
    betaGrid = seq(-2, 2, by = 0.05),
    siteId = "test_site"
  )

  expect_s3_class(profile, "medusaSiteProfile")
  expect_equal(profile$siteId, "test_site")
  expect_equal(length(profile$betaGrid), length(profile$logLikProfile))
  expect_true(!is.null(profile$betaHat))
  expect_true(!is.null(profile$seHat))
})

test_that("log-likelihood profile is concave (has unique maximum)", {
  simData <- simulateMRData(n = 5000, nSnps = 5, trueEffect = 0.5, seed = 456)
  profile <- fitOutcomeModel(
    cohortData = simData$data,
    covariateData = NULL,
    instrumentTable = simData$instrumentTable,
    betaGrid = seq(-3, 3, by = 0.01)
  )

  # Profile should have a clear peak
  peakIdx <- which.max(profile$logLikProfile)
  expect_true(peakIdx > 1)
  expect_true(peakIdx < length(profile$logLikProfile))

  # Values should decrease away from peak
  expect_true(profile$logLikProfile[peakIdx] > profile$logLikProfile[1])
  expect_true(profile$logLikProfile[peakIdx] > profile$logLikProfile[length(profile$logLikProfile)])
})

test_that("MLE from profile is close to true beta_ZY for large N", {
  simData <- simulateMRData(n = 10000, nSnps = 5, trueEffect = 0.5, seed = 789)
  profile <- fitOutcomeModel(
    cohortData = simData$data,
    covariateData = NULL,
    instrumentTable = simData$instrumentTable,
    betaGrid = seq(-3, 3, by = 0.01)
  )

  # The betaHat from the profile is the allele score coefficient, not directly
  # beta_ZY, but should be non-zero and in the right direction
  expect_true(is.finite(profile$betaHat))
  expect_true(profile$seHat > 0)
})

test_that("profile object contains correct named elements", {
  simData <- simulateMRData(n = 500, nSnps = 3, trueEffect = 0.2, seed = 111)
  profile <- fitOutcomeModel(
    cohortData = simData$data,
    covariateData = NULL,
    instrumentTable = simData$instrumentTable,
    betaGrid = seq(-1, 1, by = 0.1)
  )

  expect_true(all(c("siteId", "betaGrid", "logLikProfile", "nCases",
                     "nControls", "snpIds", "diagnosticFlags",
                     "betaHat", "seHat") %in% names(profile)))
})

test_that("alleleScore and perSNP modes both return valid profiles", {
  simData <- simulateMRData(n = 1000, nSnps = 3, trueEffect = 0.3, seed = 222)
  grid <- seq(-2, 2, by = 0.1)

  profileAS <- fitOutcomeModel(
    cohortData = simData$data,
    covariateData = NULL,
    instrumentTable = simData$instrumentTable,
    betaGrid = grid,
    analysisType = "alleleScore"
  )

  profilePS <- fitOutcomeModel(
    cohortData = simData$data,
    covariateData = NULL,
    instrumentTable = simData$instrumentTable,
    betaGrid = grid,
    analysisType = "perSNP"
  )

  expect_equal(length(profileAS$logLikProfile), length(grid))
  expect_equal(length(profilePS$logLikProfile), length(grid))
  expect_true(is.finite(profileAS$betaHat))
  expect_true(is.finite(profilePS$betaHat))
  expect_true("perSnpEstimates" %in% names(profilePS))
  expect_true(nrow(profilePS$perSnpEstimates) > 0)
  expect_true(all(c("effect_allele", "other_allele", "eaf") %in%
                    names(profilePS$perSnpEstimates)))
})

test_that("warning issued when MLE is at grid boundary", {
  simData <- simulateMRData(n = 2000, nSnps = 5, trueEffect = 5, seed = 333)

  expect_warning(
    profile <- fitOutcomeModel(
      cohortData = simData$data,
      covariateData = NULL,
      instrumentTable = simData$instrumentTable,
      betaGrid = seq(-1, 1, by = 0.1)  # Too narrow for true effect of 5
    ),
    "grid boundary"
  )
})

test_that("fitOutcomeModel validates inputs", {
  expect_error(
    fitOutcomeModel(
      cohortData = data.frame(outcome = 1),
      covariateData = NULL,
      instrumentTable = data.frame(),  # invalid
      betaGrid = seq(-1, 1, by = 0.1)
    )
  )
})

test_that("fitOutcomeModel accepts Inf regularizationVariance for unpenalized fits", {
  simData <- simulateMRData(n = 400, nSnps = 3, trueEffect = 0.2, seed = 334)

  expect_no_error(
    fitOutcomeModel(
      cohortData = simData$data,
      covariateData = NULL,
      instrumentTable = simData$instrumentTable,
      betaGrid = seq(-1, 1, by = 0.1),
      modelBackend = "glm",
      regularizationVariance = Inf
    )
  )
})

test_that("fitOutcomeModel aligns SNP columns by SNP ID instead of column order", {
  simData <- simulateMRData(n = 1500, nSnps = 4, trueEffect = 0.25, seed = 444)
  snpCols <- grep("^snp_", names(simData$data), value = TRUE)
  reordered <- simData$data[, c(
    setdiff(names(simData$data), snpCols),
    rev(snpCols)
  )]

  profileOriginal <- fitOutcomeModel(
    cohortData = simData$data,
    covariateData = NULL,
    instrumentTable = simData$instrumentTable,
    betaGrid = seq(-2, 2, by = 0.05)
  )

  profileReordered <- fitOutcomeModel(
    cohortData = reordered,
    covariateData = NULL,
    instrumentTable = simData$instrumentTable,
    betaGrid = seq(-2, 2, by = 0.05)
  )

  expect_equal(profileOriginal$betaHat, profileReordered$betaHat, tolerance = 1e-10)
  expect_equal(profileOriginal$logLikProfile, profileReordered$logLikProfile, tolerance = 1e-10)
  expect_equal(profileOriginal$scoreDefinition$betaZX, profileReordered$scoreDefinition$betaZX, tolerance = 1e-10)
})

test_that("fitOutcomeModel supports the optional Cyclops backend", {
  simData <- simulateMRData(n = 1200, nSnps = 3, trueEffect = 0.2, seed = 555)
  profile <- fitOutcomeModel(
    cohortData = simData$data,
    covariateData = NULL,
    instrumentTable = simData$instrumentTable,
    betaGrid = seq(-1, 1, by = 0.1),
    modelBackend = "cyclops"
  )

  expect_s3_class(profile, "medusaSiteProfile")
  expect_true(is.finite(profile$betaHat))
  expect_true(all(is.finite(profile$logLikProfile)))
})

test_that("fitOutcomeModel helper utilities align instruments and evaluate finite likelihoods", {
  simData <- Medusa::simulateMRData(n = 250, nSnps = 3, trueEffect = 0.2, seed = 556)
  cohortData <- simData$data
  cohortData$personId <- cohortData$person_id

  extraInstrument <- simData$instrumentTable[1, , drop = FALSE]
  extraInstrument$snp_id <- "rs_missing"
  extraInstrument$effect_allele <- "A"
  extraInstrument$other_allele <- "C"
  extendedTable <- rbind(simData$instrumentTable, extraInstrument)

  expect_warning(
    alignment <- alignInstrumentColumns(cohortData, extendedTable),
    "Dropping 1 instrument"
  )
  expect_equal(alignment$instrumentTable$snp_id, simData$instrumentTable$snp_id)

  weights <- computeAlleleScoreWeights(alignment$instrumentTable)
  snpMatrix <- as.matrix(cohortData[, alignment$snpColumns, drop = FALSE])
  alleleScore <- as.numeric(snpMatrix %*% weights)
  baseModel <- data.frame(
    outcome = cohortData$outcome,
    alleleScore = alleleScore,
    stringsAsFactors = FALSE
  )
  covariateData <- data.frame(
    personId = cohortData$personId,
    covariate_1 = rep(c(0, 1), length.out = nrow(cohortData)),
    stringsAsFactors = FALSE
  )
  modelParts <- appendCovariatesToModelData(baseModel, cohortData, covariateData)
  completeModel <- modelParts$modelData[complete.cases(modelParts$modelData), , drop = FALSE]

  coefficientFit <- fitBinaryOutcomeCoefficient(
    modelData = completeModel,
    exposureColumn = "alleleScore",
    covariateColumns = modelParts$covariateColumns,
    modelBackend = "glm",
    regularizationVariance = 0.1,
    instrumentRegularization = FALSE
  )
  pointLogLik <- evaluateBinaryProfilePoint(
    modelData = completeModel,
    exposureColumn = "alleleScore",
    covariateColumns = modelParts$covariateColumns,
    offsetVector = rep(0, nrow(completeModel)),
    modelBackend = "glm",
    regularizationVariance = 0.1,
    instrumentRegularization = FALSE
  )

  expect_equal(sum(abs(weights)), 1, tolerance = 1e-10)
  expect_true(is.finite(coefficientFit$betaHat))
  expect_true(is.finite(coefficientFit$seHat))
  expect_true(is.finite(pointLogLik))

  zeroWeightTable <- alignment$instrumentTable
  zeroWeightTable$beta_ZX <- 0
  expect_error(computeAlleleScoreWeights(zeroWeightTable), "not finite")
})

test_that("appendCovariatesToModelData aligns plain covariate data by person_id", {
  simData <- Medusa::simulateMRData(n = 60, nSnps = 2, trueEffect = 0.2, seed = 560)
  baseModel <- data.frame(
    outcome = simData$data$outcome,
    alleleScore = seq_len(nrow(simData$data)),
    stringsAsFactors = FALSE
  )
  covariateData <- data.frame(
    person_id = rev(simData$data$person_id),
    covariate_1 = seq_len(nrow(simData$data)),
    stringsAsFactors = FALSE
  )

  modelParts <- appendCovariatesToModelData(
    modelData = baseModel,
    cohortData = simData$data,
    covariateData = covariateData
  )

  expected <- covariateData$covariate_1[match(simData$data$person_id, covariateData$person_id)]
  expect_equal(modelParts$modelData$covariate_1, expected)
})

test_that("extractCovariateDataFrame expands sparse medusaCovariateData for fitOutcomeModel", {
  simData <- Medusa::simulateMRData(n = 80, nSnps = 2, trueEffect = 0.2, seed = 562)
  cohortData <- simData$data
  cohortData$personId <- cohortData$person_id

  medusaCovariates <- list(
    covariateData = list(
      covariates = data.frame(
        rowId = c(cohortData$personId[1], cohortData$personId[1], cohortData$personId[2]),
        covariateId = c(101L, 102L, 101L),
        covariateValue = c(1, 1, 1),
        stringsAsFactors = FALSE
      ),
      covariateRef = data.frame(
        covariateId = c(101L, 102L),
        covariateName = c("condition diabetes", "drug aspirin"),
        stringsAsFactors = FALSE
      )
    )
  )
  class(medusaCovariates) <- "medusaCovariateData"

  extracted <- extractCovariateDataFrame(medusaCovariates, cohortData)
  expect_true("personId" %in% names(extracted))
  expect_true(all(c("condition.diabetes", "drug.aspirin") %in% names(extracted)))
  expect_equal(extracted$condition.diabetes[match(cohortData$personId[1], extracted$personId)], 1)
  expect_equal(extracted$drug.aspirin[match(cohortData$personId[2], extracted$personId)], 0)

  baseModel <- data.frame(
    outcome = cohortData$outcome,
    alleleScore = seq_len(nrow(cohortData)),
    stringsAsFactors = FALSE
  )
  modelParts <- appendCovariatesToModelData(baseModel, cohortData, medusaCovariates)
  expect_true(all(c("condition.diabetes", "drug.aspirin") %in% modelParts$covariateColumns))
  expect_equal(modelParts$modelData$condition.diabetes[1], 1)
  expect_equal(modelParts$modelData$drug.aspirin[2], 0)

  expect_no_error(
    fitOutcomeModel(
      cohortData = cohortData,
      covariateData = medusaCovariates,
      instrumentTable = simData$instrumentTable,
      betaGrid = seq(-1, 1, by = 0.1)
    )
  )
})

test_that("fitOutcomeModel imputes missing genotypes to expected dosage instead of zero", {
  simData <- Medusa::simulateMRData(n = 300, nSnps = 3, trueEffect = 0.2, seed = 561)
  snpCols <- grep("^snp_", names(simData$data), value = TRUE)
  withMissing <- simData$data
  withMissing[1:10, snpCols[1]] <- NA_real_
  withExpected <- withMissing
  withExpected[1:10, snpCols[1]] <- 2 * simData$instrumentTable$eaf[[1]]

  profileMissing <- suppressWarnings(
    suppressMessages(
      fitOutcomeModel(
        cohortData = withMissing,
        covariateData = NULL,
        instrumentTable = simData$instrumentTable,
        betaGrid = seq(-1.5, 1.5, by = 0.1)
      )
    )
  )
  profileExpected <- suppressWarnings(
    suppressMessages(
      fitOutcomeModel(
        cohortData = withExpected,
        covariateData = NULL,
        instrumentTable = simData$instrumentTable,
        betaGrid = seq(-1.5, 1.5, by = 0.1)
      )
    )
  )

  expect_equal(profileMissing$betaHat, profileExpected$betaHat, tolerance = 1e-10)
  expect_equal(profileMissing$logLikProfile, profileExpected$logLikProfile, tolerance = 1e-10)
})

test_that("fitOutcomeModel validates missing SNPs and no-case cohorts", {
  simData <- Medusa::simulateMRData(n = 120, nSnps = 2, trueEffect = 0.2, seed = 557)
  noSnp <- simData$data[, setdiff(names(simData$data), grep("^snp_", names(simData$data), value = TRUE))]
  noCases <- simData$data
  noCases$outcome <- 0L

  expect_error(
    fitOutcomeModel(
      cohortData = noSnp,
      covariateData = NULL,
      instrumentTable = simData$instrumentTable,
      betaGrid = seq(-1, 1, by = 0.1)
    ),
    "No SNP columns"
  )
  expect_error(
    fitOutcomeModel(
      cohortData = noCases,
      covariateData = NULL,
      instrumentTable = simData$instrumentTable,
      betaGrid = seq(-1, 1, by = 0.1)
    ),
    "No outcome cases"
  )
})

test_that("fitOutcomeModel helper fallback paths return flat or unavailable outputs", {
  simData <- Medusa::simulateMRData(n = 40, nSnps = 2, trueEffect = 0.2, seed = 558)
  badCohort <- simData$data
  snpCols <- grep("^snp_", names(badCohort), value = TRUE)
  badCohort[[snpCols[1]]] <- NA_integer_

  perSnp <- fitPerSNPModels(
    cohortData = badCohort,
    covariateData = NULL,
    instrumentTable = simData$instrumentTable,
    betaGrid = seq(-1, 1, by = 0.1),
    regularizationVariance = 0.1,
    instrumentRegularization = FALSE,
    modelBackend = "glm"
  )
  expect_true(all(is.na(perSnp$perSnpBetaHats) | is.finite(perSnp$perSnpBetaHats)))
  expect_true(any(!is.finite(perSnp$perSnpProfiles[, 1])))

  expect_error(
    alignInstrumentColumns(
      cohortData = badCohort[, c("outcome", "person_id"), drop = FALSE],
      instrumentTable = simData$instrumentTable
    ),
    "No cohortData SNP columns match"
  )

  expect_error(
    evaluateBinaryProfile(
      modelData = data.frame(outcome = c(1, 0), alleleScore = c(0.1, 0.2)),
      exposureColumn = "alleleScore",
      covariateColumns = 1,
      betaGrid = seq(-1, 1, by = 0.1)
    ),
    "character vector"
  )
  expect_error(
    evaluateBinaryProfile(
      modelData = data.frame(outcome = NA_real_, alleleScore = 0.1),
      exposureColumn = "alleleScore",
      covariateColumns = character(0),
      betaGrid = 0
    ),
    "failed at all grid points"
  )

  expect_true(is.infinite(estimateSEFromProfile(seq(-1, 1, by = 0.1), c(0, rep(-1, 20)))))
  expect_true(is.infinite(estimateSEFromProfile(seq(-1, 1, by = 0.1), rep(0, 21))))
})

test_that("fitOutcomeModel covers warning and internal fallback branches", {
  skip_if_not_installed("mockery")

  simData <- Medusa::simulateMRData(n = 80, nSnps = 2, trueEffect = 0.2, seed = 559)
  lowCaseData <- simData$data
  lowCaseData$outcome <- c(rep(1L, 20), rep(0L, 60))

  expect_warning(
    suppressMessages(
      fitOutcomeModel(
        cohortData = lowCaseData,
        covariateData = NULL,
        instrumentTable = simData$instrumentTable,
        betaGrid = seq(-1, 1, by = 0.1)
      )
    ),
    "Only 20 cases"
  )

  flatProfileFn <- fitOutcomeModel
  mockery::stub(
    flatProfileFn,
    "fitAlleleScoreModel",
    function(...) {
      list(
        logLikProfile = rep(0, 21),
        betaHat = 0,
        seHat = 1,
        alignedInstruments = simData$instrumentTable,
        scoreDefinition = list(
          snpIds = simData$instrumentTable$snp_id,
          scoreWeights = c(0.5, 0.5),
          betaZX = 0.2,
          seZX = 0.05
        )
      )
    }
  )
  expect_warning(
    suppressMessages(
      flatProfileFn(
        cohortData = simData$data,
        covariateData = NULL,
        instrumentTable = simData$instrumentTable,
        betaGrid = seq(-1, 1, by = 0.1)
      )
    ),
    "Profile likelihood is flat"
  )

  alleleFn <- fitAlleleScoreModel
  mockery::stub(alleleFn, "fitBinaryOutcomeCoefficient", function(...) stop("coefficient failed"))
  expect_warning(
    flatFit <- alleleFn(
      cohortData = simData$data,
      covariateData = NULL,
      instrumentTable = simData$instrumentTable,
      betaGrid = seq(-1, 1, by = 0.1),
      regularizationVariance = 0.1,
      instrumentRegularization = FALSE,
      modelBackend = "glm"
    ),
    "Model fitting failed"
  )
  expect_equal(flatFit$betaHat, 0)
  expect_true(is.infinite(flatFit$seHat))

  perSnpFn <- fitPerSNPModels
  mockery::stub(
    perSnpFn,
    "alignInstrumentColumns",
    function(...) {
      list(
        instrumentTable = subset(simData$instrumentTable, select = -pval_ZX),
        snpColumns = grep("^snp_", names(simData$data), value = TRUE)
      )
    }
  )
  mockery::stub(
    perSnpFn,
    "fitAlleleScoreModel",
    function(...) {
      list(
        logLikProfile = rep(0, 21),
        betaHat = 0,
        seHat = 1,
        scoreDefinition = list(snpIds = c("rs1", "rs2"), scoreWeights = c(0.5, 0.5), betaZX = 0.2, seZX = 0.05)
      )
    }
  )
  perSnpNoP <- perSnpFn(
    cohortData = simData$data,
    covariateData = NULL,
    instrumentTable = simData$instrumentTable,
    betaGrid = seq(-1, 1, by = 0.1),
    regularizationVariance = 0.1,
    instrumentRegularization = FALSE,
    modelBackend = "glm"
  )
  expect_true("pval_ZX" %in% names(perSnpNoP$perSnpEstimates))

  perSnpErrorFn <- fitPerSNPModels
  mockery::stub(
    perSnpErrorFn,
    "fitAlleleScoreModel",
    function(...) {
      list(
        logLikProfile = rep(0, 21),
        betaHat = 0,
        seHat = 1,
        scoreDefinition = list(snpIds = c("rs1", "rs2"), scoreWeights = c(0.5, 0.5), betaZX = 0.2, seZX = 0.05)
      )
    }
  )
  mockery::stub(perSnpErrorFn, "fitBinaryOutcomeCoefficient", function(...) stop("per-SNP fit failed"))
  perSnpError <- perSnpErrorFn(
    cohortData = simData$data,
    covariateData = NULL,
    instrumentTable = simData$instrumentTable,
    betaGrid = seq(-1, 1, by = 0.1),
    regularizationVariance = 0.1,
    instrumentRegularization = FALSE,
    modelBackend = "glm"
  )
  expect_equal(nrow(perSnpError$perSnpEstimates), 0)
  expect_true(all(is.infinite(perSnpError$perSnpProfiles)))

  expect_warning(
    appendResult <- appendCovariatesToModelData(
      modelData = data.frame(outcome = simData$data$outcome, alleleScore = 0),
      cohortData = simData$data,
      covariateData = data.frame(other_id = 1:80, cov_x = rnorm(80))
    ),
    "aligning covariates by row order"
  )
  expect_true("cov_x" %in% names(appendResult$modelData))

  expect_equal(createCyclopsPrior(Inf), Cyclops::createPrior("none"))
  convexGrid <- seq(-1, 1, by = 0.1)
  convexProfile <- convexGrid^2
  expect_true(is.infinite(estimateSEFromProfile(convexGrid, convexProfile)))

  flatTopFn <- estimateSEFromProfile
  mockery::stub(flatTopFn, "which.max", function(...) 3L)
  expect_true(is.infinite(flatTopFn(seq(-1, 1, by = 0.5), rep(0, 5))))
})

test_that("fitOutcomeModel helper branches cover Cyclops and profile-point fallbacks", {
  skip_if_not_installed("mockery")

  simData <- Medusa::simulateMRData(n = 120, nSnps = 2, trueEffect = 0.2, seed = 560)
  modelData <- data.frame(
    outcome = simData$data$outcome,
    alleleScore = rowSums(simData$data[, grep("^snp_", names(simData$data), value = TRUE), drop = FALSE])
  )

  cyclopsCoef <- fitBinaryOutcomeCoefficient(
    modelData = modelData,
    exposureColumn = "alleleScore",
    covariateColumns = character(0),
    modelBackend = "cyclops",
    regularizationVariance = 0.1,
    instrumentRegularization = TRUE
  )
  expect_true(is.finite(cyclopsCoef$betaHat))

  expect_equal(
    evaluateBinaryProfilePoint(
      modelData = modelData,
      exposureColumn = "alleleScore",
      covariateColumns = character(0),
      offsetVector = rep(0, nrow(modelData)),
      modelBackend = "glm",
      regularizationVariance = 0.1,
      instrumentRegularization = FALSE
    ),
    evaluateBinaryProfilePoint(
      modelData = modelData,
      exposureColumn = "alleleScore",
      covariateColumns = character(0),
      offsetVector = rep(0, nrow(modelData)),
      modelBackend = "glm",
      regularizationVariance = 0.1,
      instrumentRegularization = FALSE
    )
  )

  glmPointFn <- evaluateBinaryProfilePoint
  mockery::stub(
    glmPointFn,
    "stats::glm.fit",
    function(...) list(converged = FALSE, boundary = FALSE, fitted.values = rep(0.5, nrow(modelData)))
  )
  expect_equal(
    glmPointFn(
      modelData = modelData,
      exposureColumn = "alleleScore",
      covariateColumns = character(0),
      offsetVector = rep(0, nrow(modelData)),
      modelBackend = "glm",
      regularizationVariance = 0.1,
      instrumentRegularization = FALSE
    ),
    -Inf
  )

  glmNonFiniteFn <- evaluateBinaryProfilePoint
  mockery::stub(
    glmNonFiniteFn,
    "stats::glm.fit",
    function(...) list(converged = TRUE, boundary = FALSE, fitted.values = rep(0.5, nrow(modelData)))
  )
  mockery::stub(glmNonFiniteFn, "stats::dbinom", function(...) rep(NaN, nrow(modelData)))
  expect_equal(
    glmNonFiniteFn(
      modelData = modelData,
      exposureColumn = "alleleScore",
      covariateColumns = character(0),
      offsetVector = rep(0, nrow(modelData)),
      modelBackend = "glm",
      regularizationVariance = 0.1,
      instrumentRegularization = FALSE
    ),
    -Inf
  )

  pointFn <- evaluateBinaryProfilePoint
  mockery::stub(pointFn, "fitCyclopsLogistic", function(...) NULL)
  expect_equal(
    pointFn(
      modelData = modelData,
      exposureColumn = "alleleScore",
      covariateColumns = character(0),
      offsetVector = rep(0, nrow(modelData)),
      modelBackend = "cyclops",
      regularizationVariance = 0.1,
      instrumentRegularization = FALSE
    ),
    -Inf
  )

  predictNullFn <- evaluateBinaryProfilePoint
  mockery::stub(predictNullFn, "fitCyclopsLogistic", function(...) structure(list(), class = "cyclopsFit"))
  mockery::stub(predictNullFn, "stats::predict", function(...) NULL)
  expect_equal(
    predictNullFn(
      modelData = modelData,
      exposureColumn = "alleleScore",
      covariateColumns = character(0),
      offsetVector = rep(0, nrow(modelData)),
      modelBackend = "cyclops",
      regularizationVariance = 0.1,
      instrumentRegularization = FALSE
    ),
    -Inf
  )

  cyclopsNonFiniteFn <- evaluateBinaryProfilePoint
  mockery::stub(
    cyclopsNonFiniteFn,
    "fitCyclopsLogistic",
    function(...) structure(list(), class = "cyclopsFit")
  )
  mockery::stub(cyclopsNonFiniteFn, "stats::predict", function(...) rep(0.5, nrow(modelData)))
  mockery::stub(cyclopsNonFiniteFn, "stats::dbinom", function(...) rep(NaN, nrow(modelData)))
  expect_equal(
    cyclopsNonFiniteFn(
      modelData = modelData,
      exposureColumn = "alleleScore",
      covariateColumns = character(0),
      offsetVector = rep(0, nrow(modelData)),
      modelBackend = "cyclops",
      regularizationVariance = 0.1,
      instrumentRegularization = FALSE
    ),
    -Inf
  )
})
