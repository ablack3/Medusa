test_that("createDefaultMRCovariateSettings returns a FeatureExtraction settings object", {
  settings <- createDefaultMRCovariateSettings()

  expect_false(is.null(settings))
  expect_true(inherits(settings, "covariateSettings"))
})

test_that("buildMRCovariates assembles covariates with stubbed database calls", {
  skip_if_not_installed("mockery")

  buildFn <- buildMRCovariates
  fakeConnection <- structure(list(), class = "fakeConnection")
  connectionDetails <- structure(list(dbms = "postgresql"), class = "connectionDetails")

  mockery::stub(buildFn, "DatabaseConnector::connect", function(...) fakeConnection)
  mockery::stub(buildFn, "DatabaseConnector::disconnect", function(...) invisible(NULL))
  mockery::stub(
    buildFn,
    "FeatureExtraction::getDbCovariateData",
    function(...) {
      list(
        covariateRef = data.frame(
          covariateId = c(1L, 2L),
          stringsAsFactors = FALSE
        ),
        covariates = data.frame(
          rowId = c(1L, 2L, 1L),
          covariateId = c(1L, 1L, 2L),
          covariateValue = c(1, 0, 1),
          stringsAsFactors = FALSE
        )
      )
    }
  )
  mockery::stub(
    buildFn,
    "SqlRender::translateSql",
    function(...) list(sql = "SELECT person_id, pc1, pc2 FROM ancestry")
  )
  mockery::stub(
    buildFn,
    "DatabaseConnector::querySql",
    function(...) {
      data.frame(
        personId = c(1L, 2L),
        pc1 = c(0.1, -0.2),
        pc2 = c(0.3, 0.4),
        stringsAsFactors = FALSE
      )
    }
  )

  covariates <- suppressMessages(
    buildFn(
      connectionDetails = connectionDetails,
      cdmDatabaseSchema = "cdm",
      cohortDatabaseSchema = "results",
      cohortTable = "cohort",
      outcomeCohortId = 1L,
      covariateSettings = structure(list(dummy = TRUE), class = "covariateSettings"),
      ancestryPCsTable = "ancestry",
      ancestryPCsSchema = "scratch",
      numAncestryPCs = 2L
    )
  )

  expect_s3_class(covariates, "medusaCovariateData")
  expect_equal(nrow(covariates$covariateData$covariates), 3)
  expect_equal(nrow(covariates$ancestryPCs), 2)
})

test_that("buildMRCovariates can create default settings and skip ancestry extraction", {
  skip_if_not_installed("mockery")

  buildFn <- buildMRCovariates
  fakeConnection <- structure(list(), class = "fakeConnection")
  connectionDetails <- structure(list(dbms = "postgresql"), class = "connectionDetails")
  defaultSettings <- structure(list(default = TRUE), class = "covariateSettings")

  mockery::stub(buildFn, "createDefaultMRCovariateSettings", function(...) defaultSettings)
  mockery::stub(buildFn, "DatabaseConnector::connect", function(...) fakeConnection)
  mockery::stub(buildFn, "DatabaseConnector::disconnect", function(...) invisible(NULL))
  mockery::stub(
    buildFn,
    "FeatureExtraction::getDbCovariateData",
    function(...) {
      list(
        covariateRef = data.frame(covariateId = 1L, stringsAsFactors = FALSE),
        covariates = data.frame(rowId = 1L, covariateId = 1L, covariateValue = 1, stringsAsFactors = FALSE)
      )
    }
  )

  covariates <- suppressMessages(
    buildFn(
      connectionDetails = connectionDetails,
      cdmDatabaseSchema = "cdm",
      cohortDatabaseSchema = "results",
      cohortTable = "cohort",
      outcomeCohortId = 1L
    )
  )

  expect_s3_class(covariates, "medusaCovariateData")
  expect_null(covariates$ancestryPCs)
  expect_true(identical(covariates$settings, defaultSettings))
})
