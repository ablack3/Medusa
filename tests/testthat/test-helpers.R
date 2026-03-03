test_that("helper validation and formatting utilities behave as expected", {
  instrumentTable <- data.frame(
    snp_id = c("rs1", "rs2"),
    effect_allele = c("A", "G"),
    other_allele = c("C", "T"),
    beta_ZX = c(0.5, 0.3),
    se_ZX = c(0.05, 0.08),
    pval_ZX = c(1e-10, 1e-5),
    eaf = c(0.3, 0.45),
    stringsAsFactors = FALSE
  )
  siteProfile <- list(
    siteId = "site_1",
    betaGrid = seq(-1, 1, by = 0.1),
    logLikProfile = rep(0, 21),
    nCases = 10,
    nControls = 90,
    snpIds = c("rs1", "rs2"),
    scoreDefinition = list(
      snpIds = c("rs1", "rs2"),
      scoreWeights = c(0.6, 0.4),
      betaZX = 0.42,
      seZX = 0.05
    )
  )

  expect_true(isTRUE(validateInstrumentTable(instrumentTable)))
  expect_true(isTRUE(validateBetaGrid(seq(-1, 1, by = 0.1))))
  expect_true(isTRUE(validateSiteProfile(siteProfile)))
  expect_equal(computeApproxFStatistic(c(0.5, 0.3), c(0.05, 0.1)), c(100, 9))
  expect_equal(makeSnpColumnName(c("rs1", "rs-2")), c("snp_rs1", "snp_rs_2"))
  expect_equal(isStrandAmbiguous(c("A", "G", "A"), c("T", "C", "C")), c(TRUE, TRUE, FALSE))

  expect_error(validateInstrumentTable(data.frame(snp_id = "rs1")), "missing required columns")
  expect_error(validateBetaGrid(c(seq(-1, 0.8, by = 0.2), -0.9)), "sorted")
  badProfile <- siteProfile
  badProfile$logLikProfile <- rep(0, 20)
  expect_error(validateSiteProfile(badProfile), "same length")
  badProfile2 <- siteProfile
  badProfile2$scoreDefinition$scoreWeights <- 0.6
  expect_error(validateSiteProfile(badProfile2), "must have the same length")
})

test_that("mrTheme returns a usable object", {
  themeObject <- mrTheme()

  expect_s3_class(themeObject, "theme")
})

test_that("loadRenderTranslateSql renders and translates package SQL", {
  skip_if_not_installed("SqlRender")

  sql <- loadRenderTranslateSql(
    sqlFileName = "extractOutcomeCohort.sql",
    dbms = "postgresql",
    cdm_database_schema = "cdm",
    cohort_database_schema = "results",
    cohort_table = "cohort",
    outcome_cohort_id = 1L,
    washout_days = 365L,
    exclude_prior_outcome = 1L
  )

  expect_type(sql, "character")
  expect_true(nchar(sql) > 0)
  expect_match(sql, "cdm")
  expect_false(grepl("@cdm_database_schema", sql, fixed = TRUE))
})

test_that("loadRenderTranslateSql validates dependencies and SQL file presence", {
  skip_if_not_installed("mockery")

  helperFn <- loadRenderTranslateSql
  mockery::stub(helperFn, "requireNamespace", function(package, ...) {
    if (identical(package, "SqlRender")) {
      FALSE
    } else {
      base::requireNamespace(package, quietly = TRUE)
    }
  })

  expect_error(
    helperFn("extractOutcomeCohort.sql", dbms = "postgresql"),
    "Package 'SqlRender' is required"
  )

  expect_error(
    loadRenderTranslateSql("definitely_missing.sql", dbms = "postgresql"),
    "SQL file not found"
  )
})
