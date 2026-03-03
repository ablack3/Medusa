test_that("allele harmonization flips genotypes when alleles are swapped", {
  # Test that when effect allele doesn't match, genotypes get flipped
  instruments <- data.frame(
    snp_id = "rs1",
    effect_allele = "A",
    other_allele = "C",
    beta_ZX = 0.5,
    se_ZX = 0.05,
    pval_ZX = 1e-10,
    eaf = 0.3,
    stringsAsFactors = FALSE
  )
  genoAlleles <- data.frame(
    snp_id = "rs1",
    allele_coded = "C",
    allele_noncoded = "A",
    stringsAsFactors = FALSE
  )

  result <- harmonizeAlleles(instruments, genoAlleles)
  # After harmonization, beta should be flipped
 expect_equal(result$instrumentTable$beta_ZX, -0.5)
  expect_equal(result$instrumentTable$eaf, 0.7)
})

test_that("missing genotypes are coded as NA not 0", {
  genotypeData <- data.frame(
    personId = c(1, 1, 2),
    snpId = c("rs1", "rs2", "rs1"),
    genotype = c(2, 1, 0),
    stringsAsFactors = FALSE
  )
  instrumentTable <- data.frame(
    snp_id = c("rs1", "rs2"),
    effect_allele = c("A", "G"),
    other_allele = c("C", "T"),
    beta_ZX = c(0.5, 0.3),
    se_ZX = c(0.05, 0.08),
    pval_ZX = c(1e-10, 1e-5),
    eaf = c(0.3, 0.5),
    stringsAsFactors = FALSE
  )

  result <- reshapeGenotypes(genotypeData, instrumentTable)

  # Person 2 is missing rs2, should be NA
  expect_true(is.na(result$snp_rs2[result$personId == 2]))
  # Person 1 should have values
  expect_equal(result$snp_rs1[result$personId == 1], 2L)
  expect_equal(result$snp_rs2[result$personId == 1], 1L)
})

test_that("reshapeGenotypes validates genotype values", {
  genotypeData <- data.frame(
    personId = c(1, 1),
    snpId = c("rs1", "rs2"),
    genotype = c(3, 1),  # 3 is invalid
    stringsAsFactors = FALSE
  )
  instrumentTable <- data.frame(
    snp_id = c("rs1", "rs2"),
    effect_allele = c("A", "G"),
    other_allele = c("C", "T"),
    beta_ZX = c(0.5, 0.3),
    se_ZX = c(0.05, 0.08),
    pval_ZX = c(1e-10, 1e-5),
    eaf = c(0.3, 0.5),
    stringsAsFactors = FALSE
  )

  expect_warning(
    result <- reshapeGenotypes(genotypeData, instrumentTable),
    "invalid genotype values"
  )
  expect_true(is.na(result$snp_rs1[1]))
})

test_that("reshapeGenotypes produces correct wide format", {
  genotypeData <- data.frame(
    personId = c(1, 1, 2, 2, 3, 3),
    snpId = c("rs1", "rs2", "rs1", "rs2", "rs1", "rs2"),
    genotype = c(0, 1, 2, 0, 1, 2),
    stringsAsFactors = FALSE
  )
  instrumentTable <- data.frame(
    snp_id = c("rs1", "rs2"),
    effect_allele = c("A", "G"),
    other_allele = c("C", "T"),
    beta_ZX = c(0.5, 0.3),
    se_ZX = c(0.05, 0.08),
    pval_ZX = c(1e-10, 1e-5),
    eaf = c(0.3, 0.5),
    stringsAsFactors = FALSE
  )

  result <- reshapeGenotypes(genotypeData, instrumentTable)

  expect_equal(nrow(result), 3)
  expect_true("snp_rs1" %in% names(result))
  expect_true("snp_rs2" %in% names(result))
  expect_equal(result$snp_rs1, c(0L, 2L, 1L))
  expect_equal(result$snp_rs2, c(1L, 0L, 2L))
})

test_that("convertGenotypeString handles VCF-style genotypes", {
  expect_equal(
    convertGenotypeString(c("0/0", "0/1", "1/1", "0|0", "0|1", "1|1")),
    c(0L, 1L, 2L, 0L, 1L, 2L)
  )
  expect_equal(
    convertGenotypeString(c("1/0", "1|0")),
    c(1L, 1L)
  )
})

test_that("convertGenotypeString handles plain integer strings", {
  expect_equal(
    convertGenotypeString(c("0", "1", "2")),
    c(0L, 1L, 2L)
  )
})

test_that("convertGenotypeString warns on unrecognized values", {
  expect_warning(
    result <- convertGenotypeString(c("0/0", "unknown", "1/1")),
    "unrecognized genotype values"
  )
  expect_equal(result, c(0L, NA_integer_, 2L))
})

test_that("convertGenotypeString handles NA values", {
  result <- convertGenotypeString(c("0/1", NA_character_, "1/1"))
  expect_equal(result, c(1L, NA_integer_, 2L))
})

test_that("computeCohortAlleleFrequencies computes correct frequencies", {
  genotypeData <- data.frame(
    snpId = c("rs1", "rs1", "rs1", "rs1", "rs2", "rs2"),
    genotype = c(0L, 1L, 2L, 1L, 0L, 2L),
    stringsAsFactors = FALSE
  )
  freqs <- computeCohortAlleleFrequencies(genotypeData)
  expect_equal(freqs[["rs1"]], mean(c(0, 1, 2, 1)) / 2)
  expect_equal(freqs[["rs2"]], mean(c(0, 2)) / 2)
})

test_that("buildMRCohort validates inputs", {
  expect_error(
    buildMRCohort(
      connectionDetails = "not_a_connection",
      cdmDatabaseSchema = "cdm",
      cohortDatabaseSchema = "results",
      cohortTable = "cohort",
      outcomeCohortId = 1,
      instrumentTable = data.frame(),
      genomicDatabaseSchema = "genomics"
    ),
    "connectionDetails"
  )
})

test_that("buildMRCohort fails cleanly when required dependencies are missing", {
  skip_if_not_installed("mockery")

  dbMissingFn <- buildMRCohort
  mockery::stub(dbMissingFn, "requireNamespace", function(package, ...) {
    if (identical(package, "DatabaseConnector")) {
      FALSE
    } else {
      TRUE
    }
  })

  expect_error(
    dbMissingFn(
      connectionDetails = structure(list(dbms = "postgresql"), class = "connectionDetails"),
      cdmDatabaseSchema = "cdm",
      cohortDatabaseSchema = "results",
      cohortTable = "cohort",
      outcomeCohortId = 1L,
      instrumentTable = data.frame(
        snp_id = "rs1",
        effect_allele = "A",
        other_allele = "C",
        beta_ZX = 0.5,
        se_ZX = 0.05,
        pval_ZX = 1e-10,
        eaf = 0.3,
        stringsAsFactors = FALSE
      ),
      genomicDatabaseSchema = "genomics"
    ),
    "Package 'DatabaseConnector' is required"
  )

  sqlMissingFn <- buildMRCohort
  mockery::stub(sqlMissingFn, "requireNamespace", function(package, ...) {
    if (identical(package, "SqlRender")) {
      FALSE
    } else {
      TRUE
    }
  })

  expect_error(
    sqlMissingFn(
      connectionDetails = structure(list(dbms = "postgresql"), class = "connectionDetails"),
      cdmDatabaseSchema = "cdm",
      cohortDatabaseSchema = "results",
      cohortTable = "cohort",
      outcomeCohortId = 1L,
      instrumentTable = data.frame(
        snp_id = "rs1",
        effect_allele = "A",
        other_allele = "C",
        beta_ZX = 0.5,
        se_ZX = 0.05,
        pval_ZX = 1e-10,
        eaf = 0.3,
        stringsAsFactors = FALSE
      ),
      genomicDatabaseSchema = "genomics"
    ),
    "Package 'SqlRender' is required"
  )
})

test_that("buildMRCohort executes the main extraction flow with stubbed database calls", {
  skip_if_not_installed("mockery")

  buildFn <- buildMRCohort
  connectionDetails <- structure(list(dbms = "postgresql"), class = "connectionDetails")
  fakeConnection <- structure(list(), class = "fakeConnection")
  instrumentTable <- data.frame(
    snp_id = c("rs1", "rs2"),
    effect_allele = c("A", "G"),
    other_allele = c("C", "T"),
    beta_ZX = c(0.5, 0.3),
    se_ZX = c(0.05, 0.08),
    pval_ZX = c(1e-10, 1e-5),
    eaf = c(0.3, 0.5),
    stringsAsFactors = FALSE
  )
  cohortFrame <- data.frame(
    personId = seq_len(100),
    outcome = c(rep(1L, 60), rep(0L, 40)),
    stringsAsFactors = FALSE
  )
  # VARIANT_OCCURRENCE-style genotype data with allele info
  genotypeFrame <- rbind(
    data.frame(
      personId = seq_len(100),
      snpId = "rs1",
      genotypeRaw = rep(c("0/0", "0/1", "1/1", "0/1"), length.out = 100),
      referenceAllele = "C",
      alternateAllele = "A",
      stringsAsFactors = FALSE
    ),
    data.frame(
      personId = seq_len(80),
      snpId = "rs2",
      genotypeRaw = rep(c("0/0", "1/1"), length.out = 80),
      referenceAllele = "T",
      alternateAllele = "G",
      stringsAsFactors = FALSE
    )
  )
  queryCount <- 0L

  mockery::stub(buildFn, "DatabaseConnector::connect", function(...) fakeConnection)
  mockery::stub(buildFn, "DatabaseConnector::disconnect", function(...) invisible(NULL))
  mockery::stub(buildFn, "loadRenderTranslateSql", function(...) "SELECT 1")
  mockery::stub(buildFn, "DatabaseConnector::executeSql", function(...) invisible(NULL))
  mockery::stub(
    buildFn,
    "DatabaseConnector::querySql",
    function(...) {
      queryCount <<- queryCount + 1L
      if (queryCount == 1L) {
        return(cohortFrame)
      }
      genotypeFrame
    }
  )

  expect_warning(
    result <- suppressMessages(
      buildFn(
        connectionDetails = connectionDetails,
        cdmDatabaseSchema = "cdm",
        cohortDatabaseSchema = "results",
        cohortTable = "cohort",
        outcomeCohortId = 1L,
        instrumentTable = instrumentTable,
        genomicDatabaseSchema = "genomics"
      )
    ),
    "20.0% missing genotypes"
  )

  expect_equal(nrow(result), 100)
  expect_true(all(c("snp_rs1", "snp_rs2") %in% names(result)))
  expect_true(any(is.na(result$snp_rs2)))
})

test_that("buildMRCohort stops when the cohort query returns no persons", {
  skip_if_not_installed("mockery")

  buildFn <- buildMRCohort
  connectionDetails <- structure(list(dbms = "postgresql"), class = "connectionDetails")
  instrumentTable <- data.frame(
    snp_id = "rs1",
    effect_allele = "A",
    other_allele = "C",
    beta_ZX = 0.5,
    se_ZX = 0.05,
    pval_ZX = 1e-10,
    eaf = 0.3,
    stringsAsFactors = FALSE
  )

  mockery::stub(buildFn, "DatabaseConnector::connect", function(...) structure(list(), class = "fakeConnection"))
  mockery::stub(buildFn, "DatabaseConnector::disconnect", function(...) invisible(NULL))
  mockery::stub(buildFn, "loadRenderTranslateSql", function(...) "SELECT 1")
  mockery::stub(buildFn, "DatabaseConnector::executeSql", function(...) invisible(NULL))
  mockery::stub(buildFn, "DatabaseConnector::querySql", function(...) data.frame())

  expect_error(
    suppressMessages(
      buildFn(
        connectionDetails = connectionDetails,
        cdmDatabaseSchema = "cdm",
        cohortDatabaseSchema = "results",
        cohortTable = "cohort",
        outcomeCohortId = 1L,
        instrumentTable = instrumentTable,
        genomicDatabaseSchema = "genomics"
      )
    ),
    "No persons found"
  )
})

test_that("buildMRCohort stops when genotype data are unavailable", {
  skip_if_not_installed("mockery")

  buildFn <- buildMRCohort
  connectionDetails <- structure(list(dbms = "postgresql"), class = "connectionDetails")
  instrumentTable <- data.frame(
    snp_id = "rs1",
    effect_allele = "A",
    other_allele = "C",
    beta_ZX = 0.5,
    se_ZX = 0.05,
    pval_ZX = 1e-10,
    eaf = 0.3,
    stringsAsFactors = FALSE
  )
  queryCount <- 0L

  mockery::stub(buildFn, "DatabaseConnector::connect", function(...) structure(list(), class = "fakeConnection"))
  mockery::stub(buildFn, "DatabaseConnector::disconnect", function(...) invisible(NULL))
  mockery::stub(buildFn, "loadRenderTranslateSql", function(...) "SELECT 1")
  mockery::stub(buildFn, "DatabaseConnector::executeSql", function(...) invisible(NULL))
  mockery::stub(
    buildFn,
    "DatabaseConnector::querySql",
    function(...) {
      queryCount <<- queryCount + 1L
      if (queryCount == 1L) {
        return(data.frame(personId = 1:5, outcome = c(1L, 1L, 1L, 0L, 0L)))
      }
      data.frame()
    }
  )

  expect_error(
    suppressMessages(
      suppressWarnings(
        buildFn(
          connectionDetails = connectionDetails,
          cdmDatabaseSchema = "cdm",
          cohortDatabaseSchema = "results",
          cohortTable = "cohort",
          outcomeCohortId = 1L,
          instrumentTable = instrumentTable,
          genomicDatabaseSchema = "genomics"
        )
      )
    ),
    "No persons in cohort have genotype data"
  )
})
