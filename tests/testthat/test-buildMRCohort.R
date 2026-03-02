test_that("allele harmonization flips genotypes when alleles are swapped", {
  # Test that when effect allele doesn't match, genotypes get flipped (2 - genotype)
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

test_that("buildMRCohort validates inputs", {
  expect_error(
    buildMRCohort(
      connectionDetails = "not_a_connection",
      cdmDatabaseSchema = "cdm",
      cohortDatabaseSchema = "results",
      cohortTable = "cohort",
      outcomeCohortId = 1,
      instrumentTable = data.frame(),
      genomicLinkageSchema = "genomics",
      genomicLinkageTable = "genotype"
    ),
    "connectionDetails"
  )
})
