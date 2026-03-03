test_that("non-ambiguous SNP is correctly harmonized without changes", {
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
    allele_coded = "A",
    allele_noncoded = "C",
    stringsAsFactors = FALSE
  )
  result <- harmonizeAlleles(instruments, genoAlleles)
  expect_equal(result$instrumentTable$beta_ZX, 0.5)
  expect_equal(result$instrumentTable$eaf, 0.3)
  expect_false(result$instrumentTable$flipped)
  expect_equal(nrow(result$removedSnps), 0)
})

test_that("allele flip correctly inverts effect direction and EAF", {
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
  expect_equal(result$instrumentTable$beta_ZX, -0.5)
  expect_equal(result$instrumentTable$eaf, 0.7)
  expect_true(result$instrumentTable$flipped)
})

test_that("A/T palindromic SNP is flagged as strand-ambiguous", {
  instruments <- data.frame(
    snp_id = "rs1",
    effect_allele = "A",
    other_allele = "T",
    beta_ZX = 0.5,
    se_ZX = 0.05,
    pval_ZX = 1e-10,
    eaf = 0.48,
    stringsAsFactors = FALSE
  )
  genoAlleles <- data.frame(
    snp_id = "rs1",
    allele_coded = "A",
    allele_noncoded = "T",
    stringsAsFactors = FALSE
  )
  result <- suppressMessages(
    harmonizeAlleles(instruments, genoAlleles, eafPalindromicThreshold = 0.08)
  )
  # EAF 0.48 is within 0.08 of 0.5, so should be removed
  expect_equal(nrow(result$instrumentTable), 0)
  expect_equal(nrow(result$removedSnps), 1)
  expect_equal(result$removedSnps$reason, "palindromic_ambiguous_eaf")
})

test_that("palindromic SNP with EAF far from 0.5 is retained", {
  instruments <- data.frame(
    snp_id = "rs1",
    effect_allele = "A",
    other_allele = "T",
    beta_ZX = 0.5,
    se_ZX = 0.05,
    pval_ZX = 1e-10,
    eaf = 0.2,
    stringsAsFactors = FALSE
  )
  genoAlleles <- data.frame(
    snp_id = "rs1",
    allele_coded = "A",
    allele_noncoded = "T",
    stringsAsFactors = FALSE
  )
  result <- harmonizeAlleles(instruments, genoAlleles, eafPalindromicThreshold = 0.08)
  # EAF 0.2 is far from 0.5, should be kept
  expect_equal(nrow(result$instrumentTable), 1)
  expect_equal(nrow(result$removedSnps), 0)
})

test_that("G/C palindromic SNP with ambiguous EAF is dropped with warning", {
  instruments <- data.frame(
    snp_id = "rs1",
    effect_allele = "G",
    other_allele = "C",
    beta_ZX = 0.3,
    se_ZX = 0.04,
    pval_ZX = 1e-8,
    eaf = 0.50,
    stringsAsFactors = FALSE
  )
  genoAlleles <- data.frame(
    snp_id = "rs1",
    allele_coded = "G",
    allele_noncoded = "C",
    stringsAsFactors = FALSE
  )
  expect_message(
    result <- harmonizeAlleles(instruments, genoAlleles),
    "palindromic"
  )
  expect_equal(nrow(result$instrumentTable), 0)
  expect_equal(result$removedSnps$reason, "palindromic_ambiguous_eaf")
})

test_that("SNP not found in genotype data is removed", {
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
    snp_id = "rs999",
    allele_coded = "A",
    allele_noncoded = "C",
    stringsAsFactors = FALSE
  )
  result <- suppressMessages(harmonizeAlleles(instruments, genoAlleles))
  expect_equal(nrow(result$instrumentTable), 0)
  expect_equal(result$removedSnps$reason, "not_in_genotype_data")
})

test_that("completely mismatched alleles are removed", {
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
    allele_coded = "G",
    allele_noncoded = "T",
    stringsAsFactors = FALSE
  )
  # A/C vs G/T — complement of A is T, complement of C is G
  # So complement match: T==G? No. T==T? Let me think...
  # A->T, C->G, so complement pair is T/G which matches G/T (swapped)
  # This should actually be harmonized via complement + flip
  # Let me use truly mismatched alleles instead
  genoAlleles2 <- data.frame(
    snp_id = "rs1",
    allele_coded = "A",
    allele_noncoded = "A",  # Invalid — same allele
    stringsAsFactors = FALSE
  )
  result <- suppressMessages(harmonizeAlleles(instruments, genoAlleles2))
  expect_equal(nrow(result$instrumentTable), 0)
})

test_that("multiple SNPs harmonized correctly in batch", {
  instruments <- data.frame(
    snp_id = c("rs1", "rs2", "rs3"),
    effect_allele = c("A", "G", "A"),
    other_allele = c("C", "T", "G"),
    beta_ZX = c(0.5, -0.3, 0.2),
    se_ZX = c(0.05, 0.08, 0.06),
    pval_ZX = c(1e-10, 1e-5, 1e-8),
    eaf = c(0.3, 0.45, 0.6),
    stringsAsFactors = FALSE
  )
  genoAlleles <- data.frame(
    snp_id = c("rs1", "rs2", "rs3"),
    allele_coded = c("C", "G", "A"),
    allele_noncoded = c("A", "T", "G"),
    stringsAsFactors = FALSE
  )

  result <- harmonizeAlleles(instruments, genoAlleles)

  # rs1: effect=A, other=C, geno coded=C noncoded=A => flip
  expect_equal(result$instrumentTable$beta_ZX[result$instrumentTable$snp_id == "rs1"], -0.5)

  # rs2: G/T — palindromic? No, G/T is not palindromic. effect=G, geno coded=G => match
  expect_equal(result$instrumentTable$beta_ZX[result$instrumentTable$snp_id == "rs2"], -0.3)

  # rs3: A/G, geno coded=A noncoded=G => direct match
  expect_equal(result$instrumentTable$beta_ZX[result$instrumentTable$snp_id == "rs3"], 0.2)
})

test_that("complement allele match is retained without flipping", {
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
    allele_coded = "T",
    allele_noncoded = "G",
    stringsAsFactors = FALSE
  )

  result <- harmonizeAlleles(instruments, genoAlleles)

  expect_equal(result$instrumentTable$beta_ZX, 0.5)
  expect_equal(result$instrumentTable$eaf, 0.3)
  expect_false(result$instrumentTable$flipped)
  expect_equal(nrow(result$removedSnps), 0)
})

test_that("complement allele swap flips effect direction and eaf", {
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
    allele_coded = "G",
    allele_noncoded = "T",
    stringsAsFactors = FALSE
  )

  result <- harmonizeAlleles(instruments, genoAlleles)

  expect_equal(result$instrumentTable$beta_ZX, -0.5)
  expect_equal(result$instrumentTable$eaf, 0.7)
  expect_equal(result$instrumentTable$effect_allele, "G")
  expect_equal(result$instrumentTable$other_allele, "T")
  expect_true(result$instrumentTable$flipped)
  expect_equal(nrow(result$removedSnps), 0)
})

test_that("harmonizeAlleles is idempotent after SNPs are aligned to genotype coding", {
  instruments <- data.frame(
    snp_id = c("rs1", "rs2", "rs3"),
    effect_allele = c("A", "G", "A"),
    other_allele = c("C", "T", "G"),
    beta_ZX = c(0.5, -0.3, 0.2),
    se_ZX = c(0.05, 0.08, 0.06),
    pval_ZX = c(1e-10, 1e-5, 1e-8),
    eaf = c(0.3, 0.45, 0.6),
    stringsAsFactors = FALSE
  )
  genoAlleles <- data.frame(
    snp_id = c("rs1", "rs2", "rs3"),
    allele_coded = c("C", "G", "A"),
    allele_noncoded = c("A", "T", "G"),
    stringsAsFactors = FALSE
  )

  firstPass <- suppressMessages(harmonizeAlleles(instruments, genoAlleles))
  secondPass <- suppressMessages(harmonizeAlleles(firstPass$instrumentTable, genoAlleles))

  stableCols <- c(
    "snp_id", "effect_allele", "other_allele",
    "beta_ZX", "se_ZX", "pval_ZX", "eaf"
  )
  expect_equal(
    secondPass$instrumentTable[, stableCols, drop = FALSE],
    firstPass$instrumentTable[, stableCols, drop = FALSE]
  )
  expect_true(all(!secondPass$instrumentTable$flipped))
  expect_equal(nrow(secondPass$removedSnps), 0)
})

test_that("complementAllele returns correct complements", {
  expect_equal(complementAllele("A"), "T")
  expect_equal(complementAllele("T"), "A")
  expect_equal(complementAllele("G"), "C")
  expect_equal(complementAllele("C"), "G")
})

test_that("complementAllele handles lowercase", {
  expect_equal(complementAllele("a"), "T")
  expect_equal(complementAllele("g"), "C")
})

test_that("complementAllele errors on invalid allele", {
  expect_error(complementAllele("X"), "Invalid allele")
})

test_that("isStrandAmbiguous correctly identifies AT and GC pairs", {
  expect_true(isStrandAmbiguous("A", "T"))
  expect_true(isStrandAmbiguous("T", "A"))
  expect_true(isStrandAmbiguous("G", "C"))
  expect_true(isStrandAmbiguous("C", "G"))
  expect_false(isStrandAmbiguous("A", "C"))
  expect_false(isStrandAmbiguous("A", "G"))
  expect_false(isStrandAmbiguous("T", "C"))
  expect_false(isStrandAmbiguous("T", "G"))
})
