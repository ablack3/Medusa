test_that("getMRInstruments returns expected columns with mocked API", {
  mockAssociations <- data.frame(
    rsid = c("rs1", "rs2", "rs3"),
    ea = c("A", "G", "T"),
    nea = c("G", "C", "A"),
    beta = c(0.5, 0.3, 0.2),
    se = c(0.05, 0.08, 0.03),
    p = c(1e-20, 1e-10, 1e-15),
    eaf = c(0.3, 0.5, 0.2),
    id = rep("ieu-a-1119", 3),
    stringsAsFactors = FALSE
  )

  mockClumped <- data.frame(
    rsid = c("rs1", "rs3"),
    pval = c(1e-20, 1e-15),
    id = rep("ieu-a-1119", 2),
    stringsAsFactors = FALSE
  )

  local_mocked_bindings(
    associations = function(...) mockAssociations,
    ld_clump = function(...) mockClumped,
    .package = "ieugwasr"
  )

  result <- getMRInstruments("ieu-a-1119")

  expectedCols <- c("snp_id", "effect_allele", "other_allele",
                     "beta_ZX", "se_ZX", "pval_ZX", "eaf",
                     "gene_region", "fStatistic", "strandAmbiguous")
  expect_true(all(expectedCols %in% names(result)))
  expect_equal(nrow(result), 2)
})

test_that("strand-ambiguous SNPs are flagged correctly", {
  mockAssociations <- data.frame(
    rsid = c("rs1", "rs2", "rs3"),
    ea = c("A", "G", "A"),
    nea = c("T", "C", "G"),
    beta = c(0.5, 0.3, 0.2),
    se = c(0.05, 0.08, 0.03),
    p = c(1e-20, 1e-10, 1e-15),
    eaf = c(0.3, 0.5, 0.2),
    id = rep("trait1", 3),
    stringsAsFactors = FALSE
  )

  local_mocked_bindings(
    associations = function(...) mockAssociations,
    ld_clump = function(...) {
      data.frame(rsid = c("rs1", "rs2", "rs3"),
                 pval = c(1e-20, 1e-10, 1e-15),
                 id = rep("trait1", 3),
                 stringsAsFactors = FALSE)
    },
    .package = "ieugwasr"
  )

  expect_warning(
    result <- getMRInstruments("trait1"),
    "strand-ambiguous"
  )

  # rs1 is A/T (ambiguous), rs2 is G/C (ambiguous), rs3 is A/G (not ambiguous)
  expect_equal(result$strandAmbiguous[result$snp_id == "rs1"], TRUE)
  expect_equal(result$strandAmbiguous[result$snp_id == "rs2"], TRUE)
  expect_equal(result$strandAmbiguous[result$snp_id == "rs3"], FALSE)
})

test_that("function stops with informative error when zero SNPs returned", {
  local_mocked_bindings(
    associations = function(...) {
      data.frame(
        rsid = character(0), ea = character(0), nea = character(0),
        beta = numeric(0), se = numeric(0), p = numeric(0),
        eaf = numeric(0), id = character(0),
        stringsAsFactors = FALSE
      )
    },
    .package = "ieugwasr"
  )

  expect_error(
    getMRInstruments("bad-trait-id"),
    "No genome-wide significant SNPs found"
  )
})

test_that("function handles API failure gracefully", {
  local_mocked_bindings(
    associations = function(...) stop("Connection timeout"),
    .package = "ieugwasr"
  )

  expect_error(
    getMRInstruments("ieu-a-1119"),
    "Failed to query OpenGWAS"
  )
})

test_that("F-statistic approximation is computed correctly", {
  mockAssociations <- data.frame(
    rsid = c("rs1", "rs2"),
    ea = c("A", "G"),
    nea = c("C", "T"),
    beta = c(0.5, 0.1),
    se = c(0.05, 0.1),
    p = c(1e-20, 1e-10),
    eaf = c(0.3, 0.5),
    id = rep("trait1", 2),
    stringsAsFactors = FALSE
  )

  local_mocked_bindings(
    associations = function(...) mockAssociations,
    ld_clump = function(...) {
      data.frame(rsid = c("rs1", "rs2"),
                 pval = c(1e-20, 1e-10),
                 id = rep("trait1", 2),
                 stringsAsFactors = FALSE)
    },
    .package = "ieugwasr"
  )

  expect_warning(
    result <- getMRInstruments("trait1"),
    "F-statistic < 10"
  )

  # F = (beta/se)^2
  expect_equal(result$fStatistic[result$snp_id == "rs1"], (0.5 / 0.05)^2)
  expect_equal(result$fStatistic[result$snp_id == "rs2"], (0.1 / 0.1)^2)
})

test_that("warning issued when fewer than 3 instruments", {
  mockAssociations <- data.frame(
    rsid = c("rs1", "rs2"),
    ea = c("A", "G"),
    nea = c("C", "T"),
    beta = c(0.5, 0.3),
    se = c(0.05, 0.08),
    p = c(1e-20, 1e-10),
    eaf = c(0.3, 0.5),
    id = rep("trait1", 2),
    stringsAsFactors = FALSE
  )

  local_mocked_bindings(
    associations = function(...) mockAssociations,
    ld_clump = function(...) {
      data.frame(rsid = c("rs1", "rs2"),
                 pval = c(1e-20, 1e-10),
                 id = rep("trait1", 2),
                 stringsAsFactors = FALSE)
    },
    .package = "ieugwasr"
  )

  expect_warning(
    result <- getMRInstruments("trait1"),
    "Only 2 instruments available"
  )
})

test_that("createInstrumentTable produces valid instrument table", {
  result <- createInstrumentTable(
    snpId = c("rs1", "rs2"),
    effectAllele = c("A", "G"),
    otherAllele = c("C", "T"),
    betaZX = c(0.5, 0.3),
    seZX = c(0.05, 0.08),
    pvalZX = c(1e-10, 1e-5),
    eaf = c(0.3, 0.45)
  )

  expect_equal(nrow(result), 2)
  expect_true("snp_id" %in% names(result))
  expect_true("fStatistic" %in% names(result))
  expect_true("strandAmbiguous" %in% names(result))
  expect_equal(result$fStatistic[1], (0.5 / 0.05)^2)
  expect_false(is.null(attr(result, "retrievalTimestamp")))
})

test_that("createInstrumentTable validates inputs", {
  expect_error(
    createInstrumentTable(
      snpId = c("rs1"),
      effectAllele = c("A", "G"),  # length mismatch
      otherAllele = c("C"),
      betaZX = c(0.5),
      seZX = c(0.05),
      pvalZX = c(1e-10),
      eaf = c(0.3)
    )
  )
})

test_that("instrument table carries metadata attributes", {
  mockAssociations <- data.frame(
    rsid = c("rs1", "rs2", "rs3"),
    ea = c("A", "G", "C"),
    nea = c("C", "T", "A"),
    beta = c(0.5, 0.3, 0.4),
    se = c(0.05, 0.08, 0.06),
    p = c(1e-20, 1e-10, 1e-15),
    eaf = c(0.3, 0.5, 0.2),
    id = rep("ieu-a-1119", 3),
    stringsAsFactors = FALSE
  )

  local_mocked_bindings(
    associations = function(...) mockAssociations,
    ld_clump = function(...) {
      data.frame(rsid = c("rs1", "rs2", "rs3"),
                 pval = c(1e-20, 1e-10, 1e-15),
                 id = rep("ieu-a-1119", 3),
                 stringsAsFactors = FALSE)
    },
    .package = "ieugwasr"
  )

  result <- getMRInstruments("ieu-a-1119")

  expect_false(is.null(attr(result, "retrievalTimestamp")))
  expect_equal(attr(result, "exposureTraitId"), "ieu-a-1119")
  expect_equal(attr(result, "parameters")$pThreshold, 5e-8)
})

test_that("additionalSnps are appended when available from OpenGWAS", {
  mainAssociations <- data.frame(
    rsid = c("rs1", "rs2"),
    ea = c("A", "G"),
    nea = c("C", "T"),
    beta = c(0.5, 0.3),
    se = c(0.05, 0.08),
    p = c(1e-20, 1e-10),
    eaf = c(0.3, 0.5),
    id = rep("trait1", 2),
    stringsAsFactors = FALSE
  )
  additionalAssociations <- data.frame(
    rsid = "rs_extra",
    ea = "C",
    nea = "A",
    beta = 0.25,
    se = 0.07,
    p = 1e-6,
    eaf = 0.4,
    id = "trait1",
    gene = "GENE_EXTRA",
    stringsAsFactors = FALSE
  )

  local_mocked_bindings(
    associations = function(variants = NULL, ...) {
      if (is.null(variants)) {
        return(mainAssociations)
      }
      additionalAssociations
    },
    ld_clump = function(...) {
      data.frame(
        rsid = "rs1",
        pval = 1e-20,
        id = "trait1",
        stringsAsFactors = FALSE
      )
    },
    .package = "ieugwasr"
  )

  expect_warning(
    result <- getMRInstruments("trait1", additionalSnps = "rs_extra"),
    "Only 2 instruments available"
  )

  expect_true(all(c("rs1", "rs_extra") %in% result$snp_id))
  expect_equal(result$gene_region[result$snp_id == "rs_extra"], "GENE_EXTRA")
})

test_that("LD clumping failures surface an informative error", {
  mockAssociations <- data.frame(
    rsid = c("rs1", "rs2"),
    ea = c("A", "G"),
    nea = c("C", "T"),
    beta = c(0.5, 0.3),
    se = c(0.05, 0.08),
    p = c(1e-20, 1e-10),
    eaf = c(0.3, 0.5),
    id = rep("trait1", 2),
    stringsAsFactors = FALSE
  )

  local_mocked_bindings(
    associations = function(...) mockAssociations,
    ld_clump = function(...) stop("reference panel unavailable"),
    .package = "ieugwasr"
  )

  expect_error(
    getMRInstruments("trait1"),
    "LD clumping failed"
  )
})

test_that("getMRInstruments handles package-missing, post-filter empty, and post-clumping empty branches", {
  skip_if_not_installed("mockery")

  noPkgFn <- getMRInstruments
  mockery::stub(noPkgFn, "requireNamespace", function(package, ...) {
    if (identical(package, "ieugwasr")) {
      FALSE
    } else {
      TRUE
    }
  })
  expect_error(
    noPkgFn("trait1"),
    "Package 'ieugwasr' is required"
  )

  highPAssociations <- data.frame(
    rsid = c("rs1", "rs2"),
    ea = c("A", "G"),
    nea = c("C", "T"),
    beta = c(0.5, 0.3),
    se = c(0.05, 0.08),
    p = c(1e-3, 1e-4),
    eaf = c(0.3, 0.5),
    id = rep("trait1", 2),
    stringsAsFactors = FALSE
  )
  local_mocked_bindings(
    associations = function(...) highPAssociations,
    .package = "ieugwasr"
  )
  expect_error(
    getMRInstruments("trait1", pThreshold = 1e-8),
    "at p <"
  )

  local_mocked_bindings(
    associations = function(...) {
      data.frame(
        rsid = c("rs1", "rs2"),
        ea = c("A", "G"),
        nea = c("C", "T"),
        beta = c(0.5, 0.3),
        se = c(0.05, 0.08),
        p = c(1e-20, 1e-10),
        eaf = c(0.3, 0.5),
        id = rep("trait1", 2),
        stringsAsFactors = FALSE
      )
    },
    ld_clump = function(...) data.frame(rsid = character(0), pval = numeric(0), id = character(0)),
    .package = "ieugwasr"
  )
  expect_error(
    getMRInstruments("trait1"),
    "No SNPs remained after LD clumping"
  )
})

test_that("getMRInstruments warns when forced-include retrieval fails and createInstrumentTable accepts geneRegion", {
  mainAssociations <- data.frame(
    rsid = c("rs1", "rs2"),
    ea = c("A", "G"),
    nea = c("C", "T"),
    beta = c(0.5, 0.3),
    se = c(0.05, 0.08),
    p = c(1e-20, 1e-10),
    eaf = c(0.3, 0.5),
    id = rep("trait1", 2),
    gene = c("GENE1", "GENE2"),
    stringsAsFactors = FALSE
  )

  local_mocked_bindings(
    associations = function(variants = NULL, ...) {
      if (is.null(variants)) {
        return(mainAssociations)
      }
      stop("api down")
    },
    ld_clump = function(...) {
      data.frame(
        rsid = "rs1",
        pval = 1e-20,
        id = "trait1",
        stringsAsFactors = FALSE
      )
    },
    .package = "ieugwasr"
  )

  expect_warning(
    result <- getMRInstruments("trait1", additionalSnps = "rs_extra"),
    "Failed to retrieve additional SNPs"
  )
  expect_equal(result$gene_region[1], "GENE1")

  created <- createInstrumentTable(
    snpId = c("rs1", "rs2"),
    effectAllele = c("a", "g"),
    otherAllele = c("c", "t"),
    betaZX = c(0.5, 0.3),
    seZX = c(0.05, 0.08),
    pvalZX = c(1e-10, 1e-5),
    eaf = c(0.3, 0.45),
    geneRegion = c("GENE1", "GENE2")
  )

  expect_equal(created$gene_region, c("GENE1", "GENE2"))
  expect_equal(created$effect_allele, c("A", "G"))
})
