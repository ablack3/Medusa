# Simulation helpers for Medusa test suite
# These functions generate synthetic data for testing without requiring
# live database connections or GWAS API access.

#' Simulate MR Data with Known Causal Effect
#'
#' Generates a complete dataset with SNP genotypes, an exposure, confounders,
#' and a binary outcome where the true causal effect of the exposure on the
#' outcome is known and recoverable.
#'
#' @param n Number of individuals.
#' @param nSnps Number of SNPs (instruments).
#' @param trueEffect True causal effect of exposure on outcome (log-OR scale).
#' @param confoundingStrength Strength of confounding (coefficient of confounder
#'   on both exposure and outcome).
#' @param snpEffectRange Range of SNP-exposure effect sizes.
#' @param seed Random seed for reproducibility.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{data}{Data frame with person_id, outcome (0/1), snp_1..snp_K,
#'       confounder_1, confounder_2, exposure.}
#'     \item{instrumentTable}{Data frame mimicking getMRInstruments() output.}
#'     \item{trueEffect}{The true causal effect used in simulation.}
#'   }
simulateMRData <- function(n = 5000,
                           nSnps = 10,
                           trueEffect = 0.5,
                           confoundingStrength = 0.3,
                           snpEffectRange = c(0.1, 0.5),
                           seed = 42) {
  set.seed(seed)

  # Generate confounders
  confounder1 <- rnorm(n)
  confounder2 <- rbinom(n, 1, 0.5)

  # Generate SNP genotypes (0, 1, 2) with random MAFs
  mafs <- runif(nSnps, 0.1, 0.4)
  snpMatrix <- matrix(NA_integer_, nrow = n, ncol = nSnps)
  for (j in seq_len(nSnps)) {
    snpMatrix[, j] <- rbinom(n, 2, mafs[j])
  }
  colnames(snpMatrix) <- paste0("snp_", seq_len(nSnps))

  # True SNP-exposure effects
  betaZX <- runif(nSnps, snpEffectRange[1], snpEffectRange[2])
  seZX <- abs(betaZX) / runif(nSnps, 4, 8)  # Reasonable SE

  # Generate exposure: function of SNPs + confounders + noise
  snpComponent <- as.numeric(snpMatrix %*% betaZX)
  exposure <- snpComponent +
    confoundingStrength * confounder1 +
    confoundingStrength * confounder2 +
    rnorm(n)

  # Generate binary outcome: function of exposure + confounders
  logOdds <- trueEffect * exposure +
    confoundingStrength * confounder1 +
    confoundingStrength * confounder2
  probOutcome <- 1 / (1 + exp(-logOdds))
  outcome <- rbinom(n, 1, probOutcome)

  # Build data frame
  df <- data.frame(
    person_id = seq_len(n),
    outcome = outcome,
    snpMatrix,
    confounder_1 = confounder1,
    confounder_2 = confounder2,
    exposure = exposure,
    stringsAsFactors = FALSE
  )

  # Build instrument table
  alleles <- c("A", "C", "G", "T")
  instrumentTable <- data.frame(
    snp_id = paste0("rs", seq_len(nSnps)),
    effect_allele = sample(alleles, nSnps, replace = TRUE),
    other_allele = sample(alleles, nSnps, replace = TRUE),
    beta_ZX = betaZX,
    se_ZX = seZX,
    pval_ZX = 2 * pnorm(-abs(betaZX / seZX)),
    eaf = mafs,
    gene_region = paste0("GENE", seq_len(nSnps)),
    stringsAsFactors = FALSE
  )
  # Ensure effect and other alleles differ
  for (i in seq_len(nSnps)) {
    while (instrumentTable$effect_allele[i] == instrumentTable$other_allele[i]) {
      instrumentTable$other_allele[i] <- sample(alleles, 1)
    }
  }

  list(
    data = df,
    instrumentTable = instrumentTable,
    trueEffect = trueEffect
  )
}


#' Simulate an Instrument Table
#'
#' Creates a synthetic instrument table mimicking the output of getMRInstruments().
#'
#' @param nSnps Number of SNPs to simulate.
#' @param seed Random seed for reproducibility.
#'
#' @return Data frame with columns: snp_id, effect_allele, other_allele,
#'   beta_ZX, se_ZX, pval_ZX, eaf, gene_region.
simulateInstrumentTable <- function(nSnps = 10, seed = 42) {
  set.seed(seed)
  betaZX <- rnorm(nSnps, mean = 0, sd = 0.3)
  seZX <- runif(nSnps, 0.02, 0.08)

  nonAmbiguousPairs <- list(
    c("A", "C"), c("A", "G"), c("T", "C"), c("T", "G"),
    c("C", "A"), c("G", "A"), c("C", "T"), c("G", "T")
  )

  allelePairs <- nonAmbiguousPairs[sample(length(nonAmbiguousPairs),
                                          nSnps, replace = TRUE)]

  data.frame(
    snp_id = paste0("rs", sample(1e6:9e6, nSnps)),
    effect_allele = vapply(allelePairs, `[`, character(1), 1),
    other_allele = vapply(allelePairs, `[`, character(1), 2),
    beta_ZX = betaZX,
    se_ZX = seZX,
    pval_ZX = 2 * pnorm(-abs(betaZX / seZX)),
    eaf = runif(nSnps, 0.05, 0.95),
    gene_region = paste0("GENE", seq_len(nSnps)),
    stringsAsFactors = FALSE
  )
}


#' Simulate Site Profile Likelihood Objects
#'
#' Generates a list of site profile objects as would be returned by
#' fitOutcomeModel() at multiple sites. Log-likelihood profiles are
#' quadratic (normal approximation) centered near the true beta_ZY value.
#'
#' @param nSites Number of sites to simulate.
#' @param betaGrid Numeric vector of grid points.
#' @param trueBeta True beta_ZY value. Profiles are centered near this.
#' @param nPerSite Approximate sample size per site (affects profile width).
#' @param seed Random seed for reproducibility.
#'
#' @return Named list of site profile objects, each with elements:
#'   siteId, betaGrid, logLikProfile, nCases, nControls, snpIds, diagnosticFlags.
simulateSiteProfiles <- function(nSites = 3,
                                 betaGrid = seq(-3, 3, by = 0.01),
                                 trueBeta = 0.5,
                                 nPerSite = 2000,
                                 seed = 42) {
  set.seed(seed)
  profiles <- list()

  for (s in seq_len(nSites)) {
    # Each site has a slightly different MLE due to sampling variability
    siteBeta <- trueBeta + rnorm(1, 0, 0.1)
    # Information proportional to sample size
    info <- nPerSite * runif(1, 0.8, 1.2) * 0.01
    # Quadratic log-likelihood profile
    logLik <- -0.5 * info * (betaGrid - siteBeta)^2
    # Normalize so max is 0
    logLik <- logLik - max(logLik)

    nCases <- round(nPerSite * runif(1, 0.05, 0.15))
    nControls <- nPerSite - nCases

    profiles[[paste0("site_", LETTERS[s])]] <- list(
      siteId = paste0("site_", LETTERS[s]),
      betaGrid = betaGrid,
      logLikProfile = logLik,
      nCases = nCases,
      nControls = nControls,
      snpIds = paste0("rs", 1:5),
      diagnosticFlags = list(
        weakInstruments = FALSE,
        lowCaseCount = nCases < 50,
        gridBoundaryMLE = FALSE
      )
    )
  }

  profiles
}


#' Simulate Covariate Data
#'
#' Generates a covariate matrix mimicking FeatureExtraction output.
#'
#' @param n Number of individuals.
#' @param nCovariates Number of covariates to generate.
#' @param seed Random seed for reproducibility.
#'
#' @return A data frame with person_id and nCovariates columns of binary
#'   or continuous covariates.
simulateCovariateData <- function(n = 1000,
                                  nCovariates = 50,
                                  seed = 42) {
  set.seed(seed)

  df <- data.frame(person_id = seq_len(n))

  for (j in seq_len(nCovariates)) {
    if (j <= nCovariates / 2) {
      # Binary covariates (conditions, drugs)
      prevalence <- runif(1, 0.01, 0.3)
      df[[paste0("covariate_", j)]] <- rbinom(n, 1, prevalence)
    } else {
      # Continuous covariates (measurements, age)
      df[[paste0("covariate_", j)]] <- rnorm(n)
    }
  }

  df
}
