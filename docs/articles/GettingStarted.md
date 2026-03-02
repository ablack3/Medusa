# Getting Started with Medusa

## What is Medusa?

Medusa (**M**endelian **E**stimation in **D**istributed **S**tandardized
**A**nalytics) is an R package that implements two-sample Mendelian
Randomization (MR) within the OHDSI ecosystem using the OMOP Common Data
Model.

**Mendelian Randomization** uses genetic variants as instrumental
variables to estimate causal effects of exposures on outcomes. Because
genetic variants are randomly allocated at conception (like a natural
randomized trial), MR can provide evidence for causality from
observational data.

**Two-sample MR** separates the analysis into two parts:

1.  **SNP-exposure associations** (beta_ZX): Obtained from published
    GWAS (typically from large biobank studies)
2.  **SNP-outcome associations** (beta_ZY): Estimated from your own
    cohort data

The causal estimate is the **Wald ratio**: beta_MR = beta_ZY / beta_ZX.

## Installation

``` r
# Install from GitHub
remotes::install_github("OHDSI/Medusa")

# This will also install OHDSI dependencies:
# DatabaseConnector, SqlRender, Cyclops, FeatureExtraction
```

## The Federated Approach

Medusa’s key innovation is **one-shot federated pooling**. Here’s how it
works:

1.  The **coordinator** assembles genetic instruments from public GWAS
    data
2.  Each **site** (hospital/biobank with OMOP CDM + genomic data):
    - Extracts the outcome cohort and genotype data locally
    - Fits a regression model locally
    - Computes a **profile log-likelihood curve** (a vector of numbers)
    - Shares ONLY this vector — no individual-level data leaves the site
3.  The **coordinator** sums the log-likelihood vectors across sites
4.  The causal estimate is extracted from the combined curve

This is privacy-preserving by design: the shared data is a smooth
numeric vector that cannot be reverse-engineered to identify individual
patients.

## Quick Start with Simulated Data

No database connection needed for this example.

``` r
library(Medusa)

# Simulate MR data with a known causal effect of 0.5 (log-OR)
simData <- simulateMRData(
  n = 5000,
  nSnps = 10,
  trueEffect = 0.5,
  seed = 42
)

# Look at the instrument table
head(simData$instrumentTable[, c("snp_id", "beta_ZX", "se_ZX", "eaf")])
#>   snp_id   beta_ZX      se_ZX       eaf
#> 1    rs1 0.1784231 0.04386823 0.2846324
#> 2    rs2 0.3862010 0.07320910 0.3445011
#> 3    rs3 0.2940426 0.05558554 0.2029028
#> 4    rs4 0.1500947 0.02021748 0.2164829
#> 5    rs5 0.4452723 0.05744303 0.1043854
#> 6    rs6 0.4793527 0.06068079 0.1473405
```

### Fit the Outcome Model

``` r
# Fit the outcome model and get the profile likelihood
profile <- fitOutcomeModel(
  cohortData = simData$data,
  covariateData = NULL,
  instrumentTable = simData$instrumentTable,
  betaGrid = seq(-2, 2, by = 0.05),
  siteId = "simulated_site"
)
#> Fitting outcome model at site 'simulated_site' (3452 cases, 1548 controls)...
#> Site 'simulated_site': beta_ZY_hat = 0.6137 (SE = 0.1410).

# The profile contains no individual-level data
names(profile)
#> [1] "siteId"          "betaGrid"        "logLikProfile"   "nCases"         
#> [5] "nControls"       "snpIds"          "diagnosticFlags" "betaHat"        
#> [9] "seHat"
cat(sprintf("beta_ZY estimate: %.3f (SE: %.3f)\n", profile$betaHat, profile$seHat))
#> beta_ZY estimate: 0.614 (SE: 0.141)
```

### Pool Profiles (Simulating Federation)

``` r
# Simulate 3 sites with similar data
profiles <- simulateSiteProfiles(
  nSites = 3,
  betaGrid = seq(-2, 2, by = 0.05),
  trueBeta = 0.5,
  nPerSite = 3000,
  seed = 42
)

# Pool the profiles
combined <- poolLikelihoodProfiles(profiles)
#> Pooling profile likelihoods from 3 site(s)...
#> Pooling complete: 3 sites, 955 total cases, 8045 total controls.
cat(sprintf("Pooled %d sites: %d total cases, %d total controls\n",
            combined$nSites, combined$totalCases, combined$totalControls))
#> Pooled 3 sites: 955 total cases, 8045 total controls
```

### Compute the MR Estimate

``` r
instruments <- simulateInstrumentTable(nSnps = 5)
estimate <- computeMREstimate(combined, instruments)
#> MR estimate: beta = 3.1823 (95% CI: 2.3144, 4.3395), p = 1.33e-05
#> Odds ratio: 24.102 (95% CI: 10.119, 76.667)

cat(sprintf("Causal estimate (beta_MR): %.3f\n", estimate$betaMR))
#> Causal estimate (beta_MR): 3.182
cat(sprintf("95%% CI: [%.3f, %.3f]\n", estimate$ciLower, estimate$ciUpper))
#> 95% CI: [2.314, 4.339]
cat(sprintf("P-value: %.2e\n", estimate$pValue))
#> P-value: 1.33e-05
cat(sprintf("Odds ratio: %.3f\n", estimate$oddsRatio))
#> Odds ratio: 24.102
```

### Run Sensitivity Analyses

``` r
# Create per-SNP estimates for sensitivity analyses
set.seed(42)
nSnps <- 10
betaZX <- rnorm(nSnps, 0.3, 0.05)
perSnp <- data.frame(
  snp_id = paste0("rs", 1:nSnps),
  beta_ZY = 0.5 * betaZX + rnorm(nSnps, 0, 0.02),
  se_ZY = rep(0.02, nSnps),
  beta_ZX = betaZX,
  se_ZX = rep(0.05, nSnps)
)

results <- runSensitivityAnalyses(perSnp)
#> Running sensitivity analyses with 10 SNPs...
#>   IVW...
#>   MR-Egger...
#>   Weighted Median...
#>   Steiger filtering...
#>     Steiger filter removed 8 of 10 SNPs.
#>   Leave-One-Out...
#> Sensitivity analyses complete.
print(results$summary)
#>                 method   beta_MR      se_MR   ci_lower  ci_upper          pval
#> 1                  IVW 0.4859306 0.01917957  0.4483387 0.5235226 1.287302e-141
#> 2             MR-Egger 0.2072935 0.25588407 -0.2942393 0.7088263  4.413080e-01
#> 3      Weighted Median 0.4848631 0.04409626  0.3984344 0.5712918  4.014147e-28
#> 4 Steiger-filtered IVW 0.3577932 0.04016741  0.2790651 0.4365213  5.217304e-19
```

## Understanding the Report

The [`generateMRReport()`](../reference/generateMRReport.md) function
creates a self-contained HTML report with:

1.  **Executive Summary** — Plain-language interpretation with strength
    rating
2.  **Instrument Summary** — Table of all SNPs with F-statistics
3.  **Likelihood Profile** — Visual of site and combined likelihood
    curves
4.  **Main Result** — Odds ratio with confidence interval
5.  **Sensitivity Analyses** — Forest plot comparing methods
6.  **Diagnostics** — PheWAS plot, negative controls, missingness
7.  **Methods Section** — Auto-generated text for manuscripts

## FAQ

**What genomic data linkage is required?** Your OMOP CDM site needs a
genomic linkage table mapping person_id to SNP genotypes (coded as 0/1/2
for allele count). This is increasingly common in biobank-linked health
systems.

**What if my site has no genomic data?** Medusa requires at least one
site with genomic data linked to OMOP CDM. Sites without genomic data
cannot participate in the federated analysis.

**How many sites are needed?** A single site is sufficient. More sites
increase statistical power and allow for cross-site validation. The
federated approach adds value with 2+ sites.

**What OMOP CDM version is required?** Medusa supports OMOP CDM v5.3 and
v5.4.
