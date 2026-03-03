<div id="main" class="col-md-9" role="main">

# Getting Started with Medusa

<div class="section level2">

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

1.  **SNP-exposure associations** (beta\_ZX): Obtained from published
    GWAS (typically from large biobank studies)
2.  **SNP-outcome associations** (beta\_ZY): Estimated from your own
    cohort data

The causal estimate is the **Wald ratio**: beta\_MR = beta\_ZY /
beta\_ZX.

This vignette is a runnable sandbox: every example below uses synthetic
data, so the goal is to show the workflow and the shape of the outputs
rather than to make a scientific claim.

</div>

<div class="section level2">

## Installation

<div id="cb1" class="sourceCode">

``` r
# Install from GitHub
remotes::install_github("OHDSI/Medusa")

# This will also install OHDSI dependencies:
# DatabaseConnector, SqlRender, Cyclops, FeatureExtraction
```

</div>

</div>

<div class="section level2">

## The Federated Approach

Medusa’s key innovation is **one-shot federated pooling**. Here’s how it
works:

1.  The **coordinator** assembles genetic instruments from public GWAS
    data
2.  Each **site** (hospital/biobank with OMOP CDM + genomic data):
    -   Extracts the outcome cohort and genotype data locally
    -   Fits a regression model locally
    -   Computes a **profile log-likelihood curve** (a vector of
        numbers)
    -   Shares ONLY this vector — no individual-level data leaves the
        site
3.  The **coordinator** sums the log-likelihood vectors across sites
4.  The causal estimate is extracted from the combined curve

This is privacy-preserving by design: the shared data is a smooth
numeric vector that cannot be reverse-engineered to identify individual
patients.

</div>

<div class="section level2">

## Quick Start with Simulated Data

No database connection needed for this example.

<div id="cb2" class="sourceCode">

``` r
library(Medusa)

betaGrid <- seq(-2, 2, by = 0.05)

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

</div>

<div class="section level3">

### Fit the Outcome Model

<div id="cb3" class="sourceCode">

``` r
# Fit the outcome model and get the profile likelihood
profileSiteA <- fitOutcomeModel(
  cohortData = simData$data,
  covariateData = NULL,
  instrumentTable = simData$instrumentTable,
  betaGrid = betaGrid,
  siteId = "site_A"
)
#> Fitting outcome model at site 'site_A' (3452 cases, 1548 controls)...
#> Site 'site_A': beta_ZY_hat = 0.6137 (SE = 0.1410).

# The profile contains no individual-level data
names(profileSiteA)
#>  [1] "siteId"          "betaGrid"        "logLikProfile"   "nCases"         
#>  [5] "nControls"       "snpIds"          "diagnosticFlags" "betaHat"        
#>  [9] "seHat"           "scoreDefinition"
cat(sprintf("beta_ZY estimate: %.3f (SE: %.3f)\n",
            profileSiteA$betaHat, profileSiteA$seHat))
#> beta_ZY estimate: 0.614 (SE: 0.141)
```

</div>

</div>

<div class="section level3">

### Pool Profiles (Simulating Federation)

<div id="cb4" class="sourceCode">

``` r
# Create synthetic partner sites that reuse the same allele-score definition.
# This keeps the pooled likelihood compatible with the denominator used later
# in computeMREstimate().
makePartnerProfile <- function(template, siteId, shift, infoScale, nCases, nControls) {
  partner <- template
  partner$siteId <- siteId
  partner$betaHat <- template$betaHat + shift
  partner$seHat <- template$seHat / sqrt(infoScale)
  partner$logLikProfile <- -0.5 * infoScale *
    ((template$betaGrid - partner$betaHat) / template$seHat)^2
  partner$logLikProfile <- partner$logLikProfile - max(partner$logLikProfile)
  partner$nCases <- nCases
  partner$nControls <- nControls
  partner
}

profiles <- list(
  site_A = profileSiteA,
  site_B = makePartnerProfile(profileSiteA, "site_B", shift = 0.06,
                              infoScale = 1.1, nCases = 295, nControls = 2705),
  site_C = makePartnerProfile(profileSiteA, "site_C", shift = -0.04,
                              infoScale = 0.9, nCases = 260, nControls = 2740)
)

# Pool the profiles
combined <- poolLikelihoodProfiles(profiles)
#> Pooling profile likelihoods from 3 site(s)...
#> Pooling complete: 3 sites, 4007 total cases, 6993 total controls.
cat(sprintf("Pooled %d sites: %d total cases, %d total controls\n",
            combined$nSites, combined$totalCases, combined$totalControls))
#> Pooled 3 sites: 4007 total cases, 6993 total controls
```

</div>

</div>

<div class="section level3">

### Compute the MR Estimate

<div id="cb5" class="sourceCode">

``` r
estimate <- computeMREstimate(combined, simData$instrumentTable)
#> MR estimate: beta = 2.1158 (95% CI: 1.7631, 2.6447), p = 3.99e-12
#> Odds ratio: 8.296 (95% CI: 5.831, 14.079)

cat(sprintf("Causal estimate (beta_MR): %.3f\n", estimate$betaMR))
#> Causal estimate (beta_MR): 2.116
cat(sprintf("95%% CI: [%.3f, %.3f]\n", estimate$ciLower, estimate$ciUpper))
#> 95% CI: [1.763, 2.645]
cat(sprintf("P-value: %.2e\n", estimate$pValue))
#> P-value: 3.99e-12
cat(sprintf("Odds ratio: %.3f\n", estimate$oddsRatio))
#> Odds ratio: 8.296
```

</div>

The key detail is that the pooled profile carries the exact allele-score
definition used at each site. Passing the same instrument table back
into `computeMREstimate()` keeps the Wald ratio denominator aligned with
the fitted score rather than with an unrelated simulated instrument set.

</div>

<div class="section level3">

### Run Sensitivity Analyses

<div id="cb6" class="sourceCode">

``` r
# Create per-SNP estimates for sensitivity analyses
set.seed(42)
nSnps <- 10
instrumentSummary <- simulateInstrumentTable(nSnps = nSnps, seed = 42)
perSnp <- data.frame(
  snp_id = instrumentSummary$snp_id,
  beta_ZY = 0.5 * instrumentSummary$beta_ZX + rnorm(nSnps, 0, 0.02),
  se_ZY = rep(0.02, nSnps),
  beta_ZX = instrumentSummary$beta_ZX,
  se_ZX = instrumentSummary$se_ZX,
  effect_allele = instrumentSummary$effect_allele,
  other_allele = instrumentSummary$other_allele,
  eaf = instrumentSummary$eaf
)

results <- runSensitivityAnalyses(
  perSnp,
  methods = c("IVW", "MREgger", "WeightedMedian", "LeaveOneOut"),
  engine = "internal"
)
#> Running sensitivity analyses with 10 SNPs...
#>   Engine: internal
#>   IVW...
#>   MR-Egger...
#>   Weighted Median...
#>   Leave-One-Out...
#> Sensitivity analyses complete.
print(results$summary)
#>            method   beta_MR      se_MR  ci_lower  ci_upper         pval
#> 1             IVW 0.4586540 0.02188781 0.4157539 0.5015541 1.697680e-97
#> 2        MR-Egger 0.4634544 0.03576486 0.3933553 0.5335535 1.191296e-06
#> 3 Weighted Median 0.4764280 0.02650498 0.4244783 0.5283778 3.056401e-72
```

</div>

</div>

</div>

<div class="section level2">

## Understanding the Report

The `generateMRReport()` function creates a self-contained HTML report
with:

1.  **Executive Summary** — Plain-language interpretation with strength
    rating
2.  **Instrument Summary** — Table of all SNPs with F-statistics
3.  **Likelihood Profile** — Visual of site and combined likelihood
    curves
4.  **Main Result** — Odds ratio with confidence interval
5.  **Sensitivity Analyses** — Forest plot comparing methods
6.  **Diagnostics** — instrument PheWAS, allele-frequency checks,
    missingness
7.  **Methods Section** — Auto-generated text for manuscripts

</div>

<div class="section level2">

## FAQ

**What genomic data linkage is required?** Your OMOP CDM site needs a
genomic linkage table mapping person\_id to SNP genotypes (coded as
0/1/2 for allele count). This is increasingly common in biobank-linked
health systems.

**What if my site has no genomic data?** Medusa requires at least one
site with genomic data linked to OMOP CDM. Sites without genomic data
cannot participate in the federated analysis.

**How many sites are needed?** A single site is sufficient. More sites
increase statistical power and allow for cross-site validation. The
federated approach adds value with 2+ sites.

**What OMOP CDM version is required?** Medusa supports OMOP CDM v5.3 and
v5.4.

</div>

</div>
