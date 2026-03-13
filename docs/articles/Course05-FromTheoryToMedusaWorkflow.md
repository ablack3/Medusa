<div id="main" class="col-md-9" role="main">

# Course 5: From Theory to the Medusa Workflow

<div class="section level2">

## What you’ll learn

By the end of this chapter, you should be able to:

-   map textbook two-sample MR quantities onto Medusa objects,
-   explain the roles of the coordinator and local sites,
-   explain why Medusa shares profile likelihoods instead of
    person-level data,
-   walk through the main Medusa analysis functions on simulated data.

</div>

<div class="section level2">

## Why this matters

The earlier chapters explained two-sample MR as a summarized-data
method. Medusa implements that same causal logic in a privacy-preserving
federated workflow.

The exposure side still comes from external summary statistics. The main
difference is how the outcome side is estimated and shared.

</div>

<div class="section level2">

## The Medusa division of labor

<div class="section level3">

### Coordinator tasks

The coordinator typically:

-   defines the scientific question,
-   assembles or imports the instrument table,
-   distributes the instrument table to participating sites,
-   receives site-level summary artifacts,
-   pools them and computes the final MR estimate.

</div>

<div class="section level3">

### Site tasks

Each site typically:

-   extracts the outcome cohort,
-   harmonizes genotype coding,
-   fits the outcome model locally,
-   exports a profile likelihood summary,
-   optionally exports per-SNP summary estimates for sensitivity
    analyses.

**Key takeaway:** Medusa separates the exposure side, the local
outcome-model fitting, and the final pooled MR estimation.

</div>

</div>

<div class="section level2">

## How Medusa represents the MR ingredients

In textbook notation:

-   `beta_ZX`: SNP-exposure associations from external GWAS,
-   `beta_ZY`: SNP-outcome association estimated from the outcome
    sample.

In Medusa:

-   the instrument table stores the SNP-exposure side (`beta_ZX`,
    `se_ZX`, alleles, effect allele frequency),
-   `fitOutcomeModel()` estimates the outcome-side signal locally,
-   `poolLikelihoodProfiles()` combines site summaries,
-   `computeMREstimate()` converts the pooled profile into the final MR
    estimate.

</div>

<div class="section level2">

## Why profile likelihood is shared

Medusa’s primary federated estimate does not ask each site to export all
person-level outcomes or raw genotypes. Instead, each site evaluates a
profile log-likelihood over a grid of candidate SNP-outcome effect
values.

The coordinator receives:

-   the grid,
-   the log-likelihood values,
-   metadata such as case counts,
-   the allele-score definition.

Because log-likelihoods from independent sites add, the coordinator can
pool the profiles exactly by summation.

This is what makes the approach both federated and statistically
principled.

</div>

<div class="section level2">

## A runnable Medusa example

<div id="cb1" class="sourceCode">

``` r
library(Medusa)

beta_grid <- seq(-2, 2, by = 0.05)

sim_data <- simulateMRData(
  n = 4000,
  nSnps = 8,
  trueEffect = 0.40,
  seed = 1005
)

head(sim_data$instrumentTable[, c("snp_id", "beta_ZX", "se_ZX", "eaf")])
#>   snp_id   beta_ZX      se_ZX       eaf
#> 1    rs1 0.4936651 0.06309374 0.1736507
#> 2    rs2 0.1734556 0.02556424 0.2347729
#> 3    rs3 0.1970980 0.04546400 0.1961966
#> 4    rs4 0.3080739 0.05459666 0.1726639
#> 5    rs5 0.4419107 0.05693616 0.3474797
#> 6    rs6 0.4473812 0.07394154 0.2378871
```

</div>

<div class="section level3">

### Step 1: fit the local outcome model

<div id="cb2" class="sourceCode">

``` r
profile_site_a <- fitOutcomeModel(
  cohortData = sim_data$data,
  covariateData = NULL,
  instrumentTable = sim_data$instrumentTable,
  betaGrid = beta_grid,
  siteId = "site_A"
)
#> Fitting outcome model at site 'site_A' (2621 cases, 1379 controls)...
#> Site 'site_A': beta_ZY_hat = 0.4640 (SE = 0.1438).

names(profile_site_a)
#>  [1] "siteId"          "betaGrid"        "logLikProfile"   "nCases"         
#>  [5] "nControls"       "snpIds"          "diagnosticFlags" "betaHat"        
#>  [9] "seHat"           "scoreDefinition"
```

</div>

This local profile is Medusa’s site-level representation of the
SNP-outcome evidence.

</div>

<div class="section level3">

### Step 2: create partner site summaries

To simulate federation, we create two compatible partner profiles that
reuse the same allele-score definition.

<div id="cb3" class="sourceCode">

``` r
make_partner_profile <- function(template, site_id, shift, info_scale, n_cases, n_controls) {
  partner <- template
  partner$siteId <- site_id
  partner$betaHat <- template$betaHat + shift
  partner$seHat <- template$seHat / sqrt(info_scale)
  partner$logLikProfile <- -0.5 * info_scale *
    ((template$betaGrid - partner$betaHat) / template$seHat)^2
  partner$logLikProfile <- partner$logLikProfile - max(partner$logLikProfile)
  partner$nCases <- n_cases
  partner$nControls <- n_controls
  partner
}

site_profiles <- list(
  site_A = profile_site_a,
  site_B = make_partner_profile(profile_site_a, "site_B", shift = 0.04,
                                info_scale = 1.10, n_cases = 230, n_controls = 1770),
  site_C = make_partner_profile(profile_site_a, "site_C", shift = -0.03,
                                info_scale = 0.95, n_cases = 210, n_controls = 1790)
)
```

</div>

</div>

<div class="section level3">

### Step 3: pool the profiles

<div id="cb4" class="sourceCode">

``` r
combined_profile <- poolLikelihoodProfiles(site_profiles)
#> Pooling profile likelihoods from 3 site(s)...
#> Pooling complete: 3 sites, 3061 total cases, 4939 total controls.

combined_profile$nSites
#> [1] 3
combined_profile$totalCases
#> [1] 3061
combined_profile$totalControls
#> [1] 4939
```

</div>

</div>

<div class="section level3">

### Step 4: recover the MR estimate

<div id="cb5" class="sourceCode">

``` r
mr_estimate <- computeMREstimate(combined_profile, sim_data$instrumentTable)
#> MR estimate: beta = 1.4779 (95% CI: 1.1495, 1.9705), p = 1.68e-07
#> Odds ratio: 4.384 (95% CI: 3.157, 7.175)

mr_estimate$betaMR
#> [1] 1.477905
mr_estimate$ciLower
#> [1] 1.149481
mr_estimate$ciUpper
#> [1] 1.97054
mr_estimate$oddsRatio
#> [1] 4.383751
```

</div>

</div>

</div>

<div class="section level2">

## Mapping objects back to theory

-   `sim_data$instrumentTable` corresponds to the SNP-exposure side.
-   `profile_site_a$betaHat` is the local site’s fitted score-to-outcome
    effect.
-   `combined_profile$logLikProfile` is the pooled evidence about the
    outcome-side parameter.
-   `mr_estimate$betaMR` is the causal estimate after dividing the
    pooled SNP-outcome signal by the exposure effect of the same fitted
    score.

This last point is important: Medusa’s main federated estimator is built
around the exact allele score used locally, not a generic average SNP
effect.

**Key takeaway:** Medusa preserves the logic of two-sample MR while
changing the way outcome information is communicated across sites.

</div>

<div class="section level2">

## Optional per-SNP mode for sensitivity analyses

<div id="cb6" class="sourceCode">

``` r
profile_per_snp <- fitOutcomeModel(
  cohortData = sim_data$data,
  covariateData = NULL,
  instrumentTable = sim_data$instrumentTable,
  betaGrid = beta_grid,
  siteId = "site_A",
  analysisType = "perSNP"
)
#> Fitting outcome model at site 'site_A' (2621 cases, 1379 controls)...
#> Site 'site_A': beta_ZY_hat = 0.4640 (SE = 0.1438).

head(profile_per_snp$perSnpEstimates)
#>   snp_id effect_allele other_allele       eaf     beta_ZY      se_ZY   beta_ZX
#> 1    rs1             C            A 0.1736507  0.15532721 0.06590682 0.4936651
#> 2    rs2             C            T 0.2347729 -0.03576511 0.05656689 0.1734556
#> 3    rs3             T            C 0.1961966  0.05255082 0.06120696 0.1970980
#> 4    rs4             T            G 0.1726639  0.10805974 0.06489736 0.3080739
#> 5    rs5             C            A 0.3474797  0.22455152 0.05161634 0.4419107
#> 6    rs6             C            T 0.2378871  0.19933685 0.05721299 0.4473812
#>        se_ZX      pval_ZX      pval_ZY
#> 1 0.06309374 5.104392e-15 0.0184346747
#> 2 0.02556424 1.160163e-11 0.5272155342
#> 3 0.04546400 1.455910e-05 0.3905745821
#> 4 0.05459666 1.673798e-08 0.0958954331
#> 5 0.05693616 8.392247e-15 0.0000135892
#> 6 0.07394154 1.444222e-09 0.0004937608
```

</div>

Those `perSnpEstimates` are the package’s bridge back to standard
summarized-data MR sensitivity analyses.

</div>

<div class="section level2">

## Exercise set 1

<div class="section level3">

### Exercise 1

In Medusa, where does the SNP-exposure information live?

Show solution

It lives in the instrument table, which carries quantities such as
`beta_ZX`, `se_ZX`, allele labels, and effect allele frequencies.

</div>

<div class="section level3">

### Exercise 2

Why does Medusa pool log-likelihood profiles instead of person-level
records?

Show solution

Because the site-level log-likelihood profile is a privacy-preserving
summary that can be added across independent sites to recover the pooled
evidence without sharing person-level genotypes or outcomes.

</div>

<div class="section level3">

### Exercise 3

What does `analysisType = "perSNP"` add beyond the main federated
estimate?

Show solution

It provides per-SNP outcome associations that can be used for
summarized-data sensitivity analyses such as IVW, MR-Egger, weighted
median, and leave-one-out.

</div>

</div>

<div class="section level2">

## Exercise set 2: inspect Medusa objects

Use the objects created above and answer:

1.  Which object stores the number of pooled sites?
2.  Which object stores the final odds ratio?

Show solution

`combined_profile$nSites` stores the number of pooled sites, and
`mr_estimate$oddsRatio` stores the final odds ratio.

</div>

<div class="section level2">

## Chapter summary

-   Medusa keeps the exposure side in the instrument table and estimates
    the outcome side locally.
-   Sites share profile-likelihood summaries rather than person-level
    records.
-   The coordinator pools those profiles and computes the MR estimate
    from the same fitted allele score used at the sites.
-   Per-SNP output is available when you want standard summarized-data
    sensitivity analyses.

</div>

<div class="section level2">

## Next chapter

Next: [Course 6: Doing and critiquing a two-sample MR
study](Course06-DoingAndCritiquingATwoSampleMRStudy.md)

For more package detail after this chapter, see [Getting Started with
Medusa](GettingStarted.md) and [Federated Analysis
Guide](FederatedAnalysisGuide.md).

</div>

</div>
