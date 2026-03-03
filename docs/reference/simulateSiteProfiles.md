<div id="main" class="col-md-9" role="main">

# Simulate Site Profile Likelihood Objects

<div class="ref-description section level2">

Generates a list of site profile objects as would be returned by
`fitOutcomeModel` at multiple sites. Log-likelihood profiles are
quadratic (normal approximation) centered near the true beta\_ZY value.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
simulateSiteProfiles(
  nSites = 3,
  betaGrid = seq(-3, 3, by = 0.01),
  trueBeta = 0.5,
  nPerSite = 2000,
  seed = 42
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   nSites:

    Number of sites to simulate.

-   betaGrid:

    Numeric vector of grid points.

-   trueBeta:

    True beta\_ZY value. Profiles are centered near this.

-   nPerSite:

    Approximate sample size per site (affects profile width).

-   seed:

    Random seed for reproducibility.

</div>

<div class="section level2">

## Value

Named list of site profile objects, each with elements: siteId,
betaGrid, logLikProfile, nCases, nControls, snpIds, diagnosticFlags.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
profiles <- simulateSiteProfiles(nSites = 3, trueBeta = 0.5)
names(profiles)
#> [1] "site_A" "site_B" "site_C"
```

</div>

</div>

</div>
