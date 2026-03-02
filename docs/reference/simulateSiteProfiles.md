# Simulate Site Profile Likelihood Objects

Generates a list of site profile objects as would be returned by
[`fitOutcomeModel`](fitOutcomeModel.md) at multiple sites.
Log-likelihood profiles are quadratic (normal approximation) centered
near the true beta_ZY value.

## Usage

``` r
simulateSiteProfiles(
  nSites = 3,
  betaGrid = seq(-3, 3, by = 0.01),
  trueBeta = 0.5,
  nPerSite = 2000,
  seed = 42
)
```

## Arguments

- nSites:

  Number of sites to simulate.

- betaGrid:

  Numeric vector of grid points.

- trueBeta:

  True beta_ZY value. Profiles are centered near this.

- nPerSite:

  Approximate sample size per site (affects profile width).

- seed:

  Random seed for reproducibility.

## Value

Named list of site profile objects, each with elements: siteId,
betaGrid, logLikProfile, nCases, nControls, snpIds, diagnosticFlags.

## Examples

``` r
profiles <- simulateSiteProfiles(nSites = 3, trueBeta = 0.5)
names(profiles)
#> [1] "site_A" "site_B" "site_C"
```
