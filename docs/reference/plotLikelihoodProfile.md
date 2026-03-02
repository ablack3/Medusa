# Likelihood profile visualization

Creates a ggplot2 visualization of the profile log-likelihood curve(s),
showing individual site profiles (if provided) and the combined profile
with MLE and confidence interval.

## Usage

``` r
plotLikelihoodProfile(
  combinedProfile,
  siteProfileList = NULL,
  mrEstimate = NULL,
  title = "Profile Log-Likelihood"
)
```

## Arguments

- combinedProfile:

  Output of [`poolLikelihoodProfiles`](poolLikelihoodProfiles.md).

- siteProfileList:

  Optional named list of site profile objects.

- mrEstimate:

  Optional output of [`computeMREstimate`](computeMREstimate.md) for
  annotating the MLE and CI.

- title:

  Plot title. Default is "Profile Log-Likelihood".

## Value

A ggplot2 object.

## Details

Plot Profile Likelihood

## Examples

``` r
profiles <- simulateSiteProfiles(nSites = 3, trueBeta = 0.5)
combined <- poolLikelihoodProfiles(profiles)
#> Pooling profile likelihoods from 3 site(s)...
#> Pooling complete: 3 sites, 637 total cases, 5363 total controls.
plotLikelihoodProfile(combined, profiles)

```
