<div id="main" class="col-md-9" role="main">

# Likelihood profile visualization

<div class="ref-description section level2">

Creates a ggplot2 visualization of the profile log-likelihood curve(s),
showing individual site profiles (if provided) and the combined profile
with MLE and confidence interval.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
plotLikelihoodProfile(
  combinedProfile,
  siteProfileList = NULL,
  mrEstimate = NULL,
  title = "Profile Log-Likelihood"
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   combinedProfile:

    Output of `poolLikelihoodProfiles`.

-   siteProfileList:

    Optional named list of site profile objects.

-   mrEstimate:

    Optional output of `computeMREstimate` for annotating the MLE
    and CI.

-   title:

    Plot title. Default is "Profile Log-Likelihood".

</div>

<div class="section level2">

## Value

A ggplot2 object.

</div>

<div class="section level2">

## Details

Plot Profile Likelihood

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
profiles <- simulateSiteProfiles(nSites = 3, trueBeta = 0.5)
combined <- poolLikelihoodProfiles(profiles)
#> Pooling profile likelihoods from 3 site(s)...
#> Pooling complete: 3 sites, 637 total cases, 5363 total controls.
plotLikelihoodProfile(combined, profiles)

```

</div>

</div>

</div>
