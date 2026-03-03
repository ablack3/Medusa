<div id="main" class="col-md-9" role="main">

# Estimate Standard Error from Profile Log-Likelihood Curvature

<div class="ref-description section level2">

Uses the second derivative of the profile log-likelihood at the MLE to
estimate the standard error.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
estimateSEFromProfile(betaGrid, logLikProfile)
```

</div>

</div>

<div class="section level2">

## Arguments

-   betaGrid:

    Numeric vector of grid points.

-   logLikProfile:

    Numeric vector of log-likelihood values.

</div>

<div class="section level2">

## Value

Numeric standard error estimate.

</div>

</div>
