# Estimate Standard Error from Profile Log-Likelihood Curvature

Uses the second derivative of the profile log-likelihood at the MLE to
estimate the standard error.

## Usage

``` r
estimateSEFromProfile(betaGrid, logLikProfile)
```

## Arguments

- betaGrid:

  Numeric vector of grid points.

- logLikProfile:

  Numeric vector of log-likelihood values.

## Value

Numeric standard error estimate.
