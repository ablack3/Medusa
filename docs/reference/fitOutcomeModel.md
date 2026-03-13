<div id="main" class="col-md-9" role="main">

# Outcome model with exact grid likelihood evaluation

<div class="ref-description section level2">

Fits a logistic regression of the binary outcome on SNP genotype(s) plus
covariates, then evaluates the profile log-likelihood across a
pre-specified grid of beta\_ZY values. This is the core methodological
function for federated MR: each site runs this locally and shares only
the resulting log-likelihood profile vector.

Two analysis modes are supported: `alleleScore` fits a single model
using a weighted allele score as the genetic exposure variable, while
`perSNP` fits separate models for each SNP (needed for multi-SNP
sensitivity analyses).

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
fitOutcomeModel(
  cohortData,
  covariateData = NULL,
  instrumentTable,
  betaGrid = seq(-3, 3, by = 0.01),
  regularizationVariance = 0.1,
  instrumentRegularization = FALSE,
  outcomeType = "binary",
  analysisType = "alleleScore",
  siteId = "site_1",
  modelBackend = "glm"
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   cohortData:

    Data frame. Output of `buildMRCohort`.

-   covariateData:

    Covariate data. Output of `buildMRCovariates` or a data frame with
    person\_id and covariate columns.

-   instrumentTable:

    Data frame. Output of `getMRInstruments`.

-   betaGrid:

    Numeric vector. Grid of beta\_ZY values at which to evaluate the
    profile log-likelihood. Default is `seq(-3, 3, by = 0.01)` (601 grid
    points).

-   regularizationVariance:

    Numeric. When `modelBackend = "cyclops"`, this is the variance of
    the normal prior applied to nuisance coefficients. Smaller values
    imply stronger shrinkage; use `Inf` for an unpenalized Cyclops fit.
    Ignored by the `"glm"` backend. Default is 0.1.

-   instrumentRegularization:

    Logical. When `modelBackend = "cyclops"`, whether the exposure
    coefficient is included in the Cyclops prior. Default is FALSE so
    the profiled exposure coefficient remains unpenalized. Ignored by
    the `"glm"` backend.

-   outcomeType:

    Character. Type of outcome: "binary" for logistic regression.
    Default is "binary".

-   analysisType:

    Character. "alleleScore" for single weighted score model or "perSNP"
    for separate per-SNP models. Default is "alleleScore".

-   siteId:

    Character. Identifier for this site. Included in the returned
    profile object. Default is "site\_1".

-   modelBackend:

    Character. Outcome-model fitting backend: `"glm"` uses base R
    logistic regression, while `"cyclops"` uses Cyclops for scalable
    logistic regression with optional Gaussian shrinkage on nuisance
    covariates. Default is `"glm"`.

</div>

<div class="section level2">

## Value

A list with class "medusaSiteProfile" containing:

-   siteId:

    Character site identifier.

-   betaGrid:

    Numeric vector of grid points (same as input).

-   logLikProfile:

    Numeric vector of profile log-likelihood values, same length as
    betaGrid. In `perSNP` mode this remains the valid allele-score
    profile used for pooling.

-   nCases:

    Number of outcome cases.

-   nControls:

    Number of controls.

-   snpIds:

    Character vector of SNP IDs used.

-   diagnosticFlags:

    List of diagnostic flags from the model fit.

-   betaHat:

    Point estimate of beta\_ZY (MLE from profile).

-   seHat:

    Approximate standard error from profile curvature.

-   perSnpEstimates:

    When `analysisType = "perSNP"`, a data frame of per-SNP beta\_ZY /
    se\_ZY estimates for summarized-data sensitivity analyses.

This object contains no individual-level data and is safe to share.

</div>

<div class="section level2">

## Details

Fit Outcome Model and Evaluate Profile Log-Likelihood

The profile log-likelihood evaluation works as follows:

1.  Fit the unconstrained model to obtain the MLE and its Wald SE.

2.  At each grid point, fix beta\_ZY to that value and re-fit the
    nuisance parameters using the fixed term as an offset.

3.  Record the maximized constrained log-likelihood at each grid point.
    This is the exact profile likelihood for the coefficient in a
    logistic generalized linear model.

When `modelBackend = "cyclops"` and `regularizationVariance` is finite,
the nuisance parameters are estimated under a Gaussian prior at both the
unconstrained optimum and every grid point. In that configuration the
returned objective is a penalized profile, not an unpenalized MLE
profile.

When using the allele score, weights are \\(w\_j = \\gamma\_j /
\\mathrm{SE}(\\gamma\_j)^2\\), normalized by \\(\\sum\_j \|w\_j\|\\).
The same weights are reused downstream so that the MR denominator
matches the exact score fitted at each site.

Missing SNP dosages in the allele-score model are imputed to the
expected dosage \\(2 \\times \\mathrm{EAF}\_j\\) from the instrument
table, rather than being treated as homozygous reference. This avoids a
systematic downward bias in the score when genotype missingness is
present.

</div>

<div class="section level2">

## References

Suchard, M. A., et al. (2013). Massive parallelization of serial
inference algorithms for a complex generalized linear model. *ACM
Transactions on Modeling and Computer Simulation*, 23(1), 1-17.

</div>

<div class="section level2">

## See also

<div class="dont-index">

`poolLikelihoodProfiles`, `computeMREstimate`, `buildMRCohort`

</div>

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
simData <- simulateMRData(n = 2000, nSnps = 5, trueEffect = 0.3)
profile <- fitOutcomeModel(
  cohortData = simData$data,
  covariateData = NULL,
  instrumentTable = simData$instrumentTable,
  betaGrid = seq(-2, 2, by = 0.05)
)
#> Fitting outcome model at site 'site_1' (1171 cases, 829 controls)...
#> Site 'site_1': beta_ZY_hat = 0.3273 (SE = 0.1919).
plot(profile$betaGrid, profile$logLikProfile, type = "l",
     xlab = "beta_ZY", ylab = "Profile log-likelihood")

```

</div>

</div>

</div>
