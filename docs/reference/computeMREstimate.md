<div id="main" class="col-md-9" role="main">

# Wald ratio MR estimate from pooled profile likelihood

<div class="ref-description section level2">

Takes the pooled profile log-likelihood from `poolLikelihoodProfiles`
and the instrument table to compute the MR causal estimate via the Wald
ratio. The SNP-outcome association (beta\_ZY) is estimated from the
profile likelihood peak, and the causal estimate is beta\_MR = beta\_ZY
/ beta\_ZX. For multi-SNP analyses, beta\_ZX is the association of the
exact weighted allele score used at each site, not a simple mean of the
component SNP effects. Uncertainty is propagated using the delta method.

Confidence intervals are computed both from the profile likelihood
(likelihood ratio method) and via the delta method (normal
approximation). The likelihood- based CI is preferred as it does not
assume normality.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
computeMREstimate(combinedProfile, instrumentTable, ciLevel = 0.95)
```

</div>

</div>

<div class="section level2">

## Arguments

-   combinedProfile:

    Output of `poolLikelihoodProfiles`.

-   instrumentTable:

    Data frame. Output of `getMRInstruments`. Used as a fallback for
    beta\_ZX and se\_ZX values when the pooled profile does not carry a
    stored score definition.

-   ciLevel:

    Numeric. Confidence interval level. Default is 0.95.

</div>

<div class="section level2">

## Value

A list with class "medusaMREstimate" containing:

-   betaZY:

    Point estimate of SNP-outcome association (profile MLE).

-   seZY:

    Standard error of betaZY from profile curvature.

-   betaMR:

    Causal effect estimate: betaZY / betaZX.

-   seMR:

    Standard error of betaMR via delta method.

-   ciLower:

    Lower confidence limit for betaMR.

-   ciUpper:

    Upper confidence limit for betaMR.

-   pValue:

    Two-sided p-value for betaMR.

-   oddsRatio:

    Exp(betaMR) — causal odds ratio for binary outcomes.

-   orCiLower:

    Lower CI for odds ratio.

-   orCiUpper:

    Upper CI for odds ratio.

-   ciLevel:

    The confidence level used.

-   betaZX:

    The exposure effect of the fitted allele score.

-   seZX:

    The SE of betaZX used.

-   nInstruments:

    Number of instruments used.

-   combinedProfile:

    The full combined profile for plotting.

</div>

<div class="section level2">

## Details

Compute Mendelian Randomization Causal Estimate

**Wald ratio**: For a single instrument, beta\_MR = beta\_ZY / beta\_ZX.
When multiple instruments are pooled into a single allele score, the
denominator must be the exposure effect of that same score:
$$\\beta\_{ZX, score} = \\sum\_j w\_j \\gamma\_j$$ where \\(w\_j\\) are
the fixed score weights used at the sites and \\(\\gamma\_j\\) are the
SNP-exposure coefficients.

**Delta method SE**: se\_MR = sqrt( (se\_ZY/beta\_ZX)^2 + (beta\_ZY \*
se\_ZX / beta\_ZX^2)^2 ).

**Likelihood-based CI**: The region of betaZY values where the
log-likelihood does not drop below peak - qchisq(ciLevel, df=1)/2. This
CI is then transformed to the MR scale via division by beta\_ZX.

</div>

<div class="section level2">

## References

Burgess, S., Butterworth, A., & Thompson, S. G. (2013). Mendelian
randomization analysis with multiple genetic variants using summarized
data. *Genetic Epidemiology*, 37(7), 658-665. doi:10.1002/gepi.21758.
Open access: https://pmc.ncbi.nlm.nih.gov/articles/PMC4377079/

</div>

<div class="section level2">

## See also

<div class="dont-index">

`poolLikelihoodProfiles`, `runSensitivityAnalyses`, `generateMRReport`

</div>

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
profiles <- simulateSiteProfiles(nSites = 3, trueBeta = 0.5)
combined <- poolLikelihoodProfiles(profiles)
#> Pooling profile likelihoods from 3 site(s)...
#> Pooling complete: 3 sites, 402 total cases, 5598 total controls.
instruments <- simulateInstrumentTable(nSnps = 5)
estimate <- computeMREstimate(combined, instruments)
#> MR estimate: beta = 1.3561 (95% CI: 0.7195, 1.9926), p = 6.42e-05
#> Odds ratio: 3.881 (95% CI: 2.054, 7.335)
print(estimate$betaMR)
#> [1] 1.356072
```

</div>

</div>

</div>
