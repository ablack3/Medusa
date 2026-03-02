# Wald ratio MR estimate from pooled profile likelihood

Takes the pooled profile log-likelihood from
[`poolLikelihoodProfiles`](poolLikelihoodProfiles.md) and the instrument
table to compute the MR causal estimate via the Wald ratio. The
SNP-outcome association (beta_ZY) is estimated from the profile
likelihood peak, and the causal estimate is beta_MR = beta_ZY / beta_ZX.
Uncertainty is propagated using the delta method.

Confidence intervals are computed both from the profile likelihood
(likelihood ratio method) and via the delta method (normal
approximation). The likelihood- based CI is preferred as it does not
assume normality.

## Usage

``` r
computeMREstimate(combinedProfile, instrumentTable, ciLevel = 0.95)
```

## Arguments

- combinedProfile:

  Output of [`poolLikelihoodProfiles`](poolLikelihoodProfiles.md).

- instrumentTable:

  Data frame. Output of [`getMRInstruments`](getMRInstruments.md). Used
  for beta_ZX and se_ZX values.

- ciLevel:

  Numeric. Confidence interval level. Default is 0.95.

## Value

A list with class "medusaMREstimate" containing:

- betaZY:

  Point estimate of SNP-outcome association (profile MLE).

- seZY:

  Standard error of betaZY from profile curvature.

- betaMR:

  Causal effect estimate: betaZY / betaZX.

- seMR:

  Standard error of betaMR via delta method.

- ciLower:

  Lower confidence limit for betaMR.

- ciUpper:

  Upper confidence limit for betaMR.

- pValue:

  Two-sided p-value for betaMR.

- oddsRatio:

  Exp(betaMR) — causal odds ratio for binary outcomes.

- orCiLower:

  Lower CI for odds ratio.

- orCiUpper:

  Upper CI for odds ratio.

- ciLevel:

  The confidence level used.

- betaZX:

  The SNP-exposure effect used (mean across instruments).

- seZX:

  The SE of betaZX used.

- nInstruments:

  Number of instruments used.

- combinedProfile:

  The full combined profile for plotting.

## Details

Compute Mendelian Randomization Causal Estimate

**Wald ratio**: For a single instrument, beta_MR = beta_ZY / beta_ZX.
When multiple instruments are pooled into a single allele score, the
effective beta_ZX is the precision-weighted mean.

**Delta method SE**: se_MR = sqrt( (se_ZY/beta_ZX)^2 + (beta_ZY \* se_ZX
/ beta_ZX^2)^2 ).

**Likelihood-based CI**: The region of betaZY values where the
log-likelihood does not drop below peak - qchisq(ciLevel, df=1)/2. This
CI is then transformed to the MR scale via division by beta_ZX.

## References

Burgess, S., Butterworth, A., & Thompson, S. G. (2013). Mendelian
randomization analysis with multiple genetic variants using summarized
data. *Genetic Epidemiology*, 37(7), 658-665.

## See also

[`poolLikelihoodProfiles`](poolLikelihoodProfiles.md),
[`runSensitivityAnalyses`](runSensitivityAnalyses.md),
[`generateMRReport`](generateMRReport.md)

## Examples

``` r
profiles <- simulateSiteProfiles(nSites = 3, trueBeta = 0.5)
combined <- poolLikelihoodProfiles(profiles)
#> Pooling profile likelihoods from 3 site(s)...
#> Pooling complete: 3 sites, 637 total cases, 5363 total controls.
instruments <- simulateInstrumentTable(nSnps = 5)
estimate <- computeMREstimate(combined, instruments)
#> MR estimate: beta = 3.2980 (95% CI: 1.8515, 4.7445), p = 1.16e-04
#> Odds ratio: 27.058 (95% CI: 6.369, 114.948)
print(estimate$betaMR)
#> [1] 3.297994
```
