# Standard Mendelian Randomization sensitivity analyses

Implements multiple MR estimation methods that are robust to different
patterns of instrument invalidity. Requires per-SNP beta_ZY estimates
(available when `analysisType = "perSNP"` in
[`fitOutcomeModel`](fitOutcomeModel.md)).

Methods include: Inverse Variance Weighted (IVW), MR-Egger regression,
Weighted Median, Steiger directional filtering, and Leave-One-Out
analysis. Concordance across methods strengthens causal evidence.

## Usage

``` r
runSensitivityAnalyses(
  perSnpEstimates,
  methods = c("IVW", "MREgger", "WeightedMedian", "Steiger", "LeaveOneOut"),
  outcomeSampleSize = NULL,
  exposureSampleSize = NULL
)
```

## Arguments

- perSnpEstimates:

  Data frame with per-SNP summary statistics. Required columns: snp_id,
  beta_ZY, se_ZY, beta_ZX, se_ZX.

- methods:

  Character vector of methods to run. Default is
  `c("IVW", "MREgger", "WeightedMedian", "Steiger", "LeaveOneOut")`.

- outcomeSampleSize:

  Optional integer. Total sample size for outcome GWAS/analysis. Needed
  for Steiger filtering.

- exposureSampleSize:

  Optional integer. Total sample size for exposure GWAS. Needed for
  Steiger filtering.

## Value

A named list with class "medusaSensitivity" containing:

- ivw:

  Data frame with method, beta_MR, se_MR, ci_lower, ci_upper, pval.

- mrEgger:

  Data frame with beta_MR, se_MR, pval, plus intercept, intercept_se,
  intercept_pval.

- weightedMedian:

  Data frame with beta_MR, se_MR, ci_lower, ci_upper, pval.

- steiger:

  Data frame with IVW results after Steiger filtering, plus n_removed
  (number of SNPs removed).

- leaveOneOut:

  Data frame with snp_removed, beta_MR, se_MR, pval for each SNP
  dropped.

- summary:

  Data frame comparing all methods side by side.

## Details

Run MR Sensitivity Analyses

**IVW**: Weighted regression of beta_ZY on beta_ZX through the origin,
with weights 1/se_ZY^2. This is the primary MR estimate assuming all
instruments are valid.

**MR-Egger**: Weighted regression of beta_ZY on beta_ZX with an
intercept. A non-zero intercept indicates directional pleiotropy. The
slope is a consistent causal estimate even under directional pleiotropy,
provided the InSIDE assumption holds.

**Weighted Median**: Bootstrap-based weighted median of per-SNP Wald
ratio estimates. Produces a consistent estimate if fewer than 50 weight
comes from invalid instruments.

**Steiger**: Tests whether each SNP explains more variance in the
exposure than the outcome (the expected causal direction). SNPs failing
this test are removed and IVW is re-run.

**Leave-One-Out**: Drops each SNP in turn and recomputes the IVW
estimate. Identifies influential outlier instruments.

## References

Bowden, J., Davey Smith, G., & Burgess, S. (2015). Mendelian
randomization with invalid instruments: effect estimation and bias
detection through Egger regression. *International Journal of
Epidemiology*, 44(2), 512-525.

Bowden, J., et al. (2016). Consistent estimation in Mendelian
randomization with some invalid instruments using a weighted median
estimator. *Genetic Epidemiology*, 40(4), 304-314.

Hemani, G., Tilling, K., & Davey Smith, G. (2017). Orienting the causal
relationship between imprecisely measured traits using GWAS summary
data. *PLoS Genetics*, 13(11), e1007081.

## See also

[`computeMREstimate`](computeMREstimate.md),
[`fitOutcomeModel`](fitOutcomeModel.md),
[`generateMRReport`](generateMRReport.md)

## Examples

``` r
# Simulate per-SNP estimates
set.seed(42)
nSnps <- 10
betaZX <- rnorm(nSnps, 0.3, 0.05)
betaZY <- 0.5 * betaZX + rnorm(nSnps, 0, 0.02)
perSnp <- data.frame(
  snp_id = paste0("rs", 1:nSnps),
  beta_ZY = betaZY,
  se_ZY = rep(0.02, nSnps),
  beta_ZX = betaZX,
  se_ZX = rep(0.05, nSnps)
)
results <- runSensitivityAnalyses(perSnp)
#> Running sensitivity analyses with 10 SNPs...
#>   IVW...
#>   MR-Egger...
#>   Weighted Median...
#>   Steiger filtering...
#>     Steiger filter removed 8 of 10 SNPs.
#>   Leave-One-Out...
#> Sensitivity analyses complete.
results$summary
#>                 method   beta_MR      se_MR   ci_lower  ci_upper          pval
#> 1                  IVW 0.4859306 0.01917957  0.4483387 0.5235226 1.287302e-141
#> 2             MR-Egger 0.2072935 0.25588407 -0.2942393 0.7088263  4.413080e-01
#> 3      Weighted Median 0.4848631 0.04409626  0.3984344 0.5712918  4.014147e-28
#> 4 Steiger-filtered IVW 0.3577932 0.04016741  0.2790651 0.4365213  5.217304e-19
```
