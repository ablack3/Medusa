# Forest plot comparing MR methods

Creates a horizontal forest plot showing estimates from all sensitivity
analysis methods side by side with confidence intervals.

## Usage

``` r
plotSensitivityForest(sensitivityResults)
```

## Arguments

- sensitivityResults:

  Output of [`runSensitivityAnalyses`](runSensitivityAnalyses.md).

## Value

A ggplot2 object.

## Details

Plot Sensitivity Analysis Forest Plot

## Examples

``` r
set.seed(42)
nSnps <- 10
perSnp <- data.frame(
  snp_id = paste0("rs", 1:nSnps),
  beta_ZY = rnorm(nSnps, 0.15, 0.02),
  se_ZY = rep(0.02, nSnps),
  beta_ZX = rnorm(nSnps, 0.3, 0.05),
  se_ZX = rep(0.05, nSnps)
)
results <- runSensitivityAnalyses(perSnp)
#> Running sensitivity analyses with 10 SNPs...
#>   IVW...
#>   MR-Egger...
#>   Weighted Median...
#>   Steiger filtering...
#>     Steiger filter removed 9 of 10 SNPs.
#>   Leave-One-Out...
#> Sensitivity analyses complete.
plotSensitivityForest(results)
#> `height` was translated to `width`.

```
