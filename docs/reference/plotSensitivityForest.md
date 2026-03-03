<div id="main" class="col-md-9" role="main">

# Forest plot comparing MR methods

<div class="ref-description section level2">

Creates a horizontal forest plot showing estimates from all sensitivity
analysis methods side by side with confidence intervals.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
plotSensitivityForest(sensitivityResults)
```

</div>

</div>

<div class="section level2">

## Arguments

-   sensitivityResults:

    Output of `runSensitivityAnalyses`.

</div>

<div class="section level2">

## Value

A ggplot2 object.

</div>

<div class="section level2">

## Details

Plot Sensitivity Analysis Forest Plot

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

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
#> Error in validatePerSnpSummaryData(perSnpEstimates): Assertion on 'requiredCols' failed: Must be a subset of {'snp_id','beta_ZY','se_ZY','beta_ZX','se_ZX'}, but has additional elements {'effect_allele','other_allele','eaf'}.
plotSensitivityForest(results)
#> Error: object 'results' not found
```

</div>

</div>

</div>
