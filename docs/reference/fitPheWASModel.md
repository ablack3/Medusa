<div id="main" class="col-md-9" role="main">

# Fit a single PheWAS model for one covariate against all SNPs using Cyclops

<div class="ref-description section level2">

Fit a single PheWAS model for one covariate against all SNPs using
Cyclops

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
fitPheWASModel(modelDf, snpCols, covariateName, pValueThreshold)
```

</div>

</div>

<div class="section level2">

## Arguments

-   modelDf:

    Data frame with column `y` (covariate values) and one column per
    SNP.

-   snpCols:

    Character vector of SNP column names in `modelDf`.

-   covariateName:

    Character. Name of the covariate being tested.

-   pValueThreshold:

    Numeric. Bonferroni-corrected significance threshold.

</div>

<div class="section level2">

## Value

A data frame with one row per SNP (columns: snp\_id, covariate\_name,
beta, se, pval, significant), or NULL if the model cannot be fit.

</div>

</div>
