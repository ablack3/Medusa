<div id="main" class="col-md-9" role="main">

# Compute Approximate F-statistic from GWAS Summary Statistics

<div class="ref-description section level2">

Computes the approximation F = (beta\_ZX / se\_ZX)^2 for each SNP. This
is the standard approximation used when individual-level data is not
available.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
computeApproxFStatistic(betaZX, seZX)
```

</div>

</div>

<div class="section level2">

## Arguments

-   betaZX:

    Numeric vector of SNP-exposure effect estimates.

-   seZX:

    Numeric vector of standard errors for SNP-exposure effects.

</div>

<div class="section level2">

## Value

Numeric vector of approximate F-statistics.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
computeApproxFStatistic(c(0.5, 0.3), c(0.05, 0.1))
#> [1] 100   9
```

</div>

</div>

</div>
