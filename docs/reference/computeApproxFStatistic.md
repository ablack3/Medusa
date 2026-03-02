# Compute Approximate F-statistic from GWAS Summary Statistics

Computes the approximation F = (beta_ZX / se_ZX)^2 for each SNP. This is
the standard approximation used when individual-level data is not
available.

## Usage

``` r
computeApproxFStatistic(betaZX, seZX)
```

## Arguments

- betaZX:

  Numeric vector of SNP-exposure effect estimates.

- seZX:

  Numeric vector of standard errors for SNP-exposure effects.

## Value

Numeric vector of approximate F-statistics.

## Examples

``` r
computeApproxFStatistic(c(0.5, 0.3), c(0.05, 0.1))
#> [1] 100   9
```
