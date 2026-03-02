# Simulate MR Data with Known Causal Effect

Generates a complete dataset with SNP genotypes, an exposure,
confounders, and a binary outcome where the true causal effect of the
exposure on the outcome is known and recoverable. Useful for testing and
vignette examples without requiring a database connection.

## Usage

``` r
simulateMRData(
  n = 5000,
  nSnps = 10,
  trueEffect = 0.5,
  confoundingStrength = 0.3,
  snpEffectRange = c(0.1, 0.5),
  seed = 42
)
```

## Arguments

- n:

  Number of individuals.

- nSnps:

  Number of SNPs (instruments).

- trueEffect:

  True causal effect of exposure on outcome (log-OR scale).

- confoundingStrength:

  Strength of confounding (coefficient of confounder on both exposure
  and outcome).

- snpEffectRange:

  Range of SNP-exposure effect sizes.

- seed:

  Random seed for reproducibility.

## Value

A list with elements:

- data:

  Data frame with person_id, outcome (0/1), snp_1..snp_K, confounder_1,
  confounder_2, exposure.

- instrumentTable:

  Data frame mimicking getMRInstruments() output.

- trueEffect:

  The true causal effect used in simulation.

## Examples

``` r
simData <- simulateMRData(n = 1000, nSnps = 5, trueEffect = 0.3)
head(simData$data)
#>   person_id outcome snp_1 snp_2 snp_3 snp_4 snp_5 confounder_1 confounder_2
#> 1         1       1     1     1     0     0     1    1.3709584            1
#> 2         2       0     1     0     1     0     0   -0.5646982            0
#> 3         3       1     1     0     2     1     0    0.3631284            1
#> 4         4       1     0     0     1     1     0    0.6328626            1
#> 5         5       1     0     2     0     0     0    0.4042683            1
#> 6         6       1     0     1     0     0     0   -0.1061245            1
#>     exposure
#> 1  3.3068319
#> 2 -0.4567887
#> 3  1.8067824
#> 4  0.7664246
#> 5  3.0807866
#> 6 -0.1501653
```
