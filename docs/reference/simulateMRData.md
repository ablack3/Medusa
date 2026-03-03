<div id="main" class="col-md-9" role="main">

# Simulate MR Data with Known Causal Effect

<div class="ref-description section level2">

Generates a complete dataset with SNP genotypes, an exposure,
confounders, and a binary outcome where the true causal effect of the
exposure on the outcome is known and recoverable. Useful for testing and
vignette examples without requiring a database connection.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

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

</div>

</div>

<div class="section level2">

## Arguments

-   n:

    Number of individuals.

-   nSnps:

    Number of SNPs (instruments).

-   trueEffect:

    True causal effect of exposure on outcome (log-OR scale).

-   confoundingStrength:

    Strength of confounding (coefficient of confounder on both exposure
    and outcome).

-   snpEffectRange:

    Range of SNP-exposure effect sizes.

-   seed:

    Random seed for reproducibility.

</div>

<div class="section level2">

## Value

A list with elements:

-   data:

    Data frame with person\_id, outcome (0/1), `snp_<sanitized rsID>`
    columns, confounder\_1, confounder\_2, exposure.

-   instrumentTable:

    Data frame mimicking getMRInstruments() output.

-   trueEffect:

    The true causal effect used in simulation.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
simData <- simulateMRData(n = 1000, nSnps = 5, trueEffect = 0.3)
head(simData$data)
#>   person_id outcome snp_rs1 snp_rs2 snp_rs3 snp_rs4 snp_rs5 confounder_1
#> 1         1       1       1       1       0       0       1    1.3709584
#> 2         2       0       1       0       1       0       0   -0.5646982
#> 3         3       1       1       0       2       1       0    0.3631284
#> 4         4       1       0       0       1       1       0    0.6328626
#> 5         5       1       0       2       0       0       0    0.4042683
#> 6         6       1       0       1       0       0       0   -0.1061245
#>   confounder_2   exposure
#> 1            1  3.3068319
#> 2            0 -0.4567887
#> 3            1  1.8067824
#> 4            1  0.7664246
#> 5            1  3.0807866
#> 6            1 -0.1501653
```

</div>

</div>

</div>
