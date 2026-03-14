<div id="main" class="col-md-9" role="main">

# Simulate Negative Control Outcome Columns

<div class="ref-description section level2">

Adds negative control outcome columns to existing cohort data. These
outcomes are generated independently of genotype (true null effects)
with varying prevalence.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
simulateNegativeControlOutcomes(cohortData, nControls = 20, seed = 99)
```

</div>

</div>

<div class="section level2">

## Arguments

-   cohortData:

    Data frame with at least a `person_id` column.

-   nControls:

    Integer. Number of negative control outcomes to simulate.

-   seed:

    Random seed for reproducibility.

</div>

<div class="section level2">

## Value

The input `cohortData` with appended `nc_outcome_<i>` columns (binary
0/1).

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
simData <- simulateMRData(n = 1000, nSnps = 5)
withNC <- simulateNegativeControlOutcomes(simData$data, nControls = 10)
head(withNC)
#>   person_id outcome snp_rs1 snp_rs2 snp_rs3 snp_rs4 snp_rs5 confounder_1
#> 1         1       1       1       1       0       0       1    1.3709584
#> 2         2       0       1       0       1       0       0   -0.5646982
#> 3         3       1       1       0       2       1       0    0.3631284
#> 4         4       1       0       0       1       1       0    0.6328626
#> 5         5       1       0       2       0       0       0    0.4042683
#> 6         6       1       0       1       0       0       0   -0.1061245
#>   confounder_2   exposure nc_outcome_1 nc_outcome_2 nc_outcome_3 nc_outcome_4
#> 1            1  3.3068319            0            1            0            0
#> 2            0 -0.4567887            0            0            0            0
#> 3            1  1.8067824            1            1            1            0
#> 4            1  0.7664246            0            1            0            0
#> 5            1  3.0807866            1            0            0            0
#> 6            1 -0.1501653            0            0            0            0
#>   nc_outcome_5 nc_outcome_6 nc_outcome_7 nc_outcome_8 nc_outcome_9
#> 1            1            0            0            0            1
#> 2            0            0            0            0            1
#> 3            0            0            0            0            0
#> 4            0            0            0            0            0
#> 5            0            0            0            0            0
#> 6            0            0            0            0            0
#>   nc_outcome_10
#> 1             0
#> 2             1
#> 3             0
#> 4             0
#> 5             0
#> 6             0
```

</div>

</div>

</div>
