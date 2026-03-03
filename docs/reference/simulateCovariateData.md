<div id="main" class="col-md-9" role="main">

# Simulate Covariate Data

<div class="ref-description section level2">

Generates a covariate matrix mimicking FeatureExtraction output.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
simulateCovariateData(n = 1000, nCovariates = 50, seed = 42)
```

</div>

</div>

<div class="section level2">

## Arguments

-   n:

    Number of individuals.

-   nCovariates:

    Number of covariates to generate.

-   seed:

    Random seed for reproducibility.

</div>

<div class="section level2">

## Value

A data frame with person\_id and nCovariates columns of binary or
continuous covariates.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
covData <- simulateCovariateData(n = 100, nCovariates = 10)
head(covData)
#>   person_id covariate_1 covariate_2 covariate_3 covariate_4 covariate_5
#> 1         1           1           0           0           0           0
#> 2         2           0           0           0           1           1
#> 3         3           1           1           0           0           0
#> 4         4           0           1           1           0           0
#> 5         5           0           0           0           0           0
#> 6         6           1           0           0           0           0
#>   covariate_6 covariate_7 covariate_8 covariate_9 covariate_10
#> 1   1.2449329   0.5428314   1.5977432   0.3625802    1.2499199
#> 2  -1.0512478  -0.5304505  -0.6919057  -0.2747531    1.1463134
#> 3  -0.6650712   0.1811529  -1.1865409   1.1425879   -0.0569313
#> 4   0.7340706   1.1606279  -0.4467962  -0.3645218    0.5824822
#> 5  -0.5089338   0.6364922   1.3443038   0.6579949   -0.4493756
#> 6  -0.4272471   0.5952489   0.7001235  -1.2361387   -0.4379546
```

</div>

</div>

</div>
