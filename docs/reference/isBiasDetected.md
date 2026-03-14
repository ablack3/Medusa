<div id="main" class="col-md-9" role="main">

# Detect systematic bias from negative control estimates

<div class="ref-description section level2">

Tests whether the mean of negative control MR estimates is significantly
different from zero using a one-sample t-test.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
isBiasDetected(ncEstimates)
```

</div>

</div>

<div class="section level2">

## Arguments

-   ncEstimates:

    Data frame of negative control estimates with columns `beta_MR` and
    `se_MR`.

</div>

<div class="section level2">

## Value

Logical. TRUE if bias is detected (p &lt; 0.05).

</div>

</div>
