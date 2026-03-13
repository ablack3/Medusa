<div id="main" class="col-md-9" role="main">

# Compute Cohort Allele Frequencies from Genotype Data

<div class="ref-description section level2">

Computes the coded (alternate) allele frequency for each SNP from the
extracted genotype data. Used for palindromic SNP resolution during
allele harmonization.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
computeCohortAlleleFrequencies(genotypeData)
```

</div>

</div>

<div class="section level2">

## Arguments

-   genotypeData:

    Data frame with columns snpId and genotype (integer).

</div>

<div class="section level2">

## Value

Named numeric vector of allele frequencies, keyed by SNP ID.

</div>

</div>
