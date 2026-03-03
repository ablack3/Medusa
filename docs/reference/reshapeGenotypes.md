<div id="main" class="col-md-9" role="main">

# Reshape Long-Format Genotype Data to Wide Format

<div class="ref-description section level2">

Converts genotype data from long format (person\_id, snp\_id, genotype)
to wide format with one column per SNP. Missing genotypes are coded as
NA.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
reshapeGenotypes(genotypeData, instrumentTable)
```

</div>

</div>

<div class="section level2">

## Arguments

-   genotypeData:

    Data frame in long format with columns personId, snpId, genotype.

-   instrumentTable:

    Instrument table for SNP ordering.

</div>

<div class="section level2">

## Value

Data frame in wide format with personId and one column per SNP.

</div>

</div>
