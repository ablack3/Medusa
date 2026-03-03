<div id="main" class="col-md-9" role="main">

# Check for Strand-Ambiguous SNPs

<div class="ref-description section level2">

Identifies SNPs with ambiguous strand (A/T or G/C allele pairs).

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
isStrandAmbiguous(effectAllele, otherAllele)
```

</div>

</div>

<div class="section level2">

## Arguments

-   effectAllele:

    Character vector of effect alleles.

-   otherAllele:

    Character vector of other alleles.

</div>

<div class="section level2">

## Value

Logical vector indicating which SNPs are strand-ambiguous.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
isStrandAmbiguous(c("A", "G", "A"), c("T", "C", "C"))
#> [1]  TRUE  TRUE FALSE
```

</div>

</div>

</div>
