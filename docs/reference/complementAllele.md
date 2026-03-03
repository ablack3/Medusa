<div id="main" class="col-md-9" role="main">

# Get Complement of a DNA Allele

<div class="ref-description section level2">

Returns the Watson-Crick complement: A&lt;-&gt;T, G&lt;-&gt;C.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
complementAllele(allele)
```

</div>

</div>

<div class="section level2">

## Arguments

-   allele:

    Character string representing an allele (A, T, G, or C).

</div>

<div class="section level2">

## Value

Character string with the complement allele.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
complementAllele("A")  # returns "T"
#> [1] "T"
complementAllele("G")  # returns "C"
#> [1] "C"
```

</div>

</div>

</div>
