# Get Complement of a DNA Allele

Returns the Watson-Crick complement: A\<-\>T, G\<-\>C.

## Usage

``` r
complementAllele(allele)
```

## Arguments

- allele:

  Character string representing an allele (A, T, G, or C).

## Value

Character string with the complement allele.

## Examples

``` r
complementAllele("A")  # returns "T"
#> [1] "T"
complementAllele("G")  # returns "C"
#> [1] "C"
```
