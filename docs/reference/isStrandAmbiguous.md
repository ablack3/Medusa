# Check for Strand-Ambiguous SNPs

Identifies SNPs with ambiguous strand (A/T or G/C allele pairs).

## Usage

``` r
isStrandAmbiguous(effectAllele, otherAllele)
```

## Arguments

- effectAllele:

  Character vector of effect alleles.

- otherAllele:

  Character vector of other alleles.

## Value

Logical vector indicating which SNPs are strand-ambiguous.

## Examples

``` r
isStrandAmbiguous(c("A", "G", "A"), c("T", "C", "C"))
#> [1]  TRUE  TRUE FALSE
```
