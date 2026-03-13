<div id="main" class="col-md-9" role="main">

# Convert VCF-Style Genotype Strings to Integer Allele Dosage

<div class="ref-description section level2">

Converts genotype strings from the OMOP Genomic Extension's
VARIANT\_OCCURRENCE table to integer allele dosage values (0, 1, 2).
Handles VCF-style genotypes ("0/0", "0/1", "1/1", and phased equivalents
"0\|0", "0\|1", "1\|1") as well as plain integer strings ("0", "1",
"2"). Unrecognized values are converted to NA with a warning.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
convertGenotypeString(genotypeRaw)
```

</div>

</div>

<div class="section level2">

## Arguments

-   genotypeRaw:

    Character vector of raw genotype strings.

</div>

<div class="section level2">

## Value

Integer vector of allele dosage values (0, 1, or 2). Unrecognized values
are set to NA.

</div>

<div class="section level2">

## References

OHDSI Genomic CDM: <https://github.com/OHDSI/Genomic-CDM>

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
convertGenotypeString(c("0/0", "0/1", "1/1", "0|1", "2"))
#> [1] 0 1 2 1 2
```

</div>

</div>

</div>
