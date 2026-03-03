<div id="main" class="col-md-9" role="main">

# Build instrument table from user-provided data

<div class="ref-description section level2">

Creates a properly formatted instrument table from user-provided GWAS
summary statistics, bypassing the OpenGWAS API. Useful when working
offline or with custom GWAS results not available in OpenGWAS.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
createInstrumentTable(
  snpId,
  effectAllele,
  otherAllele,
  betaZX,
  seZX,
  pvalZX,
  eaf,
  geneRegion = NULL
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   snpId:

    Character vector of SNP rsIDs.

-   effectAllele:

    Character vector of effect alleles.

-   otherAllele:

    Character vector of other alleles.

-   betaZX:

    Numeric vector of SNP-exposure effect estimates.

-   seZX:

    Numeric vector of standard errors.

-   pvalZX:

    Numeric vector of p-values.

-   eaf:

    Numeric vector of effect allele frequencies.

-   geneRegion:

    Optional character vector of gene/region annotations.

</div>

<div class="section level2">

## Value

A data frame with the same structure as `getMRInstruments` output.

</div>

<div class="section level2">

## Details

Create Instrument Table from Local Data

</div>

<div class="section level2">

## See also

<div class="dont-index">

`getMRInstruments`

</div>

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
instruments <- createInstrumentTable(
  snpId = c("rs1234", "rs5678"),
  effectAllele = c("A", "G"),
  otherAllele = c("G", "C"),
  betaZX = c(0.5, 0.3),
  seZX = c(0.05, 0.08),
  pvalZX = c(1e-10, 1e-5),
  eaf = c(0.3, 0.45)
)
```

</div>

</div>

</div>
