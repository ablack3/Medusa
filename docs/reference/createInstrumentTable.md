# Build instrument table from user-provided data

Creates a properly formatted instrument table from user-provided GWAS
summary statistics, bypassing the OpenGWAS API. Useful when working
offline or with custom GWAS results not available in OpenGWAS.

## Usage

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

## Arguments

- snpId:

  Character vector of SNP rsIDs.

- effectAllele:

  Character vector of effect alleles.

- otherAllele:

  Character vector of other alleles.

- betaZX:

  Numeric vector of SNP-exposure effect estimates.

- seZX:

  Numeric vector of standard errors.

- pvalZX:

  Numeric vector of p-values.

- eaf:

  Numeric vector of effect allele frequencies.

- geneRegion:

  Optional character vector of gene/region annotations.

## Value

A data frame with the same structure as
[`getMRInstruments`](getMRInstruments.md) output.

## Details

Create Instrument Table from Local Data

## See also

[`getMRInstruments`](getMRInstruments.md)

## Examples

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
