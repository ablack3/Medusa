# Reshape Long-Format Genotype Data to Wide Format

Converts genotype data from long format (person_id, snp_id, genotype) to
wide format with one column per SNP. Missing genotypes are coded as NA.

## Usage

``` r
reshapeGenotypes(genotypeData, instrumentTable)
```

## Arguments

- genotypeData:

  Data frame in long format with columns personId, snpId, genotype.

- instrumentTable:

  Instrument table for SNP ordering.

## Value

Data frame in wide format with personId and one column per SNP.
