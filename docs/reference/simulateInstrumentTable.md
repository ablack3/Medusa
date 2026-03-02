# Simulate an Instrument Table

Creates a synthetic instrument table mimicking the output of
[`getMRInstruments`](getMRInstruments.md).

## Usage

``` r
simulateInstrumentTable(nSnps = 10, seed = 42)
```

## Arguments

- nSnps:

  Number of SNPs to simulate.

- seed:

  Random seed for reproducibility.

## Value

Data frame with columns: snp_id, effect_allele, other_allele, beta_ZX,
se_ZX, pval_ZX, eaf, gene_region.

## Examples

``` r
instruments <- simulateInstrumentTable(nSnps = 5)
head(instruments)
#>      snp_id effect_allele other_allele    beta_ZX      se_ZX      pval_ZX
#> 1 rs8218050             C            A  0.4112875 0.04746451 4.508801e-18
#> 2 rs3683480             G            A -0.1694095 0.06314674 7.301074e-03
#> 3 rs2737635             T            G  0.1089385 0.07608033 1.521759e-01
#> 4 rs7612958             A            G  0.1898588 0.03532573 7.678766e-08
#> 5 rs8132676             A            G  0.1212805 0.04773757 1.106729e-02
#>         eaf gene_region
#> 1 0.7138361       GENE1
#> 2 0.7799496       GENE2
#> 3 0.3992975       GENE3
#> 4 0.6666528       GENE4
#> 5 0.0535535       GENE5
```
