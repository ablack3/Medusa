<div id="main" class="col-md-9" role="main">

# Simulate an Instrument Table

<div class="ref-description section level2">

Creates a synthetic instrument table mimicking the output of
`getMRInstruments`.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
simulateInstrumentTable(nSnps = 10, seed = 42)
```

</div>

</div>

<div class="section level2">

## Arguments

-   nSnps:

    Number of SNPs to simulate.

-   seed:

    Random seed for reproducibility.

</div>

<div class="section level2">

## Value

Data frame with columns: snp\_id, effect\_allele, other\_allele,
beta\_ZX, se\_ZX, pval\_ZX, eaf, gene\_region.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
instruments <- simulateInstrumentTable(nSnps = 5)
head(instruments)
#>      snp_id effect_allele other_allele    beta_ZX      se_ZX      pval_ZX
#> 1 rs8218050             C            A  0.4659224 0.04746451 9.586362e-23
#> 2 rs3683480             G            A  0.4748302 0.06314674 5.499750e-14
#> 3 rs2737635             T            G  0.2144558 0.07608033 4.820243e-03
#> 4 rs7612958             A            G -0.4321791 0.03532573 2.043187e-34
#> 5 rs8132676             A            G  0.3566982 0.04773757 7.894616e-14
#>         eaf gene_region
#> 1 0.7138361       GENE1
#> 2 0.7799496       GENE2
#> 3 0.3992975       GENE3
#> 4 0.6666528       GENE4
#> 5 0.0535535       GENE5
```

</div>

</div>

</div>
