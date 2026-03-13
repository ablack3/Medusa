<div id="main" class="col-md-9" role="main">

# Allele harmonization for Mendelian Randomization

<div class="ref-description section level2">

Ensures that the effect allele coding in the instrument table (from
GWAS) matches the allele coding in the genotype data at each site.
Handles allele flipping (when alleles are swapped) and detects
strand-ambiguous SNPs (A/T and G/C pairs) that cannot be reliably
harmonized.

When alleles need to be flipped, beta\_ZX is multiplied by -1 and EAF is
replaced by 1 - EAF. Palindromic SNPs with EAF near 0.5 (within
`eafPalindromicThreshold` of 0.5) are dropped because strand cannot be
reliably inferred from allele frequency.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
harmonizeAlleles(
  instrumentTable,
  genotypeAlleles,
  eafPalindromicThreshold = 0.08,
  cohortAlleleFrequencies = NULL
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   instrumentTable:

    Data frame with columns: snp\_id, effect\_allele, other\_allele,
    beta\_ZX, se\_ZX, eaf. Typically the output of `getMRInstruments`.

-   genotypeAlleles:

    Data frame with columns: snp\_id, allele\_coded, allele\_noncoded.
    Describes the allele coding used in the genotype data at a
    particular site.

-   eafPalindromicThreshold:

    Numeric threshold for removing palindromic SNPs. Palindromic SNPs
    with EAF between (0.5 - threshold) and (0.5 + threshold) are
    removed. Default is 0.08 (i.e., EAF between 0.42 and 0.58).

-   cohortAlleleFrequencies:

    Optional named numeric vector of cohort allele frequencies, keyed by
    SNP ID. The frequency should be for the coded allele in the genotype
    data (i.e., mean(genotype)/2). Used to resolve palindromic SNPs when
    EAF is far enough from 0.5.

</div>

<div class="section level2">

## Value

A list with two elements:

-   instrumentTable:

    The harmonized instrument table with updated effect\_allele,
    other\_allele, beta\_ZX, and eaf where allele flips were applied. A
    logical column `flipped` is added.

-   removedSnps:

    Data frame of SNPs that were removed, with a column `reason`
    indicating why (e.g., "palindromic\_ambiguous\_eaf",
    "allele\_mismatch").

</div>

<div class="section level2">

## Details

Harmonize Alleles Between Instrument Table and Genotype Data

The harmonization algorithm proceeds as follows for each SNP:

1.  If the instrument effect\_allele matches the genotype allele\_coded:
    no action needed.

2.  If the instrument effect\_allele matches the genotype
    allele\_noncoded: flip the coding (beta\_ZX \*= -1, eaf = 1 - eaf).

3.  If neither direct match is found, try complement alleles
    (A&lt;-&gt;T, G&lt;-&gt;C):

    -   If the SNP is palindromic (A/T or G/C) and EAF is near 0.5,
        remove the SNP.

    -   If palindromic but EAF clearly differs from 0.5, use EAF to
        infer strand orientation.

4.  If no match is found even with complements: remove the SNP with
    reason "allele\_mismatch".

</div>

<div class="section level2">

## References

Hartwig, F. P., Davies, N. M., Hemani, G., & Davey Smith, G. (2016).
Two-sample Mendelian randomization: avoiding the downsides of a
powerful, widely applicable but potentially fallible technique.
*International Journal of Epidemiology*, 45(6), 1717-1726.

</div>

<div class="section level2">

## See also

<div class="dont-index">

`isStrandAmbiguous`, `getMRInstruments`, `buildMRCohort`

</div>

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
instruments <- data.frame(
  snp_id = c("rs1", "rs2", "rs3"),
  effect_allele = c("A", "G", "A"),
  other_allele = c("G", "T", "C"),
  beta_ZX = c(0.5, -0.3, 0.2),
  se_ZX = c(0.05, 0.08, 0.06),
  pval_ZX = c(1e-10, 1e-5, 1e-8),
  eaf = c(0.3, 0.45, 0.6),
  stringsAsFactors = FALSE
)
genotypeAlleles <- data.frame(
  snp_id = c("rs1", "rs2", "rs3"),
  allele_coded = c("G", "G", "A"),
  allele_noncoded = c("A", "T", "C"),
  stringsAsFactors = FALSE
)
result <- harmonizeAlleles(instruments, genotypeAlleles)
```

</div>

</div>

</div>
