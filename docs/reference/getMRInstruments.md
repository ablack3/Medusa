# Instrument assembly from OpenGWAS

Queries the IEU OpenGWAS database via the `ieugwasr` package to retrieve
GWAS summary statistics for a specified exposure trait, applies LD
clumping to obtain independent instruments, computes approximate
F-statistics, flags strand-ambiguous SNPs, and returns a clean
instrument table ready for distribution to federated analysis sites.

This function runs at the coordinator node. The returned instrument
table is serialized to disk and distributed to all sites unchanged.

## Usage

``` r
getMRInstruments(
  exposureTraitId,
  pThreshold = 5e-08,
  r2Threshold = 0.001,
  kb = 10000,
  ancestryPopulation = "EUR",
  additionalSnps = NULL
)
```

## Arguments

- exposureTraitId:

  Character string. IEU OpenGWAS trait ID (e.g., "ieu-a-1119" for IL-6
  receptor levels).

- pThreshold:

  Numeric. Genome-wide significance p-value threshold for selecting
  SNPs. Default is 5e-8.

- r2Threshold:

  Numeric. LD clumping r-squared threshold. SNP pairs with r-squared
  above this value will be pruned, keeping the more significant SNP.
  Default is 0.001.

- kb:

  Numeric. LD clumping window in kilobases. Default is 10000.

- ancestryPopulation:

  Character. Reference panel ancestry population for LD clumping. One of
  "EUR", "EAS", "AFR", "SAS", "AMR". Default is "EUR".

- additionalSnps:

  Optional character vector of SNP rsIDs to force-include in the
  instrument set (added after clumping, not subject to LD pruning).
  Default is NULL.

## Value

A data frame with one row per independent instrument SNP and columns:

- snp_id:

  rsID of the SNP (e.g., "rs2228145").

- effect_allele:

  The allele associated with increased exposure level.

- other_allele:

  The non-effect allele.

- beta_ZX:

  Effect size of the SNP on the exposure (log scale).

- se_ZX:

  Standard error of beta_ZX.

- pval_ZX:

  P-value for the SNP-exposure association.

- eaf:

  Effect allele frequency in the GWAS reference population.

- gene_region:

  Nearest gene or genomic region annotation.

- fStatistic:

  Approximate F-statistic: (beta_ZX / se_ZX)^2.

- strandAmbiguous:

  Logical. TRUE if the SNP is strand-ambiguous (A/T or G/C allele pair).

The data frame also carries the following attributes:

- retrievalTimestamp:

  POSIXct timestamp of when instruments were retrieved.

- exposureTraitId:

  The trait ID used for retrieval.

- parameters:

  List of all parameter values used.

## Details

Retrieve and Clump Genetic Instruments for Mendelian Randomization

The function queries the IEU OpenGWAS API for all SNP associations with
the specified trait below `pThreshold`, then applies LD clumping using
the specified reference panel to retain only independent instruments.
The approximate F-statistic is computed as (beta / SE)^2 for each SNP.
SNPs with F \< 10 are flagged as potentially weak instruments via a
warning message but are not removed automatically.

If the OpenGWAS API is unavailable, the function throws an informative
error suggesting the user provide a cached instrument table instead.

## References

Hemani, G., et al. (2018). The MR-Base platform supports systematic
causal inference across the human phenome. *eLife*, 7, e34408.

## See also

[`harmonizeAlleles`](harmonizeAlleles.md),
[`buildMRCohort`](buildMRCohort.md),
[`computeApproxFStatistic`](computeApproxFStatistic.md)

## Examples

``` r
# Using simulated data (no API call)
instruments <- simulateInstrumentTable(nSnps = 10)
head(instruments)
#>      snp_id effect_allele other_allele     beta_ZX      se_ZX      pval_ZX
#> 1 rs3781070             T            G  0.41128753 0.07424188 3.027626e-08
#> 2 rs2440057             A            G -0.16940945 0.02832261 2.211995e-09
#> 3 rs1521507             T            G  0.10893852 0.07933350 1.696990e-01
#> 4 rs5432853             G            T  0.18985878 0.07680009 1.343157e-02
#> 5 rs2154668             T            C  0.12128050 0.02494625 1.163985e-06
#> 6 rs6986500             A            G -0.03183735 0.05085271 5.312690e-01
#>         eaf gene_region
#> 1 0.6580465       GENE1
#> 2 0.9345355       GENE2
#> 3 0.7335898       GENE3
#> 4 0.5598396       GENE4
#> 5 0.8147207       GENE5
#> 6 0.2205265       GENE6

if (FALSE) { # \dontrun{
# Real API call (requires internet)
instruments <- getMRInstruments(
  exposureTraitId = "ieu-a-1119",
  pThreshold = 5e-8,
  r2Threshold = 0.001,
  kb = 10000,
  ancestryPopulation = "EUR"
)
} # }
```
