<div id="main" class="col-md-9" role="main">

# Standard Mendelian Randomization sensitivity analyses

<div class="ref-description section level2">

Implements multiple MR estimation methods that are robust to different
patterns of instrument invalidity. Requires per-SNP beta\_ZY estimates
(available when `analysisType = "perSNP"` in `fitOutcomeModel`).

Methods include: Inverse Variance Weighted (IVW), MR-Egger regression,
Weighted Median, Steiger directional filtering, and Leave-One-Out
analysis. Concordance across methods strengthens causal evidence.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
runSensitivityAnalyses(
  perSnpEstimates,
  methods = c("IVW", "MREgger", "WeightedMedian", "Steiger", "LeaveOneOut"),
  outcomeSampleSize = NULL,
  exposureSampleSize = NULL,
  outcomeType = "binary",
  engine = "auto",
  bootstrapSeed = 42
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   perSnpEstimates:

    Data frame with per-SNP summary statistics. Required columns:
    `snp_id`, `beta_ZY`, `se_ZY`, `beta_ZX`, `se_ZX`, `effect_allele`,
    `other_allele`, and `eaf`. These allele fields make the object
    harmonisation-ready for delegation to TwoSampleMR. Optional columns
    `pval_ZX` and `pval_ZY` are used when present and otherwise derived
    from the reported estimates.

-   methods:

    Character vector of methods to run. Default is
    `c("IVW", "MREgger", "WeightedMedian", "Steiger", "LeaveOneOut")`.

-   outcomeSampleSize:

    Optional integer. Total sample size for outcome GWAS/analysis.
    Needed for Steiger filtering.

-   exposureSampleSize:

    Optional integer. Total sample size for exposure GWAS. Needed for
    Steiger filtering.

-   outcomeType:

    Character. Type of outcome summary statistics represented in
    `perSnpEstimates`. Medusa's internal Steiger implementation is
    available only for `"continuous"` outcomes, while
    `engine = "TwoSampleMR"` delegates to TwoSampleMR's Steiger
    implementation on the harmonised data. Default is `"binary"`.

-   engine:

    Character. `"auto"` prefers TwoSampleMR when installed and otherwise
    falls back to Medusa's internal implementations. Use `"TwoSampleMR"`
    to require delegation, or `"internal"` to force Medusa's built-in
    implementations. Default is `"auto"`.

-   bootstrapSeed:

    Integer or NULL. Seed for the parametric bootstrap used by the
    internal Weighted Median SE estimator. Default is 42, which gives
    reproducible results. Set to NULL to use the ambient RNG state.

</div>

<div class="section level2">

## Value

A named list with class "medusaSensitivity" containing:

-   ivw:

    Data frame with method, beta\_MR, se\_MR, ci\_lower, ci\_upper,
    pval.

-   mrEgger:

    Data frame with beta\_MR, se\_MR, pval, plus intercept,
    intercept\_se, intercept\_pval.

-   weightedMedian:

    Data frame with beta\_MR, se\_MR, ci\_lower, ci\_upper, pval.

-   steiger:

    Data frame with IVW results after Steiger filtering, plus n\_removed
    (number of SNPs removed).

-   leaveOneOut:

    Data frame with snp\_removed, beta\_MR, se\_MR, pval for each SNP
    dropped.

-   summary:

    Data frame comparing all methods side by side.

</div>

<div class="section level2">

## Details

Run MR Sensitivity Analyses

**IVW**: Weighted regression of beta\_ZY on beta\_ZX through the origin,
with weights 1/se\_ZY^2. This is the primary MR estimate assuming all
instruments are valid.

**MR-Egger**: Weighted regression of beta\_ZY on beta\_ZX with an
intercept, after orienting all SNPs so that beta\_ZX is positive. A
non-zero intercept indicates directional pleiotropy. The slope is a
consistent causal estimate even under directional pleiotropy, provided
the InSIDE assumption holds.

**Weighted Median**: Weighted median of per-SNP Wald ratio estimates
using first-order inverse-variance weights for the ratio estimates
(\\(w\_j = \\beta\_{Xj}^2 / \\mathrm{SE}(\\beta\_{Yj})^2\\)). Produces a
consistent estimate if fewer than 50

**Steiger**: Tests whether each SNP explains more variance in the
exposure than the outcome (the expected causal direction). With
`engine = "internal"`, Medusa applies the continuous-trait
correlation-based approximation described by Hemani et al. (2017). With
`engine = "TwoSampleMR"`, Medusa delegates to
`TwoSampleMR::steiger_filtering()` on the harmonised summary data.

**Leave-One-Out**: Drops each SNP in turn and recomputes the IVW
estimate. Identifies influential outlier instruments.

</div>

<div class="section level2">

## References

Bowden, J., Davey Smith, G., & Burgess, S. (2015). Mendelian
randomization with invalid instruments: effect estimation and bias
detection through Egger regression. *International Journal of
Epidemiology*, 44(2), 512-525. doi:10.1093/ije/dyv080. Abstract:
https://pubmed.ncbi.nlm.nih.gov/26050253/ Full text (if available to
you):
https://academic.oup.com/ije/article-pdf/44/2/512/18631079/dyv080.pdf

Bowden, J., et al. (2016). Consistent estimation in Mendelian
randomization with some invalid instruments using a weighted median
estimator. *Genetic Epidemiology*, 40(4), 304-314.
doi:10.1002/gepi.21965. Open access:
https://pmc.ncbi.nlm.nih.gov/articles/PMC4849733/

Hemani, G., Tilling, K., & Davey Smith, G. (2017). Orienting the causal
relationship between imprecisely measured traits using GWAS summary
data. *PLoS Genetics*, 13(11), e1007081.
doi:10.1371/journal.pgen.1007081. Open access:
https://pmc.ncbi.nlm.nih.gov/articles/PMC5711033/

</div>

<div class="section level2">

## See also

<div class="dont-index">

`computeMREstimate`, `fitOutcomeModel`, `generateMRReport`

</div>

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
# Simulate per-SNP estimates
set.seed(42)
nSnps <- 10
betaZX <- rnorm(nSnps, 0.3, 0.05)
betaZY <- 0.5 * betaZX + rnorm(nSnps, 0, 0.02)
perSnp <- data.frame(
  snp_id = paste0("rs", 1:nSnps),
  effect_allele = c("A", "C", "G", "T", "A", "C", "G", "T", "A", "C"),
  other_allele = c("C", "G", "T", "A", "C", "G", "T", "A", "C", "G"),
  eaf = seq(0.1, 0.55, length.out = nSnps),
  beta_ZY = betaZY,
  se_ZY = rep(0.02, nSnps),
  beta_ZX = betaZX,
  se_ZX = rep(0.05, nSnps)
)
results <- runSensitivityAnalyses(
  perSnp,
  outcomeSampleSize = 10000,
  exposureSampleSize = 10000,
  outcomeType = "continuous"
)
#> Running sensitivity analyses with 10 SNPs...
#>   Engine: TwoSampleMR
#> Harmonising exposure (medusa_exposure) and outcome (medusa_outcome)
#> Removing the following SNPs for being palindromic with intermediate allele frequencies:
#> rs10, rs8
#>   IVW...
#> Analysing 'medusa_exposure' on 'medusa_outcome'
#>   MR-Egger...
#> Analysing 'medusa_exposure' on 'medusa_outcome'
#>   Weighted Median...
#> Analysing 'medusa_exposure' on 'medusa_outcome'
#>   Steiger filtering...
#> Analysing 'medusa_exposure' on 'medusa_outcome'
#>   Leave-One-Out...
#> Sensitivity analyses complete.
results$summary
#>                 method    beta_MR      se_MR   ci_lower  ci_upper         pval
#> 1                  IVW 0.49183243 0.03090498  0.4312587 0.5524062 5.039535e-57
#> 2             MR-Egger 0.08879307 0.22182774 -0.3459893 0.5235754 7.028132e-01
#> 3      Weighted Median 0.48415523 0.04343351  0.3990256 0.5692849 7.402441e-29
#> 4 Steiger-filtered IVW 0.37825702 0.04988512  0.2804822 0.4760319 3.387952e-14
```

</div>

</div>

</div>
