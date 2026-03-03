<div id="main" class="col-md-9" role="main">

# Single-site cross-check against standard summary-data MR

<div class="ref-description section level2">

Runs Medusa's single-site profile-likelihood estimate and compares it
with TwoSampleMR estimates computed from the same site's per-SNP summary
statistics. This is a calibration check: agreement is expected when the
profile is close to quadratic and the allele-score approximation aligns
with standard summarized-data MR assumptions.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
validateAgainstTwoSampleMR(
  siteProfile,
  instrumentTable,
  methods = c("IVW", "MREgger", "WeightedMedian"),
  outcomeSampleSize = NULL,
  exposureSampleSize = NULL,
  outcomeType = "binary",
  ciLevel = 0.95
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   siteProfile:

    A single `medusaSiteProfile` object produced by `fitOutcomeModel`
    with `analysisType = "perSNP"`.

-   instrumentTable:

    Data frame. Output of `getMRInstruments`. Used both for the Medusa
    Wald ratio and, if needed, to backfill allele columns missing from
    older `perSnpEstimates` objects.

-   methods:

    Character vector of summary-data methods to compare. Default is
    `c("IVW", "MREgger", "WeightedMedian")`.

-   outcomeSampleSize:

    Optional integer. Passed to `runSensitivityAnalyses`.

-   exposureSampleSize:

    Optional integer. Passed to `runSensitivityAnalyses`.

-   outcomeType:

    Character. Outcome summary-statistic scale represented by
    `siteProfile$perSnpEstimates`. Default is `"binary"`.

-   ciLevel:

    Numeric. Confidence interval level for the Medusa estimate. Default
    is 0.95.

</div>

<div class="section level2">

## Value

A data frame with one row for the Medusa profile-likelihood estimate and
one row for each requested TwoSampleMR method. The `delta_vs_medusa`
column reports each summary-data estimate minus the Medusa estimate.

</div>

<div class="section level2">

## Details

Validate a Single-Site Medusa Estimate Against TwoSampleMR

</div>

<div class="section level2">

## See also

<div class="dont-index">

`fitOutcomeModel`, `computeMREstimate`, `runSensitivityAnalyses`

</div>

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
simData <- simulateMRData(n = 3000, nSnps = 5, trueEffect = 0.3)
siteProfile <- fitOutcomeModel(
  cohortData = simData$data,
  instrumentTable = simData$instrumentTable,
  analysisType = "perSNP",
  betaGrid = seq(-2, 2, by = 0.05)
)
#> Fitting outcome model at site 'site_1' (1725 cases, 1275 controls)...
#> Site 'site_1': beta_ZY_hat = 0.0426 (SE = 0.1182).
comparison <- validateAgainstTwoSampleMR(siteProfile, simData$instrumentTable)
#> Pooling profile likelihoods from 1 site(s)...
#> Pooling complete: 1 sites, 1725 total cases, 1275 total controls.
#> MR estimate: beta = 0.2369 (95% CI: -0.7107, 1.1845), p = 6.72e-01
#> Odds ratio: 1.267 (95% CI: 0.491, 3.269)
#> Running sensitivity analyses with 5 SNPs...
#>   Engine: TwoSampleMR
#> Harmonising exposure (medusa_exposure) and outcome (medusa_outcome)
#>   IVW...
#> Analysing 'medusa_exposure' on 'medusa_outcome'
#>   MR-Egger...
#> Analysing 'medusa_exposure' on 'medusa_outcome'
#>   Weighted Median...
#> Analysing 'medusa_exposure' on 'medusa_outcome'
#> Sensitivity analyses complete.
comparison
#>        source             method   beta_MR     se_MR     ci_lower  ci_upper
#> 1      Medusa Profile likelihood 0.2369024 0.5603015 -0.710707109 1.1845118
#> 2 TwoSampleMR                IVW 0.2194567 0.1167514 -0.009376005 0.4482894
#> 3 TwoSampleMR           MR-Egger 0.6304650 0.2779381  0.085706258 1.1752237
#> 4 TwoSampleMR    Weighted Median 0.2426033 0.1391244 -0.030080435 0.5152870
#>         pval delta_vs_medusa
#> 1 0.67243220     0.000000000
#> 2 0.06014999    -0.017445678
#> 3 0.10808801     0.393562607
#> 4 0.08119616     0.005700936
```

</div>

</div>

</div>
