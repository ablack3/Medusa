<div id="main" class="col-md-9" role="main">

# Empirical method validation using negative control outcomes

<div class="ref-description section level2">

Implements OHDSI-style empirical calibration for MR estimates. For each
negative control outcome (an outcome with no expected causal
relationship to the exposure), the same allele score model is fitted to
obtain a null MR estimate. The distribution of these null estimates is
then used to assess systematic bias and, optionally, to calibrate the
primary estimate's p-value and confidence interval via the
EmpiricalCalibration package.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
runNegativeControlAnalysis(
  cohortData,
  instrumentTable,
  covariateData = NULL,
  negativeControlColumns = NULL,
  primaryEstimate = NULL,
  modelBackend = "glm"
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   cohortData:

    Data frame. Output of `buildMRCohort`, expected to contain negative
    control outcome columns (named `nc_outcome_<id>` or as specified in
    `negativeControlColumns`).

-   instrumentTable:

    Data frame. Output of `getMRInstruments`.

-   covariateData:

    Optional covariate data for adjustment. Same format as accepted by
    `fitOutcomeModel`.

-   negativeControlColumns:

    Optional character vector of column names in `cohortData` that
    contain negative control outcomes. If NULL, auto-detects columns
    matching `nc_outcome_*`.

-   primaryEstimate:

    Optional list. Output of `computeMREstimate`. If provided, the
    primary estimate is calibrated using the negative control
    distribution.

-   modelBackend:

    Character. Model fitting backend: `"glm"` or `"cyclops"`. Default is
    `"glm"`.

</div>

<div class="section level2">

## Value

A list with class `"medusaNegativeControls"` containing:

-   ncEstimates:

    Data frame with one row per negative control outcome: `outcome_id`,
    `beta_ZY`, `se_ZY`, `beta_MR`, `se_MR`, `pval`, `log_rr`,
    `se_log_rr`.

-   calibration:

    Output of `EmpiricalCalibration::fitSystematicErrorModel()` if the
    package is available, otherwise NULL.

-   calibratedPrimary:

    Named list with `calibratedP`, `calibratedCiLower`,
    `calibratedCiUpper` for the primary estimate, or NULL if
    `primaryEstimate` was not provided or calibration is unavailable.

-   biasDetected:

    Logical. TRUE if the systematic error model indicates the null
    distribution mean is significantly different from zero.

</div>

<div class="section level2">

## Details

Run Negative Control Outcome Analysis for Empirical Calibration

Negative control outcomes should be conditions that are biologically
implausible as consequences of the exposure. A well-calibrated MR
analysis should produce null estimates for these outcomes. Systematic
deviation from the null suggests residual pleiotropy, population
stratification, or other bias.

When EmpiricalCalibration is installed, the function fits a systematic
error model to the negative control estimates and uses it to produce
calibrated p-values and confidence intervals for the primary estimate.
This approach is described in Schuemie et al. (2014, 2018) and is
standard practice in OHDSI observational studies.

</div>

<div class="section level2">

## References

Schuemie, M. J., et al. (2014). Interpreting observational studies: why
empirical calibration is needed to correct p-values. *Statistics in
Medicine*, 33(2), 209-218.

Schuemie, M. J., et al. (2018). Empirical confidence interval
calibration for population-level effect estimation studies in
observational healthcare data. *PNAS*, 115(11), 2571-2577.

</div>

<div class="section level2">

## See also

<div class="dont-index">

`runInstrumentDiagnostics`, `computeMREstimate`, `buildMRCohort`

</div>

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
simData <- simulateMRData(n = 2000, nSnps = 5, seed = 123)
cohortData <- simulateNegativeControlOutcomes(simData$data, nControls = 10)
ncResults <- runNegativeControlAnalysis(
  cohortData = cohortData,
  instrumentTable = simData$instrumentTable
)
#> Running negative control analysis for 10 outcomes...
#>   10 of 10 negative control models fitted successfully.
#>   Package 'EmpiricalCalibration' not installed. Skipping systematic error model. Install with: remotes::install_github('ohdsi/EmpiricalCalibration')
#> Negative control analysis complete.
ncResults$ncEstimates
#>    outcome_id     beta_ZY     se_ZY     beta_MR     se_MR      pval      log_rr
#> 1           1 -0.15699438 0.2728961 -0.49566773 0.8624621 0.5654858 -0.49566773
#> 2           2  0.08931716 0.2343064  0.28199503 0.7400863 0.7031811  0.28199503
#> 3           3 -0.73834441 0.4406275 -2.33112483 1.4029781 0.0966021 -2.33112483
#> 4           4 -0.07535182 0.3249300 -0.23790321 1.0260469 0.8166437 -0.23790321
#> 5           5 -0.45467701 0.4086319 -1.43552094 1.2949883 0.2676371 -1.43552094
#> 6           6 -0.43102268 0.2793367 -1.36083873 0.8882859 0.1255267 -1.36083873
#> 7           7 -0.05332100 0.2544409 -0.16834679 0.8034363 0.8340318 -0.16834679
#> 8           8  0.08689478 0.3269849  0.27434702 1.0325888 0.7904791  0.27434702
#> 9           9 -0.22835265 0.2461892 -0.72096237 0.7793050 0.3548956 -0.72096237
#> 10         10  0.01321082 0.3272555  0.04170962 1.0332269 0.9677995  0.04170962
#>    se_log_rr
#> 1  0.8624621
#> 2  0.7400863
#> 3  1.4029781
#> 4  1.0260469
#> 5  1.2949883
#> 6  0.8882859
#> 7  0.8034363
#> 8  1.0325888
#> 9  0.7793050
#> 10 1.0332269
```

</div>

</div>

</div>
