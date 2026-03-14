<div id="main" class="col-md-9" role="main">

# Shiny dashboard for MR diagnostics and results

<div class="ref-description section level2">

Opens an interactive bslib dashboard displaying instrument diagnostics
(F-statistics, PheWAS, allele frequencies, missingness, heterogeneity)
and MR results (primary estimate, sensitivity analyses, scatter/funnel/
leave-one-out plots). Styled to match the Medusa package website.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
launchResultsExplorer(
  mrEstimate,
  sensitivityResults = NULL,
  diagnosticResults = NULL,
  combinedProfile = NULL,
  siteProfileList = NULL,
  instrumentTable = NULL,
  perSnpEstimates = NULL,
  negativeControlResults = NULL,
  launch.browser = TRUE,
  port = NULL
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   mrEstimate:

    Output of `computeMREstimate`.

-   sensitivityResults:

    Output of `runSensitivityAnalyses`. Can be NULL.

-   diagnosticResults:

    Output of `runInstrumentDiagnostics`. Can be NULL.

-   combinedProfile:

    Output of `poolLikelihoodProfiles`.

-   siteProfileList:

    Optional named list of site profile objects from `fitOutcomeModel`.

-   instrumentTable:

    Output of `getMRInstruments`.

-   perSnpEstimates:

    Optional data frame of per-SNP summary statistics (output of
    `fitOutcomeModel` with `analysisType = "perSNP"`). Needed for
    scatter and funnel plots.

-   negativeControlResults:

    Output of `runNegativeControlAnalysis`. Can be NULL.

-   launch.browser:

    Logical. Whether to open a browser window. Default is TRUE.

-   port:

    Integer or NULL. Port for the Shiny server. NULL uses a random
    available port.

</div>

<div class="section level2">

## Value

Invisible NULL. Launches a Shiny app.

</div>

<div class="section level2">

## Details

Launch Interactive Results Explorer

</div>

<div class="section level2">

## See also

<div class="dont-index">

`computeMREstimate`, `runSensitivityAnalyses`,
`runInstrumentDiagnostics`, `generateMRReport`

</div>

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
if (FALSE) { # \dontrun{
launchResultsExplorer(
  mrEstimate = estimate,
  sensitivityResults = sensitivity,
  diagnosticResults = diagnostics,
  combinedProfile = combined,
  instrumentTable = instruments,
  perSnpEstimates = perSnpEst
)
} # }
```

</div>

</div>

</div>
