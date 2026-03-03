<div id="main" class="col-md-9" role="main">

# Self-contained HTML report for MR analysis

<div class="ref-description section level2">

Produces a self-contained HTML report (single file, no external
dependencies) summarizing the entire MR analysis including instrument
summary, likelihood profile plots, main MR result, sensitivity analyses,
PheWAS diagnostics, and site contributions. Suitable for sharing with
non-technical stakeholders.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
generateMRReport(
  mrEstimate,
  sensitivityResults = NULL,
  diagnosticResults = NULL,
  combinedProfile,
  siteProfileList = NULL,
  instrumentTable = NULL,
  exposureLabel = "Exposure",
  outcomeLabel = "Outcome",
  outputPath = "./Medusa_report.html"
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   mrEstimate:

    Output of `computeMREstimate`.

-   sensitivityResults:

    Output of `runSensitivityAnalyses`. Can be NULL if no sensitivity
    analyses were run.

-   diagnosticResults:

    Output of `runInstrumentDiagnostics`. Can be NULL if no diagnostics
    were run.

-   combinedProfile:

    Output of `poolLikelihoodProfiles`.

-   siteProfileList:

    Named list of site profile objects from `fitOutcomeModel`.

-   instrumentTable:

    Output of `getMRInstruments`.

-   exposureLabel:

    Character. Human-readable name for the exposure. Default is
    "Exposure".

-   outcomeLabel:

    Character. Human-readable name for the outcome. Default is
    "Outcome".

-   outputPath:

    Character. Path for the output HTML file. Default is
    "./Medusa\_report.html".

</div>

<div class="section level2">

## Value

Character string with the path to the generated report (invisibly).

</div>

<div class="section level2">

## Details

Generate Mendelian Randomization Analysis Report

</div>

<div class="section level2">

## See also

<div class="dont-index">

`computeMREstimate`, `runSensitivityAnalyses`,
`runInstrumentDiagnostics`

</div>

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
if (FALSE) { # \dontrun{
generateMRReport(
  mrEstimate = estimate,
  sensitivityResults = sensitivity,
  diagnosticResults = diagnostics,
  combinedProfile = combined,
  siteProfileList = siteProfiles,
  instrumentTable = instruments,
  exposureLabel = "IL-6 receptor levels",
  outcomeLabel = "Colorectal cancer"
)
} # }
```

</div>

</div>

</div>
