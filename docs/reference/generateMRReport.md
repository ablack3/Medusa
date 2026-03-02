# Self-contained HTML report for MR analysis

Produces a self-contained HTML report (single file, no external
dependencies) summarizing the entire MR analysis including instrument
summary, likelihood profile plots, main MR result, sensitivity analyses,
PheWAS diagnostics, and site contributions. Suitable for sharing with
non-technical stakeholders.

## Usage

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

## Arguments

- mrEstimate:

  Output of [`computeMREstimate`](computeMREstimate.md).

- sensitivityResults:

  Output of [`runSensitivityAnalyses`](runSensitivityAnalyses.md). Can
  be NULL if no sensitivity analyses were run.

- diagnosticResults:

  Output of [`runInstrumentDiagnostics`](runInstrumentDiagnostics.md).
  Can be NULL if no diagnostics were run.

- combinedProfile:

  Output of [`poolLikelihoodProfiles`](poolLikelihoodProfiles.md).

- siteProfileList:

  Named list of site profile objects from
  [`fitOutcomeModel`](fitOutcomeModel.md).

- instrumentTable:

  Output of [`getMRInstruments`](getMRInstruments.md).

- exposureLabel:

  Character. Human-readable name for the exposure. Default is
  "Exposure".

- outcomeLabel:

  Character. Human-readable name for the outcome. Default is "Outcome".

- outputPath:

  Character. Path for the output HTML file. Default is
  "./Medusa_report.html".

## Value

Character string with the path to the generated report (invisibly).

## Details

Generate Mendelian Randomization Analysis Report

## See also

[`computeMREstimate`](computeMREstimate.md),
[`runSensitivityAnalyses`](runSensitivityAnalyses.md),
[`runInstrumentDiagnostics`](runInstrumentDiagnostics.md)

## Examples

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
