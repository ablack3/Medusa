<div id="main" class="col-md-9" role="main">

# Export a Site Profile to CSV

<div class="ref-description section level2">

Writes a site profile (output of `fitOutcomeModel`) to two
human-readable CSV files: one containing the profile log-likelihood grid
and one containing site metadata. These CSV files are the artifacts
shared between sites and the coordinator in a federated Medusa analysis.

CSV is used instead of binary formats so that every value leaving a site
is human-readable and auditable.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
exportSiteProfile(profile, outputDir = ".", prefix = "medusa")
```

</div>

</div>

<div class="section level2">

## Arguments

-   profile:

    A site profile object (output of `fitOutcomeModel`).

-   outputDir:

    Character. Directory to write files to. Default is current working
    directory.

-   prefix:

    Character. Filename prefix. Default is "medusa".

</div>

<div class="section level2">

## Value

A named character vector with the paths to the written files
(invisibly). Names are "profile" and "metadata".

</div>

<div class="section level2">

## See also

<div class="dont-index">

`importSiteProfile`, `fitOutcomeModel`

</div>

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
simData <- simulateMRData(n = 500, nSnps = 3, trueEffect = 0.3)
profile <- fitOutcomeModel(
  cohortData = simData$data,
  covariateData = NULL,
  instrumentTable = simData$instrumentTable,
  betaGrid = seq(-2, 2, by = 0.1),
  siteId = "example_site"
)
#> Fitting outcome model at site 'example_site' (281 cases, 219 controls)...
#> Site 'example_site': beta_ZY_hat = 0.3486 (SE = 0.2158).
if (FALSE) { # \dontrun{
exportSiteProfile(profile, outputDir = tempdir())
} # }
```

</div>

</div>

</div>
