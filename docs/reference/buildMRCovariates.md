<div id="main" class="col-md-9" role="main">

# Covariate assembly using FeatureExtraction

<div class="ref-description section level2">

Leverages the OHDSI FeatureExtraction package to assemble a rich
covariate matrix from OMOP CDM data. This covariate matrix serves two
purposes: (1) adjustment covariates in the Cyclops outcome model to
control for confounders and population stratification, and (2) the full
phenome for the instrument PheWAS diagnostic that checks for pleiotropic
associations.

If ancestry principal components (PCs) are available, they are merged
into the covariate matrix. Ancestry PCs are critical for controlling
population stratification in genetic association studies.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
buildMRCovariates(
  connectionDetails,
  cdmDatabaseSchema,
  cohortDatabaseSchema,
  cohortTable,
  outcomeCohortId,
  covariateSettings = NULL,
  ancestryPCsTable = NULL,
  ancestryPCsSchema = NULL,
  numAncestryPCs = 10
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   connectionDetails:

    A `DatabaseConnector::connectionDetails` object.

-   cdmDatabaseSchema:

    Character. Schema containing OMOP CDM tables.

-   cohortDatabaseSchema:

    Character. Schema containing the cohort table.

-   cohortTable:

    Character. Name of the cohort table.

-   outcomeCohortId:

    Integer. Cohort definition ID for the outcome.

-   covariateSettings:

    A FeatureExtraction covariate settings object. If NULL (default),
    uses `createDefaultMRCovariateSettings`.

-   ancestryPCsTable:

    Character or NULL. Name of a table containing person\_id and
    ancestry principal components (PC1 through PC\_K). If NULL, ancestry
    PCs are not included.

-   ancestryPCsSchema:

    Character or NULL. Schema containing the ancestry PCs table.

-   numAncestryPCs:

    Integer. Number of ancestry PCs to include (1 through this value).
    Default is 10.

</div>

<div class="section level2">

## Value

A list with class "medusaCovariateData" containing:

-   covariates:

    Data frame with person\_id and covariate columns.

-   covariateRef:

    Data frame mapping covariate IDs to names and domains.

-   ancestryPCs:

    Data frame of ancestry PCs if provided, NULL otherwise.

-   settings:

    The covariate settings object used.

</div>

<div class="section level2">

## Details

Assemble Covariate Matrix for Mendelian Randomization

Default covariate settings include: conditions in 365-day lookback
(binary), drug exposures in 365-day lookback (binary), most recent
measurement values, and demographics (age group, sex, index year). These
are assembled using standard FeatureExtraction covariate setting
objects.

</div>

<div class="section level2">

## See also

<div class="dont-index">

`createDefaultMRCovariateSettings`, `buildMRCohort`,
`runInstrumentDiagnostics`

</div>

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
if (FALSE) { # \dontrun{
covData <- buildMRCovariates(
  connectionDetails = connDetails,
  cdmDatabaseSchema = "cdm",
  cohortDatabaseSchema = "results",
  cohortTable = "cohort",
  outcomeCohortId = 1234
)
} # }
```

</div>

</div>

</div>
