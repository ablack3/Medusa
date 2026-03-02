# Covariate assembly using FeatureExtraction

Leverages the OHDSI FeatureExtraction package to assemble a rich
covariate matrix from OMOP CDM data. This covariate matrix serves two
purposes: (1) adjustment covariates in the Cyclops outcome model to
control for confounders and population stratification, and (2) the full
phenome for the instrument PheWAS diagnostic that checks for pleiotropic
associations.

If ancestry principal components (PCs) are available, they are merged
into the covariate matrix. Ancestry PCs are critical for controlling
population stratification in genetic association studies.

## Usage

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

## Arguments

- connectionDetails:

  A `DatabaseConnector::connectionDetails` object.

- cdmDatabaseSchema:

  Character. Schema containing OMOP CDM tables.

- cohortDatabaseSchema:

  Character. Schema containing the cohort table.

- cohortTable:

  Character. Name of the cohort table.

- outcomeCohortId:

  Integer. Cohort definition ID for the outcome.

- covariateSettings:

  A FeatureExtraction covariate settings object. If NULL (default), uses
  [`createDefaultMRCovariateSettings`](createDefaultMRCovariateSettings.md).

- ancestryPCsTable:

  Character or NULL. Name of a table containing person_id and ancestry
  principal components (PC1 through PC_K). If NULL, ancestry PCs are not
  included.

- ancestryPCsSchema:

  Character or NULL. Schema containing the ancestry PCs table.

- numAncestryPCs:

  Integer. Number of ancestry PCs to include (1 through this value).
  Default is 10.

## Value

A list with class "medusaCovariateData" containing:

- covariates:

  Data frame with person_id and covariate columns.

- covariateRef:

  Data frame mapping covariate IDs to names and domains.

- ancestryPCs:

  Data frame of ancestry PCs if provided, NULL otherwise.

- settings:

  The covariate settings object used.

## Details

Assemble Covariate Matrix for Mendelian Randomization

Default covariate settings include: conditions in 365-day lookback
(binary), drug exposures in 365-day lookback (binary), most recent
measurement values, and demographics (age group, sex, index year). These
are assembled using standard FeatureExtraction covariate setting
objects.

## See also

[`createDefaultMRCovariateSettings`](createDefaultMRCovariateSettings.md),
[`buildMRCohort`](buildMRCohort.md),
[`runInstrumentDiagnostics`](runInstrumentDiagnostics.md)

## Examples

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
