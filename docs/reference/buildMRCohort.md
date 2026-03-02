# Site-level cohort extraction from OMOP CDM

Runs locally at each OMOP CDM site. Executes parameterized SQL to
extract outcome cohort data and genotype data, performs allele
harmonization to ensure genotype coding matches the instrument table,
and returns a local R data frame suitable for downstream analysis. The
returned data never leaves the site.

The function queries: PERSON (age, sex), CONDITION_OCCURRENCE (outcome),
OBSERVATION_PERIOD (eligibility), and the genomic linkage table
(genotypes). All SQL is rendered and translated via SqlRender for
cross-dialect compatibility.

## Usage

``` r
buildMRCohort(
  connectionDetails,
  cdmDatabaseSchema,
  cohortDatabaseSchema,
  cohortTable,
  outcomeCohortId,
  instrumentTable,
  genomicLinkageSchema,
  genomicLinkageTable,
  indexDateOffset = 0,
  washoutPeriod = 365,
  excludePriorOutcome = TRUE,
  genomicPersonIdColumn = "person_id"
)
```

## Arguments

- connectionDetails:

  A `DatabaseConnector::connectionDetails` object specifying the
  database connection.

- cdmDatabaseSchema:

  Character. Schema containing the OMOP CDM tables.

- cohortDatabaseSchema:

  Character. Schema containing the cohort table.

- cohortTable:

  Character. Name of the cohort table.

- outcomeCohortId:

  Integer. Cohort definition ID for the outcome of interest (e.g.,
  incident colorectal cancer).

- instrumentTable:

  Data frame. Output of [`getMRInstruments`](getMRInstruments.md) or
  [`createInstrumentTable`](createInstrumentTable.md) containing the
  instrument SNPs.

- genomicLinkageSchema:

  Character. Schema containing the genomic linkage table.

- genomicLinkageTable:

  Character. Name of the table linking person_id to SNP genotypes. Must
  contain columns for person identifier, snp_id, and genotype (coded as
  0/1/2 count of effect alleles).

- indexDateOffset:

  Integer. Days offset from cohort start date for defining the index
  date. Default is 0.

- washoutPeriod:

  Integer. Minimum days of prior observation required before index date.
  Default is 365.

- excludePriorOutcome:

  Logical. If TRUE, persons with the outcome before their index date are
  excluded. Default is TRUE.

- genomicPersonIdColumn:

  Character. Name of the person identifier column in the genomic linkage
  table if it differs from "person_id". Default is "person_id".

## Value

A data frame with one row per person and columns:

- person_id:

  Integer person identifier.

- outcome:

  Integer 0/1 outcome status.

- age_at_index:

  Age in years at index date.

- gender_concept_id:

  OMOP concept ID for gender.

- index_date:

  Date of cohort entry.

- snp_1 through snp_K:

  Integer genotype values (0, 1, or 2) for each instrument SNP, coded as
  count of the effect allele after harmonization. Missing genotypes are
  NA, not 0.

## Details

Build Study Cohort for Mendelian Randomization at a Single Site

Genotype coding: genotypes are coded as 0, 1, 2 representing the count
of effect alleles. The function performs allele harmonization by
comparing the effect allele in the instrument table to the allele coding
in the genotype data. If alleles are swapped, genotypes are flipped (2 -
genotype).

## References

Hripcsak, G., et al. (2015). Observational Health Data Sciences and
Informatics (OHDSI): Opportunities for Observational Researchers.
*Studies in Health Technology and Informatics*, 216, 574-578.

## See also

[`getMRInstruments`](getMRInstruments.md),
[`harmonizeAlleles`](harmonizeAlleles.md),
[`buildMRCovariates`](buildMRCovariates.md)

## Examples

``` r
if (FALSE) { # \dontrun{
connectionDetails <- DatabaseConnector::createConnectionDetails(
  dbms = "postgresql",
  server = "localhost/ohdsi",
  user = "user",
  password = "password"
)
instruments <- getMRInstruments("ieu-a-1119")
cohort <- buildMRCohort(
  connectionDetails = connectionDetails,
  cdmDatabaseSchema = "cdm",
  cohortDatabaseSchema = "results",
  cohortTable = "cohort",
  outcomeCohortId = 1234,
  instrumentTable = instruments,
  genomicLinkageSchema = "genomics",
  genomicLinkageTable = "genotype_data"
)
} # }
```
