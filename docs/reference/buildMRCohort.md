<div id="main" class="col-md-9" role="main">

# Site-level cohort extraction from OMOP CDM

<div class="ref-description section level2">

Runs locally at each OMOP CDM site. Executes parameterized SQL to
extract outcome cohort data and genotype data, performs allele
harmonization to ensure genotype coding matches the instrument table,
and returns a local R data frame suitable for downstream analysis. The
returned data never leaves the site.

Genotype data is extracted from the **VARIANT\_OCCURRENCE** table
defined by the OMOP CDM Genomic Extension. The minimal required columns
from that table are:

-   `person_id` — links variants to persons

-   `rs_id` — dbSNP rs identifier for the variant

-   `genotype` — genotype call (VCF-style "0/0", "0/1", "1/1" or plain
    integer "0", "1", "2")

Additionally, `reference_allele` and `alternate_allele` are used for
allele harmonization when available.

The function also queries: PERSON (age, sex), CONDITION\_OCCURRENCE
(outcome), and OBSERVATION\_PERIOD (eligibility). All SQL is rendered
and translated via SqlRender for cross-dialect compatibility.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
buildMRCohort(
  connectionDetails,
  cdmDatabaseSchema,
  cohortDatabaseSchema,
  cohortTable,
  outcomeCohortId,
  instrumentTable,
  genomicDatabaseSchema = cdmDatabaseSchema,
  indexDateOffset = 0,
  washoutPeriod = 365,
  excludePriorOutcome = TRUE
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   connectionDetails:

    A `DatabaseConnector::connectionDetails` object specifying the
    database connection.

-   cdmDatabaseSchema:

    Character. Schema containing the OMOP CDM tables.

-   cohortDatabaseSchema:

    Character. Schema containing the cohort table.

-   cohortTable:

    Character. Name of the cohort table.

-   outcomeCohortId:

    Integer. Cohort definition ID for the outcome of interest (e.g.,
    incident colorectal cancer).

-   instrumentTable:

    Data frame. Output of `getMRInstruments` or `createInstrumentTable`
    containing the instrument SNPs.

-   genomicDatabaseSchema:

    Character. Schema containing the VARIANT\_OCCURRENCE table from the
    OMOP Genomic Extension. Defaults to `cdmDatabaseSchema`.

-   indexDateOffset:

    Integer. Days offset from cohort start date for defining the index
    date. Default is 0.

-   washoutPeriod:

    Integer. Minimum days of prior observation required before index
    date. Default is 365.

-   excludePriorOutcome:

    Logical. If TRUE, persons with the outcome before their index date
    are excluded. Default is TRUE.

</div>

<div class="section level2">

## Value

A data frame with one row per person and columns:

-   person\_id:

    Integer person identifier.

-   outcome:

    Integer 0/1 outcome status.

-   age\_at\_index:

    Age in years at index date.

-   gender\_concept\_id:

    OMOP concept ID for gender.

-   index\_date:

    Date of cohort entry.

-   snp\_&lt;sanitized rsID&gt;:

    Integer genotype values (0, 1, or 2) for each instrument SNP, coded
    as count of the effect allele after harmonization. Missing genotypes
    are NA, not 0.

</div>

<div class="section level2">

## Details

Build Study Cohort for Mendelian Randomization at a Single Site

Genotype coding: genotypes are coded as 0, 1, 2 representing the count
of effect alleles. The function performs allele harmonization by
comparing the effect allele in the instrument table to the allele coding
in the genotype data (alternate allele from VARIANT\_OCCURRENCE). If
alleles are swapped, genotypes are flipped (2 - genotype) and instrument
beta\_ZX is negated.

</div>

<div class="section level2">

## References

OHDSI Genomic CDM: <https://github.com/OHDSI/Genomic-CDM>

Hripcsak, G., et al. (2015). Observational Health Data Sciences and
Informatics (OHDSI): Opportunities for Observational Researchers.
*Studies in Health Technology and Informatics*, 216, 574-578.

</div>

<div class="section level2">

## See also

<div class="dont-index">

`getMRInstruments`, `harmonizeAlleles`, `buildMRCovariates`

</div>

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

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
  genomicDatabaseSchema = "genomics"
)
} # }
```

</div>

</div>

</div>
