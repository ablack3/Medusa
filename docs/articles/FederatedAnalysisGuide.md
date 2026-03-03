<div id="main" class="col-md-9" role="main">

# Federated Analysis Guide for Network Coordinators

<div class="section level2">

## Overview

This guide walks OHDSI network coordinators through setting up and
running a federated Medusa analysis across multiple sites. The key
principle: **only profile log-likelihood vectors leave each site** — no
individual-level data is shared.

</div>

<div class="section level2">

## What Data Leaves Each Site

| Shared with coordinator                                | NOT shared            |
|--------------------------------------------------------|-----------------------|
| Log-likelihood profile (numeric vector, \~600 numbers) | Individual genotypes  |
| Number of cases and controls                           | Person-level outcomes |
| Diagnostic flags (logical values)                      | Covariate values      |
| Site identifier                                        | Demographics          |

The profile vector is a smooth curve that represents the aggregate
statistical evidence at a site. It cannot be reverse-engineered to
identify individuals.

</div>

<div class="section level2">

## Site Setup Requirements

<div class="section level3">

### OMOP CDM Requirements

-   OMOP CDM version 5.3 or 5.4
-   Standard tables: PERSON, OBSERVATION\_PERIOD, CONDITION\_OCCURRENCE
-   A cohort table with the outcome cohort pre-defined (e.g., via ATLAS)

</div>

<div class="section level3">

### Genomic Linkage Table

Each site needs a table linking person\_id to SNP genotypes:

<div id="cb1" class="sourceCode">

``` sql
CREATE TABLE genomics.genotype_data (
  person_id BIGINT NOT NULL,
  snp_id VARCHAR(50) NOT NULL,
  genotype INTEGER NOT NULL  -- 0, 1, or 2 (count of effect alleles)
);
```

</div>

If the person identifier column has a different name (e.g.,
`subject_id`), specify it via the `genomicPersonIdColumn` parameter.

</div>

<div class="section level3">

### R Package Dependencies

<div id="cb2" class="sourceCode">

``` r
# Required at each site
remotes::install_github("OHDSI/Medusa")

# This installs: DatabaseConnector, SqlRender, Cyclops, FeatureExtraction
```

</div>

</div>

<div class="section level3">

### Ancestry Principal Components (Optional but Recommended)

A table with ancestry PCs controls for population stratification:

<div id="cb3" class="sourceCode">

``` sql
CREATE TABLE genomics.ancestry_pcs (
  person_id BIGINT NOT NULL,
  pc1 FLOAT, pc2 FLOAT, pc3 FLOAT, pc4 FLOAT, pc5 FLOAT,
  pc6 FLOAT, pc7 FLOAT, pc8 FLOAT, pc9 FLOAT, pc10 FLOAT
);
```

</div>

</div>

</div>

<div class="section level2">

## Template Site Analysis Script

Send this script to each participating site, customized with their local
connection details:

<div id="cb4" class="sourceCode">

``` r
# ============================================================
# Medusa Site Analysis Script
# Study: [YOUR STUDY NAME]
# Site: [SITE NAME]
# Date: [DATE]
# ============================================================

library(Medusa)

# ----- Site-specific configuration -----
connectionDetails <- DatabaseConnector::createConnectionDetails(
  dbms = "postgresql",              # Change to your DBMS
  server = "localhost/ohdsi",       # Change to your server
  user = "ohdsi_user",             # Change to your user
  password = keyring::key_get("ohdsi_db")  # Use secure credential storage
)

cdmDatabaseSchema <- "cdm"           # Change to your CDM schema
cohortDatabaseSchema <- "results"     # Change to your results schema
cohortTable <- "cohort"               # Change if different
outcomeCohortId <- 1234               # Provided by coordinator
genomicLinkageSchema <- "genomics"    # Change to your genomics schema
genomicLinkageTable <- "genotype_data"
siteId <- "site_A"                    # Unique identifier for this site

# ----- Load instrument table (provided by coordinator) -----
instrumentTable <- read.csv("instruments.csv", stringsAsFactors = FALSE)

# ----- Step 1: Build cohort -----
cohortData <- buildMRCohort(
  connectionDetails = connectionDetails,
  cdmDatabaseSchema = cdmDatabaseSchema,
  cohortDatabaseSchema = cohortDatabaseSchema,
  cohortTable = cohortTable,
  outcomeCohortId = outcomeCohortId,
  instrumentTable = instrumentTable,
  genomicLinkageSchema = genomicLinkageSchema,
  genomicLinkageTable = genomicLinkageTable,
  washoutPeriod = 365,
  excludePriorOutcome = TRUE
)

# ----- Step 2: Build covariates -----
covariateData <- buildMRCovariates(
  connectionDetails = connectionDetails,
  cdmDatabaseSchema = cdmDatabaseSchema,
  cohortDatabaseSchema = cohortDatabaseSchema,
  cohortTable = cohortTable,
  outcomeCohortId = outcomeCohortId
)

# ----- Step 3: Run diagnostics -----
diagnostics <- runInstrumentDiagnostics(
  cohortData = cohortData,
  covariateData = covariateData,
  instrumentTable = instrumentTable
)

# ----- Step 4: Fit outcome model -----
# Use the SAME betaGrid as all other sites
betaGrid <- seq(-3, 3, by = 0.01)

profile <- fitOutcomeModel(
  cohortData = cohortData,
  covariateData = covariateData,
  instrumentTable = instrumentTable,
  betaGrid = betaGrid,
  siteId = siteId,
  analysisType = "alleleScore"
)

# Also fit per-SNP models for sensitivity analyses
profilePerSnp <- fitOutcomeModel(
  cohortData = cohortData,
  covariateData = covariateData,
  instrumentTable = instrumentTable,
  betaGrid = betaGrid,
  siteId = siteId,
  analysisType = "perSNP"
)

# ----- Step 5: Export and share -----
# ONLY these CSV files leave the site:
exportSiteProfile(profile, outputDir = ".", prefix = "medusa")

message("Analysis complete. Share the CSV files with the coordinator.")
```

</div>

</div>

<div class="section level2">

## Coordinator Pooling Script

After collecting profile CSV files from all sites:

<div id="cb5" class="sourceCode">

``` r
library(Medusa)

# Load instrument table
instrumentTable <- read.csv("instruments.csv", stringsAsFactors = FALSE)

# Import site profiles from CSV files
siteProfiles <- list(
  site_A = importSiteProfile("medusa_profile_site_A.csv"),
  site_B = importSiteProfile("medusa_profile_site_B.csv"),
  site_C = importSiteProfile("medusa_profile_site_C.csv")
)

# Pool
combined <- poolLikelihoodProfiles(siteProfiles)

# Estimate
estimate <- computeMREstimate(combined, instrumentTable)

# Sensitivity analyses (if per-SNP profiles available)
# ... construct perSnpEstimates from per-SNP profiles ...

# Report
generateMRReport(
  mrEstimate = estimate,
  combinedProfile = combined,
  siteProfileList = siteProfiles,
  instrumentTable = instrumentTable,
  exposureLabel = "IL-6 signaling",
  outcomeLabel = "Colorectal cancer"
)
```

</div>

</div>

<div class="section level2">

## Secure File Transfer

Profile CSV files should be transferred securely between sites and
coordinator. Options include:

-   SFTP with encrypted credentials
-   Institutional secure file sharing (e.g., Box, SharePoint with
    encryption)
-   OHDSI network file transfer protocols (if available)

The CSV files are human-readable and typically very small (&lt; 100 KB)
as they contain only numeric grid values and summary statistics. Using
CSV ensures that every value leaving a site can be inspected and audited
before transfer.

</div>

<div class="section level2">

## Handling Different OMOP Versions

Sites running OMOP CDM v5.3 vs v5.4 can participate together. The SQL
templates in Medusa use only core CDM tables that are consistent across
versions. If a site uses non-standard table names, these can be
configured via function parameters.

</div>

<div class="section level2">

## Troubleshooting

<div class="section level3">

### “No persons have genotype data”

-   Verify the genomic linkage table exists and has data
-   Check that `person_id` column names match (use
    `genomicPersonIdColumn` parameter)
-   Ensure the cohort and genotype table share person\_id values

</div>

<div class="section level3">

### “Profile likelihood is flat”

-   Instruments may be too weak at this site
-   Check F-statistics in diagnostics
-   Verify that genotype data is coded correctly (0/1/2)

</div>

<div class="section level3">

### “MLE is at grid boundary”

-   Expand the `betaGrid` range (e.g., `seq(-5, 5, by = 0.01)`)
-   Ensure the same `betaGrid` is used at all sites

</div>

<div class="section level3">

### Weak instrument warning (F &lt; 10)

-   The instrument explains very little outcome variance at this site
-   Consider excluding weak instruments or using more instruments
-   Results may be biased toward the null

</div>

</div>

</div>
