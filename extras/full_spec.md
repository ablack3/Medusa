# Medusa R Package — AI Coder Specification

**Federated Mendelian Randomization on OMOP CDM**  
Version 1.0 | March 2026

---

## Table of Contents

1. [Project Overview](#1-project-overview)
2. [Package Architecture](#2-package-architecture)
3. [Package File Structure](#3-package-file-structure)
4. [Package Dependencies](#4-package-dependencies)
5. [DESCRIPTION File](#5-description-file)
6. [Testing Requirements](#6-testing-requirements)
7. [Vignette Requirements](#7-vignette-requirements)
8. [Documentation Requirements](#8-documentation-requirements)
9. [Coding Standards and Style](#9-coding-standards-and-style)
10. [Error Handling and Edge Cases](#10-error-handling-and-edge-cases)
11. [SQL Template Requirements](#11-sql-template-requirements)
12. [Scientific Validation](#12-scientific-validation)
13. [Deliverables Checklist](#13-deliverables-checklist)
14. [Implementation Notes for the AI Coder](#14-implementation-notes-for-the-ai-coder)

---

## 1. Project Overview

You are building **Medusa**, an open-source R package that implements two-sample Mendelian Randomization (MR) natively within the OHDSI ecosystem, using the OMOP Common Data Model (CDM) as the data substrate. This package enables federated causal inference across distributed health networks without requiring individual-level data to leave any site.

The package bridges two previously separate worlds: the OHDSI federated network infrastructure and the Mendelian Randomization methodology used in drug target validation. It is intended for use by epidemiologists, statisticians, and pharmacoepidemiology researchers working in pharma, biotech, and academic settings.

The core methodological innovation is **one-shot federated pooling via profile likelihood aggregation** — each site computes a log-likelihood profile across a grid of parameter values and shares only that vector of numbers. The coordinator sums profiles across sites to obtain a pooled estimate without any iterative communication protocol and without any individual-level data leaving any site.

### 1.1 Key Scientific Context

Mendelian Randomization uses genetic variants as instrumental variables to estimate causal effects of exposures on outcomes, analogous to a natural randomized trial. In two-sample MR, the SNP-exposure association (`beta_ZX`) is obtained from a published GWAS, and the SNP-outcome association (`beta_ZY`) is estimated from the analyst's own cohort. The causal estimate is the Wald ratio:

```
beta_MR = beta_ZY / beta_ZX
```

The OMOP CDM provides standardized clinical data (diagnoses, drugs, labs, procedures) across thousands of health systems globally. Genomic data is increasingly being linked to OMOP records via `person_id` in biobank-linked networks such as All of Us, Vanderbilt BioVU, Michigan Genomics Initiative, and VA Million Veteran Program. Medusa targets these linked sites for the outcome model while drawing exposure-side GWAS summary statistics from public databases.

### 1.2 Overall Data Flow

```
PUBLIC GWAS DATABASE              EACH OMOP SITE
(OpenGWAS API)                    (runs locally)
      |                                  |
      |                                  |
  beta_ZX, se_ZX              1. Cohort extraction
  SNP-exposure                2. Covariate assembly
  association                 3. Cyclops outcome model
      |                       4. Grid likelihood profile
      |                                  |
      └─────────────┬────────────────────┘
                    |
           COORDINATOR NODE
           (analyst machine)
                    |
       5. Pointwise log-likelihood summation
       6. Combined profile → pooled estimate + CI
       7. Wald ratio MR estimate
       8. Diagnostics and sensitivity analyses
       9. Output HTML report
```

---

## 2. Package Architecture

The package is organized into eight functional modules, each corresponding to a stage of the analysis pipeline. All modules must be implemented as documented, exported R functions with full roxygen2 documentation.

---

### 2.1 Module 1 — Instrument Assembly (`getMRInstruments`)

This module runs at the coordinator node. It queries the IEU OpenGWAS database via the `ieugwasr` package to retrieve GWAS summary statistics for a specified exposure trait, applies LD clumping to obtain independent instruments, and returns a clean instrument table.

**Function signature:**

```r
getMRInstruments(
  exposureTraitId,          # IEU OpenGWAS trait ID (e.g. 'ieu-a-1119' for IL-6)
  pThreshold       = 5e-8,  # Genome-wide significance threshold
  r2Threshold      = 0.001, # LD clumping r-squared threshold
  kb               = 10000, # LD clumping window in kb
  ancestryPopulation = 'EUR',
  additionalSnps   = NULL   # Optional vector of SNP IDs to force-include
)
```

**Returns:** A data frame with columns: `snp_id`, `effect_allele`, `other_allele`, `beta_ZX`, `se_ZX`, `pval_ZX`, `eaf` (effect allele frequency), `gene_region`. This instrument table is serialized to disk and distributed to all sites unchanged. Attach a timestamp and OpenGWAS version as attributes for reproducibility.

**Validation checks to implement:**

- Minimum instrument count warning if fewer than 3 SNPs returned
- F-statistic approximation using `beta_ZX` and `se_ZX` to flag potentially weak instruments
- Flag strand-ambiguous SNPs (A/T or G/C) and warn user
- Stop with informative error if zero SNPs returned

---

### 2.2 Module 2 — Site Cohort Extraction (`buildMRCohort`)

Runs locally at each OMOP site. Distributes parameterized SQL that queries standard OMOP CDM tables. Must be compatible with all major SQL dialects supported by `DatabaseConnector`: SQL Server, PostgreSQL, RedShift, Oracle, Snowflake, DuckDB.

**Function signature:**

```r
buildMRCohort(
  connectionDetails,          # DatabaseConnector connection details object
  cdmDatabaseSchema,          # Schema containing CDM tables
  cohortDatabaseSchema,       # Schema for cohort tables
  cohortTable,                # Table name for cohort
  outcomeCohortId,            # Atlas cohort ID for outcome (e.g. incident colorectal cancer)
  instrumentTable,            # Output of getMRInstruments()
  genomicLinkageSchema,       # Schema containing genomic linkage table
  genomicLinkageTable,        # Table with person_id, snp_id, genotype (0/1/2)
  indexDateOffset    = 0,     # Days offset from cohort start for index date
  washoutPeriod      = 365,   # Days of prior observation required
  excludePriorOutcome = TRUE
)
```

**Queries these OMOP tables:** `PERSON` (age, sex, race), `CONDITION_OCCURRENCE` (outcome status, exclusion criteria), `OBSERVATION_PERIOD` (eligibility, follow-up), and the genomic linkage table (SNP genotypes by `person_id`).

**Returns:** A local R data frame — never written back to the network or database.

**Genotype coding:** Genotypes must be coded as 0, 1, 2 (count of effect alleles). The function must handle allele harmonization — checking that the effect allele in the instrument table matches the coding in the genomic linkage table, and flipping if necessary (`2 - genotype`). Strand-ambiguous SNPs must be flagged. Missing genotypes must be coded as `NA`, not `0`.

---

### 2.3 Module 3 — Covariate Assembly (`buildMRCovariates`)

Leverages the existing OHDSI `FeatureExtraction` package to assemble a rich covariate matrix from OMOP. This is critical for both confounding adjustment in the outcome model and for the instrument diagnostic PheWAS.

**Function signature:**

```r
buildMRCovariates(
  connectionDetails,
  cdmDatabaseSchema,
  cohortDatabaseSchema,
  cohortTable,
  outcomeCohortId,
  covariateSettings = NULL,  # If NULL, uses default MR covariate settings
  ancestryPCsTable  = NULL,  # Table with person_id and PC1-PC10
  numAncestryPCs    = 10
)
```

**Default covariate settings must include:**

- Conditions in 365-day lookback window (binary)
- Drug exposures in 365-day lookback window (binary)
- Measurements — most recent value
- Demographics (age group, sex, index year)

Ancestry principal components (PC1 through PC10) must be merged from `ancestryPCsTable` if provided — these are mandatory for controlling population stratification and must never be regularized in the outcome model.

The covariate matrix serves double duty: (1) adjustment covariates in the Cyclops outcome model, and (2) the full phenome for the instrument PheWAS diagnostic. Document this dual use explicitly in roxygen2.

---

### 2.4 Module 4 — Instrument Diagnostics (`runInstrumentDiagnostics`)

Performs instrument validation using OMOP covariate data. This is a major differentiating feature of this package. Produces a structured diagnostic report object that feeds into the final HTML report.

**Function signature:**

```r
runInstrumentDiagnostics(
  cohortData,                   # Output of buildMRCohort()
  covariateData,                # Output of buildMRCovariates()
  instrumentTable,              # Output of getMRInstruments()
  exposureProxyConceptIds = NULL,   # OMOP concept IDs for exposure measurement proxy
  negativeControlOutcomeIds = NULL, # Cohort IDs for negative control outcomes
  pValueThreshold = 0.05 / nrow(covariateData)  # Bonferroni correction default
)
```

**Diagnostics to implement:**

- **First stage F-statistic:** Regress exposure proxy measurement on each SNP. Report F-statistic per SNP. Flag if F < 10. If no proxy measurement is available, compute approximate F from `beta_ZX` and `se_ZX` from GWAS.
- **Instrument PheWAS:** Regress each SNP against every covariate in the covariate matrix using logistic (binary) or linear (continuous) regression. Apply Bonferroni correction. Flag associations surviving correction. Output a Manhattan-style plot with covariates on x-axis ordered by domain, -log10(p) on y-axis.
- **Negative control outcome test:** Test SNP association with each negative control outcome. Significant associations suggest pleiotropy. Plot results.
- **Allele frequency check:** Compare effect allele frequency in the cohort to the GWAS reference. Flag discrepancy > 0.1 as possible strand flip or population mismatch.
- **Missing genotype summary:** Report missingness per SNP. Flag SNPs with > 10% missing genotypes.

**Returns:** A named list:

```r
list(
  fStatistics           = ...,  # Data frame: snp_id, f_statistic, flag
  phewasResults         = ...,  # Data frame: snp_id, covariate_id, beta, pval, flag
  negativeControlResults = ..., # Data frame: snp_id, outcome_id, beta, pval
  afComparison          = ...,  # Data frame: snp_id, eaf_cohort, eaf_gwas, discrepancy
  missingnessReport     = ...,  # Data frame: snp_id, n_missing, pct_missing
  diagnosticFlags       = ...   # Named logical vector: which checks failed
)
```

---

### 2.5 Module 5 — Outcome Model and Grid Likelihood (`fitOutcomeModel`)

The methodological core. Fits a regularized logistic regression of the outcome on SNP genotype(s) plus covariates using the `Cyclops` package, then evaluates the profile log-likelihood across a pre-specified grid of `beta_ZY` values.

**Function signature:**

```r
fitOutcomeModel(
  cohortData,
  covariateData,
  instrumentTable,
  betaGrid              = seq(-3, 3, by = 0.01),  # 600 grid points default
  regularizationVariance = 0.1,    # Cyclops prior variance for covariates
  instrumentRegularization = FALSE, # Do NOT regularize SNP coefficients
  outcomeType           = 'binary', # 'binary' or 'survival'
  analysisType          = 'alleleScore' # 'alleleScore' or 'perSNP'
)
```

**Analysis types:**

- `alleleScore`: Fits a single model using a weighted allele score (weights = `beta_ZX / se_ZX^2`, normalized) as the exposure variable. Primary mode.
- `perSNP`: Fits separate models for each SNP and returns a profile per SNP for use in multi-SNP sensitivity analyses.

**Grid likelihood evaluation:** After fitting the unconstrained Cyclops model, the profile log-likelihood at each grid point is evaluated by fixing `beta_ZY` to that value and optimizing over all other parameters. Use warm starting from adjacent grid points — do not refit from scratch at each grid point. SNP coefficients (`instrumentRegularization = FALSE`) must never be penalized regardless of the regularization setting for other covariates.

**Return object (shared with coordinator):**

```r
list(
  siteId        = 'site_A',
  betaGrid      = betaGrid,
  logLikProfile = logLikProfile,  # numeric vector, same length as betaGrid
  nCases        = nCases,
  nControls     = nControls,
  snpIds        = instrumentTable$snp_id,
  diagnosticFlags = diagnosticFlags
)
```

This object contains no individual-level data. It must be serializable to `.rds` format for secure file transfer between sites and coordinator.

---

### 2.6 Module 6 — Coordinator Pooling (`poolLikelihoodProfiles` + `computeMREstimate`)

Runs at the analyst's coordinator machine after receiving profile objects from all sites.

**Pooling function:**

```r
poolLikelihoodProfiles(
  siteProfileList,            # Named list of site profile objects
  validateGridAlignment = TRUE # Check all sites used same grid
)
```

Validate that all sites used identical `betaGrid` vectors. If not, interpolate to a common grid using spline interpolation and warn the user. The pooling operation is a pointwise sum of log-likelihood vectors:

```r
log_lik_combined <- Reduce("+", lapply(siteProfileList, function(x) x$logLikProfile))
```

**MR estimation function:**

```r
computeMREstimate(
  combinedProfile,   # Output of poolLikelihoodProfiles()
  instrumentTable,   # For beta_ZX and se_ZX
  ciLevel = 0.95
)
```

**Implementation steps:**

1. Find peak of combined log-likelihood → `beta_ZY_hat`
2. Compute likelihood-based CI: region where log-likelihood does not drop below `peak - qchisq(ciLevel, df=1)/2` (= 1.92 for 95% CI)
3. Apply Wald ratio: `beta_MR = beta_ZY_hat / beta_ZX`
4. Propagate uncertainty via delta method: `se_MR = sqrt((se_ZY/beta_ZX)^2 + (beta_ZY_hat * se_ZX / beta_ZX^2)^2)`
5. Compute two-sided p-value
6. Warn if MLE is at a grid boundary

With multiple SNPs, use inverse-variance weighted combination across per-SNP estimates.

---

### 2.7 Module 7 — Sensitivity Analyses (`runSensitivityAnalyses`)

Implements standard MR sensitivity analyses. Requires per-SNP `beta_ZY` estimates, available when `analysisType = 'perSNP'` in Module 5.

**Function signature:**

```r
runSensitivityAnalyses(
  perSnpEstimates,   # Data frame: snp_id, beta_ZY, se_ZY, beta_ZX, se_ZX
  methods = c('IVW', 'MREgger', 'WeightedMedian', 'Steiger', 'LeaveOneOut')
)
```

**Methods to implement:**

| Method | Description | Minimum SNPs |
|--------|-------------|--------------|
| **IVW** | Weighted regression of `beta_ZY` on `beta_ZX` through origin, weights = `1/se_ZY^2`. Primary estimate. | 1 |
| **MR-Egger** | Weighted regression with intercept. Non-zero intercept = directional pleiotropy. Use `TwoSampleMR::mr_egger_regression()`. | 3 |
| **Weighted Median** | Bootstrap weighted median of per-SNP ratio estimates. Valid if <50% of instruments invalid. Use `TwoSampleMR::mr_weighted_median()`. | 3 |
| **Steiger Filtering** | Test whether each SNP explains more variance in exposure than outcome. Remove failing SNPs, re-run IVW. | 2 |
| **Leave-One-Out** | Drop each SNP in turn, recompute IVW. Returns K estimates for forest plot. | 3 |

**Returns:** Named list of result data frames, one per method. Each data frame contains: `method`, `beta_MR`, `se_MR`, `ci_lower`, `ci_upper`, `pval`, plus method-specific columns (e.g., `egger_intercept`, `egger_intercept_pval`).

---

### 2.8 Module 8 — Reporting (`generateMRReport`)

Produces a self-contained HTML report using R Markdown. Must be a single HTML file with no external dependencies, suitable for sharing with non-technical stakeholders.

**Function signature:**

```r
generateMRReport(
  mrEstimate,          # Output of computeMREstimate()
  sensitivityResults,  # Output of runSensitivityAnalyses()
  diagnosticResults,   # Output of runInstrumentDiagnostics()
  combinedProfile,     # For likelihood profile plot
  siteProfileList,     # For per-site profile visualization
  exposureLabel = 'Exposure',
  outcomeLabel  = 'Outcome',
  outputPath    = './Medusa_report.html'
)
```

**Report sections (all required):**

1. **Executive Summary** — Plain-language interpretation. Evidence strength rating: Strong / Moderate / Weak / Inconclusive, based on sensitivity analysis concordance.
2. **Instrument Summary Table** — All SNPs: `snp_id`, `beta_ZX`, `se_ZX`, p-value, gene region, F-statistic per site.
3. **Likelihood Profile Plot** — Individual site curves in light blue, combined curve in dark blue, MLE marked, CI shaded. Built with ggplot2.
4. **Main MR Result** — Odds ratio (or beta for continuous outcome), 95% CI, p-value.
5. **Sensitivity Analysis Plot** — Horizontal forest plot: IVW, MR-Egger, Weighted Median, Steiger, color-coded by agreement with main estimate.
6. **PheWAS Plot** — Manhattan-style plot of SNP associations across OMOP covariate domains. Red horizontal line at Bonferroni threshold.
7. **Negative Control Plot** — Point estimates for SNP–negative control outcome associations, should cluster near null.
8. **Site Contribution Table** — N cases, N controls, F-statistic, diagnostic flags per site.
9. **Auto-generated Methods Text** — Manuscript-ready methods paragraph populated with actual analysis parameters.

---

## 3. Package File Structure

Implement the following directory and file structure exactly:

```
Medusa/
├── DESCRIPTION
├── NAMESPACE
├── LICENSE
├── README.md
├── R/
│   ├── getMRInstruments.R
│   ├── buildMRCohort.R
│   ├── buildMRCovariates.R
│   ├── runInstrumentDiagnostics.R
│   ├── fitOutcomeModel.R
│   ├── poolLikelihoodProfiles.R
│   ├── computeMREstimate.R
│   ├── runSensitivityAnalyses.R
│   ├── generateMRReport.R
│   ├── helpers.R               # Internal utility functions, mrTheme()
│   ├── harmonizeAlleles.R      # Allele harmonization utilities
│   └── Medusa-package.R       # Package-level documentation
├── inst/
│   ├── sql/
│   │   ├── extractOutcomeCohort.sql
│   │   ├── extractGenotypes.sql
│   │   └── extractNegativeControls.sql
│   ├── rmd/
│   │   └── MRReport.Rmd        # Report template
│   └── extdata/
│       └── defaultNegativeControlOutcomes.csv
├── tests/
│   ├── testthat/
│   │   ├── helper-simulateData.R
│   │   ├── test-getMRInstruments.R
│   │   ├── test-buildMRCohort.R
│   │   ├── test-fitOutcomeModel.R
│   │   ├── test-poolLikelihoodProfiles.R
│   │   ├── test-computeMREstimate.R
│   │   ├── test-runSensitivityAnalyses.R
│   │   ├── test-harmonizeAlleles.R
│   │   └── test-scientificValidation.R
│   └── testthat.R
└── vignettes/
    ├── GettingStarted.Rmd
    ├── ColorectalCancerIL6Example.Rmd
    └── FederatedAnalysisGuide.Rmd
```

---

## 4. Package Dependencies

| Package | Source | Purpose | Type |
|---------|--------|---------|------|
| Cyclops | OHDSI GitHub | Regularized regression + likelihood evaluation | Imports |
| FeatureExtraction | OHDSI GitHub | OMOP covariate assembly | Imports |
| DatabaseConnector | OHDSI GitHub | Multi-dialect SQL execution | Imports |
| SqlRender | OHDSI GitHub | SQL parameterization and dialect translation | Imports |
| ieugwasr | CRAN | OpenGWAS API queries | Imports |
| TwoSampleMR | GitHub (mrcieu) | MR-Egger, weighted median sensitivity analyses | Imports |
| ggplot2 | CRAN | All visualizations | Imports |
| dplyr | CRAN | Data manipulation | Imports |
| rmarkdown | CRAN | HTML report generation | Imports |
| knitr | CRAN | Report rendering | Imports |
| patchwork | CRAN | Plot composition | Imports |
| survival | CRAN | Survival outcome support | Suggests |
| testthat | CRAN | Unit testing (>= 3.0.0) | Suggests |
| withr | CRAN | Test isolation | Suggests |
| mockery | CRAN | Mocking external calls in tests | Suggests |

OHDSI packages are not on CRAN. The `DESCRIPTION` file must include:

```
Remotes:
  ohdsi/Cyclops,
  ohdsi/FeatureExtraction,
  ohdsi/DatabaseConnector,
  ohdsi/SqlRender,
  mrcieu/TwoSampleMR
```

---

## 5. DESCRIPTION File

```
Package: Medusa
Title: Federated Mendelian Randomization on OMOP Common Data Model
Version: 0.1.0
Authors@R: person('First', 'Last', email = 'author@institution.edu',
           role = c('aut', 'cre'))
Description: Implements two-sample Mendelian Randomization within the OHDSI
  ecosystem using the OMOP Common Data Model. Supports federated analysis
  across distributed health networks via one-shot profile likelihood pooling.
  Includes instrument validation diagnostics leveraging OMOP covariate richness,
  integration with public GWAS databases via OpenGWAS, and comprehensive
  sensitivity analyses. Designed for drug target validation in
  pharmacoepidemiology and oncology research.
License: Apache License 2.0
Encoding: UTF-8
LazyData: true
RoxygenNote: 7.3.0
Suggests:
    testthat (>= 3.0.0),
    withr,
    mockery,
    survival
Config/testthat/edition: 3
Remotes:
  ohdsi/Cyclops,
  ohdsi/FeatureExtraction,
  ohdsi/DatabaseConnector,
  ohdsi/SqlRender,
  mrcieu/TwoSampleMR
```

---

## 6. Testing Requirements

All tests use the `testthat` framework (edition 3). A simulation helper file must be created at `tests/testthat/helper-simulateData.R` that generates synthetic data for testing without requiring a live database connection or GWAS API access.

### 6.1 Simulation Helper (`helper-simulateData.R`)

Implement the following functions, available to all test files:

**`simulateMRData(n, nSnps, trueEffect, confoundingStrength)`**  
Generates a data frame with columns `person_id`, `outcome`, `snp_1` through `snp_K` (coded 0/1/2), plus confounders. The true causal effect parameter must be recoverable from this data. Used to validate that the Wald ratio estimator is approximately unbiased.

**`simulateInstrumentTable(nSnps)`**  
Returns a data frame mimicking `getMRInstruments()` output with plausible `beta_ZX` and `se_ZX` values.

**`simulateSiteProfiles(nSites, betaGrid, trueBeta)`**  
Generates a list of site profile objects as if returned by `fitOutcomeModel()` at multiple sites, with log-likelihood profiles consistent with a known true `beta_ZY` value. Used to validate pooling.

**`simulateCovariateData(n, nCovariates)`**  
Returns a sparse matrix mimicking `FeatureExtraction` output for PheWAS testing.

---

### 6.2 Test File Specifications

#### `test-getMRInstruments.R`

- Mock `ieugwasr::associations()` and `ieugwasr::ld_clump()` — never make live API calls in tests
- Test that function returns expected column names
- Test that strand-ambiguous SNPs (A/T and G/C) are flagged correctly
- Test that function stops with informative error when zero SNPs are returned
- Test that F-statistic approximation is computed correctly from `beta_ZX` and `se_ZX`

#### `test-buildMRCohort.R`

- Mock `DatabaseConnector::connect()` and `DatabaseConnector::querySql()`
- Test allele harmonization: if effect allele in instrument table doesn't match genomic table, genotypes must be flipped (`2 - genotype`) and beta direction inverted
- Test that missing genotypes are coded as `NA`, not `0`
- Test that washout period filtering removes persons with insufficient prior observation
- Test that persons with prior outcome are excluded when `excludePriorOutcome = TRUE`

#### `test-fitOutcomeModel.R`

- Use `simulateMRData()` to generate data with a known true effect
- Test that the log-likelihood profile is concave (unique maximum)
- Test that the MLE from the profile is within 0.1 of the true `beta_ZY` for n = 10,000
- Test that the profile object contains all required named elements
- Test that `instrumentRegularization = FALSE` means SNP coefficient is not penalized
- Test that both `alleleScore` and `perSNP` modes return valid profile objects

#### `test-poolLikelihoodProfiles.R`

- Use `simulateSiteProfiles()` to generate profiles from 3 simulated sites
- Test that pointwise sum equals manual sum of vectors
- Test that pooled MLE recovers true `beta` within tolerance for large simulated N
- Test that misaligned grids trigger a warning and spline interpolation
- Test that single-site pooling returns the same result as the un-pooled profile

#### `test-computeMREstimate.R`

- Test Wald ratio formula: `beta_MR = beta_ZY / beta_ZX`
- Test delta method standard error formula
- Test that 95% CI from likelihood profile contains true parameter in simulation
- Test that `ciLevel` parameter correctly changes the chi-squared threshold
- Test that a profile with multiple local maxima raises a warning

#### `test-runSensitivityAnalyses.R`

- Test IVW estimate matches manual weighted regression calculation
- Test MR-Egger returns intercept estimate and its p-value
- Test leave-one-out returns K estimates for K SNPs
- Test Steiger filtering removes SNPs with wrong variance direction
- Test that all methods handle the single-SNP case gracefully (some require >= 3 SNPs — must skip with warning, not error)

#### `test-harmonizeAlleles.R`

- Test that an A/C SNP is correctly harmonized (no ambiguity)
- Test that an A/T SNP is flagged as strand-ambiguous
- Test that a palindromic SNP with EAF near 0.5 is dropped with a warning
- Test that an allele flip correctly inverts effect direction (`beta * -1`) and EAF (`1 - EAF`)

#### `test-scientificValidation.R`

See [Section 12](#12-scientific-validation). All tests in this file may use `skip_on_cran()` as they are computationally intensive.

---

## 7. Vignette Requirements

Three vignettes are required. All must be fully executable end-to-end using only synthetic data or public data — no proprietary database connections should be required to build any vignette. Use the simulation helpers or publicly available cached data for all demonstrations.

### 7.1 `GettingStarted.Rmd`

Target audience: epidemiologist new to the package.

- Installation instructions including OHDSI package remotes
- Conceptual overview of two-sample MR in plain language
- Explanation of the federated profile likelihood approach and why it preserves privacy
- Quick start example using fully simulated data (no database connection)
- Description of the output report and how to interpret each section
- FAQ section covering:
  - What genomic data linkage is required at each site?
  - What if my site has no genomic data?
  - How many sites are needed for adequate power?
  - What OMOP CDM version is supported?

### 7.2 `ColorectalCancerIL6Example.Rmd`

Target audience: researcher wanting to run a real oncology analysis. This is the primary scientific vignette.

- Scientific motivation: IL-6 signaling and colorectal cancer risk
- Step-by-step walkthrough of `getMRInstruments()` using real IL6R region SNPs from OpenGWAS — wrap the API call in `tryCatch` with cached fallback data so vignette builds offline
- Cohort definition guidance: how to define incident colorectal cancer in OMOP, SNOMED/ICD codes, exclusion criteria
- Running `buildMRCovariates()` and interpreting the covariate summary
- Running `runInstrumentDiagnostics()` and interpreting the PheWAS plot — show an example with a problematic SNP and explain how to handle it
- Running `fitOutcomeModel()` and explaining the profile likelihood output
- Demonstrating `poolLikelihoodProfiles()` with three simulated site profiles
- Running all sensitivity analyses and interpreting concordance
- Generating and walking through the HTML report
- Interpretation section: translating the result into drug development language (target validation evidence strength)

### 7.3 `FederatedAnalysisGuide.Rmd`

Target audience: OHDSI network coordinator or data engineer.

- Detailed explanation of what data leaves each site (only the log-likelihood profile vector) and what does not (individual genotypes, outcomes, covariates)
- Required OMOP CDM version and genomic linkage table schema specification
- Template site analysis script, parameterized for site-specific connection details
- Instructions for securely transferring `.rds` profile files from sites to coordinator
- How to handle sites with different OMOP versions or incomplete genomic coverage
- Troubleshooting guide: common errors, diagnostic flag meanings, handling weak instruments

---

## 8. Documentation Requirements

Every exported function must have complete roxygen2 documentation including:

- `@title` — one-line title
- `@description` — paragraph describing what the function does and when to use it
- `@param` — every parameter documented with type, allowed values, and explanation of default
- `@return` — complete description of return object structure including all named elements
- `@details` — technical details of the algorithm (especially for Modules 5 and 6)
- `@references` — relevant statistical literature in APA format
- `@examples` — at least one runnable example using simulated data, executable without database or API access
- `@seealso` — links to related functions in the package
- `@export` — all user-facing functions must be exported

**Package-level documentation** (`Medusa-package.R`): Describe the overall package, its scientific purpose, the federated architecture, and cite key MR methods literature:

- Davey Smith & Hemani (2014) — foundational MR reference
- Bowden et al. (2015) — MR-Egger
- Bowden et al. (2016) — Weighted median estimator
- Burgess & Thompson (2015) — Mendelian randomization methods textbook

**`README.md`** must include:

- Installation instructions
- Minimal working example using simulated data
- ASCII diagram of the data flow architecture
- Table of main functions and their purpose
- Links to vignettes on pkgdown site

---

## 9. Coding Standards and Style

- Follow the OHDSI R package style guide: `snake_case` for function names, `camelCase` for variable names within functions, `SCREAMING_SNAKE_CASE` for constants
- No hard-coded database schema names or table names — always accept as function parameters
- All SQL must go through `SqlRender::renderSql()` and `SqlRender::translateSql()` — never write dialect-specific SQL directly
- All database queries must use `DatabaseConnector` — never use `DBI` or `RODBC` directly
- Wrap all external API calls (`ieugwasr`) in `tryCatch()` with informative error messages
- Use `message()` not `print()` for informative output — users can suppress with `suppressMessages()`
- Use `checkmate` or equivalent for input validation at the start of each function — validate types, ranges, and required columns before any computation
- Define `mrTheme()` in `helpers.R` and apply it to all ggplot2 plots for consistent styling
- Never use `<<-` for assignment — all state must be passed explicitly through function arguments and return values
- The package must pass `R CMD CHECK` with **zero errors and zero warnings**
- Study existing OHDSI packages (`CohortMethod`, `SelfControlledCaseSeries`) and follow their structural and naming conventions — OHDSI reviewers will expect familiar patterns

---

## 10. Error Handling and Edge Cases

The following situations must be handled with informative, actionable error messages:

| Situation | Required Behavior |
|-----------|-------------------|
| Zero SNPs from `getMRInstruments()` | `stop()` with: *"No genome-wide significant SNPs found for trait [id]. Consider relaxing pThreshold or checking trait ID."* |
| Fewer than 3 SNPs after LD clumping | `warning()`: *"Only [n] instruments available. Sensitivity analyses requiring >=3 SNPs will be skipped."* |
| No persons with genotype data | `stop()` with: *"No persons in cohort have genotype data in [table]. Verify genomic linkage table schema and person_id join."* |
| Fewer than 50 outcome cases | `warning()`: *"Only [n] cases in outcome cohort. Results may be unstable."* |
| Flat log-likelihood profile | `warning()`: *"Profile likelihood is flat — instrument may be too weak to estimate beta_ZY. Check F-statistic."* |
| MLE at grid boundary | `warning()`: *"MLE is at grid boundary [value]. Expand betaGrid range. Current range: [min, max]."* |
| Misaligned grids across sites | `warning()` + interpolate: *"Site [id] used different betaGrid. Interpolating to common grid."* |
| All SNPs fail Steiger filter | `warning()`: *"All SNPs failed Steiger filter. Investigate instrument validity."* |
| OpenGWAS API unavailable | Informative error suggesting use of `cachedInstrumentTable` parameter |

---

## 11. SQL Template Requirements

All SQL files in `inst/sql/` must use `SqlRender` parameterization syntax (`@paramName`) and must be compatible with: SQL Server, PostgreSQL, RedShift, Oracle, Snowflake, DuckDB. Test dialect translation with `SqlRender::translateSql()` for each target dialect.

### 11.1 `extractOutcomeCohort.sql`

Extracts: `person_id`, `outcome` (0/1), `index_date`, `age_at_index`, `gender_concept_id`, `observation_period_start`, `observation_period_end`.

Joins: `PERSON`, `CONDITION_OCCURRENCE`, `OBSERVATION_PERIOD` using standard OMOP column names.

Parameters: `@cdmDatabaseSchema`, `@cohortDatabaseSchema`, `@cohortTable`, `@outcomeCohortId`, `@washoutDays`, `@startDate`, `@endDate`.

### 11.2 `extractGenotypes.sql`

Extracts: `person_id`, `snp_id`, `genotype` (integer 0/1/2).

Joins the genomic linkage table to the cohort table on `person_id`. Handles the case where the genomic linkage table uses a different person identifier column name.

Parameters: `@genomicLinkageSchema`, `@genomicLinkageTable`, `@cohortDatabaseSchema`, `@cohortTable`, `@outcomeCohortId`, `@genomicPersonIdColumn`.

### 11.3 `extractNegativeControls.sql`

For each negative control outcome cohort ID, extracts whether each person in the main cohort has the outcome.

Returns long format: `person_id`, `outcome_cohort_id`, `has_outcome` (0/1).

Parameters: `@cdmDatabaseSchema`, `@cohortDatabaseSchema`, `@cohortTable`, `@outcomeCohortId`, `@negativeControlCohortIds` (comma-separated list, rendered by SqlRender).

---

## 12. Scientific Validation

The following validation tests must be implemented in `tests/testthat/test-scientificValidation.R`. All may use `skip_on_cran()`.

| Test | Criterion |
|------|-----------|
| **Estimator bias check** | Using `simulateMRData(n=50000, trueEffect=0.5)`, the IVW estimate must be within 0.05 of 0.5 averaged over 100 simulations with `set.seed()` |
| **Coverage check** | Using `simulateMRData(n=10000)`, the 95% likelihood-based CI must contain the true parameter in at least 93 of 100 simulations |
| **Pooling consistency** | Using `simulateSiteProfiles(nSites=5)` each with n=2000, the pooled estimate must be within 0.05 of the estimate from all 10,000 persons analyzed jointly |
| **Type I error check** | Using `simulateMRData(trueEffect=0.0)`, p-value must be < 0.05 in no more than 7 of 100 simulations (empirical type I error < 0.07) |

---

## 13. Deliverables Checklist

| Item | Required |
|------|----------|
| All 9 R function files with complete roxygen2 documentation | ✅ |
| `helpers.R` with `mrTheme()` and internal utilities | ✅ |
| `harmonizeAlleles.R` with allele harmonization functions | ✅ |
| `Medusa-package.R` with package-level documentation | ✅ |
| 3 SQL template files in `inst/sql/` | ✅ |
| HTML report R Markdown template in `inst/rmd/` | ✅ |
| Default negative control outcomes CSV in `inst/extdata/` | ✅ |
| Simulation helper `helper-simulateData.R` | ✅ |
| 7 unit test files | ✅ |
| Scientific validation test file | ✅ |
| 3 vignettes, fully executable without live database | ✅ |
| `DESCRIPTION` with correct dependencies and `Remotes` | ✅ |
| `README.md` with installation, quickstart, architecture diagram | ✅ |
| `R CMD CHECK` passes with 0 errors, 0 warnings | ✅ |
| Scientific validation tests passing | ✅ |

---

## 14. Implementation Notes for the AI Coder

### Start with the simulation helpers

Everything else depends on being able to generate synthetic test data. Build `helper-simulateData.R` first, then build and test each module in pipeline order (1 through 8). Do not move to the next module until the current module's tests pass.

### Cyclops profile likelihood is the hardest part

The grid likelihood evaluation in `fitOutcomeModel()` is the most technically challenging component. Study the Cyclops documentation carefully — specifically `Cyclops::fitCyclopsModel()` with the `control` argument and offset/constraint capabilities.

If profile likelihood via constrained optimization proves intractable within Cyclops, a fallback is the quadratic approximation:

```
logLik(beta) ≈ logLik_max - 0.5 * (beta - beta_hat)^2 * observed_information
```

The true profile is preferred and must be attempted first. If falling back to the quadratic approximation, document this clearly in the function's `@details` section.

### Mock all external dependencies in tests

Use `mockery::stub()` or `testthat::local_mocked_bindings()` to mock `ieugwasr` functions and `DatabaseConnector` functions. Tests must never make real network calls or database connections. API call tests must work offline.

### SQL dialect testing

Use `SqlRender::translateSql()` with `targetDialect` set to each supported dialect in your SQL tests to verify parameterization works correctly. You do not need a live database — just verify the SQL renders without syntax errors.

### The vignettes are the scientific face of the package

Invest in making `ColorectalCancerIL6Example.Rmd` read well scientifically, not just technically. A pharmacoepidemiology scientist at a pharma company should be able to read it and understand both the method and its value for drug target validation without needing to read the function documentation.

### Follow OHDSI conventions closely

Study the source code of `CohortMethod` and `SelfControlledCaseSeries` on GitHub before writing any R code. OHDSI community reviewers will expect familiar patterns for database connection handling, cohort table management, covariate data objects, and result object structure. Diverging from these conventions without good reason will slow down community adoption.

### Version the instrument table

Attach a timestamp and OpenGWAS version as attributes to the instrument table returned by `getMRInstruments()`. GWAS results change as studies update — reproducibility requires knowing exactly which data was used at analysis time.

---

*End of Medusa Package Specification v1.0*
