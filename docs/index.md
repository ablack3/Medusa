<div id="main" class="col-md-9" role="main">

# Medusa

<div class="section level1">

**Mendelian Estimation in Distributed Standardized Analytics**

Federated two-sample Mendelian Randomization on the OMOP Common Data
Model.

**PACKAGE UNDER CONSTRUCTION: DO NOT USE**

<div class="section level2">

## Overview

Medusa implements two-sample Mendelian Randomization (MR) natively
within the OHDSI ecosystem. It enables causal inference across
distributed health networks without requiring individual-level data to
leave any site.

The core innovation is **one-shot federated pooling** via profile
likelihood aggregation: each site computes a log-likelihood profile and
shares only site-level summary files centered on that numeric vector.
The coordinator sums profiles across sites to obtain a pooled estimate
with no iterative communication protocol.

</div>

<div class="section level2">

## Architecture

        ┌──────────────────────────────────────────────────────────┐
        │                    COORDINATOR NODE                      │
        │                                                          │
        │  getMRInstruments()  ──────►  Instrument Table           │
        │       │                       (from OpenGWAS)            │
        │       │ distribute                                       │
        │       ▼                                                  │
        │  ┌─────────┐  ┌─────────┐  ┌─────────┐                   │
        │  │ Site A  │  │ Site B  │  │ Site C  │   OMOP CDM        │
        │  │─────────│  │─────────│  │─────────│   Sites           │
        │  │buildMR  │  │buildMR  │  │buildMR  │                   │
        │  │Cohort() │  │Cohort() │  │Cohort() │                   │
        │  │    ▼    │  │    ▼    │  │    ▼    │                   │
        │  │fitOut   │  │fitOut   │  │fitOut   │                   │
        │  │come     │  │come     │  │come     │                   │
        │  │Model()  │  │Model()  │  │Model()  │                   │
        │  │    │    │  │    │    │  │    │    │                   │
        │  └────┼────┘  └────┼────┘  └────┼────┘                   │
        │       │            │            │                        │
        │       ▼            ▼            ▼                        │
        │    Profile       Profile      Profile   ◄── Only these   │
        │    vectors       vectors      vectors       leave sites  │
        │       │            │            │                        │
        │       └────────────┼────────────┘                        │
        │                    ▼                                     │
        │         poolLikelihoodProfiles()                         │
        │                    │                                     │
        │                    ▼                                     │
        │           computeMREstimate()                            │
        │           runSensitivityAnalyses()                       │
        │           generateMRReport()                             │
        └──────────────────────────────────────────────────────────┘

</div>

<div class="section level2">

## Installation

<div id="cb2" class="sourceCode">

``` r
# Install from GitHub
remotes::install_github("OHDSI/Medusa")

# Or install with all dependencies
remotes::install_github("OHDSI/Medusa", dependencies = TRUE)
```

</div>

</div>

<div class="section level2">

## Quick Start

<div id="cb3" class="sourceCode">

``` r
library(Medusa)

# 1. Assemble instruments from OpenGWAS (coordinator)
instruments <- getMRInstruments(
  exposureTraitId = "ieu-a-1119",  # IL-6 receptor
  pThreshold = 5e-8,
  r2Threshold = 0.001
)

# 2. At each site: build cohort and fit outcome model
cohort <- buildMRCohort(
  connectionDetails = connectionDetails,
  cdmDatabaseSchema = "cdm",
  cohortDatabaseSchema = "results",
  cohortTable = "cohort",
  outcomeCohortId = 1234,
  instrumentTable = instruments,
  genomicDatabaseSchema = "genomics"  # Schema with VARIANT_OCCURRENCE
)

profile <- fitOutcomeModel(
  cohortData = cohort,
  instrumentTable = instruments,
  siteId = "my_site"
)

# 3. At coordinator: pool profiles and estimate
combined <- poolLikelihoodProfiles(list(siteA = profileA, siteB = profileB))
estimate <- computeMREstimate(combined, instruments)

# 4. Generate report
generateMRReport(
  mrEstimate = estimate,
  combinedProfile = combined,
  exposureLabel = "IL-6 signaling",
  outcomeLabel = "Colorectal cancer"
)
```

</div>

</div>

<div class="section level2">

## Main Functions

| Function                     | Module              | Description                                              |
|------------------------------|---------------------|----------------------------------------------------------|
| `getMRInstruments()`         | Instrument Assembly | Query OpenGWAS for instruments, LD clump                 |
| `createInstrumentTable()`    | Instrument Assembly | Build instruments from local data                        |
| `buildMRCohort()`            | Cohort Extraction   | Extract cohort + genotypes from OMOP CDM                 |
| `buildMRCovariates()`        | Covariate Assembly  | Assemble covariates via FeatureExtraction                |
| `runInstrumentDiagnostics()` | Diagnostics         | F-stats, PheWAS, allele-frequency and missingness checks |
| `fitOutcomeModel()`          | Outcome Model       | Fit outcome model, evaluate profile likelihood           |
| `poolLikelihoodProfiles()`   | Pooling             | Sum log-likelihood profiles across sites                 |
| `computeMREstimate()`        | Estimation          | Wald ratio with delta method SE                          |
| `runSensitivityAnalyses()`   | Sensitivity         | IVW, MR-Egger, weighted median, etc.                     |
| `generateMRReport()`         | Reporting           | Self-contained HTML report                               |

</div>

<div class="section level2">

## Vignettes

-   **Getting Started** — Installation, concepts, quick example with
    simulated data
-   **IL-6 and Colorectal Cancer** — Complete scientific walkthrough
-   **Federated Analysis Guide** — Network coordinator instructions

</div>

<div class="section level2">

## Requirements

-   R &gt;= 4.1.0
-   OHDSI packages: DatabaseConnector, SqlRender, Cyclops,
    FeatureExtraction
-   OMOP CDM database with [OMOP Genomic
    CDM](https://github.com/OHDSI/Genomic-CDM) (VARIANT\_OCCURRENCE
    table)
-   For instrument retrieval: internet access to OpenGWAS API

</div>

<div class="section level2">

## License

Apache License 2.0

</div>

<div class="section level2">

## References

-   Davey Smith & Hemani (2014). Mendelian randomization: genetic
    anchors for causal inference. *Human Molecular Genetics*.
-   Bowden et al. (2015). Mendelian randomization with invalid
    instruments: MR-Egger. *IJE*.
-   Bowden et al. (2016). Weighted median estimator. *Genetic
    Epidemiology*.
-   Suchard et al. (2013). Cyclops: massive parallelization of serial
    inference. *ACM TOMACS*.

</div>

</div>

</div>
