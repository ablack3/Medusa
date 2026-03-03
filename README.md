# Medusa <img src="man/figures/logo.png" align="right" height="80" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/ablack3/Medusa/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ablack3/Medusa/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/ablack3/Medusa/graph/badge.svg)](https://app.codecov.io/gh/ablack3/Medusa)
<!-- badges: end -->
  
**Mendelian Estimation in Distributed Standardized Analytics**

Federated two-sample Mendelian Randomization on the OMOP Common Data Model.

## Overview

Medusa implements two-sample Mendelian Randomization (MR) natively within the OHDSI ecosystem. It enables causal inference across distributed health networks without requiring individual-level data to leave any site.

The core innovation is **one-shot federated pooling** via profile likelihood aggregation: each site computes a log-likelihood profile and shares only that numeric vector. The coordinator sums profiles across sites to obtain a pooled estimate — no iterative communication protocol needed.

## Architecture

```
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
```

## Installation

```r
# Install from GitHub
remotes::install_github("OHDSI/Medusa")

# Or install with all dependencies
remotes::install_github("OHDSI/Medusa", dependencies = TRUE)
```

## Quick Start

```r
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
  genomicLinkageSchema = "genomics",
  genomicLinkageTable = "genotype_data"
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

## Main Functions

| Function | Module | Description |
|----------|--------|-------------|
| `getMRInstruments()` | Instrument Assembly | Query OpenGWAS for instruments, LD clump |
| `createInstrumentTable()` | Instrument Assembly | Build instruments from local data |
| `buildMRCohort()` | Cohort Extraction | Extract cohort + genotypes from OMOP CDM |
| `buildMRCovariates()` | Covariate Assembly | Assemble covariates via FeatureExtraction |
| `runInstrumentDiagnostics()` | Diagnostics | F-stats, PheWAS, negative controls |
| `fitOutcomeModel()` | Outcome Model | Fit model, evaluate profile likelihood |
| `poolLikelihoodProfiles()` | Pooling | Sum log-likelihood profiles across sites |
| `computeMREstimate()` | Estimation | Wald ratio with delta method SE |
| `runSensitivityAnalyses()` | Sensitivity | IVW, MR-Egger, weighted median, etc. |
| `generateMRReport()` | Reporting | Self-contained HTML report |

## Vignettes

- **Getting Started** — Installation, concepts, quick example with simulated data
- **IL-6 and Colorectal Cancer** — Complete scientific walkthrough
- **Federated Analysis Guide** — Network coordinator instructions

## Requirements

- R >= 4.1.0
- OHDSI packages: DatabaseConnector, SqlRender, Cyclops, FeatureExtraction
- OMOP CDM database with genomic linkage table
- For instrument retrieval: internet access to OpenGWAS API

## License

Apache License 2.0

## References

- Davey Smith & Hemani (2014). Mendelian randomization: genetic anchors for causal inference. *Human Molecular Genetics*.
- Bowden et al. (2015). Mendelian randomization with invalid instruments: MR-Egger. *IJE*.
- Bowden et al. (2016). Weighted median estimator. *Genetic Epidemiology*.
- Suchard et al. (2013). Cyclops: massive parallelization of serial inference. *ACM TOMACS*.
