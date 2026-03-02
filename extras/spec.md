
**`Medusa`** — keeps it simple and discoverable in the OHDSI ecosystem.

---

## Overall Data Flow

```
PUBLIC GWAS DATABASE          EACH OMOP SITE
(OpenGWAS API)                (runs locally)
      |                              |
      |                              |
  β_ZX, SE_ZX              1. Cohort extraction
  SNP-exposure              2. Covariate assembly
  association               3. Cyclops outcome model
      |                     4. Grid likelihood profile
      |                              |
      └──────────┬───────────────────┘
                 |
         COORDINATOR NODE
         (analyst machine)
                 |
    5. Pointwise log-likelihood summation
    6. Combined profile → pooled estimate + CI
    7. Wald ratio MR estimate
    8. Diagnostics and sensitivity analyses
    9. Output report
```

---

## Module 1: Exposure Side — GWAS Instrument Assembly

This runs at the coordinator, not at sites. It queries public databases to assemble the instrument.

**Function: `getMRInstruments()`**

Inputs:
- Exposure trait (e.g. "IL-6" or an IEU OpenGWAS trait ID)
- P-value threshold for SNP inclusion (default 5e-8, genome-wide significance)
- R-squared threshold for LD clumping (to ensure SNPs are independent of each other)
- Ancestry of reference panel for LD clumping

What it does:
- Queries OpenGWAS API for SNPs associated with the exposure trait
- Performs LD clumping to get a set of independent instruments
- Returns a clean instrument table

Output table — one row per SNP:

| snp_id | effect_allele | other_allele | beta_ZX | se_ZX | pval | eaf |
|--------|--------------|--------------|---------|-------|------|-----|
| rs2228145 | A | C | 0.34 | 0.02 | 3e-12 | 0.41 |
| rs4537545 | T | G | 0.21 | 0.03 | 8e-9 | 0.28 |

This instrument table gets distributed to all sites. It tells each site exactly which SNPs to look for in their genomic data and which allele coding to use.

---

## Module 2: Site-Level Cohort Extraction

This runs locally at each OMOP site. The package distributes parameterized SQL that runs against standard OMOP tables.

**Function: `buildMRCohort()`**

Inputs:
- Target cohort definition (outcome — e.g. incident colorectal cancer) as an ATLAS cohort ID or parameterized SQL
- Instrument table from Module 1
- Covariate settings object (see below)
- Genomic linkage table name and schema (site-specific configuration)

What it does against OMOP tables:

**From PERSON:** age, sex, race/ethnicity

**From CONDITION_OCCURRENCE:** outcome status (cancer yes/no), exclusion criteria (prior cancer), negative control outcomes

**From DRUG_EXPOSURE:** medication covariates for confounder adjustment and PheWAS diagnostics

**From MEASUREMENT:** relevant lab values (e.g. CRP as IL-6 proxy for instrument validation)

**From OBSERVATION_PERIOD:** follow-up time, index date, eligibility windows

**From genomic linkage table:** genotype at each instrument SNP, coded as 0/1/2 (number of effect alleles), matched to person_id

Output: a local analysis-ready dataset — never leaves the site.

---

## Module 3: Covariate Assembly via FeatureExtraction

Rather than writing custom covariate SQL, leverage the existing OHDSI `FeatureExtraction` package which already knows how to pull thousands of covariates from OMOP in a standardized way.

**Function: `buildMRCovariates()`**

Covariate groups to include:

- **Ancestry principal components** — pulled from genomic linkage table, essential for population stratification control
- **Demographics** — age, sex, index year
- **Condition covariates** — conditions in 365-day pre-index window (FeatureExtraction handles this)
- **Drug covariates** — drug exposures in pre-index window
- **Measurement covariates** — recent lab values

For the outcome model, Cyclops with LASSO regularization handles the high dimensionality — you can include thousands of covariates without overfitting because Cyclops regularizes automatically. This is the key advantage of using Cyclops here rather than standard logistic regression.

---

## Module 4: Instrument Diagnostics — PheWAS and Validation

This is where OMOP's covariate richness pays off. Runs locally at each site before the main analysis.

**Function: `runInstrumentDiagnostics()`**

**First stage check:**
Regress each SNP on the exposure proxy (e.g. CRP from MEASUREMENT table) to confirm the instrument actually predicts the exposure in this population. Reports F-statistic — convention is F > 10 indicates adequate instrument strength. Weak instruments bias MR estimates toward the confounded observational association, so this is critical.

**Instrument PheWAS:**
Regress each SNP against every available covariate from FeatureExtraction — hundreds of conditions, drugs, measurements. Under the null (clean instrument), SNPs should associate with nothing except the exposure pathway. Flag any covariate associations that survive multiple testing correction. This is the balance table equivalent for genetic instruments — analogous to checking covariate balance after propensity score matching, a concept you know well from observational work.

**Negative control outcomes:**
Test instrument association with a set of pre-specified negative control outcomes — conditions that IL-6 biologically should not cause. Significant associations suggest pleiotropy. The package would ship with a default negative control outcome set for common exposures, overridable by the user.

Output: a diagnostic report per site, flagging any instrument validity concerns before the main analysis proceeds.

---

## Module 5: Outcome Model via Cyclops — Grid Likelihood Estimation

This is the methodological core. Runs locally at each site.

**Function: `fitOutcomeModel()`**

For each SNP (or for a combined allele score), fit a logistic regression of cancer outcome on SNP genotype, adjusting for ancestry PCs and selected covariates, using Cyclops for regularized estimation.

**Then the key step — grid likelihood profile:**

Rather than just extracting the point estimate and standard error, the package evaluates the profile log-likelihood at each point on a pre-specified grid of β_ZY values.

```r
# Conceptually what happens inside the function
beta_grid <- seq(-2, 2, by = 0.01)  # 400 grid points
log_lik_profile <- numeric(length(beta_grid))

for (i in seq_along(beta_grid)) {
  # Evaluate Cyclops log-likelihood at this fixed beta value
  log_lik_profile[i] <- evaluateCyclopsLogLik(
    cyclopsModel, 
    beta_ZY = beta_grid[i]
  )
}
```

In practice Cyclops exposes likelihood evaluation functions so this doesn't require refitting the model 400 times — it evaluates the likelihood at fixed parameter values efficiently.

**What each site shares back to coordinator:**

```r
list(
  site_id = "site_A",
  beta_grid = beta_grid,          # the grid (same at all sites)
  log_lik_profile = log_lik_profile,  # 400 numbers
  n_cases = n_cases,              # for reporting
  n_controls = n_controls,
  f_statistic = f_stat,           # from first stage
  diagnostic_flags = diag_flags   # from Module 4
)
```

That's all that leaves the site. A vector of ~400 numbers plus metadata.

---

## Module 6: Coordinator — Pooling and MR Estimation

Runs at the analyst's machine after receiving profiles from all sites.

**Function: `poolLikelihoodProfiles()`**

Step 1 — sum log-likelihoods pointwise across sites:

```r
log_lik_combined <- Reduce("+", lapply(site_results, function(x) x$log_lik_profile))
```

Step 2 — find the peak (maximum likelihood estimate of β_ZY):

```r
beta_ZY_hat <- beta_grid[which.max(log_lik_combined)]
```

Step 3 — likelihood-based confidence interval. Find the region where the log-likelihood doesn't drop more than 1.92 below the peak (this gives the 95% CI without assuming normality):

```r
ci_threshold <- max(log_lik_combined) - 1.92
ci_region <- beta_grid[log_lik_combined >= ci_threshold]
ci_lower <- min(ci_region)
ci_upper <- max(ci_region)
```

**Function: `computeMREstimate()`**

Apply the Wald ratio using the pooled β_ZY and the β_ZX from GWAS:

```r
beta_MR <- beta_ZY_hat / beta_ZX
se_MR <- sqrt((se_ZY / beta_ZX)^2 + (beta_ZY_hat * se_ZX / beta_ZX^2)^2)
```

With multiple SNPs, use inverse-variance weighted combination across SNPs.

---

## Module 7: Sensitivity Analyses

**Function: `runSensitivityAnalyses()`**

All of these operate on the pooled estimate and the multi-SNP results:

**MR-Egger** — fits a weighted regression of SNP-outcome associations on SNP-exposure associations, allowing a non-zero intercept. A significant intercept indicates directional pleiotropy. The slope gives a pleiotropy-robust causal estimate.

**Weighted median** — sorts SNPs by their ratio estimates and takes the weighted median. Valid even if up to 50% of instruments are invalid. Computationally simple once you have per-SNP β_ZY estimates from each site.

**Leave-one-SNP-out** — drops each SNP in turn and recomputes the pooled estimate. Checks whether any single SNP is driving the result. Requires re-running the pooling step K times (once per SNP) but this is fast since you already have the likelihood profiles.

**Steiger filtering** — tests whether each SNP explains more variance in the exposure than the outcome, as expected if the causal direction is exposure → outcome rather than reverse. Flags SNPs that fail this test.

---

## Module 8: Output and Reporting

**Function: `generateMRReport()`**

Produces a structured output object and an HTML report containing:

**Instrument summary table** — SNPs used, effect sizes, F-statistics per site

**Likelihood profile plot** — individual site curves in light color, combined curve bold, peak and CI marked. This is the signature visualization of the package and makes the pooling transparent and intuitive.

**Main MR result** — pooled causal estimate, CI, p-value, interpretation

**Sensitivity analysis plot** — showing main IVW estimate alongside MR-Egger and weighted median estimates. If they all agree, confidence in the result is high.

**PheWAS plot** — Manhattan-style plot of instrument associations across all OMOP covariates. Clean instruments show nothing significant except the exposure pathway.

**Negative control plot** — instrument associations with negative control outcomes, should cluster around null.

**Site-level diagnostics table** — F-statistics, sample sizes, any flags per site.

---

## Package Structure Summary

```
MRohdsi/
├── R/
│   ├── getMRInstruments.R       # Module 1 - GWAS query
│   ├── buildMRCohort.R          # Module 2 - OMOP extraction
│   ├── buildMRCovariates.R      # Module 3 - FeatureExtraction
│   ├── runInstrumentDiagnostics.R  # Module 4 - PheWAS + validation
│   ├── fitOutcomeModel.R        # Module 5 - Cyclops + grid likelihood
│   ├── poolLikelihoodProfiles.R # Module 6 - Coordinator pooling
│   ├── computeMREstimate.R      # Module 6 - Wald ratio + IVW
│   ├── runSensitivityAnalyses.R # Module 7 - Egger, median, Steiger
│   └── generateMRReport.R       # Module 8 - Output
├── inst/
│   └── sql/                     # Parameterized OMOP SQL
├── vignettes/
│   └── colorectal_IL6_example.Rmd  # The example we walked through
└── DESCRIPTION
```

---

## Dependencies

| Package | Purpose |
|---------|---------|
| Cyclops | Regularized outcome regression + likelihood evaluation |
| FeatureExtraction | OMOP covariate assembly |
| DatabaseConnector | OMOP database connections |
| TwoSampleMR | Sensitivity analyses (MR-Egger, weighted median) |
| ieugwasr | OpenGWAS API queries |
| ggplot2 | Visualization |
| survival | Time-to-event outcome support |

---

## What Makes This Novel

To summarize the genuine contributions:

1. **First MR implementation native to OMOP CDM** — standardized cohort extraction and covariate assembly using existing OHDSI infrastructure
2. **One-shot federated pooling via profile likelihood** — no iterative protocol, one communication round, exact likelihood-based inference
3. **Instrument diagnostics using OMOP covariate richness** — PheWAS of instruments, negative control outcomes, automated balance checking at scale
4. **Federated design** — individual data never leaves sites, consistent with OHDSI network principles

This is a methods paper and a software paper. The IL-6 colorectal cancer example could be the applied demonstration in both.
