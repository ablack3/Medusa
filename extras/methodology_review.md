# Medusa Methodology Review

Date: 2026-03-03

Thorough review of statistical methodology, cross-module consistency,
and correctness across all R source files in the Medusa package.

---

## Critical Issues

### 1. `buildMRCovariates` output silently ignored by `fitOutcomeModel`

**Files**: `R/buildMRCovariates.R`, `R/fitOutcomeModel.R`

`buildMRCovariates()` returns a list with class `"medusaCovariateData"`.
`fitOutcomeModel()` -> `appendCovariatesToModelData()` checks
`is.data.frame(covariateData)`, which is FALSE for the list object, so the
entire covariate block is skipped. Covariates from the real OMOP extraction
pipeline will be **silently dropped** from the outcome model.

`runInstrumentDiagnostics()` handles this correctly with
`inherits(covariateData, "medusaCovariateData")`, so only the outcome model
is affected.

**Fix**: Add `medusaCovariateData` handling in `appendCovariatesToModelData`.

### 2. `scoreDefinition` lost in CSV round-trip

**File**: `R/profileIO.R`

`exportSiteProfile()` does not write `scoreDefinition` (the allele-score
weights, betaZX, seZX) to any CSV file. `importSiteProfile()` does not
reconstruct it. Imported profiles have `scoreDefinition = NULL`.

Consequences:
- The `identicalScoreDefinition` check in `poolLikelihoodProfiles()` silently
  passes (both NULL -> TRUE), so the safety guard that all sites used the same
  allele score is bypassed.
- `computeMREstimate()` falls back to recomputing from the `instrumentTable`
  argument, which is correct only if the instrument table exactly matches what
  was used at the sites.

**Fix**: Export `scoreDefinition` as a third CSV file (e.g.,
`{prefix}_score_{siteId}.csv`), and import/reconstruct it.

### 3. Palindromic SNPs with EAF away from 0.5 not actually resolved using EAF

**File**: `R/harmonizeAlleles.R`

The docstring says "If palindromic but EAF clearly differs from 0.5, use EAF
to infer strand orientation." However, the code does not implement this. For
palindromic SNPs (e.g., A/T), the complement alleles are identical to the
original alleles, so Cases 1-4 in the allele-matching logic cannot
distinguish a strand flip from a direct match using alleles alone. The code
would need to compare the GWAS EAF against the cohort allele frequency to
resolve the ambiguity. If the GWAS is on the opposite strand, the effect
allele could be assigned incorrectly.

**Fix**: When a palindromic SNP passes the EAF threshold (EAF far from 0.5),
compare GWAS EAF against cohort allele frequency. If they are discordant
(e.g., GWAS EAF = 0.8, cohort AF = 0.2), flip beta and EAF.

---

## Medium-Severity Issues

### 4. SNP column alignment uses positional indexing

**File**: `R/runInstrumentDiagnostics.R` (allele frequency comparison and
genotype missingness)

`snpCols` from `grep("^snp_", names(cohortData))` may not be in the same
order as rows in `instrumentTable`. The comparison
`if (i <= length(snpCols)) { geno <- cohortData[[snpCols[i]]] }` associates
SNP i in the instrument table with column i from the grep results. If the
column order differs, the wrong allele frequencies will be compared.

`fitOutcomeModel.R` handles this correctly via `alignInstrumentColumns()`,
which matches by name using `makeSnpColumnName()`.

**Fix**: Use `makeSnpColumnName(instrumentTable$snp_id[i])` to look up the
correct column instead of positional indexing.

### 5. `perSnpEstimates` not exported in CSV

**File**: `R/profileIO.R`

When `analysisType = "perSNP"`, the profile object contains
`perSnpEstimates` (a data frame of per-SNP beta_ZY/se_ZY). This is not
written by `exportSiteProfile()`. After CSV import, sensitivity analyses
(IVW, MR-Egger, weighted median) cannot be run from the imported profiles
alone.

**Fix**: Optionally export `perSnpEstimates` as a separate CSV. The
federated vignette already shows this as a manual step -- could be integrated
into `exportSiteProfile()`.

### 6. Steiger filtering uses point estimate only

**File**: `R/runSensitivityAnalyses.R`

The Steiger filter criterion is `abs(rExposure) > abs(rOutcome)` (point
estimate comparison). The Fisher z-test p-value is computed but not used for
filtering. Hemani et al. (2017) recommend using the p-value threshold to
establish statistical confidence in the causal direction. A SNP could pass
the filter with a trivially small difference in correlations.

**Fix**: Consider adding `steigerPass & steigerPval < 0.05` as the filtering
criterion, or at minimum document the choice.

### 7. `simulateSiteProfiles` missing `scoreDefinition`

**File**: `R/simulate.R`

Simulated profiles lack `scoreDefinition`, `betaHat`, and `seHat`. This
means the primary code path in `poolLikelihoodProfiles()` (score-definition
validation) and `computeMREstimate()` (using scoreDefinition for betaZX) is
not exercised by simulation-based tests and vignettes.

**Fix**: Add `scoreDefinition`, `betaHat`, and `seHat` to the simulated
profile objects.

### 8. `estimateSEFromProfile` fragile for non-uniform grids

**File**: `R/fitOutcomeModel.R`

The three-point finite difference for the second derivative at the profile
peak assumes uniform grid spacing. After interpolation in
`poolLikelihoodProfiles()`, the grid may not be perfectly uniform. The
asymmetric formula should be used, or a local quadratic fit over more points.

Low practical impact with the default grid `seq(-3, 3, by = 0.01)`.

### 9. Cyclops `predict()` return type unverified

**File**: `R/fitOutcomeModel.R`

In `evaluateBinaryProfilePoint()`, the Cyclops backend calls
`stats::predict(fit)` and uses the result as probabilities in
`dbinom(..., prob = fittedProb, ...)`. If Cyclops returns linear predictors
(log-odds) instead of response-scale probabilities, the log-likelihood
computation would be wrong. Should use `predict(fit, type = "response")` if
the Cyclops API supports it.

### 10. Strand-flip does not update allele labels

**File**: `R/harmonizeAlleles.R`

In the strand-flip direct-match case (Case 3, complement matches coded
allele), no allele update is performed. `effect_allele` and `other_allele`
remain the original alleles (e.g., "A") even though the genotype coded allele
is "T" (the complement). Effect estimates are correct but allele labels are
misleading after harmonization.

---

## Low-Severity / Design Notes

### 11. MR-Egger uses random-effects SE model

**File**: `R/runSensitivityAnalyses.R`

The internal MR-Egger uses `lm()` with weights, which incorporates the
residual variance (sigma-hat-squared) into the SE. This is the
"multiplicative random effects" version. The fixed-effect version would
assume sigma^2 = 1. The random-effects version is commonly used (and matches
TwoSampleMR default) but should be documented.

### 12. Weighted median bootstrap seed is fixed

**File**: `R/runSensitivityAnalyses.R`

`set.seed(42)` makes bootstrap SEs deterministic. Fine for reproducibility
but prevents assessing bootstrap variability.

### 13. Missing genotype imputation is crude

**File**: `R/fitOutcomeModel.R`

Missing genotypes are imputed as 0 (homozygous reference). Mean imputation
(2 * EAF) would reduce bias toward the null for common variants.

### 14. `testNegativeControls` is a stub

**File**: `R/runInstrumentDiagnostics.R`

Always returns an empty data frame. The `negativeControlFailure` flag can
never be TRUE. Should be implemented or documented as planned.

### 15. Cyclops `excludeIndices` assumes fixed column ordering

**File**: `R/fitOutcomeModel.R`

`excludeIndices = c(1L, 2L)` assumes intercept is at index 1 and the
exposure is at index 2 in the Cyclops design matrix. This matches R's
formula convention but is fragile if Cyclops reorders columns.

### 16. `buildMRCovariates` FeatureExtraction Andromeda access

**File**: `R/buildMRCovariates.R`

`length(unique(covariateData$covariateRef$covariateId))` may not work
directly on Andromeda-backed tables without `collect()`.

### 17. `computeAlleleScoreWeights` weighting convention

**File**: `R/fitOutcomeModel.R`

Weights are `beta_ZX / se_ZX^2` normalized by `sum(|w|)`. This is correct
per Burgess et al. (2013). The `sum(abs())` normalization preserves sign
structure for negative betas. Not a bug, but worth noting the convention.

### 18. Force-included SNPs bypass LD clumping

**File**: `R/getMRInstruments.R`

SNPs added via `forceIncludeSnps` are not checked for LD with clumped
instruments. Could introduce correlated instruments that violate MR
assumptions. A warning would be appropriate.

### 19. `simulateInstrumentTable` generates mean-zero betaZX

**File**: `R/simulate.R`

`betaZX <- rnorm(nSnps, mean = 0, sd = 0.3)` produces instruments with
near-zero effects, causing numerical issues in ratio computation.
`simulateMRData` avoids this with `runif(nSnps, 0.1, 0.5)`.

### 20. Delta method SE assumes zero covariance (two-sample)

**File**: `R/computeMREstimate.R`

The Wald ratio SE omits the covariance term between beta_ZX and beta_ZY.
This is correct for two-sample MR by design but should be documented as an
explicit assumption.

---

## What Is Correct

The following components are methodologically sound and well-implemented:

- **Profile likelihood via offset** (fitOutcomeModel.R): Standard GLM profile
  likelihood. Fixing the exposure as an offset and re-fitting nuisance
  parameters is the textbook approach.
- **Pointwise log-likelihood summation** (poolLikelihoodProfiles.R): Exact
  under independence across sites. The normalization by max(logLik) is
  numerically correct.
- **IVW estimator** (runSensitivityAnalyses.R): Matches Burgess et al. (2013)
  fixed-effect formula.
- **MR-Egger** (runSensitivityAnalyses.R): Correct orientation step, correct
  WLS regression with intercept.
- **Weighted median** (runSensitivityAnalyses.R): Matches Bowden et al. (2016)
  with parametric bootstrap SE.
- **Cochran's Q** (runSensitivityAnalyses.R): Standard heterogeneity
  statistic, correctly implemented.
- **Wald ratio and delta-method SE** (computeMREstimate.R): Correct for
  two-sample MR.
- **Likelihood-ratio CI** (computeMREstimate.R): Standard chi-squared(1)/2
  threshold, CI reordering for negative betaZX.
- **Allele harmonization cases 1-4** (harmonizeAlleles.R): All four
  allele-matching scenarios are handled correctly.
- **Allele score weights** (fitOutcomeModel.R): Standard IVW weighting with
  sign preservation.
- **Per-SNP as auxiliary output** (fitOutcomeModel.R): Correctly avoids summing
  incompatible per-SNP profiles; uses the allele-score profile as the poolable
  object.
- **Score definition validation** (poolLikelihoodProfiles.R): Critical guard
  that all sites used the same allele score is well-implemented.
- **Grid boundary and flatness diagnostics** (fitOutcomeModel.R): Sensible
  warnings for degenerate profiles.
