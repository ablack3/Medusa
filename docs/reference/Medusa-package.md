<div id="main" class="col-md-9" role="main">

# Medusa: Mendelian Estimation in Distributed Standardized Analytics

<div class="ref-description section level2">

Medusa implements two-sample Mendelian Randomization (MR) within the
OHDSI ecosystem using the OMOP Common Data Model as the data substrate.
The package enables federated causal inference across distributed health
networks without requiring individual-level data to leave any site.

The core methodological innovation is one-shot federated pooling via
profile likelihood aggregation: each site computes a log-likelihood
profile across a grid of parameter values and shares only that vector of
numbers. The coordinator sums profiles across sites to obtain a pooled
estimate without any iterative communication protocol.

</div>

<div class="section level2">

## Details

The analysis pipeline consists of eight modules:

1.  **Instrument Assembly** (`getMRInstruments`): Query OpenGWAS for
    GWAS summary statistics and apply LD clumping.

2.  **Cohort Extraction** (`buildMRCohort`): Extract outcome cohorts and
    genotype data from OMOP CDM sites.

3.  **Covariate Assembly** (`buildMRCovariates`): Assemble covariates
    via FeatureExtraction for adjustment and diagnostics.

4.  **Instrument Diagnostics** (`runInstrumentDiagnostics`): Validate
    instruments via F-statistics, PheWAS, and negative controls.

5.  **Outcome Model** (`fitOutcomeModel`): Fit regularized outcome model
    and evaluate profile log-likelihood on a grid.

6.  **Likelihood Pooling** (`poolLikelihoodProfiles`): Aggregate
    site-level log-likelihood profiles via pointwise summation.

7.  **MR Estimation** (`computeMREstimate`): Compute Wald ratio estimate
    with delta method standard errors.

8.  **Sensitivity Analyses** (`runSensitivityAnalyses`): IVW, MR-Egger,
    weighted median, Steiger filtering, leave-one-out.

9.  **Reporting** (`generateMRReport`): Generate self-contained HTML
    report with all results and diagnostics.

</div>

<div class="section level2">

## References

Davey Smith, G., & Hemani, G. (2014). Mendelian randomization: genetic
anchors for causal inference in epidemiological studies. *Human
Molecular Genetics*, 23(R1), R89-R98.

Bowden, J., Davey Smith, G., & Burgess, S. (2015). Mendelian
randomization with invalid instruments: effect estimation and bias
detection through Egger regression. *International Journal of
Epidemiology*, 44(2), 512-525.

Bowden, J., Davey Smith, G., Haycock, P. C., & Burgess, S. (2016).
Consistent estimation in Mendelian randomization with some invalid
instruments using a weighted median estimator. *Genetic Epidemiology*,
40(4), 304-314.

</div>

<div class="section level2">

## Author

**Maintainer**: Adam Black <black@ohdsi.org>

</div>

</div>
