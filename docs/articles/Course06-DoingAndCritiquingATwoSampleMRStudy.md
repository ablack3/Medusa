<div id="main" class="col-md-9" role="main">

# Course 6: Doing and Critiquing a Two-Sample MR Study

<div class="section level2">

## What you’ll learn

By the end of this chapter, you should be able to:

-   describe an end-to-end two-sample MR study workflow,
-   identify major design and interpretation checkpoints,
-   use Medusa outputs to critique an analysis,
-   explain when an MR estimate should be treated cautiously.

</div>

<div class="section level2">

## Why this matters

By this point, you know the logic of MR and the mechanics of estimation.
The last step is judgment: can you design a credible study, read one
critically, and know when the output supports a cautious causal
interpretation?

</div>

<div class="section level2">

## A practical study checklist

Before running an MR study, ask:

1.  Is the exposure well defined and biologically coherent?
2.  Are the instruments strong enough?
3.  Are the outcome data appropriate for the question?
4.  Is allele harmonization trustworthy?
5.  Are the sensitivity analyses informative?
6.  Are the limitations stated honestly?

This checklist is often more important than the final p-value.

</div>

<div class="section level2">

## A small end-to-end synthetic analysis

<div id="cb1" class="sourceCode">

``` r
library(Medusa)

set.seed(1006)
beta_grid <- seq(-1.5, 1.5, by = 0.05)

sim_data <- simulateMRData(
  n = 5000,
  nSnps = 6,
  trueEffect = 0.30,
  seed = 1006
)

profile <- fitOutcomeModel(
  cohortData = sim_data$data,
  covariateData = NULL,
  instrumentTable = sim_data$instrumentTable,
  betaGrid = beta_grid,
  siteId = "site_A",
  analysisType = "perSNP"
)
#> Fitting outcome model at site 'site_A' (3141 cases, 1859 controls)...
#> Site 'site_A': beta_ZY_hat = 0.7014 (SE = 0.1192).

combined <- poolLikelihoodProfiles(list(site_A = profile))
#> Pooling profile likelihoods from 1 site(s)...
#> Pooling complete: 1 sites, 3141 total cases, 1859 total controls.
estimate <- computeMREstimate(combined, sim_data$instrumentTable)
#> MR estimate: beta = 1.9712 (95% CI: 1.4080, 2.5344), p = 3.52e-08
#> Odds ratio: 7.179 (95% CI: 4.088, 12.609)
sensitivity <- runSensitivityAnalyses(
  profile$perSnpEstimates,
  methods = c("IVW", "MREgger", "WeightedMedian", "LeaveOneOut"),
  engine = "internal"
)
#> Running sensitivity analyses with 6 SNPs...
#>   Engine: internal
#>   IVW...
#>   MR-Egger...
#>   Weighted Median...
#>   Leave-One-Out...
#> Sensitivity analyses complete.
```

</div>

<div class="section level3">

### Main estimate

<div id="cb2" class="sourceCode">

``` r
estimate$betaMR
#> [1] 1.971226
estimate$ciLower
#> [1] 1.408019
estimate$ciUpper
#> [1] 2.534433
estimate$oddsRatio
#> [1] 7.179473
```

</div>

</div>

<div class="section level3">

### Sensitivity summary

<div id="cb3" class="sourceCode">

``` r
sensitivity$summary
#>            method   beta_MR      se_MR   ci_lower  ci_upper         pval
#> 1             IVW 0.3500159 0.05251378  0.2470889 0.4529429 2.642698e-11
#> 2        MR-Egger 0.4332913 0.30920997 -0.1727602 1.0393429 2.337452e-01
#> 3 Weighted Median 0.4095011 0.07367214  0.2651037 0.5538985 2.722192e-08
```

</div>

</div>

</div>

<div class="section level2">

## How to critique the analysis

<div class="section level3">

### 1. Instrument quality

Check:

-   approximate F-statistics,
-   whether instruments come from a biologically sensible region,
-   whether the number of SNPs is enough to support meaningful
    sensitivity analyses.

</div>

<div class="section level3">

### 2. Harmonization and data integrity

Check:

-   allele labels,
-   palindromic SNP handling,
-   allele-frequency consistency,
-   genotype missingness.

</div>

<div class="section level3">

### 3. Estimator agreement

Check:

-   whether IVW, weighted median, and MR-Egger point in the same general
    direction,
-   whether one SNP dominates leave-one-out results,
-   whether uncertainty is very large.

</div>

<div class="section level3">

### 4. Interpretation discipline

Even a well-conducted MR result should be described carefully.

Good language:

-   “consistent with a causal effect,”
-   “supports a likely causal role,”
-   “evidence is compatible with…”

Overstated language:

-   “proves causation,”
-   “settles the question,”
-   “guarantees a drug will work.”

**Key takeaway:** the output of MR is evidence, not certainty.

</div>

</div>

<div class="section level2">

## Common reasons not to trust an MR estimate

-   weak instruments,
-   obvious horizontal pleiotropy,
-   severe disagreement across methods,
-   poor harmonization,
-   implausible biology,
-   strong chance of sample selection artifacts,
-   outcome definitions that do not match the scientific question.

</div>

<div class="section level2">

## Exercise set 1

<div class="section level3">

### Exercise 1

If IVW is positive, weighted median is near zero, and leave-one-out
changes sign when one SNP is removed, what should you conclude?

Show solution

You should conclude that the result is unstable and should not be
treated as robust causal evidence without deeper investigation. One SNP
may be overly influential, or the methods may be reacting to invalid
instruments.

</div>

<div class="section level3">

### Exercise 2

Why is “MR proves causation” usually too strong?

Show solution

Because MR depends on assumptions that can fail, especially exclusion
restriction and independence. It can strengthen causal inference
substantially, but it does not remove uncertainty completely.

</div>

<div class="section level3">

### Exercise 3

Name two things you would inspect before trusting a surprising MR
result.

Show solution

Reasonable answers include instrument strength, allele harmonization,
population-structure concerns, sensitivity-analysis agreement,
leave-one-out results, and biological plausibility.

</div>

</div>

<div class="section level2">

## Exercise set 2: critique a mock abstract

Mock claim:

> We found a statistically significant MR estimate linking biomarker X
> to disease Y, therefore biomarker X is definitively causal and should
> be targeted therapeutically.

What is missing from this claim?

Show solution

The claim ignores instrument validity, pleiotropy, sensitivity analyses,
harmonization quality, uncertainty, and the difference between causal
evidence and therapeutic success. A responsible abstract would discuss
the assumptions and supporting diagnostics rather than making a
definitive statement.

</div>

<div class="section level2">

## Where to go next

You now have the conceptual foundation to:

-   read ordinary two-sample MR papers critically,
-   understand why sensitivity analyses matter,
-   understand how Medusa’s federated workflow connects to textbook MR
    logic.

For next steps inside the package:

-   [Getting Started with Medusa](GettingStarted.md) for a runnable
    workflow,
-   [IL-6 Signaling and Colorectal
    Cancer](ColorectalCancerIL6Example.md) for a scientific walkthrough,
-   [Federated Analysis Guide](FederatedAnalysisGuide.md) for
    network-level deployment.

</div>

<div class="section level2">

## Chapter summary

-   A good MR study is judged by design quality and robustness, not only
    by the headline estimate.
-   You should inspect instrument strength, harmonization, diagnostics,
    and sensitivity analyses together.
-   Medusa provides tools that support this critical workflow, but the
    final judgment remains scientific and epidemiologic.

</div>

<div class="section level2">

## End of course

You have now moved from observational causal problems, to instrumental
variables, to two-sample MR mechanics, to sensitivity analyses, and
finally to Medusa’s federated implementation and critical
interpretation.

</div>

</div>
