<div id="main" class="col-md-9" role="main">

# Course 4: Estimation and Sensitivity Analysis

<div class="section level2">

## What you’ll learn

By the end of this chapter, you should be able to:

-   explain what IVW, MR-Egger, weighted median, leave-one-out, and
    Steiger filtering do,
-   describe what each method protects against,
-   explain why sensitivity analyses are necessary in MR,
-   run Medusa’s summarized-data sensitivity analysis workflow on
    synthetic data.

</div>

<div class="section level2">

## Why this matters

A single MR estimate is rarely enough. If your result changes
dramatically when one SNP is removed, or if methods with different
assumptions disagree sharply, you need to know that before interpreting
the analysis as causal evidence.

Sensitivity analyses exist because MR assumptions can fail in different
ways.

</div>

<div class="section level2">

## The main summarized-data estimators

<div class="section level3">

### IVW: inverse-variance weighted MR

IVW combines SNP-specific information into one overall estimate, giving
more weight to precise SNP-outcome associations.

It is often treated as the default summarized-data MR estimator when
most or all instruments are believed to be valid.

**What it is good at:** efficient estimation when the instruments are
broadly valid.

**What it is vulnerable to:** directional pleiotropy and influential
invalid instruments.

</div>

<div class="section level3">

### MR-Egger

MR-Egger allows a non-zero intercept in the weighted regression of
`beta_ZY` on `beta_ZX`.

Intuition:

-   if the intercept is far from zero, the SNPs may be pushing the
    outcome in a systematic way outside the exposure pathway.

**What it is good at:** detecting and sometimes adjusting for
directional pleiotropy under additional assumptions.

**What it is vulnerable to:** low precision and sensitivity to
weak-instrument problems.

</div>

<div class="section level3">

### Weighted median

The weighted median estimator takes the weighted median of the
SNP-specific ratio estimates.

**What it is good at:** giving a consistent estimate when less than half
of the total weight comes from invalid instruments.

**What it is vulnerable to:** strong violations concentrated in many
SNPs or complex dependence structures.

</div>

<div class="section level3">

### Leave-one-out

Leave-one-out removes one SNP at a time and re-runs the IVW estimate.

**What it is good at:** detecting whether one SNP dominates the result.

</div>

<div class="section level3">

### Steiger filtering

Steiger filtering asks whether a SNP appears to explain more variance in
the exposure than the outcome. If not, the proposed causal direction
becomes less convincing.

In Medusa:

-   the internal Steiger implementation is for continuous outcomes,
-   binary outcomes require different assumptions and are not handled by
    the internal version.

</div>

</div>

<div class="section level2">

## A small synthetic summary dataset

<div id="cb1" class="sourceCode">

``` r
library(Medusa)

set.seed(1004)
n_snps <- 8
beta_zx <- c(0.22, 0.18, 0.25, 0.20, 0.16, 0.24, 0.19, 0.21)
true_effect <- 0.45
beta_zy <- true_effect * beta_zx + rnorm(n_snps, sd = 0.015)

per_snp <- data.frame(
  snp_id = paste0("rs", seq_len(n_snps)),
  effect_allele = c("A", "C", "G", "T", "A", "C", "G", "A"),
  other_allele = c("G", "T", "A", "C", "C", "T", "T", "G"),
  eaf = seq(0.15, 0.50, length.out = n_snps),
  beta_ZY = beta_zy,
  se_ZY = rep(0.03, n_snps),
  beta_ZX = beta_zx,
  se_ZX = rep(0.04, n_snps),
  stringsAsFactors = FALSE
)

per_snp
#>   snp_id effect_allele other_allele  eaf    beta_ZY se_ZY beta_ZX se_ZX
#> 1    rs1             A            G 0.15 0.08988239  0.03    0.22  0.04
#> 2    rs2             C            T 0.20 0.09251444  0.03    0.18  0.04
#> 3    rs3             G            A 0.25 0.11003592  0.03    0.25  0.04
#> 4    rs4             T            C 0.30 0.08956835  0.03    0.20  0.04
#> 5    rs5             A            C 0.35 0.07220337  0.03    0.16  0.04
#> 6    rs6             C            T 0.40 0.09726425  0.03    0.24  0.04
#> 7    rs7             G            T 0.45 0.06333035  0.03    0.19  0.04
#> 8    rs8             A            G 0.50 0.10692022  0.03    0.21  0.04
```

</div>

</div>

<div class="section level2">

## Running Medusa sensitivity analyses

<div id="cb2" class="sourceCode">

``` r
sensitivity <- runSensitivityAnalyses(
  per_snp,
  methods = c("IVW", "MREgger", "WeightedMedian", "Steiger", "LeaveOneOut"),
  outcomeSampleSize = 15000,
  exposureSampleSize = 15000,
  outcomeType = "continuous",
  engine = "internal"
)
#> Running sensitivity analyses with 8 SNPs...
#>   Engine: internal
#>   IVW...
#>   MR-Egger...
#>   Weighted Median...
#>   Steiger filtering...
#> Warning in computeSteiger(perSnpEstimates, outcomeSampleSize,
#> exposureSampleSize, : All SNPs failed Steiger filter. Steiger-filtered estimate
#> not available. Investigate instrument validity.
#>   Leave-One-Out...
#> Sensitivity analyses complete.

sensitivity$summary
#>            method   beta_MR      se_MR   ci_lower  ci_upper         pval
#> 1             IVW 0.4362023 0.05095002 0.33634029 0.5360644 1.115270e-17
#> 2        MR-Egger 0.3720929 0.15241266 0.07336409 0.6708217 5.037909e-02
#> 3 Weighted Median 0.4401437 0.06266758 0.31731524 0.5629721 2.164288e-12
```

</div>

</div>

<div class="section level2">

## How to read the results

Questions to ask:

-   Are the IVW and weighted median estimates in the same rough
    direction?
-   Is the MR-Egger slope wildly different or very imprecise?
-   Is the MR-Egger intercept suggestive of directional pleiotropy?
-   Does leave-one-out show a single highly influential SNP?
-   Does Steiger filtering remove many SNPs?

Agreement across methods does not prove causality, but it strengthens
the case. Strong disagreement tells you to slow down and inspect the
assumptions.

**Key takeaway:** sensitivity analyses are not optional add-ons. They
are part of the causal argument.

</div>

<div class="section level2">

## Common misunderstanding

There is no single “bias-proof” MR estimator. Each method protects
against some problems under specific assumptions. The purpose of running
several methods is to see whether your inference is robust across
different failure scenarios.

</div>

<div class="section level2">

## Exercise set 1

<div class="section level3">

### Exercise 1

Which method is most directly designed to show whether one SNP is
dominating the overall result?

Show solution

Leave-one-out analysis, because it removes each SNP in turn and
recomputes the estimate.

</div>

<div class="section level3">

### Exercise 2

If the MR-Egger intercept is clearly non-zero, what concern does that
raise?

Show solution

It raises concern about directional pleiotropy, meaning the SNPs may be
affecting the outcome through pathways other than the exposure of
interest.

</div>

<div class="section level3">

### Exercise 3

Why might IVW and weighted median differ, even when both are computed
correctly?

Show solution

They make different assumptions about invalid instruments. Weighted
median can remain reliable when a minority of the weight comes from
invalid SNPs, whereas IVW can be more sensitive to those SNPs.

</div>

</div>

<div class="section level2">

## Exercise set 2: inspect a result table

Use the summary table above and answer:

1.  Which methods produce point estimates in the same general direction?
2.  Is there evidence that one method is much less precise than the
    others?

Show solution

You should usually see IVW, MR-Egger, and weighted median pointing in
the same general direction in this synthetic example. MR-Egger is often
less precise than IVW because it estimates both an intercept and a
slope.

</div>

<div class="section level2">

## Chapter summary

-   IVW is the standard efficient summarized-data estimator under
    broadly valid instruments.
-   MR-Egger checks for and can adjust for some directional pleiotropy.
-   Weighted median is robust when less than half of the total weight is
    invalid.
-   Leave-one-out identifies influential SNPs.
-   Steiger filtering checks whether the proposed causal direction looks
    plausible.

</div>

<div class="section level2">

## References

-   Burgess, S., Butterworth, A., & Thompson, S.G. (2013). Mendelian
    randomization analysis with multiple genetic variants using
    summarized data. *Genetic Epidemiology*, 37(7), 658–665.
-   Bowden, J., Davey Smith, G., & Burgess, S. (2015). Mendelian
    randomization with invalid instruments: effect estimation and bias
    detection through Egger regression. *International Journal of
    Epidemiology*, 44(2), 512–525.
-   Bowden, J., Davey Smith, G., Haycock, P.C., & Burgess, S. (2016).
    Consistent estimation in Mendelian randomization with some invalid
    instruments using a weighted median estimator. *Genetic
    Epidemiology*, 40(4), 304–314.
-   Hemani, G., Tilling, K., & Davey Smith, G. (2017). Orienting the
    causal relationship between imprecisely measured traits using GWAS
    summary data. *PLOS Genetics*, 13(11), e1007081.

</div>

<div class="section level2">

## Next chapter

Next: [Course 5: From textbook two-sample MR to the Medusa
workflow](Course05-FromTheoryToMedusaWorkflow.md)

For a more complete package example after this chapter, see [Getting
Started with Medusa](GettingStarted.md).

</div>

</div>
