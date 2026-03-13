<div id="main" class="col-md-9" role="main">

# Course 2: Instrumental Variables and MR Assumptions

<div class="section level2">

## What you’ll learn

By the end of this chapter, you should be able to:

-   define an instrumental variable in plain language,
-   state the three core MR assumptions,
-   explain weak instruments, pleiotropy, and population stratification,
-   recognize common threats to a valid MR design.

</div>

<div class="section level2">

## Why this matters

MR is persuasive only when the genetic variants really behave like
useful instruments. That is why serious MR work spends a lot of time
checking assumptions, not just calculating estimates.

</div>

<div class="section level2">

## What is an instrumental variable?

An **instrumental variable (IV)** is a variable that changes the
exposure in a way that can help us learn about the causal effect of the
exposure on the outcome.

For MR, the instrument is usually one or more genetic variants.

To be useful, the variant needs to do three things:

1.  predict the exposure,
2.  avoid tracking major confounders,
3.  affect the outcome mainly through the exposure.

Those are the three classic IV assumptions.

</div>

<div class="section level2">

## Assumption 1: Relevance

The instrument must genuinely affect the exposure.

If a variant barely changes LDL-C, it contains very little information
about the LDL-C to outcome pathway. In practice, such variants are
called **weak instruments**.

Weak instruments matter because they can make MR estimates unstable and
biased.

<div class="section level3">

### A simple strength check

MR studies often summarize instrument strength using the approximate
F-statistic:

\[ F ()^2 \]

Here:

-   (\_{ZX}) is the SNP-exposure association,
-   ((\_{ZX})) is its standard error.

A common rule of thumb is that (F &lt; 10) suggests a weak instrument.

<div id="cb1" class="sourceCode">

``` r
beta_zx <- c(0.20, 0.06, 0.30)
se_zx <- c(0.04, 0.05, 0.06)
data.frame(
  snp = c("rsA", "rsB", "rsC"),
  beta_zx = beta_zx,
  se_zx = se_zx,
  approx_F = (beta_zx / se_zx)^2
)
#>   snp beta_zx se_zx approx_F
#> 1 rsA    0.20  0.04    25.00
#> 2 rsB    0.06  0.05     1.44
#> 3 rsC    0.30  0.06    25.00
```

</div>

**Key takeaway:** if the instrument does not move the exposure, it
cannot tell you much about the exposure’s causal effect.

</div>

</div>

<div class="section level2">

## Assumption 2: Independence

The instrument should be independent of major confounders of the
exposure-outcome relationship.

Why might this fail?

-   **Population stratification:** ancestry differences can create
    associations between variants and environment.
-   **Selection bias:** the analytic sample may depend on both genotype
    and outcome-related factors.
-   **Dynastic or family effects:** parental genotype can shape the
    environment of offspring.

In practice, this assumption is motivated by the biology of inheritance,
but it is not automatically guaranteed in every dataset.

</div>

<div class="section level2">

## Assumption 3: Exclusion restriction

The instrument should affect the outcome mainly through the exposure of
interest.

This is usually the hardest assumption to defend.

The main threat is **pleiotropy**, which means a genetic variant affects
more than one trait.

There are two broad cases:

-   **Vertical pleiotropy:** the variant affects an upstream trait,
    which then affects downstream traits in the same causal pathway.
    This is usually not a problem for MR.
-   **Horizontal pleiotropy:** the variant affects the outcome through a
    separate pathway outside the exposure of interest. This can bias MR
    estimates.

<div class="section level3">

### Example

Suppose a variant changes LDL-C and also independently changes
inflammation. If inflammation affects CAD risk, then the variant has a
second route to the outcome. That violates exclusion restriction for an
LDL-C MR analysis.

**Key takeaway:** the instrument must not have an important back door
into the outcome outside the intended exposure pathway.

</div>

</div>

<div class="section level2">

## A picture in words

For a valid MR design, we want:

-   `genetic variant -> exposure -> outcome`

and we do **not** want:

-   `confounder -> genetic variant`
-   `genetic variant -> other pathway -> outcome`

That verbal graph is often enough to orient your thinking, even before
drawing a formal DAG.

</div>

<div class="section level2">

## A small simulation: strong vs weak instruments

<div id="cb2" class="sourceCode">

``` r
set.seed(1002)
n <- 5000
strong_snp <- rbinom(n, 2, 0.3)
weak_snp <- rbinom(n, 2, 0.3)

exposure_strong <- 0.40 * strong_snp + rnorm(n)
exposure_weak <- 0.04 * weak_snp + rnorm(n)

summary(lm(exposure_strong ~ strong_snp))$coefficients["strong_snp", ]
#>     Estimate   Std. Error      t value     Pr(>|t|) 
#> 4.025216e-01 2.185009e-02 1.842196e+01 2.242399e-73
summary(lm(exposure_weak ~ weak_snp))$coefficients["weak_snp", ]
#>    Estimate  Std. Error     t value    Pr(>|t|) 
#> 0.003518649 0.022374208 0.157263618 0.875043445
```

</div>

The strong variant produces a clearer exposure shift. The weak variant
produces a noisier, less informative relationship.

</div>

<div class="section level2">

## Exercise set 1

<div class="section level3">

### Exercise 1

Which MR assumption is threatened most directly by a variant that barely
changes the exposure?

Show solution

The relevance assumption is threatened. If the variant barely affects
the exposure, it is a weak instrument.

</div>

<div class="section level3">

### Exercise 2

A variant affects body mass index and, independently, smoking behavior.
You are using it to study body mass index and lung disease. Which
assumption is at risk?

Show solution

The exclusion restriction is at risk, because the variant may affect
lung disease through smoking behavior, not only through body mass index.

</div>

<div class="section level3">

### Exercise 3

Why can population stratification threaten MR, even though genes are
inherited at conception?

Show solution

If ancestry groups differ in both allele frequencies and environmental
or health-related factors, genotype can become associated with
confounders at the sample level. That weakens the independence
assumption.

</div>

</div>

<div class="section level2">

## How Medusa helps with assumptions

Medusa cannot prove the assumptions, but it gives you tools that support
assumption checking:

-   `runInstrumentDiagnostics()` reports approximate F-statistics,
-   the instrument PheWAS looks for associations with measured
    covariates,
-   allele-frequency comparisons can flag harmonization or population
    issues,
-   genotype missingness summaries can reveal data-quality concerns.

Those diagnostics do not replace scientific judgment, but they improve
transparency.

</div>

<div class="section level2">

## Common misunderstanding

It is tempting to say “because genes are randomized, MR assumptions are
automatically true.” That is too strong. Inheritance helps, but weak
instruments, pleiotropy, population structure, and sample selection can
still cause problems.

</div>

<div class="section level2">

## Chapter summary

-   MR is an instrumental variable method using genetic variants as
    instruments.
-   The three core assumptions are relevance, independence, and
    exclusion restriction.
-   Weak instruments threaten relevance.
-   Pleiotropy most often threatens exclusion restriction.
-   Population stratification and selection can threaten independence.

</div>

<div class="section level2">

## Next chapter

Next: [Course 3: How two-sample MR works in
practice](Course03-TwoSampleMRMechanics.md)

For a package-oriented preview of diagnostics, see [Getting Started with
Medusa](GettingStarted.md) and [IL-6 Signaling and Colorectal
Cancer](ColorectalCancerIL6Example.md).

</div>

</div>
