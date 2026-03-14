<div id="main" class="col-md-9" role="main">

# Course 3: The Mechanics of Two-Sample MR

<div class="section level2">

## What you’ll learn

By the end of this chapter, you should be able to:

-   explain what “two-sample” means,
-   interpret `beta_ZX` and `beta_ZY`,
-   compute a Wald ratio by hand,
-   explain why allele harmonization matters,
-   understand how multiple-SNP summarized-data MR is organized.

</div>

<div class="section level2">

## Why this matters

This chapter is where MR becomes operational. The main conceptual move
in two-sample MR is simple:

-   one dataset tells us how strongly each SNP predicts the exposure,
-   another dataset tells us how strongly the same SNP predicts the
    outcome.

Those two pieces are combined into a causal estimate.

</div>

<div class="section level2">

## What does “two-sample” mean?

In **two-sample MR**, the SNP-exposure and SNP-outcome associations are
estimated in different samples.

Typical pattern:

-   exposure sample: a large external GWAS,
-   outcome sample: your cohort, biobank, or federated network.

This does **not** mean that two hospitals are required. It means the
exposure and outcome associations come from distinct datasets.

</div>

<div class="section level2">

## The two key quantities

For a given SNP (Z):

-   (\_{ZX}): association of SNP (Z) with exposure (X),
-   (\_{ZY}): association of SNP (Z) with outcome (Y).

If the SNP affects the outcome only through the exposure, then a simple
causal estimate is:

\[ \_{MR} = \]

This is the **Wald ratio**.

<div class="section level3">

### Interpreting the ratio

If:

-   a SNP increases LDL-C by 0.20 units, and
-   the same SNP increases log-odds of CAD by 0.10,

then the Wald ratio is (0.10 / 0.20 = 0.50). That means:

-   the implied causal effect is 0.50 log-odds units of CAD per 1-unit
    increase in LDL-C.

</div>

</div>

<div class="section level2">

## One-SNP worked example

<div id="cb1" class="sourceCode">

``` r
beta_zx <- 0.20
beta_zy <- 0.10
beta_mr <- beta_zy / beta_zx
beta_mr
#> [1] 0.5
```

</div>

That is the most basic MR estimate possible.

</div>

<div class="section level2">

## Why the SNP signs must match

Suppose one dataset defines the effect allele as allele `A`, but the
other defines the effect allele as allele `G` at the same SNP. Then the
reported effect directions can point in opposite directions.

If you do not harmonize them first, the ratio can be wrong even when
both input studies were internally correct.

<div class="section level3">

### Mini harmonization example

<div id="cb2" class="sourceCode">

``` r
harmonization_example <- data.frame(
  snp_id = c("rs1", "rs2", "rs3"),
  effect_allele = c("A", "C", "T"),
  other_allele = c("G", "T", "C"),
  beta_ZX = c(0.20, -0.10, 0.15),
  beta_ZY = c(0.08, -0.06, 0.09),
  stringsAsFactors = FALSE
)

harmonization_example
#>   snp_id effect_allele other_allele beta_ZX beta_ZY
#> 1    rs1             A            G    0.20    0.08
#> 2    rs2             C            T   -0.10   -0.06
#> 3    rs3             T            C    0.15    0.09
```

</div>

If one SNP is coded with the opposite effect allele in one dataset, both
`beta_ZX` and `beta_ZY` for that SNP must be flipped together before the
ratio is computed.

**Key takeaway:** harmonization is not a cosmetic step. It is necessary
to keep the numerator and denominator talking about the same allele.

</div>

</div>

<div class="section level2">

## Multiple SNPs

Most MR analyses use more than one SNP.

For each SNP (j), we have:

\[ \_{MR,j} = \]

Those SNP-specific ratios can then be combined, often using
inverse-variance weighting.

Here is a simple summarized-data example:

<div id="cb3" class="sourceCode">

``` r
summary_data <- data.frame(
  snp_id = c("rsA", "rsB", "rsC"),
  beta_ZX = c(0.20, 0.12, 0.30),
  se_ZX = c(0.04, 0.03, 0.06),
  beta_ZY = c(0.10, 0.05, 0.17),
  se_ZY = c(0.03, 0.02, 0.05),
  effect_allele = c("A", "C", "T"),
  other_allele = c("G", "T", "C"),
  eaf = c(0.30, 0.45, 0.25),
  stringsAsFactors = FALSE
)

summary_data$wald_ratio <- summary_data$beta_ZY / summary_data$beta_ZX
summary_data
#>   snp_id beta_ZX se_ZX beta_ZY se_ZY effect_allele other_allele  eaf wald_ratio
#> 1    rsA    0.20  0.04    0.10  0.03             A            G 0.30  0.5000000
#> 2    rsB    0.12  0.03    0.05  0.02             C            T 0.45  0.4166667
#> 3    rsC    0.30  0.06    0.17  0.05             T            C 0.25  0.5666667
```

</div>

</div>

<div class="section level2">

## Why Medusa is a little different

Textbook summarized-data MR often combines SNP-specific ratio estimates
directly.

Medusa supports those methods in `runSensitivityAnalyses()`, but its
**primary federated estimator** is different:

-   sites estimate a pooled SNP-outcome signal through a weighted allele
    score,
-   they share a profile log-likelihood rather than person-level data,
-   the coordinator reconstructs the MR estimate from the pooled profile
    and the score’s SNP-exposure association.

So Medusa uses standard summarized-data MR ideas for sensitivity
analysis, but a profile-likelihood route for the main federated
estimate.

</div>

<div class="section level2">

## Exercise set 1

<div class="section level3">

### Exercise 1

If `beta_ZX = 0.25` and `beta_ZY = 0.05`, what is the Wald ratio?

Show solution

The Wald ratio is `0.05 / 0.25 = 0.20`.

</div>

<div class="section level3">

### Exercise 2

Why must `beta_ZX` and `beta_ZY` be aligned to the same effect allele
before forming a ratio?

Show solution

Because the ratio only makes sense when the numerator and denominator
describe the effect of the same allele. If one association is for allele
`A` and the other is for allele `G`, the directions can be inconsistent
and the ratio can be wrong.

</div>

<div class="section level3">

### Exercise 3

In one sentence, what does “two-sample” refer to in two-sample MR?

Show solution

It means the SNP-exposure associations and SNP-outcome associations are
estimated in different samples or datasets.

</div>

</div>

<div class="section level2">

## Exercise set 2: small coding practice

Use the small table below to compute SNP-specific Wald ratios in R.

<div id="cb4" class="sourceCode">

``` r
practice_df <- data.frame(
  snp = c("rs1", "rs2"),
  beta_ZX = c(0.10, 0.25),
  beta_ZY = c(0.03, 0.09)
)

practice_df$beta_MR <- practice_df$beta_ZY / practice_df$beta_ZX
practice_df
#>   snp beta_ZX beta_ZY beta_MR
#> 1 rs1    0.10    0.03    0.30
#> 2 rs2    0.25    0.09    0.36
```

</div>

Show solution

The key code is:

<div id="cb5" class="sourceCode">

``` r
practice_df$beta_MR <- practice_df$beta_ZY / practice_df$beta_ZX
practice_df
```

</div>

The resulting ratio estimates are 0.30 and 0.36. They are similar but
not identical, which is typical in multi-SNP MR because each
SNP-specific estimate contains sampling noise and may also be affected
by different biases.

</div>

<div class="section level2">

## Chapter summary

-   Two-sample MR combines SNP-exposure and SNP-outcome associations
    from different datasets.
-   `beta_ZX` describes SNP to exposure association.
-   `beta_ZY` describes SNP to outcome association.
-   The simplest causal estimate is the Wald ratio `beta_ZY / beta_ZX`.
-   Allele harmonization is essential.
-   Medusa uses conventional summarized-data logic for sensitivity
    analyses, but a federated profile-likelihood approach for its main
    pooled estimate.

</div>

<div class="section level2">

## References

-   Pierce, B.L. & Burgess, S. (2013). Efficient design for Mendelian
    randomization studies: subsample and 2-sample instrumental variable
    estimators. *American Journal of Epidemiology*, 178(7), 1177–1184.
-   Hartwig, F.P., Davies, N.M., Hemani, G., & Davey Smith, G. (2016).
    Two-sample Mendelian randomization: avoiding the downsides of a
    powerful, widely applicable but potentially fallible technique.
    *International Journal of Epidemiology*, 45(6), 1717–1726.
-   Hemani, G., Zheng, J., Elsworth, B., et al. (2018). The MR-Base
    platform supports systematic causal inference across the human
    phenome. *eLife*, 7, e34408.

</div>

<div class="section level2">

## Next chapter

Next: [Course 4: Estimation and sensitivity
analysis](Course04-EstimationAndSensitivityAnalysis.md)

If you want to see the package workflow in action before moving on, see
[Getting Started with Medusa](GettingStarted.md).

</div>

</div>
