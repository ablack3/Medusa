<div id="main" class="col-md-9" role="main">

# Course 1: Why Mendelian Randomization?

<div class="section level2">

## What you’ll learn

By the end of this chapter, you should be able to:

-   explain why observational associations do not automatically imply
    causation,
-   describe confounding and reverse causation in plain language,
-   explain why randomization helps with causal inference,
-   describe the core idea of Mendelian randomization (MR) without using
    formulas.

</div>

<div class="section level2">

## Why this matters

Epidemiology often begins with a practical question: if one factor
changes, does disease risk really change too?

For example, people with higher low-density lipoprotein cholesterol
(LDL-C) tend to have more coronary artery disease (CAD). That
observation is useful, but it leaves open at least three possibilities:

1.  LDL-C causes CAD.
2.  Some third factor raises both LDL-C and CAD.
3.  Early disease processes change LDL-C, so the arrow points in the
    opposite direction.

MR is valuable because it tries to separate those possibilities using
genetic variation as a source of quasi-random assignment.

</div>

<div class="section level2">

## Association is not causation

Suppose we observe that people who carry umbrellas are more likely to be
wet. Umbrellas do not cause wetness. Rain causes both umbrella use and
wet clothes.

That is the basic causal challenge in observational data:

-   the exposure and outcome can move together,
-   but they may do so for reasons other than a direct causal effect.

<div class="section level3">

### A health example

Higher C-reactive protein (CRP) is often associated with worse
cardiovascular health. But what does that mean?

-   Maybe CRP is itself part of the causal pathway.
-   Maybe smoking, adiposity, or chronic illness raise CRP and
    independently raise cardiovascular risk.
-   Maybe developing disease raises inflammation, so CRP is partly an
    effect rather than a cause.

These are exactly the situations where traditional regression adjustment
helps, but may not fully solve the problem.

**Key takeaway:** strong association can still reflect confounding,
reverse causation, or both.

</div>

<div class="section level3">

### Common misunderstanding

It is easy to think that “we adjusted for confounders” means “we solved
confounding.” It does not. Adjustment only helps for confounders that
were measured well enough and modeled well enough.

</div>

</div>

<div class="section level2">

## Why randomized trials are so persuasive

Randomized controlled trials are powerful because treatment assignment
is independent of most patient characteristics. In expectation, the
treated and untreated groups differ mainly in the treatment itself.

That changes the causal question from:

-   “Are exposed people different from unexposed people?”

to:

-   “What happened after a random allocation changed exposure?”

This is the benchmark MR is trying to imitate.

</div>

<div class="section level2">

## The core idea of Mendelian randomization

Genes are allocated at conception. For many research questions, this
allocation happens before the disease develops and before most
later-life behaviors or clinical decisions occur.

If a genetic variant reliably changes an exposure, then the variant can
act as a kind of natural experiment:

-   some people inherit slightly higher lifelong exposure,
-   some inherit slightly lower lifelong exposure,
-   and that difference is not usually chosen by the individual.

This is the intuition behind Mendelian randomization.

In MR, a genetic variant is used as an **instrument**:

-   it nudges the exposure,
-   and we ask whether the outcome changes in a compatible way.

</div>

<div class="section level2">

## A simple trial analogy

Think of an LDL-raising variant.

-   In a trial, a drug assignment changes LDL-C.
-   In MR, inheritance of a variant changes LDL-C.

If the variant raises LDL-C and people carrying it also have higher CAD
risk, that pattern supports a causal role for LDL-C, provided the
variant is not operating through some other pathway.

The point is not that genes are literally a clinical trial. The point is
that they can sometimes create exposure differences that are much less
entangled with behavioral and social confounding than ordinary
observational measurements.

**Key takeaway:** MR uses genetic variants as instrumental variables to
mimic some of the causal leverage we get from randomization.

</div>

<div class="section level2">

## A concrete simulated example

The code below creates a small observational dataset where a confounder
affects both LDL-C and CAD risk. Then it adds a genetic variant that
also shifts LDL-C.

<div id="cb1" class="sourceCode">

``` r
set.seed(1001)
n <- 2000

confounder <- rnorm(n)
genetic_score <- rbinom(n, 2, 0.35)

ldl <- 0.8 * confounder + 0.4 * genetic_score + rnorm(n, sd = 0.8)
cad_linear_predictor <- 0.6 * ldl + 0.8 * confounder
cad <- rbinom(n, 1, plogis(-1 + cad_linear_predictor))

sim_df <- data.frame(
  confounder = confounder,
  genetic_score = genetic_score,
  ldl = ldl,
  cad = cad
)

round(cor(sim_df[, c("confounder", "genetic_score", "ldl", "cad")]), 2)
#>               confounder genetic_score  ldl  cad
#> confounder          1.00          0.01 0.70 0.45
#> genetic_score       0.01          1.00 0.25 0.09
#> ldl                 0.70          0.25 1.00 0.46
#> cad                 0.45          0.09 0.46 1.00
```

</div>

You should see:

-   the confounder correlated with both LDL-C and CAD,
-   the genetic score correlated with LDL-C,
-   the genetic score also related to CAD because it moves LDL-C.

This does not prove a valid MR design on its own, but it shows the basic
logic: we are looking for exposure variation that comes from inheritance
rather than from ordinary confounding structure.

</div>

<div class="section level2">

## Exercise set 1

<div class="section level3">

### Exercise 1

A study finds that people who drink more coffee have higher rates of
heartburn. Name one possible causal interpretation and one possible
non-causal interpretation.

Show solution

A causal interpretation is that coffee directly increases heartburn
risk. A non-causal interpretation is that stress or irregular meal
patterns cause both higher coffee intake and more heartburn. The
observation alone cannot separate those explanations.

</div>

<div class="section level3">

### Exercise 2

Why is reverse causation a problem in observational epidemiology?

Show solution

Reverse causation means the outcome or early disease process changes the
exposure, rather than the exposure changing the outcome. If that
happens, an observed association can look causal in the wrong direction.

</div>

<div class="section level3">

### Exercise 3

In one sentence, explain why randomization helps with causal inference.

Show solution

Randomization helps because it makes treatment assignment independent of
most patient characteristics, reducing systematic confounding between
groups.

</div>

</div>

<div class="section level2">

## From randomized trials to instrumental variables

MR is part of a broader family of methods called **instrumental
variable** methods. The instrument is something that changes the
exposure without being a plain confounder.

In the next chapter, we will make that idea precise and learn the three
core assumptions behind any valid MR analysis.

</div>

<div class="section level2">

## Chapter summary

-   Observational associations are not enough for causal inference.
-   Confounding and reverse causation can both create misleading
    associations.
-   Randomization is powerful because it changes exposure independently
    of many other determinants of the outcome.
-   MR uses genetic variants as instrumental variables to approximate
    that kind of leverage.

</div>

<div class="section level2">

## References

-   Davey Smith, G. & Ebrahim, S. (2003). ‘Mendelian randomization’: can
    genetic epidemiology contribute to understanding environmental
    determinants of disease? *International Journal of Epidemiology*,
    32(1), 1–22.
-   Lawlor, D.A., Harbord, R.M., Sterne, J.A.C., Timpson, N., & Davey
    Smith, G. (2008). Mendelian randomization: using genes as
    instruments for making causal inferences in epidemiology.
    *Statistics in Medicine*, 27(8), 1133–1163.
-   Davey Smith, G. & Hemani, G. (2014). Mendelian randomization:
    genetic anchors for causal inference in epidemiological studies.
    *Human Molecular Genetics*, 23(R1), R89–R98.

</div>

<div class="section level2">

## Next chapter

Next: [Course 2: Instrumental variables and MR
assumptions](Course02-InstrumentalVariablesAndMRAssumptions.md)

If you want a workflow-first overview after this chapter, see [Getting
Started with Medusa](GettingStarted.md).

</div>

</div>
