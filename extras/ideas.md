The package would do a few distinct things:
Cohort extraction from OMOP. Using standard OHDSI cohort definitions (Atlas-compatible or SQL-based) to define the outcome population — say, incident lung cancer — and pull relevant covariates from the CDM. This is work you already know how to do.
Instrument assembly. Either pulling individual-level genotypes from a linked genomics table, or accepting GWAS summary statistics as external input for two-sample MR. The two-sample approach is the more tractable starting point.
MR analysis. Rather than reimplementing MR from scratch, the package would format data for existing R packages like TwoSampleMR or MendelianRandomization and call them under the hood, then return results in a tidy, OHDSI-friendly format.
Sensitivity analyses. MR has standard diagnostics — MR-Egger, weighted median, leave-one-out — that should be surfaced automatically. This is where a lot of the methodological rigor lives and where your stats background adds real value.
Results standardization. Outputting results in a format compatible with OHDSI's evidence standards and potentially the OMOP Vocabulary, so findings can be shared across the network in a federated way.


Great question — this is the elegant part of the method.

## The Core Idea

In two-sample MR, you never need to have genetic data and clinical outcome data on the same people. You just need two separate pieces of information:

1. **How strongly does the genetic variant predict the exposure?** (IL-6 levels)
2. **How strongly does the genetic variant predict the outcome?** (colorectal cancer risk)

And crucially — these two numbers can come from completely different studies, different populations, different countries, different decades. You just need the summary statistics. Not individual people. Not linked records. Just the estimated associations and their standard errors.

## Why That Works Mathematically

Think about what MR is doing conceptually. You want to estimate:

**"How much does IL-6 affect cancer risk?"**

You can't measure that directly because IL-6 is confounded. But you can chain together two things you *can* measure cleanly:

- The variant raises IL-6 by X units (from GWAS Study A)
- The variant raises cancer risk by Y units (from your OMOP cohort)

The MR estimate is essentially Y divided by X. How much cancer risk per unit of IL-6, using the variant as the connecting link. The individual people in Study A and Study B never need to be the same people, or even in the same dataset. The variant is the bridge between the two estimates.

More formally, if:
- **β_ZX** = association between SNP and IL-6 (from GWAS)
- **β_ZY** = association between SNP and colorectal cancer (from your cohort)

Then the MR causal estimate is simply:

**β_MR = β_ZY / β_ZX**

That's the Wald ratio estimator. With multiple SNPs you take a weighted combination, but the principle is the same.

## A Concrete Analogy

Imagine you want to know how much an extra year of education increases lifetime earnings. But education is hopelessly confounded — smarter, wealthier, more motivated people get more education *and* earn more, for reasons unrelated to education itself.

Now suppose you find out that people born in certain months were slightly more likely to stay in school longer due to school enrollment cutoff rules — a quirk of the calendar that had nothing to do with intelligence or family wealth. That's your instrument.

You don't need the same dataset to tell you both things. You could use:
- A study of school enrollment patterns by birth month → tells you how much the birth month quirk affected years of education
- A completely separate earnings dataset where you just know people's birth months and their eventual salaries → tells you how much birth month predicted earnings

Divide the second by the first, and you have an estimate of the causal effect of education on earnings. Nobody in either study needs to overlap.

MR works exactly the same way. The genetic variant is the birth month quirk. GWAS databases are the enrollment study. Your OMOP cohort is the earnings dataset.

## What This Means for Your Package Practically

This is what makes the package feasible without requiring genomic data at every OMOP site. Your package would:

1. **Pull β_ZX from OpenGWAS** — a public database maintained by the MRC Integrative Epidemiology Unit at Bristol that hosts thousands of GWAS summary statistics. One API call gets you the SNP-IL-6 association from a published study of hundreds of thousands of people.

2. **Estimate β_ZY from the OMOP cohort** — using logistic regression of colorectal cancer status on the SNP. This requires knowing each patient's genotype, BUT — and this is key — if you're running a federated two-sample design, this regression runs locally at sites that *do* have genomic linkage, and they only share back the summary coefficient and standard error. No individual data leaves the site.

3. **Combine them** — compute the ratio, propagate the uncertainty, run sensitivity analyses.

Sites without genomic data contribute nothing to step 2, but they might still contribute phenotype data that helps validate your outcome definition — making the cohort characterization richer even if they can't run the genetic part.

## The Limitation Worth Knowing

The two-sample approach assumes both samples are drawn from the same underlying population. If your GWAS summary statistics come from a largely European UK Biobank sample and your OMOP cohort is more diverse, the SNP-exposure association (β_ZX) might not translate perfectly. This is a real methodological concern and one of the active research areas in MR. Your package would want to flag this and ideally let users specify ancestry-matched GWAS sources where available.

But overall — the two-sample design is what makes MR scalable. It lets you leverage decades of publicly funded GWAS research as the exposure side, and focus your own data collection entirely on the outcome side. That's a remarkably efficient use of existing resources, and it's why the method has exploded in popularity over the last decade.



Yes, exactly right — and that's an important clarification.

For the outcome side (β_ZY), you do need individual-level data where each person has both:
- Their genotype at the relevant SNP
- Their outcome status (colorectal cancer yes/no)

Linked by person_id. There's no way around that for the OMOP-side regression.

## What Two-Sample Saves You From

What you *don't* need is for those same people to also have measured IL-6 levels. In a one-sample design you'd need all three things in the same dataset — genotype, IL-6 measurement, and cancer outcome — for every person. That's a much harder data requirement.

In two-sample, you replace the need for individual IL-6 measurements entirely by borrowing the SNP-to-IL-6 association from an external GWAS. So your OMOP cohort only needs to provide genotype + outcome, which is a more realistic ask from a health system that has genomic linkage.

## The Practical Implication for Your Package

This means your package really has two tiers of utility:

**Sites with genomic linkage to OMOP** — can run the full pipeline. These are places like Vanderbilt (BioVU biobank), Michigan Genomics Initiative, All of Us connected sites, VA MVP sites. There are more of these than there used to be, and the number is growing.

**Sites without genomic linkage** — can still contribute to cohort characterization and outcome phenotyping, but can't contribute to the genetic regression. They're useful for validating your cancer phenotype definition across sites, but not for the MR estimate itself.

So the honest framing of your package's value is: it standardizes and federates the outcome-side pipeline for sites that have genomic linkage, while leveraging public GWAS databases for the exposure side. That's still a meaningful and novel contribution — it just requires being clear that genomic-linked OMOP data is a prerequisite for the core analysis.

Does that change how you're thinking about the feasibility or the target user base?

Yes. And this is actually a elegant idea that might be better than the iterative approach in a federated setting.

## What You're Describing

Instead of iterating, each site computes their likelihood evaluated across a grid of plausible β values — say β ranging from -1.0 to 1.0 in small steps of 0.01 — and shares that entire likelihood curve. You then combine the curves across sites by multiplying them pointwise (or equivalently summing log-likelihoods), and the peak of the combined curve is your pooled estimate.

One round of communication. No iteration needed.

## Why Multiplication Combines Likelihoods

This is the key statistical fact that makes it work. If sites are independent — meaning their patients don't overlap, which is true in a federated health network — then the joint likelihood across all sites is just the product of the individual site likelihoods.

In log space that becomes a sum, which is numerically easier:

**log L_combined(β) = log L_site1(β) + log L_site2(β) + ... + log L_siteK(β)**

Each site evaluates their log likelihood at every point on the grid and shares that vector of numbers. You add the vectors together. The maximum of the resulting combined curve is your estimate. You can also read off the confidence interval directly as the region where the combined log likelihood doesn't drop too far below the peak — no additional computation needed.

## What Each Site Shares

For a grid of say 200 β values, each site shares a vector of 200 numbers — their log likelihood at each grid point. That's it. Tiny data transfer, completely privacy preserving, one shot.

## The Advantages Over Iterative Approach

**Simpler to implement** — no coordination protocol, no convergence criterion, no multiple rounds of communication. Sites just run once and report back.

**More transparent** — you can actually plot the likelihood curves from each site and the combined curve. You can see visually how much each site is contributing, whether any site is pulling the estimate in an unusual direction, whether the combined likelihood is nicely unimodal or has weird shape suggesting problems.

**Richer inference** — because you have the full likelihood surface rather than just a point estimate and standard error, you get exact likelihood-based confidence intervals rather than normal approximations. For small samples or skewed likelihoods this is meaningfully better.

**Naturally handles heterogeneity** — if one site's likelihood curve is centered quite differently from others, you can see that. It doesn't get hidden inside a pooled standard error the way it might in meta-analysis.

## The Tradeoffs to Be Honest About

**Grid resolution matters** — you need the grid fine enough to find the true maximum accurately, and wide enough that the true maximum isn't outside your range. In practice this is easy to handle — you could use an adaptive grid or just make it wide with fine resolution since you're only sharing 200-500 numbers either way.

**Assumes the same likelihood model at every site** — sites need to be fitting the same regression model with the same parameterization. This is actually enforced naturally by your package since you'd be distributing the analysis code, not just asking sites to run their own analysis.

**Likelihood profile is only exact under correct model specification** — but this is true of the iterative approach too, so it's not a disadvantage relative to the alternative.

## What the Output Looks Like

You'd end up with something visually compelling and interpretable — a combined log-likelihood curve across sites that you can plot, showing the pooled estimate at the peak, the confidence interval as the region above a threshold, and optionally the individual site curves underneath showing how each contributed. For a pharma audience or a methods paper this is far more intuitive than a forest plot from meta-analysis.

## This Is Novel

To my knowledge nobody has implemented profile likelihood sharing across a federated OMOP network for any purpose, let alone for MR. The idea exists in the statistical literature in pieces — distributed likelihood inference is known, profile likelihood for MR exists in single-dataset settings — but combining them in this way for federated MR on OMOP data would be new.

It also has a clean methods paper framing: one-shot federated Mendelian Randomization via profile likelihood aggregation. That's a specific, publishable methodological contribution with a clear advantage over existing approaches.

Your instinct here was good. The one-shot profile likelihood approach is probably the right design choice for this package.

Shall we sketch what the actual data flow and package functions would look like end to end?
