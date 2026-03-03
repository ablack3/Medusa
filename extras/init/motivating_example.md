Let's walk through a realistic one. This is based on actual published research directions in the field.

## The Question

**Does elevated IL-6 signaling causally increase the risk of colorectal cancer?**

IL-6 is an inflammatory cytokine — a signaling protein involved in immune response. It's been observed repeatedly in epidemiological studies that people with high circulating IL-6 tend to have higher rates of colorectal cancer. But is that causal, or is it reverse causation (early tumors causing inflammation) or confounding (obese people have high IL-6 *and* high cancer risk, but obesity is the real driver)?

This matters enormously for drug development because there are already approved IL-6 inhibitors — tocilizumab, siltuximab — used in autoimmune disease. If IL-6 signaling genuinely causes colorectal cancer, these drugs or drugs like them might be repurposed for cancer prevention or treatment. That's a billion-dollar question and a relatively cheap one to answer with MR compared to running a prevention trial.

## Setting Up the Instrument

Genome-wide association studies (GWAS) have identified several genetic variants — SNPs — that reliably influence circulating IL-6 levels or IL-6 receptor signaling. A well-known one is a variant in the *IL6R* gene (rs2228145) that affects how well the IL-6 receptor binds IL-6. People who inherit this variant have measurably different IL-6 signaling their entire lives, and it was assigned randomly at conception.

This variant becomes your instrument. It's not perfectly clean — IL6R variants affect other inflammatory pathways too, which is a pleiotropy concern — but it's been used in published MR studies and is considered reasonably valid.

So your three conditions:
1. rs2228145 affects IL-6 signaling ✓ (well established in GWAS)
2. It affects colorectal cancer risk only through IL-6 signaling (assumed, tested with sensitivity analyses)
3. It's independent of obesity, diet, smoking (approximately true because it's inherited randomly)

## What You'd Need in the CDM

This is where your package comes in. Here's concretely what you'd be pulling:

**From OMOP CONDITION_OCCURRENCE:** Incident colorectal cancer cases, using standard SNOMED or ICD codes mapped to OMOP concepts. You'd define a careful phenotype — first diagnosis, excluding people with prior colorectal cancer, probably requiring two occurrences to reduce miscoding.

**From OMOP MEASUREMENT:** Serum IL-6 measurements if available, or CRP as a proxy (C-reactive protein is downstream of IL-6 signaling and more commonly measured). This would let you validate that the genetic instrument actually correlates with IL-6 in your population — called the first stage.

**From OMOP PERSON and OBSERVATION_PERIOD:** Age, sex, ancestry, follow-up time. Ancestry is critical because genetic variant frequencies differ by population and you need to account for that.

**Genomic linkage table:** rs2228145 genotype for each person_id. This is the piece that requires the CDM genomics extension. In a two-sample MR design, you could skip individual-level genotypes entirely and just use published GWAS summary statistics for IL-6, which makes the whole thing feasible even at sites without genomic data.

## Running the Analysis

In two-sample MR, you'd do roughly this:

**Sample 1** — from a large published GWAS: the association between rs2228145 and IL-6 levels. This is your exposure estimate. You get this from a public database like IEU OpenGWAS, not from your CDM at all.

**Sample 2** — from your OMOP network: the association between rs2228145 genotype and colorectal cancer diagnosis in your population. This is your outcome estimate.

You then combine them. The MR estimate is essentially: how much does a genetically-predicted unit increase in IL-6 signaling change the risk of colorectal cancer? If the answer is a significant odds ratio above 1.0, that's causal evidence.

Your package would handle extracting Sample 2 from the CDM, formatting it correctly, pulling Sample 1 from OpenGWAS via API, and passing both to TwoSampleMR in R to compute the estimate and run sensitivity analyses.

## What the Sensitivity Analyses Would Check

Because IL6R variants might affect other pathways, you'd run:

**MR-Egger regression** — detects directional pleiotropy. If the intercept is significantly non-zero, your instrument is affecting the outcome through pathways other than IL-6.

**Weighted median estimator** — gives a valid estimate even if up to 50% of your instruments are invalid. You'd use multiple IL6R-region SNPs as instruments, not just one.

**Leave-one-out analysis** — drops each SNP one at a time and checks whether the result is driven by a single outlier variant.

If all three sensitivity analyses point in the same direction as your main estimate, you have reasonably robust evidence.

## What the Result Would Mean

If your MR finds a significant positive association — higher genetically-predicted IL-6 signaling, higher colorectal cancer risk — that's human genetic evidence that IL-6 is causally involved. That would:

- Support further investigation of IL-6 inhibitors in colorectal cancer
- Potentially justify a repurposing trial of tocilizumab in high-risk populations
- Give a pharma company developing a novel IL-6 inhibitor additional confidence that their target is valid in oncology

If the result is null or negative, that's also valuable — it suggests the epidemiological association is likely confounded, and saves resources that would have gone into a doomed drug program.

## Why This Example Is Realistic

This isn't hypothetical — there are actually published MR studies on IL-6 and colorectal cancer risk using UK Biobank data. The results have been mixed, which is itself interesting and illustrates why having larger, more diverse populations through an OHDSI federated network would add real value. Running this across ten health systems with genomic linkage instead of one biobank would substantially increase power and generalizability.

That's exactly the gap your package would fill.

Does this make the end-to-end workflow feel concrete enough to start sketching what the package would actually need to do, step by step?
