# Study Specification: GLP-1 Receptor Agonism and Cancer Risk
## A Drug-Target Mendelian Randomization Study Using Medusa

*Version 1.0 | 2026-03-14*

---

## 1. Study Objective

To estimate the causal effect of genetically proxied GLP-1 receptor (GLP-1R)
activation on risk of obesity-related cancers using drug-target Mendelian
Randomization, with replication in a diverse population (All of Us) and
empirical calibration via negative control outcomes.

### Primary Question

Does genetically proxied GLP-1R activation reduce or increase risk of breast,
colorectal, endometrial, hepatocellular, or pancreatic cancer?

### Secondary Questions

1. Are effects consistent across ancestries (European, African, Hispanic)?
2. Is the cancer effect mediated by BMI reduction or direct GLP-1R signaling?
3. Does the effect differ by cancer subtype (e.g., ER+/ER- breast cancer)?

---

## 2. Genetic Instruments

### Source

cis-eQTLs for the *GLP1R* gene from the **eQTLGen Consortium** (blood samples
from 31,684 European-ancestry individuals).

### Gene Region

- Gene: *GLP1R*
- Chromosome: 6
- Position: 39,016,557 -- 39,059,079 (GRCh38)
- Instrument window: +/- 100 kb (or +/- 1 Mb for broader window)

### Instrument Selection Criteria

1. Significantly associated with *GLP1R* expression in blood (p < 5 x 10^-8)
2. Within +/- 100 kb of the *GLP1R* gene
3. LD clumped (r^2 < 0.3 within 250 kb window, using 1000 Genomes EUR reference)
4. F-statistic > 10 for each instrument
5. Minor allele frequency > 1%

### Expected Instruments

Based on published studies, approximately **12--22 independent cis-eQTL SNPs**
are expected after clumping. These instruments have been validated in prior
drug-target MR studies with:
- Positive control: reduced T2DM risk (OR ~ 0.82)
- Positive control: reduced BMI (OR ~ 0.95)

### Alternative/Supplementary Instruments

- **cis-pQTLs** for GLP1R protein from the INTERVAL study (n = 3,301)
- **HbA1c-associated variants** near *GLP1R* from UK Biobank (used in
  multivariable MR to separate glycemic from weight-loss effects)

### Key Pharmacogenomic Variants (for context, not as instruments)

- rs6923761: GLP1R missense variant (MAF ~ 0.37), studied for semaglutide response
- rs3765467: GLP1R missense variant, associated with early-onset T2DM risk

### Instrument Retrieval in Medusa

```r
# Option 1: From OpenGWAS/eQTLGen via Medusa
instruments <- getMRInstruments(
  traitId = "eqtl-a-ENSG00000112164",  # GLP1R eQTL from eQTLGen
  pThreshold = 5e-8,
  clumpR2 = 0.3,
  clumpKb = 250
)

# Option 2: Manual instrument table from published supplementary data
# (preferred for reproducibility)
instruments <- data.frame(
  snp_id = c("rs...", "rs...", ...),    # From eQTLGen cis-eQTL results
  effect_allele = c(...),
  other_allele = c(...),
  beta_ZX = c(...),                      # Effect on GLP1R expression
  se_ZX = c(...),
  pval_ZX = c(...),
  eaf = c(...),
  gene_region = rep("GLP1R", n)
)
```

---

## 3. Outcome Definitions (OMOP Concept Sets)

### Primary Cancer Outcomes

| Cancer | OMOP Concept ID | Concept Name | SNOMED Code | Record Count |
|--------|----------------|--------------|-------------|--------------|
| Breast cancer | 4112853 | Malignant tumor of breast | 254837009 | 1,755,760 |
| Colorectal cancer (colon) | 4180790 | Malignant tumor of colon | 363406005 | 282,420 |
| Colorectal cancer (rectum) | 443390 | Malignant neoplasm of rectum | 363351006 | 96,420 |
| Endometrial cancer | 4110871 | Endometrial carcinoma | 254878006 | 95,210 |
| Hepatocellular carcinoma | 4001171 | Liver cell carcinoma | 109841003 | 4,459,600 |
| Pancreatic cancer | 4180793 | Malignant tumor of pancreas | 363418001 | 64,210 |

### Secondary/Exploratory Cancer Outcomes

| Cancer | OMOP Concept ID | Concept Name | Record Count |
|--------|----------------|--------------|--------------|
| Ovarian cancer | 4181351 | Malignant neoplasm of ovary | 100,020 |
| Renal cell carcinoma | 45765451 | Renal cell carcinoma | 609,960 |
| Thyroid cancer | 40488900 | Carcinoma of thyroid | 2,220 |
| Basal cell carcinoma | 4112752 | Basal cell carcinoma of skin | 5,097,970 |

### Cancer Subtype Concepts (for secondary analyses)

| Subtype | OMOP Concept ID | Concept Name | Record Count |
|---------|----------------|--------------|--------------|
| Adenocarcinoma of colon | 42872396 | Primary adenocarcinoma of colon | 18,520 |
| Adenocarcinoma of pancreas | 45763891 | Adenocarcinoma of pancreas | 167,150 |
| Adenocarcinoma of endometrium | 4048226 | Adenocarcinoma of endometrium | 151,410 |
| Papillary thyroid carcinoma | 4116228 | Papillary thyroid carcinoma | 412,970 |
| Clear cell carcinoma of kidney | 4110880 | Clear cell carcinoma of kidney | 154,370 |

### Cohort Definition Strategy

For each cancer outcome, the cohort definition should:

1. Use the parent SNOMED concept + all descendant concepts to capture all
   coding variations
2. Require at least 1 condition occurrence in the `condition_occurrence` table
3. Use the **first** occurrence date as the index date
4. Exclude prevalent cases (require a minimum washout period of 365 days)
5. Consider requiring a second confirmatory diagnosis within 90 days for
   specificity (sensitivity analysis)

---

## 4. Positive Control Outcomes

These outcomes have established causal relationships with GLP-1R activation
and serve to validate instrument strength.

| Outcome | Expected Direction | OMOP Concept ID | Concept Name |
|---------|--------------------|----------------|--------------|
| Type 2 diabetes | Protective (OR ~ 0.82) | 201826 | Type 2 diabetes mellitus |
| Obesity | Protective (OR ~ 0.95 for BMI) | 4215968 | Obesity |

---

## 5. Negative Control Outcomes

Outcomes with no plausible causal pathway from GLP-1R activation. These are
used for empirical calibration via Medusa's `runNegativeControlAnalysis()`.

**Selection criteria**: No known biological mechanism linking GLP-1R signaling,
glucose metabolism, insulin sensitivity, body weight, or GI motility to the
outcome. Common enough to have sufficient cases. Recordable in EHR data.

| # | Outcome | OMOP Concept ID | Concept Name | Rationale |
|---|---------|----------------|--------------|-----------|
| 1 | Fracture of bone | 75053 | Fracture of bone | Mechanical/trauma, not metabolic |
| 2 | Osteoarthritis | 80180 | Osteoarthritis | Degenerative joint disease |
| 3 | Cataract | 375545 | Cataract | Lens opacity, age-related |
| 4 | Appendicitis | 440448 | Appendicitis | Infectious/anatomical |
| 5 | Inguinal hernia | 4288544 | Inguinal hernia | Anatomical/mechanical |
| 6 | Ingrowing nail | 139099 | Ingrowing nail | Mechanical/local |
| 7 | Sensorineural hearing loss | 374366 | Sensorineural hearing loss | Cochlear degeneration |
| 8 | Carpal tunnel syndrome | 380094 | Carpal tunnel syndrome | Nerve compression |
| 9 | Allergic rhinitis | 257007 | Allergic rhinitis | IgE-mediated allergy |
| 10 | Seborrheic dermatitis | 137053 | Seborrheic dermatitis | Skin/fungal |
| 11 | Kidney stones | 4148260 | Calculus of kidney and ureter | Mineral precipitation |
| 12 | Varicose veins | 4141501 | Simple varicose veins | Venous valve incompetence |
| 13 | Dental caries | 133228 | Dental caries | Bacterial/dietary |
| 14 | Rotator cuff tear | 4172970 | Nontraumatic rotator cuff tear | Tendon degeneration |
| 15 | Migraine | 318736 | Migraine | Neurovascular |
| 16 | BPPV | 81878 | Benign paroxysmal positional vertigo | Otolith displacement |

**Exclusions from negative controls** (outcomes potentially on GLP-1R pathway):
- Cardiovascular outcomes (GLP-1R has known cardioprotective effects)
- Pancreatitis (GLP-1RA label warning)
- Gallbladder disease (associated with rapid weight loss)
- Depression/anxiety (GLP-1R expressed in brain; active MR research area)
- Diabetic retinopathy, neuropathy (downstream of glucose effects)

---

## 6. Analysis Plan

### 6.1 Primary Analysis: Allele Score MR

Using Medusa's profile likelihood approach:

```r
library(Medusa)

# Step 1: Get instruments
instruments <- getMRInstruments(...)

# Step 2: Build MR cohort (at each site)
cohortData <- buildMRCohort(
  connectionDetails = connectionDetails,
  cdmDatabaseSchema = "cdm_schema",
  cohortDatabaseSchema = "results_schema",
  cohortTable = "cohort",
  outcomeCohortId = 4112853,  # e.g., breast cancer
  instrumentTable = instruments,
  negativeControlCohortIds = c(75053, 80180, 375545, 440448,
                                4288544, 139099, 374366, 380094,
                                257007, 137053, 4148260, 4141501,
                                133228, 4172970, 318736, 81878)
)

# Step 3: Build covariates
covariateData <- buildMRCovariates(
  connectionDetails = connectionDetails,
  cdmDatabaseSchema = "cdm_schema",
  cohortData = cohortData
)

# Step 4: Run diagnostics
diagnostics <- runInstrumentDiagnostics(
  cohortData = cohortData,
  covariateData = covariateData,
  instrumentTable = instruments
)

# Step 5: Fit outcome model (profile likelihood)
siteProfile <- fitOutcomeModel(
  cohortData = cohortData,
  instrumentTable = instruments,
  covariateData = covariateData,
  betaGrid = seq(-3, 3, by = 0.01)
)

# Step 6: Export profile for federated pooling
exportSiteProfile(siteProfile, outputDir = "export/")
```

### 6.2 Federated Pooling (Coordinator)

```r
# Pool profiles from All of Us, MVP, eMERGE, etc.
siteProfiles <- importSiteProfiles("all_site_exports/")

combined <- poolLikelihoodProfiles(siteProfiles)

estimate <- computeMREstimate(combined, instruments)

sensitivity <- runSensitivityAnalyses(
  combinedProfile = combined,
  instrumentTable = instruments,
  perSnpEstimates = perSnpResults
)
```

### 6.3 Negative Control Calibration

```r
ncResults <- runNegativeControlAnalysis(
  cohortData = cohortData,
  instrumentTable = instruments,
  primaryEstimate = list(
    betaMR = estimate$betaMR,
    seMR = estimate$seMR,
    pValue = estimate$pValue
  )
)

# Calibrated p-value and CI
ncResults$calibratedPrimary
ncResults$biasDetected
```

### 6.4 Sensitivity Analyses

1. **MR-Egger**: Test for directional pleiotropy (intercept test)
2. **Weighted median**: Robust to up to 50% invalid instruments
3. **MR-PRESSO**: Outlier detection and correction
4. **Cochran's Q / I-squared**: Instrument heterogeneity
5. **Leave-one-out**: Identify influential individual SNPs
6. **Multivariable MR**: Separate GLP-1R effect from BMI-mediated effect
   using HbA1c-associated and BMI-associated variants in the GLP-1R region

### 6.5 Subgroup Analyses

1. **By ancestry**: European, African, Hispanic (All of Us enables this)
2. **By sex**: Particularly for breast and endometrial cancer
3. **By age**: <60 vs >=60 at cohort entry
4. **By cancer subtype**: Where OMOP coding allows (e.g., adenocarcinoma vs.
   other histologies)

---

## 7. Expected Sample Sizes

### All of Us (spring 2026 release, n ~ 535,000 WGS)

| Cancer | Estimated prevalence | Expected cases |
|--------|---------------------|----------------|
| Breast cancer | 1.5-2.0% | 5,000-8,000 |
| Colorectal cancer | 0.5-1.0% | 2,500-5,000 |
| Endometrial cancer | 0.2-0.4% | 1,000-2,000 |
| Hepatocellular carcinoma | 0.05-0.1% | 250-500 |
| Pancreatic cancer | 0.05-0.1% | 250-500 |
| Basal cell carcinoma | 2-4% | 10,000-20,000 |
| Ovarian cancer | 0.1-0.2% | 500-1,000 |
| Renal cell carcinoma | 0.2-0.4% | 1,000-2,000 |

**Note**: All of Us is a "sicker" population than general US demographics
(per program director Josh Denny), so cancer prevalence may be higher.
Cancer registry linkage (planned) will improve case ascertainment.

### Federated (All of Us + MVP + eMERGE, n ~ 1.5M)

Case counts roughly triple, making even pancreatic cancer and HCC well-powered.

---

## 8. Power Considerations

For a binary outcome with allele score MR:

- With 5,000 breast cancer cases, 12-22 instruments explaining ~1% of GLP1R
  expression variance, power is >80% to detect OR of 0.85 or 1.15 at alpha = 0.05
- For rarer cancers (500 cases), minimum detectable OR is ~0.70 or 1.30
- Federated pooling substantially improves power for rare outcomes

---

## 9. Prior MR Evidence to Compare Against

| Cancer | Sun et al. 2024 MR finding | Direction | Our hypothesis |
|--------|---------------------------|-----------|----------------|
| Breast | Decreased risk | Protective | Replicate in diverse pop. |
| Colorectal | Increased risk | Harmful | Replicate; test by subtype |
| BCC | Decreased risk | Protective | Replicate |
| Ovarian | No association | Null | Confirm null |
| Lung | No association | Null | Confirm null |
| Thyroid | No association | Null | Confirm null (RCT concerns) |
| Endometrial | Not tested | Unknown | **Novel** |
| HCC | Not tested | Unknown | **Novel** |
| Pancreatic | Not tested | Unknown | **Novel** |

---

## 10. Timeline and Deliverables

| Phase | Activity | Timeline |
|-------|----------|----------|
| 1 | All of Us workspace setup, IRB, DUA | Month 1-2 |
| 2 | Cohort definitions (ATLAS/Capr), instrument validation | Month 2-3 |
| 3 | Pilot: Breast cancer + BCC (highest power) | Month 3-4 |
| 4 | All primary outcomes + negative controls | Month 4-6 |
| 5 | Sensitivity analyses, subgroup analyses | Month 6-7 |
| 6 | Federated analysis with partner sites | Month 7-9 |
| 7 | Manuscript preparation, preprint | Month 9-10 |

### Key Deliverables

1. Medusa study package (reusable R code for any OMOP CDM site)
2. Manuscript for high-impact journal (e.g., JAMA Oncology, Nature Medicine)
3. Interactive Shiny dashboard (via `launchResultsExplorer()`)
4. Open-source code and cohort definitions for community replication

---

## 11. References

1. Sun H, et al. (2024). Association of GLP-1 receptor agonists with risk of
   cancers -- evidence from a drug-target MR and clinical trials. *International
   Journal of Surgery*.
2. Silverii GA, et al. (2025). GLP-1 receptor agonists and the risk for cancer:
   a meta-analysis of RCTs. *Diabetes, Obesity and Metabolism*.
3. All of Us Research Program (2024). Genomic data in the All of Us Research
   Program. *Nature*.
4. Zheng J, et al. (2020). Phenome-wide MR mapping the influence of the plasma
   proteome on complex diseases. *Nature Genetics*.
5. Schuemie MJ, et al. (2014). Interpreting observational studies: why empirical
   calibration is needed. *Statistics in Medicine*.
6. eQTLGen Consortium (2019). eQTLGen: large-scale cis- and trans-eQTL
   analyses identify thousands of genetic loci. *Nature Genetics*.
