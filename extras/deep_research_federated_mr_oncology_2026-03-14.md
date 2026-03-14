# Federated Mendelian Randomization for Oncology: Why Medusa's Architecture Is Uniquely Valuable

Date: 2026-03-14
Repository: `Medusa` (OHDSI)

**Note**: This synthesis draws on published literature through early 2025 and the Medusa codebase. Web search tools were unavailable during generation; citations reflect the author's knowledge of the field and should be verified against current databases.

---

## Executive Summary

Federated individual-level Mendelian randomization (MR), as implemented by Medusa, addresses five fundamental limitations that block conventional summary-statistic MR from reaching its potential in oncology:

1. **Sample size for rare cancers and subtypes** -- most cancer GWAS are insufficient for subtype-specific MR; federated pooling across biobank networks solves this without centralizing data.
2. **Multi-ancestry generalizability** -- health system biobanks are ancestrally diverse; individual-level data allows proper stratification, ancestry-PC adjustment, and trans-ethnic instrument validation.
3. **Nonlinear and subgroup analyses** -- individual-level data enables dose-response modeling, gene-environment interaction testing, and clinically meaningful subgroup MR that summary statistics cannot support.
4. **Covariate-rich adjustment** -- OMOP CDM provides structured clinical phenotypes (comorbidities, medications, labs) for population stratification control and negative-control analyses impossible with public GWAS summaries.
5. **Privacy-preserving scale** -- Medusa's one-shot profile likelihood pooling is mathematically lossless while sharing only ~600 numbers per site, making it governance-compatible across institutional boundaries.

---

## 1. Why Federated/Distributed MR Matters for Oncology

### 1.1 The Fundamental Problem: Cancer Is Statistically Sparse

Cancer incidence rates create a severe power problem for MR. Even common cancers like colorectal cancer (CRC, ~4-5% lifetime risk) yield only hundreds to low thousands of cases in a typical health system biobank of 50,000-200,000 genotyped individuals. For MR, where effect sizes are modest (typical ORs of 1.05-1.30 per genetically proxied unit change), statistical power demands are large.

The problem is dramatically worse for:

- **Rare cancers**: Pancreatic cancer (~1.5% lifetime risk), esophageal cancer (~0.5%), gallbladder cancer, mesothelioma, and most pediatric cancers. A single site with 100,000 genotyped patients may have only 50-200 incident cases of pancreatic cancer -- far too few for any MR analysis.
- **Cancer subtypes**: Triple-negative breast cancer (~15% of breast cancers), MSI-high colorectal cancer (~15% of CRC), small-cell lung cancer (~15% of lung cancers), clear cell renal cell carcinoma. Subtype-level analyses typically cut case counts by 3-10x.
- **Molecular subtypes**: HER2+ breast cancer, EGFR-mutant NSCLC, BRAF-mutant melanoma -- the categories that actually drive therapeutic decision-making.

### 1.2 Which Cancers Have Sufficient Cases Across Biobank Networks?

The OHDSI network, combined with institutions like the eMERGE network, UK Biobank, Million Veteran Program (MVP), All of Us, and BioVU, collectively represents millions of genotyped individuals. Realistic case count estimates across a 10-site federated network (assuming ~1M total genotyped individuals):

| Cancer type | Estimated cases (single 100K site) | Estimated cases (10-site, 1M total) | MR feasibility |
|---|---|---|---|
| Breast cancer (all) | 1,000-2,000 | 10,000-20,000 | Excellent |
| Colorectal cancer (all) | 400-800 | 4,000-8,000 | Good |
| Prostate cancer | 800-1,500 | 8,000-15,000 | Excellent |
| Lung cancer (all) | 500-1,000 | 5,000-10,000 | Good |
| Melanoma | 200-500 | 2,000-5,000 | Moderate |
| Pancreatic cancer | 50-150 | 500-1,500 | Marginal -> feasible |
| Bladder cancer | 150-300 | 1,500-3,000 | Moderate |
| Kidney cancer | 100-250 | 1,000-2,500 | Moderate |
| Triple-negative breast cancer | 150-300 | 1,500-3,000 | Moderate (federated only) |
| MSI-high CRC | 60-120 | 600-1,200 | Marginal -> feasible (federated only) |
| Esophageal adenocarcinoma | 30-80 | 300-800 | Marginal |
| Hepatocellular carcinoma | 50-150 | 500-1,500 | Marginal -> feasible |

**Key insight**: Federation moves several clinically important cancer types from "impossible" to "feasible" for MR. This is not merely an efficiency gain -- it enables entirely new research questions.

### 1.3 The Governance Advantage of One-Shot Pooling

Medusa's profile likelihood pooling (`poolLikelihoodProfiles()`) is specifically designed for the governance constraints that dominate multi-site genomic research:

- Each site exports only a numeric vector of ~600 log-likelihood values, case/control counts, and optional per-SNP summary statistics.
- No individual-level genotypes, phenotypes, or demographics leave any site.
- The pooled likelihood is mathematically equivalent to fitting the model on the combined data (under standard regularity conditions), making the approach lossless.
- This satisfies IRB, HIPAA, and data use agreement (DUA) requirements that typically prohibit sharing individual-level genomic data.

This is qualitatively different from meta-analysis of summary statistics: profile likelihood pooling preserves the full shape of the likelihood surface, not just point estimates and standard errors. This matters when likelihoods are skewed, multimodal, or when sites have very different case counts.

### 1.4 Comparison With Alternative Approaches

| Approach | Data shared | Statistical efficiency | Governance burden | Supports nonlinear MR |
|---|---|---|---|---|
| Centralized individual-level pooling | All patient data | Optimal | Very high (often infeasible) | Yes |
| **Medusa federated likelihood** | ~600 numbers per site | Lossless (equivalent to centralized) | Low | Possible with extensions |
| Fixed-effect meta-analysis of site estimates | Point estimate + SE per site | Efficient but loses likelihood shape | Low | No |
| Summary-statistic MR (e.g., from public GWAS) | Published GWAS summaries | Depends on GWAS size | None | No |

---

## 2. Multi-Ancestry MR Advantages

### 2.1 How Ancestral Diversity Strengthens MR Inference

Multi-ancestry MR provides several methodological advantages that are particularly relevant for oncology:

**Pleiotropy detection through cross-ancestry consistency**. If a genetic instrument affects a cancer outcome through the intended exposure pathway (e.g., IL-6 signaling), the causal estimate should be consistent across ancestries after accounting for allele frequency and LD differences. Inconsistency across ancestries is a signal of horizontal pleiotropy or population-specific confounding. This is one of the most powerful pleiotropy checks available.

**LD pattern differences resolve confounding by LD**. A major concern in cis-MR (common in drug-target validation) is that the instrumental SNP may be in LD with a nearby causal variant that affects the outcome through a different pathway. LD patterns differ across ancestries, so if the same causal estimate is obtained using instruments selected in European, East Asian, and African-ancestry populations, the probability that LD confounding explains the result drops substantially.

**Allele frequency variation provides natural dose-response information**. SNPs that are common in one ancestry but rare in another create natural variation in instrument strength, allowing implicit dose-response assessment across populations.

**Generalizability of drug targets**. A drug target validated by MR only in European populations may not generalize. Multi-ancestry MR is essential for assessing whether a target is relevant across the patient populations who will actually receive the therapy.

### 2.2 Cancer-Target Associations Studied in Non-European Populations

Published trans-ethnic or non-European MR studies in oncology include:

- **BMI and breast cancer**: Multiple MR studies across European and East Asian populations show consistent positive causal effects of adiposity on postmenopausal breast cancer, but the relationship is inverse for premenopausal breast cancer. East Asian GWAS (e.g., BCAC Asian component, BioBank Japan) provide independent replication.

- **Circulating testosterone and prostate cancer**: MR analyses using instruments from multi-ancestry GWAS of testosterone levels have been conducted. Results are largely consistent across European and African-ancestry populations, though instrument strength varies.

- **Height and cancer risk**: Genetically predicted height shows consistent positive associations with colorectal cancer, breast cancer, and melanoma risk across European and East Asian MR analyses. This is one of the most robustly replicated trans-ethnic MR findings in oncology.

- **Lipid fractions and cancer**: HDL-cholesterol and LDL-cholesterol MR studies for various cancers have been conducted using instruments from GLGC (multi-ethnic) and BioBank Japan. Results for LDL-breast cancer and HDL-ovarian cancer have been examined across ancestries.

- **Type 2 diabetes and cancer**: T2D instruments from DIAGRAM (European) and East Asian GWAS have been used for MR studies of pancreatic, liver, and endometrial cancer risk, with generally consistent findings across ancestries.

- **Inflammation markers (CRP, IL-6)**: CRP MR studies have been conducted across European and East Asian populations for lung and colorectal cancer, with some evidence of ancestry-specific differences.

### 2.3 Instrument Portability Across Ancestries

The portability of MR instruments across ancestries depends on the exposure being studied:

**Generally portable** (effect direction conserved, instruments often transferable):
- Lipid-pathway variants (HMGCR, PCSK9, NPC1L1 regions) -- well-validated across multiple ancestries
- IL6R rs2228145 -- functionally characterized missense variant, broadly portable
- BMI instruments from large GWAS (GIANT, multi-ethnic) -- partially portable but require re-weighting

**Partially portable** (some instruments transfer, others do not):
- Most polygenic instruments for complex traits -- LD pruning must be redone per ancestry, effect sizes need re-estimation
- eQTL-based instruments -- gene regulation is partially ancestry-specific due to regulatory variant frequency differences

**Poorly portable** (require ancestry-specific instrument construction):
- Instruments for traits where the genetic architecture differs substantially across ancestries
- Variants that are monomorphic or very rare in some populations (e.g., some European-common variants are absent in East Asian populations)

**Medusa's advantage here**: Because Medusa works with individual-level data at each site, it can incorporate ancestry principal components (the `buildMRCovariates()` function supports an `ancestryPCsTable` parameter with up to 10 PCs). This allows proper within-site ancestry stratification or adjustment, which is impossible with pre-computed GWAS summary statistics. Sites with diverse populations can contribute without introducing population stratification bias.

Li and Morrison (2025/2026) showed that even intra-continental ancestry mismatch between exposure and outcome samples causes attenuation bias in two-sample MR. Medusa's individual-level design, with per-site ancestry PC adjustment, directly mitigates this.

---

## 3. Individual-Level vs Summary-Level MR

### 3.1 Capabilities Enabled by Individual-Level Data

Access to individual-level data in OMOP CDM enables analyses that are fundamentally impossible with GWAS summary statistics alone:

#### 3.1.1 Nonlinear MR

Standard summary-statistic MR assumes a linear (log-linear for binary outcomes) relationship between exposure and outcome. Individual-level data allows:

- **Stratified MR / localized average causal response (LACR)**: Dividing individuals into quantiles of the genetic instrument (or exposure residualized on covariates) and estimating the causal effect within each stratum. This reveals U-shaped, threshold, or saturation effects. For oncology, this is critical: the dose-response relationship between a drug target and cancer risk may be nonlinear. For example, very low LDL-cholesterol might have different cancer implications than moderate reductions.

- **Fractional polynomial MR**: Methods like those of Staley and Burgess (2017) require individual-level data to estimate flexible dose-response curves.

- **Doubly-ranked method**: The approach of Tian et al. (2023) for nonlinear MR requires individual-level instrument and exposure data.

**Oncology relevance**: Drug development needs to know not just "does the target affect cancer risk?" but "at what dose/level of modulation does the effect become clinically meaningful?" Nonlinear MR with individual-level data can inform this.

#### 3.1.2 Subgroup Analyses

Individual-level data allows MR within clinically defined subgroups:

- **By sex**: Essential for hormone-related cancers (breast, prostate, endometrial). Effect modification by sex is a real biological phenomenon.
- **By age group**: Cancer biology changes with age; MR effects may differ between early-onset and late-onset disease.
- **By comorbidity status**: Does the causal effect of adiposity on CRC differ in patients with vs without type 2 diabetes? Individual-level data from OMOP CDM makes this testable.
- **By treatment exposure**: Does the effect of genetically proxied inflammation on cancer progression differ in patients who have vs have not received immunotherapy? This requires linked treatment data unavailable in GWAS summaries.
- **By cancer stage/grade**: With OMOP oncology extensions, MR analyses can be stratified by clinical features of the cancer itself.

#### 3.1.3 Gene-Environment Interaction MR

The most powerful application unique to individual-level MR:

- **GxE interaction testing**: Does the causal effect of a genetically proxied exposure differ by an environmental modifier? For example, does genetically proxied alcohol metabolism have a different effect on esophageal cancer risk in smokers vs non-smokers?
- **Factorial MR designs**: Combining genetic instruments for two exposures to test interaction effects on cancer outcomes. This can inform combination therapy strategies.
- **Time-varying exposures**: With longitudinal OMOP data, the relationship between genetically proxied exposure trajectories and cancer incidence can be modeled.

#### 3.1.4 Covariate Adjustment for Population Stratification

This is where Medusa's OMOP CDM foundation provides a decisive advantage:

- **Ancestry PC adjustment**: Medusa's `buildMRCovariates()` function incorporates ancestry PCs directly into the outcome model. Summary-statistic MR relies on the original GWAS having done this correctly -- and many older GWAS had suboptimal stratification control.
- **Clinical covariate adjustment**: OMOP CDM provides structured access to comorbidities, medications, laboratory values, and procedures. These can be included as covariates to reduce residual confounding in the outcome model.
- **Collider bias detection**: Individual-level data allows testing whether conditioning on covariates induces collider bias, which is invisible in summary-statistic analyses.

### 3.2 Summary: What Individual-Level Federated MR Adds

| Analysis type | Summary-stat MR | Individual-level MR (Medusa) |
|---|---|---|
| Linear causal effect estimation | Yes | Yes |
| Nonlinear dose-response | No | Yes |
| Sex-stratified MR | Only if sex-stratified GWAS exist | Yes (flexibly) |
| Age-stratified MR | Rarely available | Yes |
| Comorbidity-stratified MR | No | Yes |
| Gene-environment interactions | Very limited | Yes |
| Ancestry PC adjustment at analysis time | No (must trust source GWAS) | Yes |
| Covariate-adjusted outcome models | No | Yes |
| Time-to-event MR | No | Yes (with longitudinal OMOP data) |
| Negative control outcome testing | Limited | Yes (rich phenome available) |

---

## 4. The OHDSI Network for Oncology

### 4.1 OHDSI Sites With Genomic Data Linked to OMOP CDM

Several OHDSI-participating institutions have or are developing genomic data linked to their OMOP CDM instances:

- **Columbia University Irving Medical Center**: One of the founding OHDSI sites, with biobank genotyping linked to their OMOP CDM. Active in pharmacogenomics research.

- **Vanderbilt University Medical Center (BioVU)**: BioVU is one of the largest academic biobanks (~300,000+ genotyped individuals as of 2024), with data mapped to OMOP CDM. Strong oncology data through the Vanderbilt-Ingram Cancer Center. Has the VARIANT_OCCURRENCE table or equivalent genotype storage.

- **Mount Sinai (BioMe)**: The BioMe biobank (~50,000+ participants) is ancestrally diverse (large Hispanic, African-American representation) and linked to OMOP CDM. Mount Sinai is active in OHDSI.

- **University of Colorado (UCHealth)**: Has biobank data being mapped to OMOP. Part of the OHDSI network.

- **Veterans Affairs (MVP)**: The Million Veteran Program has ~900,000+ enrolled participants with genotyping, and VA data is being mapped to OMOP CDM through the VA OHDSI node. MVP has particular strength in prostate cancer, lung cancer, and liver cancer due to VA demographics.

- **All of Us Research Program**: NIH's All of Us has ~250,000+ participants with whole-genome sequencing, mapped to an OMOP-based data model. Highly diverse (>50% from underrepresented populations). Cancer cases are accumulating as the cohort matures.

- **UK Biobank**: While not a formal OHDSI member, UK Biobank data (~500,000 genotyped) has been mapped to OMOP CDM by several groups. Cancer ascertainment through linked cancer registry data is comprehensive.

- **eMERGE Network sites**: Several eMERGE sites (Northwestern, Geisinger, Marshfield, Partners/Mass General Brigham) participate in OHDSI and have genotyped biobanks. eMERGE III and IV have focused on clinical implementation of genomics.

- **Ajou University, Samsung Medical Center, and other Korean OHDSI sites**: Korean institutions active in OHDSI have large-scale genotyping efforts, providing East Asian ancestry representation.

- **IQVIA, Optum, and other data partners**: While these commercial partners have large OMOP CDM instances, genomic linkage is more limited. Some are exploring integration.

### 4.2 OMOP Oncology Extensions

The OHDSI community has developed several oncology-specific extensions to the OMOP CDM:

#### 4.2.1 OMOP Oncology CDM Module

Developed by the OHDSI Oncology Working Group, this extension adds structured representation for:

- **Cancer diagnosis**: Integration with tumor registries, including ICD-O-3 morphology and topography coding mapped to OMOP concepts.
- **Cancer staging**: TNM staging (AJCC/UICC) represented through the MEASUREMENT table with specific staging concepts. This allows MR analyses stratified by cancer stage at diagnosis.
- **Cancer-specific treatments**: Structured representation of chemotherapy regimens, radiation therapy episodes, and surgical procedures mapped to standard OMOP drug/procedure concepts.
- **Treatment response**: RECIST criteria, PSA response, pathologic response categories stored as MEASUREMENT or OBSERVATION records.

#### 4.2.2 OMOP Genomic CDM Extension

The extension that Medusa directly leverages (via the `VARIANT_OCCURRENCE` table):

- **VARIANT_OCCURRENCE**: Stores per-person variant calls with `person_id`, `rs_id`, `genotype`, `reference_allele`, `alternate_allele`. This is what Medusa's `buildMRCohort()` queries.
- **STRUCTURAL_VARIANT_OCCURRENCE**: For copy number variants and structural variants relevant to cancer genomics.
- **SEQUENCE_OCCURRENCE**: For gene expression or other sequence-level data.

The Genomic CDM was designed with pharmacogenomics and genetic epidemiology in mind. Its integration with the clinical OMOP CDM is what makes Medusa-style federated MR possible.

#### 4.2.3 Tumor Registry Data in OMOP

Several OHDSI sites have mapped tumor registry data (e.g., SEER-linked data, institutional cancer registries) to OMOP:

- Cancer-specific diagnosis dates with histology information
- Stage at diagnosis (localized, regional, distant)
- First course of treatment
- Vital status and cause of death

This data, when combined with genomic data, enables cancer-specific MR analyses that account for tumor characteristics.

### 4.3 What This Means for Medusa

The OHDSI oncology infrastructure means Medusa can leverage:

1. **Standardized cancer phenotyping**: Outcome cohorts defined using OMOP concept IDs are portable across all OHDSI sites. A colorectal cancer cohort defined in ATLAS at one site can be instantiated identically at every other site.

2. **Rich covariate space**: Through `buildMRCovariates()` and FeatureExtraction, Medusa accesses the full clinical phenome for each patient -- not just cancer status, but all conditions, medications, lab values, and procedures. This enables the negative-control testing and covariate-adjusted MR that distinguishes individual-level from summary-statistic approaches.

3. **Cancer subtype resolution**: With oncology extensions, outcome cohorts can be defined at the subtype level (e.g., ER+/HER2- breast cancer, MSI-high CRC) using structured OMOP concepts rather than just ICD codes.

4. **Treatment-aware analyses**: MR analyses can be conditioned on or stratified by prior treatment exposure -- critical for understanding whether a genetic proxy for a drug target predicts outcomes in treatment-naive vs treated patients.

---

## 5. Specific Oncology Applications Where Federated Individual-Level MR Is Uniquely Valuable

### 5.1 Drug Target Validation for Cancer

The strongest near-term application. Minikel et al. (2024) showed that drug mechanisms with human genetic support have ~2.6x higher probability of clinical success. For oncology targets, federated MR via Medusa enables:

- **Validating targets for rare cancers** that no single biobank can study (e.g., is genetically proxied VEGF-A causally related to renal cell carcinoma risk? Each site may have only 100-200 RCC cases, but 10 sites together have 1,000-2,000).
- **Subtype-specific target validation** (e.g., does genetically proxied PD-L1 expression have different effects on MSI-high vs MSI-stable CRC? Requires subtype-level case counts only achievable through federation).
- **Multi-ancestry target validation** to assess generalizability before clinical trials enroll diverse populations.

### 5.2 Cancer Chemoprevention Candidate Identification

MR can identify modifiable risk factors for cancer prevention. Individual-level federated MR extends this by:

- **Nonlinear dose-response for chemopreventive agents**: Does the protective effect of genetically proxied aspirin-like COX-2 inhibition on CRC show a threshold effect?
- **Subgroup identification for chemoprevention**: Does genetically proxied statin exposure reduce cancer risk equally across sexes, age groups, and comorbidity profiles?
- **Time-varying effects**: Using longitudinal OMOP data, does the "MR-predicted" chemopreventive effect depend on duration of exposure?

### 5.3 Understanding Cancer Health Disparities

Multi-ancestry federated MR is uniquely positioned to address disparities:

- **Do known risk factors have different causal magnitudes across ancestries?** For example, does genetically proxied obesity have a stronger causal effect on endometrial cancer in African-American vs European-American women?
- **Are therapeutic targets equally relevant across populations?** Critical for ensuring equitable benefit from targeted therapies.
- **Population-specific confounding**: Individual-level ancestry PC adjustment, site by site, prevents spurious cross-ancestry associations that plague summary-statistic meta-analyses.

### 5.4 Immune-Oncology Target Discovery

The immune-oncology revolution creates new MR opportunities:

- **Cytokine-cancer relationships**: Genetically proxied levels of IL-6, TNF, IL-10, interferon-gamma, and dozens of other cytokines can be tested against cancer incidence. The IL-6/CRC example in Medusa's vignette is a prototype for this.
- **Immune cell composition**: Using instruments from blood cell trait GWAS to test whether genetically determined immune cell proportions (e.g., CD8+ T cell counts) causally affect cancer risk or progression.
- **Checkpoint molecule expression**: cis-MR using instruments near PD-1, PD-L1, CTLA-4, LAG-3, and TIGIT genes to test whether genetically proxied expression levels predict cancer outcomes.

### 5.5 Cancer Survivorship and Treatment Response

A more speculative but high-value application:

- **Pharmacogenomic MR**: Does genetically proxied drug metabolism (e.g., CYP2D6 activity for tamoxifen) causally affect breast cancer recurrence? This requires individual-level treatment and outcome data.
- **Adverse event prediction**: MR for genetically proxied drug exposure and specific treatment toxicities.
- **Second primary cancer risk**: Does genetically proxied exposure to a risk factor predict second cancer risk in cancer survivors?

---

## 6. Methodological Considerations Specific to Oncology MR

### 6.1 Challenges

- **Competing risks**: Cancer patients face competing mortality, which complicates MR interpretation. Individual-level data allows competing-risk models; summary statistics do not.
- **Survival bias / incidence-prevalence bias**: Using prevalent rather than incident cancer cases can bias MR. OMOP CDM's longitudinal design with observation periods helps identify truly incident cases.
- **Reverse causation windows**: Cancer may take decades to develop, but genetic instruments act from conception. The long latency of most cancers actually strengthens the MR design against reverse causation compared to acute diseases.
- **Tumor heterogeneity**: The same ICD code can represent biologically distinct diseases. OMOP oncology extensions with histology and molecular subtyping improve phenotype precision.

### 6.2 Medusa-Specific Strengths for Oncology

1. **Profile likelihood preserves information**: For rare cancers with small case counts at individual sites, the full likelihood profile captures uncertainty better than point estimates. Pooling profiles is more informative than meta-analyzing point estimates when site-level likelihoods are skewed or poorly approximated by normal distributions.

2. **One-shot communication**: Oncology data is among the most sensitive clinical data. Medusa's one-shot design (no iterative communication) minimizes governance friction and attack surface compared to iterative federated learning approaches.

3. **OMOP CDM standardization**: Cancer cohort definitions are portable. A coordinator defines the outcome cohort (e.g., incident colorectal cancer with specific SNOMED concepts, washout period, exclusion criteria), and every OHDSI site can instantiate it identically using ATLAS or the same SQL.

4. **FeatureExtraction integration**: Medusa's `buildMRCovariates()` leverages the OHDSI FeatureExtraction package to automatically assemble a rich covariate matrix from the full OMOP phenome. This provides the covariate adjustment that makes individual-level MR superior to summary-statistic approaches.

---

## 7. Key References

### Federated and Distributed Analysis Methods
- Luo et al. (2022). dPQL: lossless distributed GLMM. *Biostatistics*. PMID: 35579348
- Zhang et al. (2025). COLA-GLM-H: one-shot lossless distributed GLM. PMID: 41259033
- Duan et al. (2020). ODAL: one-shot distributed algorithm for logistic regression. *Journal of Biomedical Informatics*
- Li et al. (2023). Federated causal inference. Various distributed causal inference methodologies

### Individual-Level MR Methods
- Burgess, Butterworth, Thompson (2013). Multi-variant summarized-data MR (IVW framework). *Genetic Epidemiology*
- Staley and Burgess (2017). Semiparametric MR for nonlinear effects. *Int J Epidemiol*
- Tian et al. (2023). Doubly-ranked method for nonlinear MR
- Spiller et al. (2019). Detecting and accommodating GxE interactions in MR

### Multi-Ancestry MR
- Li, Morrison (2025/2026). Population mismatch bias in two-sample MR. *AJHG*
- Burgess et al. (2022). Guidelines for MR across ancestry groups
- Mahajan et al. (2022). Multi-ancestry GWAS of type 2 diabetes -- instruments used in trans-ethnic MR
- Chen et al. (2023). Multi-ancestry Mendelian randomization studies across various cancer types

### Cancer MR Applications
- Cornish et al. (2020). Modifiable risk factors for cancer identified by MR. *JNCI*
- Dimitrakopoulou et al. (2022). Circulating interleukins and cancer risk (MR)
- Larsson et al. (2020). Body mass index and cancer risk: comprehensive MR study
- Bull et al. (2020). Adiposity and cancer risk: systematic review of MR studies

### OHDSI Oncology
- Belenkaya et al. (2021). OMOP Oncology CDM extension. *JCO Clinical Cancer Informatics*
- Warner et al. (2023). Oncology module for OMOP CDM
- OHDSI Genomic CDM: https://github.com/OHDSI/Genomic-CDM
- Hripcsak et al. (2015). OHDSI collaborative approach. *JAMIA*

### Drug Target Validation
- Minikel et al. (2024). Genetic evidence and drug development success (~2.6x improvement). *Nature*. PMID: 38632401
- Schmidt et al. (2020). Genetic drug target validation framework. PMID: 32591531
- Zheng et al. (2020). Plasma proteome MR + colocalization. *Nature Genetics*

---

## 8. Conclusion: The Case for Medusa in Oncology

Federated individual-level MR, as implemented by Medusa, is not merely an incremental improvement over summary-statistic MR for oncology applications. It enables qualitatively different analyses:

1. **It makes rare cancer and subtype MR possible** by aggregating cases across sites while preserving privacy.
2. **It makes multi-ancestry MR rigorous** by allowing per-site ancestry adjustment rather than trusting external GWAS stratification.
3. **It enables nonlinear, subgroup, and interaction analyses** that are the difference between "this target affects cancer risk" and "this target affects cancer risk in this patient population at this level of modulation."
4. **It leverages the full clinical richness of OMOP CDM** for covariate adjustment, negative controls, and cancer-specific phenotyping that no public GWAS can match.
5. **It does all of this in a governance-compatible one-shot design** that requires sharing only aggregate statistical summaries.

The closest analogues in the literature are DataSHIELD-based federated analyses and consortium meta-analyses of individual participant data. Medusa's profile likelihood approach is more statistically efficient than DataSHIELD's iterative distributed regression (which requires multiple communication rounds and is lossy), and more privacy-preserving than traditional IPD meta-analysis (which requires data pooling). For the specific application of two-sample MR in oncology, Medusa occupies a unique and valuable methodological niche.
