# Medusa Oncology Research Program
## Two-Sample Mendelian Randomization for Drug Target Validation

*Generated 2026-03-14 | Medusa (Federated MR on OMOP CDM)*

---

## Executive Summary

This document outlines a program of oncology-focused Mendelian Randomization
(MR) studies using Medusa on the All of Us dataset and other OMOP CDM databases
with linked genomic data. The studies are prioritized by (1) clinical impact
potential, (2) feasibility given available data, and (3) novelty relative to
existing literature. We focus on hypotheses where individual-level, federated MR
on diverse populations adds value beyond what published summary-level MR has
already established.

---

## Data Landscape

### All of Us Research Program

- **535,000 short-read whole genome sequences** (spring 2026 release), up from 414,000
- **77% from historically underrepresented communities**, 46% racial/ethnic minorities
- EHR data harmonized to OMOP CDM from ~50 health care organizations
- Upcoming linkages to **cancer registries**, mortality data, and claims
- New proteomic (10,000 participants) and metabolomic (5,500) datasets in 2026
- Compute via the Researcher Workbench (Jupyter/R in a controlled environment)

### Complementary Datasets

| Dataset | N (genomic) | Strengths |
|---------|-------------|-----------|
| UK Biobank | ~500,000 | Mature, deeply phenotyped, large cancer cohort |
| Million Veteran Program (MVP) | ~900,000 | Large, diverse, VA health system |
| eMERGE Network | ~100,000 | Multi-site, linked EHR, OMOP conversion ongoing |
| BioVU (Vanderbilt) | ~300,000 | Longstanding biobank, strong oncology registry |
| GERA (Kaiser) | ~100,000 | Integrated health system, cancer registry linked |

### Medusa's Unique Advantages for This Program

1. **Federated pooling**: Combine All of Us + MVP + eMERGE without sharing patient-level data
2. **Multi-ancestry instruments**: Diverse All of Us population enables trans-ethnic MR, strengthening causal inference
3. **Individual-level data**: Enables nonlinear MR, subgroup analyses, gene-environment interaction
4. **OMOP CDM**: Standardized phenotyping across sites; reusable cohort definitions
5. **Rich covariate adjustment**: OMOP's covariate richness enables better instrument validation (PheWAS, covariate balance)
6. **Negative control calibration**: Built-in empirical validation guards against systematic bias

---

## Tier 1: High-Impact, High-Feasibility Studies

These studies target well-powered questions where Medusa adds clear value over
existing summary-level MR and where clinical translation is plausible.

### 1.1 GLP-1 Receptor Agonism and Cancer Risk

**Hypothesis**: Genetically proxied GLP-1 receptor activation reduces risk of
obesity-related cancers (breast, colorectal, endometrial, hepatocellular,
pancreatic) through metabolic pathways independent of weight loss.

**Background**: GLP-1 receptor agonists (semaglutide, tirzepatide) are the
fastest-growing drug class globally. A drug-target MR study (Sun et al., 2024)
found GLP-1RA may decrease breast cancer and BCC risk but increase colorectal
cancer risk. RCT meta-analyses show neutral-to-protective overall cancer effects.
The conflicting colorectal finding needs replication in diverse populations.

**Why Medusa adds value**:
- All of Us diversity enables trans-ethnic replication of the colorectal cancer
  signal, which was derived from European GWAS
- Individual-level data allows testing whether the effect is mediated by BMI
  reduction (mediation MR) or through direct GLP-1R signaling
- Federated analysis with MVP provides power for rarer cancers (pancreatic,
  hepatocellular)
- Can examine cancer subtypes (ER+/ER- breast, MSI/MSS colorectal) using OMOP
  phenotyping

**Instruments**: cis-eQTLs near *GLP1R* (7 SNPs from deCODE/UKB pQTL data)

**Outcomes**: Breast, colorectal, endometrial, hepatocellular, pancreatic cancer

**Expected sample sizes in All of Us** (estimated from US prevalence):
- Breast cancer: ~8,000-12,000 cases
- Colorectal cancer: ~3,000-5,000 cases
- Pancreatic cancer: ~500-1,000 cases

**Clinical impact**: Could inform or accelerate cancer prevention trials of
GLP-1 agonists. Given the massive at-risk population already on semaglutide,
even small absolute risk reductions would be transformative.

---

### 1.2 IL-6 Pathway Inhibition and Cancer Risk Across Ancestries

**Hypothesis**: Genetically proxied IL-6 receptor inhibition (mimicking
tocilizumab/sarilumab) reduces risk of endometrial cancer and potentially
other inflammation-driven cancers.

**Background**: The comprehensive Yarmolinsky et al. (2024) MR study in
eBioMedicine found that most inflammatory markers showed no causal association
with cancer, challenging conventional wisdom. However, IL-6 showed a positive
association with endometrial cancer, and IL-1RL1 showed a protective association
with triple-negative breast cancer. Prior work has validated IL-6R inhibition
for multiple diseases via drug-target MR.

**Why Medusa adds value**:
- IL-6R MR has been done almost exclusively in European populations; All of Us
  enables the first large-scale trans-ethnic replication
- Individual-level data allows testing IL-6R × BMI interaction (obesity drives
  both IL-6 and endometrial cancer risk)
- Can use negative controls to calibrate estimates and rule out pleiotropy
- OMOP phenotyping captures cancer subtypes and histology

**Instruments**: cis-pQTLs near *IL6R* (rs2228145 and flanking variants)

**Outcomes**: Endometrial, colorectal, breast (overall and TNBC), lung,
ovarian, prostate cancer

**Clinical impact**: Tocilizumab is already FDA-approved for RA and COVID-19.
If MR supports a causal role in endometrial cancer prevention, this is a
directly actionable drug repurposing opportunity.

---

### 1.3 PCSK9 Inhibition, LDL Cholesterol, and Gastrointestinal Cancers

**Hypothesis**: Genetically proxied PCSK9 inhibition increases hepatocellular
carcinoma risk through LDL receptor-mediated hepatocyte proliferation, while
lowering breast and lung cancer risk.

**Background**: MR studies show conflicting signals. One study found PCSK9
inhibition associated with increased hepatic cancer risk (OR = 1.99) and
decreased breast/lung cancer risk. Preclinical work shows PCSK9 directly
promotes KRAS-mutant colorectal cancer via cholesterol biosynthesis pathways.
PCSK9 inhibitors (evolocumab, alirocumab) are widely prescribed.

**Why Medusa adds value**:
- The hepatocellular carcinoma signal is critical to resolve given widespread
  PCSK9 inhibitor use; All of Us diverse population provides independent
  replication
- Individual-level MR can test whether the liver cancer signal is driven by
  underlying liver disease (confounding) or is a direct PCSK9 effect
- Can separate PCSK9 pathway from HMGCR (statin) pathway effects on same
  cancers using multivariable MR at the individual level
- Federated analysis across sites increases HCC case count (rare in single sites)

**Instruments**: cis-variants near *PCSK9* (rs11591147 and LD proxies)

**Outcomes**: Hepatocellular carcinoma, colorectal cancer (KRAS-mutant subtype),
breast cancer, lung cancer

**Clinical impact**: Safety signal monitoring for one of the most widely
prescribed new drug classes. Could inform FDA post-marketing surveillance.

---

### 1.4 Insulin/IGF-1 Axis and Pancreatic Cancer Prevention

**Hypothesis**: Genetically proxied fasting insulin reduction (mimicking
insulin-sensitizing drugs) causally reduces pancreatic cancer risk.

**Background**: MR evidence supports fasting insulin as a causal risk factor
for pancreatic cancer (OR = 1.66 per SD increase). The insulin → IGF-1 →
PI3K/AKT/mTOR pathway is well-characterized in pancreatic tumorigenesis.
Metformin was hypothesized as protective, but drug-target MR found it may
*increase* prostate cancer risk (OR = 1.55).

**Why Medusa adds value**:
- Pancreatic cancer is rare (~500-1,000 cases in All of Us); federated pooling
  across All of Us + MVP + eMERGE is essential for power
- Individual-level data enables stratification by diabetes status (pre-diabetic
  vs. diabetic vs. normal glycemia)
- Can decompose the pathway: instrument BMI, fasting insulin, IGF-1, and
  IGFBP-3 separately to identify the specific causal mediator
- All of Us proteomic data (2026 release) enables direct measurement of
  circulating IGF-1 for validation

**Instruments**: Variants near *INSR*, *IGF1R*, *IRS1*, *FTO*, *MC4R*

**Outcomes**: Pancreatic cancer, with secondary outcomes of colorectal, breast,
endometrial cancer

**Clinical impact**: Could redirect the metformin-cancer prevention field toward
more specific insulin-sensitizing targets (e.g., PPAR-γ agonists, GLP-1 agonists).

---

## Tier 2: Novel Hypotheses with High Potential

These studies explore less-established territory where positive results would
be particularly impactful.

### 2.1 Next-Generation Immune Checkpoint Targets

**Hypothesis**: Genetically proxied inhibition of LAG-3, TIGIT, or CD47 is
causally associated with altered cancer risk, providing human genetic validation
(or refutation) of these immuno-oncology targets.

**Background**: LAG-3 inhibition (relatlimab) is FDA-approved for melanoma.
TIGIT trials have shown mixed results. CD47 (magrolimab) trials have failed.
**No published MR study has evaluated these targets using drug-target MR** ---
this is a clear gap in the literature.

**Why Medusa adds value**:
- This has never been done; first-mover advantage
- Individual-level MR allows testing of checkpoint expression × tumor
  microenvironment interactions
- All of Us WGS provides full coverage of cis-regulatory variants near
  *LAG3*, *TIGIT*, *HAVCR2* (TIM-3), and *CD47*
- OMOP phenotyping can capture specific cancer types where these targets are
  being tested in clinical trials

**Instruments**: cis-eQTLs/pQTLs near *LAG3* (12p13.31), *TIGIT* (3q13.31),
*HAVCR2* (5q33.3), *CD47* (3q13.12)

**Outcomes**: Melanoma, NSCLC, AML, colorectal cancer, breast cancer

**Clinical impact**: Could validate or refute billion-dollar drug development
programs. If MR shows no causal link between genetically proxied TIGIT
inhibition and cancer outcomes, it provides a genetic explanation for mixed
clinical trial results.

---

### 2.2 Complement System and Colorectal Cancer

**Hypothesis**: Genetically proxied complement C1QB levels causally influence
colorectal cancer risk, suggesting complement modulation as a therapeutic
strategy.

**Background**: A recent drug-target MR study (Frontiers in Genetics, 2024)
identified complement C1QB as a novel causal factor for colorectal cancer.
Complement activation in the tumor microenvironment is increasingly recognized
as a driver of cancer progression. Multiple complement-targeting drugs are in
development for other indications (e.g., eculizumab for PNH).

**Why Medusa adds value**:
- Replication in diverse ancestry populations
- Individual-level data allows testing whether complement effects differ by
  MSI status (microsatellite instability affects immune microenvironment)
- Can use PheWAS to validate instrument specificity for complement pathway

**Instruments**: cis-pQTLs near *C1QB*, *C3*, *C5*, *CFH*

**Outcomes**: Colorectal cancer (all, MSI-H, MSS), with exploratory analyses
of gastric, ovarian, bladder cancer

---

### 2.3 CHIP/Clonal Hematopoiesis and Therapy-Related Cancer

**Hypothesis**: Genetically determined susceptibility to clonal hematopoiesis
(via telomere length, DNA repair variants) causally increases risk of
therapy-related myeloid neoplasms (t-MN) following cancer treatment.

**Background**: CHIP occurs in >10% of adults over 70 and dramatically
increases risk of therapy-related myeloid neoplasms. Mendelian randomization
has shown bidirectional causality between telomere length and CHIP. The
clinical question is whether baseline germline risk can predict who will
develop t-MN after chemotherapy.

**Why Medusa adds value**:
- Unique capability: Medusa can analyze patients who *received cancer treatment*
  and track subsequent t-MN as an outcome, using OMOP drug exposure data to
  identify the at-risk cohort
- All of Us WGS provides direct measurement of germline variants affecting
  DNA repair (*TP53*, *PPM1D*, *DNMT3A*, *TET2*)
- Federated analysis essential: t-MN is very rare (1-3% of treated patients)
- Individual-level data enables stratification by treatment type (alkylating
  agents, PARP inhibitors, platinum-based)

**Instruments**: Telomere length PRS, variants near *TERT*, *DNMT3A*, *TET2*

**Outcomes**: Therapy-related MDS/AML, secondary leukemia, second primary
malignancy

**Clinical impact**: Could identify patients who should avoid specific
chemotherapy regimens or receive enhanced monitoring. Directly translatable
to precision oncology decision-making.

---

### 2.4 Proteome-Wide MR Scan for Druggable Cancer Targets

**Hypothesis**: Systematic screen of ~3,000 circulating proteins using pQTL
instruments to identify novel druggable targets for the 10 most common cancers.

**Background**: Multiple proteome-wide MR studies have been published (Zheng
et al. 2023, lung cancer 2025, lymphoma 2025, gastric cancer 2024) but these
all used European-only GWAS summary statistics. None has leveraged
individual-level data from a diverse population.

**Why Medusa adds value**:
- All of Us 2026 proteomic data on 10,000 participants enables validation of
  pQTL instruments directly in the study population
- Trans-ethnic MR dramatically reduces false positives from population
  stratification
- Can apply negative control calibration (Medusa's empirical calibration
  pipeline) to the full proteome-wide scan, a methodological innovation not
  used in prior proteome-wide cancer MR studies
- Individual-level data enables colocalization-equivalent analyses within the
  MR framework

**Instruments**: cis-pQTLs from deCODE (n=35,559), UKB-PPP (n=54,306)

**Outcomes**: Breast, colorectal, lung, prostate, endometrial, ovarian,
pancreatic, hepatocellular, renal, melanoma

**Clinical impact**: Systematic target discovery pipeline; results feed directly
into pharma drug development prioritization.

---

## Tier 3: Speculative but Potentially Transformative

### 3.1 Gut Microbiome Metabolites and Colorectal Cancer

**Hypothesis**: Genetically determined levels of microbial metabolites (TMAO,
secondary bile acids, short-chain fatty acids) causally influence colorectal
cancer risk.

**Instruments**: Host genetic variants influencing circulating TMAO (*FMO3*),
bile acid composition (*CYP7A1*, *NR1H4*), SCFA levels

**Why speculative**: Instrument strength for microbiome metabolites is often
weak; biology is complex. But if valid instruments exist, this connects the
massive microbiome-cancer literature to causal inference for the first time.

---

### 3.2 Aspirin/COX-2 Pathway and Pan-Cancer Prevention

**Hypothesis**: Genetically proxied COX-2 inhibition (via *PTGS2* eQTLs)
reduces risk of colorectal, breast, and melanoma through prostaglandin-mediated
inflammation.

**Background**: Observational and RCT evidence supports aspirin for colorectal
cancer prevention. However, the mechanism (COX-dependent vs. COX-independent)
is debated. A drug-target MR approach can isolate the COX-2-specific signal.

**Instruments**: cis-eQTLs near *PTGS2*

**Outcomes**: Colorectal, breast, melanoma, ovarian cancer

---

### 3.3 Sex Hormone Pathways and Hormone-Sensitive Cancers

**Hypothesis**: Genetically proxied aromatase inhibition (via *CYP19A1*
variants) has differential cancer effects across ancestries, explaining
disparities in breast and endometrial cancer incidence.

**Why All of Us matters**: Estrogen metabolism genetics differ substantially
across ancestries; European-only MR may miss or mischaracterize effects in
African American and Hispanic populations.

---

## Implementation Roadmap

### Phase 1: Foundation (Months 1-3)

| Study | Target | Primary Cancer | Est. Cases |
|-------|--------|---------------|------------|
| 1.1 | GLP-1R | Breast, Colorectal | 8,000+ |
| 1.2 | IL-6R | Endometrial | 2,000+ |
| 1.3 | PCSK9 | Hepatocellular | 500+ |

**Deliverables**:
- All of Us workspace setup and IRB approval
- Cohort definitions for 10 primary cancer phenotypes in OMOP (reuse OHDSI
  community phenotype library where available)
- Instrument validation: F-statistics, PheWAS, allele frequency comparison
  in All of Us vs. GWAS source populations
- Negative control outcome selection (50 outcomes per study)
- Pilot results for GLP-1R and breast cancer

### Phase 2: Core Studies (Months 4-9)

| Study | Target | Primary Cancer | Est. Cases |
|-------|--------|---------------|------------|
| 1.4 | Insulin/IGF-1 | Pancreatic | 500+ |
| 2.1 | LAG3/TIGIT/CD47 | Melanoma, NSCLC | 3,000+ |
| 2.2 | Complement | Colorectal | 4,000+ |

**Deliverables**:
- Tier 1 studies completed, manuscripts in preparation
- Federated analysis protocol established (Medusa package deployed to MVP or
  eMERGE partner site)
- Tier 2 studies initiated with preliminary results

### Phase 3: Discovery (Months 10-15)

| Study | Target | Primary Cancer | Est. Cases |
|-------|--------|---------------|------------|
| 2.3 | CHIP/Telomeres | Therapy-related MN | 200+ |
| 2.4 | Proteome scan | Pan-cancer | Varies |
| 3.1-3.3 | Exploratory | Various | Varies |

**Deliverables**:
- Proteome-wide scan with empirical calibration (methodological innovation)
- Complete program results for 8+ studies
- OHDSI community study packages for replication at other sites

### Phase 4: Translation (Months 16-18)

- Synthesize results across all studies into a target prioritization framework
- Identify top 3-5 candidates for clinical trial prioritization
- Publish program overview paper in high-impact journal
- Open-source all analysis code and Medusa study packages

---

## Methodological Innovations Enabled by Medusa

1. **Empirically calibrated proteome-wide MR**: First application of OHDSI-style
   negative control calibration to a proteome-wide MR drug target scan
2. **Trans-ethnic federated MR**: Pooled individual-level MR across ancestrally
   diverse populations without sharing patient data
3. **Treatment-stratified MR**: Using OMOP drug exposure tables to stratify MR
   analyses by treatment history (e.g., CHIP → t-MN conditional on chemo type)
4. **Cancer subtype MR**: Leveraging OMOP phenotyping to define molecular
   subtypes (ER+/ER-, MSI-H/MSS, KRAS-mutant) as distinct outcomes
5. **Instrument validation via PheWAS**: Using OMOP covariate richness to
   comprehensively test instrument exclusion restriction

---

## Key References

- Yarmolinsky J, et al. (2024). Association between circulating inflammatory markers and adult cancer risk: a Mendelian randomization analysis. *eBioMedicine*.
- Sun H, et al. (2024). Association of glucagon-like peptide-1 receptor agonists with risk of cancers. *International Journal of Surgery*.
- Zheng J, et al. (2023). Proteome-wide Mendelian randomization implicates therapeutic targets in common cancers. *Journal of Translational Medicine*.
- Feng J, et al. (2024). Identifying genetically-supported drug repurposing targets for NSCLC through MR of the druggable genome. *Translational Lung Cancer Research*.
- Tsilidis KK, et al. (2022). Circulating inflammatory cytokines and risk of five cancers. *BMC Medicine*.
- All of Us Research Program (2024). Genomic data in the All of Us Research Program. *Nature*.
- British Journal of Cancer (2024). Blood lipids and LDL cholesterol lowering drug-targets with colorectal cancer risk.
- MDPI Genes (2024). PCSK9 inhibitors and malignant tumors: a Mendelian randomization study.
- Silverii GA, et al. (2025). GLP-1 receptor agonists and the risk for cancer: a meta-analysis of RCTs. *Diabetes, Obesity and Metabolism*.
