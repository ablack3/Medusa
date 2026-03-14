# Deep Research Review: Two-sample Mendelian Randomization in Observational Databases for Drug Development

Date: 2026-03-11
Repository: `Medusa` (OHDSI)

## 1. Scope and question
This review addresses:
1. The research landscape for two-sample Mendelian randomization (MR), especially for drug-target discovery/validation.
2. Whether the methods currently implemented in Medusa are methodologically correct.
3. What should be improved to align with current (through 2026) best practice for drug-development-grade MR.

## 2. Executive verdict
### Short answer
- **Core Medusa estimator is statistically correct** for a federated two-sample MR setup under its stated assumptions.
- **Sensitivity methods implemented are correct and standard** (IVW, MR-Egger, weighted median, Steiger, leave-one-out).
- **For drug-target MR, Medusa is not yet state-of-the-art** because key modern components are missing: robust cis-MR with correlated variants, colocalization-native workflows, overlap/winner's-curse corrections, and richer pleiotropy-robust estimators.

### Bottom-line assessment
- **Are the current methods correct?** Yes, for baseline inference.
- **Are they sufficient for high-confidence drug target prioritization?** Not yet.

## 3. What Medusa currently does (code audit)
The package currently implements:
- Federated pooling by summing profile log-likelihoods across sites (one-shot communication) in [`/Users/ablack/ohdsi/Medusa/R/poolLikelihoodProfiles.R:17`](/Users/ablack/ohdsi/Medusa/R/poolLikelihoodProfiles.R:17).
- Site-level outcome modeling via exact profile likelihood (offset-based profiling in logistic regression) in [`/Users/ablack/ohdsi/Medusa/R/fitOutcomeModel.R:283`](/Users/ablack/ohdsi/Medusa/R/fitOutcomeModel.R:283).
- MR effect computed as Wald ratio using the exact fitted score denominator and delta-method SE in [`/Users/ablack/ohdsi/Medusa/R/computeMREstimate.R:150`](/Users/ablack/ohdsi/Medusa/R/computeMREstimate.R:150).
- Sensitivity analyses (IVW, MR-Egger, weighted median, Steiger, leave-one-out; optional TwoSampleMR delegation) in [`/Users/ablack/ohdsi/Medusa/R/runSensitivityAnalyses.R:17`](/Users/ablack/ohdsi/Medusa/R/runSensitivityAnalyses.R:17).
- Instrument assembly from OpenGWAS + LD clumping in [`/Users/ablack/ohdsi/Medusa/R/getMRInstruments.R:168`](/Users/ablack/ohdsi/Medusa/R/getMRInstruments.R:168).
- Allele harmonization with palindromic handling + optional cohort AF disambiguation in [`/Users/ablack/ohdsi/Medusa/R/harmonizeAlleles.R:163`](/Users/ablack/ohdsi/Medusa/R/harmonizeAlleles.R:163).
- Diagnostics (F-statistics, PheWAS, AF discrepancy, missingness) in [`/Users/ablack/ohdsi/Medusa/R/runInstrumentDiagnostics.R:17`](/Users/ablack/ohdsi/Medusa/R/runInstrumentDiagnostics.R:17).
- Cross-checking against TwoSampleMR in [`/Users/ablack/ohdsi/Medusa/R/validateAgainstTwoSampleMR.R:17`](/Users/ablack/ohdsi/Medusa/R/validateAgainstTwoSampleMR.R:17).

Current gaps visible directly in code:
- Negative-control testing is a stub: [`/Users/ablack/ohdsi/Medusa/R/runInstrumentDiagnostics.R:419`](/Users/ablack/ohdsi/Medusa/R/runInstrumentDiagnostics.R:419).
- Internal Steiger unavailable for binary outcomes: [`/Users/ablack/ohdsi/Medusa/R/runSensitivityAnalyses.R:796`](/Users/ablack/ohdsi/Medusa/R/runSensitivityAnalyses.R:796).
- No built-in colocalization/fine-mapping module.
- No explicit sample-overlap/winner's-curse correction.
- No robust cis-MR engine for correlated cis-variants beyond clumping + standard estimators.

## 4. Annotated bibliography (drug-target and two-sample MR focused)

### Foundations and core estimators
1. **Burgess, Butterworth, Thompson (2013)**. *Mendelian randomization analysis with multiple genetic variants using summarized data*.
   - Link: https://pmc.ncbi.nlm.nih.gov/articles/PMC4377079/
   - Contribution: Formalized multi-variant summarized-data MR (IVW framework).
   - Why it matters here: Medusa’s IVW logic and weighted score framing are aligned with this foundation.

2. **Bowden, Davey Smith, Burgess (2015)**. *MR-Egger regression with invalid instruments*.
   - Link: https://pubmed.ncbi.nlm.nih.gov/26050253/
   - Contribution: Intercept-based pleiotropy detection and slope estimate under InSIDE.
   - Relevance: Medusa’s MR-Egger implementation matches this class of sensitivity analysis.

3. **Bowden et al. (2016)**. *Weighted median estimator*.
   - Link: https://pmc.ncbi.nlm.nih.gov/articles/PMC4849733/
   - Contribution: Consistent estimate if <50% of instrument weight is invalid.
   - Relevance: Included in Medusa; good baseline robustness check.

4. **Hemani, Tilling, Davey Smith (2017)**. *Steiger directionality*.
   - Link: https://pmc.ncbi.nlm.nih.gov/articles/PMC5711033/
   - Contribution: Direction-of-causation checks using variance explained.
   - Relevance: Medusa includes Steiger, but internal binary-outcome branch is intentionally unavailable.

### Bias, assumptions, and reporting
5. **Hartwig et al. (2016)**. *Two-sample MR: avoiding downsides*.
   - Link: https://pubmed.ncbi.nlm.nih.gov/27616674/
   - Contribution: Practical pitfalls and assumption stress-testing in two-sample MR.
   - Relevance: Supports Medusa’s sensitivity-analysis-first design philosophy.

6. **Burgess, Davies, Thompson (2016)**. *Bias due to participant overlap*.
   - Link: https://pubmed.ncbi.nlm.nih.gov/27625185/
   - Contribution: Quantifies overlap-induced weak-instrument bias behavior.
   - Relevance: Medusa currently does not model this bias explicitly.

7. **Mounier, Kutalik (2023)**. *MRlap overlap/winner's curse correction*.
   - Link: https://pubmed.ncbi.nlm.nih.gov/37036286/
   - Contribution: Analytic correction for overlap + winner’s curse + weak instruments.
   - Relevance: Strong candidate for optional post-hoc correction in Medusa coordinator workflows.

8. **Li, Morrison (2026 AJHG; preprint 2025)**. *Population mismatch in two-sample MR*.
   - Links: https://www.sciencedirect.com/science/article/pii/S0002929726000662 and https://pmc.ncbi.nlm.nih.gov/articles/PMC12324648/
   - Contribution: Shows attenuation bias toward zero under ancestry/population mismatch, including intra-continental mismatch.
   - Relevance: Directly relevant to multi-site observational networks using heterogeneous ancestry mixes.

9. **Skrivankova et al. (2021)**. *STROBE-MR statement*.
   - Links: https://pubmed.ncbi.nlm.nih.gov/34698778/ and https://pubmed.ncbi.nlm.nih.gov/34702754/
   - Contribution: 20-item reporting framework for transparent MR studies.
   - Relevance: Medusa report generation can be upgraded to be STROBE-MR-compliant by default.

10. **Burgess et al. (2023 update)**. *Guidelines for performing MR investigations*.
    - Link: https://pubmed.ncbi.nlm.nih.gov/32760811/
    - Contribution: Practical guidance covering instruments, harmonization, primary/sensitivity analyses, interpretation.
    - Relevance: Useful blueprint for default Medusa analysis templates and QC gates.

### Pleiotropy-robust and modern methods not yet in Medusa
11. **Verbanck et al. (2018)**. *MR-PRESSO*.
    - Link: https://pubmed.ncbi.nlm.nih.gov/29686387/
    - Contribution: Outlier detection/correction under horizontal pleiotropy.
    - Relevance: Not currently implemented; would materially strengthen robustness.

12. **Zhao et al. (2020)**. *MR-RAPS*.
    - Link: https://doi.org/10.1214/19-AOS1866
    - Contribution: Robust adjusted profile score, designed for many weak instruments.
    - Relevance: Important for weak-instrument settings common in observational cohorts.

13. **Morrison et al. (2020)**. *CAUSE*.
    - Link: https://pubmed.ncbi.nlm.nih.gov/32451458/
    - Contribution: Models correlated and uncorrelated pleiotropy.
    - Relevance: Not in Medusa; high-value optional method for difficult causal claims.

### Drug-target and cis-MR specific evidence
14. **Schmidt et al. (2020)**. *Genetic drug target validation using MR*.
    - Link: https://pubmed.ncbi.nlm.nih.gov/32591531/
    - Contribution: Drug-target MR framework; emphasizes cis-region design and robustness strategy.
    - Relevance: Conceptual template for Medusa’s drug-development positioning.

15. **Gordillo-Marañón et al. (2021)**. *Validation of lipid-related therapeutic targets*.
    - Link: https://pubmed.ncbi.nlm.nih.gov/34675202/
    - Contribution: End-to-end target triage: MR + replication + colocalization + external evidence.
    - Relevance: Demonstrates the workflow Medusa should support natively for drug target work.

16. **Zheng et al. (2020)**. *Phenome-wide MR of plasma proteome*.
    - Link: https://www.nature.com/articles/s41588-020-0682-6
    - Contribution: Showed many naïve MR hits fail colocalization; MR+coloc strongly improves target plausibility.
    - Relevance: The single biggest methodological gap in current Medusa is no native colocalization module.

17. **Lin, Pan (2024)**. *cisMR-cML robust cis-MR*.
    - Link: https://pubmed.ncbi.nlm.nih.gov/39025905/
    - Contribution: Robust cis-MR with correlated SNPs, conditional effects, invalid-IV tolerance.
    - Relevance: Highly relevant for protein/gene drug targets where cis instruments are correlated.

18. **van der Graaf et al. (2025)**. *MR-link-2*.
    - Link: https://www.nature.com/articles/s41467-025-60868-1
    - Contribution: Pleiotropy-robust cis-MR with improved type I error behavior in limited-region settings.
    - Relevance: Supports adding a modern cis-MR backend beyond clumped-IVW.

19. **Gkatzionis, Burgess, Newcombe (2023)**. *Statistical methods for cis-MR*.
    - Link: https://pubmed.ncbi.nlm.nih.gov/36273411/
    - Contribution: Technical review of numerical instability, LD conditioning, correlated-variant handling.
    - Relevance: Directly informs redesign of Medusa’s drug-target mode.

20. **Giambartolomei et al. (2014)** and **Wallace (2020)**. *Coloc framework and extensions*.
    - Links: https://pubmed.ncbi.nlm.nih.gov/24830394/ and https://pubmed.ncbi.nlm.nih.gov/32310995/
    - Contribution: Bayesian colocalization and handling multiple causal variants/priors.
    - Relevance: Essential companions to MR for drug-target claims.

21. **Minikel et al. (2024)**. *Genetic evidence and clinical success*.
    - Link: https://pubmed.ncbi.nlm.nih.gov/38632401/
    - Contribution: Drug mechanisms with genetic support had ~2.6x higher probability of success.
    - Relevance: Strong strategic justification for investment in high-quality MR target validation.

22. **Burgess et al. (2023 AJHG review)**. *Using genetic association data for drug discovery*.
    - Link: https://pubmed.ncbi.nlm.nih.gov/36736292/
    - Contribution: Integrative review of MR and related genetics tools in drug development.
    - Relevance: Good high-level map for Medusa roadmap prioritization.

### Federated/distributed methodology context (for Medusa architecture)
23. **Luo et al. (2022)**. *dPQL lossless distributed GLMM*.
    - Link: https://pubmed.ncbi.nlm.nih.gov/35579348/
    - Contribution: Demonstrates lossless distributed inference using institution-level summaries.
    - Relevance: Supports Medusa’s federated likelihood-pooling design rationale.

24. **Zhang et al. (2025)**. *COLA-GLM-H one-shot lossless distributed GLM*.
    - Link: https://pubmed.ncbi.nlm.nih.gov/41259033/
    - Contribution: One-round lossless heterogeneous distributed GLM estimation.
    - Relevance: Reinforces feasibility of Medusa’s one-shot federated model philosophy.

## 5. Comparison: Medusa methods vs current research consensus

| Method area | Medusa status | Literature expectation (2023-2026) | Assessment |
|---|---|---|---|
| Core two-sample causal estimator | Wald-ratio from profile MLE + delta SE | Valid under two-sample independence assumptions | **Correct** |
| Multi-site federation | Sum of site log-likelihood profiles | Lossless summary-level federation is acceptable if model assumptions hold | **Correct and innovative** |
| Standard sensitivity methods | IVW, MR-Egger, weighted median, Steiger, LOO | Minimum expected baseline | **Correct but minimal** |
| Drug-target cis-MR | Clumping + standard MR methods | Correlated cis-variant modeling, conditional effects, robust cis estimators | **Behind current best practice** |
| Colocalization | Not built in | MR + colocalization is now expected for target validation | **Major gap** |
| Overlap/winner's curse | No explicit correction | Growing consensus to evaluate/correct (MRlap and related) | **Gap** |
| Population mismatch robustness | No explicit modeling | Explicit checks/stratified analyses now strongly advised | **Gap** |
| Negative controls | Stub implementation | Expected in observational translational pipelines | **Gap** |
| Reporting standardization | HTML report exists | STROBE-MR-level transparency expected | **Partial** |

## 6. Are Medusa’s current methods correct?
### What is correct
1. **Primary estimator math is coherent**: profile-based beta_ZY, score-consistent beta_ZX, Wald ratio, and delta-method propagation are implemented correctly for two-sample assumptions.
2. **Federated pooling is statistically principled** under site independence and comparable model definitions.
3. **Sensitivity estimators are implemented in standard form** and optionally cross-checked via TwoSampleMR.
4. **Scientific validation tests exist and pass locally** when package is loaded (`NOT_CRAN=true`).

Executed checks in this workspace:
- `tests/testthat/test-runSensitivityAnalyses.R`: passed (76 assertions).
- `tests/testthat/test-scientificValidation.R`: passed (4 assertions).
- `tests/testthat/test-validateAgainstTwoSampleMR.R`: passed (7 assertions with full validation enabled).

### What remains methodologically incomplete for drug development
- No native colocalization/fine-mapping workflow.
- No modern cis-MR robust estimators for correlated IVs (e.g., cisMR-cML/MR-link-2 class).
- No explicit overlap/winner's-curse correction path.
- No explicit handling/reporting of population mismatch bias.
- Negative-control testing unfinished.

## 7. Prioritized improvement roadmap for Medusa

### Priority 1 (high impact, drug-target readiness)
1. **Add coloc integration** (`coloc`/PWCoCo-style support, posterior reporting, QC thresholds).
2. **Add robust cis-MR backend** (at least one of: cisMR-cML class, MR-link-2 style wrapper, or equivalent).
3. **Add overlap/winner's-curse diagnostics + optional correction** (MRlap-compatible pipeline).

### Priority 2 (robustness and interpretability)
4. Implement optional **MR-PRESSO / MR-RAPS / CAUSE** wrappers for pleiotropy and weak instruments.
5. Add **ancestry/population mismatch diagnostics** for exposure vs outcome summary sources and stratified analyses.
6. Complete **negative-control outcome testing** in diagnostics.

### Priority 3 (reporting quality and reproducibility)
7. Generate a **STROBE-MR appendix** in `generateMRReport()` automatically.
8. Add a required analysis metadata block: sample overlap declaration, population matching statement, colocalization evidence, and per-method decision log.

## 8. Practical interpretation for this package today
If used now for research prototyping, Medusa is defensible for:
- Federated baseline two-sample MR with standard sensitivity checks.

If used for drug-target go/no-go decisions, current evidence standards suggest adding (outside Medusa or as new modules):
- Colocalization,
- Robust cis-MR with correlated variants,
- Overlap/population-bias diagnostics,
- Extended pleiotropy-robust estimators.

Without these, target-prioritization claims risk being overstated even when the core MR estimate is technically correct.

## 9. Source links (primary)
- Burgess 2013: https://pmc.ncbi.nlm.nih.gov/articles/PMC4377079/
- Bowden 2015: https://pubmed.ncbi.nlm.nih.gov/26050253/
- Bowden 2016: https://pmc.ncbi.nlm.nih.gov/articles/PMC4849733/
- Hemani 2017: https://pmc.ncbi.nlm.nih.gov/articles/PMC5711033/
- Hartwig 2016: https://pubmed.ncbi.nlm.nih.gov/27616674/
- Burgess overlap 2016: https://pubmed.ncbi.nlm.nih.gov/27625185/
- Verbanck 2018: https://pubmed.ncbi.nlm.nih.gov/29686387/
- Zhao (MR-RAPS) 2020: https://doi.org/10.1214/19-AOS1866
- Morrison (CAUSE) 2020: https://pubmed.ncbi.nlm.nih.gov/32451458/
- Schmidt 2020: https://pubmed.ncbi.nlm.nih.gov/32591531/
- Gordillo-Marañón 2021: https://pubmed.ncbi.nlm.nih.gov/34675202/
- Zheng 2020: https://www.nature.com/articles/s41588-020-0682-6
- Gkatzionis 2023: https://pubmed.ncbi.nlm.nih.gov/36273411/
- Lin & Pan 2024: https://pubmed.ncbi.nlm.nih.gov/39025905/
- van der Graaf 2025: https://www.nature.com/articles/s41467-025-60868-1
- Li & Morrison 2025/2026: https://pmc.ncbi.nlm.nih.gov/articles/PMC12324648/ and https://www.sciencedirect.com/science/article/pii/S0002929726000662
- Giambartolomei 2014: https://pubmed.ncbi.nlm.nih.gov/24830394/
- Wallace 2020: https://pubmed.ncbi.nlm.nih.gov/32310995/
- STROBE-MR 2021: https://pubmed.ncbi.nlm.nih.gov/34698778/ and https://pubmed.ncbi.nlm.nih.gov/34702754/
- MR guidelines update 2023: https://pubmed.ncbi.nlm.nih.gov/32760811/
- Minikel 2024: https://pubmed.ncbi.nlm.nih.gov/38632401/
- AJHG review 2023: https://pubmed.ncbi.nlm.nih.gov/36736292/
- dPQL 2022: https://pubmed.ncbi.nlm.nih.gov/35579348/
- COLA-GLM-H 2025: https://pubmed.ncbi.nlm.nih.gov/41259033/
