# Medusa Oncology Upgrade Gaps

This document records the known gaps remaining after the first implementation
pass of the oncology target-validation roadmap.

## Methodology gaps

- `runCisMr()` currently supports GLS IVW and GLS Egger only.
  It does not yet implement stronger cis-MR methods for invalid or correlated
  instruments such as constrained ML, MR-RAPS-class approaches, or fine-mapped
  multi-signal estimators.
- `runColocalization()` is a native single-signal approximate Bayes factor
  implementation.
  It is not yet a full multi-signal coloc/fine-mapping workflow and does not
  expose posterior decomposition across multiple causal variants.
- `runOverlapDiagnostics()` provides heuristics and an approximate correction
  factor.
  It is not an MRlap-grade overlap and winner's-curse correction engine.
- The current decision engine is deterministic threshold logic.
  It does not yet support richer evidence synthesis across multiple exposure
  sources, datasets, or external drug-development priors.

## Data and workflow gaps

- The external dosage path supports standardized R data frames only.
  It does not yet include direct readers for PLINK, PGEN, BGEN, Hail, or VCF.
- The All of Us support is a documented adapter workflow, not a turnkey
  Workbench extraction pipeline.
  There is no packaged code yet for pulling dosage assets directly from All of
  Us infrastructure.
- The oncology outcome library is still recipe-level metadata.
  It does not yet ship full concept sets, cohort SQL, or ATLAS-ready
  definitions for the cancer endpoints in the research program.
- `runStudyManifest()` expects a caller-supplied executor.
  There is not yet a built-in end-to-end portfolio runner that fetches summary
  statistics, executes cohorts, runs all diagnostics, and assembles final
  reports automatically.

## Reporting and package gaps

- `generateTargetValidationReport()` is in place, but the package does not yet
  generate a full STROBE-MR appendix with all checklist items populated from
  structured provenance.
- Roxygen/man/pkgdown outputs were not regenerated in this increment.
  The new exported functions are available in code and `NAMESPACE`, but the
  manual pages and site artifacts still need a documentation rebuild.
- The current report templates focus on HTML output only.
  There is no compact machine-readable export for downstream portfolio tracking
  or triage dashboards.

## Validation gaps

- The new APIs are covered by unit and integration-style tests, and the full
  local test suite passes.
  However, there is not yet a gold-standard oncology validation package that
  reproduces known positive and negative target-control pairs from external
  benchmarks.
- The implementation has not yet been validated against real All of Us,
  UK Biobank, or FinnGen extraction pipelines in this repository.

## Scope intentionally deferred

- Time-to-event and survival modeling.
- Somatic tumor genomics and treatment-response modeling.
- Full clinico-genomic platform behavior for resources such as AACR GENIE BPC
  or Flatiron CGDB.
