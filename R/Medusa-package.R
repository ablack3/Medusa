# Copyright 2026 Observational Health Data Sciences and Informatics
#
# This file is part of Medusa
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' @title Medusa: Mendelian Estimation in Distributed Standardized Analytics
#'
#' @description
#' Medusa implements two-sample Mendelian Randomization (MR) within the OHDSI
#' ecosystem using the OMOP Common Data Model as the data substrate. The package
#' enables federated causal inference across distributed health networks without
#' requiring individual-level data to leave any site.
#'
#' The core methodological innovation is one-shot federated pooling via profile
#' likelihood aggregation: each site computes a log-likelihood profile across a
#' grid of parameter values and shares only that vector of numbers. The coordinator
#' sums profiles across sites to obtain a pooled estimate without any iterative
#' communication protocol.
#'
#' @details
#' The analysis pipeline consists of nine modules:
#'
#' \enumerate{
#'   \item \strong{Instrument Assembly} (\code{\link{getMRInstruments}}):
#'     Query OpenGWAS for GWAS summary statistics and apply LD clumping.
#'   \item \strong{Cohort Extraction} (\code{\link{buildMRCohort}}):
#'     Extract outcome cohorts and genotype data from OMOP CDM sites.
#'   \item \strong{Covariate Assembly} (\code{\link{buildMRCovariates}}):
#'     Assemble covariates via FeatureExtraction for adjustment and diagnostics.
#'   \item \strong{Instrument Diagnostics} (\code{\link{runInstrumentDiagnostics}}):
#'     Validate instruments via F-statistics, PheWAS, allele-frequency checks,
#'     missingness summaries, and a placeholder negative-control interface.
#'   \item \strong{Outcome Model} (\code{\link{fitOutcomeModel}}):
#'     Fit the binary outcome model and evaluate an exact or penalized profile
#'     log-likelihood on a grid.
#'   \item \strong{Likelihood Pooling} (\code{\link{poolLikelihoodProfiles}}):
#'     Aggregate site-level log-likelihood profiles via pointwise summation.
#'   \item \strong{MR Estimation} (\code{\link{computeMREstimate}}):
#'     Compute Wald ratio estimate with delta method standard errors.
#'   \item \strong{Sensitivity Analyses} (\code{\link{runSensitivityAnalyses}}):
#'     IVW, MR-Egger, weighted median, Steiger filtering, leave-one-out.
#'   \item \strong{Reporting} (\code{\link{generateMRReport}}):
#'     Generate self-contained HTML report with all results and diagnostics.
#' }
#'
#' @references
#' Davey Smith, G., & Hemani, G. (2014). Mendelian randomization: genetic anchors
#' for causal inference in epidemiological studies. \emph{Human Molecular Genetics},
#' 23(R1), R89-R98.
#'
#' Bowden, J., Davey Smith, G., & Burgess, S. (2015). Mendelian randomization with
#' invalid instruments: effect estimation and bias detection through Egger regression.
#' \emph{International Journal of Epidemiology}, 44(2), 512-525.
#'
#' Bowden, J., Davey Smith, G., Haycock, P. C., & Burgess, S. (2016). Consistent
#' estimation in Mendelian randomization with some invalid instruments using a
#' weighted median estimator. \emph{Genetic Epidemiology}, 40(4), 304-314.
#'
#' @importFrom stats coef confint lm glm binomial gaussian pchisq pnorm qchisq
#'   qnorm quantile var sd weighted.mean setNames spline approx median
#'   as.formula complete.cases rnorm rbinom runif
#' @importFrom rlang .data %||%
#' @importFrom checkmate assertCharacter assertNumeric assertLogical assertDataFrame
#'   assertSubset assertCount assertNumber assertString assertList
"_PACKAGE"
