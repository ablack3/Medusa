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

#' Generate Mendelian Randomization Analysis Report
#'
#' @title Self-contained HTML report for MR analysis
#'
#' @description
#' Produces a self-contained HTML report (single file, no external dependencies)
#' summarizing the entire MR analysis including instrument summary, likelihood
#' profile plots, main MR result, sensitivity analyses, PheWAS diagnostics,
#' and site contributions. Suitable for sharing with non-technical stakeholders.
#'
#' @param mrEstimate Output of \code{\link{computeMREstimate}}.
#' @param sensitivityResults Output of \code{\link{runSensitivityAnalyses}}.
#'   Can be NULL if no sensitivity analyses were run.
#' @param diagnosticResults Output of \code{\link{runInstrumentDiagnostics}}.
#'   Can be NULL if no diagnostics were run.
#' @param combinedProfile Output of \code{\link{poolLikelihoodProfiles}}.
#' @param siteProfileList Named list of site profile objects from
#'   \code{\link{fitOutcomeModel}}.
#' @param instrumentTable Output of \code{\link{getMRInstruments}}.
#' @param exposureLabel Character. Human-readable name for the exposure.
#'   Default is "Exposure".
#' @param outcomeLabel Character. Human-readable name for the outcome.
#'   Default is "Outcome".
#' @param outputPath Character. Path for the output HTML file. Default is
#'   "./Medusa_report.html".
#'
#' @return Character string with the path to the generated report (invisibly).
#'
#' @examples
#' \dontrun{
#' generateMRReport(
#'   mrEstimate = estimate,
#'   sensitivityResults = sensitivity,
#'   diagnosticResults = diagnostics,
#'   combinedProfile = combined,
#'   siteProfileList = siteProfiles,
#'   instrumentTable = instruments,
#'   exposureLabel = "IL-6 receptor levels",
#'   outcomeLabel = "Colorectal cancer"
#' )
#' }
#'
#' @seealso \code{\link{computeMREstimate}}, \code{\link{runSensitivityAnalyses}},
#'   \code{\link{runInstrumentDiagnostics}}
#'
#' @export
generateMRReport <- function(mrEstimate,
                             sensitivityResults = NULL,
                             diagnosticResults = NULL,
                             combinedProfile,
                             siteProfileList = NULL,
                             instrumentTable = NULL,
                             exposureLabel = "Exposure",
                             outcomeLabel = "Outcome",
                             outputPath = "./Medusa_report.html") {
  # Input validation
  checkmate::assertList(mrEstimate)
  checkmate::assertList(combinedProfile)
  checkmate::assertString(exposureLabel)
  checkmate::assertString(outcomeLabel)
  checkmate::assertString(outputPath)

  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("Package 'rmarkdown' is required for generateMRReport(). ",
         "Install it with: install.packages('rmarkdown')",
         call. = FALSE)
  }
  if (!requireNamespace("knitr", quietly = TRUE)) {
    stop("Package 'knitr' is required for generateMRReport(). ",
         "Install it with: install.packages('knitr')",
         call. = FALSE)
  }

  message("Generating Medusa MR analysis report...")

  # Find RMarkdown template
  rmdTemplate <- system.file("rmd", "MRReport.Rmd", package = "Medusa")
  if (rmdTemplate == "") {
    # Fallback for development
    rmdTemplate <- file.path(system.file(package = "Medusa"), "rmd", "MRReport.Rmd")
    if (!file.exists(rmdTemplate)) {
      stop("Cannot find report template MRReport.Rmd. Ensure the package is properly installed.")
    }
  }

  # Create output directory if needed
  outputDir <- dirname(outputPath)
  if (!dir.exists(outputDir)) {
    dir.create(outputDir, recursive = TRUE)
  }

  # Copy template to temp location for rendering
  tempRmd <- tempfile(fileext = ".Rmd")
  file.copy(rmdTemplate, tempRmd, overwrite = TRUE)

  # Render with parameters
  tryCatch(
    {
      rmarkdown::render(
        input = tempRmd,
        output_file = basename(outputPath),
        output_dir = outputDir,
        params = list(
          mrEstimate = mrEstimate,
          sensitivityResults = sensitivityResults,
          diagnosticResults = diagnosticResults,
          combinedProfile = combinedProfile,
          siteProfileList = siteProfileList,
          instrumentTable = instrumentTable,
          exposureLabel = exposureLabel,
          outcomeLabel = outcomeLabel
        ),
        quiet = TRUE
      )
    },
    error = function(e) {
      stop(sprintf("Report generation failed: %s", conditionMessage(e)))
    }
  )

  # Clean up temp file
  unlink(tempRmd)

  message(sprintf("Report saved to: %s", outputPath))
  invisible(outputPath)
}


#' Plot Profile Likelihood
#'
#' @title Likelihood profile visualization
#'
#' @description
#' Creates a ggplot2 visualization of the profile log-likelihood curve(s),
#' showing individual site profiles (if provided) and the combined profile
#' with MLE and confidence interval.
#'
#' @param combinedProfile Output of \code{\link{poolLikelihoodProfiles}}.
#' @param siteProfileList Optional named list of site profile objects.
#' @param mrEstimate Optional output of \code{\link{computeMREstimate}} for
#'   annotating the MLE and CI.
#' @param title Plot title. Default is "Profile Log-Likelihood".
#'
#' @return A ggplot2 object.
#'
#' @examples
#' profiles <- simulateSiteProfiles(nSites = 3, trueBeta = 0.5)
#' combined <- poolLikelihoodProfiles(profiles)
#' plotLikelihoodProfile(combined, profiles)
#'
#' @export
plotLikelihoodProfile <- function(combinedProfile,
                                  siteProfileList = NULL,
                                  mrEstimate = NULL,
                                  title = "Profile Log-Likelihood") {
  plotData <- data.frame(
    beta = combinedProfile$betaGrid,
    logLik = combinedProfile$logLikProfile,
    source = "Combined"
  )

  # Add site profiles
  if (!is.null(siteProfileList)) {
    for (siteName in names(siteProfileList)) {
      profile <- siteProfileList[[siteName]]
      siteData <- data.frame(
        beta = profile$betaGrid,
        logLik = profile$logLikProfile - max(profile$logLikProfile),
        source = profile$siteId
      )
      plotData <- rbind(plotData, siteData)
    }
  }

  p <- ggplot2::ggplot(plotData, ggplot2::aes(x = .data$beta, y = .data$logLik)) +
    ggplot2::geom_line(
      data = plotData[plotData$source != "Combined", ],
      ggplot2::aes(group = .data$source),
      color = MR_COLORS$site, alpha = 0.6, linewidth = 0.5
    ) +
    ggplot2::geom_line(
      data = plotData[plotData$source == "Combined", ],
      color = MR_COLORS$combined, linewidth = 1.2
    ) +
    ggplot2::labs(
      title = title,
      x = expression(beta[ZY]),
      y = "Profile log-likelihood"
    ) +
    mrTheme()

  # Add MLE and CI annotation
  if (!is.null(mrEstimate)) {
    p <- p +
      ggplot2::geom_vline(
        xintercept = mrEstimate$betaZY,
        linetype = "dashed", color = MR_COLORS$primary
      ) +
      ggplot2::geom_hline(
        yintercept = -qchisq(mrEstimate$ciLevel, df = 1) / 2,
        linetype = "dotted", color = MR_COLORS$warning, alpha = 0.5
      )
  }

  p
}


#' Plot Sensitivity Analysis Forest Plot
#'
#' @title Forest plot comparing MR methods
#'
#' @description
#' Creates a horizontal forest plot showing estimates from all sensitivity
#' analysis methods side by side with confidence intervals.
#'
#' @param sensitivityResults Output of \code{\link{runSensitivityAnalyses}}.
#'
#' @return A ggplot2 object.
#'
#' @examples
#' set.seed(42)
#' nSnps <- 10
#' perSnp <- data.frame(
#'   snp_id = paste0("rs", 1:nSnps),
#'   effect_allele = rep(c("A", "C", "G", "T", "A"), length.out = nSnps),
#'   other_allele = rep(c("C", "G", "T", "A", "C"), length.out = nSnps),
#'   eaf = seq(0.1, 0.7, length.out = nSnps),
#'   beta_ZY = rnorm(nSnps, 0.15, 0.02),
#'   se_ZY = rep(0.02, nSnps),
#'   beta_ZX = rnorm(nSnps, 0.3, 0.05),
#'   se_ZX = rep(0.05, nSnps)
#' )
#' results <- runSensitivityAnalyses(perSnp, engine = "internal")
#' plotSensitivityForest(results)
#'
#' @export
plotSensitivityForest <- function(sensitivityResults) {
  if (is.null(sensitivityResults$summary) || nrow(sensitivityResults$summary) == 0) {
    message("No sensitivity analysis results to plot.")
    return(ggplot2::ggplot() + ggplot2::theme_void())
  }

  plotData <- sensitivityResults$summary
  plotData$method <- factor(plotData$method, levels = rev(plotData$method))

  ggplot2::ggplot(plotData, ggplot2::aes(
    x = .data$beta_MR, y = .data$method,
    xmin = .data$ci_lower, xmax = .data$ci_upper
  )) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_errorbarh(height = 0.2, linewidth = 0.8,
                             color = MR_COLORS$primary) +
    ggplot2::geom_point(size = 3, color = MR_COLORS$primary) +
    ggplot2::labs(
      title = "Sensitivity Analysis Comparison",
      x = expression(beta[MR]),
      y = NULL
    ) +
    mrTheme()
}
