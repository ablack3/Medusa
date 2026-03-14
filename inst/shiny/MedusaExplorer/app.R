# Medusa Results Explorer
# Interactive bslib dashboard for MR diagnostics and results
#
# Launched via Medusa::launchResultsExplorer()

library(shiny)
library(bslib)
library(ggplot2)

# ---------------------------------------------------------------------------
# Retrieve results from shinyOptions (set by launchResultsExplorer)
# ---------------------------------------------------------------------------
getResults <- function() {
  res <- getShinyOption("medusaResults")
  if (is.null(res)) {
    stop("No results found. Launch this app via Medusa::launchResultsExplorer().")
  }
  res
}

# ---------------------------------------------------------------------------
# Theme & palette (matches pkgdown website)
# ---------------------------------------------------------------------------
MEDUSA_COLORS <- list(
  navyBlue   = "#244268",
  ohdsiBlue  = "#264553",
  blueLight  = "#3d6573",
  orange     = "#E69528",
  text       = "#333333",
  background = "#ffffff",
  codeBg     = "#f6f8fa",
  border     = "#e8ecef",
  success    = "#27AE60",
  warning    = "#E74C3C",
  neutral    = "#95A5A6",
  primary    = "#1B4F72",
  secondary  = "#85C1E9",
  methods    = c(
    IVW             = "#1B4F72",
    `MR-Egger`      = "#E74C3C",
    `Weighted Median` = "#27AE60",
    `Steiger-filtered IVW` = "#F39C12",
    `Leave-One-Out` = "#8E44AD"
  )
)

medusaTheme <- bs_theme(
  version  = 5,
  bootswatch = "default",
  primary  = MEDUSA_COLORS$navyBlue,
  success  = MEDUSA_COLORS$success,
  danger   = MEDUSA_COLORS$warning,
  info     = MEDUSA_COLORS$ohdsiBlue,
  warning  = MEDUSA_COLORS$orange,
  fg       = MEDUSA_COLORS$text,
  bg       = MEDUSA_COLORS$background,
  base_font = font_google("Source Sans 3"),
  heading_font = font_google("Source Sans 3"),
  "navbar-bg" = MEDUSA_COLORS$navyBlue,
  "navbar-dark-color" = "#ffffff",
  "navbar-dark-hover-color" = MEDUSA_COLORS$orange,
  "navbar-dark-active-color" = MEDUSA_COLORS$orange
)

ggTheme <- function(baseSize = 13) {
  theme_minimal(base_size = baseSize, base_family = "Source Sans 3") +
    theme(
      plot.title    = element_text(face = "bold", size = baseSize + 2,
                                    color = MEDUSA_COLORS$ohdsiBlue),
      plot.subtitle = element_text(color = "grey40"),
      axis.title    = element_text(face = "bold"),
      axis.text     = element_text(color = "grey30"),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
      strip.text = element_text(face = "bold")
    )
}

alleleLabelLayer <- function() {
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    return(ggrepel::geom_text_repel(
      aes(label = .data$snp_id),
      size = 3,
      max.overlaps = 15
    ))
  }

  geom_text(
    aes(label = .data$snp_id),
    size = 3,
    check_overlap = TRUE,
    nudge_y = 0.02
  )
}

# ---------------------------------------------------------------------------
# Helper: value box card
# ---------------------------------------------------------------------------
statusCard <- function(title, value, theme = "primary", icon = NULL) {
  value_box(
    title = title,
    value = value,
    showcase = if (!is.null(icon)) bsicons::bs_icon(icon) else NULL,
    theme = theme
  )
}

# ---------------------------------------------------------------------------
# Helper: format p-value
# ---------------------------------------------------------------------------
fmtP <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return(sprintf("%.2e", p))
  sprintf("%.4f", p)
}

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------
ui <- page_navbar(
  title = tags$span(
    tags$strong("Medusa"),
    tags$span("Results Explorer", style = "font-weight: 300; margin-left: 4px;")
  ),
  id = "nav",
  theme = medusaTheme,
  bg = MEDUSA_COLORS$navyBlue,

  # ---- Overview ----
  nav_panel(
    title = "Overview",
    icon = bsicons::bs_icon("speedometer2"),
    layout_columns(
      fill = FALSE,
      col_widths = c(3, 3, 3, 3),
      uiOutput("overviewCards")
    ),
    layout_columns(
      col_widths = c(8, 4),
      card(
        card_header("Sensitivity Analysis Comparison"),
        plotOutput("overviewForest", height = "350px")
      ),
      card(
        card_header("Diagnostic Flags"),
        uiOutput("diagnosticFlagCards")
      )
    )
  ),

  # ---- Instrument Diagnostics ----
  nav_panel(
    title = "Instrument Diagnostics",
    icon = bsicons::bs_icon("clipboard2-check"),
    navset_card_tab(
      title = "MR Assumption Checks",
      nav_panel(
        "F-Statistics (Relevance)",
        p(class = "text-muted mt-2",
          tags$strong("IV Assumption 1 (Relevance):"),
          " Instruments must be strongly associated with the exposure. F > 10 indicates adequate strength."
        ),
        layout_columns(
          col_widths = c(7, 5),
          card(
            card_header("F-Statistic by SNP"),
            plotOutput("fStatPlot", height = "350px")
          ),
          card(
            card_header("F-Statistics Table"),
            tableOutput("fStatTable")
          )
        )
      ),
      nav_panel(
        "PheWAS (Independence)",
        p(class = "text-muted mt-2",
          tags$strong("IV Assumption 2 (Independence):"),
          " Instruments should not be associated with confounders.",
          " Significant SNP-covariate associations (after Bonferroni correction) may indicate violation."
        ),
        layout_columns(
          col_widths = c(7, 5),
          card(
            card_header("PheWAS Results"),
            plotOutput("phewasPlot", height = "400px")
          ),
          card(
            card_header("Significant Associations"),
            tableOutput("phewasTable")
          )
        )
      ),
      nav_panel(
        "Allele Frequencies",
        p(class = "text-muted mt-2",
          "Compares GWAS effect allele frequency to cohort-observed frequency.",
          " Discrepancy > 0.10 may indicate strand issues or population mismatch."
        ),
        layout_columns(
          col_widths = c(7, 5),
          card(
            card_header("GWAS vs Cohort EAF"),
            plotOutput("afPlot", height = "350px")
          ),
          card(
            card_header("Allele Frequency Table"),
            tableOutput("afTable")
          )
        )
      ),
      nav_panel(
        "Genotype Missingness",
        p(class = "text-muted mt-2",
          "SNPs with > 10% missingness are flagged. High missingness reduces power",
          " and may introduce bias if missing-not-at-random."
        ),
        card(
          card_header("Missingness by SNP"),
          layout_columns(
            col_widths = c(7, 5),
            plotOutput("missingnessPlot", height = "300px"),
            tableOutput("missingnessTable")
          )
        )
      ),
      nav_panel(
        "Heterogeneity",
        p(class = "text-muted mt-2",
          tags$strong("IV Assumption 3 (Exclusion Restriction):"),
          " If all instruments are valid, per-SNP estimates should be homogeneous.",
          " Significant Cochran's Q or high I\u00B2 suggests potential pleiotropy."
        ),
        layout_columns(
          col_widths = c(6, 6),
          card(
            card_header("Heterogeneity Statistics"),
            tableOutput("heterogeneityTable")
          ),
          card(
            card_header("MR-Egger Intercept Test"),
            p(class = "text-muted",
              "A non-zero intercept indicates directional pleiotropy",
              " (violation of the exclusion restriction)."
            ),
            tableOutput("eggerInterceptTable")
          )
        )
      )
    )
  ),

  # ---- Results ----
  nav_panel(
    title = "Results",
    icon = bsicons::bs_icon("graph-up"),
    layout_columns(
      fill = FALSE,
      col_widths = c(4, 4, 4),
      uiOutput("resultCards")
    ),
    layout_columns(
      col_widths = c(7, 5),
      card(
        card_header("Profile Log-Likelihood"),
        plotOutput("profilePlot", height = "400px")
      ),
      card(
        card_header("Primary Estimate"),
        tableOutput("primaryEstimateTable")
      )
    )
  ),

  # ---- Sensitivity Analyses ----
  nav_panel(
    title = "Sensitivity Analyses",
    icon = bsicons::bs_icon("shield-check"),
    navset_card_tab(
      title = "Robustness Checks",
      nav_panel(
        "Method Comparison",
        p(class = "text-muted mt-2",
          "Concordance of estimates across methods robust to different",
          " patterns of instrument invalidity strengthens causal evidence."
        ),
        layout_columns(
          col_widths = c(7, 5),
          card(plotOutput("sensitivityForest", height = "350px")),
          card(
            card_header("Method Summary"),
            tableOutput("sensitivitySummaryTable")
          )
        )
      ),
      nav_panel(
        "Scatter Plot",
        p(class = "text-muted mt-2",
          "Per-SNP exposure effects (\u03B2_ZX) vs outcome effects (\u03B2_ZY).",
          " Regression lines show fitted slopes from each MR method."
        ),
        card(plotOutput("scatterPlot", height = "450px"))
      ),
      nav_panel(
        "Funnel Plot",
        p(class = "text-muted mt-2",
          "Per-SNP Wald ratio estimates vs precision (1/SE).",
          " Asymmetry suggests directional pleiotropy."
        ),
        card(plotOutput("funnelPlot", height = "400px"))
      ),
      nav_panel(
        "Leave-One-Out",
        p(class = "text-muted mt-2",
          "IVW estimate recomputed after dropping each SNP.",
          " Identifies influential outlier instruments."
        ),
        card(plotOutput("looPlot", height = "450px"))
      )
    )
  ),

  # ---- Negative Controls ----
  nav_panel(
    title = "Negative Controls",
    icon = bsicons::bs_icon("bullseye"),
    p(class = "text-muted mt-2",
      tags$strong("Empirical Calibration:"),
      " Negative control outcomes are conditions with no expected causal effect.",
      " The MR pipeline is run on each NC outcome to assess systematic bias.",
      " A well-calibrated analysis should produce null estimates for these outcomes."
    ),
    layout_columns(
      fill = FALSE,
      col_widths = c(4, 4, 4),
      uiOutput("ncSummaryCards")
    ),
    layout_columns(
      col_widths = c(7, 5),
      card(
        card_header("Negative Control Effect Estimates"),
        plotOutput("ncEffectPlot", height = "400px")
      ),
      card(
        card_header("Calibrated Primary Estimate"),
        uiOutput("ncCalibratedCard"),
        tags$hr(),
        card_header("NC Estimates Table"),
        tableOutput("ncEstimatesTable")
      )
    )
  ),

  # ---- About ----
  nav_panel(
    title = "About",
    icon = bsicons::bs_icon("info-circle"),
    card(
      card_header("About Medusa"),
      tags$div(
        class = "p-3",
        tags$h5("Mendelian Estimation in Distributed Standardized Analytics"),
        tags$p("Medusa implements two-sample Mendelian Randomization within the OHDSI",
               "ecosystem using the OMOP Common Data Model. It supports federated analysis",
               "across distributed health networks via one-shot profile likelihood pooling."),
        tags$h6("MR Assumptions Validated"),
        tags$ol(
          tags$li(tags$strong("Relevance:"), " F-statistics for instrument strength (F > 10)"),
          tags$li(tags$strong("Independence:"), " PheWAS tests for SNP-confounder associations"),
          tags$li(tags$strong("Exclusion Restriction:"), " Multiple sensitivity methods",
                  " (MR-Egger intercept, heterogeneity, Steiger directionality)")
        ),
        tags$h6("Diagnostic Checks"),
        tags$ul(
          tags$li("Per-SNP F-statistics with weak instrument threshold"),
          tags$li("Instrument PheWAS across OMOP covariate domains (Bonferroni-corrected)"),
          tags$li("Allele frequency GWAS-vs-cohort comparison"),
          tags$li("Genotype missingness summary"),
          tags$li("Cochran's Q and I\u00B2 heterogeneity statistics"),
          tags$li("MR-Egger intercept test for directional pleiotropy"),
          tags$li("Steiger directionality filtering"),
          tags$li("Leave-one-out influence analysis"),
          tags$li("Negative control outcome empirical calibration")
        ),
        tags$hr(),
        tags$p(class = "text-muted",
               "Medusa v0.1.0 | Apache License 2.0 | OHDSI")
      )
    )
  )
)


# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------
server <- function(input, output, session) {
  results <- getResults()

  mrEst    <- reactive(results$mrEstimate)
  sens     <- reactive(results$sensitivityResults)
  diag     <- reactive(results$diagnosticResults)
  combined <- reactive(results$combinedProfile)
  sites    <- reactive(results$siteProfileList)
  inst     <- reactive(results$instrumentTable)
  perSnp   <- reactive(results$perSnpEstimates)
  ncRes    <- reactive(results$negativeControlResults)

  # ========================================================================
  # OVERVIEW
  # ========================================================================
  output$overviewCards <- renderUI({
    est <- mrEst()
    if (is.null(est)) return(NULL)
    tagList(
      statusCard("Causal Estimate (\u03B2)",
                 sprintf("%.4f", est$betaMR), "primary", "graph-up-arrow"),
      statusCard("Odds Ratio",
                 sprintf("%.3f", est$oddsRatio), "primary", "calculator"),
      statusCard("P-value",
                 fmtP(est$pValue),
                 if (est$pValue < 0.05) "success" else "secondary",
                 "check-circle"),
      statusCard("N Instruments",
                 as.character(est$nInstruments), "info", "list-ol")
    )
  })

  output$overviewForest <- renderPlot({
    s <- sens()
    if (is.null(s) || is.null(s$summary) || nrow(s$summary) == 0) return(NULL)
    plotData <- s$summary
    plotData$method <- factor(plotData$method, levels = rev(plotData$method))
    methodColors <- MEDUSA_COLORS$methods[as.character(plotData$method)]
    methodColors[is.na(methodColors)] <- MEDUSA_COLORS$primary

    ggplot(plotData, aes(x = .data$beta_MR, y = .data$method,
                         xmin = .data$ci_lower, xmax = .data$ci_upper)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      geom_errorbarh(height = 0.2, linewidth = 0.8,
                     color = methodColors) +
      geom_point(size = 3.5, color = methodColors) +
      labs(x = expression(beta[MR]), y = NULL,
           title = "Sensitivity Analysis Comparison") +
      ggTheme()
  })

  output$diagnosticFlagCards <- renderUI({
    d <- diag()
    if (is.null(d)) {
      return(tags$p(class = "text-muted", "No diagnostics available."))
    }
    flags <- d$diagnosticFlags
    flagItems <- lapply(names(flags), function(nm) {
      label <- switch(nm,
        weakInstruments       = "Weak Instruments (F < 10)",
        phewasSignificant     = "PheWAS Significant Associations",
        negativeControlFailure = "Negative Control Failure",
        alleleFreqDiscrepancy = "Allele Frequency Discrepancy",
        highMissingness       = "High Genotype Missingness",
        nm
      )
      icon <- if (flags[nm]) bsicons::bs_icon("exclamation-triangle-fill") else bsicons::bs_icon("check-circle-fill")
      color <- if (flags[nm]) MEDUSA_COLORS$warning else MEDUSA_COLORS$success
      tags$div(
        class = "d-flex align-items-center mb-2 p-2 rounded",
        style = sprintf("border-left: 4px solid %s; background-color: %s10;",
                        color, color),
        tags$span(icon, style = sprintf("color: %s; margin-right: 8px;", color)),
        tags$span(label, style = "font-weight: 500;")
      )
    })
    tagList(flagItems)
  })

  # ========================================================================
  # INSTRUMENT DIAGNOSTICS
  # ========================================================================

  # -- F-Statistics --
  output$fStatPlot <- renderPlot({
    d <- diag()
    if (is.null(d) || is.null(d$fStatistics)) return(NULL)
    fs <- d$fStatistics
    fs$snp_id <- factor(fs$snp_id, levels = fs$snp_id[order(fs$fStatistic)])
    fs$color <- ifelse(fs$weakFlag, MEDUSA_COLORS$warning, MEDUSA_COLORS$ohdsiBlue)

    ggplot(fs, aes(x = .data$fStatistic, y = .data$snp_id)) +
      geom_vline(xintercept = 10, linetype = "dashed", color = MEDUSA_COLORS$orange,
                 linewidth = 0.8) +
      geom_segment(aes(x = 0, xend = .data$fStatistic,
                       y = .data$snp_id, yend = .data$snp_id),
                   color = "grey70", linewidth = 0.4) +
      geom_point(size = 3, color = fs$color) +
      annotate("text", x = 10, y = 0.5, label = "F = 10",
               hjust = -0.1, color = MEDUSA_COLORS$orange, fontface = "italic",
               size = 3.5) +
      labs(x = "F-Statistic", y = NULL, title = "Instrument Strength") +
      ggTheme()
  })

  output$fStatTable <- renderTable({
    d <- diag()
    if (is.null(d) || is.null(d$fStatistics)) return(NULL)
    fs <- d$fStatistics
    fs$Status <- ifelse(fs$weakFlag, "Weak", "OK")
    fs$fStatistic <- round(fs$fStatistic, 1)
    fs[, c("snp_id", "fStatistic", "source", "Status")]
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  # -- PheWAS --
  output$phewasPlot <- renderPlot({
    d <- diag()
    if (is.null(d) || is.null(d$phewasResults) || nrow(d$phewasResults) == 0) {
      return(NULL)
    }
    pw <- d$phewasResults
    pw$negLogP <- -log10(pw$pval)
    threshold <- if (any(pw$significant)) {
      min(pw$negLogP[pw$significant])
    } else {
      -log10(0.05 / nrow(pw))
    }

    ggplot(pw, aes(x = .data$covariate_name, y = .data$negLogP,
                   color = .data$significant)) +
      geom_hline(yintercept = threshold, linetype = "dashed",
                 color = MEDUSA_COLORS$orange) +
      geom_point(size = 2, alpha = 0.7) +
      scale_color_manual(values = c("FALSE" = MEDUSA_COLORS$neutral,
                                     "TRUE"  = MEDUSA_COLORS$warning),
                         labels = c("Not Significant", "Significant"),
                         name = NULL) +
      labs(x = NULL, y = expression(-log[10](p)),
           title = "Instrument PheWAS") +
      coord_flip() +
      ggTheme() +
      theme(legend.position = "bottom")
  })

  output$phewasTable <- renderTable({
    d <- diag()
    if (is.null(d) || is.null(d$phewasResults)) return(NULL)
    pw <- d$phewasResults
    if (nrow(pw) == 0) return(data.frame(Message = "No associations tested."))
    sig <- pw[pw$significant, , drop = FALSE]
    if (nrow(sig) == 0) return(data.frame(Message = "No significant associations found."))
    sig$pval <- sapply(sig$pval, fmtP)
    sig$beta <- round(sig$beta, 4)
    sig$se <- round(sig$se, 4)
    sig[, c("snp_id", "covariate_name", "beta", "se", "pval")]
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  # -- Allele Frequencies --
  output$afPlot <- renderPlot({
    d <- diag()
    if (is.null(d) || is.null(d$afComparison)) return(NULL)
    af <- d$afComparison
    af <- af[!is.na(af$eaf_cohort), , drop = FALSE]
    if (nrow(af) == 0) return(NULL)

    ggplot(af, aes(x = .data$eaf_gwas, y = .data$eaf_cohort)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      geom_abline(slope = 1, intercept = 0.1, linetype = "dotted",
                  color = MEDUSA_COLORS$orange, alpha = 0.5) +
      geom_abline(slope = 1, intercept = -0.1, linetype = "dotted",
                  color = MEDUSA_COLORS$orange, alpha = 0.5) +
      geom_point(aes(color = .data$discrepancyFlag), size = 3) +
      alleleLabelLayer() +
      scale_color_manual(values = c("FALSE" = MEDUSA_COLORS$ohdsiBlue,
                                     "TRUE"  = MEDUSA_COLORS$warning),
                         labels = c("OK", "Discrepancy > 0.1"),
                         name = NULL) +
      labs(x = "EAF (GWAS)", y = "EAF (Cohort)",
           title = "Allele Frequency Comparison") +
      coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
      ggTheme()
  })

  output$afTable <- renderTable({
    d <- diag()
    if (is.null(d) || is.null(d$afComparison)) return(NULL)
    af <- d$afComparison
    af$eaf_gwas <- round(af$eaf_gwas, 3)
    af$eaf_cohort <- round(af$eaf_cohort, 3)
    af$eaf_diff <- round(af$eaf_diff, 3)
    af$Status <- ifelse(af$discrepancyFlag, "Flag", "OK")
    af[, c("snp_id", "eaf_gwas", "eaf_cohort", "eaf_diff", "Status")]
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  # -- Missingness --
  output$missingnessPlot <- renderPlot({
    d <- diag()
    if (is.null(d) || is.null(d$missingnessReport)) return(NULL)
    mr <- d$missingnessReport
    mr$snp_id <- factor(mr$snp_id, levels = mr$snp_id[order(mr$pct_missing)])
    mr$color <- ifelse(mr$highMissingFlag, MEDUSA_COLORS$warning, MEDUSA_COLORS$ohdsiBlue)

    ggplot(mr, aes(x = .data$pct_missing, y = .data$snp_id)) +
      geom_vline(xintercept = 10, linetype = "dashed",
                 color = MEDUSA_COLORS$orange) +
      geom_segment(aes(x = 0, xend = .data$pct_missing,
                       y = .data$snp_id, yend = .data$snp_id),
                   color = "grey70", linewidth = 0.4) +
      geom_point(size = 3, color = mr$color) +
      labs(x = "% Missing", y = NULL,
           title = "Genotype Missingness") +
      ggTheme()
  })

  output$missingnessTable <- renderTable({
    d <- diag()
    if (is.null(d) || is.null(d$missingnessReport)) return(NULL)
    mr <- d$missingnessReport
    mr$pct_missing <- round(mr$pct_missing, 1)
    mr$Status <- ifelse(mr$highMissingFlag, "High", "OK")
    mr[, c("snp_id", "n_total", "n_missing", "pct_missing", "Status")]
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  # -- Heterogeneity --
  output$heterogeneityTable <- renderTable({
    s <- sens()
    if (is.null(s) || is.null(s$ivw)) return(NULL)
    ivw <- s$ivw
    df <- data.frame(
      Statistic = c("Cochran's Q", "Q p-value", "I\u00B2 (%)"),
      Value = c(
        if (!is.na(ivw$cochran_Q)) sprintf("%.2f", ivw$cochran_Q) else "NA",
        if (!is.na(ivw$cochran_Q_pval)) fmtP(ivw$cochran_Q_pval) else "NA",
        if ("i_squared" %in% names(ivw) && !is.na(ivw$i_squared)) {
          sprintf("%.1f%%", ivw$i_squared)
        } else {
          "NA"
        }
      ),
      Interpretation = c(
        if (!is.na(ivw$cochran_Q_pval) && ivw$cochran_Q_pval < 0.05) {
          "Significant heterogeneity"
        } else {
          "No significant heterogeneity"
        },
        "",
        if ("i_squared" %in% names(ivw) && !is.na(ivw$i_squared)) {
          if (ivw$i_squared > 75) "High heterogeneity"
          else if (ivw$i_squared > 50) "Moderate heterogeneity"
          else if (ivw$i_squared > 25) "Low heterogeneity"
          else "Very low heterogeneity"
        } else {
          ""
        }
      ),
      stringsAsFactors = FALSE
    )
    df
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  output$eggerInterceptTable <- renderTable({
    s <- sens()
    if (is.null(s) || is.null(s$mrEgger)) {
      return(data.frame(Message = "MR-Egger not available (requires >= 4 SNPs)."))
    }
    eg <- s$mrEgger
    df <- data.frame(
      Statistic = c("Intercept", "Intercept SE", "Intercept p-value"),
      Value = c(
        sprintf("%.4f", eg$intercept),
        sprintf("%.4f", eg$intercept_se),
        fmtP(eg$intercept_pval)
      ),
      Interpretation = c(
        "",
        "",
        if (!is.na(eg$intercept_pval) && eg$intercept_pval < 0.05) {
          "Evidence of directional pleiotropy"
        } else {
          "No evidence of directional pleiotropy"
        }
      ),
      stringsAsFactors = FALSE
    )
    df
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  # ========================================================================
  # RESULTS
  # ========================================================================
  output$resultCards <- renderUI({
    est <- mrEst()
    if (is.null(est)) return(NULL)
    tagList(
      statusCard(
        "Causal \u03B2 (95% CI)",
        sprintf("%.4f (%.4f, %.4f)", est$betaMR, est$ciLower, est$ciUpper),
        "primary", "graph-up-arrow"
      ),
      statusCard(
        "Odds Ratio (95% CI)",
        sprintf("%.3f (%.3f, %.3f)", est$oddsRatio, est$orCiLower, est$orCiUpper),
        "primary", "calculator"
      ),
      statusCard(
        "P-value",
        fmtP(est$pValue),
        if (est$pValue < 0.05) "success" else "secondary",
        "check-circle"
      )
    )
  })

  output$profilePlot <- renderPlot({
    comb <- combined()
    est <- mrEst()
    if (is.null(comb)) return(NULL)

    plotData <- data.frame(
      beta = comb$betaGrid,
      logLik = comb$logLikProfile
    )

    p <- ggplot(plotData, aes(x = .data$beta, y = .data$logLik)) +
      geom_line(color = MEDUSA_COLORS$primary, linewidth = 1.2) +
      labs(x = expression(beta[ZY]),
           y = "Profile log-likelihood",
           title = "Combined Profile Log-Likelihood") +
      ggTheme()

    # Add site profiles
    sl <- sites()
    if (!is.null(sl) && length(sl) > 0) {
      for (nm in names(sl)) {
        sp <- sl[[nm]]
        siteData <- data.frame(
          beta = sp$betaGrid,
          logLik = sp$logLikProfile - max(sp$logLikProfile)
        )
        p <- p + geom_line(data = siteData, color = MEDUSA_COLORS$secondary,
                           alpha = 0.5, linewidth = 0.5)
      }
    }

    if (!is.null(est)) {
      p <- p +
        geom_vline(xintercept = est$betaZY, linetype = "dashed",
                   color = MEDUSA_COLORS$primary) +
        geom_hline(yintercept = max(plotData$logLik) - qchisq(0.95, 1) / 2,
                   linetype = "dotted", color = MEDUSA_COLORS$warning, alpha = 0.5)
    }
    p
  })

  output$primaryEstimateTable <- renderTable({
    est <- mrEst()
    if (is.null(est)) return(NULL)
    data.frame(
      Parameter = c(
        "\u03B2_ZY (SNP-Outcome)", "SE(\u03B2_ZY)",
        "\u03B2_ZX (SNP-Exposure)", "SE(\u03B2_ZX)",
        "\u03B2_MR (Causal)", "SE(\u03B2_MR)",
        "95% CI", "P-value",
        "Odds Ratio", "OR 95% CI",
        "N Instruments"
      ),
      Value = c(
        sprintf("%.4f", est$betaZY), sprintf("%.4f", est$seZY),
        sprintf("%.4f", est$betaZX), sprintf("%.4f", est$seZX),
        sprintf("%.4f", est$betaMR), sprintf("%.4f", est$seMR),
        sprintf("(%.4f, %.4f)", est$ciLower, est$ciUpper),
        fmtP(est$pValue),
        sprintf("%.3f", est$oddsRatio),
        sprintf("(%.3f, %.3f)", est$orCiLower, est$orCiUpper),
        as.character(est$nInstruments)
      ),
      stringsAsFactors = FALSE
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  # ========================================================================
  # SENSITIVITY ANALYSES
  # ========================================================================
  output$sensitivityForest <- renderPlot({
    s <- sens()
    if (is.null(s) || is.null(s$summary) || nrow(s$summary) == 0) return(NULL)
    plotData <- s$summary
    plotData$method <- factor(plotData$method, levels = rev(plotData$method))
    methodColors <- MEDUSA_COLORS$methods[as.character(plotData$method)]
    methodColors[is.na(methodColors)] <- MEDUSA_COLORS$primary

    ggplot(plotData, aes(x = .data$beta_MR, y = .data$method,
                         xmin = .data$ci_lower, xmax = .data$ci_upper)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      geom_errorbarh(height = 0.2, linewidth = 0.8, color = methodColors) +
      geom_point(size = 3.5, color = methodColors) +
      labs(x = expression(beta[MR]), y = NULL,
           title = "MR Methods Comparison") +
      ggTheme()
  })

  output$sensitivitySummaryTable <- renderTable({
    s <- sens()
    if (is.null(s) || is.null(s$summary)) return(NULL)
    df <- s$summary
    df$beta_MR <- round(df$beta_MR, 4)
    df$se_MR <- round(df$se_MR, 4)
    df$ci_lower <- round(df$ci_lower, 4)
    df$ci_upper <- round(df$ci_upper, 4)
    df$pval <- sapply(df$pval, fmtP)
    df$OR <- round(exp(s$summary$beta_MR), 3)
    df
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  # -- Scatter Plot --
  output$scatterPlot <- renderPlot({
    ps <- perSnp()
    s <- sens()
    if (is.null(ps)) return(NULL)

    p <- ggplot(ps, aes(x = .data$beta_ZX, y = .data$beta_ZY)) +
      geom_hline(yintercept = 0, linetype = "dotted", color = "grey70") +
      geom_vline(xintercept = 0, linetype = "dotted", color = "grey70") +
      geom_errorbar(aes(ymin = .data$beta_ZY - 1.96 * .data$se_ZY,
                        ymax = .data$beta_ZY + 1.96 * .data$se_ZY),
                    width = 0, color = "grey60", alpha = 0.6) +
      geom_errorbarh(aes(xmin = .data$beta_ZX - 1.96 * .data$se_ZX,
                         xmax = .data$beta_ZX + 1.96 * .data$se_ZX),
                     height = 0, color = "grey60", alpha = 0.6) +
      geom_point(size = 3, color = MEDUSA_COLORS$ohdsiBlue) +
      labs(x = expression(beta[ZX]~"(SNP-Exposure)"),
           y = expression(beta[ZY]~"(SNP-Outcome)"),
           title = "Per-SNP Effect Estimates") +
      ggTheme()

    # Add method slopes
    if (!is.null(s) && !is.null(s$summary)) {
      for (i in seq_len(nrow(s$summary))) {
        meth <- s$summary$method[i]
        slope <- s$summary$beta_MR[i]
        if (is.na(slope)) next
        clr <- MEDUSA_COLORS$methods[meth]
        if (is.na(clr)) clr <- MEDUSA_COLORS$neutral

        if (meth == "MR-Egger" && !is.null(s$mrEgger)) {
          p <- p + geom_abline(intercept = s$mrEgger$intercept,
                               slope = slope, color = clr,
                               linetype = "dashed", linewidth = 0.7)
        } else {
          p <- p + geom_abline(intercept = 0, slope = slope,
                               color = clr, linetype = "dashed",
                               linewidth = 0.7)
        }
      }
    }

    p
  })

  # -- Funnel Plot --
  output$funnelPlot <- renderPlot({
    ps <- perSnp()
    s <- sens()
    if (is.null(ps)) return(NULL)

    waldRatio <- ps$beta_ZY / ps$beta_ZX
    waldSE <- ps$se_ZY / abs(ps$beta_ZX)
    precision <- 1 / waldSE

    funnelData <- data.frame(
      ratio = waldRatio,
      precision = precision,
      snp_id = ps$snp_id
    )

    p <- ggplot(funnelData, aes(x = .data$ratio, y = .data$precision)) +
      geom_point(size = 3, color = MEDUSA_COLORS$ohdsiBlue) +
      labs(x = expression(beta[MR]~"(Wald ratio)"),
           y = "Precision (1/SE)",
           title = "Funnel Plot") +
      ggTheme()

    # Add IVW estimate line
    if (!is.null(s) && !is.null(s$ivw) && !is.na(s$ivw$beta_MR)) {
      p <- p + geom_vline(xintercept = s$ivw$beta_MR,
                          linetype = "dashed", color = MEDUSA_COLORS$primary)
    }
    p
  })

  # -- Leave-One-Out --
  output$looPlot <- renderPlot({
    s <- sens()
    if (is.null(s) || is.null(s$leaveOneOut)) return(NULL)
    loo <- s$leaveOneOut
    loo$snp_removed <- factor(loo$snp_removed,
                               levels = rev(loo$snp_removed))
    loo$ci_lower <- loo$beta_MR - 1.96 * loo$se_MR
    loo$ci_upper <- loo$beta_MR + 1.96 * loo$se_MR

    p <- ggplot(loo, aes(x = .data$beta_MR, y = .data$snp_removed,
                          xmin = .data$ci_lower, xmax = .data$ci_upper)) +
      geom_vline(xintercept = 0, linetype = "dotted", color = "grey70") +
      geom_errorbarh(height = 0.2, linewidth = 0.6,
                     color = MEDUSA_COLORS$ohdsiBlue) +
      geom_point(size = 2.5, color = MEDUSA_COLORS$ohdsiBlue) +
      labs(x = expression(beta[MR]~"(IVW excluding SNP)"),
           y = "SNP Removed",
           title = "Leave-One-Out Analysis") +
      ggTheme()

    # Add full IVW estimate
    if (!is.null(s$ivw) && !is.na(s$ivw$beta_MR)) {
      p <- p + geom_vline(xintercept = s$ivw$beta_MR,
                          linetype = "dashed", color = MEDUSA_COLORS$primary,
                          linewidth = 0.7)
    }
    p
  })

  # ========================================================================
  # NEGATIVE CONTROLS
  # ========================================================================
  output$ncSummaryCards <- renderUI({
    nc <- ncRes()
    if (is.null(nc)) {
      return(tags$p(class = "text-muted",
                    "No negative control results available.",
                    " Provide negativeControlResults to launchResultsExplorer()."))
    }
    nNC <- nrow(nc$ncEstimates)
    biasLabel <- if (nc$biasDetected) "Bias Detected" else "No Bias Detected"
    biasTheme <- if (nc$biasDetected) "danger" else "success"
    biasIcon <- if (nc$biasDetected) "exclamation-triangle" else "check-circle"

    calLabel <- if (!is.null(nc$calibratedPrimary)) {
      fmtP(nc$calibratedPrimary$calibratedP)
    } else {
      "Not available"
    }

    tagList(
      statusCard("NC Outcomes Tested", as.character(nNC), "info", "list-ol"),
      statusCard("Systematic Bias", biasLabel, biasTheme, biasIcon),
      statusCard("Calibrated P-value", calLabel, "primary", "bullseye")
    )
  })

  output$ncEffectPlot <- renderPlot({
    nc <- ncRes()
    if (is.null(nc) || is.null(nc$ncEstimates) || nrow(nc$ncEstimates) == 0) {
      return(NULL)
    }
    ncEst <- nc$ncEstimates

    p <- ggplot(ncEst, aes(x = .data$beta_MR, y = seq_len(nrow(ncEst)))) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      geom_errorbarh(aes(xmin = .data$beta_MR - 1.96 * .data$se_MR,
                          xmax = .data$beta_MR + 1.96 * .data$se_MR),
                     height = 0.3, color = MEDUSA_COLORS$neutral) +
      geom_point(size = 2.5, color = MEDUSA_COLORS$ohdsiBlue) +
      scale_y_continuous(breaks = seq_len(nrow(ncEst)),
                         labels = ncEst$outcome_id) +
      labs(x = expression(beta[MR]~"(Wald ratio)"),
           y = "Negative Control Outcome",
           title = "NC Effect Estimates (Expected: Null)") +
      ggTheme()

    # Add primary estimate if available
    est <- mrEst()
    if (!is.null(est)) {
      p <- p + geom_vline(xintercept = est$betaMR,
                          linetype = "dotted", color = MEDUSA_COLORS$orange,
                          linewidth = 0.8)
    }
    p
  })

  output$ncCalibratedCard <- renderUI({
    nc <- ncRes()
    est <- mrEst()
    if (is.null(nc)) return(NULL)

    items <- list()

    if (nc$biasDetected) {
      items <- c(items, list(
        tags$div(class = "alert alert-warning mb-2",
                 bsicons::bs_icon("exclamation-triangle-fill"),
                 " Systematic bias detected in negative control distribution.",
                 " Interpret primary results with caution.")
      ))
    } else {
      items <- c(items, list(
        tags$div(class = "alert alert-success mb-2",
                 bsicons::bs_icon("check-circle-fill"),
                 " No systematic bias detected in negative controls.")
      ))
    }

    if (!is.null(nc$calibratedPrimary) && !is.null(est)) {
      cal <- nc$calibratedPrimary
      items <- c(items, list(
        tags$table(class = "table table-sm table-bordered mt-2",
          tags$thead(tags$tr(
            tags$th(""), tags$th("Uncalibrated"), tags$th("Calibrated")
          )),
          tags$tbody(
            tags$tr(
              tags$td("P-value"),
              tags$td(fmtP(est$pValue)),
              tags$td(fmtP(cal$calibratedP))
            ),
            tags$tr(
              tags$td("95% CI Lower"),
              tags$td(sprintf("%.4f", est$ciLower)),
              tags$td(if (!is.null(cal$calibratedCiLower))
                sprintf("%.4f", cal$calibratedCiLower) else "NA")
            ),
            tags$tr(
              tags$td("95% CI Upper"),
              tags$td(sprintf("%.4f", est$ciUpper)),
              tags$td(if (!is.null(cal$calibratedCiUpper))
                sprintf("%.4f", cal$calibratedCiUpper) else "NA")
            )
          )
        )
      ))
    }

    tagList(items)
  })

  output$ncEstimatesTable <- renderTable({
    nc <- ncRes()
    if (is.null(nc) || is.null(nc$ncEstimates)) return(NULL)
    df <- nc$ncEstimates
    df$beta_MR <- round(df$beta_MR, 4)
    df$se_MR <- round(df$se_MR, 4)
    df$pval <- sapply(df$pval, fmtP)
    df$OR <- round(exp(nc$ncEstimates$beta_MR), 3)
    df[, c("outcome_id", "beta_MR", "se_MR", "pval", "OR")]
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
}


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
shinyApp(ui = ui, server = server)
