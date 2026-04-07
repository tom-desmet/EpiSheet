# mod_study_size.R — Study Size (precision-based) module
# Precision-based sample size: works backward from desired CI width to required N.
# Method: SE(ln OR/RR) = sqrt(1/a + 1/b + 1/c + 1/d); Upper CI = exp(ln(est) + z*SE);
#   Pr = pnorm((ln(bound) - ln(est)) / SE) — probability upper CI ≤ desired bound.
# Reference: Rothman KJ, Modern Epidemiology.

studySizeUI <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(class = "calc-layout",

      # ── Input panel ────────────────────────────────────────────────────────────
      tags$div(class = "panel",
        tags$div(class = "panel-head ph-teal", "Input"),
        tags$div(class = "panel-body",

          tags$p(class = "panel-hint",
            "Precision-based study size: find N such that the upper confidence ",
            "limit is likely to fall below a desired bound."),

          tags$div(class = "input-group",
            selectInput(ns("study_type"), "Study type",
                        choices = c("Cohort" = "cohort",
                                    "Case-control" = "cc"),
                        selected = "cohort")
          ),

          tags$div(class = "input-group",
            numericInput(ns("conf_level"), "Confidence level (%)",
                         value = 95, min = 80, max = 99.9, step = 1)
          ),

          tags$div(class = "input-group",
            numericInput(ns("bound"), "Upper CI bound < (OR/RR)",
                         value = 2.0, min = 1.001, max = 1000, step = 0.1)
          ),

          tags$div(class = "input-group",
            numericInput(ns("est"), "Expected OR/RR >",
                         value = 1.0, min = 0.001, max = 999, step = 0.1)
          ),

          tags$div(class = "input-group",
            numericInput(ns("p0"), "Exposure prevalence (p\u2080)",
                         value = 0.3, min = 1e-4, max = 0.9999, step = 0.05)
          ),

          tags$div(class = "input-group",
            numericInput(ns("alloc"), "Allocation ratio (N\u2081/N\u2082)",
                         value = 1, min = 0.01, max = 100, step = 0.1)
          ),

          tags$div(class = "input-group",
            numericInput(ns("max_size"), "Maximum total size",
                         value = 10000, min = 100, max = 1e7, step = 100)
          ),

          tags$button(class = "calc-btn",
            onclick = sprintf("Shiny.setInputValue('%s', Math.random())",
                              ns("calc")),
            "Calculate")
        )
      ),

      # ── Results panel ──────────────────────────────────────────────────────────
      tags$div(class = "panel",
        tags$div(class = "panel-head ph-blue", "Results"),
        tags$div(class = "panel-body",
          uiOutput(ns("results"))
        )
      )
    ),

    uiOutput(ns("interpretation"))
  )
}

studySizeServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    computed <- eventReactive(input$calc, {
      validate(
        need(!is.na(input$conf_level) && input$conf_level > 50 && input$conf_level < 100,
             "Confidence level must be between 50 and 100."),
        need(!is.na(input$bound) && input$bound > 1,
             "Upper CI bound must be greater than 1."),
        need(!is.na(input$est) && input$est > 0,
             "Expected OR/RR must be positive."),
        need(input$est < input$bound,
             "Expected OR/RR must be less than the upper CI bound."),
        need(!is.na(input$p0) && input$p0 > 0 && input$p0 < 1,
             "Exposure prevalence must be between 0 and 1."),
        need(!is.na(input$alloc) && input$alloc > 0,
             "Allocation ratio must be positive."),
        need(!is.na(input$max_size) && input$max_size >= 100,
             "Maximum total size must be at least 100.")
      )

      alpha     <- 1 - input$conf_level / 100
      z         <- qnorm(1 - alpha / 2)
      bound     <- input$bound
      est       <- input$est
      p0        <- input$p0
      r         <- input$alloc       # N1/N2
      max_size  <- input$max_size
      study_type <- input$study_type

      ln_est   <- log(est)
      ln_bound <- log(bound)

      # Build grid of total N values — ~80 steps on log scale for smooth curve
      n_seq <- unique(round(exp(seq(log(max(10, max_size / 200)),
                                    log(max_size),
                                    length.out = 80))))
      n_seq <- n_seq[n_seq >= 4]  # need at least 4 for 2x2 cells

      compute_row <- function(N_total) {
        # Split N_total into N1 (exposed/cases) and N2 (unexposed/controls)
        N1 <- N_total * r / (1 + r)
        N2 <- N_total / (1 + r)

        if (study_type == "cohort") {
          # Cohort: p1 = RR * p0, proportion exposed developing disease = p1
          p1 <- est * p0
          if (p1 >= 1) return(NULL)
          a  <- N1 * p1          # exposed cases
          b  <- N1 * (1 - p1)   # exposed non-cases
          c  <- N2 * p0          # unexposed cases
          d  <- N2 * (1 - p0)   # unexposed non-cases
        } else {
          # Case-control: p1 derived from OR and p0
          p1 <- (est * p0) / (1 + p0 * (est - 1))
          if (p1 >= 1 || p1 <= 0) return(NULL)
          # N1 = cases, N2 = controls; proportion exposed among cases = p1, controls = p0
          a  <- N1 * p1          # cases exposed
          b  <- N1 * (1 - p1)   # cases unexposed
          c  <- N2 * p0          # controls exposed
          d  <- N2 * (1 - p0)   # controls unexposed
        }

        if (any(c(a, b, c, d) < 0.5)) return(NULL)

        se    <- sqrt(1/a + 1/b + 1/c + 1/d)
        span  <- 2 * z * se
        upper_ci <- exp(ln_est + z * se)
        # Pr: probability that a study of this size yields upper CI ≤ bound
        # Treating ln(upper CI) ~ Normal(ln(est) + z*SE, ...) is deterministic for fixed N.
        # Proper interpretation: Pr = P(Z ≤ (ln(bound) - ln(est)) / se - z)
        # = pnorm((ln_bound - ln_est) / se - z)
        # This is derived from: upper CI = exp(ln(est) + z * se_hat),
        # where se_hat is the realised SE — approximated by its expected value here.
        pr <- pnorm((ln_bound - ln_est) / se - z)

        data.frame(
          N_total  = N_total,
          N1       = round(N1),
          N2       = round(N2),
          SE_lnOR  = se,
          Span     = span,
          Upper_CI = upper_ci,
          Pr       = pr
        )
      }

      rows <- lapply(n_seq, compute_row)
      rows <- rows[!sapply(rows, is.null)]
      validate(need(length(rows) > 0,
        "No valid sample sizes could be computed. Check inputs."))

      df <- do.call(rbind, rows)

      # Find N where Pr first crosses 0.50 and 0.80 (if it does)
      idx_50 <- which(df$Pr >= 0.50)
      idx_80 <- which(df$Pr >= 0.80)
      n_50 <- if (length(idx_50) > 0) df$N_total[idx_50[1]] else NA_integer_
      n_80 <- if (length(idx_80) > 0) df$N_total[idx_80[1]] else NA_integer_

      # Display table: ~10 evenly spaced rows
      tbl_idx <- unique(round(seq(1, nrow(df), length.out = min(12, nrow(df)))))
      tbl <- df[tbl_idx, ]
      tbl$SE_lnOR  <- round(tbl$SE_lnOR, 4)
      tbl$Span     <- round(tbl$Span, 4)
      tbl$Upper_CI <- round(tbl$Upper_CI, 4)
      tbl$Pr       <- round(tbl$Pr, 4)

      col_lbl <- if (study_type == "cohort") c("N\u2081 (exposed)", "N\u2082 (unexposed)")
                 else c("N\u2081 (cases)", "N\u2082 (controls)")

      list(df = df, tbl = tbl, z = z,
           n_50 = n_50, n_80 = n_80,
           bound = bound, est = est,
           col_lbl = col_lbl,
           study_type = study_type,
           conf_level = input$conf_level)
    })

    output$results <- renderUI({
      r <- computed()

      metric_cards <- tags$div(class = "metrics",
        tags$div(class = "metric m-primary wide",
          tags$div(class = "metric-label", "N for Pr \u2265 50%"),
          tags$div(class = "metric-value",
                   if (!is.na(r$n_50)) formatC(r$n_50, format = "d", big.mark = ",")
                   else "\u2014"),
          tags$div(class = "metric-sub",
                   paste0("Upper CI < ", fmt4(r$bound), " with \u226550% probability"))
        ),
        tags$div(class = "metric m-blue",
          tags$div(class = "metric-label", "N for Pr \u2265 80%"),
          tags$div(class = "metric-value",
                   if (!is.na(r$n_80)) formatC(r$n_80, format = "d", big.mark = ",")
                   else "\u2014"),
          tags$div(class = "metric-sub", "80% probability threshold")
        ),
        tags$div(class = "metric m-teal",
          tags$div(class = "metric-label", "Confidence level"),
          tags$div(class = "metric-value", paste0(r$conf_level, "%")),
          tags$div(class = "metric-sub",
                   paste0("z = ", fmt4(r$z)))
        )
      )

      tagList(
        metric_cards,
        tags$div(style = "margin-top: 16px;",
          plotOutput(ns("size_chart"), height = "300px")
        ),
        tags$h4(style = "margin-top: 14px; font-size: 12px; font-weight: 500; color: #4a6070;",
                "Sample size table"),
        tableOutput(ns("size_table"))
      )
    })

    output$size_chart <- renderPlot({
      r <- computed()
      df <- r$df

      ggplot2::ggplot(df, ggplot2::aes(x = N_total, y = Pr)) +
        ggplot2::geom_line(colour = "#005f6a", linewidth = 1.3) +
        ggplot2::geom_hline(yintercept = 0.80, linetype = "dashed",
                            colour = "#165c7d", linewidth = 0.7, alpha = 0.7) +
        ggplot2::geom_hline(yintercept = 0.50, linetype = "dotted",
                            colour = "#4a9aa5", linewidth = 0.6, alpha = 0.7) +
        ggplot2::annotate("text", x = max(df$N_total) * 0.98, y = 0.82,
                          label = "80%", hjust = 1, size = 3.2,
                          colour = "#165c7d") +
        ggplot2::annotate("text", x = max(df$N_total) * 0.98, y = 0.52,
                          label = "50%", hjust = 1, size = 3.2,
                          colour = "#4a9aa5") +
        ggplot2::scale_x_continuous(
          labels = scales::comma,
          expand = ggplot2::expansion(mult = c(0, 0.02))
        ) +
        ggplot2::scale_y_continuous(
          limits = c(0, 1),
          breaks = seq(0, 1, 0.2),
          labels = scales::percent
        ) +
        ggplot2::labs(
          title = paste0("Precision: Pr(upper CI \u2264 ", fmt4(r$bound), ")"),
          x     = "Total study size",
          y     = paste0("Pr(upper CI \u2264 bound)")
        ) +
        ggplot2::theme_bw(base_size = 13) +
        ggplot2::theme(
          panel.grid.minor  = ggplot2::element_blank(),
          plot.title        = ggplot2::element_text(size = 12, colour = "#005f6a"),
          axis.title        = ggplot2::element_text(size = 11)
        )
    })

    output$size_table <- renderTable({
      r    <- computed()
      tbl  <- r$tbl
      lbls <- r$col_lbl
      # Rename columns for display
      names(tbl) <- c("Total N", lbls[1], lbls[2],
                      "SE(ln OR)", "Span (2\u00b7z\u00b7SE)", "Upper CI", "Pr")
      tbl
    }, digits = 4, striped = TRUE, hover = TRUE, bordered = TRUE)

    output$interpretation <- renderUI({
      r    <- computed()
      type_str <- if (r$study_type == "cohort") "cohort" else "case-control"
      n50_str  <- if (!is.na(r$n_50))
        formatC(r$n_50, format = "d", big.mark = ",") else "not reached"
      n80_str  <- if (!is.na(r$n_80))
        formatC(r$n_80, format = "d", big.mark = ",") else "not reached"

      tags$div(class = "conclusion",
        tags$div(class = "conc-head", "Interpretation"),
        tags$p(class = "conc-text",
          "For a ", tags$strong(type_str), " study with an expected ",
          "OR/RR of ", tags$strong(fmt4(r$est)),
          " and a desired upper ", tags$strong(paste0(r$conf_level, "%")),
          " confidence limit below ", tags$strong(fmt4(r$bound)),
          ": a total sample size of approximately ", tags$strong(n50_str),
          " is needed for a 50% probability of achieving this precision",
          " (the \u2018median\u2019 required size). ",
          "To reach 80% probability, approximately ", tags$strong(n80_str),
          " subjects are needed. ",
          "These are precision-based estimates, not power calculations — ",
          "they reflect how large a study must be to observe an upper confidence ",
          "limit no greater than ", tags$strong(fmt4(r$bound)),
          " with the specified probability."
        )
      )
    })

    outputOptions(output, "results",        suspendWhenHidden = FALSE)
    outputOptions(output, "interpretation", suspendWhenHidden = FALSE)
  })
}
