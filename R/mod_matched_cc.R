# mod_matched_cc.R — Matched Case-Control Data module
#
# Focuses on 1:1 McNemar analysis.  The user enters the 2×2 discordant-pairs
# table: n11 (both exposed), n10 (case exposed / control unexposed),
#         n01 (case unexposed / control exposed), n00 (both unexposed).
#
# McNemar OR  = n10 / n01
# SE(ln OR)   = sqrt(1/n10 + 1/n01)
# 95% CI      = exp(ln(OR) ± 1.96 × SE)
# 99% CI      = exp(ln(OR) ± 2.576 × SE)
# chi²        = (|n10 - n01| - 1)² / (n10 + n01)  [with continuity correction]
# p-value     = pchisq(chi², df = 1, lower.tail = FALSE)
#
# Crude OR uses marginal totals derived from the 2×2 pairs table.
#
# Test values (CLAUDE.md, matching-ratio 1):
#   Case Exposed=41/23/16, Case Unexposed=31/8/0 → RRmh=2.0615
#   → these represent 1:1 pairs where columns are 0/1/2+ controls exposed.
#   For MVP the user enters the standard discordant pairs (n10, n01, n11, n00).
#
# Reference: Modern Epidemiology 3rd ed., Ch. 8 & 15.

`%||%` <- function(a, b) if (!is.null(a)) a else b

matchedCcUI <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(class = "calc-layout",

      # ── Input panel ──────────────────────────────────────────────────────────
      tags$div(class = "panel",
        tags$div(class = "panel-head ph-teal", "Input"),
        tags$div(class = "panel-body",

          tags$p(
            style = "font-size: 11px; color: #4a6070; margin-bottom: 12px;",
            "Enter the 2\u00d72 table of matched pairs (1:1 matching). Each cell ",
            "represents a pair type defined by whether the case and the matched ",
            "control were exposed."
          ),

          # 2×2 matched-pairs table
          tags$table(class = "tbl",
            style = "margin-bottom: 14px;",
            tags$thead(
              tags$tr(
                tags$th(class = "rh", ""),
                tags$th("Control Exposed"),
                tags$th("Control Unexposed"),
                tags$th("Total")
              )
            ),
            tags$tbody(
              # Row 1: Case Exposed
              tags$tr(
                tags$td(class = "cl", "Case Exposed"),
                tags$td(
                  style = "position: relative;",
                  tags$span(class = "cell-letter", "n\u2081\u2081"),
                  numericInput(ns("n11"), label = NULL, value = NULL, min = 0, step = 1)
                ),
                tags$td(
                  style = "position: relative;",
                  tags$span(class = "cell-letter", "n\u2081\u2080"),
                  numericInput(ns("n10"), label = NULL, value = NULL, min = 0, step = 1)
                ),
                tags$td(class = "tot", uiOutput(ns("row1_tot")))
              ),
              # Row 2: Case Unexposed
              tags$tr(
                tags$td(class = "cl", "Case Unexposed"),
                tags$td(
                  style = "position: relative;",
                  tags$span(class = "cell-letter", "n\u2080\u2081"),
                  numericInput(ns("n01"), label = NULL, value = NULL, min = 0, step = 1)
                ),
                tags$td(
                  style = "position: relative;",
                  tags$span(class = "cell-letter", "n\u2080\u2080"),
                  numericInput(ns("n00"), label = NULL, value = NULL, min = 0, step = 1)
                ),
                tags$td(class = "tot", uiOutput(ns("row2_tot")))
              ),
              # Total row
              tags$tr(
                tags$td(class = "cl", "Total pairs"),
                tags$td(class = "tot", uiOutput(ns("col1_tot"))),
                tags$td(class = "tot", uiOutput(ns("col2_tot"))),
                tags$td(class = "tot", uiOutput(ns("grand_tot")))
              )
            )
          ),

          # Cell annotation key
          tags$p(
            style = "font-size: 10px; color: #4a6070; margin-bottom: 14px;",
            tags$strong("n\u2081\u2081"), ": both exposed; ",
            tags$strong("n\u2081\u2080"), ": case exposed, control unexposed (discordant); ",
            tags$strong("n\u2080\u2081"), ": case unexposed, control exposed (discordant); ",
            tags$strong("n\u2080\u2080"), ": both unexposed."
          ),

          tags$button(
            class   = "calc-btn",
            onclick = sprintf("Shiny.setInputValue('%s', Math.random())", ns("calc")),
            "Calculate"
          )
        )
      ),

      # ── Results panel ─────────────────────────────────────────────────────────
      tags$div(class = "panel",
        tags$div(class = "panel-head ph-blue", "Results"),
        tags$div(class = "panel-body",
          uiOutput(ns("results"))
        )
      )
    ),

    # ── Interpretation ─────────────────────────────────────────────────────────
    uiOutput(ns("interpretation"))
  )
}

matchedCcServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ── Live margin totals ─────────────────────────────────────────────────────
    get_cell <- function(name) {
      v <- input[[name]]
      if (is.null(v) || is.na(v)) 0L else as.integer(v)
    }

    output$row1_tot  <- renderUI({ n11 <- get_cell("n11"); n10 <- get_cell("n10"); tags$span(n11 + n10) })
    output$row2_tot  <- renderUI({ n01 <- get_cell("n01"); n00 <- get_cell("n00"); tags$span(n01 + n00) })
    output$col1_tot  <- renderUI({ n11 <- get_cell("n11"); n01 <- get_cell("n01"); tags$span(n11 + n01) })
    output$col2_tot  <- renderUI({ n10 <- get_cell("n10"); n00 <- get_cell("n00"); tags$span(n10 + n00) })
    output$grand_tot <- renderUI({
      tags$span(get_cell("n11") + get_cell("n10") + get_cell("n01") + get_cell("n00"))
    })

    # ── Core computation ───────────────────────────────────────────────────────
    computed <- eventReactive(input$calc, {
      n11 <- input$n11; n10 <- input$n10
      n01 <- input$n01; n00 <- input$n00

      validate(
        need(!is.null(n11) && !is.null(n10) && !is.null(n01) && !is.null(n00),
             "Enter values for all four cells of the matched-pairs table."),
        need(!anyNA(c(n11, n10, n01, n00)),
             "All four cell counts must be numeric."),
        need(all(c(n11, n10, n01, n00) >= 0),
             "Cell counts must be non-negative."),
        need(n10 + n01 > 0,
             "Both discordant cells (n10 and n01) cannot be zero simultaneously.")
      )

      n11 <- as.numeric(n11); n10 <- as.numeric(n10)
      n01 <- as.numeric(n01); n00 <- as.numeric(n00)

      N <- n11 + n10 + n01 + n00

      # ── McNemar odds ratio ────────────────────────────────────────────────────
      validate(need(n10 > 0 && n01 > 0,
        "Both discordant cells must be > 0 to compute the McNemar OR and its CI. ",
        "If only one is zero, an exact method is needed."))

      OR_mc  <- n10 / n01
      lnOR   <- log(OR_mc)
      SE_ln  <- sqrt(1 / n10 + 1 / n01)

      ci95   <- exp(lnOR + c(-1, 1) * 1.96  * SE_ln)
      ci99   <- exp(lnOR + c(-1, 1) * 2.576 * SE_ln)
      ci90   <- exp(lnOR + c(-1, 1) * 1.645 * SE_ln)

      # ── McNemar chi-squared (with continuity correction) ──────────────────────
      chi2  <- (abs(n10 - n01) - 1)^2 / (n10 + n01)
      p_val <- pchisq(chi2, df = 1, lower.tail = FALSE)

      # McNemar chi-squared without continuity correction (for display)
      chi2_nc  <- (n10 - n01)^2 / (n10 + n01)
      p_val_nc <- pchisq(chi2_nc, df = 1, lower.tail = FALSE)

      # ── Crude OR from marginal 2×2 reconstructed from pair data ───────────────
      # From matched pairs: a = exposed cases = n11 + n10
      #                     b = unexposed cases = n01 + n00
      #                     c = exposed controls = n11 + n01
      #                     d = unexposed controls = n10 + n00
      a <- n11 + n10   # cases exposed
      b <- n01 + n00   # cases unexposed
      c <- n11 + n01   # controls exposed
      d <- n10 + n00   # controls unexposed

      crude_OR <- NA_real_
      crude_ci95 <- c(NA_real_, NA_real_)
      if (b > 0 && c > 0 && a > 0 && d > 0) {
        crude_OR  <- (a * d) / (b * c)
        SE_crude  <- sqrt(1/a + 1/b + 1/c + 1/d)
        crude_ci95 <- exp(log(crude_OR) + c(-1, 1) * 1.96 * SE_crude)
      }

      list(
        n11 = n11, n10 = n10, n01 = n01, n00 = n00, N = N,
        OR_mc    = OR_mc,
        ci90     = ci90,
        ci95     = ci95,
        ci99     = ci99,
        chi2     = chi2,
        chi2_nc  = chi2_nc,
        p_val    = p_val,
        p_val_nc = p_val_nc,
        crude_OR   = crude_OR,
        crude_ci95 = crude_ci95
      )
    })

    # ── Results output ─────────────────────────────────────────────────────────
    output$results <- renderUI({
      r <- computed()

      crude_txt <- if (!is.na(r$crude_OR)) {
        paste0(
          fmt4(r$crude_OR),
          "  (95% CI: ", fmt4(r$crude_ci95[1]), "\u2013", fmt4(r$crude_ci95[2]), ")"
        )
      } else {
        "\u2014"
      }

      chi2_txt <- paste0(
        "Continuity-corrected: \u03c7\u00b2 = ", fmt4(r$chi2), ",  p = ", fmt4(r$p_val),
        "\nUncorrected: \u03c7\u00b2 = ", fmt4(r$chi2_nc), ",  p = ", fmt4(r$p_val_nc)
      )

      tags$div(class = "metrics",
        metric_card(
          "McNemar Odds Ratio",
          fmt4(r$OR_mc),
          paste0("n10 = ", r$n10, ", n01 = ", r$n01,
                 " (", r$n10 + r$n01, " discordant pairs)"),
          "m-primary wide",
          tip = "McNemar odds ratio from discordant pairs: n\u2081\u2080 / n\u2080\u2081. Only discordant pairs contribute information about the exposure-disease association."
        ),
        metric_card(
          "95% Confidence Interval",
          paste0(fmt4(r$ci95[1]), "\u2013", fmt4(r$ci95[2])),
          paste0("90% CI: ", fmt4(r$ci90[1]), "\u2013", fmt4(r$ci90[2])),
          "m-blue",
          tip = "Normal-approximation confidence interval for the McNemar OR, based on SE(ln OR) = sqrt(1/n10 + 1/n01)."
        ),
        metric_card(
          "99% Confidence Interval",
          paste0(fmt4(r$ci99[1]), "\u2013", fmt4(r$ci99[2])),
          "Normal approximation on ln(OR)",
          "m-teal",
          tip = "99% normal-approximation confidence interval for the McNemar OR."
        ),
        metric_card(
          "McNemar P-value",
          fmt4(r$p_val),
          paste0("\u03c7\u00b2 (corrected) = ", fmt4(r$chi2), " (df = 1)"),
          "m-purple",
          tip = "McNemar chi-square p-value (with continuity correction) testing the null hypothesis that the OR = 1."
        ),
        metric_card(
          "Crude OR (marginal)",
          crude_txt,
          "From reconstructed marginal 2\u00d72",
          "m-blue",
          tip = "Unadjusted odds ratio from the marginal 2\u00d72 table reconstructed from the matched-pairs counts, ignoring the matching."
        ),
        metric_card(
          "Total pairs (N)",
          as.character(r$N),
          paste0("Concordant: ", r$n11 + r$n00,
                 " \u2502 Discordant: ", r$n10 + r$n01),
          "m-teal",
          tip = "Total number of matched pairs. Only discordant pairs (n10, n01) contribute to the McNemar estimate."
        )
      )
    })

    # ── Interpretation ─────────────────────────────────────────────────────────
    output$interpretation <- renderUI({
      r   <- computed()
      sig <- r$p_val < 0.05
      dir <- if (r$OR_mc > 1) "higher" else "lower"

      tags$div(class = "conclusion",
        tags$div(class = "conc-head", "Interpretation"),
        tags$p(class = "conc-text",
          "In this 1:1 matched case-control study with ",
          tags$strong(as.character(r$N)), " pairs (",
          tags$strong(as.character(as.integer(r$n10 + r$n01))),
          " discordant), the McNemar odds ratio is ",
          tags$strong(fmt4(r$OR_mc)),
          " (95% CI: ", fmt4(r$ci95[1]), "\u2013", fmt4(r$ci95[2]),
          "; 99% CI: ", fmt4(r$ci99[1]), "\u2013", fmt4(r$ci99[2]), "). ",
          "The odds of exposure among cases is ", tags$strong(dir), " than among controls. ",
          "The McNemar test gives \u03c7\u00b2\u202f=\u202f", tags$strong(fmt4(r$chi2)),
          " (continuity-corrected, p\u202f=\u202f", tags$strong(fmt4(r$p_val)), "). ",
          if (sig)
            "There is a statistically significant association between exposure and case status (p < 0.05)."
          else
            "There is no statistically significant association between exposure and case status (p \u2265 0.05).",
          if (!is.na(r$crude_OR)) {
            tagList(
              " The crude OR from the marginal totals is ",
              tags$strong(fmt4(r$crude_OR)),
              " (95% CI: ", fmt4(r$crude_ci95[1]), "\u2013", fmt4(r$crude_ci95[2]), ")."
            )
          }
        )
      )
    })

    outputOptions(output, "results",        suspendWhenHidden = FALSE)
    outputOptions(output, "interpretation", suspendWhenHidden = FALSE)
  })
}
