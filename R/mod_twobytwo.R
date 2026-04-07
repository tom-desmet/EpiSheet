# mod_twobytwo.R — 2×2 Table module

twobytwoUI <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(class = "calc-layout",

      # Input panel
      tags$div(class = "panel",
        tags$div(class = "panel-head ph-teal", "Input"),
        tags$div(class = "panel-body",
          tags$table(class = "tbl",
            tags$thead(
              tags$tr(
                tags$th(class = "rh", ""),
                tags$th("Cases"),
                tags$th("Non-cases"),
                tags$th("Total")
              )
            ),
            tags$tbody(
              tags$tr(
                tags$td(class = "cl", "Exposed"),
                tags$td(numericInput(ns("ca"), tags$span(class = "cell-lbl", "a"), value = NA, min = 0, width = "80px")),
                tags$td(numericInput(ns("cb"), tags$span(class = "cell-lbl", "b"), value = NA, min = 0, width = "80px")),
                tags$td(class = "tot", uiOutput(ns("n1"), inline = TRUE))
              ),
              tags$tr(
                tags$td(class = "cl", "Unexposed"),
                tags$td(numericInput(ns("cc"), tags$span(class = "cell-lbl", "c"), value = NA, min = 0, width = "80px")),
                tags$td(numericInput(ns("cd"), tags$span(class = "cell-lbl", "d"), value = NA, min = 0, width = "80px")),
                tags$td(class = "tot", uiOutput(ns("n0"), inline = TRUE))
              ),
              tags$tr(
                tags$td(class = "cl", "Total"),
                tags$td(class = "tot", uiOutput(ns("m1"), inline = TRUE)),
                tags$td(class = "tot", uiOutput(ns("m0"), inline = TRUE)),
                tags$td(class = "tot", uiOutput(ns("N"),  inline = TRUE))
              )
            )
          ),
          tags$button(
            class   = "calc-btn",
            onclick = sprintf("Shiny.setInputValue('%s', Math.random())", ns("calc")),
            "Calculate"
          )
        )
      ),

      # Results panel
      tags$div(class = "panel",
        tags$div(class = "panel-head ph-blue", "Results"),
        tags$div(class = "panel-body",
          uiOutput(ns("results"))
        )
      )
    ),

    # Interpretation block
    uiOutput(ns("interpretation"))
  )
}

twobytwoServer <- function(id) {
  moduleServer(id, function(input, output, session) {

    # Marginal totals
    output$n1 <- renderUI({
      if (!is.na(input$ca) && !is.na(input$cb)) input$ca + input$cb else "\u2014"
    })
    output$n0 <- renderUI({
      if (!is.na(input$cc) && !is.na(input$cd)) input$cc + input$cd else "\u2014"
    })
    output$m1 <- renderUI({
      if (!is.na(input$ca) && !is.na(input$cc)) input$ca + input$cc else "\u2014"
    })
    output$m0 <- renderUI({
      if (!is.na(input$cb) && !is.na(input$cd)) input$cb + input$cd else "\u2014"
    })
    output$N <- renderUI({
      vals <- c(input$ca, input$cb, input$cc, input$cd)
      if (!any(is.na(vals))) sum(vals) else "\u2014"
    })

    # Core computation triggered by Calculate button
    computed <- eventReactive(input$calc, {
      ca <- input$ca; cb <- input$cb
      cc <- input$cc; cd <- input$cd

      validate(
        need(!any(is.na(c(ca, cb, cc, cd))), "Enter all four cell counts."),
        need(all(c(ca, cb, cc, cd) >= 0),    "All counts must be >= 0.")
      )

      n1 <- ca + cb; n0 <- cc + cd
      m1 <- ca + cc; m0 <- cb + cd
      N  <- n1 + n0

      validate(need(n1 > 0 && n0 > 0, "Row totals cannot be zero."))
      validate(need(cb > 0 && cc > 0,  "Cannot compute: one or more cells is zero."))

      p1 <- ca / n1
      p0 <- cc / n0
      RR <- p1 / p0
      RD <- p1 - p0
      OR <- (ca * cd) / (cb * cc)

      # Taylor CIs
      se_ln_rr <- sqrt(cb / (ca * n1) + cd / (cc * n0))
      rr95 <- exp(log(RR) + c(-1, 1) * 1.96  * se_ln_rr)
      rr99 <- exp(log(RR) + c(-1, 1) * 2.576 * se_ln_rr)

      se_rd <- sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0)
      rd95 <- RD + c(-1, 1) * 1.96  * se_rd
      rd99 <- RD + c(-1, 1) * 2.576 * se_rd

      se_ln_or <- sqrt(1/ca + 1/cb + 1/cc + 1/cd)
      or95 <- exp(log(OR) + c(-1, 1) * 1.96  * se_ln_or)
      or99 <- exp(log(OR) + c(-1, 1) * 2.576 * se_ln_or)

      # Chi-square
      chi2    <- N * (ca*cd - cb*cc)^2 / (n1 * n0 * m1 * m0)
      chi2_cc <- N * (max(0, abs(ca*cd - cb*cc) - N/2))^2 / (n1 * n0 * m1 * m0)
      p_chi2    <- pchisq(chi2,    1, lower.tail = FALSE)
      p_chi2_cc <- pchisq(chi2_cc, 1, lower.tail = FALSE)
      p_fisher  <- fisher.test(matrix(c(ca, cc, cb, cd), 2, 2))$p.value

      # Attributable risk / NNT
      ar_exp <- (p1 - p0) / p1
      ar_pop <- ((ca + cc) / N - p0) / ((ca + cc) / N)
      nnt    <- 1 / abs(RD)

      list(
        RR = RR, rr95 = rr95, rr99 = rr99,
        RD = RD, rd95 = rd95, rd99 = rd99,
        OR = OR, or95 = or95, or99 = or99,
        chi2 = chi2, chi2_cc = chi2_cc,
        p_chi2 = p_chi2, p_chi2_cc = p_chi2_cc,
        p_fisher = p_fisher,
        ar_exp = ar_exp, ar_pop = ar_pop,
        nnt = nnt, RD_sign = sign(RD)
      )
    })

    fmt_p <- function(x) {
      if (is.na(x) || is.nan(x)) return("\u2014")
      if (x < 0.0001) return("< 0.0001")
      formatC(x, digits = 4, format = "f")
    }

    output$results <- renderUI({
      r <- computed()
      nnt_label <- ifelse(r$RD_sign > 0, "NNH", "NNT")
      tags$div(class = "metrics",
        metric_card("Risk Ratio (RR)", fmt2(r$RR),
           paste0("95% CI: ", fmt2(r$rr95[1]), "\u2013", fmt2(r$rr95[2])),
           "m-primary wide",
           "Ratio of disease risk in the exposed vs unexposed group. RR > 1 indicates increased risk."),
        metric_card("Odds Ratio (OR)", fmt2(r$OR),
           paste0("95% CI: ", fmt2(r$or95[1]), "\u2013", fmt2(r$or95[2])),
           "m-blue",
           "Ratio of the odds of disease in exposed vs unexposed. Approximates RR when disease is rare."),
        metric_card("Risk Difference", fmt2(r$RD),
           paste0("95% CI: ", fmt2(r$rd95[1]), "\u2013", fmt2(r$rd95[2])),
           "m-blue",
           "Absolute difference in disease probability between exposed and unexposed groups."),
        metric_card("\u03c7\u00b2 (no cc)", fmt2(r$chi2),    paste0("p = ", fmt_p(r$p_chi2)),    "m-teal",
           "Pearson chi-square test of association without continuity correction."),
        metric_card("\u03c7\u00b2 (cc)",    fmt2(r$chi2_cc), paste0("p = ", fmt_p(r$p_chi2_cc)), "m-teal",
           "Chi-square with Yates continuity correction. More conservative than uncorrected."),
        metric_card("Fisher exact p", fmt_p(r$p_fisher), cls = "m-purple",
           tip = "Exact p-value based on the hypergeometric distribution. Preferred for small cell counts."),
        metric_card("AR (exposed)",   fmt2(r$ar_exp),    cls = "m-purple",
           tip = "Attributable risk among the exposed: proportion of disease in the exposed group due to the exposure."),
        metric_card("AR (population)", fmt2(r$ar_pop),   cls = "m-blue",
           tip = "Population attributable risk: proportion of total disease in the population attributable to the exposure."),
        metric_card(nnt_label, fmt2(r$nnt), cls = "m-blue",
           tip = "Number needed to harm (NNH) or treat (NNT): how many people must be exposed to cause one additional case.")
      )
    })

    output$interpretation <- renderUI({
      r <- computed()
      dir   <- ifelse(r$RR > 1, "higher", "lower")
      pct   <- abs(round((r$RR - 1) * 100))
      p_ok  <- !is.na(r$p_chi2) && !is.nan(r$p_chi2)
      assoc <- if (p_ok && r$p_chi2 < 0.05) "statistically significant (p < 0.05)" else
               if (p_ok) "not statistically significant (p \u2265 0.05)" else
               "unknown (check inputs)"
      tags$div(class = "conclusion",
        tags$div(class = "conc-head", "Interpretation"),
        tags$p(class = "conc-text",
          "The risk of disease was ", tags$strong(paste0(pct, "% ", dir)),
          " in the exposed group (RR = ", tags$strong(fmt4(r$RR)),
          ", 95% CI: ", fmt4(r$rr95[1]), "\u2013", fmt4(r$rr95[2]), "). ",
          "The odds ratio was ", tags$strong(fmt4(r$OR)),
          " and the absolute risk difference was ", tags$strong(fmt4(r$RD)), ". ",
          "The association is ", tags$strong(assoc), "."
        )
      )
    })

    outputOptions(output, "results",        suspendWhenHidden = FALSE)
    outputOptions(output, "interpretation", suspendWhenHidden = FALSE)
  })
}
