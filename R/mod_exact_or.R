# mod_exact_or.R — O.R. Exact module
# Exact OR confidence limits via hypergeometric bisection (port of VBA Module2)
# Test: input 4/386/8/1636 → OR=3.2383; 95% Mid-P CI: 0.7267–14.4052; Fisher: 0.5997–17.4565

exactOrUI <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(class = "calc-layout",

      # ── Input panel ────────────────────────────────────────────────────────────
      tags$div(class = "panel",
        tags$div(class = "panel-head ph-teal", "Input"),
        tags$div(class = "panel-body",

          tags$p(class = "panel-hint",
            "Enter cell counts for a single 2\u00d72 table."),

          tags$table(class = "tbl",
            tags$thead(tags$tr(
              tags$th(class = "rh", ""),
              tags$th("Disease"),
              tags$th("No Disease"),
              tags$th("Total")
            )),
            tags$tbody(
              tags$tr(
                tags$td(class = "cl", "Exposed"),
                tags$td(numericInput(ns("ce"),  tags$span(class = "cell-lbl", "a"), value = NA, min = 0, width = "70px")),
                tags$td(numericInput(ns("nce"), tags$span(class = "cell-lbl", "b"), value = NA, min = 0, width = "70px")),
                tags$td(class = "tot", textOutput(ns("n1_tot"), inline = TRUE))
              ),
              tags$tr(
                tags$td(class = "cl", "Unexposed"),
                tags$td(numericInput(ns("cu"),  tags$span(class = "cell-lbl", "c"), value = NA, min = 0, width = "70px")),
                tags$td(numericInput(ns("ncu"), tags$span(class = "cell-lbl", "d"), value = NA, min = 0, width = "70px")),
                tags$td(class = "tot", textOutput(ns("n2_tot"), inline = TRUE))
              ),
              tags$tr(
                tags$td(class = "cl", "Total"),
                tags$td(class = "tot", textOutput(ns("m1_tot"), inline = TRUE)),
                tags$td(class = "tot", textOutput(ns("m0_tot"), inline = TRUE)),
                tags$td(class = "tot", textOutput(ns("N_tot"),  inline = TRUE))
              )
            )
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

exactOrServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    vals <- reactive({
      list(ce = input$ce, nce = input$nce, cu = input$cu, ncu = input$ncu)
    })

    # Marginal totals for table display
    output$n1_tot <- renderText({
      v <- vals()
      if (!is.na(v$ce) && !is.na(v$nce)) v$ce + v$nce else "\u2014"
    })
    output$n2_tot <- renderText({
      v <- vals()
      if (!is.na(v$cu) && !is.na(v$ncu)) v$cu + v$ncu else "\u2014"
    })
    output$m1_tot <- renderText({
      v <- vals()
      if (!is.na(v$ce) && !is.na(v$cu)) v$ce + v$cu else "\u2014"
    })
    output$m0_tot <- renderText({
      v <- vals()
      if (!is.na(v$nce) && !is.na(v$ncu)) v$nce + v$ncu else "\u2014"
    })
    output$N_tot <- renderText({
      v <- vals()
      if (!any(sapply(v, is.na))) Reduce("+", v) else "\u2014"
    })

    computed <- eventReactive(input$calc, {
      v <- vals()
      validate(
        need(!any(sapply(v, is.na)), "Enter all four cell counts."),
        need(all(sapply(v, function(x) x >= 0)), "All counts must be \u2265 0.")
      )

      x1 <- v$ce   # exposed cases
      x2 <- v$cu   # unexposed cases
      n1 <- v$ce  + v$nce   # exposed total
      n2 <- v$cu  + v$ncu   # unexposed total

      validate(
        need(n1 + n2 <= 20000, "Limit: n1 + n2 must be \u2264 20,000."),
        need(n1 > 0 && n2 > 0, "Row totals cannot be zero."),
        need(x1 + x2 > 0, "At least one case is required.")
      )

      # Point estimate
      OR_est <- if (x1 == 0 || (n1 - x1) == 0 || x2 == 0 || (n2 - x2) == 0) {
        NA_real_
      } else {
        (x1 * (n2 - x2)) / ((n1 - x1) * x2)
      }

      # 95% CI
      ci95_midp   <- limitcalc(x1, x2, n1, n2, level = 0.025, fishermid = FALSE)
      ci95_fisher <- limitcalc(x1, x2, n1, n2, level = 0.025, fishermid = TRUE)

      # 90% CI
      ci90_midp   <- limitcalc(x1, x2, n1, n2, level = 0.05,  fishermid = FALSE)
      ci90_fisher <- limitcalc(x1, x2, n1, n2, level = 0.05,  fishermid = TRUE)

      # P-values
      # Fisher exact (two-sided) via fisher.test
      ft <- fisher.test(matrix(c(x1, x2, n1 - x1, n2 - x2), nrow = 2))
      p_fisher <- ft$p.value

      # Mid-P: P(X >= x1) - 0.5*P(X = x1) under OR=1 (hypergeometric)
      # Under null OR=1, X ~ Hypergeometric(N, n1, x1+x2)
      N  <- n1 + n2
      m1 <- x1 + x2   # total cases
      # Two-sided mid-P
      x_min <- max(0L, m1 - n2)
      x_max <- min(m1, n1)
      probs_null <- dhyper(x_min:x_max, n1, n2, m1)
      p_obs <- probs_null[x1 - x_min + 1]
      # One-sided upper: P(X >= x1) - 0.5 * P(X = x1)
      p_upper_midp <- sum(probs_null[x_min:x_max >= x1]) - 0.5 * p_obs
      # One-sided lower: P(X <= x1) - 0.5 * P(X = x1)
      p_lower_midp <- sum(probs_null[x_min:x_max <= x1]) - 0.5 * p_obs
      p_midp <- 2 * min(p_upper_midp, p_lower_midp)

      list(
        OR       = OR_est,
        ci95_midp   = ci95_midp,
        ci95_fisher = ci95_fisher,
        ci90_midp   = ci90_midp,
        ci90_fisher = ci90_fisher,
        p_midp   = p_midp,
        p_fisher = p_fisher
      )
    })

    metric <- function(label, value, sub = NULL, cls = "m-blue") {
      tags$div(class = paste("metric", cls),
        tags$div(class = "metric-label", label),
        tags$div(class = "metric-value", value),
        if (!is.null(sub)) tags$div(class = "metric-sub", sub)
      )
    }

    ci_str <- function(lo, hi) {
      paste0(fmt4(lo), " \u2013 ", fmt4(hi))
    }

    output$results <- renderUI({
      r <- computed()
      tags$div(class = "metrics",
        metric("Odds Ratio", fmt4(r$OR),
               "Point estimate", "m-primary wide"),

        metric("95% Mid-P CI", ci_str(r$ci95_midp$lower, r$ci95_midp$upper),
               "Exact mid-p method", "m-blue"),
        metric("95% Fisher CI", ci_str(r$ci95_fisher$lower, r$ci95_fisher$upper),
               "Fisher exact method", "m-blue"),

        metric("90% Mid-P CI", ci_str(r$ci90_midp$lower, r$ci90_midp$upper),
               "Exact mid-p method", "m-teal"),
        metric("90% Fisher CI", ci_str(r$ci90_fisher$lower, r$ci90_fisher$upper),
               "Fisher exact method", "m-teal"),

        metric("Mid-P exact p", fmt4(r$p_midp),
               "Two-sided", "m-purple"),
        metric("Fisher exact p", fmt4(r$p_fisher),
               "Two-sided", "m-purple")
      )
    })

    output$interpretation <- renderUI({
      r <- computed()
      sig    <- !is.na(r$p_fisher) && r$p_fisher < 0.05
      dir    <- if (!is.na(r$OR) && r$OR > 1) "higher" else "lower"
      assoc  <- if (sig) "statistically significant (p < 0.05)"
                else "not statistically significant (p \u2265 0.05)"

      tags$div(class = "conclusion",
        tags$div(class = "conc-head", "Interpretation"),
        tags$p(class = "conc-text",
          "The estimated odds ratio is ", tags$strong(fmt4(r$OR)),
          " (95% Mid-P CI: ", ci_str(r$ci95_midp$lower, r$ci95_midp$upper),
          "; 95% Fisher CI: ",  ci_str(r$ci95_fisher$lower, r$ci95_fisher$upper), "). ",
          "The exposed group had ", tags$strong(dir), " odds of disease. ",
          "The Fisher exact two-sided p-value is ", tags$strong(fmt4(r$p_fisher)),
          "; the association is ", tags$strong(assoc), "."
        )
      )
    })

    outputOptions(output, "results",        suspendWhenHidden = FALSE)
    outputOptions(output, "interpretation", suspendWhenHidden = FALSE)
  })
}
