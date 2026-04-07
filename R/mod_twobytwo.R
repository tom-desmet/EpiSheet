# mod_twobytwo.R — 2×2 Table module

twobytwoUI <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(class = "calc-layout",

      # ── Input panel ──────────────────────────────────────────────────────────
      tags$div(class = "panel",
        tags$div(class = "panel-head ph-teal", "Input"),
        tags$div(class = "panel-body",
          tags$table(class = "tbl",
            tags$thead(tags$tr(
              tags$th(class = "rh", ""),
              tags$th("Cases"),
              tags$th("Non-cases"),
              tags$th("Total")
            )),
            tags$tbody(
              tags$tr(
                tags$td(class = "cl", "Exposed"),
                tags$td(tags$div(class = "cell-wrap",
                  tags$span(class = "cell-letter", "a"),
                  tags$input(id = ns("a"), type = "number", min = "0", placeholder = "0")
                )),
                tags$td(tags$div(class = "cell-wrap",
                  tags$span(class = "cell-letter", "b"),
                  tags$input(id = ns("b"), type = "number", min = "0", placeholder = "0")
                )),
                tags$td(class = "tot", textOutput(ns("n1"), inline = TRUE))
              ),
              tags$tr(
                tags$td(class = "cl", "Unexposed"),
                tags$td(tags$div(class = "cell-wrap",
                  tags$span(class = "cell-letter", "c"),
                  tags$input(id = ns("c"), type = "number", min = "0", placeholder = "0")
                )),
                tags$td(tags$div(class = "cell-wrap",
                  tags$span(class = "cell-letter", "d"),
                  tags$input(id = ns("d"), type = "number", min = "0", placeholder = "0")
                )),
                tags$td(class = "tot", textOutput(ns("n0"), inline = TRUE))
              ),
              tags$tr(
                tags$td(class = "cl", "Total"),
                tags$td(class = "tot", textOutput(ns("m1"), inline = TRUE)),
                tags$td(class = "tot", textOutput(ns("m0"), inline = TRUE)),
                tags$td(class = "tot", textOutput(ns("N"),  inline = TRUE))
              )
            )
          ),
          # Wire raw inputs from HTML input elements via JS
          tags$script(HTML(sprintf("
            ['%s','%s','%s','%s'].forEach(function(id) {
              document.getElementById(id).addEventListener('input', function() {
                Shiny.setInputValue(id, this.value === '' ? null : +this.value, {priority: 'event'});
              });
            });
          ", ns("a"), ns("b"), ns("c"), ns("d"))),
          tags$button(class = "calc-btn",
            onclick = sprintf("Shiny.setInputValue('%s', Math.random())", ns("calc")),
            "Calculate")
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

    # ── Interpretation block ────────────────────────────────────────────────────
    uiOutput(ns("interpretation"))
  )
}

twobytwoServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Reactive cell values (from JS-wired inputs)
    vals <- reactive({
      list(
        a = input$a, b = input$b,
        c = input$c, d = input$d
      )
    })

    # Marginal totals
    output$n1 <- renderText({ v <- vals(); if (!is.null(v$a) && !is.null(v$b)) v$a + v$b else "—" })
    output$n0 <- renderText({ v <- vals(); if (!is.null(v$c) && !is.null(v$d)) v$c + v$d else "—" })
    output$m1 <- renderText({ v <- vals(); if (!is.null(v$a) && !is.null(v$c)) v$a + v$c else "—" })
    output$m0 <- renderText({ v <- vals(); if (!is.null(v$b) && !is.null(v$d)) v$b + v$d else "—" })
    output$N  <- renderText({ v <- vals(); if (!any(sapply(v, is.null))) Reduce("+", v) else "—" })

    # Core computation, triggered by calculate button
    computed <- eventReactive(input$calc, {
      v <- vals()
      a <- v$a; b <- v$b; c <- v$c; d <- v$d
      validate(
        need(!any(sapply(v, is.null)), "Enter all four cell counts."),
        need(all(sapply(v, function(x) x >= 0)), "All counts must be \u2265 0.")
      )

      n1 <- a + b; n0 <- c + d
      m1 <- a + c; m0 <- b + d
      N  <- n1 + n0

      validate(need(n1 > 0 && n0 > 0, "Row totals cannot be zero."))
      validate(need(b > 0 && c > 0, "Cannot compute ratios: one or more cells are zero."))

      p1 <- a / n1; p0 <- c / n0
      RR <- p1 / p0
      RD <- p1 - p0
      OR <- (a * d) / (b * c)

      # CI (Taylor)
      se_ln_rr <- sqrt(b/(a*n1) + d/(c*n0))
      rr95 <- exp(log(RR) + c(-1,1) * 1.96  * se_ln_rr)
      rr99 <- exp(log(RR) + c(-1,1) * 2.576 * se_ln_rr)

      se_rd <- sqrt(p1*(1-p1)/n1 + p0*(1-p0)/n0)
      rd95 <- RD + c(-1,1) * 1.96  * se_rd
      rd99 <- RD + c(-1,1) * 2.576 * se_rd

      se_ln_or <- sqrt(1/a + 1/b + 1/c + 1/d)
      or95 <- exp(log(OR) + c(-1,1) * 1.96  * se_ln_or)
      or99 <- exp(log(OR) + c(-1,1) * 2.576 * se_ln_or)

      # Chi-square
      chi2    <- N * (a*d - b*c)^2 / (n1 * n0 * m1 * m0)
      chi2_cc <- N * (max(0, abs(a*d - b*c) - N/2))^2 / (n1 * n0 * m1 * m0)
      p_chi2    <- pchisq(chi2,    1, lower.tail = FALSE)
      p_chi2_cc <- pchisq(chi2_cc, 1, lower.tail = FALSE)

      p_fisher <- fisher.test(matrix(c(a, c, b, d), 2, 2))$p.value

      # Attributable risk
      ar_exp <- (p1 - p0) / p1
      ar_pop <- (p1 * n1/N - p0 * n0/N) / (p1 * n1/N + p0 * n0/N)  # simplified population AR
      nnt    <- 1 / abs(RD)

      list(RR=RR, rr95=rr95, rr99=rr99,
           RD=RD, rd95=rd95, rd99=rd99,
           OR=OR, or95=or95, or99=or99,
           chi2=chi2, chi2_cc=chi2_cc, p_chi2=p_chi2, p_chi2_cc=p_chi2_cc,
           p_fisher=p_fisher, ar_exp=ar_exp, ar_pop=ar_pop, nnt=nnt,
           p1=p1, p0=p0)
    })

    metric <- function(label, value, sub = NULL, cls = "m-blue") {
      tags$div(class = paste("metric", cls),
        tags$div(class = "metric-label", label),
        tags$div(class = "metric-value", value),
        if (!is.null(sub)) tags$div(class = "metric-sub", sub)
      )
    }

    output$results <- renderUI({
      r <- computed()
      tags$div(class = "metrics",
        # Primary: RR
        metric("Risk Ratio (RR)", fmt4(r$RR),
               paste0("95% CI: ", fmt4(r$rr95[1]), "\u2013", fmt4(r$rr95[2])),
               "m-primary wide"),
        # Row 2
        metric("Odds Ratio (OR)", fmt4(r$OR),
               paste0("95% CI: ", fmt4(r$or95[1]), "\u2013", fmt4(r$or95[2])),
               "m-blue"),
        metric("Risk Difference", fmt4(r$RD),
               paste0("95% CI: ", fmt4(r$rd95[1]), "\u2013", fmt4(r$rd95[2])),
               "m-blue"),
        # Row 3
        metric("\u03c7\u00b2 (no cc)", paste0(fmt4(r$chi2), "  p=", fmt4(r$p_chi2)),   cls = "m-teal"),
        metric("\u03c7\u00b2 (cc)",     paste0(fmt4(r$chi2_cc), "  p=", fmt4(r$p_chi2_cc)), cls = "m-teal"),
        # Row 4
        metric("Fisher exact p", fmt4(r$p_fisher), cls = "m-purple"),
        metric("AR (exposed)",   fmt4(r$ar_exp),   cls = "m-purple"),
        # Row 5
        metric("AR (population)", fmt4(r$ar_pop), cls = "m-blue"),
        metric(if (r$RD > 0) "NNH" else "NNT", fmt2(r$nnt), cls = "m-blue")
      )
    })

    output$interpretation <- renderUI({
      r <- computed()
      dir  <- if (r$RR > 1) "higher" else "lower"
      pct  <- abs(round((r$RR - 1) * 100))
      sig  <- r$p_chi2 < 0.05
      assoc <- if (sig) "statistically significant (p < 0.05)" else "not statistically significant (p \u2265 0.05)"

      tags$div(class = "conclusion",
        tags$div(class = "conc-head", "Interpretation"),
        tags$p(class = "conc-text",
          "The risk of disease was ", tags$strong(paste0(pct, "% ", dir)),
          " in the exposed group (RR = ", tags$strong(fmt4(r$RR)), ", 95% CI: ",
          fmt4(r$rr95[1]), "\u2013", fmt4(r$rr95[2]), "). ",
          "The odds ratio was ", tags$strong(fmt4(r$OR)),
          " and the absolute risk difference was ", tags$strong(fmt4(r$RD)), ". ",
          "The association is ", tags$strong(assoc), "."
        )
      )
    })

    outputOptions(output, "results", suspendWhenHidden = FALSE)
    outputOptions(output, "interpretation", suspendWhenHidden = FALSE)
  })
}
