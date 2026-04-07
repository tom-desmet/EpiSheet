# mod_rate_data.R — Mantel-Haenszel Rate Data module

rateDataUI <- function(id) {
  ns <- NS(id)

  # Helper: one stratum block
  # Columns: Cases / Person-time; Rows: Exposed / Unexposed
  stratum_inputs <- function(i) {
    lbl <- paste0("Stratum ", i)
    pfx <- paste0("s", i, "_")
    tags$div(style = "margin-bottom: 12px;",
      tags$div(class = "strata-label", lbl),
      tags$table(class = "tbl",
        tags$thead(tags$tr(
          tags$th(class = "rh", ""),
          tags$th("Cases"),
          tags$th("Person-time")
        )),
        tags$tbody(
          tags$tr(
            tags$td(class = "cl", "Exposed"),
            tags$td(numericInput(ns(paste0(pfx, "a")),  tags$span(class = "cell-lbl", "a"),  value = NA, min = 0, width = "70px")),
            tags$td(numericInput(ns(paste0(pfx, "T1")), tags$span(class = "cell-lbl", "T\u2081"), value = NA, min = 0, width = "70px"))
          ),
          tags$tr(
            tags$td(class = "cl", "Unexposed"),
            tags$td(numericInput(ns(paste0(pfx, "b")),  tags$span(class = "cell-lbl", "b"),  value = NA, min = 0, width = "70px")),
            tags$td(numericInput(ns(paste0(pfx, "T0")), tags$span(class = "cell-lbl", "T\u2080"), value = NA, min = 0, width = "70px"))
          )
        )
      )
    )
  }

  tagList(
    tags$div(class = "calc-layout",

      # ── Input panel ──────────────────────────────────────────────────────────
      tags$div(class = "panel",
        tags$div(class = "panel-head ph-teal", "Input"),
        tags$div(class = "panel-body",

          stratum_inputs(1),
          stratum_inputs(2),
          stratum_inputs(3),

          tags$hr(),
          tags$div(style = "margin-bottom: 12px;",
            tags$div(class = "strata-label", "Crude Data"),
            tags$table(class = "tbl",
              tags$thead(tags$tr(
                tags$th(class = "rh", ""),
                tags$th("Cases"),
                tags$th("Person-time")
              )),
              tags$tbody(
                tags$tr(
                  tags$td(class = "cl", "Exposed"),
                  tags$td(numericInput(ns("cr_a"),  tags$span(class = "cell-lbl", "a"),  value = NA, min = 0, width = "70px")),
                  tags$td(numericInput(ns("cr_T1"), tags$span(class = "cell-lbl", "T\u2081"), value = NA, min = 0, width = "70px"))
                ),
                tags$tr(
                  tags$td(class = "cl", "Unexposed"),
                  tags$td(numericInput(ns("cr_b"),  tags$span(class = "cell-lbl", "b"),  value = NA, min = 0, width = "70px")),
                  tags$td(numericInput(ns("cr_T0"), tags$span(class = "cell-lbl", "T\u2080"), value = NA, min = 0, width = "70px"))
                )
              )
            )
          ),

          tags$br(),
          tags$button(class = "calc-btn",
            onclick = sprintf("Shiny.setInputValue('%s', Math.random())", ns("calc")),
            "Calculate")
        )
      ),

      # ── Results panel ────────────────────────────────────────────────────────
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

rateDataServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    metric <- function(label, value, sub = NULL, cls = "m-blue") {
      tags$div(class = paste("metric", cls),
        tags$div(class = "metric-label", label),
        tags$div(class = "metric-value", value),
        if (!is.null(sub)) tags$div(class = "metric-sub", sub)
      )
    }

    # Collect stratum data — returns list of complete strata (all 4 fields non-NA and valid)
    get_strata <- reactive({
      strata <- list()
      for (i in 1:3) {
        pfx <- paste0("s", i, "_")
        a  <- input[[paste0(pfx, "a")]]
        T1 <- input[[paste0(pfx, "T1")]]
        b  <- input[[paste0(pfx, "b")]]
        T0 <- input[[paste0(pfx, "T0")]]
        if (!is.na(a) && !is.na(T1) && !is.na(b) && !is.na(T0) &&
            T1 > 0 && T0 > 0) {
          strata[[length(strata) + 1]] <- list(a = a, T1 = T1, b = b, T0 = T0)
        }
      }
      strata
    })

    computed <- eventReactive(input$calc, {
      strata <- get_strata()
      validate(need(length(strata) >= 1, "Enter data for at least one stratum (all four fields required)."))

      # MH Rate Ratio components
      # R_MH = sum(a_i * T0_i / T_i) / sum(b_i * T1_i / T_i)
      R_num <- 0; R_den <- 0
      # Taylor variance numerator accumulator
      var_num <- 0
      # Homogeneity: stratum-specific RR for Woolf chi-sq
      ln_rr_i <- numeric(length(strata))
      w_i     <- numeric(length(strata))

      for (i in seq_along(strata)) {
        s  <- strata[[i]]
        Ti <- s$T1 + s$T0
        R_num <- R_num + s$a * s$T0 / Ti
        R_den <- R_den + s$b * s$T1 / Ti
        # Taylor variance per stratum:
        # numerator term = (T1*T0*(a+b) - a*b*T) / T^2
        var_num <- var_num + (s$T1 * s$T0 * (s$a + s$b) - s$a * s$b * Ti) / Ti^2
        # Stratum RR and Woolf weight
        if (s$b > 0 && s$T1 > 0 && s$a > 0 && s$T0 > 0) {
          rr_i       <- (s$a / s$T1) / (s$b / s$T0)
          ln_rr_i[i] <- log(rr_i)
          w_i[i]     <- 1 / (1/s$a + 1/s$b)
        } else {
          ln_rr_i[i] <- NA_real_
          w_i[i]     <- NA_real_
        }
      }

      validate(need(R_den > 0, "Denominator of MH estimator is zero — check unexposed person-time and cases."))

      RR_mh  <- R_num / R_den
      var_ln <- var_num / (R_num * R_den)

      # Confidence intervals
      z95  <- qnorm(0.975); z99 <- qnorm(0.995)
      ci95 <- exp(log(RR_mh) + c(-1, 1) * z95 * sqrt(var_ln))
      ci99 <- exp(log(RR_mh) + c(-1, 1) * z99 * sqrt(var_ln))

      # Woolf chi-sq for homogeneity
      ok <- !is.na(w_i) & !is.na(ln_rr_i) & w_i > 0
      chi2_hom <- NA_real_
      p_hom    <- NA_real_
      df_hom   <- sum(ok) - 1
      if (sum(ok) >= 2) {
        ln_rr_w   <- sum(w_i[ok] * ln_rr_i[ok]) / sum(w_i[ok])
        chi2_hom  <- sum(w_i[ok] * (ln_rr_i[ok] - ln_rr_w)^2)
        p_hom     <- pchisq(chi2_hom, df = df_hom, lower.tail = FALSE)
      }

      # Crude rate ratio
      crude_RR  <- NA_real_
      crude_ci95 <- c(NA_real_, NA_real_)
      cr_a  <- input$cr_a
      cr_T1 <- input$cr_T1
      cr_b  <- input$cr_b
      cr_T0 <- input$cr_T0
      if (!is.na(cr_a) && !is.na(cr_T1) && !is.na(cr_b) && !is.na(cr_T0) &&
          cr_T1 > 0 && cr_T0 > 0 && cr_b > 0 && cr_a > 0) {
        crude_RR   <- (cr_a / cr_T1) / (cr_b / cr_T0)
        se_cr      <- sqrt(1/cr_a + 1/cr_b)
        crude_ci95 <- exp(log(crude_RR) + c(-1, 1) * z95 * se_cr)
      }

      list(
        RR_mh = RR_mh, var_ln = var_ln,
        ci95 = ci95, ci99 = ci99,
        chi2_hom = chi2_hom, p_hom = p_hom, df_hom = df_hom,
        crude_RR = crude_RR, crude_ci95 = crude_ci95,
        n_strata = length(strata)
      )
    })

    output$results <- renderUI({
      r <- computed()

      hom_text <- if (!is.na(r$chi2_hom)) {
        paste0("\u03c7\u00b2(", r$df_hom, ") = ", fmt4(r$chi2_hom), ",  p = ", fmt4(r$p_hom))
      } else {
        "Need \u2265 2 complete strata"
      }

      crude_text <- if (!is.na(r$crude_RR)) {
        paste0(fmt4(r$crude_RR), "  (95% CI: ", fmt4(r$crude_ci95[1]), "\u2013", fmt4(r$crude_ci95[2]), ")")
      } else {
        "\u2014"
      }

      tags$div(class = "metrics",
        metric("MH Rate Ratio",
               fmt4(r$RR_mh),
               paste0(r$n_strata, " strat", if (r$n_strata == 1) "um" else "a", " pooled"),
               "m-primary wide"),
        metric("95% Confidence Interval",
               paste0(fmt4(r$ci95[1]), " \u2013 ", fmt4(r$ci95[2])),
               cls = "m-blue"),
        metric("99% Confidence Interval",
               paste0(fmt4(r$ci99[1]), " \u2013 ", fmt4(r$ci99[2])),
               cls = "m-blue"),
        metric("Variance of ln(RR)",
               fmt4(r$var_ln),
               cls = "m-teal"),
        metric("P for Homogeneity",
               hom_text,
               cls = "m-purple"),
        metric("Crude Rate Ratio",
               crude_text,
               cls = "m-teal")
      )
    })

    output$interpretation <- renderUI({
      r <- computed()
      sig_hom <- !is.na(r$p_hom) && r$p_hom < 0.05
      dir     <- if (r$RR_mh > 1) "higher" else "lower"

      tags$div(class = "conclusion",
        tags$div(class = "conc-head", "Interpretation"),
        tags$p(class = "conc-text",
          "The Mantel-Haenszel rate ratio pooled across ", r$n_strata,
          " strat", if (r$n_strata == 1) "um" else "a", " is ",
          tags$strong(fmt4(r$RR_mh)),
          " (95% CI: ", fmt4(r$ci95[1]), "\u2013", fmt4(r$ci95[2]),
          "; 99% CI: ", fmt4(r$ci99[1]), "\u2013", fmt4(r$ci99[2]), "). ",
          "The exposed group had a ", tags$strong(paste0(dir, " rate")),
          " than the unexposed. ",
          if (!is.na(r$p_hom)) {
            paste0(
              "The test for homogeneity across strata gave \u03c7\u00b2 = ",
              fmt4(r$chi2_hom), " (p = ", fmt4(r$p_hom), "); ",
              if (sig_hom)
                "there is significant heterogeneity between strata (p < 0.05), so the pooled estimate should be interpreted with caution."
              else
                "there is no significant heterogeneity between strata (p \u2265 0.05)."
            )
          }
        )
      )
    })

    outputOptions(output, "results",        suspendWhenHidden = FALSE)
    outputOptions(output, "interpretation", suspendWhenHidden = FALSE)
  })
}
