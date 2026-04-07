# mod_risk_data.R — Mantel-Haenszel Risk Data module

riskDataUI <- function(id) {
  ns <- NS(id)

  # Helper: one stratum block
  stratum_inputs <- function(i) {
    lbl <- paste0("Stratum ", i)
    pfx <- paste0("s", i, "_")
    tags$div(class = "stratum-block",
      tags$div(class = "stratum-title", lbl),
      tags$div(class = "stratum-grid",
        tags$div(
          tags$label("Cases (exposed)"),
          numericInput(ns(paste0(pfx, "a")), NULL, value = NULL, min = 0, step = 1)
        ),
        tags$div(
          tags$label("Total exposed (N\u2081)"),
          numericInput(ns(paste0(pfx, "N1")), NULL, value = NULL, min = 1, step = 1)
        ),
        tags$div(
          tags$label("Cases (unexposed)"),
          numericInput(ns(paste0(pfx, "c")), NULL, value = NULL, min = 0, step = 1)
        ),
        tags$div(
          tags$label("Total unexposed (N\u2080)"),
          numericInput(ns(paste0(pfx, "N0")), NULL, value = NULL, min = 1, step = 1)
        )
      )
    )
  }

  tagList(
    tags$div(class = "calc-layout",

      # ── Input panel ────────────────────────────────────────────────────────────
      tags$div(class = "panel",
        tags$div(class = "panel-head ph-teal", "Input"),
        tags$div(class = "panel-body",

          stratum_inputs(1),
          stratum_inputs(2),
          stratum_inputs(3),

          tags$hr(),
          tags$div(class = "stratum-title", "Crude Data"),
          tags$div(class = "stratum-grid",
            tags$div(
              tags$label("Cases (exposed)"),
              numericInput(ns("cr_a"), NULL, value = NULL, min = 0, step = 1)
            ),
            tags$div(
              tags$label("Total exposed"),
              numericInput(ns("cr_N1"), NULL, value = NULL, min = 1, step = 1)
            ),
            tags$div(
              tags$label("Cases (unexposed)"),
              numericInput(ns("cr_c"), NULL, value = NULL, min = 0, step = 1)
            ),
            tags$div(
              tags$label("Total unexposed"),
              numericInput(ns("cr_N0"), NULL, value = NULL, min = 1, step = 1)
            )
          ),

          tags$br(),
          tags$button(class = "calc-btn",
            onclick = sprintf("Shiny.setInputValue('%s', Math.random())", ns("calc")),
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

riskDataServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    metric <- function(label, value, sub = NULL, cls = "m-blue") {
      tags$div(class = paste("metric", cls),
        tags$div(class = "metric-label", label),
        tags$div(class = "metric-value", value),
        if (!is.null(sub)) tags$div(class = "metric-sub", sub)
      )
    }

    # Collect complete strata
    get_strata <- reactive({
      strata <- list()
      for (i in 1:3) {
        pfx <- paste0("s", i, "_")
        a  <- input[[paste0(pfx, "a")]]
        N1 <- input[[paste0(pfx, "N1")]]
        c  <- input[[paste0(pfx, "c")]]
        N0 <- input[[paste0(pfx, "N0")]]
        if (!is.null(a) && !is.null(N1) && !is.null(c) && !is.null(N0) &&
            !is.na(a) && !is.na(N1) && !is.na(c) && !is.na(N0) &&
            N1 > 0 && N0 > 0 && a <= N1 && c <= N0) {
          strata[[length(strata) + 1]] <- list(a = a, N1 = N1, c = c, N0 = N0)
        }
      }
      strata
    })

    computed <- eventReactive(input$calc, {
      strata <- get_strata()
      validate(need(length(strata) >= 1,
        "Enter data for at least one stratum (all four fields required, cases \u2264 total)."))

      # MH Risk Ratio:
      # RR_MH = sum(a_i * N0_i / N_i) / sum(c_i * N1_i / N_i)
      # N_i = N1_i + N0_i
      R_num <- 0; R_den <- 0

      # Pooled Risk Difference:
      # RD_pool = sum(w_i * RD_i) / sum(w_i)
      # w_i = N1_i * N0_i / N_i  (MH weights)
      # Var(RD_pool) = sum(w_i^2 * Var_i) / sum(w_i)^2
      # Var_i = p1_i*(1-p1_i)/N1_i + p0_i*(1-p0_i)/N0_i
      sum_w   <- 0; sum_wRD <- 0; sum_w2var <- 0

      # Taylor variance for ln(RR_MH):
      # Var(ln RR_MH) = sum(P_i * R_i + Q_i * S_i) / (2 * R_MH_num * R_MH_den)
      # where P_i = a_i/N1_i, Q_i = 1-P_i, R_i = a_i*N0_i/N_i, S_i = c_i*N1_i/N_i
      # Using Greenland-Robins formula:
      # = [sum((N1_i*N0_i*(a_i+c_i) - a_i*c_i*N_i)/N_i^2)] / (R_num * R_den)
      var_num <- 0

      # Stratum-specific estimates for homogeneity (Woolf)
      ln_rr_i <- numeric(length(strata))
      w_i     <- numeric(length(strata))

      for (i in seq_along(strata)) {
        s  <- strata[[i]]
        Ni <- s$N1 + s$N0
        p1 <- s$a / s$N1
        p0 <- s$c / s$N0
        RD_i <- p1 - p0

        # MH components
        ri <- s$a * s$N0 / Ni   # contributes to numerator
        si <- s$c * s$N1 / Ni   # contributes to denominator
        R_num <- R_num + ri
        R_den <- R_den + si

        # Greenland-Robins Taylor variance term for ln(RR)
        var_num <- var_num + (s$N1 * s$N0 * (s$a + s$c) - s$a * s$c * Ni) / Ni^2

        # Pooled RD accumulation
        wi       <- s$N1 * s$N0 / Ni
        var_RD_i <- p1 * (1 - p1) / s$N1 + p0 * (1 - p0) / s$N0
        sum_w    <- sum_w  + wi
        sum_wRD  <- sum_wRD + wi * RD_i
        sum_w2var <- sum_w2var + wi^2 * var_RD_i

        # Woolf weight for homogeneity
        b  <- s$N1 - s$a    # non-cases exposed
        d  <- s$N0 - s$c    # non-cases unexposed
        if (s$a > 0 && b > 0 && s$c > 0 && d > 0) {
          rr_i       <- p1 / p0
          ln_rr_i[i] <- log(rr_i)
          # Var(ln RR_i) approx = 1/a - 1/N1 + 1/c - 1/N0
          w_i[i]     <- 1 / (1/s$a - 1/s$N1 + 1/s$c - 1/s$N0)
        } else {
          ln_rr_i[i] <- NA_real_
          w_i[i]     <- NA_real_
        }
      }

      validate(need(R_den > 0,
        "Denominator of MH estimator is zero \u2014 check that unexposed group has cases."))

      RR_mh  <- R_num / R_den
      var_ln <- if (R_num > 0 && R_den > 0) var_num / (R_num * R_den) else NA_real_

      z95 <- qnorm(0.975)
      z99 <- qnorm(0.995)

      ci_rr95 <- exp(log(RR_mh) + c(-1, 1) * z95 * sqrt(var_ln))
      ci_rr99 <- exp(log(RR_mh) + c(-1, 1) * z99 * sqrt(var_ln))

      # Pooled RD and its CI
      RD_pool <- sum_wRD / sum_w
      se_RD   <- sqrt(sum_w2var) / sum_w
      ci_rd95 <- RD_pool + c(-1, 1) * z95 * se_RD
      ci_rd99 <- RD_pool + c(-1, 1) * z99 * se_RD

      # Woolf homogeneity chi-sq
      ok <- !is.na(w_i) & !is.na(ln_rr_i) & w_i > 0
      chi2_hom <- NA_real_
      p_hom    <- NA_real_
      df_hom   <- sum(ok) - 1
      if (sum(ok) >= 2) {
        ln_rr_w  <- sum(w_i[ok] * ln_rr_i[ok]) / sum(w_i[ok])
        chi2_hom <- sum(w_i[ok] * (ln_rr_i[ok] - ln_rr_w)^2)
        p_hom    <- pchisq(chi2_hom, df = df_hom, lower.tail = FALSE)
      }

      # Crude RR
      crude_RR   <- NA_real_
      crude_ci95 <- c(NA_real_, NA_real_)
      cr_a  <- input$cr_a
      cr_N1 <- input$cr_N1
      cr_c  <- input$cr_c
      cr_N0 <- input$cr_N0
      if (!is.null(cr_a) && !is.null(cr_N1) && !is.null(cr_c) && !is.null(cr_N0) &&
          !is.na(cr_a) && !is.na(cr_N1) && !is.na(cr_c) && !is.na(cr_N0) &&
          cr_N1 > 0 && cr_N0 > 0 && cr_a > 0 && cr_c > 0 &&
          cr_a <= cr_N1 && cr_c <= cr_N0) {
        p1_cr      <- cr_a / cr_N1
        p0_cr      <- cr_c / cr_N0
        crude_RR   <- p1_cr / p0_cr
        b_cr       <- cr_N1 - cr_a
        d_cr       <- cr_N0 - cr_c
        se_cr      <- sqrt(b_cr / (cr_a * cr_N1) + d_cr / (cr_c * cr_N0))
        crude_ci95 <- exp(log(crude_RR) + c(-1, 1) * z95 * se_cr)
      }

      list(
        RR_mh = RR_mh, var_ln = var_ln,
        ci_rr95 = ci_rr95, ci_rr99 = ci_rr99,
        RD_pool = RD_pool,
        ci_rd95 = ci_rd95, ci_rd99 = ci_rd99,
        chi2_hom = chi2_hom, p_hom = p_hom, df_hom = df_hom,
        crude_RR = crude_RR, crude_ci95 = crude_ci95,
        n_strata = length(strata)
      )
    })

    output$results <- renderUI({
      r <- computed()

      hom_text <- if (!is.na(r$chi2_hom)) {
        paste0("\u03c7\u00b2(", r$df_hom, ") = ", fmt4(r$chi2_hom),
               ",  p = ", fmt4(r$p_hom))
      } else {
        "Need \u2265 2 complete strata"
      }

      crude_text <- if (!is.na(r$crude_RR)) {
        paste0(fmt4(r$crude_RR), "  (95% CI: ",
               fmt4(r$crude_ci95[1]), "\u2013", fmt4(r$crude_ci95[2]), ")")
      } else {
        "\u2014"
      }

      tags$div(class = "metrics",
        metric("MH Risk Ratio",
               fmt4(r$RR_mh),
               paste0(r$n_strata, " strat",
                      if (r$n_strata == 1) "um" else "a", " pooled"),
               "m-primary wide"),
        metric("95% CI for Risk Ratio",
               paste0(fmt4(r$ci_rr95[1]), " \u2013 ", fmt4(r$ci_rr95[2])),
               cls = "m-blue"),
        metric("99% CI for Risk Ratio",
               paste0(fmt4(r$ci_rr99[1]), " \u2013 ", fmt4(r$ci_rr99[2])),
               cls = "m-blue"),
        metric("Pooled Risk Difference",
               fmt4(r$RD_pool),
               paste0("95% CI: ", fmt4(r$ci_rd95[1]), "\u2013", fmt4(r$ci_rd95[2])),
               cls = "m-teal"),
        metric("99% CI for Risk Difference",
               paste0(fmt4(r$ci_rd99[1]), " \u2013 ", fmt4(r$ci_rd99[2])),
               cls = "m-teal"),
        metric("P for Homogeneity",
               hom_text,
               cls = "m-purple"),
        metric("Crude Risk Ratio",
               crude_text,
               cls = "m-blue")
      )
    })

    output$interpretation <- renderUI({
      r <- computed()
      sig_hom <- !is.na(r$p_hom) && r$p_hom < 0.05
      dir_rr  <- if (r$RR_mh > 1) "higher" else "lower"
      dir_rd  <- if (r$RD_pool > 0) "excess" else "deficit"

      tags$div(class = "conclusion",
        tags$div(class = "conc-head", "Interpretation"),
        tags$p(class = "conc-text",
          "The Mantel-Haenszel risk ratio pooled across ", r$n_strata,
          " strat", if (r$n_strata == 1) "um" else "a", " is ",
          tags$strong(fmt4(r$RR_mh)),
          " (95% CI: ", fmt4(r$ci_rr95[1]), "\u2013", fmt4(r$ci_rr95[2]), "). ",
          "The exposed group had a ", tags$strong(paste0(dir_rr, " risk")),
          " of the outcome. ",
          "The pooled risk difference is ",
          tags$strong(fmt4(r$RD_pool)),
          " (95% CI: ", fmt4(r$ci_rd95[1]), "\u2013", fmt4(r$ci_rd95[2]), "), ",
          "representing an absolute ", dir_rd, " in risk. ",
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
