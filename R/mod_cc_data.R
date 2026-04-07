# mod_cc_data.R — Mantel-Haenszel Case-Control Data module
# Test values (CLAUDE.md): Stratum 1 Cases 104/1059, Controls 12/1163 → crude OR ≈ 7.333

ccDataUI <- function(id) {
  ns <- NS(id)

  # Helper: one stratum 2×2 block
  # Columns: Cases Exposed (a), Cases Unexposed (b), Controls Exposed (c), Controls Unexposed (d)
  stratum_inputs <- function(i) {
    lbl <- paste0("Stratum ", i)
    pfx <- paste0("s", i, "_")
    tags$div(class = "stratum-block",
      tags$div(class = "stratum-title", lbl),
      tags$div(class = "stratum-grid",
        tags$div(
          tags$label("Cases exposed (a)"),
          numericInput(ns(paste0(pfx, "a")), NULL, value = NULL, min = 0, step = 1)
        ),
        tags$div(
          tags$label("Cases unexposed (b)"),
          numericInput(ns(paste0(pfx, "b")), NULL, value = NULL, min = 0, step = 1)
        ),
        tags$div(
          tags$label("Controls exposed (c)"),
          numericInput(ns(paste0(pfx, "c")), NULL, value = NULL, min = 0, step = 1)
        ),
        tags$div(
          tags$label("Controls unexposed (d)"),
          numericInput(ns(paste0(pfx, "d")), NULL, value = NULL, min = 0, step = 1)
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
              tags$label("Cases exposed (a)"),
              numericInput(ns("cr_a"), NULL, value = NULL, min = 0, step = 1)
            ),
            tags$div(
              tags$label("Cases unexposed (b)"),
              numericInput(ns("cr_b"), NULL, value = NULL, min = 0, step = 1)
            ),
            tags$div(
              tags$label("Controls exposed (c)"),
              numericInput(ns("cr_c"), NULL, value = NULL, min = 0, step = 1)
            ),
            tags$div(
              tags$label("Controls unexposed (d)"),
              numericInput(ns("cr_d"), NULL, value = NULL, min = 0, step = 1)
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

ccDataServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    metric <- function(label, value, sub = NULL, cls = "m-blue") {
      tags$div(class = paste("metric", cls),
        tags$div(class = "metric-label", label),
        tags$div(class = "metric-value", value),
        if (!is.null(sub)) tags$div(class = "metric-sub", sub)
      )
    }

    # Collect complete strata — all four cells non-NA and >= 0, at least one cell > 0 per row
    get_strata <- reactive({
      strata <- list()
      for (i in 1:3) {
        pfx <- paste0("s", i, "_")
        a <- input[[paste0(pfx, "a")]]
        b <- input[[paste0(pfx, "b")]]
        c <- input[[paste0(pfx, "c")]]
        d <- input[[paste0(pfx, "d")]]
        if (!is.null(a) && !is.null(b) && !is.null(c) && !is.null(d) &&
            !is.na(a) && !is.na(b) && !is.na(c) && !is.na(d) &&
            a >= 0 && b >= 0 && c >= 0 && d >= 0 &&
            (a + b) > 0 && (c + d) > 0) {
          strata[[length(strata) + 1]] <- list(a = a, b = b, c = c, d = d)
        }
      }
      strata
    })

    computed <- eventReactive(input$calc, {
      strata <- get_strata()
      validate(need(length(strata) >= 1,
        "Enter data for at least one stratum (all four cells required; cases and controls must be > 0)."))

      # MH Odds Ratio:
      # OR_MH = sum(a_i * d_i / N_i) / sum(b_i * c_i / N_i)
      # N_i = a_i + b_i + c_i + d_i
      mh_num <- 0; mh_den <- 0

      # Woolf CI for OR_MH:
      # SE(ln OR_MH) = 1 / sqrt(sum(w_i))
      # where w_i = 1 / (1/a_i + 1/b_i + 1/c_i + 1/d_i)  [Woolf weights]
      # Note: w_i = a_i*b_i*c_i*d_i / (a_i*b_i + a_i*c_i + ... ) — use harmonic form
      sum_w_woolf <- 0

      # Stratum-level for homogeneity test
      ln_or_i <- numeric(length(strata))
      w_i     <- numeric(length(strata))

      for (i in seq_along(strata)) {
        s  <- strata[[i]]
        Ni <- s$a + s$b + s$c + s$d

        # MH components
        mh_num <- mh_num + s$a * s$d / Ni
        mh_den <- mh_den + s$b * s$c / Ni

        # Woolf weight for this stratum (requires all four cells > 0)
        if (s$a > 0 && s$b > 0 && s$c > 0 && s$d > 0) {
          or_i       <- (s$a * s$d) / (s$b * s$c)
          ln_or_i[i] <- log(or_i)
          # Woolf variance of ln(OR_i) = 1/a + 1/b + 1/c + 1/d
          var_ln_i   <- 1/s$a + 1/s$b + 1/s$c + 1/s$d
          w_i[i]     <- 1 / var_ln_i
          sum_w_woolf <- sum_w_woolf + w_i[i]
        } else {
          ln_or_i[i] <- NA_real_
          w_i[i]     <- NA_real_
        }
      }

      validate(
        need(mh_den > 0,
             "Denominator of MH-OR is zero \u2014 check that there are exposed controls in every stratum."),
        need(sum_w_woolf > 0,
             "Cannot compute Woolf CI: all cells must be > 0 in at least one stratum.")
      )

      OR_mh  <- mh_num / mh_den
      ln_or  <- log(OR_mh)
      se_woolf <- 1 / sqrt(sum_w_woolf)

      z95 <- qnorm(0.975)
      z99 <- qnorm(0.995)

      ci_woolf95 <- exp(ln_or + c(-1, 1) * z95 * se_woolf)
      ci_woolf99 <- exp(ln_or + c(-1, 1) * z99 * se_woolf)

      # Woolf homogeneity chi-sq:
      # Q = sum(w_i * (ln OR_i - ln OR_MH)^2)  with df = k - 1
      ok <- !is.na(w_i) & !is.na(ln_or_i) & w_i > 0
      chi2_hom <- NA_real_
      p_hom    <- NA_real_
      df_hom   <- sum(ok) - 1
      if (sum(ok) >= 2) {
        chi2_hom <- sum(w_i[ok] * (ln_or_i[ok] - ln_or)^2)
        p_hom    <- pchisq(chi2_hom, df = df_hom, lower.tail = FALSE)
      }

      # Crude OR
      crude_OR   <- NA_real_
      crude_ci95 <- c(NA_real_, NA_real_)
      cr_a <- input$cr_a; cr_b <- input$cr_b
      cr_c <- input$cr_c; cr_d <- input$cr_d
      if (!is.null(cr_a) && !is.null(cr_b) && !is.null(cr_c) && !is.null(cr_d) &&
          !is.na(cr_a) && !is.na(cr_b) && !is.na(cr_c) && !is.na(cr_d) &&
          cr_a > 0 && cr_b > 0 && cr_c > 0 && cr_d > 0) {
        crude_OR   <- (cr_a * cr_d) / (cr_b * cr_c)
        se_cr      <- sqrt(1/cr_a + 1/cr_b + 1/cr_c + 1/cr_d)
        crude_ci95 <- exp(log(crude_OR) + c(-1, 1) * z95 * se_cr)
      }

      list(
        OR_mh = OR_mh, se_woolf = se_woolf,
        ci_woolf95 = ci_woolf95, ci_woolf99 = ci_woolf99,
        chi2_hom = chi2_hom, p_hom = p_hom, df_hom = df_hom,
        crude_OR = crude_OR, crude_ci95 = crude_ci95,
        n_strata = length(strata)
      )
    })

    output$results <- renderUI({
      r <- computed()

      hom_text <- if (!is.na(r$chi2_hom)) {
        paste0("\u03c7\u00b2(", r$df_hom, ") = ", fmt4(r$chi2_hom),
               ",  p = ", fmt4(r$p_hom))
      } else {
        "Need \u2265 2 complete strata with all cells > 0"
      }

      crude_text <- if (!is.na(r$crude_OR)) {
        paste0(fmt4(r$crude_OR), "  (95% CI: ",
               fmt4(r$crude_ci95[1]), "\u2013", fmt4(r$crude_ci95[2]), ")")
      } else {
        "\u2014"
      }

      tags$div(class = "metrics",
        metric("MH Odds Ratio",
               fmt4(r$OR_mh),
               paste0(r$n_strata, " strat",
                      if (r$n_strata == 1) "um" else "a", " pooled"),
               "m-primary wide"),
        metric("95% CI (Woolf)",
               paste0(fmt4(r$ci_woolf95[1]), " \u2013 ", fmt4(r$ci_woolf95[2])),
               cls = "m-blue"),
        metric("99% CI (Woolf)",
               paste0(fmt4(r$ci_woolf99[1]), " \u2013 ", fmt4(r$ci_woolf99[2])),
               cls = "m-blue"),
        metric("Crude Odds Ratio",
               crude_text,
               cls = "m-teal"),
        metric("P for Homogeneity",
               hom_text,
               cls = "m-purple")
      )
    })

    output$interpretation <- renderUI({
      r <- computed()
      sig_hom <- !is.na(r$p_hom) && r$p_hom < 0.05
      dir     <- if (r$OR_mh > 1) "higher odds" else "lower odds"
      conf_str <- paste0("95% Woolf CI: ",
                         fmt4(r$ci_woolf95[1]), "\u2013", fmt4(r$ci_woolf95[2]))

      tags$div(class = "conclusion",
        tags$div(class = "conc-head", "Interpretation"),
        tags$p(class = "conc-text",
          "The Mantel-Haenszel odds ratio pooled across ", r$n_strata,
          " strat", if (r$n_strata == 1) "um" else "a", " is ",
          tags$strong(fmt4(r$OR_mh)),
          " (", conf_str, "; 99% CI: ",
          fmt4(r$ci_woolf99[1]), "\u2013", fmt4(r$ci_woolf99[2]), "). ",
          "Cases had ", tags$strong(dir),
          " of exposure than controls. ",
          if (!is.na(r$crude_OR)) {
            paste0("The crude OR was ", fmt4(r$crude_OR), ". ")
          },
          if (!is.na(r$p_hom)) {
            paste0(
              "The Woolf test for homogeneity across strata gave \u03c7\u00b2 = ",
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
