# mod_pvalue_functions.R — P-value Functions module
#
# For each of two curves the user supplies CI bounds at three levels (90/95/99%).
# The module recovers a point estimate and SE from the 95% bounds, then plots
# the p-value function:
#
#   p(RR_0) = 2 * (1 - pnorm( |ln(est) - ln(RR_0)| / SE ))
#
# Reference: Poole C. Beyond the confidence interval.
#            Am J Public Health 1987;77(2):195-199.
# Adapted from episheet::pvalueplot() (CRAN).

`%||%` <- function(a, b) if (!is.null(a)) a else b

pvalueFunctionsUI <- function(id) {
  ns <- NS(id)

  # Helper to build one curve's input block
  curve_inputs <- function(curve_num, color_hex) {
    pfx   <- paste0("c", curve_num, "_")
    title <- paste0("Curve ", curve_num)
    border_style <- paste0(
      "border-left: 3px solid ", color_hex, "; ",
      "padding: 10px 12px; margin-bottom: 10px; ",
      "background: #fafcfd; border-radius: 0 8px 8px 0;"
    )
    tags$div(
      style = border_style,
      tags$p(
        style = paste0(
          "font-size: 11px; font-weight: 500; color: ", color_hex, "; ",
          "margin: 0 0 8px;"
        ),
        title
      ),
      tags$div(
        style = "display: grid; grid-template-columns: 1fr 1fr; gap: 6px;",
        numericInput(ns(paste0(pfx, "ll90")), "Lower 90% CI", value = NULL, min = 0.001, step = 0.01),
        numericInput(ns(paste0(pfx, "ul90")), "Upper 90% CI", value = NULL, min = 0.001, step = 0.01),
        numericInput(ns(paste0(pfx, "ll95")), "Lower 95% CI", value = NULL, min = 0.001, step = 0.01),
        numericInput(ns(paste0(pfx, "ul95")), "Upper 95% CI", value = NULL, min = 0.001, step = 0.01),
        numericInput(ns(paste0(pfx, "ll99")), "Lower 99% CI", value = NULL, min = 0.001, step = 0.01),
        numericInput(ns(paste0(pfx, "ul99")), "Upper 99% CI", value = NULL, min = 0.001, step = 0.01)
      )
    )
  }

  tagList(
    tags$div(class = "calc-layout",

      # ── Input panel ──────────────────────────────────────────────────────────
      tags$div(class = "panel",
        tags$div(class = "panel-head ph-teal", "Input"),
        tags$div(class = "panel-body",

          tags$p(
            style = "font-size: 11px; color: #4a6070; margin-bottom: 12px;",
            "Enter the lower and upper CI bounds at three confidence levels for each curve. ",
            "The 95% bounds are used to recover the point estimate and SE; ",
            "the 90% and 99% bounds are cross-checked if provided."
          ),

          curve_inputs(1, "#005f6a"),
          curve_inputs(2, "#165c7d"),

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

    # ── Plot (full-width below the two-column layout) ─────────────────────────
    tags$div(class = "panel", style = "margin-top: 16px;",
      tags$div(class = "panel-head ph-blue", "P-value Function Plot"),
      tags$div(class = "panel-body",
        plotOutput(ns("pvalue_plot"), height = "380px")
      )
    ),

    # ── Interpretation ─────────────────────────────────────────────────────────
    uiOutput(ns("interpretation"))
  )
}

pvalueFunctionsServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ── Parse one curve's inputs; returns NULL if 95% bounds not valid ─────────
    parse_curve <- function(pfx) {
      ll95 <- input[[paste0(pfx, "ll95")]]
      ul95 <- input[[paste0(pfx, "ul95")]]
      if (is.null(ll95) || is.null(ul95))   return(NULL)
      if (anyNA(c(ll95, ul95)))              return(NULL)
      if (ll95 <= 0 || ul95 <= 0)           return(NULL)
      if (ul95 <= ll95)                      return(NULL)

      # Recover point estimate and SE from 95% CI
      ln_est <- (log(ul95) + log(ll95)) / 2
      SE     <- (log(ul95) - log(ll95)) / (2 * 1.96)
      est    <- exp(ln_est)

      # Optional bounds for annotation / cross-check
      ll90 <- input[[paste0(pfx, "ll90")]]; ul90 <- input[[paste0(pfx, "ul90")]]
      ll99 <- input[[paste0(pfx, "ll99")]]; ul99 <- input[[paste0(pfx, "ul99")]]

      list(
        est    = est,
        SE     = SE,
        ln_est = ln_est,
        ll95   = ll95,
        ul95   = ul95,
        ll90   = ll90, ul90 = ul90,
        ll99   = ll99, ul99 = ul99
      )
    }

    # ── P-value function for a single RR_0 ────────────────────────────────────
    pval_fn <- function(rr0, ln_est, SE) {
      z <- abs(ln_est - log(rr0)) / SE
      2 * (1 - pnorm(z))
    }

    # ── Core computation ───────────────────────────────────────────────────────
    computed <- eventReactive(input$calc, {
      c1 <- parse_curve("c1_")
      c2 <- parse_curve("c2_")

      validate(
        need(!is.null(c1),
          "Curve 1: enter valid lower and upper 95% CI bounds (both must be positive and LL < UL)."),
        need(!is.null(c2),
          "Curve 2: enter valid lower and upper 95% CI bounds (both must be positive and LL < UL).")
      )

      # Grid of RR_0 values on log scale from 0.01 to 100
      rr_seq <- exp(seq(log(0.01), log(100), length.out = 300))

      p1 <- sapply(rr_seq, pval_fn, ln_est = c1$ln_est, SE = c1$SE)
      p2 <- sapply(rr_seq, pval_fn, ln_est = c2$ln_est, SE = c2$SE)

      df <- data.frame(
        rr   = rep(rr_seq, 2),
        pval = c(p1, p2),
        curve = rep(c("Curve 1", "Curve 2"), each = length(rr_seq))
      )

      list(c1 = c1, c2 = c2, df = df)
    })

    # ── Results summary cards ──────────────────────────────────────────────────
    output$results <- renderUI({
      r  <- computed()
      c1 <- r$c1
      c2 <- r$c2

      tags$div(class = "metrics",
        metric_card(
          "Curve 1 — Point Estimate",
          fmt4(c1$est),
          paste0("SE(ln RR) = ", fmt4(c1$SE),
                 " \u2502 95% CI: ", fmt4(c1$ll95), "\u2013", fmt4(c1$ul95)),
          "m-primary wide",
          tip = "Point estimate recovered from the 95% CI midpoint on the log scale: exp((ln(UL) + ln(LL)) / 2)."
        ),
        metric_card(
          "Curve 2 — Point Estimate",
          fmt4(c2$est),
          paste0("SE(ln RR) = ", fmt4(c2$SE),
                 " \u2502 95% CI: ", fmt4(c2$ll95), "\u2013", fmt4(c2$ul95)),
          "m-blue",
          tip = "Point estimate recovered from the 95% CI midpoint on the log scale: exp((ln(UL) + ln(LL)) / 2)."
        ),
        metric_card(
          "Curve 1 — P at RR = 1",
          fmt4(2 * (1 - pnorm(abs(c1$ln_est) / c1$SE))),
          "two-tailed p for null hypothesis",
          "m-teal",
          tip = "Two-sided p-value at the null value RR = 1. Derived from SE(ln RR) = (ln(UL) - ln(LL)) / (2 \u00d7 1.96)."
        ),
        metric_card(
          "Curve 2 — P at RR = 1",
          fmt4(2 * (1 - pnorm(abs(c2$ln_est) / c2$SE))),
          "two-tailed p for null hypothesis",
          "m-purple",
          tip = "Two-sided p-value at the null value RR = 1. Derived from SE(ln RR) = (ln(UL) - ln(LL)) / (2 \u00d7 1.96)."
        )
      )
    })

    # ── P-value function plot ──────────────────────────────────────────────────
    output$pvalue_plot <- renderPlot({
      r <- computed()

      # Colour mapping: Curve 1 = petrol (#005f6a), Curve 2 = blue (#165c7d)
      colour_map <- c("Curve 1" = "#005f6a", "Curve 2" = "#165c7d")

      ggplot2::ggplot(r$df, ggplot2::aes(x = rr, y = pval, colour = curve)) +
        ggplot2::geom_line(linewidth = 1.2) +
        # Null line at RR = 1
        ggplot2::geom_vline(
          xintercept = 1,
          linetype   = "solid",
          colour     = "#cccccc",
          linewidth  = 0.8
        ) +
        # Horizontal reference lines at key p-values
        ggplot2::geom_hline(yintercept = 0.05, linetype = "dotted", colour = "#aaaaaa", linewidth = 0.5) +
        ggplot2::geom_hline(yintercept = 0.10, linetype = "dotted", colour = "#aaaaaa", linewidth = 0.5) +
        ggplot2::geom_hline(yintercept = 0.01, linetype = "dotted", colour = "#aaaaaa", linewidth = 0.5) +
        # Point estimates as vertical tick
        ggplot2::geom_vline(
          data = data.frame(
            rr    = c(r$c1$est, r$c2$est),
            curve = c("Curve 1", "Curve 2")
          ),
          ggplot2::aes(xintercept = rr, colour = curve),
          linetype  = "dashed",
          linewidth = 0.5,
          alpha     = 0.6
        ) +
        ggplot2::scale_colour_manual(values = colour_map, name = NULL) +
        ggplot2::scale_x_log10(
          limits = c(0.01, 100),
          breaks = c(0.01, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100),
          labels = c("0.01", "0.05", "0.10", "0.20", "0.50",
                     "1.00", "2.00", "5.00", "10", "20", "50", "100")
        ) +
        ggplot2::scale_y_continuous(
          limits = c(0, 1),
          breaks = seq(0, 1, 0.1),
          expand = ggplot2::expansion(mult = c(0, 0.02))
        ) +
        ggplot2::labs(
          x = "Relative Risk (log scale)",
          y = "P-value",
          title = "P-value Functions"
        ) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(
          panel.grid.minor   = ggplot2::element_blank(),
          panel.grid.major.x = ggplot2::element_line(colour = "#e8eef2", linewidth = 0.4),
          panel.grid.major.y = ggplot2::element_line(colour = "#e8eef2", linewidth = 0.4),
          legend.position    = "top",
          legend.text        = ggplot2::element_text(size = 10, colour = "#4a6070"),
          axis.text.x        = ggplot2::element_text(size = 9,  colour = "#4a6070"),
          axis.text.y        = ggplot2::element_text(size = 9,  colour = "#4a6070"),
          axis.title         = ggplot2::element_text(size = 10, colour = "#4a6070"),
          plot.title         = ggplot2::element_text(size = 12, colour = "#005f6a",
                                                     face = "plain"),
          plot.background    = ggplot2::element_rect(fill = "#ffffff", colour = NA),
          panel.background   = ggplot2::element_rect(fill = "#ffffff", colour = NA)
        )
    })

    # ── Interpretation ─────────────────────────────────────────────────────────
    output$interpretation <- renderUI({
      r  <- computed()
      c1 <- r$c1
      c2 <- r$c2

      p1_null <- 2 * (1 - pnorm(abs(c1$ln_est) / c1$SE))
      p2_null <- 2 * (1 - pnorm(abs(c2$ln_est) / c2$SE))

      # Compatibility at RR = 1 (p ≥ 0.05 means compatible with null)
      compat1 <- if (p1_null >= 0.05) "compatible" else "incompatible"
      compat2 <- if (p2_null >= 0.05) "compatible" else "incompatible"

      tags$div(class = "conclusion",
        tags$div(class = "conc-head", "Interpretation"),
        tags$p(class = "conc-text",
          "Each curve shows the two-sided p-value as a function of a hypothesised relative risk. ",
          "Values on the curve where p \u2265 0.05 form the 95% compatibility interval. ",
          "The peak of each curve occurs at the point estimate. ",
          tags$strong("Curve 1"), ": estimate\u202f=\u202f",
          tags$strong(fmt4(c1$est)),
          " (95% CI: ", fmt4(c1$ll95), "\u2013", fmt4(c1$ul95),
          "); p at RR\u202f=\u202f1 is ", tags$strong(fmt4(p1_null)),
          " — the data are ", compat1, " with the null hypothesis. ",
          tags$strong("Curve 2"), ": estimate\u202f=\u202f",
          tags$strong(fmt4(c2$est)),
          " (95% CI: ", fmt4(c2$ll95), "\u2013", fmt4(c2$ul95),
          "); p at RR\u202f=\u202f1 is ", tags$strong(fmt4(p2_null)),
          " — the data are ", compat2, " with the null hypothesis."
        )
      )
    })

    outputOptions(output, "results",        suspendWhenHidden = FALSE)
    outputOptions(output, "interpretation", suspendWhenHidden = FALSE)
  })
}
