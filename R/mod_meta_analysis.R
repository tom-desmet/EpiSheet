# mod_meta_analysis.R — Meta-Analysis module
# Methods: Inverse-variance (fixed effects) and DerSimonian-Laird (random effects)
# VBA ports used: sum_w, sum_wlnrr, sum_chi2, var_pop, sum_wrlnrr,
#                 sum_wvlnrr, sum_wrvlnrr, inf_w, inf_wr  (all in helpers.R)
#
# Test values: Studies A-D with RR 1.2/2.0/3.0/1.2, CIs spanning 3/7/7/4 points
# → Pooled (fixed) ≈ 1.687, P homogeneity ≈ 0.496

# NULL-coalescing helper (defined once; loaded before server runs)
`%||%` <- function(a, b) if (!is.null(a)) a else b

metaAnalysisUI <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(class = "calc-layout",

      # ── Input panel ──────────────────────────────────────────────────────────
      tags$div(class = "panel",
        tags$div(class = "panel-head ph-teal", "Input"),
        tags$div(class = "panel-body",

          textInput(ns("title"), "Analysis title",
            placeholder = "e.g. Effect of X on Y"),

          tags$div(
            style = "display: grid; grid-template-columns: 1fr 1fr; gap: 10px; margin-bottom: 10px;",
            numericInput(ns("n_studies"), "Number of studies",
              value = 4, min = 2, max = 20, step = 1),
            tags$div(
              tags$label(
                style = "font-size: 11px; font-weight: 500; color: #4a6070; display: block; margin-bottom: 4px;",
                "Model"
              ),
              radioButtons(ns("model"), label = NULL,
                choices  = c("Fixed effects" = "fixed", "Random effects" = "random"),
                selected = "fixed",
                inline   = TRUE)
            )
          ),

          tags$p(
            style = "font-size: 11px; color: #4a6070; margin: 8px 0 6px;",
            "Enter RR and 95% CI limits for each study:"
          ),

          uiOutput(ns("study_inputs")),

          tags$br(),
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

    # ── Forest plot — full-width below the two-column layout ──────────────────
    tags$div(class = "panel", style = "margin-top: 16px;",
      tags$div(class = "panel-head ph-blue", "Forest Plot"),
      tags$div(class = "panel-body",
        plotOutput(ns("forest_plot"), height = "380px")
      )
    ),

    # ── Interpretation ─────────────────────────────────────────────────────────
    uiOutput(ns("interpretation"))
  )
}

metaAnalysisServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ── Dynamic study-input rows (re-rendered when n_studies changes) ──────────
    output$study_inputs <- renderUI({
      k <- max(2L, min(20L, as.integer(input$n_studies %||% 4)))
      tagList(
        lapply(seq_len(k), function(i) {
          tags$div(
            style = paste0(
              "display: grid; ",
              "grid-template-columns: 1.4fr 0.8fr 0.8fr 0.8fr; ",
              "gap: 6px; margin-bottom: 6px; align-items: end;"
            ),
            # Study label
            tags$div(
              if (i == 1) tags$label(
                style = "font-size: 10px; color: #4a6070;", "Study"),
              textInput(ns(paste0("study_", i)), label = NULL,
                placeholder = paste0("Study ", LETTERS[min(i, 26L)]))
            ),
            # RR
            tags$div(
              if (i == 1) tags$label(
                style = "font-size: 10px; color: #4a6070;", "RR"),
              numericInput(ns(paste0("rr_", i)), label = NULL,
                value = NULL, min = 0.001, step = 0.01)
            ),
            # Lower 95% CI
            tags$div(
              if (i == 1) tags$label(
                style = "font-size: 10px; color: #4a6070;", "Lower 95% CI"),
              numericInput(ns(paste0("ll_", i)), label = NULL,
                value = NULL, min = 0.001, step = 0.01)
            ),
            # Upper 95% CI
            tags$div(
              if (i == 1) tags$label(
                style = "font-size: 10px; color: #4a6070;", "Upper 95% CI"),
              numericInput(ns(paste0("ul_", i)), label = NULL,
                value = NULL, min = 0.001, step = 0.01)
            )
          )
        })
      )
    })

    # ── Collect valid study data from inputs ───────────────────────────────────
    get_studies <- function() {
      k <- max(2L, min(20L, as.integer(input$n_studies %||% 4)))
      studies <- lapply(seq_len(k), function(i) {
        rr <- input[[paste0("rr_", i)]]
        ll <- input[[paste0("ll_", i)]]
        ul <- input[[paste0("ul_", i)]]
        lb <- input[[paste0("study_", i)]]

        if (is.null(rr) || is.null(ll) || is.null(ul)) return(NULL)
        if (anyNA(c(rr, ll, ul)))                       return(NULL)
        if (rr <= 0 || ll <= 0 || ul <= 0)              return(NULL)
        if (ul <= ll)                                    return(NULL)

        lnRR <- log(rr)
        v    <- ((log(ul) - log(ll)) / (2 * 1.96))^2

        list(
          label = if (!is.null(lb) && nchar(trimws(lb)) > 0)
                    trimws(lb) else paste0("Study ", i),
          RR    = rr,
          lnRR  = lnRR,
          v     = v,
          ll    = ll,
          ul    = ul
        )
      })
      Filter(Negate(is.null), studies)
    }

    # ── Core computation ───────────────────────────────────────────────────────
    computed <- eventReactive(input$calc, {
      studies <- get_studies()
      validate(need(length(studies) >= 2,
        "Enter a valid RR and 95% CI (LL < UL) for at least 2 studies."))

      lnRR  <- sapply(studies, `[[`, "lnRR")
      v     <- sapply(studies, `[[`, "v")
      k     <- length(studies)
      model <- input$model %||% "fixed"

      # ── Fixed-effects (inverse-variance) pooled estimate
      W_fix     <- sum_w(v)
      lnRR_fix  <- sum_wlnrr(lnRR, v) / W_fix
      var_fix   <- sum_wvlnrr(lnRR, v)          # = 1 / sum_w
      rr_fix    <- exp(lnRR_fix)
      ci_fix_95 <- exp(lnRR_fix + c(-1, 1) * 1.96  * sqrt(var_fix))
      ci_fix_99 <- exp(lnRR_fix + c(-1, 1) * 2.576 * sqrt(var_fix))

      # ── Cochran Q heterogeneity statistic
      Q   <- sum_chi2(lnRR, v)
      p_Q <- pchisq(Q, df = k - 1L, lower.tail = FALSE)
      I2  <- max(0, (Q - (k - 1)) / Q) * 100

      # ── DerSimonian-Laird tau²
      tau2      <- var_pop(lnRR, v)

      # ── Random-effects pooled estimate
      vr          <- v + tau2
      W_rand      <- sum(1 / vr, na.rm = TRUE)
      lnRR_rand   <- sum_wrlnrr(lnRR, v) / W_rand
      var_rand    <- sum_wrvlnrr(lnRR, v)
      rr_rand     <- exp(lnRR_rand)
      ci_rand_95  <- exp(lnRR_rand + c(-1, 1) * 1.96  * sqrt(var_rand))
      ci_rand_99  <- exp(lnRR_rand + c(-1, 1) * 2.576 * sqrt(var_rand))

      # ── Per-study weights (%)
      wt_fix  <- inf_w(v)
      wt_rand <- inf_wr(lnRR, v)

      # ── Select active pooled based on chosen model
      if (model == "fixed") {
        rr_pool   <- rr_fix
        ci95_pool <- ci_fix_95
        ci99_pool <- ci_fix_99
        wt_pool   <- wt_fix
      } else {
        rr_pool   <- rr_rand
        ci95_pool <- ci_rand_95
        ci99_pool <- ci_rand_99
        wt_pool   <- wt_rand
      }

      list(
        studies    = studies,
        k          = k,
        model      = model,
        # Fixed
        rr_fix     = rr_fix,
        ci_fix_95  = ci_fix_95,
        ci_fix_99  = ci_fix_99,
        wt_fix     = wt_fix,
        # Random
        rr_rand    = rr_rand,
        ci_rand_95 = ci_rand_95,
        ci_rand_99 = ci_rand_99,
        wt_rand    = wt_rand,
        tau2       = tau2,
        # Heterogeneity
        Q          = Q,
        p_Q        = p_Q,
        I2         = I2,
        # Active pooled
        rr_pool    = rr_pool,
        ci95_pool  = ci95_pool,
        ci99_pool  = ci99_pool,
        wt_pool    = wt_pool
      )
    })

    # ── Results panel output ───────────────────────────────────────────────────
    output$results <- renderUI({
      r           <- computed()
      model_label <- if (r$model == "fixed") "Fixed effects" else "Random effects (DL)"

      # Per-study weights table
      wt_tbl <- tags$table(
        style = "width: 100%; border-collapse: collapse; font-size: 11px; margin-bottom: 12px;",
        tags$thead(
          tags$tr(
            tags$th(style = "text-align: left;  padding: 4px 6px; color: #4a6070; font-weight: 500; border-bottom: 1px solid #d0dde3;", "Study"),
            tags$th(style = "text-align: right; padding: 4px 6px; color: #4a6070; font-weight: 500; border-bottom: 1px solid #d0dde3;", "RR"),
            tags$th(style = "text-align: right; padding: 4px 6px; color: #4a6070; font-weight: 500; border-bottom: 1px solid #d0dde3;", "95% CI"),
            tags$th(style = "text-align: right; padding: 4px 6px; color: #4a6070; font-weight: 500; border-bottom: 1px solid #d0dde3;", "Weight %"),
            tags$th(style = "text-align: right; padding: 4px 6px; color: #4a6070; font-weight: 500; border-bottom: 1px solid #d0dde3;", "ln(RR)"),
            tags$th(style = "text-align: right; padding: 4px 6px; color: #4a6070; font-weight: 500; border-bottom: 1px solid #d0dde3;", "V(ln RR)")
          )
        ),
        tags$tbody(
          lapply(seq_len(r$k), function(i) {
            s <- r$studies[[i]]
            tags$tr(
              tags$td(style = "padding: 3px 6px;",                          s$label),
              tags$td(style = "text-align: right; padding: 3px 6px;",       fmt4(s$RR)),
              tags$td(style = "text-align: right; padding: 3px 6px;",
                paste0(fmt4(s$ll), "\u2013", fmt4(s$ul))),
              tags$td(style = "text-align: right; padding: 3px 6px;",       fmt2(r$wt_pool[i])),
              tags$td(style = "text-align: right; padding: 3px 6px;",       fmt4(s$lnRR)),
              tags$td(style = "text-align: right; padding: 3px 6px;",       fmt4(s$v))
            )
          })
        )
      )

      tagList(
        wt_tbl,
        tags$div(class = "metrics",
          metric(
            paste0("Pooled RR (", model_label, ")"),
            fmt4(r$rr_pool),
            paste0("95% CI: ", fmt4(r$ci95_pool[1]), "\u2013", fmt4(r$ci95_pool[2])),
            "m-primary wide"
          ),
          metric(
            "99% Confidence Interval",
            paste0(fmt4(r$ci99_pool[1]), "\u2013", fmt4(r$ci99_pool[2])),
            cls = "m-blue"
          ),
          metric(
            "Cochran Q (heterogeneity)",
            paste0("Q = ", fmt2(r$Q)),
            paste0("p = ", fmt4(r$p_Q), "  (df = ", r$k - 1L, ")"),
            cls = "m-teal"
          ),
          metric(
            paste0("I\u00b2 (inconsistency)"),
            paste0(fmt2(r$I2), "%"),
            "% variation due to heterogeneity",
            cls = "m-teal"
          ),
          metric(
            "\u03c4\u00b2 (between-study variance)",
            fmt4(r$tau2),
            "DerSimonian-Laird estimate",
            cls = "m-purple"
          )
        )
      )
    })

    # ── Forest plot ────────────────────────────────────────────────────────────
    output$forest_plot <- renderPlot({
      r <- computed()

      k      <- r$k
      labels <- sapply(r$studies, `[[`, "label")
      rr_v   <- sapply(r$studies, `[[`, "RR")
      ll_v   <- sapply(r$studies, `[[`, "ll")
      ul_v   <- sapply(r$studies, `[[`, "ul")
      wt_v   <- r$wt_pool

      # Studies: y = k down to 1 (top to bottom)
      y_studies <- seq(k, 1L, by = -1L)

      df_studies <- data.frame(
        y     = y_studies,
        label = labels,
        rr    = rr_v,
        ll    = ll_v,
        ul    = ul_v,
        wt    = wt_v,
        type  = "study",
        stringsAsFactors = FALSE
      )

      pooled_label <- paste0("Pooled (", if (r$model == "fixed") "FE" else "RE", ")")
      df_pool <- data.frame(
        y     = 0,
        label = pooled_label,
        rr    = r$rr_pool,
        ll    = r$ci95_pool[1],
        ul    = r$ci95_pool[2],
        wt    = NA_real_,
        type  = "pooled",
        stringsAsFactors = FALSE
      )

      df <- rbind(df_studies, df_pool)

      y_breaks <- c(y_studies, 0L)
      y_labels <- c(as.character(labels), pooled_label)

      # Compute x-axis limits with a little padding
      x_lo <- min(0.1,  min(df$ll, na.rm = TRUE) * 0.75)
      x_hi <- max(10.0, max(df$ul, na.rm = TRUE) * 1.25)

      ggplot2::ggplot(df, ggplot2::aes(x = rr, y = y)) +
        # Null line at RR = 1
        ggplot2::geom_vline(
          xintercept = 1, linetype = "dashed",
          colour = "#888888", linewidth = 0.5
        ) +
        # Separator line above pooled row
        ggplot2::geom_hline(
          yintercept = 0.5, colour = "#d0dde3", linewidth = 0.5
        ) +
        # CI lines for studies
        ggplot2::geom_errorbarh(
          data = df[df$type == "study", ],
          ggplot2::aes(xmin = ll, xmax = ul),
          height = 0.22, colour = "#165c7d", linewidth = 0.7
        ) +
        # Point estimates for studies (size ∝ weight)
        ggplot2::geom_point(
          data = df[df$type == "study", ],
          ggplot2::aes(size = wt),
          colour = "#165c7d", shape = 15
        ) +
        # CI line for pooled
        ggplot2::geom_errorbarh(
          data = df[df$type == "pooled", ],
          ggplot2::aes(xmin = ll, xmax = ul),
          height = 0, colour = "#005f6a", linewidth = 1.2
        ) +
        # Pooled estimate diamond
        ggplot2::geom_point(
          data = df[df$type == "pooled", ],
          colour = "#005f6a", size = 4, shape = 18
        ) +
        ggplot2::scale_x_log10(
          limits = c(x_lo, x_hi),
          breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10),
          labels = c("0.10", "0.20", "0.50", "1.00", "2.00", "5.00", "10.00")
        ) +
        ggplot2::scale_y_continuous(
          breaks = y_breaks,
          labels = y_labels,
          expand = ggplot2::expansion(add = 0.6)
        ) +
        ggplot2::scale_size_continuous(range = c(2, 6), guide = "none") +
        ggplot2::labs(
          x     = "Relative Risk (log scale)",
          y     = NULL,
          title = if (!is.null(input$title) && nchar(trimws(input$title)) > 0)
                    trimws(input$title) else "Forest Plot"
        ) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(
          panel.grid.major.y = ggplot2::element_blank(),
          panel.grid.minor   = ggplot2::element_blank(),
          panel.grid.major.x = ggplot2::element_line(colour = "#e8eef2", linewidth = 0.4),
          axis.text.y        = ggplot2::element_text(size = 10, colour = "#1a1a1a"),
          axis.text.x        = ggplot2::element_text(size = 9,  colour = "#4a6070"),
          axis.title.x       = ggplot2::element_text(size = 10, colour = "#4a6070"),
          plot.title         = ggplot2::element_text(size = 12, colour = "#005f6a",
                                                     face = "plain"),
          plot.background    = ggplot2::element_rect(fill = "#ffffff", colour = NA),
          panel.background   = ggplot2::element_rect(fill = "#ffffff", colour = NA)
        )
    })

    # ── Interpretation ─────────────────────────────────────────────────────────
    output$interpretation <- renderUI({
      r           <- computed()
      model_label <- if (r$model == "fixed") "fixed-effects" else "random-effects (DerSimonian-Laird)"
      sig_het     <- r$p_Q < 0.10

      tags$div(class = "conclusion",
        tags$div(class = "conc-head", "Interpretation"),
        tags$p(class = "conc-text",
          "The ", model_label, " pooled relative risk across ",
          tags$strong(r$k), " studies is ",
          tags$strong(fmt4(r$rr_pool)),
          " (95% CI: ", fmt4(r$ci95_pool[1]), "\u2013", fmt4(r$ci95_pool[2]),
          "; 99% CI: ", fmt4(r$ci99_pool[1]), "\u2013", fmt4(r$ci99_pool[2]), "). ",
          "Cochran Q\u202f=\u202f", tags$strong(fmt2(r$Q)),
          " (p\u202f=\u202f", fmt4(r$p_Q), ", df\u202f=\u202f", r$k - 1L, "); ",
          "I\u00b2\u202f=\u202f", tags$strong(paste0(fmt2(r$I2), "%")), "; ",
          "\u03c4\u00b2\u202f=\u202f", tags$strong(fmt4(r$tau2)), ". ",
          if (sig_het)
            "There is evidence of statistical heterogeneity across studies (p < 0.10); consider the random-effects model or investigate sources of heterogeneity."
          else
            "There is no significant statistical heterogeneity across studies (p \u2265 0.10)."
        )
      )
    })

    outputOptions(output, "results",        suspendWhenHidden = FALSE)
    outputOptions(output, "interpretation", suspendWhenHidden = FALSE)
  })
}
