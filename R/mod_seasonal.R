# mod_seasonal.R — Seasonal Analysis module
# Walter-Elwood vector method for testing periodicity in monthly count data
# References:
#   Walter SD, Elwood JM. Ann Hum Gen 1975;25:83-87
#   Brookhart MA, Rothman KJ. BMC Med Res Methodol 2008;8:67
#
# Test values (CLAUDE.md): Peak/Low Ratio=1.819; Peak=Dec 17 (346°); 95% CI: 1.153–2.872

seasonalUI <- function(id) {
  ns <- NS(id)

  month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                   "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

  tagList(
    tags$div(class = "calc-layout",

      # ── Input panel ────────────────────────────────────────────────────────────
      tags$div(class = "panel",
        tags$div(class = "panel-head ph-teal", "Input"),
        tags$div(class = "panel-body",

          tags$p(class = "panel-hint",
            "Enter monthly case counts. Population denominators are optional: ",
            "if provided, rates are used instead of raw counts."),

          tags$div(style = "overflow-x: auto;",
            tags$table(
              class = "input-table",
              style = "width: 100%; border-collapse: collapse; font-size: 0.88em;",

              tags$thead(
                tags$tr(
                  tags$th(style = "padding: 4px 8px; text-align: left; min-width: 44px;",
                          "Month"),
                  tags$th(style = "padding: 4px 8px; text-align: right; min-width: 90px;",
                          "Cases"),
                  tags$th(style = "padding: 4px 8px; text-align: right; min-width: 100px;",
                          "Population (opt.)")
                )
              ),

              tags$tbody(
                lapply(1:12, function(i) {
                  tags$tr(
                    # Month label — fixed, not user-editable
                    tags$td(
                      style = "padding: 2px 8px; font-weight: 600; color: #2a7f8f;",
                      month_names[i]
                    ),
                    # Cases
                    tags$td(style = "padding: 2px 4px;",
                      numericInput(ns(paste0("cases_", i)),
                                   label = NULL,
                                   value = NULL, min = 0, step = 1,
                                   width = "90px")
                    ),
                    # Population (optional)
                    tags$td(style = "padding: 2px 4px;",
                      numericInput(ns(paste0("pop_", i)),
                                   label = NULL,
                                   value = NULL, min = 0, step = 1,
                                   width = "100px")
                    )
                  )
                })
              )
            )
          ),

          tags$br(),
          tags$button(
            class   = "calc-btn",
            onclick = sprintf("Shiny.setInputValue('%s', Math.random())", ns("calc")),
            "Calculate"
          )
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

    # ── Full-width output: chart ────────────────────────────────────────────────
    uiOutput(ns("full_width_outputs"))
  )
}


seasonalServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                     "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

    # Helper: single metric card
    metric <- function(label, value, sub = NULL, cls = "m-blue") {
      tags$div(
        class = paste("metric", cls),
        tags$div(class = "metric-label", label),
        tags$div(class = "metric-value", value),
        if (!is.null(sub)) tags$div(class = "metric-sub", sub)
      )
    }

    # ── Helper: convert fractional month to approximate calendar date ──────────
    # peak_month_frac is 1-based (1 = Jan, 12 = Dec)
    frac_month_to_label <- function(frac) {
      # Clamp to [1, 13)
      frac <- ((frac - 1) %% 12) + 1
      month_idx  <- floor(frac)
      month_idx  <- max(1L, min(12L, as.integer(month_idx)))
      # day within month (approximate: 28-31 days, use 30 for simplicity)
      day_frac   <- frac - floor(frac)
      day_approx <- round(day_frac * 30) + 1
      day_approx <- max(1L, min(30L, as.integer(day_approx)))
      sprintf("%s %d", month_names[month_idx], day_approx)
    }

    # ── Core computation ────────────────────────────────────────────────────────
    computed <- eventReactive(input$calc, {

      cases_vec <- vapply(1:12, function(i) {
        v <- input[[paste0("cases_", i)]]
        if (is.null(v) || is.na(v)) NA_real_ else as.numeric(v)
      }, numeric(1))

      pop_vec <- vapply(1:12, function(i) {
        v <- input[[paste0("pop_", i)]]
        if (is.null(v) || is.na(v)) NA_real_ else as.numeric(v)
      }, numeric(1))

      validate(
        need(!all(is.na(cases_vec)),
             "Enter case counts for at least some months."),
        need(sum(!is.na(cases_vec)) >= 6,
             "At least 6 months of case data are required."),
        need(all(is.na(cases_vec) | cases_vec >= 0),
             "Case counts cannot be negative.")
      )

      # Replace missing cases with 0 (absent months contribute no weight)
      cases_vec[is.na(cases_vec)] <- 0

      # Decide whether to use rates or counts
      # Use rates only when ALL 12 population values are present and positive
      use_rates <- !any(is.na(pop_vec)) && all(pop_vec > 0)

      if (use_rates) {
        validate(
          need(all(pop_vec >= cases_vec, na.rm = TRUE),
               "Population cannot be smaller than the case count for any month.")
        )
        n_k <- cases_vec / pop_vec  # rates
        N_total <- sum(n_k)
        y_label <- "Rate (cases / population)"
      } else {
        n_k     <- cases_vec        # raw counts
        N_total <- sum(n_k)
        y_label <- "Case count"
      }

      validate(need(N_total > 0, "Total cases (or sum of rates) must be > 0."))

      # Mid-month angles: theta_k = 2*pi*(k - 0.5) / 12
      months   <- 1:12
      theta_k  <- 2 * pi * (months - 0.5) / 12

      # Vector components
      C_comp  <- sum(n_k * cos(theta_k))
      S_comp  <- sum(n_k * sin(theta_k))

      # Mean resultant length (r)
      r_val   <- sqrt(C_comp^2 + S_comp^2) / N_total

      # Direction of peak (radians, then degrees)
      phi_rad <- atan2(S_comp, C_comp)
      phi_deg <- phi_rad * 180 / pi
      if (phi_deg < 0) phi_deg <- phi_deg + 360

      # Convert angle to fractional month
      # phi = 0 corresponds to mid-January (month 0.5 on a 0-based scale)
      # peak_month_frac (1-based) = 0.5 + phi * 12 / (2*pi)
      peak_month_frac <- 0.5 + phi_rad * 12 / (2 * pi)
      if (peak_month_frac < 1)  peak_month_frac <- peak_month_frac + 12
      if (peak_month_frac > 13) peak_month_frac <- peak_month_frac - 12

      peak_label <- frac_month_to_label(peak_month_frac)

      # Rayleigh test statistic Z = N * r^2
      z_rayleigh <- N_total * r_val^2

      # p-value: chi-squared with 2 df (Rayleigh)
      p_rayleigh <- pchisq(2 * N_total * r_val^2, df = 2, lower.tail = FALSE)

      # Peak / Low ratio
      # Fitted sinusoid: f_k = mean * (1 + 2*r*cos(theta_k - phi))
      # Maximum at theta = phi: f_peak = mean * (1 + 2*r)
      # Minimum at theta = phi + pi: f_low  = mean * (1 - 2*r)
      mean_n        <- N_total / 12
      peak_val      <- mean_n * (1 + 2 * r_val)
      low_val       <- mean_n * (1 - 2 * r_val)
      peak_low_ratio <- if (low_val > 0) peak_val / low_val else Inf

      # 95% CI for r (approximate normal)
      # SE(r) ≈ sqrt((1 - 2*r^2) / (2*N))  [Rayleigh/circular statistics]
      se_r <- sqrt(max(0, (1 - 2 * r_val^2) / (2 * N_total)))
      r_lo95 <- max(0, r_val - 1.96 * se_r)
      r_hi95 <- min(1, r_val + 1.96 * se_r)

      # Transform r CIs to Peak/Low ratio CIs
      pl_lo95 <- if (r_lo95 < 0.5) (1 + 2 * r_lo95) / (1 - 2 * r_lo95) else Inf
      pl_hi95 <- if (r_hi95 < 0.5) (1 + 2 * r_hi95) / (1 - 2 * r_hi95) else Inf

      # 90% CI for r
      r_lo90 <- max(0, r_val - 1.645 * se_r)
      r_hi90 <- min(1, r_val + 1.645 * se_r)
      pl_lo90 <- if (r_lo90 < 0.5) (1 + 2 * r_lo90) / (1 - 2 * r_lo90) else Inf
      pl_hi90 <- if (r_hi90 < 0.5) (1 + 2 * r_hi90) / (1 - 2 * r_hi90) else Inf

      # Fitted sinusoidal values
      fitted_k <- mean_n * (1 + 2 * r_val * cos(theta_k - phi_rad))

      list(
        n_k            = n_k,
        fitted_k       = fitted_k,
        r_val          = r_val,
        phi_deg        = phi_deg,
        peak_label     = peak_label,
        z_rayleigh     = z_rayleigh,
        p_rayleigh     = p_rayleigh,
        peak_low_ratio = peak_low_ratio,
        pl_lo90        = pl_lo90,
        pl_hi90        = pl_hi90,
        pl_lo95        = pl_lo95,
        pl_hi95        = pl_hi95,
        N_total        = N_total,
        use_rates      = use_rates,
        y_label        = y_label
      )
    })

    # ── Metric cards ─────────────────────────────────────────────────────────────
    output$results <- renderUI({
      r <- computed()

      # Format Peak/Low ratio card subtitle with 95% CI
      pl_sub <- if (is.finite(r$pl_lo95) && is.finite(r$pl_hi95)) {
        paste0("95% CI: ", fmt2(r$pl_lo95), "\u2013", fmt2(r$pl_hi95))
      } else {
        "95% CI: \u2014"
      }

      # p-value formatting
      p_fmt <- if (r$p_rayleigh < 0.001) {
        "< 0.001"
      } else {
        fmt4(r$p_rayleigh)
      }

      tags$div(class = "metrics",
        metric(
          label = "Peak / Low Ratio",
          value = fmt2(r$peak_low_ratio),
          sub   = pl_sub,
          cls   = "m-primary wide"
        ),
        metric(
          label = "Time of Peak",
          value = r$peak_label,
          sub   = paste0(round(r$phi_deg, 1), "\u00b0"),
          cls   = "m-blue"
        ),
        metric(
          label = "Rayleigh p-value",
          value = p_fmt,
          sub   = paste0("Z = ", fmt2(r$z_rayleigh)),
          cls   = "m-teal"
        ),
        metric(
          label = "Mean vector length (r)",
          value = fmt4(r$r_val),
          sub   = "0 = uniform; 1 = all one month",
          cls   = "m-purple"
        )
      )
    })

    # ── Full-width: chart + interpretation ────────────────────────────────────
    output$full_width_outputs <- renderUI({
      tagList(
        tags$hr(),
        tags$h4("Monthly Distribution", style = "margin-top: 16px;"),
        plotOutput(ns("seasonal_plot"), height = "380px"),
        tags$div(class = "conclusion",
          tags$div(class = "conc-head", "Interpretation"),
          uiOutput(ns("interpretation"))
        )
      )
    })

    # Petrol colour consistent with other modules
    PETROL  <- "#2a7f8f"
    OBS_COL <- "#3b6ea5"   # blue for observed

    output$seasonal_plot <- renderPlot({
      r <- computed()

      df <- data.frame(
        month   = factor(month.abb, levels = month.abb),
        month_i = 1:12,
        obs     = r$n_k,
        fitted  = r$fitted_k,
        stringsAsFactors = FALSE
      )

      ggplot2::ggplot(df, ggplot2::aes(x = month_i)) +
        # Observed line
        ggplot2::geom_line(
          ggplot2::aes(y = obs, colour = "Observed"),
          linewidth = 1.1
        ) +
        ggplot2::geom_point(
          ggplot2::aes(y = obs, colour = "Observed"),
          size = 2.5
        ) +
        # Fitted sinusoidal curve
        ggplot2::geom_line(
          ggplot2::aes(y = fitted, colour = "Fitted"),
          linewidth = 1.1,
          linetype  = "solid"
        ) +
        ggplot2::scale_colour_manual(
          name   = NULL,
          values = c("Observed" = OBS_COL, "Fitted" = PETROL),
          breaks = c("Observed", "Fitted")
        ) +
        ggplot2::scale_x_continuous(
          breaks = 1:12,
          labels = month.abb,
          name   = "Month"
        ) +
        ggplot2::scale_y_continuous(
          name   = r$y_label,
          expand = ggplot2::expansion(mult = c(0.05, 0.1))
        ) +
        ggplot2::labs(
          title    = "Seasonal Distribution: Observed vs Fitted Sinusoid",
          subtitle = sprintf(
            "r = %.4f  |  Peak: %s (%s\u00b0)  |  Peak/Low ratio = %s  |  Rayleigh p = %s",
            r$r_val,
            r$peak_label,
            round(r$phi_deg, 1),
            fmt2(r$peak_low_ratio),
            if (r$p_rayleigh < 0.001) "< 0.001" else fmt4(r$p_rayleigh)
          )
        ) +
        ggplot2::theme_bw(base_size = 13) +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank(),
          legend.position  = "top",
          legend.direction = "horizontal",
          plot.subtitle    = ggplot2::element_text(size = 9.5, colour = "grey40")
        )
    })

    output$interpretation <- renderUI({
      r <- computed()

      # Strength of seasonality narrative
      strength <- if (r$r_val < 0.1) {
        "very weak"
      } else if (r$r_val < 0.2) {
        "weak"
      } else if (r$r_val < 0.4) {
        "moderate"
      } else {
        "strong"
      }

      sig_text <- if (r$p_rayleigh < 0.05) {
        tags$span("statistically significant (p\u00a0= ",
                  if (r$p_rayleigh < 0.001) "< 0.001" else fmt4(r$p_rayleigh),
                  ")")
      } else {
        tags$span("not statistically significant (p\u00a0=\u00a0",
                  fmt4(r$p_rayleigh), ")")
      }

      pl_ci_text <- if (is.finite(r$pl_lo95) && is.finite(r$pl_hi95)) {
        tagList(
          " The estimated Peak\u2013Low ratio is\u00a0",
          tags$strong(fmt2(r$peak_low_ratio)),
          " (95%\u00a0CI:\u00a0",
          tags$strong(fmt2(r$pl_lo95)), "\u2013", tags$strong(fmt2(r$pl_hi95)),
          "); a ratio of 1.0 would indicate no seasonality."
        )
      } else {
        tagList(" The Peak\u2013Low ratio is indeterminate because r \u2265 0.5.")
      }

      tags$p(class = "conc-text",
        "The Rayleigh (Walter-Elwood) test finds ",
        strength, " seasonality that is ", sig_text, ". ",
        "The mean resultant vector length r\u00a0=\u00a0",
        tags$strong(fmt4(r$r_val)),
        " describes how concentrated the distribution is around the peak month. ",
        "The estimated peak is\u00a0", tags$strong(r$peak_label),
        " (", round(r$phi_deg, 1), "\u00b0).",
        pl_ci_text,
        " The fitted sinusoidal curve assumes a single harmonic (annual) cycle; ",
        "departures of the observed line from the fitted curve indicate ",
        "higher-order periodicities not captured by this model. ",
        tags$em("Reference: Walter & Elwood 1975; Brookhart & Rothman 2008.")
      )
    })

    outputOptions(output, "results",            suspendWhenHidden = FALSE)
    outputOptions(output, "full_width_outputs", suspendWhenHidden = FALSE)
  })
}
