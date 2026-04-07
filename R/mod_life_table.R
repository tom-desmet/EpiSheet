# mod_life_table.R — Actuarial Life Table module
# Method: Actuarial (life table) method with Greenwood variance and log-log CI
# Reference: Rothman KJ, J Chron Dis 1978;31:557-560

lifeTableUI <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(class = "calc-layout",

      # ── Input panel ──────────────────────────────────────────────────────────
      tags$div(class = "panel",
        tags$div(class = "panel-head ph-teal", "Input"),
        tags$div(class = "panel-body",

          tags$p(class = "panel-hint",
            "Enter up to 10 time periods. Leave rows blank to exclude them."),

          tags$div(style = "overflow-x: auto;",
            tags$table(
              class = "input-table",
              style = "width: 100%; border-collapse: collapse; font-size: 0.88em;",

              tags$thead(
                tags$tr(
                  tags$th(style = "padding: 4px 8px; text-align: left;",  "Period"),
                  tags$th(style = "padding: 4px 8px; text-align: right;", "At risk (n\u2093)"),
                  tags$th(style = "padding: 4px 8px; text-align: right;", "Events (d\u2093)"),
                  tags$th(style = "padding: 4px 8px; text-align: right;", "Withdrawn (w\u2093)")
                )
              ),

              tags$tbody(
                lapply(1:10, function(i) {
                  tags$tr(
                    # Period label — editable
                    tags$td(style = "padding: 2px 4px;",
                      textInput(ns(paste0("lbl_", i)),
                                label  = NULL,
                                value  = as.character(i),
                                width  = "60px")
                    ),
                    # n_x
                    tags$td(style = "padding: 2px 4px;",
                      numericInput(ns(paste0("nx_", i)),
                                   label = NULL,
                                   value = NULL, min = 0, step = 1,
                                   width = "90px")
                    ),
                    # d_x
                    tags$td(style = "padding: 2px 4px;",
                      numericInput(ns(paste0("dx_", i)),
                                   label = NULL,
                                   value = NULL, min = 0, step = 1,
                                   width = "90px")
                    ),
                    # w_x
                    tags$td(style = "padding: 2px 4px;",
                      numericInput(ns(paste0("wx_", i)),
                                   label = NULL,
                                   value = NULL, min = 0, step = 1,
                                   width = "90px")
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

      # ── Results panel ─────────────────────────────────────────────────────────
      tags$div(class = "panel",
        tags$div(class = "panel-head ph-blue", "Results"),
        tags$div(class = "panel-body",
          uiOutput(ns("results"))
        )
      )
    ),

    # ── Full-width outputs: table + chart ─────────────────────────────────────
    uiOutput(ns("full_width_outputs"))
  )
}


lifeTableServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Helper: single metric card
    metric <- function(label, value, sub = NULL, cls = "m-blue") {
      tags$div(
        class = paste("metric", cls),
        tags$div(class = "metric-label", label),
        tags$div(class = "metric-value", value),
        if (!is.null(sub)) tags$div(class = "metric-sub", sub)
      )
    }

    # ── Core computation ───────────────────────────────────────────────────────
    computed <- eventReactive(input$calc, {

      # Collect rows
      rows <- lapply(1:10, function(i) {
        nx  <- input[[paste0("nx_", i)]]
        dx  <- input[[paste0("dx_", i)]]
        wx  <- input[[paste0("wx_", i)]]
        lbl <- input[[paste0("lbl_", i)]]
        list(lbl = lbl, nx = nx, dx = dx, wx = wx)
      })

      # Keep only rows where n_x and d_x are non-missing and valid
      keep <- vapply(rows, function(r) {
        !is.null(r$nx) && !is.na(r$nx) &&
          !is.null(r$dx) && !is.na(r$dx) &&
          r$nx >= 0 && r$dx >= 0
      }, logical(1))

      validate(need(sum(keep) >= 2, "Enter at least 2 complete periods (n\u2093 and d\u2093)."))

      rows <- rows[keep]
      n_periods <- length(rows)

      lbl_vec <- vapply(rows, function(r) as.character(r$lbl), character(1))
      nx_vec  <- vapply(rows, function(r) r$nx, numeric(1))
      dx_vec  <- vapply(rows, function(r) r$dx, numeric(1))
      wx_vec  <- vapply(rows, function(r) {
        if (is.null(r$wx) || is.na(r$wx)) 0 else r$wx
      }, numeric(1))

      validate(
        need(all(dx_vec <= nx_vec), "Events cannot exceed persons at risk."),
        need(all(wx_vec >= 0),      "Withdrawals cannot be negative.")
      )

      # Effective person-time denominator (actuarial adjustment)
      n_eff <- nx_vec - wx_vec / 2

      validate(need(all(n_eff > 0), "Effective n (n\u2093 - w\u2093/2) must be > 0 for all periods."))

      # Period probability of event and survival
      qx <- dx_vec / n_eff          # conditional probability of event
      qx <- pmin(qx, 1)             # cap at 1 (handles edge cases)
      px <- 1 - qx                  # period survival probability

      # Cumulative survival S(t)
      St <- cumprod(px)

      # Greenwood variance components: d_j / (n_eff_j * (n_eff_j - d_j))
      # Guard against division by zero when n_eff == d_x (complete event)
      gw_terms <- ifelse(
        n_eff * (n_eff - dx_vec) > 0,
        dx_vec / (n_eff * (n_eff - dx_vec)),
        0
      )
      gw_cumsum <- cumsum(gw_terms)

      # Greenwood variance of S(t)
      var_St <- St^2 * gw_cumsum

      # 90% CI on log(-log(S)) scale (Rothman 1978)
      z90 <- qnorm(0.95)   # 1.645 for 90% two-sided (z for 5% each tail)

      # For periods where S == 1 or S == 0, CI is degenerate
      loglogS <- ifelse(St > 0 & St < 1, log(-log(St)), NA_real_)

      # Var of log(-log(S)):  gw_cumsum / (log(S))^2
      var_ll <- ifelse(
        !is.na(loglogS) & log(St)^2 > 0,
        gw_cumsum / log(St)^2,
        NA_real_
      )

      lo90 <- ifelse(!is.na(loglogS),
                     exp(-exp(loglogS + z90 * sqrt(var_ll))),
                     ifelse(St == 1, 1, 0))
      hi90 <- ifelse(!is.na(loglogS),
                     exp(-exp(loglogS - z90 * sqrt(var_ll))),
                     ifelse(St == 1, 1, 0))

      # Clamp CI to [0, 1]
      lo90 <- pmax(0, pmin(1, lo90))
      hi90 <- pmax(0, pmin(1, hi90))

      # Median survival: first period where S(t) <= 0.5
      median_idx <- which(St <= 0.5)
      median_surv <- if (length(median_idx) > 0) {
        lbl_vec[median_idx[1]]
      } else {
        "Not reached"
      }

      n_events   <- sum(dx_vec)
      final_surv <- St[n_periods]

      # Data frame for table and chart
      tbl <- data.frame(
        Period        = lbl_vec,
        n_x           = nx_vec,
        d_x           = dx_vec,
        w_x           = wx_vec,
        q_x           = round(qx,  4),
        p_x           = round(px,  4),
        `S(t)`        = round(St,  4),
        `Lower 90% CI` = round(lo90, 4),
        `Upper 90% CI` = round(hi90, 4),
        check.names   = FALSE,
        stringsAsFactors = FALSE
      )

      # Period index for plotting (numeric, used if labels are not numeric)
      period_num <- seq_len(n_periods)

      list(
        tbl         = tbl,
        St          = St,
        lo90        = lo90,
        hi90        = hi90,
        lbl_vec     = lbl_vec,
        period_num  = period_num,
        n_events    = n_events,
        final_surv  = final_surv,
        median_surv = median_surv,
        n_periods   = n_periods
      )
    })

    # ── Metric cards (results panel) ──────────────────────────────────────────
    output$results <- renderUI({
      r <- computed()
      tags$div(class = "metrics",
        metric(
          label = "Median Survival Period",
          value = r$median_surv,
          sub   = "S(t) first \u2264 0.5",
          cls   = "m-primary wide"
        ),
        metric(
          label = "Survival at Last Period",
          value = fmt4(r$final_surv),
          sub   = paste0("Period: ", r$lbl_vec[r$n_periods]),
          cls   = "m-blue"
        ),
        metric(
          label = "Total Events",
          value = as.character(r$n_events),
          sub   = paste0("over ", r$n_periods, " periods"),
          cls   = "m-teal"
        )
      )
    })

    # ── Full-width table + chart ───────────────────────────────────────────────
    output$full_width_outputs <- renderUI({
      tagList(
        tags$hr(),
        tags$h4("Life Table", style = "margin-top: 16px;"),
        tableOutput(ns("surv_table")),
        tags$h4("Survival Curve with 90% Confidence Band",
                style = "margin-top: 24px;"),
        plotOutput(ns("surv_plot"), height = "400px"),
        tags$div(class = "conclusion",
          tags$div(class = "conc-head", "Interpretation"),
          uiOutput(ns("interpretation"))
        )
      )
    })

    output$surv_table <- renderTable({
      computed()$tbl
    }, digits = 4, striped = TRUE, hover = TRUE, bordered = TRUE,
       na = "\u2014")

    # Petrol / teal colour used across the app: #2a7f8f
    PETROL <- "#2a7f8f"
    CI_FILL <- "#c8e6e9"  # light teal for CI band

    output$surv_plot <- renderPlot({
      r <- computed()

      # Build step-function data: add period 0 as S = 1 anchor
      time_pts <- c(0, r$period_num)
      St_pts   <- c(1, r$St)
      lo_pts   <- c(1, r$lo90)
      hi_pts   <- c(1, r$hi90)
      lbl_pts  <- c("0", r$lbl_vec)

      df_plot <- data.frame(
        t    = time_pts,
        St   = St_pts,
        lo90 = lo_pts,
        hi90 = hi_pts,
        lbl  = lbl_pts,
        stringsAsFactors = FALSE
      )

      ggplot2::ggplot(df_plot, ggplot2::aes(x = t)) +
        # CI ribbon (step)
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = lo90, ymax = hi90),
          fill = CI_FILL, alpha = 0.7
        ) +
        # Survival step function
        ggplot2::geom_step(
          ggplot2::aes(y = St),
          colour    = PETROL,
          linewidth = 1.2
        ) +
        # Reference line at S = 0.5
        ggplot2::geom_hline(
          yintercept = 0.5,
          linetype   = "dashed",
          colour     = "grey55",
          linewidth  = 0.5
        ) +
        ggplot2::scale_x_continuous(
          breaks = r$period_num,
          labels = r$lbl_vec,
          name   = "Time Period"
        ) +
        ggplot2::scale_y_continuous(
          limits = c(0, 1),
          breaks = seq(0, 1, 0.1),
          labels = scales::percent_format(accuracy = 1),
          name   = "Survival Probability S(t)"
        ) +
        ggplot2::labs(
          title    = "Actuarial Survival Curve",
          subtitle = "Shaded band = 90% confidence limits (log-log scale); dashed line = S = 0.50"
        ) +
        ggplot2::theme_bw(base_size = 13) +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank(),
          plot.subtitle    = ggplot2::element_text(size = 10, colour = "grey40")
        )
    })

    output$interpretation <- renderUI({
      r <- computed()
      tags$p(class = "conc-text",
        "The actuarial life table covers ", tags$strong(r$n_periods), " time period(s) ",
        "with a total of ", tags$strong(r$n_events), " event(s). ",
        "Cumulative survival at the final period (", tags$strong(r$lbl_vec[r$n_periods]),
        ") is ", tags$strong(fmt4(r$final_surv)), ". ",
        if (r$median_surv == "Not reached") {
          "The median survival time was not reached within the observed follow-up."
        } else {
          tagList(
            "The median survival period (first period with S(t) \u2264 0.50) is ",
            tags$strong(r$median_surv), "."
          )
        },
        " Confidence limits are 90% intervals derived using the log-log transformation ",
        "of the Greenwood variance (Rothman 1978)."
      )
    })

    outputOptions(output, "results",            suspendWhenHidden = FALSE)
    outputOptions(output, "full_width_outputs", suspendWhenHidden = FALSE)
  })
}
