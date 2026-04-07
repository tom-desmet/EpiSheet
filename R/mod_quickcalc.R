# mod_quickcalc.R — Quick Calculator module

quickcalcUI <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(class = "calc-layout",

      # ── Input panel ────────────────────────────────────────────────────────────
      tags$div(class = "panel",
        tags$div(class = "panel-head ph-teal", "Input"),
        tags$div(class = "panel-body",

          tabsetPanel(id = ns("calc_type"), type = "tabs",

            # ── Tab 1: Proportion ──────────────────────────────────────────────
            tabPanel("Proportion",
              tags$div(style = "padding-top: 12px;",
                numericInput(ns("prop_x"), "Numerator (x)", value = NULL, min = 0, step = 1),
                numericInput(ns("prop_N"), "Denominator (N)", value = NULL, min = 1, step = 1),
                numericInput(ns("prop_conf"), "Confidence Level (%)", value = 95, min = 50, max = 99.9, step = 1)
              )
            ),

            # ── Tab 2: Rate / SMR ──────────────────────────────────────────────
            tabPanel("Rate / SMR",
              tags$div(style = "padding-top: 12px;",
                numericInput(ns("rate_x"), "Numerator (events)", value = NULL, min = 0, step = 1),
                numericInput(ns("rate_N"), "Denominator (person-time or expected)", value = NULL, min = 0),
                numericInput(ns("rate_conf"), "Confidence Level (%)", value = 95, min = 50, max = 99.9, step = 1)
              )
            ),

            # ── Tab 3: Statistical Tests ───────────────────────────────────────
            tabPanel("Statistical Tests",
              tags$div(style = "padding-top: 12px;",
                selectInput(ns("test_type"), "Test type",
                  choices = c("Standard Normal (Z)" = "z", "Chi-Square" = "chisq")),
                numericInput(ns("test_val"), "Test statistic value", value = NULL),
                numericInput(ns("test_df"), "Degrees of freedom (Chi-sq only)", value = 1, min = 1, step = 1)
              )
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

quickcalcServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    metric <- function(label, value, sub = NULL, cls = "m-blue") {
      tags$div(class = paste("metric", cls),
        tags$div(class = "metric-label", label),
        tags$div(class = "metric-value", value),
        if (!is.null(sub)) tags$div(class = "metric-sub", sub)
      )
    }

    computed <- eventReactive(input$calc, {
      tab <- input$calc_type

      if (tab == "Proportion") {
        x    <- input$prop_x
        N    <- input$prop_N
        conf <- input$prop_conf
        validate(
          need(!is.null(x) && !is.null(N), "Enter numerator and denominator."),
          need(N > 0, "Denominator must be > 0."),
          need(x >= 0 && x <= N, "Numerator must be between 0 and N."),
          need(!is.null(conf) && conf > 0 && conf < 100, "Confidence level must be between 0 and 100.")
        )
        conflevel <- conf / 100
        alpha     <- 1 - conflevel

        pt_est <- x / N

        # Exact CI — Mid-P
        lo_midp <- binomialcalc(N, x, conflevel, "lower", fisher = FALSE)
        hi_midp <- binomialcalc(N, x, conflevel, "upper", fisher = FALSE)

        # Exact CI — Fisher
        lo_fish <- binomialcalc(N, x, conflevel, "lower", fisher = TRUE)
        hi_fish <- binomialcalc(N, x, conflevel, "upper", fisher = TRUE)

        # Wilson approximate CI
        z_w  <- qnorm(1 - alpha / 2)
        denom_w <- N + z_w^2
        centre  <- (x + z_w^2 / 2) / denom_w
        margin  <- z_w * sqrt(N * (x / N) * (1 - x / N) + z_w^2 / 4) / denom_w
        lo_wils <- centre - margin
        hi_wils <- centre + margin

        # Exact 2-tail P-value (H0: p = 0.5 as a default; use mid-p and Fisher exact vs p=x/N is degenerate,
        # so we follow Episheet convention: 2-sided P that observed count or more extreme under H0: p=0.5)
        # Actually the Episheet Quickcalc gives P that proportion = 0 (i.e. test p=0); we report
        # the mid-p and Fisher 2-sided p-value testing H0: p = 0 vs x observed.
        # More usefully, report the exact 2-sided p for H0: p=0 isn't sensible.
        # Standard approach: 2-sided p = 2 * min(lower tail, upper tail) for mid-p and Fisher.
        p_midp <- 2 * min(
          pbinom(x, N, 0.5) - dbinom(x, N, 0.5) / 2,
          1 - pbinom(x - 1, N, 0.5) - dbinom(x, N, 0.5) / 2
        )
        p_fish <- 2 * min(pbinom(x, N, 0.5), 1 - pbinom(x - 1, N, 0.5))

        list(
          type     = "proportion",
          pt_est   = pt_est,
          conf     = conf,
          lo_midp  = lo_midp, hi_midp = hi_midp,
          lo_fish  = lo_fish, hi_fish = hi_fish,
          lo_wils  = lo_wils, hi_wils = hi_wils,
          p_midp   = p_midp,  p_fish  = p_fish,
          x = x, N = N
        )

      } else if (tab == "Rate / SMR") {
        x    <- input$rate_x
        N    <- input$rate_N
        conf <- input$rate_conf
        validate(
          need(!is.null(x) && !is.null(N), "Enter numerator and denominator."),
          need(N > 0, "Denominator must be > 0."),
          need(x >= 0, "Numerator must be \u2265 0."),
          need(!is.null(conf) && conf > 0 && conf < 100, "Confidence level must be between 0 and 100.")
        )
        conflevel <- conf / 100
        alpha     <- 1 - conflevel

        pt_est <- x / N

        # Exact CI — Mid-P (Poisson)
        lo_midp <- poissoncalc(N, x, conflevel, "lower", fisher = FALSE)
        hi_midp <- poissoncalc(N, x, conflevel, "upper", fisher = FALSE)

        # Exact CI — Fisher/Garwood
        lo_fish <- poissoncalc(N, x, conflevel, "lower", fisher = TRUE)
        hi_fish <- poissoncalc(N, x, conflevel, "upper", fisher = TRUE)

        # Byar approximate CI  (per Rothman: (x*(1 - 1/(9x) ± z/sqrt(9x))^3) / N)
        z_b     <- qnorm(1 - alpha / 2)
        if (x > 0) {
          lo_byar <- x * (1 - 1 / (9 * x) - z_b / sqrt(9 * x))^3 / N
        } else {
          lo_byar <- 0
        }
        hi_byar <- (x + 1) * (1 - 1 / (9 * (x + 1)) + z_b / sqrt(9 * (x + 1)))^3 / N

        # 2-sided p-value (Poisson, mid-p and Fisher/exact)
        mu0 <- N * pt_est  # expected under H0 = observed rate — not useful; use H0: lambda=0 not standard
        # Follow Episheet: two-sided p from Poisson (testing lambda vs observed)
        # mid-p: P(X <= x) - P(X=x)/2
        if (N > 0 && x >= 0) {
          # Use chi-sq equivalence: Fisher 2-sided = 2 * min tail
          p_fish <- 2 * min(ppois(x, x), 1 - ppois(x - 1, x))
          p_midp <- 2 * min(
            ppois(x, x) - dpois(x, x) / 2,
            1 - ppois(x - 1, x) - dpois(x, x) / 2
          )
        } else {
          p_fish <- NA_real_
          p_midp <- NA_real_
        }

        list(
          type     = "rate",
          pt_est   = pt_est,
          conf     = conf,
          lo_midp  = lo_midp, hi_midp = hi_midp,
          lo_fish  = lo_fish, hi_fish = hi_fish,
          lo_byar  = lo_byar, hi_byar = hi_byar,
          p_midp   = p_midp,  p_fish  = p_fish,
          x = x, N = N
        )

      } else {
        # Statistical tests
        val  <- input$test_val
        df   <- input$test_df
        type <- input$test_type
        validate(need(!is.null(val), "Enter a test statistic value."))
        if (type == "z") {
          p_two  <- 2 * pnorm(-abs(val))
          p_upper <- pnorm(val, lower.tail = FALSE)
          p_lower <- pnorm(val)
          list(type = "test", test_type = "z",  val = val, df = NULL,
               p_two = p_two, p_upper = p_upper, p_lower = p_lower)
        } else {
          validate(need(!is.null(df) && df >= 1, "Enter degrees of freedom \u2265 1."))
          p_upper <- pchisq(val, df, lower.tail = FALSE)
          list(type = "test", test_type = "chisq", val = val, df = df,
               p_two = NA, p_upper = p_upper, p_lower = NULL)
        }
      }
    })

    output$results <- renderUI({
      r <- computed()

      if (r$type == "proportion") {
        tags$div(class = "metrics",
          metric("Point Estimate",
                 fmt4(r$pt_est),
                 paste0(r$x, " / ", r$N),
                 "m-primary wide"),
          metric(paste0("Exact Mid-P CI (", r$conf, "%)"),
                 paste0(fmt4(r$lo_midp), " \u2013 ", fmt4(r$hi_midp)),
                 cls = "m-blue"),
          metric(paste0("Exact Fisher CI (", r$conf, "%)"),
                 paste0(fmt4(r$lo_fish), " \u2013 ", fmt4(r$hi_fish)),
                 cls = "m-blue"),
          metric(paste0("Wilson Approx CI (", r$conf, "%)"),
                 paste0(fmt4(r$lo_wils), " \u2013 ", fmt4(r$hi_wils)),
                 cls = "m-teal"),
          metric("2-tail P (Mid-P, H\u2080: p = 0.5)", fmt4(r$p_midp), cls = "m-purple"),
          metric("2-tail P (Fisher, H\u2080: p = 0.5)", fmt4(r$p_fish), cls = "m-purple")
        )

      } else if (r$type == "rate") {
        tags$div(class = "metrics",
          metric("Rate (per unit time)",
                 fmt4(r$pt_est),
                 paste0(r$x, " events / ", r$N, " person-time"),
                 "m-primary wide"),
          metric(paste0("Exact Mid-P CI (", r$conf, "%)"),
                 paste0(fmt4(r$lo_midp), " \u2013 ", fmt4(r$hi_midp)),
                 cls = "m-blue"),
          metric(paste0("Exact Garwood CI (", r$conf, "%)"),
                 paste0(fmt4(r$lo_fish), " \u2013 ", fmt4(r$hi_fish)),
                 cls = "m-blue"),
          metric(paste0("Byar Approx CI (", r$conf, "%)"),
                 paste0(fmt4(r$lo_byar), " \u2013 ", fmt4(r$hi_byar)),
                 cls = "m-teal"),
          metric("2-tail P (Mid-P)", fmt4(r$p_midp), cls = "m-purple"),
          metric("2-tail P (Fisher/Exact)", fmt4(r$p_fish), cls = "m-purple")
        )

      } else {
        # Statistical test results
        if (r$test_type == "z") {
          tags$div(class = "metrics",
            metric("Standard Normal (Z)",
                   fmt4(r$val), "Test statistic", "m-primary wide"),
            metric("2-sided P-value",  fmt4(r$p_two),   cls = "m-blue"),
            metric("Upper-tail P",     fmt4(r$p_upper), cls = "m-teal"),
            metric("Lower-tail P",     fmt4(r$p_lower), cls = "m-purple")
          )
        } else {
          tags$div(class = "metrics",
            metric("Chi-Square",
                   fmt4(r$val),
                   paste0("df = ", r$df),
                   "m-primary wide"),
            metric("Upper-tail P-value", fmt4(r$p_upper), cls = "m-blue")
          )
        }
      }
    })

    output$interpretation <- renderUI({
      r <- computed()
      if (r$type == "proportion") {
        tags$div(class = "conclusion",
          tags$div(class = "conc-head", "Interpretation"),
          tags$p(class = "conc-text",
            "The observed proportion is ", tags$strong(fmt4(r$pt_est)),
            " (", r$x, "/", r$N, "). ",
            "The ", r$conf, "% exact Mid-P confidence interval is ",
            tags$strong(paste0(fmt4(r$lo_midp), "\u2013", fmt4(r$hi_midp))),
            " and the Fisher exact interval is ",
            tags$strong(paste0(fmt4(r$lo_fish), "\u2013", fmt4(r$hi_fish))), "."
          )
        )
      } else if (r$type == "rate") {
        tags$div(class = "conclusion",
          tags$div(class = "conc-head", "Interpretation"),
          tags$p(class = "conc-text",
            "The observed rate is ", tags$strong(fmt4(r$pt_est)), " per unit person-time. ",
            "The ", r$conf, "% exact Mid-P confidence interval is ",
            tags$strong(paste0(fmt4(r$lo_midp), "\u2013", fmt4(r$hi_midp))),
            " and the Garwood (Fisher exact) interval is ",
            tags$strong(paste0(fmt4(r$lo_fish), "\u2013", fmt4(r$hi_fish))), "."
          )
        )
      } else {
        if (r$test_type == "z") {
          sig <- r$p_two < 0.05
          tags$div(class = "conclusion",
            tags$div(class = "conc-head", "Interpretation"),
            tags$p(class = "conc-text",
              "Z = ", tags$strong(fmt4(r$val)),
              ". The two-sided p-value is ", tags$strong(fmt4(r$p_two)), ". ",
              "The result is ", tags$strong(if (sig) "statistically significant" else "not statistically significant"),
              " at the \u03b1 = 0.05 level."
            )
          )
        } else {
          sig <- r$p_upper < 0.05
          tags$div(class = "conclusion",
            tags$div(class = "conc-head", "Interpretation"),
            tags$p(class = "conc-text",
              "\u03c7\u00b2 = ", tags$strong(fmt4(r$val)), " on ", r$df,
              " degree(s) of freedom. ",
              "The upper-tail p-value is ", tags$strong(fmt4(r$p_upper)), ". ",
              "The result is ", tags$strong(if (sig) "statistically significant" else "not statistically significant"),
              " at the \u03b1 = 0.05 level."
            )
          )
        }
      }
    })

    outputOptions(output, "results",        suspendWhenHidden = FALSE)
    outputOptions(output, "interpretation", suspendWhenHidden = FALSE)
  })
}
