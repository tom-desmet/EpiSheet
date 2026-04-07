# mod_cc_power.R — Case-Control Power module
# Schlesselman power formula (port of VBA Module3: ccPowercalc)
# Structure mirrors mod_cohort_power.R

ccPowerUI <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(class = "calc-layout",

      # ── Input panel ────────────────────────────────────────────────────────────
      tags$div(class = "panel",
        tags$div(class = "panel-head ph-teal", "Input"),
        tags$div(class = "panel-body",

          tags$div(class = "input-group",
            numericInput(ns("z_alpha"), "Z-alpha (one-sided \u03b1)",
                         value = 1.645, min = 0.5, max = 4, step = 0.001)
          ),
          tags$div(class = "input-group",
            numericInput(ns("r"), "Controls per case (r)",
                         value = 1, min = 0.01, max = 100, step = 0.5)
          ),
          tags$div(class = "input-group",
            numericInput(ns("p0"), "Exposure prevalence in controls (p\u2080)",
                         value = 0.3, min = 1e-6, max = 0.999, step = 0.01)
          ),
          tags$div(class = "input-group",
            numericInput(ns("max_n"), "Maximum N cases",
                         value = 1000, min = 10, max = 1e6, step = 50)
          ),

          tags$hr(),
          tags$p(class = "panel-hint", "OR effect levels (6 curves):"),

          tags$div(class = "input-row",
            tags$div(class = "input-group",
              numericInput(ns("or1"), "OR\u2081", value = 1.5, min = 0.01, max = 100, step = 0.1)
            ),
            tags$div(class = "input-group",
              numericInput(ns("or2"), "OR\u2082", value = 2,   min = 0.01, max = 100, step = 0.1)
            ),
            tags$div(class = "input-group",
              numericInput(ns("or3"), "OR\u2083", value = 3,   min = 0.01, max = 100, step = 0.1)
            )
          ),
          tags$div(class = "input-row",
            tags$div(class = "input-group",
              numericInput(ns("or4"), "OR\u2084", value = 5,   min = 0.01, max = 100, step = 0.1)
            ),
            tags$div(class = "input-group",
              numericInput(ns("or5"), "OR\u2085", value = 7,   min = 0.01, max = 100, step = 0.1)
            ),
            tags$div(class = "input-group",
              numericInput(ns("or6"), "OR\u2086", value = 10,  min = 0.01, max = 100, step = 0.1)
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
    )
  )
}

ccPowerServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    computed <- eventReactive(input$calc, {
      validate(
        need(!is.na(input$z_alpha) && input$z_alpha > 0,
             "Z-alpha must be a positive number."),
        need(!is.na(input$r) && input$r > 0,
             "Controls per case (r) must be positive."),
        need(!is.na(input$p0) && input$p0 > 0 && input$p0 < 1,
             "p\u2080 must be between 0 and 1."),
        need(!is.na(input$max_n) && input$max_n >= 10,
             "Maximum N must be at least 10.")
      )

      z     <- input$z_alpha
      r     <- input$r
      p0    <- input$p0
      max_n <- input$max_n

      or_levels <- c(input$or1, input$or2, input$or3,
                     input$or4, input$or5, input$or6)
      or_levels <- or_levels[!is.na(or_levels) & or_levels > 0]
      validate(need(length(or_levels) > 0, "Enter at least one OR level."))

      # N sequence: 50 evenly-spaced steps from 10 to max_n
      n_seq <- unique(round(seq(10, max_n, length.out = 50)))

      # Power matrix
      pow_mat <- sapply(or_levels, function(OR) {
        sapply(n_seq, function(N) cc_powercalc(z, r, p0, N, OR))
      })
      if (!is.matrix(pow_mat)) pow_mat <- matrix(pow_mat, nrow = length(n_seq))
      colnames(pow_mat) <- paste0("OR=", or_levels)

      # Long-format for ggplot
      df_long <- do.call(rbind, lapply(seq_along(or_levels), function(j) {
        data.frame(
          N     = n_seq,
          Power = pow_mat[, j],
          Level = factor(paste0("OR = ", or_levels[j]),
                         levels = paste0("OR = ", or_levels))
        )
      }))

      # Summary table
      key_idx  <- unique(round(seq(1, length(n_seq), length.out = 8)))
      key_ns   <- n_seq[key_idx]
      tbl_data <- data.frame(N_cases = key_ns)
      for (j in seq_along(or_levels)) {
        col <- sprintf("OR=%s", or_levels[j])
        tbl_data[[col]] <- round(pow_mat[key_idx, j], 4)
      }

      list(df_long = df_long, tbl_data = tbl_data,
           or_levels = or_levels, max_n = max_n)
    })

    output$results <- renderUI({
      tagList(
        plotOutput(ns("power_chart"), height = "380px"),
        tags$h4("Power at selected sample sizes"),
        tableOutput(ns("power_table"))
      )
    })

    output$power_chart <- renderPlot({
      r <- computed()
      ggplot2::ggplot(r$df_long,
                      ggplot2::aes(x = N, y = Power, colour = Level)) +
        ggplot2::geom_line(linewidth = 1.1) +
        ggplot2::geom_hline(yintercept = 0.80, linetype = "dashed",
                            colour = "grey50", linewidth = 0.6) +
        ggplot2::scale_colour_brewer(palette = "Dark2", name = "Effect level") +
        ggplot2::scale_x_continuous(labels = scales::comma) +
        ggplot2::scale_y_continuous(limits = c(0, 1),
                                    breaks = seq(0, 1, 0.2),
                                    labels = scales::percent) +
        ggplot2::labs(
          title = "Case-control study power",
          x     = "Number of cases",
          y     = "Power"
        ) +
        ggplot2::theme_bw(base_size = 13) +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank(),
          legend.position  = "right"
        )
    })

    output$power_table <- renderTable({
      computed()$tbl_data
    }, digits = 4, striped = TRUE, hover = TRUE, bordered = TRUE)

    outputOptions(output, "results", suspendWhenHidden = FALSE)
  })
}
