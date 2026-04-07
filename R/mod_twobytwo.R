# mod_twobytwo.R — 2×2 Table module

twobytwoUI <- function(id) {
  ns <- NS(id)
  tagList(
    h3("2\u00d72 Table Analysis"),
    p("Enter cell counts for exposed/unexposed \u00d7 cases/non-cases."),
    fluidRow(
      column(4,
        tags$table(class = "table table-bordered",
          tags$thead(tags$tr(
            tags$th(""),
            tags$th("Cases"),
            tags$th("Non-cases")
          )),
          tags$tbody(
            tags$tr(
              tags$td(tags$strong("Exposed")),
              tags$td(numericInput(ns("a"), NULL, value = NA, min = 0, width = "90px")),
              tags$td(numericInput(ns("b"), NULL, value = NA, min = 0, width = "90px"))
            ),
            tags$tr(
              tags$td(tags$strong("Unexposed")),
              tags$td(numericInput(ns("c"), NULL, value = NA, min = 0, width = "90px")),
              tags$td(numericInput(ns("d"), NULL, value = NA, min = 0, width = "90px"))
            )
          )
        )
      )
    ),
    hr(),
    uiOutput(ns("results"))
  )
}

twobytwoServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    output$results <- renderUI({
      a <- input$a; b <- input$b; c <- input$c; d <- input$d
      if (any(is.na(c(a, b, c, d)))) return(p("Enter all four cell counts.", class = "text-muted"))
      if (any(c(a, b, c, d) < 0)) return(p("All counts must be \u2265 0.", class = "text-danger"))

      n1 <- a + b  # exposed total
      n0 <- c + d  # unexposed total
      m1 <- a + c  # total cases
      m0 <- b + d  # total non-cases
      N  <- a + b + c + d

      if (n1 == 0 || n0 == 0) return(p("Row totals cannot be zero.", class = "text-danger"))

      # Point estimates
      p1 <- a / n1
      p0 <- c / n0
      RR <- p1 / p0
      RD <- p1 - p0
      OR <- (a * d) / (b * c)

      # RR CI (Taylor)
      se_ln_rr <- sqrt(b / (a * n1) + d / (c * n0))
      rr95_lo <- exp(log(RR) - 1.96 * se_ln_rr)
      rr95_hi <- exp(log(RR) + 1.96 * se_ln_rr)
      rr99_lo <- exp(log(RR) - 2.576 * se_ln_rr)
      rr99_hi <- exp(log(RR) + 2.576 * se_ln_rr)

      # RD CI (Taylor)
      se_rd <- sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0)
      rd95_lo <- RD - 1.96 * se_rd
      rd95_hi <- RD + 1.96 * se_rd
      rd99_lo <- RD - 2.576 * se_rd
      rd99_hi <- RD + 2.576 * se_rd

      # OR CI (Woolf)
      se_ln_or <- sqrt(1/a + 1/b + 1/c + 1/d)
      or95_lo <- exp(log(OR) - 1.96 * se_ln_or)
      or95_hi <- exp(log(OR) + 1.96 * se_ln_or)
      or99_lo <- exp(log(OR) - 2.576 * se_ln_or)
      or99_hi <- exp(log(OR) + 2.576 * se_ln_or)

      # Chi-square
      expected_a <- n1 * m1 / N
      chi2 <- (N * (a * d - b * c)^2) / (n1 * n0 * m1 * m0)
      chi2_cc <- (N * (abs(a * d - b * c) - N/2)^2) / (n1 * n0 * m1 * m0)
      p_chi2 <- pchisq(chi2, df = 1, lower.tail = FALSE)
      p_chi2_cc <- pchisq(chi2_cc, df = 1, lower.tail = FALSE)

      # Fisher exact p-value
      p_fisher <- fisher.test(matrix(c(a, c, b, d), 2, 2))$p.value

      # Attributable risk
      ar_exp  <- (p1 - p0) / p1          # AR among exposed
      ar_pop  <- (N / m1 * p1 - p0) / (N / m1 * p1)  # population AR

      # NNT/NNH
      nnt <- if (RD != 0) 1 / abs(RD) else Inf

      tagList(
        fluidRow(
          column(6,
            h4("Relative Risk (RR)"),
            tags$table(class = "table table-sm",
              tags$tr(tags$td("RR"), tags$td(tags$span(class="result-red", fmt4(RR)))),
              tags$tr(tags$td("95% CI"), tags$td(paste0(fmt4(rr95_lo), " \u2013 ", fmt4(rr95_hi)))),
              tags$tr(tags$td("99% CI"), tags$td(paste0(fmt4(rr99_lo), " \u2013 ", fmt4(rr99_hi))))
            ),
            h4("Risk Difference (RD)"),
            tags$table(class = "table table-sm",
              tags$tr(tags$td("RD"), tags$td(tags$span(class="result-red", fmt4(RD)))),
              tags$tr(tags$td("95% CI"), tags$td(paste0(fmt4(rd95_lo), " \u2013 ", fmt4(rd95_hi)))),
              tags$tr(tags$td("99% CI"), tags$td(paste0(fmt4(rd99_lo), " \u2013 ", fmt4(rd99_hi))))
            )
          ),
          column(6,
            h4("Odds Ratio (OR)"),
            tags$table(class = "table table-sm",
              tags$tr(tags$td("OR"), tags$td(tags$span(class="result-red", fmt4(OR)))),
              tags$tr(tags$td("95% CI (Woolf)"), tags$td(paste0(fmt4(or95_lo), " \u2013 ", fmt4(or95_hi)))),
              tags$tr(tags$td("99% CI"), tags$td(paste0(fmt4(or99_lo), " \u2013 ", fmt4(or99_hi))))
            ),
            h4("Tests"),
            tags$table(class = "table table-sm",
              tags$tr(tags$td("Chi-square"), tags$td(fmt4(chi2)), tags$td(paste0("p=", fmt4(p_chi2)))),
              tags$tr(tags$td("Chi-sq (cc)"), tags$td(fmt4(chi2_cc)), tags$td(paste0("p=", fmt4(p_chi2_cc)))),
              tags$tr(tags$td("Fisher exact p"), tags$td(fmt4(p_fisher)))
            )
          )
        ),
        fluidRow(
          column(6,
            h4("Attributable Risk"),
            tags$table(class = "table table-sm",
              tags$tr(tags$td("AR (exposed)"),    tags$td(tags$span(class="result-orange", fmt4(ar_exp)))),
              tags$tr(tags$td("AR (population)"), tags$td(tags$span(class="result-orange", fmt4(ar_pop))))
            )
          ),
          column(6,
            h4("NNT / NNH"),
            tags$table(class = "table table-sm",
              tags$tr(tags$td(if (RD > 0) "NNH" else "NNT"),
                      tags$td(tags$span(class="result-orange", fmt2(nnt))))
            )
          )
        )
      )
    })
  })
}
