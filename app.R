library(shiny)
library(shinyjs)

# Source all modules and helpers
source("R/helpers.R")
source("R/mod_twobytwo.R")
source("R/mod_quickcalc.R")
source("R/mod_rate_data.R")
source("R/mod_risk_data.R")
source("R/mod_cc_data.R")
source("R/mod_exact_or.R")
source("R/mod_cohort_power.R")
source("R/mod_cc_power.R")
source("R/mod_study_size.R")
source("R/mod_standardization.R")
source("R/mod_meta_analysis.R")
source("R/mod_pvalue_functions.R")
source("R/mod_matched_cc.R")
source("R/mod_life_table.R")
source("R/mod_seasonal.R")

# ── Calculator registry ──────────────────────────────────────────────────────
calculators <- list(
  list(id = "twobytwo",     name = "2\u00d72 Table",            desc = "OR, RR, RD, exact tests, AR, NNT",           icon = "\U0001f4ca", cat = "strat", ui = twobytwoUI,     server = twobytwoServer),
  list(id = "rate_data",    name = "Rate Data",              desc = "MH rate ratio, 3-stratum stratified analysis", icon = "\u23f1",      cat = "strat", ui = rateDataUI,     server = rateDataServer),
  list(id = "risk_data",    name = "Risk Data",              desc = "MH risk ratio and risk difference",            icon = "\U0001f9ea", cat = "strat", ui = riskDataUI,     server = riskDataServer),
  list(id = "cc_data",      name = "Case-control Data",      desc = "MH odds ratio, stratified 2\u00d72",           icon = "\U0001f50d", cat = "strat", ui = ccDataUI,       server = ccDataServer),
  list(id = "exact_or",     name = "O.R. Exact",             desc = "Exact OR via hypergeometric CI",               icon = "\U0001f3af", cat = "strat", ui = exactOrUI,      server = exactOrServer),
  list(id = "matched_cc",   name = "Matched Case-control",   desc = "McNemar & generalised matched analysis",       icon = "\U0001f517", cat = "strat", ui = matchedCcUI,    server = matchedCcServer),
  list(id = "standardization", name = "Standardization",    desc = "Direct & indirect standardization (DSR/SMR)",  icon = "\u2696\ufe0f", cat = "strat", ui = standardizationUI, server = standardizationServer),
  list(id = "cohort_power", name = "Cohort Power",           desc = "Kelsey two-proportion power table",            icon = "\u26a1",     cat = "power", ui = cohortPowerUI,  server = cohortPowerServer),
  list(id = "cc_power",     name = "Case-control Power",     desc = "Schlesselman OR power table",                  icon = "\u26a1",     cat = "power", ui = ccPowerUI,      server = ccPowerServer),
  list(id = "study_size",   name = "Study Size",             desc = "Precision-based sample size (CI width)",       icon = "\U0001f4cf", cat = "power", ui = studySizeUI,    server = studySizeServer),
  list(id = "quickcalc",   name = "Quickcalc",               desc = "Proportion, rate, SMR, z, \u03c7\u00b2",       icon = "\u2699\ufe0f", cat = "tools", ui = quickcalcUI,  server = quickcalcServer),
  list(id = "meta_analysis", name = "Meta-Analysis",         desc = "Fixed/random effects forest plot",             icon = "\U0001f333", cat = "tools", ui = metaAnalysisUI, server = metaAnalysisServer),
  list(id = "pvalue_functions", name = "P-value Functions",  desc = "Two-curve p-value function overlay",           icon = "\U0001f4c8", cat = "tools", ui = pvalueFunctionsUI, server = pvalueFunctionsServer),
  list(id = "life_table",   name = "Life Table",             desc = "Actuarial survival curve with CI",             icon = "\U0001f4c5", cat = "tools", ui = lifeTableUI,    server = lifeTableServer),
  list(id = "seasonal",     name = "Seasonal Analysis",      desc = "Walter-Elwood periodicity test",               icon = "\U0001f321\ufe0f", cat = "tools", ui = seasonalUI, server = seasonalServer)
)

calc_ids <- sapply(calculators, `[[`, "id")

# ── Glossary entries ─────────────────────────────────────────────────────────
glossary <- list(
  list(term = "Odds Ratio (OR)",         def = "Ratio of the odds of exposure among cases to the odds among controls. Approximates the RR when disease is rare."),
  list(term = "Risk Ratio (RR)",         def = "Ratio of the probability of disease in the exposed group to the probability in the unexposed group."),
  list(term = "Risk Difference (RD)",    def = "Absolute difference in disease probability between exposed and unexposed. Also called attributable risk."),
  list(term = "95% Confidence Interval", def = "A range of values within which the true population parameter is expected to fall 95% of the time under repeated sampling."),
  list(term = "Mantel-Haenszel (MH)",    def = "A method for pooling stratum-specific estimates, weighting each stratum by its information content."),
  list(term = "Mid-P",                   def = "A modified exact p-value that removes half the probability of the observed value; less conservative than Fisher exact."),
  list(term = "Person-years",            def = "Total time at risk summed across all individuals; used as the denominator for incidence rates."),
  list(term = "SMR",                     def = "Standardised Mortality (or Morbidity) Ratio: ratio of observed to expected events, indirect standardization."),
  list(term = "Cochran Q",               def = "Statistic testing heterogeneity across studies in meta-analysis. Large Q (small p) suggests the studies are not estimating the same effect."),
  list(term = "DerSimonian-Laird",       def = "Method-of-moments estimator for the between-study variance (tau\u00b2) in random-effects meta-analysis."),
  list(term = "NNT / NNH",              def = "Number Needed to Treat (or Harm): reciprocal of the absolute risk difference. How many people must be exposed to produce one additional case."),
  list(term = "Attributable Risk",       def = "The proportion of disease in exposed individuals (or the population) attributable to the exposure, assuming causality.")
)

# ── UI ───────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  useShinyjs(),
  tags$head(
    tags$link(rel = "stylesheet", href = "styles.css"),
    tags$script(src = "tooltip.js"),
    tags$title("EpiSheet")
  ),

  # Top nav
  tags$div(class = "epi-nav",
    tags$span(class = "epi-nav-brand",
      tags$span(style = "font-size:18px; line-height:1;", "\U0001f9ec"),
      "EpiSheet"
    ),
    tags$span(class = "epi-nav-ref", "Based on Episheet by K. Rothman"),
    tags$button(class = "epi-nav-glossary", onclick = "Shiny.setInputValue('show_glossary', Math.random())",
      "\U0001f4d6 Glossary")
  ),

  # Main content area
  tags$div(class = "epi-app",
    uiOutput("main_content")
  ),

  # Glossary overlay
  uiOutput("glossary_panel")
)

# ── Server ───────────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  # State: which calculator is active (NULL = landing page)
  active_calc <- reactiveVal(NULL)
  show_glossary <- reactiveVal(FALSE)

  # Navigate to calculator
  observe({
    req(input$nav_to)
    active_calc(input$nav_to)
  })

  # Back to landing
  observe({
    req(input$nav_back)
    active_calc(NULL)
  })

  # Glossary toggle
  observe({
    req(input$show_glossary)
    show_glossary(!isolate(show_glossary()))
  })
  observe({
    req(input$close_glossary)
    show_glossary(FALSE)
  })

  # ── Landing page ────────────────────────────────────────────────────────────
  output$main_content <- renderUI({
    calc <- active_calc()
    if (is.null(calc)) {
      landing_ui()
    } else {
      calc_def <- calculators[[which(calc_ids == calc)]]
      tagList(
        tags$button(class = "calc-back",
          onclick = "Shiny.setInputValue('nav_back', Math.random())",
          "\u2190 All calculators"),
        tags$div(class = "calc-header",
          tags$h1(class = "calc-title", calc_def$name),
          tags$p(class = "calc-subtitle", calc_def$desc)
        ),
        calc_def$ui(calc)
      )
    }
  })

  # ── Glossary panel ──────────────────────────────────────────────────────────
  output$glossary_panel <- renderUI({
    if (!show_glossary()) return(NULL)
    tags$div(class = "glossary-overlay",
      tags$button(class = "glossary-close",
        onclick = "Shiny.setInputValue('close_glossary', Math.random())", "\u00d7"),
      tags$div(class = "glossary-title", "\U0001f4d6 Glossary"),
      tagList(lapply(glossary, function(g) {
        tagList(
          tags$div(class = "glossary-term", g$term),
          tags$p(class = "glossary-def", g$def)
        )
      }))
    )
  })

  # ── Initialise all calculator servers ───────────────────────────────────────
  for (calc_def in calculators) {
    local({
      cd <- calc_def
      cd$server(cd$id)
    })
  }
}

# ── Landing page UI helper ───────────────────────────────────────────────────
landing_ui <- function() {
  strat_calcs  <- Filter(function(x) x$cat == "strat",  calculators)
  power_calcs  <- Filter(function(x) x$cat == "power",  calculators)
  tools_calcs  <- Filter(function(x) x$cat == "tools",  calculators)

  make_card <- function(calc, icon_class) {
    tags$div(class = "lcard",
      onclick = sprintf("Shiny.setInputValue('nav_to', '%s', {priority: 'event'})", calc$id),
      tags$div(class = paste("lcard-icon", icon_class), calc$icon),
      tags$div(class = "lcard-name", calc$name),
      tags$div(class = "lcard-desc", calc$desc)
    )
  }

  tagList(
    tags$div(class = "landing-header",
      tags$h1(class = "landing-title", "EpiSheet"),
      tags$p(class = "landing-sub",
        "Epidemiologic analysis tools based on Rothman's Episheet. ",
        "Select a calculator to begin.")
    ),

    tags$div(class = "group-label", "Stratified Analysis"),
    tags$div(class = "card-grid",
      lapply(strat_calcs, make_card, icon_class = "ic-teal")
    ),

    tags$div(class = "group-label", "Sample Size & Power"),
    tags$div(class = "card-grid",
      lapply(power_calcs, make_card, icon_class = "ic-blue")
    ),

    tags$div(class = "group-label", "Tools & Advanced"),
    tags$div(class = "card-grid",
      lapply(tools_calcs, make_card, icon_class = "ic-purple")
    )
  )
}

shinyApp(ui, server)
