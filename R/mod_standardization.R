# mod_standardization.R — Direct Standardization module
# Method: DSR = sum(w_i * rate_i) / sum(w_i)
# Var(DSR) = sum(w_i^2 * rate_i / den_i) / sum(w_i)^2
# Normal-approximation 95%/99% CI

standardizationUI <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(class = "calc-layout",

      # ── Input panel ──────────────────────────────────────────────────────────
      tags$div(class = "panel",
        tags$div(class = "panel-head ph-teal", "Input"),
        tags$div(class = "panel-body",

          # Mode toggle and group count
          tags$div(
            style = "display: grid; grid-template-columns: 1fr 1fr; gap: 10px; margin-bottom: 14px;",
            tags$div(
              tags$label(
                style = "font-size: 11px; font-weight: 500; color: #4a6070; display: block; margin-bottom: 4px;",
                "Analysis type"
              ),
              radioButtons(ns("mode"), label = NULL,
                choices  = c("Rates" = "R", "Proportions" = "P"),
                selected = "R",
                inline   = TRUE)
            ),
            numericInput(ns("n_groups"), "Comparison groups",
              value = 2, min = 1, max = 5, step = 1)
          ),

          # Column header hint
          tags$p(
            style = "font-size: 11px; color: #4a6070; margin-bottom: 8px;",
            "Enter numerator, denominator, and weight for each stratum. ",
            "Rows with a positive weight and at least one valid group are included."
          ),

          # 10 stratum rows (always rendered; filled as needed)
          tags$div(
            id = ns("strata_container"),
            lapply(seq_len(10), function(i) {
              tags$div(
                class = "stratum-row",
                id    = ns(paste0("row_", i)),
                style = "margin-bottom: 8px; padding: 8px; border: 0.5px solid #d0dde3; border-radius: 8px;",
                # Stratum label + weight row
                tags$div(
                  style = "display: grid; grid-template-columns: 1.4fr 0.6fr; gap: 6px; margin-bottom: 6px;",
                  tags$div(
                    tags$label(
                      style = "font-size: 10px; color: #4a6070;",
                      paste0("Stratum ", i, " label")
                    ),
                    textInput(ns(paste0("label_", i)), label = NULL,
                      placeholder = paste0("Stratum ", i))
                  ),
                  numericInput(ns(paste0("weight_", i)),
                    label = paste0("Weight ", i),
                    value = 1, min = 0)
                ),
                # Per-group numerator / denominator inputs (rendered dynamically)
                uiOutput(ns(paste0("group_inputs_", i)))
              )
            })
          ),

          tags$hr(style = "border-color: #d0dde3; margin: 12px 0;"),

          # Reference group inputs (for SMR note)
          tags$p(
            style = "font-size: 11px; font-weight: 500; color: #4a6070; margin-bottom: 6px;",
            "Reference group (optional — used for indirect standardization / SMR)"
          ),
          tags$div(
            style = "display: grid; grid-template-columns: 1fr 1fr; gap: 8px;",
            numericInput(ns("ref_num"), "Reference numerator",   value = NULL, min = 0),
            numericInput(ns("ref_den"), "Reference denominator", value = NULL, min = 0)
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

    # ── Interpretation block (full-width below columns) ───────────────────────
    uiOutput(ns("interpretation"))
  )
}

standardizationServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ── NULL-coalescing helper ─────────────────────────────────────────────────
    `%||%` <- function(a, b) if (!is.null(a)) a else b

    # ── Metric card helper ─────────────────────────────────────────────────────
    metric <- function(label, value, sub = NULL, cls = "m-blue") {
      tags$div(class = paste("metric", cls),
        tags$div(class = "metric-label", label),
        tags$div(class = "metric-value", value),
        if (!is.null(sub)) tags$div(class = "metric-sub", sub)
      )
    }

    # ── Dynamic per-group numerator/denominator inputs for each stratum row ───
    lapply(seq_len(10), function(i) {
      output[[paste0("group_inputs_", i)]] <- renderUI({
        ng <- max(1L, min(5L, as.integer(input$n_groups %||% 2)))
        col_def <- paste0("repeat(", ng * 2L, ", 1fr)")
        tags$div(
          style = paste0("display: grid; grid-template-columns: ", col_def, "; gap: 6px;"),
          lapply(seq_len(ng), function(g) {
            tagList(
              numericInput(
                ns(paste0("num_", i, "_", g)),
                label = paste0("G", g, " num"),
                value = NULL, min = 0
              ),
              numericInput(
                ns(paste0("den_", i, "_", g)),
                label = paste0("G", g, " den"),
                value = NULL, min = 0
              )
            )
          })
        )
      })
    })

    # ── Core computation (triggered by Calculate button) ──────────────────────
    computed <- eventReactive(input$calc, {
      ng    <- max(1L, min(5L, as.integer(input$n_groups %||% 2)))
      mode  <- input$mode %||% "R"

      # Collect stratum rows — only keep rows where weight > 0 AND group 1 has
      # a valid (non-NA, positive) denominator
      rows <- lapply(seq_len(10), function(i) {
        w <- input[[paste0("weight_", i)]]
        if (is.null(w) || is.na(w) || w <= 0) return(NULL)
        nums <- sapply(seq_len(ng), function(g) {
          v <- input[[paste0("num_", i, "_", g)]]
          if (is.null(v)) NA_real_ else v
        })
        dens <- sapply(seq_len(ng), function(g) {
          v <- input[[paste0("den_", i, "_", g)]]
          if (is.null(v)) NA_real_ else v
        })
        # Row is complete if group 1 has both num & den
        if (is.na(nums[1]) || is.na(dens[1]) || dens[1] <= 0) return(NULL)
        list(w = w, nums = nums, dens = dens)
      })
      rows <- Filter(Negate(is.null), rows)

      validate(
        need(length(rows) >= 1,
             "Enter at least one complete stratum row (weight + numerator + denominator for Group 1).")
      )

      weights <- sapply(rows, `[[`, "w")

      # Per-group DSR, Var(DSR), 95% CI, 99% CI
      group_results <- lapply(seq_len(ng), function(g) {
        nums <- sapply(rows, function(r) r$nums[g])
        dens <- sapply(rows, function(r) r$dens[g])

        valid <- !is.na(nums) & !is.na(dens) & dens > 0
        if (sum(valid) == 0) {
          return(list(dsr = NA_real_, ci95 = c(NA_real_, NA_real_),
                      ci99 = c(NA_real_, NA_real_), se = NA_real_,
                      n_valid = 0L))
        }

        w_v   <- weights[valid]
        num_v <- nums[valid]
        den_v <- dens[valid]
        rate_v <- num_v / den_v

        sum_w   <- sum(w_v)
        dsr     <- sum(w_v * rate_v) / sum_w
        var_dsr <- sum(w_v^2 * rate_v / den_v) / sum_w^2
        se      <- sqrt(var_dsr)

        ci95 <- dsr + c(-1, 1) * 1.96  * se
        ci99 <- dsr + c(-1, 1) * 2.576 * se

        list(dsr = dsr, ci95 = ci95, ci99 = ci99, se = se, n_valid = sum(valid))
      })

      # Optional SMR for group 1 if reference group provided
      ref_num <- input$ref_num %||% NA_real_
      ref_den <- input$ref_den %||% NA_real_
      smr     <- NA_real_
      if (!is.null(ref_num) && !is.null(ref_den) &&
          !is.na(ref_num) && !is.na(ref_den) && ref_den > 0) {
        ref_rate <- ref_num / ref_den
        # Expected events for group 1 using reference rates
        dens_g1 <- sapply(rows, function(r) r$dens[1])
        expected <- sum(ref_rate * dens_g1, na.rm = TRUE)
        obs_g1   <- sum(sapply(rows, function(r) r$nums[1]), na.rm = TRUE)
        if (expected > 0) smr <- obs_g1 / expected
      }

      list(
        group_results = group_results,
        ng            = ng,
        mode          = mode,
        n_strata      = length(rows),
        smr           = smr
      )
    })

    # ── Results output ─────────────────────────────────────────────────────────
    output$results <- renderUI({
      r   <- computed()
      ng  <- r$ng
      type_label <- if (r$mode == "R") "Rate" else "Proportion"

      # Card colour cycling: first group gets m-primary wide, rest cycle through
      card_classes <- c("m-primary wide", "m-blue", "m-teal", "m-purple", "m-blue")

      group_cards <- lapply(seq_len(ng), function(g) {
        gr  <- r$group_results[[g]]
        cls <- card_classes[min(g, length(card_classes))]
        if (g > 1) cls <- gsub(" wide", "", cls, fixed = TRUE)

        if (is.na(gr$dsr)) {
          sub_txt <- "Insufficient data for this group"
          val_txt <- "\u2014"
        } else {
          sub_txt <- paste0(
            "95% CI: ", fmt4(gr$ci95[1]), "\u2013", fmt4(gr$ci95[2]),
            " \u2502 99% CI: ", fmt4(gr$ci99[1]), "\u2013", fmt4(gr$ci99[2])
          )
          val_txt <- fmt4(gr$dsr)
        }

        metric(
          label = paste0("Group ", g, " — Adjusted ", type_label),
          value = val_txt,
          sub   = sub_txt,
          cls   = cls
        )
      })

      smr_card <- if (!is.na(r$smr)) {
        metric(
          label = "SMR (Group 1 / Reference)",
          value = fmt4(r$smr),
          sub   = "Indirect standardization",
          cls   = "m-teal"
        )
      }

      tagList(
        tags$div(
          style = "font-size: 10px; color: #4a6070; margin-bottom: 8px;",
          paste0(
            r$n_strata, " str", if (r$n_strata == 1) "atum" else "ata",
            " \u2502 Direct standardization (DSR) \u2502 Normal-approximation CI"
          )
        ),
        tags$div(class = "metrics",
          do.call(tagList, group_cards),
          smr_card
        )
      )
    })

    # ── Interpretation ─────────────────────────────────────────────────────────
    output$interpretation <- renderUI({
      r     <- computed()
      g1    <- r$group_results[[1]]
      mode_word <- if (r$mode == "R") "rate" else "proportion"

      # Build comparison sentence if 2+ groups both valid
      comparison <- if (r$ng >= 2) {
        g2 <- r$group_results[[2]]
        if (!is.na(g1$dsr) && !is.na(g2$dsr)) {
          ratio  <- g1$dsr / g2$dsr
          dir    <- if (ratio > 1) "higher" else if (ratio < 1) "lower" else "equal"
          tagList(
            " The adjusted ", mode_word, " for Group 1 (",
            tags$strong(fmt4(g1$dsr)), ") is ",
            tags$strong(dir),
            " than that for Group 2 (",
            tags$strong(fmt4(g2$dsr)), ")."
          )
        }
      }

      tags$div(class = "conclusion",
        tags$div(class = "conc-head", "Interpretation"),
        tags$p(class = "conc-text",
          "Direct standardization was applied across ",
          tags$strong(r$n_strata),
          " str", if (r$n_strata == 1) "atum" else "ata",
          " using the supplied stratum weights. ",
          if (!is.na(g1$dsr)) {
            tagList(
              "The adjusted ", mode_word, " for Group 1 is ",
              tags$strong(fmt4(g1$dsr)),
              " (95% CI: ", fmt4(g1$ci95[1]), "\u2013", fmt4(g1$ci95[2]),
              "; 99% CI: ", fmt4(g1$ci99[1]), "\u2013", fmt4(g1$ci99[2]), ")."
            )
          },
          comparison,
          if (!is.na(r$smr)) {
            tagList(
              " The standardized mortality/morbidity ratio (SMR) for Group 1 relative to the reference is ",
              tags$strong(fmt4(r$smr)), "."
            )
          }
        )
      )
    })

    outputOptions(output, "results",        suspendWhenHidden = FALSE)
    outputOptions(output, "interpretation", suspendWhenHidden = FALSE)
  })
}
