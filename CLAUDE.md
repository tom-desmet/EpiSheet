# CLAUDE.md — Episheet Shiny Rebuild

## Project overview

Rebuild Kenneth Rothman's **Episheet** (epidemiologic spreadsheet tool) as a
Shiny app in R. Episheet is a protected XLS workbook with 13 analytical tabs
covering core epidemiologic methods (stratified analysis, meta-analysis, power,
life tables, seasonal analysis, etc.). This file is the authoritative spec.
Do not ask the user to re-explain any of this context.

---

## Repository

```
https://github.com/tom-desmet/EpiSheet
git clone https://github.com/tom-desmet/EpiSheet.git
cd EpiSheet
```

Branch strategy: work on `main` for now; create feature branches per tab
(e.g. `feat/quickcalc`, `feat/rate-data`) once the scaffold is in place.

## R package bootstrap

```r
# Install the partial R port that already exists — use it where applicable
install.packages("episheet")   # James Black's CRAN package (Rate Data, Risk Data, pvalueplot already done)

# Core dependencies
install.packages(c(
  "shiny", "bslib", "shinyjs",
  "epiR", "epitools", "exact2x2",
  "meta", "metafor",
  "survival",
  "circular",
  "pwr",
  "PHEindicatormethods",
  "ggplot2", "plotly",
  "dplyr", "tidyr"
))
```

---

## UI conventions (critical — match the Excel original)

| Excel cell colour | Meaning           | Shiny equivalent                                      |
|-------------------|-------------------|-------------------------------------------------------|
| **Yellow**        | User input        | `numericInput` / `textInput` / `selectInput`          |
| **Orange**        | Locked result     | `textOutput` / `verbatimTextOutput` styled orange     |
| **Red text**      | Key result value  | Output rendered in `color: #CC0000; font-weight: bold`|
| **Red border box**| Crude data section| Separate `wellPanel` with red border                  |

**No submit button.** All outputs update reactively on every input change,
exactly as in Excel.

CSS to include in `www/styles.css`:

```css
/* Yellow input cells */
.form-control, .selectize-input {
  background-color: #FFFACD !important;
}

/* Orange result cells */
.result-orange {
  background-color: #FFB347;
  color: #000;
  font-weight: bold;
  padding: 4px 8px;
  border-radius: 3px;
  display: inline-block;
  min-width: 80px;
  text-align: right;
}

/* Red result text */
.result-red {
  color: #CC0000;
  font-weight: bold;
}

/* Crude data box */
.crude-box {
  border: 2px solid #CC0000;
  padding: 10px;
  margin-top: 10px;
}
```

---

## App structure

```
episheet_shiny/
├── CLAUDE.md               ← this file
├── app.R                   ← or ui.R + server.R
├── www/
│   └── styles.css
├── R/
│   ├── mod_rate_data.R
│   ├── mod_risk_data.R
│   ├── mod_cc_data.R
│   ├── mod_standardization.R
│   ├── mod_meta_analysis.R
│   ├── mod_pvalue_functions.R
│   ├── mod_exact_or.R
│   ├── mod_matched_cc.R
│   ├── mod_life_table.R
│   ├── mod_seasonal.R
│   ├── mod_cohort_power.R
│   ├── mod_cc_power.R
│   ├── mod_study_size.R
│   ├── mod_quickcalc.R
│   └── helpers.R           ← ported VBA functions
└── tests/
    └── testthat/
```

Use **Shiny modules** for every tab. Top-level UI is `navbarPage()`.

---

## Tab inventory and computational spec

### 1. Rate Data
- **Inputs (yellow):** Up to 3 strata × (Cases + Person-time) for Exposed and
  Unexposed; Crude data section (red border box) at bottom
- **Outputs (orange/red):** MH Rate Ratio, 90%/95%/99% CI, Variance of pooled
  estimate, P-value for homogeneity
- **Method:** Mantel-Haenszel rate ratio; Taylor-series variance; χ² homogeneity
- **R:** Adapt `episheet::rate_data()` — already on CRAN
- **Reference:** Modern Epidemiology 3rd ed., Ch. 15

### 2. Risk Data
- **Inputs:** Up to 3 strata × (Cases + Non-cases + Total) for Exposed/Unexposed
  + crude box. Tooltip: "Number of Exposed People at Risk"
- **Outputs:** MH Risk Ratio, Risk Difference, 90%/95%/99% CI, homogeneity P
- **Method:** MH-RR, pooled RD, Taylor CI
- **R:** Adapt `episheet::risk_data()` — already on CRAN

### 3. Case-control Data
- **Inputs:** Up to 3 strata × 2×2 (Cases/Controls × Exposed/Unexposed) + crude box
- **Outputs:** OR=, 95% CI (Woolf), 99% CI, crude OR, P-value for homogeneity
- **Example data visible in screenshots:** Stratum 1: 104/1059/12/1163
- **Method:** MH-OR, Woolf CI, Gart CI
- **R:** `epiR::epi.2by2(method="cohort.count")`, custom MH weights

### 4. Standardization
- **Inputs:** Toggle button: **R** (rates) or **P** (proportions); up to 5
  comparison groups + 1 reference group; up to 48 strata rows (yellow); each
  row = stratum label + numerator + denominator + weight per group
- **Outputs:** Adjusted Estimates, 95%/99% lower/upper limits per group
- **Method:** Direct standardization (DSR) and indirect (ISR/SMR); Dobson-Kuulasmaa
  Poisson CI
- **R:** `PHEindicatormethods::phe_dsr()`, `phe_smr()`
- **Reference:** Rothman KJ and Greenland S, Modern Epidemiology 3rd ed., Ch. 15

### 5. Cohort Power
- **Inputs (yellow):** Z-alpha (default 1.645), Unexposed Exposure Rate,
  Risk in Unexposed, Max No. Exposed; RR Effect Levels (4 columns: e.g. 1.5, 2, 3, 5, 7, 10)
- **Outputs:** Power table (N Exposed × RR → power); embedded line chart
  (x = N Exposed linear, y = Power 0–1, one line per RR level)
- **Method:** Kelsey/Fleiss two-proportion formula (port VBA `Powercalc` from Module4):
  ```r
  powercalc <- function(z_alpha, p0, r, N, RR) {
    p1 <- RR * p0
    if (p1 * p0 == 0) return(NA)
    a <- p1 * (1 - p0) + p0 * (1 - p1)
    b <- (r - 1) * p0 * (1 - p0)
    zbnum <- -z_alpha * (a + b) + abs(p1 - p0) * sqrt(r * N * (a + b))
    zbden <- sqrt((a + b) * (r * a - b) - r * (p1 - p0)^2)
    pnorm(zbnum / zbden)
  }
  ```
- **R:** Use custom port above; `pwr::pwr.2p.test()` as cross-check

### 6. Case-control Power
- **Inputs:** Z-alpha (1.645), Control/Case Ratio, Exposure Prevalence (p0),
  Maximum No. Cases; OR Effect Levels (same 4-column structure)
- **Outputs:** Power table (N Cases × OR → power); line chart same format
- **Method:** Schlesselman formula (port VBA `ccPowercalc` from Module3):
  ```r
  cc_powercalc <- function(z_alpha, r, p1, p0, N, OR) {
    a <- p1 * (1 - p0) + p0 * (1 - p1)   # note: p1 derived from p0 and OR
    b <- (r - 1) * p0 * (1 - p0)
    zbnum <- -z_alpha * (a + b) + abs(p1 - p0) * sqrt(r * N * (a + b))
    zbden <- sqrt((a + b) * (r * a - b) - r * (p1 - p0)^2)
    pnorm(zbnum / zbden)
  }
  # p1 = (OR * p0) / (1 + p0 * (OR - 1))
  ```

### 7. Study Size (precision-based — NOT power-based, unique tab)
- **Inputs:** Type of study (cohort/case-control dropdown), Confidence Level (%),
  OR/RR Upper Bound <, Expected OR/RR >, Exposure Prevalence, Allocation Ratio,
  Maximum Size
- **Outputs:** Table of (Total Size, N1, N2, SE(ln OR), Span, Pr) for
  incrementing sample sizes; chart: x = Total Size (0–max), y = Pr Upper CI
  ≤ bound (0–1), single red curve
- **Method:** Works backward from desired CI width to required N using SE of
  ln(OR) = sqrt(1/a + 1/b + 1/c + 1/d)
- **Reference:** Rothman KJ, Modern Epidemiology

### 8. Quickcalc
- **Three sub-calculators side by side on one panel:**
  1. **Proportion:** Numerator / Denominator / Confidence Level (%) →
     Point Estimate + Exact CI (Mid-P and Fisher) + Approximate CI + Exact 2-tail P
  2. **Rate:** Numerator / Denominator / Confidence Level →
     same outputs (Poisson-based)
  3. **SMR:** Observed / Expected / Confidence Level →
     same outputs (Poisson-based)
  4. **Std Normal Deviate:** Value + DF → P-value
  5. **Chi-Square:** Value + DF → P-value
- **Method:** Port VBA `Binomialcalc` (Module5) and `Poissoncalc` (Module6).
  Both use iterative likelihood search (mid-p and Fisher exact variants).
  Approximate CIs use Wilson (binomial) and Byar (Poisson).
- **R:** `exact2x2::exact2x2()`, `stats::binom.test()`, `stats::poisson.test()`
  as cross-checks; port the iterative VBA for the mid-p variant

### 9. Seasonal Analysis
- **Inputs:** 12 rows: Month (Jan–Dec), Cases, Total Population (optional)
- **Outputs:** Peak/Low Ratio, CIs at 90% and 95%, Time of Peak (date + degrees),
  Hemi-amplitude, Adjusted monthly mean; intermediate trig columns (cos/sin/delta/theta)
- **Chart:** Line chart with "Adjusted Monthly Frequency" (blue) and
  "Fitted Frequency" (red sinusoidal curve); x = months Jan–Dec
- **Method:** Walter-Elwood test for periodicity (vector sum on unit circle):
  ```r
  # r = sqrt(sum(n*cos(theta))^2 + sum(n*sin(theta))^2) / sum(n)
  # theta_k = 2*pi*(k-0.5)/12 for month k
  # Rayleigh test: z = N * r^2
  ```
- **R:** `circular::rayleigh.test()`, `ggplot2` for plot
- **Reference:** Edwards 1961, Ann Hum Gen; Brookhart 2002, BMC Med Res Meth

### 10. Life Table
- **Inputs:** Up to 48 time periods; per period: n_x (at risk), d_x (deaths/events),
  w_x (withdrawals/losses to follow-up)
- **Outputs:** q_x, p_x, S(t), F(t), hazard, SE, lower CI, upper CI per period
- **Chart:** "Survival Curve with 90% Confidence Limits" — green CI bands,
  darker survival line; x = Time Periods, y = Survival Probability 0–1
- **Method:** Actuarial method:
  q_x = d_x / (n_x - w_x/2); p_x = 1 - q_x; S(t) = cumulative product of p_x;
  Greenwood variance; CI on log-log scale
- **R:** `survival::survfit()` or manual actuarial calculation
- **Reference:** Rothman KJ, J Chron Dis 1978;31:557-560

### 11. Matched Case-control Data
- **Inputs:** 10 matching ratio rows (1:1 through 1:10). Each row = Case
  Exposed and Case Unexposed frequencies across exposure levels (columns 0–10+)
- **Outputs (red box):** RRmh, 90%/95%/99% CI, P-value; Crude RR with CI
- **Method:** McNemar OR for 1:1; generalised matched analysis for higher ratios;
  exact binomial CI (port Module5)
- **R:** `stats::mcnemar.test()`, `epitools::oddsratio()`
- **Reference:** Modern Epidemiology 3rd ed.

### 12. O.R. Exact
- **Inputs:** Single 2×2 table: Disease/Non-Disease × Exposed/Nonexposed (4 cells)
- **Outputs:** OR point estimate; 90% and 95% CI (Mid-P lower/upper + Fisher lower/upper);
  Exact 2-tail P-value (Mid-P and Fisher)
- **Method:** Port VBA `limitcalc` (Module2) — Newton-bisection iterative search
  over hypergeometric distribution. Both Fisher exact and mid-p variants.
  Limit: n1 + n2 ≤ 20,000.
  ```r
  # Core structure of the port:
  limitcalc <- function(x1, x2, n1, n2, level, fishermid) {
    # Pre-compute log factorials
    # Iterate r (odds ratio) via bisection until tail probability = level
    # Returns r (the CI limit)
  }
  ```
- **R:** `exact2x2::fisher.exact()` covers this; validate output against VBA port
- **Reference:** Cornfield 1956, J Chron Dis

### 13. P-value Functions
- **Inputs:** Two curves; each defined by Lower Bound and Upper Bound at
  90%, 95%, AND 99% confidence levels (6 input cells total per curve)
- **Outputs:** Chart with two overlaid p-value function curves (red and yellow);
  x-axis = Relative Risk (log scale 0.01–100); y-axis = P-value (0–1);
  vertical null bar at RR=1; underlying data table
- **Method:** From CI limits, recover SE = (ln(UL) - ln(LL)) / (2 * 1.96);
  then p(RR_0) = 2 * (1 - Φ(|ln(estimate) - ln(RR_0)| / SE))
- **R:** Adapt `episheet::pvalueplot()` — already on CRAN; `ggplot2::stat_function()`
- **Reference:** Poole C, Am J Public Health 1987;77(2):195-199

### 14. Meta-Analysis
- **Inputs:** Title (text), Number of Studies (spinner, affects table rows),
  Model toggle (Fixed effects / Random effects); per study: Study label, RR,
  LL (lower CI), UL (upper CI), ln(RR) (optional), V(ln(RR)) (optional)
- **Outputs:** Per-study Relative Weight (%); Pooled row: pooled RR, LL, UL,
  ln(RR), V(ln(RR)); P for homogeneity (Cochran Q); forest plot
- **Forest plot:** Log y-axis (0.10–10.00); horizontal error bars per study;
  pooled estimate as distinct point labelled "Pooled"; study labels on x-axis
- **Model toggle:** Fixed effects = inverse-variance; Random effects = DerSimonian-Laird
- **Method:** Port all functions from VBA `Modul1`:
  - `SUM_W` → Σ(1/v_i)
  - `SUM_WLNRR` / `SUM_WRLNRR` → fixed/random pooled ln(RR)
  - `SUM_CHI2` → Cochran Q
  - `VAR_POP` → τ² (DL method)
  - `INF_W` / `INF_WR` → % influence weight per study
  - `SUM_WVLNRR` / `SUM_WRVLNRR` → pooled variance
- **R:** `meta::metagen()` or `metafor::rma()` as backend; validate against VBA
- **Reference:** Fleiss JL, Stat Methods Med Res 1993;2:121-145

---

## VBA functions to port (helpers.R)

These VBA functions have no exact CRAN equivalent and must be ported directly.
All are in the XLS file; full source was recovered.

### `limitcalc(x1, x2, n1, n2, level, fishermid)` — Module2
Exact OR confidence limits via hypergeometric bisection. Port to R verbatim.
Used by O.R. Exact tab. Limit: n1+n2 ≤ 20,000.

### `Powercalc(z, p0, r, N, RR)` — Module4
Kelsey cohort power formula. Port already shown above under Cohort Power tab.

### `ccPowercalc(z, r, p1, p0, a, b, N, RR)` — Module3
Schlesselman CC power formula. Port already shown above under CC Power tab.

### `Binomialcalc(N, x, conflevel, lowerupper, fisher)` — Module5
Exact binomial CI via iterative likelihood. Used by Quickcalc (proportion).
`lowerupper`: 0 = lower limit, 1 = upper limit.
`fisher`: 0 = mid-p, 1 = Fisher exact.

### `Poissoncalc(N, x, conflevel, lowerupper, fisher)` — Module6
Exact Poisson CI via iterative likelihood. Used by Quickcalc (rate, SMR).
Same argument convention as Binomialcalc.

### Meta-analysis functions — Modul1
`SUM_W`, `SUM_W2`, `SUM_WLNRR`, `SUM_CHI2`, `VAR_POP`, `SUM_WRLNRR`,
`SUM_WVLNRR`, `SUM_WRVLNRR`, `INF_W`, `INF_WR` — port all to R.
Full VBA source recovered; use `metafor::rma()` as validation baseline.

---

## Build order (recommended)

Start simple, prove the reactive pattern, then add complexity:

1. **Quickcalc** — single-number, self-contained, tests the exact CI ports
2. **Rate Data** — use existing `episheet` package, wrap in module
3. **Risk Data** — same
4. **Case-control Data** — introduces stratified MH logic
5. **O.R. Exact** — tests `limitcalc` port
6. **Cohort Power + CC Power** — tests power formula ports + charting
7. **Study Size** — precision-based, standalone logic
8. **Standardization** — most complex input grid (48 rows × 6 groups)
9. **Meta-Analysis** — most complex state management (dynamic row count)
10. **P-value Functions** — two-curve overlay chart
11. **Matched Case-control Data** — complex input grid
12. **Life Table** — tabular output + survival chart
13. **Seasonal Analysis** — circular stats + fitted curve

---

## Testing baseline values

Use these to validate ports against the original Excel output:

| Tab | Input | Expected output |
|-----|-------|-----------------|
| O.R. Exact | 4/386/8/1636 (Disease/Non × Exp/Nonexp) | OR=3.2383; 95% Mid-P CI: 0.7267–14.4052; 95% Fisher: 0.5997–17.4565 |
| Case-control Data | Stratum 1: Cases 104/1059, Controls 12/1163 | OR≈7.333 (crude); MH-OR shown in results |
| Seasonal Analysis | Jan–Dec counts with visible peak Dec | Peak/Low Ratio=1.819; Peak=Dec 17 (346°); 95% CI: 1.153–2.872 |
| Meta-Analysis | Studies A–D: RR 1.2/2.0/3.0/1.2 with CIs 3/7/7/4 | Pooled≈1.687 (fixed); P homogeneity≈0.496 |
| Cohort Power | z=1.645, p0=0.001, RR=4, N=1000 | Power≈0.535 (read from table) |

---

## Global clear

Implement a "Clear all inputs" action button that resets all `reactiveValues`
to `NA` / empty across all tabs. Equivalent to Ctrl+E (Macro1) in the original.
Place in the navbar, not inside individual tabs.

---

## References (include in About tab)

1. Walter SD, Elwood JM. A test for seasonality. Ann Hum Gen 1975;25:83-87.
2. Rothman KJ, Boice JD. Epidemiologic Analysis with a Programmable Calculator. 1979.
3. Dobson AJ et al. Confidence intervals for weighted sums of Poisson parameters. Stat Med 1991.
4. Fleiss JL. The statistical basis of meta-analysis. Stat Methods Med Res 1993;2:121-145.
5. Poole C. Beyond the confidence interval. Am J Public Health 1987;77(2):195-199.
6. Rothman KJ, Greenland S, Lash T. Modern Epidemiology, 3rd Ed. Lippincott, 2008.
7. Cornfield J. A statistical problem arising from retrospective studies. 1956.
8. Rothman KJ. Estimation of confidence limits for the cumulative probability of survival. J Chron Dis 1978;31:557-560.
9. Brookhart MA, Rothman KJ. Simple estimators of the intensity of seasonal occurrence. BMC Med Res Methodol 2008;8:67.

---

## Notes

- Repo: https://github.com/tom-desmet/EpiSheet
- Sheet protection password in original XLS: `"teela"` (not needed for Shiny)
- File authored by Ken Rothman; last saved Excel 12.0 (2007)
- Partial R port on CRAN: `episheet` by James Black (v0.4.2, Dec 2023)
  — covers Rate Data, Risk Data, pvalueplot only
- All other tabs must be implemented from scratch
- The `epijim/episheet` GitHub repo has the same coverage as CRAN
- Target fidelity: numerical output must match original XLS to 4 significant figures
- Do not invent epidemiological formulas — every method has a reference listed above;
  if uncertain about a formula, ask rather than approximate

---

## Deployment

Target: **shinyapps.io** (free tier sufficient for development).

```r
# When ready to deploy:
install.packages("rsconnect")
rsconnect::deployApp(appDir = ".", appName = "EpiSheet")
```

Add a `renv.lock` from the start so the deployment environment is reproducible:

```r
install.packages("renv")
renv::init()
# after installing all packages:
renv::snapshot()
```
