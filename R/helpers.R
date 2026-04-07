# helpers.R — VBA function ports and shared utilities

# ── Exact OR CI (Module2: limitcalc) ─────────────────────────────────────────
# Port of VBA limitcalc(): exact OR confidence limits via hypergeometric bisection
# x1 = exposed cases, x2 = unexposed cases, n1 = exposed total, n2 = unexposed total
# level = one-sided alpha (e.g. 0.025 for 95% CI)
# fishermid: TRUE = Fisher exact, FALSE = mid-p
limitcalc <- function(x1, x2, n1, n2, level, fishermid = TRUE) {
  # Pre-compute log factorials up to max needed
  maxn <- n1 + n2 + 1
  lf <- c(0, cumsum(log(seq_len(maxn))))  # lf[i+1] = log(i!)

  lncf <- function(n, k) {
    if (k < 0 || k > n) return(-Inf)
    lf[n + 1] - lf[k + 1] - lf[n - k + 1]
  }

  # Hypergeometric probability P(X = x | margins, OR = r)
  hyp_prob <- function(x, r) {
    # sum over all valid x' to normalise
    x_min <- max(0, x1 + x2 - n2)
    x_max <- min(x1 + x2, n1)
    log_r <- log(r)
    log_probs <- sapply(x_min:x_max, function(k) {
      lncf(n1, k) + lncf(n2, x1 + x2 - k) + k * log_r
    })
    # normalise
    mx <- max(log_probs)
    probs <- exp(log_probs - mx)
    probs / sum(probs)
  }

  tail_prob <- function(r, lower_tail) {
    x_min <- max(0, x1 + x2 - n2)
    x_max <- min(x1 + x2, n1)
    probs <- hyp_prob(x1, r)
    xs <- x_min:x_max
    p_x1 <- probs[xs == x1]
    if (lower_tail) {
      p <- sum(probs[xs <= x1])
    } else {
      p <- sum(probs[xs >= x1])
    }
    if (!fishermid) {
      p <- p - p_x1 / 2
    }
    p
  }

  # Determine if we want lower or upper limit
  # For lower CI limit: find r s.t. upper tail = level
  # For upper CI limit: find r s.t. lower tail = level
  # Use bisection
  bisect <- function(f, lo, hi, tol = 1e-7) {
    for (i in seq_len(100)) {
      mid <- (lo + hi) / 2
      if (f(mid) > 0) hi <- mid else lo <- mid
      if (hi - lo < tol) break
    }
    (lo + hi) / 2
  }

  # Lower limit: upper tail prob = level
  r_lo <- bisect(
    function(r) tail_prob(r, lower_tail = FALSE) - level,
    lo = 1e-7, hi = ifelse(x1 == 0, 1, 1e6)
  )

  # Upper limit: lower tail prob = level
  r_hi <- bisect(
    function(r) tail_prob(r, lower_tail = TRUE) - level,
    lo = ifelse(x1 == n1, 1, 1e-7), hi = 1e6
  )

  list(lower = r_lo, upper = r_hi)
}

# ── Exact binomial CI (Module5: Binomialcalc) ────────────────────────────────
# N = denominator, x = numerator, conflevel = e.g. 0.95
# lowerupper: "lower" or "upper"
# fisher: TRUE = Fisher exact, FALSE = mid-p
binomialcalc <- function(N, x, conflevel, lowerupper, fisher = TRUE) {
  alpha <- 1 - conflevel
  level <- alpha / 2

  if (x == 0 && lowerupper == "lower") return(0)
  if (x == N && lowerupper == "upper") return(1)

  # tail probability of Binomial(N, p) at x
  tail_p <- function(p, lower_tail) {
    if (lower_tail) {
      prob <- pbinom(x, N, p)
    } else {
      prob <- 1 - pbinom(x - 1, N, p)
    }
    if (!fisher) {
      prob <- prob - dbinom(x, N, p) / 2
    }
    prob
  }

  if (lowerupper == "lower") {
    # lower CI: find p s.t. upper tail = alpha/2
    uniroot(function(p) tail_p(p, lower_tail = FALSE) - level,
            lower = 1e-10, upper = x / N + 0.5, tol = 1e-8)$root
  } else {
    # upper CI: find p s.t. lower tail = alpha/2
    uniroot(function(p) tail_p(p, lower_tail = TRUE) - level,
            lower = x / N - 0.5, upper = 1 - 1e-10, tol = 1e-8)$root
  }
}

# ── Exact Poisson CI (Module6: Poissoncalc) ───────────────────────────────────
# N = person-time, x = count (events), conflevel = e.g. 0.95
poissoncalc <- function(N, x, conflevel, lowerupper, fisher = TRUE) {
  alpha <- 1 - conflevel
  level <- alpha / 2
  rate <- x / N  # rate per unit time

  if (x == 0 && lowerupper == "lower") return(0)

  # tail probability for Poisson(mu = lambda * N) at x
  tail_p <- function(lambda, lower_tail) {
    mu <- lambda * N
    if (lower_tail) {
      prob <- ppois(x, mu)
    } else {
      prob <- 1 - ppois(x - 1, mu)
    }
    if (!fisher) {
      prob <- prob - dpois(x, mu) / 2
    }
    prob
  }

  if (lowerupper == "lower") {
    if (x == 0) return(0)
    uniroot(function(lam) tail_p(lam, lower_tail = FALSE) - level,
            lower = 1e-10, upper = rate * 10 + 10 / N, tol = 1e-10)$root
  } else {
    uniroot(function(lam) tail_p(lam, lower_tail = TRUE) - level,
            lower = 1e-10, upper = rate * 10 + 10 / N, tol = 1e-10)$root
  }
}

# ── Cohort power (Module4: Powercalc) ────────────────────────────────────────
# z_alpha = 1.645 for one-sided 0.05; p0 = risk in unexposed; r = ratio N_exp/N_unexp
# N = number exposed; RR = relative risk
powercalc <- function(z_alpha, p0, r, N, RR) {
  p1 <- RR * p0
  if (p1 >= 1 || p0 <= 0 || p1 <= 0) return(NA_real_)
  a <- p1 * (1 - p0) + p0 * (1 - p1)
  b <- (r - 1) * p0 * (1 - p0)
  denom <- sqrt((a + b) * (r * a - b) - r * (p1 - p0)^2)
  if (is.nan(denom) || denom == 0) return(NA_real_)
  zbnum <- -z_alpha * (a + b) + abs(p1 - p0) * sqrt(r * N * (a + b))
  zbden <- denom
  pnorm(zbnum / zbden)
}

# ── Case-control power (Module3: ccPowercalc) ────────────────────────────────
# r = controls/cases ratio; p0 = exposure prevalence in controls
# N = number of cases; OR = odds ratio
cc_powercalc <- function(z_alpha, r, p0, N, OR) {
  p1 <- (OR * p0) / (1 + p0 * (OR - 1))
  if (p1 >= 1 || p0 <= 0 || p1 <= 0) return(NA_real_)
  a <- p1 * (1 - p0) + p0 * (1 - p1)
  b <- (r - 1) * p0 * (1 - p0)
  denom <- sqrt((a + b) * (r * a - b) - r * (p1 - p0)^2)
  if (is.nan(denom) || denom == 0) return(NA_real_)
  zbnum <- -z_alpha * (a + b) + abs(p1 - p0) * sqrt(r * N * (a + b))
  pnorm(zbnum / denom)
}

# ── Meta-analysis VBA ports (Modul1) ─────────────────────────────────────────
# Input: data frame with columns: lnRR (log RR), v (variance of log RR)

sum_w   <- function(v)            sum(1 / v, na.rm = TRUE)
sum_w2  <- function(v)            sum(1 / v^2, na.rm = TRUE)
sum_wlnrr  <- function(lnRR, v)  sum(lnRR / v, na.rm = TRUE)
sum_chi2   <- function(lnRR, v)  sum((lnRR - sum_wlnrr(lnRR, v) / sum_w(v))^2 / v, na.rm = TRUE)

# DerSimonian-Laird tau^2
var_pop <- function(lnRR, v) {
  k  <- sum(!is.na(lnRR))
  Q  <- sum_chi2(lnRR, v)
  W  <- sum_w(v)
  W2 <- sum_w2(v)
  tau2 <- max(0, (Q - (k - 1)) / (W - W2 / W))
  tau2
}

sum_wrlnrr <- function(lnRR, v) {
  tau2 <- var_pop(lnRR, v)
  vr <- v + tau2
  sum(lnRR / vr, na.rm = TRUE)
}

sum_wvlnrr <- function(lnRR, v) {
  1 / sum_w(v)
}

sum_wrvlnrr <- function(lnRR, v) {
  tau2 <- var_pop(lnRR, v)
  vr <- v + tau2
  1 / sum(1 / vr, na.rm = TRUE)
}

inf_w <- function(v) {
  w <- 1 / v
  100 * w / sum(w, na.rm = TRUE)
}

inf_wr <- function(lnRR, v) {
  tau2 <- var_pop(lnRR, v)
  wr <- 1 / (v + tau2)
  100 * wr / sum(wr, na.rm = TRUE)
}

# ── Utility: safe division ────────────────────────────────────────────────────
safe_div <- function(a, b) ifelse(b == 0, NA_real_, a / b)

# ── Format helper ─────────────────────────────────────────────────────────────
fmt <- function(x, digits = 4) {
  if (is.na(x) || is.nan(x) || is.infinite(x)) return("—")
  formatC(x, digits = digits, format = "f")
}

fmt4 <- function(x) fmt(x, 4)
fmt2 <- function(x) fmt(x, 2)

orange_val <- function(x, digits = 4) {
  tags$span(class = "result-orange", fmt(x, digits))
}

red_val <- function(x, digits = 4) {
  tags$span(class = "result-red", fmt(x, digits))
}
