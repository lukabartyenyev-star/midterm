# =============================================================================
# STEP 1 — Load the source file
# =============================================================================

source("carry_hmm3.R")


# =============================================================================
# STEP 2 — Load and prepare per-country yield + FX data
# =============================================================================
# Yields CSV: date, X2Y, X10Y  (semicolon-separated, comma decimal)
# FX CSV:     date, rate        (semicolon-separated, comma decimal)

countries <- c("GBP", "JPY", "CAD", "SEK")

country_data <- lapply(setNames(countries, countries), function(cc) {
  yields <- read.csv2(paste0("curves_data_csv/curve_", cc, ".csv"))
  fx     <- read.csv2(paste0("fx_data_csv/fx_",        cc, ".csv"))
  prepare_country(yields, fx)
})

# Sanity check — should show date, d_y2y, d_spread, fx_return, no NAs
lapply(country_data, function(cd) {
  data.frame(
    rows        = nrow(cd),
    from        = min(cd$date),
    to          = max(cd$date),
    na_d_y2y    = sum(is.na(cd$d_y2y)),
    na_d_spread = sum(is.na(cd$d_spread)),
    na_fx       = sum(is.na(cd$fx_return))
  )
})

# Confirm columns are numeric, not character
str(country_data[["GBP"]])

# Check for gaps in the time series
lapply(country_data, function(cd) {
  gaps <- diff(as.Date(cd$date))
  data.frame(
    max_gap_days  = max(gaps),
    gaps_over_5d  = sum(gaps > 5),
    gaps_over_20d = sum(gaps > 20)
  )
})


# =============================================================================
# STEP 2b — Align all countries to common start date
# =============================================================================

start_date <- max(sapply(country_data, function(cd) min(cd$date)))
cat("Common start date:", start_date, "\n")

country_data <- lapply(country_data, function(cd) {
  cd[cd$date >= start_date, ]
})


# =============================================================================
# STEP 3 — Load covariates
# =============================================================================
# TED CSV: date, ted  (other columns are ignored)
# VIX CSV: date, vix  (other columns are ignored)

ted <- read.csv2("TED.csv")
vix <- read.csv2("VIX.csv")


# =============================================================================
# STEP 4 — Fit models
# =============================================================================

hmm_results <- fit_all_countries(country_data, ted, vix, N = 2, cyc = 100)


# =============================================================================
# STEP 5 — Check state identification
# =============================================================================
# Unwind state should have:
#   - negative fx_return mean  (currency depreciating)
#   - positive d_y2y mean      (short rates rising, carry being unwound)

summarise_states(hmm_results)

# Override manually if a country is misidentified, e.g.:
# hmm_results[["JPY"]]$unwind_state <- 1
# hmm_results[["JPY"]]$unwind_prob  <- hmm_results[["JPY"]]$decoded$prb[1, ]


# =============================================================================
# STEP 6 — Extract regime probabilities
# =============================================================================

regime_probs <- lapply(names(hmm_results), function(cc) {
  r <- hmm_results[[cc]]
  data.frame(date = as.Date(r$dates), unwind_prob = r$unwind_prob)
})
names(regime_probs) <- names(hmm_results)

# Quick look at probabilities for each country
lapply(regime_probs, head, 10)


# =============================================================================
# STEP 7 — Plot unwind probabilities
# =============================================================================

par(mfrow = c(ceiling(length(countries) / 2), 2), mar = c(2, 3, 2, 1))

for (cc in countries) {
  r <- hmm_results[[cc]]
  plot(as.Date(r$dates), r$unwind_prob,
       type = "l", lwd = 1.2, col = "#D6604D",
       ylim = c(0, 1), main = cc,
       xlab = "", ylab = "P(unwind)")
  abline(h = 0.6, lty = 2, col = "grey50")
  grid()
}