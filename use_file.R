# =============================================================================
# STEP 1 — Load the source file
# =============================================================================

setwd("C:/Users/1/Downloads")
source("C:/Users/1/Downloads/Data Model/carry_hmm3.R")
library(dplyr)



# =============================================================================
# STEP 2 — Load and prepare per-country yield + FX data
# =============================================================================
# Yields CSV: date, X2Y, X10Y  (semicolon-separated, comma decimal)
# FX CSV:     date, rate        (semicolon-separated, comma decimal)

index_mapping <- read.csv("Data Model/mapping_index.csv", stringsAsFactors = FALSE, header = FALSE)
fx_rates <- read.csv("Data Model/fx_rates.csv") %>%
  mutate(Date = as.Date(Date, format = "%d.%m.%Y"))%>%
  rename(date = Date)


yields <- lapply(list.files("Data Model/Yield data", full.names = TRUE), read.csv2)
names(yields)<-gsub(list.files("Data Model/Yield data", full.names = F)     , pattern = "_wide.csv", replacement = "")# msp this from V2 in in
names(yields) <- index_mapping$V3[match(names(yields), index_mapping$V2)]

yields <- lapply(yields, function(df) {
  df %>%
    filter(type == "Mid Yield") %>%
    mutate(date = as.Date(date),
           across(c(-date, -type), ~ ifelse(abs(.x) < 0.00001, NA, .x)))%>%filter(date > as.Date("2006-01-03")) #consider uniform date
})


# Sanity check — should show date, d_y2y, d_spread, fx_return, no NAs


safe_na <- function(df, col) {
  if (col %in% names(df)) {
    is.na(df[[col]])
  } else {
    rep(TRUE, nrow(df))  # treat missing column as all NA
  }
}
DQI <- lapply(yields, function(cd) {
  
  na_3y  <- safe_na(cd, "X3Y")
  na_4y  <- safe_na(cd, "X4Y")
  na_6m  <- safe_na(cd, "X6M")
  na_9m  <- safe_na(cd, "X9M")
  na_1y  <- safe_na(cd, "X1Y")
  na_2y  <- safe_na(cd, "X2Y")
  
  na_8y  <- safe_na(cd, "X8Y")
  na_9y  <- safe_na(cd, "X9Y")
  na_12y <- safe_na(cd, "X12Y")
  na_15y <- safe_na(cd, "X15Y")
  na_10y <- safe_na(cd, "X10Y")
  
  x <- data.frame(
    date = cd$date,
    
    na_y2y = na_2y,
    
    na_we_are_done =
      !(
        (
          (!na_3y | !na_4y) &
            (!na_9m | !na_1y)
        ) |
          !na_2y
      ),
    
    na_y10y = na_10y,
    
    na_we_are_done_10 =
      !(
        (
          (!na_9y) &
            (!na_12y)
        ) |
          !na_10y
      )
  )
  
  merge(cd, x, by = "date") %>%
    dplyr::filter(na_y2y)
})
DQI_fitted <- lapply(yields_filled, function(cd) {
  
  na_3y  <- safe_na(cd, "X3Y")
  na_4y  <- safe_na(cd, "X4Y")
  na_6m  <- safe_na(cd, "X6M")
  na_9m  <- safe_na(cd, "X9M")
  na_1y  <- safe_na(cd, "X1Y")
  na_2y  <- safe_na(cd, "X2Y")
  
  na_8y  <- safe_na(cd, "X8Y")
  na_9y  <- safe_na(cd, "X9Y")
  na_12y <- safe_na(cd, "X12Y")
  na_15y <- safe_na(cd, "X15Y")
  na_10y <- safe_na(cd, "X10Y")
  
  x <- data.frame(
    date = cd$date,
    
    na_y2y = na_2y,
    
    na_we_are_done =
      !(
        (
          (!na_3y | !na_4y) &
            (!na_9m | !na_1y)
        ) |
          !na_2y
      ),
    
    na_y10y = na_10y,
    
    na_we_are_done_10 =
      !(
        (
          (!na_9y) &
            (!na_12y)
        ) |
          !na_10y
      )
  )
  
  merge(cd, x, by = "date") %>%
    dplyr::filter(na_we_are_done)
})


library(tidyverse)
library(zoo)
library(YieldCurve)

tenor_to_years <- function(t) {
  case_when(
    grepl("W$", t) ~ as.numeric(sub("W","",t)) / 52,
    grepl("M$", t) ~ as.numeric(sub("M","",t)) / 12,
    grepl("Y$", t) ~ as.numeric(sub("Y","",t)),
    TRUE           ~ NA_real_
  )
}

get_tenor_cols <- function(df) names(df)[grepl("^X", names(df))]
col_to_tenor   <- function(col) sub("^X", "", col)
tenor_to_years(col_to_tenor(get_tenor_cols(yield_test)))


roll_col <- function(vec, max_days = 3) {
  filled <- zoo::na.locf(vec, na.rm = FALSE, maxgap = max_days)
  zoo::na.locf(filled, fromLast = TRUE, na.rm = FALSE, maxgap = max_days)
}

# ── NS fit for one row, returns fitted value at target maturity ───────────────
library(xts)

fit_ns_at <- function(named_vec, target_years, min_obs = 4) {
  x   <- tenor_to_years(names(named_vec))
  y   <- as.numeric(named_vec)
  obs <- !is.na(y) & !is.na(x)
  
  if (sum(obs) < min_obs)  return(rep(NA_real_, length(target_years)))
  if (!any(x[obs] <= 2))   return(rep(NA_real_, length(target_years)))
  if (!any(x[obs] >= 5))   return(rep(NA_real_, length(target_years)))
  
  tryCatch({
    # Nelson.Siegel needs a 1-row xts — fake date is irrelevant
    rate_xts <- xts(matrix(y[obs], nrow = 1), order.by = as.Date("2000-01-01"))
    fit      <- Nelson.Siegel(rate = rate_xts, maturity = x[obs])
    fitted   <- NSrates(fit, maturity = target_years)
    as.numeric(fitted[1, ])
  }, error = function(e) {
    cat(sprintf("  NS failed: %s\n", conditionMessage(e)))
    rep(NA_real_, length(target_years))
  })
}




fill_targets <- function(df, target_cols, max_roll = 3) {
  
  tenor_cols   <- get_tenor_cols(df)
  avail_cols   <- intersect(tenor_cols, names(df))
  targets_here <- intersect(target_cols, names(df))
  target_years <- tenor_to_years(col_to_tenor(targets_here))
  
  sub_df <- df %>% arrange(date)
  
  # ── Step 1: roll each target up to max_roll days ──────────────────────────
  for (tc in targets_here) {
    sub_df[[tc]] <- roll_col(sub_df[[tc]], max_days = max_roll)
  }
  
  # ── Step 2: NS on unrolled support tenors ────────────────────────────────
  still_na <- sapply(targets_here, function(tc) is.na(sub_df[[tc]]))
  if (!is.matrix(still_na))
    still_na <- matrix(still_na, ncol = length(targets_here),
                       dimnames = list(NULL, targets_here))
  
  for (i in which(apply(still_na, 1, any))) {
    row_vals <- setNames(as.numeric(sub_df[i, avail_cols]),
                         col_to_tenor(avail_cols))
    fitted   <- fit_ns_at(row_vals, target_years)
    
    for (j in seq_along(targets_here)) {
      tc <- targets_here[j]
      if (still_na[i, tc] && !is.na(fitted[j]))
        sub_df[i, tc] <- fitted[j]
    }
  }
  
  # ── Step 3: roll support tenors, retry NS ────────────────────────────────
  still_na2 <- sapply(targets_here, function(tc) is.na(sub_df[[tc]]))
  if (!is.matrix(still_na2))
    still_na2 <- matrix(still_na2, ncol = length(targets_here),
                        dimnames = list(NULL, targets_here))
  
  rows_still_na <- which(apply(still_na2, 1, any))
  
  if (length(rows_still_na) > 0) {
    temp_df <- sub_df
    for (nc in setdiff(avail_cols, targets_here)) {
      temp_df[[nc]] <- roll_col(temp_df[[nc]], max_days = max_roll)
    }
    
    for (i in rows_still_na) {
      row_vals <- setNames(as.numeric(temp_df[i, avail_cols]),
                           col_to_tenor(avail_cols))
      fitted   <- fit_ns_at(row_vals, target_years)
      
      for (j in seq_along(targets_here)) {
        tc <- targets_here[j]
        if (still_na2[i, tc] && !is.na(fitted[j]))
          sub_df[i, tc] <- fitted[j]
      }
    }
  }
  
  # ── Report residual NAs ───────────────────────────────────────────────────
  for (tc in targets_here) {
    residual <- sum(is.na(sub_df[[tc]]))
    if (residual > 0)
      cat(sprintf("  [%s] %d NAs remain after all 3 steps\n", tc, residual))
  }
  
  sub_df
}




# ── Run ───────────────────────────────────────────────────────────────────────
target_cols <- c("X2Y","X3Y", "X10Y")

yields_filled <- imap(yields, function(df, ccy) {
  cat("Processing:", ccy, "\n")
  fill_targets(df, target_cols, max_roll = 3)
})

good_list<-c("AUD", "CAD", "EUR", "JPY", "SEK", "GBP", "USD","NZD","NOK", "CHF", "CZK", "PLN","MXN","BRL")
length(good_list)


checks<-lapply(yields_filled,function(cd){
  data.frame(na10=sum(is.na(cd$X10Y)),
             na2=sum(is.na(cd$X2Y)),
             na3=sum(is.na(cd$X3Y)))
  }
)


country_data<-lapply(good_list,function(cc){
  term_str <- yields_filled[[cc]]
  fx     <- fx_rates[,c("date",cc)]%>%
    rename(rate = cc)
  
  prepare_country(term_str,fx)[-1,]
}
)


names(country_data)<-good_list

# Confirm columns are numeric, not character


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
# STEP 3 — Load covariates
# =============================================================================
# TED CSV: date, ted  (other columns are ignored)
# VIX CSV: date, vix  (other columns are ignored)

library(dplyr)



ted <- read.csv2("C:/Users/1/Downloads/TED.csv")
vix <- read.csv2("C:/Users/1/Downloads/VIX.csv")


#' Fit TVTP-HMM to one country's data
#'
#' @param X     Matrix (T x 3): d_y2y, d_spread, fx_return
#' @param Z     Matrix (T x 2): standardised TED, VIX
#' @param N     Number of hidden states
#' @param cyc   Maximum EM iterations
#' @param tol   Convergence tolerance on log-likelihood change
#' @param seed  RNG seed
#'
#' @return List: Mu (N x 3), Cov (3 x 3 x N), x (N x N x 3), Pi (N), LL

# =============================================================================
# STEP 4 — Fit models
# =============================================================================

result_good_run<-fit_all_countries(country_data,ted,vix)
good_SAVE <- "C:/Users/1/Documents/GOODRES.rds"  # same folder as your other EWS files
saveRDS(result_good_run,good_SAVE)
processed_hmms<-readRDS(good_SAVE)

View(processed_hmms)

unwind_states<-unlist(lapply(processed_hmms,function(x){which.max(x$model$Cov[3,3,]) }))
names(unwind_states)<-names(processed_hmms)


transition_prob_betas<-lapply(processed_hmms,function(x){
  tuple <- x$model$x
  data.frame(
    transition  = c("1->2", "2->1"),
    intercept   = c(tuple[1, 2, 1], tuple[2, 1, 1]),
    beta_TED    = c(tuple[1, 2, 2], tuple[2, 1, 2]),
    beta_VIX    = c(tuple[1, 2, 3], tuple[2, 1, 3])
    )
})

cov_matrices<-lapply(processed_hmms,function(x){
  matr<-x$model$Cov
  colnames(matr)<-rownames(matr)<-c("d_y2y", "d_spread", "fx_return")
  matr
  })
Mus<-lapply(processed_hmms,function(x){
  mus<-x$model$Mu
  colnames(mus)<-c("d_y2y", "d_spread", "fx_return")
  mus
  })

wald_tables<-lapply(processed_hmms,function(x){
    x$model$wald_table
  }
)

state_history<-lapply(processed_hmms,function(x){
  data.frame(
    dates=as.Date(x$date),
    viterbi_path=x$model$viterbi_path$path,
    viterbi_probs=x$model$viterbi_path$path_probs,
    state_sequence=x$model$state_sequence,
    smoothed_probs1=x$model$smoothed_probs[,1],
    moothed_probs1=x$model$smoothed_probs[,2]
    )
  }
)

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
