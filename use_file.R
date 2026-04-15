# =============================================================================
# STEP 1 — Load the source file
# =============================================================================

setwd("C:/Users/1/Downloads")
source("Data Model/carry_hmm3.R")
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

View(prepare_country)
country_data<-lapply(index_mapping$V3,function(cc){
  term_str <- yields[[cc]]
  fx     <- fx_rates[,c("date",cc)]%>%
    rename(rate = cc)
  
  prepare_country(term_str,fx)
  }
  
)


names(country_data)<-paste0(index_mapping$V3," ",index_mapping$V2)
head(country_data[[1]])

View(country_data)

# Sanity check — should show date, d_y2y, d_spread, fx_return, no NAs

lapply(yields, function(cd) {
  data.frame(
    rows        = nrow(cd),
    from        = min(cd$date),
    to          = max(cd$date),
    na_y2y    = sum(is.na(cd$X2Y)),
    na_we_are_done = sum(
     !(
       (
         (!is.na(cd$X3Y) | !is.na(cd$X4Y)) &
           (!is.na(cd$X9M) | !is.na(cd$X1Y))
       ) |
         !is.na(cd$X2Y)
     )
   ),
   na_y10y    = sum(is.na(cd$X10Y)),
   na_we_are_done_10 = sum(
     !(
       (
         (!is.na(cd$X8Y) | !is.na(cd$X9Y)) &
           (!is.na(cd$X12M) | !is.na(cd$X15Y))
       ) |
         !is.na(cd$X10Y)
     )
   )
  )
})


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
    dplyr::filter(na_we_are_done)
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
lapply(DQI_fitted,nrow)
lapply(DQI,nrow)


lapply(DQI,function(x)diff(x[,"date"]))

View(DQI$NZD)
all_cols <- unique(unlist(lapply(yields, names)))





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

roll_col <- function(vec, max_days = 3) {
  filled <- zoo::na.locf(vec, na.rm = FALSE, maxgap = max_days)
  zoo::na.locf(filled, fromLast = TRUE, na.rm = FALSE, maxgap = max_days)
}

# ── NS fit for one row, returns fitted value at target maturity ───────────────
fit_ns_at <- function(named_vec, target_years, min_obs = 4) {
  x   <- tenor_to_years(names(named_vec))
  y   <- as.numeric(named_vec)
  obs <- !is.na(y) & !is.na(x)
  
  if (sum(obs) < min_obs)      return(NA_real_)
  if (!any(x[obs] <= 2))       return(NA_real_)  # need short end to identify beta1
  if (!any(x[obs] >= 5))       return(NA_real_)  # need long end to identify beta0
  
  tryCatch({
    fit    <- Nelson.Siegel(rate = y[obs], maturity = x[obs])
    fitted <- NSrates(fit, maturity = target_years)
    as.numeric(fitted[1, ])
  }, error = function(e) rep(NA_real_, length(target_years)))
}

# ── Main fill procedure ───────────────────────────────────────────────────────
fill_targets <- function(df, target_cols, max_roll = 3) {
  
  tenor_cols <- get_tenor_cols(df)
  types      <- unique(df$type)
  
  result <- map(types, function(tp) {
    sub_df <- df %>% filter(type == tp) %>% arrange(date)
    
    # Pre-compute target maturities in years (needed for NS)
    targets_here  <- intersect(target_cols, names(sub_df))
    target_years  <- tenor_to_years(col_to_tenor(targets_here))
    avail_cols    <- intersect(tenor_cols, names(sub_df))
    
    # ── Step 1: roll each target up to max_roll days ──────────────────────
    for (tc in targets_here) {
      sub_df[[tc]] <- roll_col(sub_df[[tc]], max_days = max_roll)
    }
    
    # ── Step 2: NS interpolation using currently available (unrolled) tenors
    still_na <- sapply(targets_here, function(tc) is.na(sub_df[[tc]]))
    # still_na is matrix: rows = observations, cols = targets
    if (!is.matrix(still_na)) still_na <- matrix(still_na, ncol = length(targets_here),
                                                 dimnames = list(NULL, targets_here))
    
    rows_needing_fill <- which(apply(still_na, 1, any))
    
    for (i in rows_needing_fill) {
      row_vals <- setNames(as.numeric(sub_df[i, avail_cols]),
                           col_to_tenor(avail_cols))
      fitted   <- fit_ns_at(row_vals, target_years)
      
      for (j in seq_along(targets_here)) {
        tc <- targets_here[j]
        if (still_na[i, tc] && !is.na(fitted[j])) {
          sub_df[i, tc] <- fitted[j]
        }
      }
    }
    
    # ── Step 3: still NA — roll all other tenor cols, then retry NS ───────
    still_na2 <- sapply(targets_here, function(tc) is.na(sub_df[[tc]]))
    if (!is.matrix(still_na2)) still_na2 <- matrix(still_na2, ncol = length(targets_here),
                                                   dimnames = list(NULL, targets_here))
    
    rows_still_na <- which(apply(still_na2, 1, any))
    
    if (length(rows_still_na) > 0) {
      # Roll all non-target tenor columns
      non_target_cols <- setdiff(avail_cols, targets_here)
      temp_df <- sub_df
      for (nc in non_target_cols) {
        temp_df[[nc]] <- roll_col(temp_df[[nc]], max_days = max_roll)
      }
      
      for (i in rows_still_na) {
        row_vals <- setNames(as.numeric(temp_df[i, avail_cols]),
                             col_to_tenor(avail_cols))
        fitted   <- fit_ns_at(row_vals, target_years)
        
        for (j in seq_along(targets_here)) {
          tc <- targets_here[j]
          if (still_na2[i, tc] && !is.na(fitted[j])) {
            sub_df[i, tc] <- fitted[j]
          }
        }
      }
    }
    
    # ── Report residual NAs ───────────────────────────────────────────────
    for (tc in targets_here) {
      residual <- sum(is.na(sub_df[[tc]]))
      if (residual > 0)
        cat(sprintf("  [%s | %s | %s] %d NAs remain after all 3 steps\n",
                    tp, tc, unique(sub_df$curve), residual))
    }
    
    sub_df
  })
  
  bind_rows(result) %>% arrange(date, type)
}

# ── Compute emissions ─────────────────────────────────────────────────────────
compute_emissions <- function(df, short_col = "X2Y", long_col = "X10Y") {
  df %>%
    filter(type == "Mid Yield") %>%
    arrange(date) %>%
    mutate(
      slope   = .data[[long_col]] - .data[[short_col]],
      d_short = .data[[short_col]] - lag(.data[[short_col]]),
      d_slope = slope - lag(slope)
    ) %>%
    select(date, short = all_of(short_col), long = all_of(long_col),
           slope, d_short, d_slope)
}

# ── Run ───────────────────────────────────────────────────────────────────────
target_cols <- c("X2Y", "X10Y")

yields_filled <- imap(yields, function(df, ccy) {
  cat("Processing:", ccy, "\n")
  fill_targets(df, target_cols, max_roll = 3)
})
























View(yields$MYR)
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

library(dplyr)

do_country <- function(yields, fx) {
  
  df <- merge(yields, fx, by = "date", sort = TRUE)
  
  # Forward-fill FX rate across gaps before computing returns
  for (i in seq(2, nrow(df)))
    if (is.na(df$rate[i])) df$rate[i] <- df$rate[i - 1]
  
  spread       <- df$X10Y - df$X2Y
  df$d_y2y     <- c(NA, diff(df$X2Y))
  df$d_spread  <- c(NA, diff(spread))
  df$fx_return <- c(NA, diff(log(df$rate)))

  df[, c("date", "d_y2y", "d_spread", "fx_return")]
}





curves<-list.files("midterm/curves_data_csv")
forex<-list.files("midterm/fx_data_csv")


yield <-



ted <- read.csv2("C:/Users/1/Downloads/TED.csv")
vix <- read.csv2("C:/Users/1/Downloads/VIX.csv")

critical<- as.Date("2002-04-10")



test<-do_country(yield,fx)%>%
  mutate(date = as.Date(date,format="%Y-%m-%d"))%>%
  filter(date>critical)


X<- as.matrix(test[-1])

y<-merge(ted,vix,by="date")%>%
  mutate(date = as.Date(date,format="%Y-%m-%d"))%>%
  select(-US0003M.Index,-GB3.Govt)%>%
  filter(date %in% test$date)

y_df<-merge(ted,vix,by="date")%>%
  mutate(date = as.Date(date,format="%Y-%m-%d"))%>%
  select(-US0003M.Index,-GB3.Govt)%>%
  mutate(ted_lag = lag(ted),
         vix_lag = lag(vix))%>%
  filter(date %in% test$date)

y<- as.matrix(y_df[,c(-1,-4,-5)])
y_lag<- as.matrix(y_df[,-1:-3])
View(y_lag)



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


check1<-fit_country_hmm(X,y)

reg_1<-1/(1+exp(rep(check1$x[1,2,1],nrow(y)) + check1$x[1,2,2]*y[,1] + check1$x[1,2,3]**y[,2]))

beta_df <- data.frame(
  transition  = c("1->2", "2->1"),
  intercept   = c(check1$x[1, 2, 1], check1$x[2, 1, 1]),
  beta_TED    = c(check1$x[1, 2, 2], check1$x[2, 1, 2]),
  beta_VIX    = c(check1$x[1, 2, 3], check1$x[2, 1, 3])
)
par(mfrow=c(2,1))
i<-2
plot(y_df$date,1/(1+exp(-rep(beta_df$intercept[i],nrow(y)) - beta_df$beta_TED[i]*y[,1] - beta_df$beta_VIX[i]*y[,2])),
     main=beta_df$transition[i])

cat(beta_df$transition[i])

reg_2<-check1$x[2,1,]

check1$Cov

View(check1$x)

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
