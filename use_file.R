# =============================================================================
# STEP 1 — Load the source file
# =============================================================================

source("C:/Users/WZHLUBD/Downloads/carry_hmm3.R")
library(dplyr)


# =============================================================================
# STEP 2 — Load and prepare per-country yield + FX data
# =============================================================================
# Yields CSV: date, X2Y, X10Y  (semicolon-separated, comma decimal)
# FX CSV:     date, rate        (semicolon-separated, comma decimal)

index_mapping <- read.csv("C:/Users/WZHLUBD/Downloads/mapping_index.csv", stringsAsFactors = FALSE, header = FALSE)
fx_rates <- read.csv("C:/Users/WZHLUBD/Downloads/fx_rates.csv") %>%
  mutate(Date = as.Date(Date, format = "%d.%m.%Y"))%>%
  rename(date = Date)


yields <- lapply(list.files("Yield data", full.names = TRUE), read.csv2)
names(yields)<-gsub(list.files("Yield data", full.names = F)     , pattern = "_wide.csv", replacement = "")# msp this from V2 in in
names(yields) <- index_mapping$V3[match(names(yields), index_mapping$V2)]

yields <- lapply(yields, function(df) {
  df %>%
    filter(type == "Mid Yield") %>%
    mutate(date = as.Date(date),
           across(c(-date, -type), ~ ifelse(abs(.x) < 0.00001, NA, .x)))
})

View(prepare_country)
country_data<-lapply(index_mapping$V3,function(cc){
  term_str <- yields[[cc]]
  fx     <- fx_rates[,c("date",cc)]%>%
    rename(rate = cc)
  
  prepare_country(term_str,fx)%>%
    filter(date > as.Date("2002-01-03"))
  }
  
)


names(country_data)<-paste0(index_mapping$V3," ",index_mapping$V2)
names(country_data)[13]

View(country_data)

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
