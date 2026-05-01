library(dplyr)
library(zoo)
library(mvtnorm)
library(TTR)

setwd("/Users/katagiri/Documents/Midterm Pitch SS2026/folder 2")
source("carry_hmm3.R")

cutoff <- as.Date("2019-07-01")
processed_hmm_full <- readRDS("result_nzd_train_2019.rds")
processed_hmm_full[["NZD"]]$unwind_state <- 2

ted <- read.csv2("TED.csv")
vix <- read.csv2("VIX.csv")

# fx prep (NZD per 1 ZAR)
fx_rates <- read.csv("Data Model/fx_rates.csv") %>%
  mutate(Date = as.Date(Date, format = "%d.%m.%Y")) %>% rename(date = Date)

nzd_zar_fx <- fx_rates %>% select(date, NZD, ZAR) %>%
  mutate(rate = NZD / ZAR) %>% select(date, rate)

# load overnight rates (as always check ur directories)
policy_nzd <- read.csv("Data Model/nzd_overnight_2018.csv") %>%
  mutate(date = as.Date(date), rate_nzd = rate / 100) %>% select(date, rate_nzd)

policy_zar <- read.csv("Data Model/zar_overnight_2018.csv") %>%
  mutate(date = as.Date(date), rate_zar = rate / 100) %>% select(date, rate_zar)

unwind_state <- processed_hmm_full[["NZD"]]$unwind_state
unwind_state

# rebuild emissions
yields_filled <- readRDS("yields_filled.rds")
country_data_full <- prepare_country(yields_filled[["NZD"]], nzd_zar_fx)[-1, ]
country_data_full$date <- as.Date(country_data_full$date)

all_dates <- country_data_full$date
model     <- processed_hmm_full[["NZD"]]$model
N         <- nrow(model$Mu)
T_full    <- length(all_dates)

X_full <- as.matrix(country_data_full[, c("d_y2y", "d_spread", "fx_return")])

B_full <- matrix(0, T_full, N)
for (s in seq_len(N))
  B_full[, s] <- dmvnorm(X_full, mean = model$Mu[s, ], sigma = model$Cov[,, s])


Z_full <- prepare_covariates(ted, vix, lagged = TRUE)
Z_filt <- merge(data.frame(date = all_dates), Z_full, by = "date", sort = TRUE)
Z_mat  <- as.matrix(Z_filt[, c("ted", "vix")])
A_full <- .compute_A(model$x, Z_mat, N)

Pi <- model$Pi

# rolling Viterbi from 2019
start_idx <- which(all_dates >= cutoff)[1]
rolling_viterbi <- rep(NA_integer_, T_full)

cat("Running rolling Viterbi", as.character(all_dates[start_idx]),
    "to", as.character(all_dates[T_full]), "\n")

for (i in start_idx:T_full) {
  v <- viterbi(B_full[1:i, , drop = FALSE],
               A_full[, , 1:i, drop = FALSE],
               Pi)
  rolling_viterbi[i] <- v$path[i]
  if (i %% 100 == 0) #just to be sure it works lol
    cat("  step", i - start_idx + 1, "of", T_full - start_idx + 1, "\n")
}

if (unwind_state == 2) {
  trans_to_unwind <- A_full[1, 2, ]
} else {
  trans_to_unwind <- A_full[2, 1, ]
}

state_df <- data.frame(date = all_dates, viterbi_rolling = rolling_viterbi,
                       trans_to_unwind = trans_to_unwind)

# build merged panel
panel <- state_df %>%
  inner_join(nzd_zar_fx, by = "date") %>%
  inner_join(policy_nzd, by = "date") %>%
  inner_join(policy_zar, by = "date") %>%
  arrange(date)

panel <- panel[panel$date >= cutoff, ]

panel$rate_nzd <- zoo::na.locf(panel$rate_nzd, na.rm = FALSE)
panel$rate_zar <- zoo::na.locf(panel$rate_zar, na.rm = FALSE)

# signal
panel$ema20 <- TTR::EMA(panel$trans_to_unwind, n = 20)
panel$skip_signal <- (panel$trans_to_unwind > panel$ema20) |
  (panel$viterbi_rolling == unwind_state)

# tweak for diff setups: 1 long, 0 flat, -1 short
if_signal <- 0
if_quiet  <- 1

n <- nrow(panel)
panel$position <- c(NA, ifelse(head(panel$skip_signal, -1), if_signal, if_quiet))

# detect exits: day t is the exit day if you held a different position at t-1
panel$is_exit_day <- FALSE
for (t in 2:n) {
  prev_pos <- panel$position[t - 1]
  curr_pos <- panel$position[t]
  if (!is.na(prev_pos) && !is.na(curr_pos) && 
      prev_pos != 0 && prev_pos != curr_pos) {
    panel$is_exit_day[t] <- TRUE
  }
}
if (!is.na(panel$position[n]) && panel$position[n] != 0) {
  panel$is_exit_day[n] <- TRUE
}

# walk forward: realize old position FX P&L on exit day, then new entry FX
# ps NA_reals are needed cus normal NA is logical type
panel$entry_fx    <- NA_real_
panel$exit_fx_pnl <- 0

current_entry_fx <- NA_real_
for (t in seq_len(n)) {
  pos      <- panel$position[t]
  prev_pos <- if (t == 1) 0 else (if (is.na(panel$position[t - 1])) 0 else panel$position[t - 1])
  
  # exit first: realize FX P&L at today's rate, sized by the position left
  if (panel$is_exit_day[t] && !is.na(current_entry_fx)) {
    panel$exit_fx_pnl[t] <- prev_pos * log(panel$rate[t] / current_entry_fx)
    current_entry_fx <- NA_real_
  }
  
  # entry: if today's position is nonzero and different from yesterday, today's FX is the new entry
  if (!is.na(pos) && pos != 0 && (is.na(prev_pos) || pos != prev_pos)) {
    current_entry_fx <- panel$rate[t]
  }
  
  panel$entry_fx[t] <- current_entry_fx
}

# returns
panel$carry_ret <- ifelse(is.na(panel$position) | panel$position == 0,
                          0,
                          panel$position * (panel$rate_zar - panel$rate_nzd) / 360)
panel$strategy_ret <- panel$carry_ret + panel$exit_fx_pnl

# benchmark: always long, daily-realized FX
panel$fx_ret_daily <- c(NA, diff(log(panel$rate)))
panel$bench_ret <- (panel$rate_zar - panel$rate_nzd) / 360 + panel$fx_ret_daily

bt <- panel[complete.cases(panel[, c("strategy_ret", "bench_ret")]), ]

# performance
perf <- function(r, label, pos = NULL, ann = 252) {
  r       <- r[is.finite(r)]
  ret_ann <- mean(r) * ann
  vol_ann <- sd(r)   * sqrt(ann)
  sharpe  <- ret_ann / vol_ann # ???
  cum     <- cumprod(1 + r)
  dd      <- (cum - cummax(cum)) / cummax(cum)
  max_dd  <- min(dd)
  win     <- mean(r > 0)
  
  cat(label, "\n")
  cat(sprintf("  Annualised return : %7.2f%%\n", ret_ann * 100))
  cat(sprintf("  Total return      : %7.2f%%\n", (prod(1 + r) - 1) * 100))
  cat(sprintf("  Annualised vol    : %7.2f%%\n", vol_ann * 100))
  cat(sprintf("  Sharpe            : %7.3f\n",   sharpe))
  cat(sprintf("  Max drawdown      : %7.2f%%\n", max_dd  * 100))
  cat(sprintf("  Win rate          : %7.2f%%\n", win     * 100))
  if (!is.null(pos)) {
    cat(sprintf("  Time long         : %7.2f%%\n", mean(pos == 1, na.rm = TRUE) * 100))
    cat(sprintf("  Time short        : %7.2f%%\n", mean(pos == -1, na.rm = TRUE) * 100))
  }
}

cat("OOS backtest period:", as.character(min(bt$date)),
    "to", as.character(max(bt$date)))
perf(bt$strategy_ret, "STRATEGY (deferred FX, daily carry)", pos = bt$position)
perf(bt$bench_ret,    "BENCHMARK (always long, daily FX)",   pos = rep(1, nrow(bt)))

# equity curves
plot(bt$date, cumprod(1 + bt$bench_ret),
     type = "l", col = "grey", lwd = 1.5, xlab = "Date", ylab = "Growth of 1",
     main = "NZD/ZAR Carry: strat vs always long")
lines(bt$date, cumprod(1 + bt$strategy_ret), col = "pink", lwd = 2)
legend("topleft",
       legend = c("Always long", "HMM-filtered, deferred FX"),
       col    = c("grey", "pink"),
       lwd    = c(1.5, 2), bty = "n")
grid()
