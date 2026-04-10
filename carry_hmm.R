# =============================================================================
# Carry Trade Regime-Switching Model
# TVTP-HMM: per-country, emissions = (slope, curvature, FX return)
#            transitions driven by VIX and TED spread
#
# Install dependencies once:
#   install.packages(c("mvtnorm", "ggplot2", "reshape2"))
# =============================================================================

library(mvtnorm)


# =============================================================================
# SECTION 1: DATA PREPARATION
# =============================================================================

#' Compute yield curve features and FX returns for one country
#'
#' @param yields   Data frame with columns: date, y3m, y2y, y5y, y10y
#'                 Yields in percent (e.g. 2.5 for 2.5%), daily frequency
#' @param fx       Data frame with columns: date, rate
#'                 FX as units of local currency per 1 USD
#'                 (higher value = local currency weaker)
#'
#' @return Data frame: date, slope, curvature, fx_return
#'   slope     = y10y - y3m             (carry quality / curve shape)
#'   curvature = y2y - 0.5*(y3m + y10y) (policy uncertainty / hump)
#'   fx_return = log(FX_t / FX_{t-1})   (positive = local ccy depreciated)

prepare_country <- function(yields, fx) {

  # Merge on date
  df <- merge(yields, fx, by = "date", sort = TRUE)

  # Yield curve features (already stationary — no differencing needed)
  df$slope     <- df$y10y - df$y3m
  df$curvature <- df$y2y - 0.5 * (df$y3m + df$y10y)

  # FX log return (stationary)
  df$fx_return <- c(NA, diff(log(df$rate)))

  # Drop first row (NA from diff) and any remaining NAs
  df <- df[complete.cases(df[, c("slope", "curvature", "fx_return")]), ]

  df[, c("date", "slope", "curvature", "fx_return")]
}


#' Prepare covariates (VIX and TED) aligned to a set of dates
#'
#' @param ted    Data frame: date, ted  (TED spread in percent)
#' @param vix    Data frame: date, vix  (VIX index level)
#' @param dates  Character or Date vector to align to
#'
#' @return Matrix (T x 2), standardised, columns: ted, vix

prepare_covariates <- function(ted, vix, dates) {

  dates <- as.character(dates)
  ted$date <- as.character(ted$date)
  vix$date <- as.character(vix$date)

  ted_vals <- ted$ted[match(dates, ted$date)]
  vix_vals <- vix$vix[match(dates, vix$date)]

  # Standardise: zero mean, unit variance
  Z <- cbind(
    ted = (ted_vals - mean(ted_vals, na.rm = TRUE)) / sd(ted_vals, na.rm = TRUE),
    vix = (vix_vals - mean(vix_vals, na.rm = TRUE)) / sd(vix_vals, na.rm = TRUE)
  )

  Z
}


# =============================================================================
# SECTION 2: CORE HMM BUILDING BLOCKS
# =============================================================================

# --- Split indices into K groups (used for initialisation) ---
.divide <- function(X, K, type = "sort") {
  n <- if (is.matrix(X)) nrow(X) else length(X)
  I <- if (type == "sort") order(if (is.matrix(X)) X[, 1] else X) else sample(n)
  lapply(seq_len(K), function(i) I[(floor((i-1)*n/K)+1):floor(i*n/K)])
}

# --- Compute time-varying transition matrix from TVTP parameters ---
# x: array [N x N x (1 + n_cov)]
#   x[i,j,1]   = intercept for row i, destination j
#   x[i,j,2:k] = covariate coefficients
# Diagonal of x is fixed to 0 (softmax reference for self-transition)
# Returns A: array [N x N x T]
.compute_A <- function(x, Z, N) {
  T     <- nrow(Z)
  n_cov <- ncol(Z)
  A     <- array(0, dim = c(N, N, T))
  for (t in seq_len(T)) {
    for (i in seq_len(N)) {
      lp   <- x[i,, 1] + as.vector(x[i,, 2:(n_cov+1), drop=FALSE] %*% Z[t,])
      lp   <- lp - max(lp)
      A[i,, t] <- exp(lp) / sum(exp(lp))
    }
  }
  A
}

# --- Negative expected complete-data log-likelihood for TVTP ---
.obj_tvtp <- function(x_vec, T, Z, Xi, N) {
  n_cov <- ncol(Z)
  x     <- array(x_vec, dim = c(N, N, n_cov + 1))
  Xi_p  <- aperm(Xi, c(2, 3, 1))   # N x N x (T-1)
  total <- 0
  for (t in seq_len(T - 1)) {
    for (i in seq_len(N)) {
      lp  <- x[i,, 1] + as.vector(x[i,, 2:(n_cov+1), drop=FALSE] %*% Z[t,])
      lp  <- lp - max(lp)
      lsm <- lp - log(sum(exp(lp)))
      total <- total + sum(Xi_p[i,, t] * lsm)
    }
  }
  -total
}


# =============================================================================
# SECTION 3: FIT TVTP-HMM FOR ONE COUNTRY
# =============================================================================

#' Fit TVTP-HMM to one country's data
#'
#' @param X     Matrix (T x 3): slope, curvature, fx_return
#' @param Z     Matrix (T x 2): standardised TED, VIX
#' @param N     Number of hidden states (2 or 3 recommended)
#' @param cyc   Maximum EM iterations
#' @param tol   Convergence tolerance on log-likelihood change
#' @param seed  RNG seed for reproducibility
#'
#' @return List:
#'   $Mu     (N x 3)       State mean vectors
#'   $Cov    (3 x 3 x N)   State covariance matrices
#'   $x      (N x N x 3)   TVTP parameters [intercept, TED coef, VIX coef]
#'   $Pi     (N)            Initial state distribution
#'   $LL     numeric        Log-likelihood history

fit_country_hmm <- function(X, Z, N = 2, cyc = 100, tol = 1e-4, seed = 42) {

  set.seed(seed)
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(Z)) Z <- as.matrix(Z)
  T     <- nrow(X)
  p     <- ncol(X)
  n_cov <- ncol(Z)
  stopifnot(nrow(Z) == T, T > N * 10)

  # --- Initialise emission parameters ---
  idx <- .divide(X, N, "sort")
  Mu  <- matrix(0, N, p)
  Cov <- array(0, dim = c(p, p, N))
  for (i in seq_len(N)) {
    sub      <- X[idx[[i]], , drop = FALSE]
    Cov[,,i] <- cov(sub) + diag(1e-4, p)
    Mu[i, ]  <- colMeans(sub)
  }

  # --- Initialise TVTP parameters (uniform transitions) ---
  Pi <- rep(1/N, N)
  x  <- array(0, dim = c(N, N, n_cov + 1))
  A  <- .compute_A(x, Z, N)

  lik <- -Inf
  LL  <- numeric(cyc)

  # --- EM loop ---
  for (cycle in seq_len(cyc)) {

    # Emission probabilities B: T x N
    B <- matrix(0, T, N)
    for (i in seq_len(N)) {
      B[, i] <- dmvnorm(X, mean = Mu[i, ], sigma = Cov[,, i])
    }
    B[B == 0] <- 1e-300

    # Forward pass
    alpha    <- matrix(0, T, N)
    sc       <- numeric(T)
    alpha[1,] <- Pi * B[1,]
    sc[1]     <- sum(alpha[1,]); if (sc[1] == 0) sc[1] <- 1e-300
    alpha[1,] <- alpha[1,] / sc[1]
    for (t in 2:T) {
      for (i in seq_len(N))
        alpha[t, i] <- sum(alpha[t-1,] * A[, i, t]) * B[t, i]
      sc[t]     <- sum(alpha[t,]); if (sc[t] == 0) sc[t] <- 1e-300
      alpha[t,] <- alpha[t,] / sc[t]
    }
    log_lik <- sum(log(sc))

    # Backward pass
    beta     <- matrix(0, T, N)
    beta[T,] <- 1 / sc[T]
    for (t in (T-1):1)
      for (i in seq_len(N))
        beta[t, i] <- sum(A[i,, t] * B[t+1,] * beta[t+1,]) / sc[t]

    # E-step: Xi [T-1 x N x N] and Gamma [T-1 x N]
    Xi <- array(0, dim = c(T-1, N, N))
    for (t in seq_len(T-1)) {
      for (i in seq_len(N))
        for (j in seq_len(N))
          Xi[t, i, j] <- alpha[t, i] * A[i, j, t] * B[t+1, j] * beta[t+1, j]
      s <- sum(Xi[t,,]); if (s > 0) Xi[t,,] <- Xi[t,,] / s
    }
    Gamma <- apply(Xi, c(1, 2), sum)

    # M-step: Pi
    Pi <- pmax(Gamma[1,] / sum(Gamma[1,]), 1e-10)

    # M-step: TVTP via L-BFGS-B
    # Fix diagonal of x to zero (softmax reference category)
    lb <- rep(-Inf, N * N * (n_cov + 1))
    ub <- rep( Inf, N * N * (n_cov + 1))
    for (i in seq_len(N))
      for (k in seq_len(n_cov + 1)) {
        idx_flat <- (k-1)*N*N + (i-1)*N + i
        lb[idx_flat] <- 0; ub[idx_flat] <- 0
      }

    res <- optim(
      par     = as.vector(x) + rnorm(length(x), 0, 0.05),
      fn      = .obj_tvtp,
      T = T, Z = Z, Xi = Xi, N = N,
      method  = "L-BFGS-B",
      lower   = lb, upper = ub,
      control = list(maxit = 300, factr = 1e8)
    )
    x <- array(res$par, dim = c(N, N, n_cov + 1))
    A <- .compute_A(x, Z, N)

    # M-step: Mu and Cov
    for (i in seq_len(N)) {
      w     <- Gamma[, i]; ws <- sum(w)
      if (ws < 1e-10) next
      Mu[i,]   <- colSums(X[1:(T-1),, drop=FALSE] * w) / ws
      d        <- sweep(X[1:(T-1),, drop=FALSE], 2, Mu[i,])
      Cov[,,i] <- (t(d * w) %*% d) / ws + diag(1e-4, p)
    }

    oldlik    <- lik
    lik       <- log_lik
    LL[cycle] <- lik

    if (cycle > 2 && abs(lik - oldlik) < tol * abs(oldlik)) break
  }

  list(Mu = Mu, Cov = Cov, x = x, Pi = Pi, LL = LL[seq_len(cycle)])
}


# =============================================================================
# SECTION 4: DECODE — VITERBI + SMOOTHED PROBABILITIES
# =============================================================================

#' Decode regime sequence for one country
#'
#' @param X      Matrix (T x 3): same features used in fitting
#' @param Z      Matrix (T x 2): TED and VIX covariates
#' @param model  Output from fit_country_hmm()
#'
#' @return List:
#'   $states  Integer vector (T): most likely state at each time (Viterbi)
#'   $prb     Matrix (N x T): smoothed state probabilities

decode_country <- function(X, Z, model) {

  if (!is.matrix(X)) X <- as.matrix(X)
  T <- nrow(X)
  N <- nrow(model$Mu)
  A <- .compute_A(model$x, Z, N)

  # Emission probabilities
  B <- matrix(0, N, T)
  for (i in seq_len(N))
    B[i,] <- dmvnorm(X, mean = model$Mu[i,], sigma = model$Cov[,,i])
  B[B == 0] <- 1e-300

  # Log-Viterbi
  delta <- matrix(-Inf, N, T)
  psi   <- matrix(0L,   N, T)
  prb   <- matrix(0,    N, T)
  sc    <- numeric(T)

  delta[, 1] <- log(model$Pi) + log(B[, 1])
  prb[, 1]   <- model$Pi * B[, 1]
  sc[1]      <- sum(prb[, 1]); prb[, 1] <- prb[, 1] / sc[1]

  for (t in 2:T) {
    for (j in seq_len(N)) {
      v           <- delta[, t-1] + log(pmax(A[, j, t], 1e-300))
      delta[j, t] <- max(v) + log(B[j, t])
      psi[j, t]   <- which.max(v)
      prb[j, t]   <- max(prb[, t-1] * A[, j, t]) * B[j, t]
    }
    sc[t] <- sum(prb[, t])
    if (sc[t] > 0) prb[, t] <- prb[, t] / sc[t]
  }

  states    <- integer(T)
  states[T] <- which.max(delta[, T])
  for (t in (T-1):1) states[t] <- psi[states[t+1], t]

  list(states = states, prb = prb)
}


# =============================================================================
# SECTION 5: IDENTIFY UNWIND STATE PER COUNTRY
# =============================================================================

#' Identify which state index corresponds to the carry unwind regime
#'
#' Unwind state = the state with the most negative mean FX return
#' (local currency depreciating, i.e. carry crash)
#'
#' FX return is assumed to be column 3 of the emission vector
#'
#' @param model  Output from fit_country_hmm()
#' @return Integer: index of the unwind state

identify_unwind_state <- function(model) {
  which.min(model$Mu[, 3])   # column 3 = fx_return; most negative = unwind
}


# =============================================================================
# SECTION 6: FIT ALL COUNTRIES
# =============================================================================

#' Fit TVTP-HMM for all countries and decode regimes
#'
#' @param country_data  Named list of outputs from prepare_country()
#'                      Names should be currency codes e.g. "AUD", "JPY"
#' @param ted           Data frame: date, ted
#' @param vix           Data frame: date, vix
#' @param N             Number of hidden states per country
#' @param cyc           Max EM iterations
#' @param tol           Convergence tolerance
#'
#' @return List per country, each containing:
#'   $model         Output from fit_country_hmm()
#'   $decoded       Output from decode_country()
#'   $unwind_state  Integer index of the unwind state
#'   $unwind_prob   Numeric vector (T): probability of being in unwind state
#'   $dates         Date vector

fit_all_countries <- function(country_data, ted, vix,
                               N = 2, cyc = 100, tol = 1e-4) {

  results <- vector("list", length(country_data))
  names(results) <- names(country_data)

  for (cc in names(country_data)) {
    cat(sprintf("\n--- Fitting %s ---\n", cc))

    cd    <- country_data[[cc]]
    dates <- as.character(cd$date)

    X <- as.matrix(cd[, c("slope", "curvature", "fx_return")])
    Z <- prepare_covariates(ted, vix, dates)

    # Drop rows where Z has NAs (missing VIX or TED on some dates)
    valid <- complete.cases(Z)
    X     <- X[valid, ]
    Z     <- Z[valid, ]
    dates <- dates[valid]

    model   <- fit_country_hmm(X, Z, N = N, cyc = cyc, tol = tol)
    decoded <- decode_country(X, Z, model)
    u_state <- identify_unwind_state(model)

    results[[cc]] <- list(
      model        = model,
      decoded      = decoded,
      unwind_state = u_state,
      unwind_prob  = decoded$prb[u_state, ],
      dates        = dates
    )

    cat(sprintf("  Unwind state: %d  |  Mean unwind prob: %.3f\n",
                u_state, mean(decoded$prb[u_state, ])))
    cat(sprintf("  Unwind state means -> slope: %.3f  curv: %.3f  fx_ret: %.4f\n",
                model$Mu[u_state, 1], model$Mu[u_state, 2], model$Mu[u_state, 3]))
  }

  results
}


# =============================================================================
# SECTION 7: SIGNAL CONSTRUCTION
# =============================================================================

#' Build daily carry trade signal from regime probabilities
#'
#' At each date: rank countries by unwind probability.
#' Short carry (or avoid) on high-unwind-probability countries.
#' Long carry on low-unwind-probability countries.
#'
#' Position entry uses previous day's unwind probability to avoid look-ahead.
#'
#' @param hmm_results     Output from fit_all_countries()
#' @param fx_returns      Named list of numeric vectors, daily log FX returns
#'                        per country, aligned to same dates
#' @param rate_diffs      Named list of numeric vectors, daily rate differential
#'                        (country 3m rate - USD 3m rate) / 252
#'                        These are the carry component
#' @param n_long          How many currencies to go long carry
#' @param n_short         How many currencies to short carry
#' @param prob_threshold  Minimum probability gap to take a position
#'                        Countries with unwind_prob between threshold and
#'                        (1-threshold) are excluded (ambiguous regime)
#'
#' @return List:
#'   $positions       Data frame (T x n_countries): position weights
#'   $carry_pnl       Data frame (T x n_countries): carry component PnL
#'   $fx_pnl          Data frame (T x n_countries): FX component PnL
#'   $total_pnl       Numeric vector (T): strategy daily return
#'   $dates           Date vector

build_signal <- function(hmm_results, fx_returns, rate_diffs,
                          n_long = 3, n_short = 3,
                          prob_threshold = 0.6) {

  countries <- names(hmm_results)

  # Find common dates across all countries
  all_dates <- Reduce(intersect, lapply(hmm_results, function(r) r$dates))
  all_dates <- sort(all_dates)
  T         <- length(all_dates)

  # Align unwind probabilities to common dates
  unwind_mat <- matrix(NA, T, length(countries),
                       dimnames = list(all_dates, countries))
  for (cc in countries) {
    r   <- hmm_results[[cc]]
    idx <- match(all_dates, r$dates)
    unwind_mat[, cc] <- r$unwind_prob[idx]
  }

  # Align FX returns and rate differentials
  fx_mat   <- matrix(NA, T, length(countries), dimnames = list(all_dates, countries))
  rate_mat <- matrix(NA, T, length(countries), dimnames = list(all_dates, countries))
  for (cc in countries) {
    fx_mat[, cc]   <- fx_returns[[cc]][match(all_dates, names(fx_returns[[cc]]))]
    rate_mat[, cc] <- rate_diffs[[cc]][match(all_dates, names(rate_diffs[[cc]]))]
  }

  # Build positions day by day (use lagged unwind_prob to avoid look-ahead)
  pos_mat <- matrix(0, T, length(countries), dimnames = list(all_dates, countries))

  for (t in 2:T) {
    prob_yesterday <- unwind_mat[t - 1, ]
    if (any(is.na(prob_yesterday))) next

    # Rank by unwind probability (ascending = safer carry)
    ranked <- order(prob_yesterday)

    # Long carry: lowest unwind probability, below threshold
    long_candidates  <- ranked[prob_yesterday[ranked] < (1 - prob_threshold)]
    short_candidates <- ranked[prob_yesterday[ranked] >       prob_threshold]

    n_l <- min(n_long,  length(long_candidates))
    n_s <- min(n_short, length(short_candidates))

    if (n_l > 0) pos_mat[t, long_candidates[1:n_l]]  <-  1 / n_l
    if (n_s > 0) pos_mat[t, short_candidates[(length(short_candidates) - n_s + 1):
                                               length(short_candidates)]] <- -1 / n_s
  }

  # PnL: carry trade return = FX return + carry differential
  # Note: position +1 = long foreign currency (short USD)
  #       FX return > 0 means foreign currency depreciated -> loss for long
  #       So correct sign: long position profits when FX return < 0
  #       (foreign currency appreciated against USD)
  carry_pnl_mat <- pos_mat *  rate_mat          # interest rate differential
  fx_pnl_mat    <- pos_mat * (-fx_mat)          # FX return (negative: long profits on appreciation)
  total_mat     <- carry_pnl_mat + fx_pnl_mat

  total_pnl <- rowSums(total_mat, na.rm = TRUE)

  list(
    positions  = as.data.frame(pos_mat),
    carry_pnl  = as.data.frame(carry_pnl_mat),
    fx_pnl     = as.data.frame(fx_pnl_mat),
    total_pnl  = total_pnl,
    dates      = as.Date(all_dates)
  )
}


# =============================================================================
# SECTION 8: PERFORMANCE EVALUATION
# =============================================================================

#' Compute standard performance metrics for a daily return series
#'
#' @param daily_ret  Numeric vector of daily returns
#' @param dates      Date vector (optional, for display)
#' @param ann        Trading days per year (default 252)
#'
#' @return Named list of metrics (also prints a summary)

evaluate_strategy <- function(daily_ret, dates = NULL, ann = 252) {

  r <- daily_ret[is.finite(daily_ret) & !is.na(daily_ret)]

  ret_ann  <- mean(r) * ann
  vol_ann  <- sd(r) * sqrt(ann)
  sharpe   <- ret_ann / vol_ann
  cum      <- cumprod(1 + r)
  dd       <- (cum - cummax(cum)) / cummax(cum)
  max_dd   <- min(dd)
  calmar   <- ret_ann / abs(max_dd)
  win_rate <- mean(r > 0)
  skew     <- mean(((r - mean(r)) / sd(r))^3)

  cat("=========================================\n")
  cat(sprintf("  Annualised Return  : %7.2f%%\n", ret_ann  * 100))
  cat(sprintf("  Annualised Vol     : %7.2f%%\n", vol_ann  * 100))
  cat(sprintf("  Sharpe Ratio       : %7.3f\n",   sharpe))
  cat(sprintf("  Max Drawdown       : %7.2f%%\n", max_dd   * 100))
  cat(sprintf("  Calmar Ratio       : %7.3f\n",   calmar))
  cat(sprintf("  Win Rate           : %7.2f%%\n", win_rate * 100))
  cat(sprintf("  Return Skewness    : %7.3f\n",   skew))
  cat("=========================================\n")

  invisible(list(ret_ann = ret_ann, vol_ann = vol_ann, sharpe = sharpe,
                 max_dd = max_dd, calmar = calmar, win_rate = win_rate,
                 skew = skew))
}


# =============================================================================
# SECTION 9: DIAGNOSTICS — INSPECT WHAT STATES MEAN
# =============================================================================

#' Print a readable summary of fitted states for each country
#'
#' @param hmm_results  Output from fit_all_countries()

summarise_states <- function(hmm_results) {

  cat("\n=== STATE SUMMARY ===\n")
  for (cc in names(hmm_results)) {
    r  <- hmm_results[[cc]]
    Mu <- r$model$Mu
    x  <- r$model$x      # TVTP parameters
    N  <- nrow(Mu)
    u  <- r$unwind_state

    cat(sprintf("\n%s (unwind state = %d)\n", cc, u))
    cat(sprintf("  %-8s  %8s  %9s  %9s\n", "State", "Slope", "Curvature", "FX_ret"))
    for (s in seq_len(N)) {
      marker <- if (s == u) " <-- UNWIND" else ""
      cat(sprintf("  State %d   %8.3f  %9.3f  %9.5f%s\n",
                  s, Mu[s, 1], Mu[s, 2], Mu[s, 3], marker))
    }

    # TVTP sensitivity to VIX and TED
    # x[i,j,2] = TED coefficient for transition i->j
    # x[i,j,3] = VIX coefficient for transition i->j
    cat("  Transition into unwind state:\n")
    for (i in seq_len(N)) {
      if (i == u) next
      cat(sprintf("    From state %d -> unwind: intercept=%+.3f  TED=%+.3f  VIX=%+.3f\n",
                  i, x[i, u, 1], x[i, u, 2], x[i, u, 3]))
    }
  }
}


# =============================================================================
# SECTION 10: PLOTTING
# =============================================================================

#' Plot regime probabilities and cumulative strategy return
#'
#' @param signal      Output from build_signal()
#' @param hmm_results Output from fit_all_countries()

plot_results <- function(signal, hmm_results) {

  countries <- names(hmm_results)
  n         <- length(countries)
  dates     <- signal$dates

  # --- Plot 1: Unwind probabilities across countries ---
  old_par <- par(mfrow = c(ceiling((n + 1) / 2), 2),
                 mar   = c(2, 3, 2, 1))
  on.exit(par(old_par))

  for (cc in countries) {
    r   <- hmm_results[[cc]]
    idx <- match(as.character(dates), r$dates)
    up  <- r$unwind_prob[idx]

    plot(dates, up, type = "l", lwd = 1.2, col = "#D6604D",
         ylim = c(0, 1), main = cc,
         xlab = "", ylab = "P(unwind)", cex.main = 1)
    abline(h = 0.6, lty = 2, col = "grey50")
    grid()
  }

  # --- Plot 2: Cumulative strategy return ---
  dev.new()
  cum_ret <- cumprod(1 + ifelse(is.na(signal$total_pnl), 0, signal$total_pnl))
  plot(dates, cum_ret, type = "l", lwd = 2, col = "#2166AC",
       main = "Carry Strategy — Cumulative Return",
       xlab = "Date", ylab = "Growth of $1")
  abline(h = 1, lty = 2, col = "grey50")
  grid()
}


# =============================================================================
# SECTION 11: COMPLETE USAGE EXAMPLE
# =============================================================================
#
# Step 1 — Prepare your data
#
# Each yield CSV needs columns: date, y3m, y2y, y5y, y10y  (yields in %)
# Each FX CSV needs columns:    date, rate  (local per USD)
# TED CSV needs columns:        date, ted   (TED spread in %)
# VIX CSV needs columns:        date, vix   (VIX index level)
#
#   countries <- c("AUD", "EUR", "GBP", "JPY", "CAD", "NZD", "SEK", "NOK")
#
#   country_data <- lapply(setNames(countries, countries), function(cc) {
#     yields <- read.csv(paste0("data/yields_", cc, ".csv"))
#     fx     <- read.csv(paste0("data/fx_",     cc, ".csv"))
#     prepare_country(yields, fx)
#   })
#
#   ted <- read.csv("data/ted_spread.csv")
#   vix <- read.csv("data/vix.csv")
#
# Step 2 — Fit models
#
#   hmm_results <- fit_all_countries(country_data, ted, vix, N = 2, cyc = 100)
#
# Step 3 — Inspect what states mean
#
#   summarise_states(hmm_results)
#
# Step 4 — Check identification: confirm unwind state is correctly labelled
#   For each country, the unwind state should have:
#     - negative fx_return mean (currency depreciated)
#     - flatter or negative slope mean (curve flattening)
#   If a country looks wrong, override manually:
#     hmm_results[["JPY"]]$unwind_state <- 1
#     hmm_results[["JPY"]]$unwind_prob  <- hmm_results[["JPY"]]$decoded$prb[1, ]
#
# Step 5 — Build carry signal
#   (You need fx_returns and rate_diffs as named lists of named numeric vectors)
#
#   signal <- build_signal(
#     hmm_results    = hmm_results,
#     fx_returns     = fx_returns,     # named list, each element named by date
#     rate_diffs     = rate_diffs,     # named list, each element named by date
#     n_long         = 3,
#     n_short        = 3,
#     prob_threshold = 0.6
#   )
#
# Step 6 — Evaluate
#
#   evaluate_strategy(signal$total_pnl)
#
# Step 7 — Plot
#
#   plot_results(signal, hmm_results)
