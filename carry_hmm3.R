# =============================================================================
# Carry Trade Regime-Switching Model
# TVTP-HMM: per-country, emissions = (d_y2y, d_spread, fx_return)
#   d_y2y    = first difference of 2y yield
#   d_spread = first difference of (10y - 2y) spread
#   Both are stationary by construction.
# Transitions driven by VIX and TED spread.
# =============================================================================

library(mvtnorm)


# =============================================================================
# SECTION 1: DATA PREPARATION
# =============================================================================

#' Compute yield curve features and FX returns for one country
#'
#' @param yields   Data frame with columns: date, X2Y, X10Y
#'                 Yields in percent (e.g. 2.5 for 2.5%), daily frequency
#' @param fx       Data frame with columns: date, rate
#'                 FX as units of local currency per 1 USD
#'
#' @return Data frame: date, d_y2y, d_spread, fx_return
#'   d_y2y    = diff(y2y)          (change in 2y yield)
#'   d_spread = diff(y10y - y2y)   (change in 10y-2y spread)
#'   fx_return = log(FX_t/FX_{t-1}) (positive = local ccy depreciated)

prepare_covariates <- function(ted, vix, lagged = TRUE) {
  
  y_df <- merge(ted, vix, by = "date") %>%
    mutate(date = as.Date(date, format = "%Y-%m-%d")) %>%
    arrange(date) %>%
    mutate(
      ted_lag = lag(ted),
      vix_lag = lag(vix)
    )
  
  if (lagged) {
    out <- y_df[, c("date", "ted_lag", "vix_lag")]
    colnames(out) <- c("date", "ted", "vix")
  } else {
    out <- y_df[, c("date", "ted", "vix")]
  }
  
  out
}


#' Prepare covariates (VIX and TED) aligned to a set of dates
#'
#' @param ted    Data frame: date, ted  (TED spread in percent)
#' @param vix    Data frame: date, vix  (VIX index level)
#' @param dates  Character or Date vector to align to
#'
#' @return Matrix (T x 2), standardised, columns: ted, vix



# =============================================================================
# SECTION 2: CORE HMM BUILDING BLOCKS
# =============================================================================

.divide <- function(X, K, type = "sort") {
  n <- if (is.matrix(X)) nrow(X) else length(X)
  I <- if (type == "sort") order(if (is.matrix(X)) X[, 1] else X) else sample(n)
  lapply(seq_len(K), function(i) I[(floor((i-1)*n/K)+1):floor(i*n/K)])
}

.compute_A <- function(x, Z, N) {
  T     <- nrow(Z)
  n_cov <- ncol(Z)
  A     <- array(0, dim = c(N, N, T))
  for (t in seq_len(T)) {
    for (i in seq_len(N)) {
      coef_mat <- matrix(x[i,, 2:(n_cov+1)], nrow = N, ncol = n_cov)
      lp       <- x[i,, 1] + as.vector(coef_mat %*% matrix(Z[t,], ncol = 1))
      lp      <- lp - max(lp)
      A[i,, t] <- exp(lp) / sum(exp(lp))
    }
  }
  A
}

# Negative expected complete-data log-likelihood for TVTP
.obj_tvtp <- function(x_vec, T, Z, Xi, N) {
  n_cov <- ncol(Z)
  x     <- array(x_vec, dim = c(N, N, n_cov + 1))
  Xi_p  <- aperm(Xi, c(2, 3, 1))   # N x N x (T-1)
  total <- 0
  for (t in seq_len(T - 1)) {
    for (i in seq_len(N)) {
      coef_mat <- matrix(x[i,, 2:(n_cov+1)], nrow = N, ncol = n_cov)
      lp       <- x[i,, 1] + as.vector(coef_mat %*% matrix(Z[t,], ncol = 1))
      lp       <- lp - max(lp)
      lsm      <- lp - log(sum(exp(lp)))
      contrib  <- Xi_p[i,, t] * lsm
      contrib[!is.finite(contrib)] <- 0
      total    <- total + sum(contrib)
    }
  }
  if (!is.finite(total)) return(1e10)
  -total
}

#joint most likely path calculation
viterbi <- function(B, A, Pi) {
  T      <- nrow(B); N <- ncol(B)
  delta  <- matrix(-Inf, T, N)
  psi    <- matrix(0L, T, N)
  
  delta[1, ] <- log(Pi) + log(B[1, ])
  
  for (t in 2:T) {
    for (j in seq_len(N)) {
      vals       <- delta[t-1, ] + log(A[, j, t])
      psi[t, j]  <- which.max(vals)
      delta[t, j] <- max(vals) + log(B[t, j])
    }
  }
  
  # Decode path
  path    <- integer(T)
  path[T] <- which.max(delta[T, ])
  for (t in (T-1):1) path[t] <- psi[t+1, path[t+1]]
  
  # Extract log-probability of the decoded path at each t
  path_log_probs <- delta[cbind(seq_len(T), path)]
  
  # Normalise across states at each t to get P(s_t = path[t] | Viterbi context)
  # This is exp(delta[t, path[t]]) / sum(exp(delta[t, ]))
  log_denom      <- apply(delta, 1, function(row) {
    m <- max(row); m + log(sum(exp(row - m)))   # log-sum-exp
  })
  path_probs     <- exp(path_log_probs - log_denom)
  
  list(
    path            = path,
    path_log_probs  = path_log_probs,   # log joint prob of best path through t
    path_probs      = path_probs,       # normalised: share of total prob mass on best path
    delta           = delta             # full delta matrix if needed
  )
}


# =============================================================================
# SECTION 3: FIT TVTP-HMM FOR ONE COUNTRY
# =============================================================================

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

fit_country_hmm <- function(X, Z, N = 2, cyc = 100, tol = 1e-4, seed = 42, compute_hessian = TRUE) {
  
  set.seed(seed)
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(Z)) Z <- as.matrix(Z)
  T     <- nrow(X)
  p     <- ncol(X)
  n_cov <- ncol(Z)
  stopifnot(nrow(Z) == T, T > N * 10)
  
  # Initialise emission parameters by sorting on d_y2y (col 1)
  idx <- .divide(X, N, "sort")
  Mu  <- matrix(0, N, p)
  Cov <- array(0, dim = c(p, p, N))
  for (i in seq_len(N)) {
    sub      <- X[idx[[i]], , drop = FALSE]
    Cov[,,i] <- cov(sub) + diag(1e-4, p)
    Mu[i, ]  <- colMeans(sub)
  }
  
  Pi <- rep(1/N, N)
  x  <- array(0, dim = c(N, N, n_cov + 1))
  A  <- .compute_A(x, Z, N)
  
  lik <- -Inf
  LL  <- numeric(cyc)
  
  for (cycle in seq_len(cyc)) {
    cat("processing cycle:",cycle)
    # Emission probabilities B: T x N
    B <- matrix(0, T, N)
    for (i in seq_len(N))
      B[, i] <- dmvnorm(X, mean = Mu[i, ], sigma = Cov[,, i])
    B[B == 0] <- 1e-300
    
    # Forward pass
    alpha     <- matrix(0, T, N)
    sc        <- numeric(T)
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
    beta[T,] <- 1 #Changed from beta[T,] <- 1 / sc[T]
    for (t in (T-1):1)
      for (i in seq_len(N))
        beta[t, i] <- sum(A[i,, t] * B[t+1,] * beta[t+1,]) / sc[t]
    
    # E-step
    Xi <- array(0, dim = c(T-1, N, N))
    for (t in seq_len(T-1)) {
      for (i in seq_len(N))
        for (j in seq_len(N))
          Xi[t, i, j] <- alpha[t, i] * A[i, j, t] * B[t+1, j] * beta[t+1, j]
      s <- sum(Xi[t,,]); if (s > 0) Xi[t,,] <- Xi[t,,] / s
    }
    Gamma <- apply(Xi, c(1, 2), sum)
    Gamma_full <- rbind(Gamma, alpha[T, ])
    
    # M-step: Pi
    Pi <- pmax(Gamma[1,] / sum(Gamma[1,]), 1e-10)
    
    # M-step: TVTP via L-BFGS-B
    # Only optimise off-diagonal elements; diagonal fixed at 0 (softmax reference)
    # Build index mapping between full array and free parameters
    free_idx <- which(array(sapply(seq_len(N*N), function(k) (k-1)%%N != (k-1)%/%N),
                            dim = c(N, N)) , arr.ind = FALSE)
    x_full <- x
    
    obj_free <- function(par) {
      x_try <- x_full
      for (k in seq_len(n_cov + 1))
        x_try[,, k][free_idx + (k-1)*N*N - (k-1)*N*N] <- par[((k-1)*length(free_idx)+1):(k*length(free_idx))]
      # rebuild properly
      x_try <- array(0, dim = c(N, N, n_cov + 1))
      par_idx <- 1
      for (k in seq_len(n_cov + 1))
        for (ij in free_idx) {
          i <- ((ij-1) %% N) + 1
          j <- ((ij-1) %/% N) + 1
          x_try[i, j, k] <- par[par_idx]
          par_idx <- par_idx + 1
        }
      .obj_tvtp(as.vector(x_try), T, Z, Xi, N)
    }
    
    # Extract current free parameters as starting values
    par0 <- numeric(length(free_idx) * (n_cov + 1))
    par_idx <- 1
    for (k in seq_len(n_cov + 1))
      for (ij in free_idx) {
        i <- ((ij-1) %% N) + 1
        j <- ((ij-1) %/% N) + 1
        par0[par_idx] <- x[i, j, k]
        par_idx <- par_idx + 1
      }
    
    res <- optim(
      par     = par0,
      fn      = obj_free,
      method  = "L-BFGS-B",
      control = list(maxit = 300, factr = 1e8)
    )
    
    # Reconstruct full x array from optimised free parameters
    x <- array(0, dim = c(N, N, n_cov + 1))
    par_idx <- 1
    for (k in seq_len(n_cov + 1))
      for (ij in free_idx) {
        i <- ((ij-1) %% N) + 1
        j <- ((ij-1) %/% N) + 1
        x[i, j, k] <- res$par[par_idx]
        par_idx <- par_idx + 1
      }
    A <- .compute_A(x, Z, N)
    
    # M-step: Mu and Cov
    for (i in seq_len(N)) {
      w  <- Gamma_full[, i]; ws <- sum(w)
      if (ws < 1e-10) next
      Mu[i,]   <- colSums(X[1:(T),, drop=FALSE] * w) / ws #T-1 change to T
      d        <- sweep(X[1:(T),, drop=FALSE], 2, Mu[i,]) #T-1 change to T
      Cov[,,i] <- (t(d * w) %*% d) / ws + diag(1e-4, p)
    }
    
    oldlik    <- lik
    lik       <- log_lik
    LL[cycle] <- lik
    
    if (cycle > 2 && abs(lik - oldlik) < tol * abs(oldlik)) break
  }
  
  
  
  hess <- NULL
  se   <- NULL
  if (compute_hessian) {
    # obj_free is no longer in scope after the loop, so rebuild it
    Xi_final   <- Xi   # Xi from the last E-step
    x_final    <- x
    obj_final  <- function(par) {
      x_try   <- array(0, dim = c(N, N, n_cov + 1))
      par_idx <- 1
      for (k in seq_len(n_cov + 1))
        for (ij in free_idx) {
          i_      <- ((ij - 1) %% N) + 1
          j_      <- ((ij - 1) %/% N) + 1
          x_try[i_, j_, k] <- par[par_idx]
          par_idx <- par_idx + 1
        }
      .obj_tvtp(as.vector(x_try), T, Z, Xi_final, N)
    }
    
    final_par <- numeric(length(free_idx) * (n_cov + 1))
    par_idx   <- 1
    for (k in seq_len(n_cov + 1))
      for (ij in free_idx) {
        i_ <- ((ij - 1) %% N) + 1
        j_ <- ((ij - 1) %/% N) + 1
        final_par[par_idx] <- x[i_, j_, k]
        par_idx <- par_idx + 1
      }
    
    hess <- numDeriv::hessian(obj_final, final_par)
    vcov <- tryCatch(solve(hess), error = function(e) {
      warning("Hessian singular — returning MASS::ginv instead"); MASS::ginv(hess)
    })
    se <- sqrt(pmax(diag(vcov), 0))   # pmax guards against tiny negatives
  }
  
  
  
  if (!is.null(se)) {
    param_labels <- character(length(free_idx) * (n_cov + 1))
    cov_names    <- c("intercept", colnames(Z))
    idx_lab      <- 1
    for (k in seq_len(n_cov + 1))
      for (ij in free_idx) {
        i_ <- ((ij - 1) %% N) + 1
        j_ <- ((ij - 1) %/% N) + 1
        param_labels[idx_lab] <- paste0(cov_names[k], "_", i_, "->", j_)
        idx_lab <- idx_lab + 1
      }
    wald_table <- data.frame(
      param  = param_labels,
      est    = final_par,
      se     = se,
      z      = final_par / se,
      pval   = 2 * pnorm(abs(final_par / se), lower.tail = FALSE)
    )
  }
  
  
  smoothed_probs <- Gamma_full   # T x N matrix, P(s_t = i | y_T, x_T)
  state_sequence <- apply(smoothed_probs, 1, which.max) #most likely pointwise
  viterbi_path <- viterbi(B, A, Pi)
  
  list(
    Mu             = Mu,
    Cov            = Cov,
    x              = x,
    Pi             = Pi,
    LL             = LL[seq_len(cycle)],
    smoothed_probs = smoothed_probs,    # T x N
    state_sequence = state_sequence,    # T-length pointwise MAP
    viterbi_path   = viterbi_path,      # T-length Viterbi path
    hessian        = hess,              # NULL unless compute_hessian = TRUE
    vcov           = if (!is.null(hess)) vcov else NULL,
    wald_table     = if (!is.null(se)) wald_table else NULL
  )
}


# =============================================================================
# SECTION 4: DECODE
# =============================================================================

#' Decode regime sequence for one country
#'
#' @param X      Matrix (T x 3): d_y2y, d_spread, fx_return
#' @param Z      Matrix (T x 2): TED, VIX covariates
#' @param model  Output from fit_country_hmm()
#'
#' @return List: states (integer T), prb (N x T)


#redundant


# =============================================================================
# SECTION 5: IDENTIFY UNWIND STATE
# =============================================================================

#' Unwind state = state with most negative mean FX return (col 3)
#' Verify this manually via summarise_states() — override if needed.


#redundant
identify_unwind_state <- function(model) {
  which.max(model$Cov[3,3])
}


# =============================================================================
# SECTION 6: FIT ALL COUNTRIES
# =============================================================================

#' @param country_data  Named list of outputs from prepare_country()
#' @param ted           Data frame: date, ted
#' @param vix           Data frame: date, vix
#' @param N             Number of hidden states per country
#' @param cyc           Max EM iterations
#' @param tol           Convergence tolerance

fit_all_countries <- function(country_data, ted, vix,
                              N = 2, cyc = 100, tol = 1e-4, lagged = TRUE) {
  
  # Build global Z once
  Z_full <- prepare_covariates(ted, vix, lagged = lagged)
  
  results <- vector("list", length(country_data))
  names(results) <- names(country_data)
  
  for (cc in names(country_data)) {
    cat(sprintf("\n--- Fitting %s ---\n", cc))
    
    cd <- country_data[[cc]]
    cd$date <- as.Date(cd$date)
    
    # Slice Z to country's dates
    merged <- merge(cd, Z_full, by = "date", sort = TRUE)
    if (any(c(is.na(merged$ted),is.na(merged$vix)))) {
      message(paste(cc,"Contains NAs in its covariates"))
      next
    }
    
    X <- as.matrix(merged[, c("d_y2y", "d_spread", "fx_return")])
    Z <- as.matrix(merged[, c("ted", "vix")])
    
    cat(sprintf("  T = %d observations after alignment\n", nrow(X)))
    
    model   <- fit_country_hmm(X, Z, N = N, cyc = cyc, tol = tol)
    u_state <- identify_unwind_state(model)
    
    results[[cc]] <- list(
      model        = model,
      unwind_state = u_state,
      dates        = merged$date
    )
    
    cat(sprintf("  Unwind state: %d\n", u_state))
    cat(sprintf("  Unwind state means -> d_y2y: %.4f  d_spread: %.4f  fx_ret: %.5f\n",
                model$Mu[u_state, 1], model$Mu[u_state, 2], model$Mu[u_state, 3]))
  }
  
  results
}




# =============================================================================
# SECTION 7: SIGNAL CONSTRUCTION ## OLD, WAS TRYING BACKTESTING
# =============================================================================

#' Build daily carry trade signal from regime probabilities
#'
#' @param hmm_results     Output from fit_all_countries()
#' @param fx_returns      Named list of named numeric vectors (log FX returns)
#' @param rate_diffs      Named list of named numeric vectors
#'                        (country 3m rate - USD 3m rate) / 252
#' @param n_long          Currencies to go long carry
#' @param n_short         Currencies to short carry
#' @param prob_threshold  Countries with unwind_prob in
#'                        (1-threshold, threshold) are excluded

build_signal <- function(hmm_results, fx_returns, rate_diffs,
                         n_long = 3, n_short = 3,
                         prob_threshold = 0.6) {
  
  countries <- names(hmm_results)
  all_dates <- sort(Reduce(intersect, lapply(hmm_results, function(r) r$dates)))
  T         <- length(all_dates)
  
  unwind_mat <- matrix(NA, T, length(countries),
                       dimnames = list(all_dates, countries))
  for (cc in countries) {
    r                <- hmm_results[[cc]]
    idx              <- match(all_dates, r$dates)
    unwind_mat[, cc] <- r$unwind_prob[idx]
  }
  
  fx_mat   <- matrix(NA, T, length(countries), dimnames = list(all_dates, countries))
  rate_mat <- matrix(NA, T, length(countries), dimnames = list(all_dates, countries))
  for (cc in countries) {
    fx_mat[, cc]   <- fx_returns[[cc]][match(all_dates, names(fx_returns[[cc]]))]
    rate_mat[, cc] <- rate_diffs[[cc]][match(all_dates, names(rate_diffs[[cc]]))]
  }
  
  pos_mat <- matrix(0, T, length(countries), dimnames = list(all_dates, countries))
  
  for (t in 2:T) {
    prob_yesterday <- unwind_mat[t - 1, ]
    if (any(is.na(prob_yesterday))) next
    
    ranked           <- order(prob_yesterday)
    long_candidates  <- ranked[prob_yesterday[ranked] < (1 - prob_threshold)]
    short_candidates <- ranked[prob_yesterday[ranked] >       prob_threshold]
    
    n_l <- min(n_long,  length(long_candidates))
    n_s <- min(n_short, length(short_candidates))
    
    if (n_l > 0) pos_mat[t, long_candidates[1:n_l]]  <-  1 / n_l
    if (n_s > 0) pos_mat[t, short_candidates[(length(short_candidates) - n_s + 1):
                                               length(short_candidates)]] <- -1 / n_s
  }
  
  # Long foreign ccy profits when it appreciates (FX rate falls) -> -fx_mat
  carry_pnl_mat <- pos_mat *  rate_mat
  fx_pnl_mat    <- pos_mat * (-fx_mat)
  total_pnl     <- rowSums(carry_pnl_mat + fx_pnl_mat, na.rm = TRUE)
  
  list(
    positions  = as.data.frame(pos_mat),
    carry_pnl  = as.data.frame(carry_pnl_mat),
    fx_pnl     = as.data.frame(fx_pnl_mat),
    total_pnl  = total_pnl,
    dates      = as.Date(all_dates)
  )
}


# =============================================================================
# SECTION 8: PERFORMANCE EVALUATION ## OLD, WAS FOR BACKTEST
# =============================================================================

evaluate_strategy <- function(daily_ret, ann = 252) {
  
  r        <- daily_ret[is.finite(daily_ret) & !is.na(daily_ret)]
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
# SECTION 9: DIAGNOSTICS
# =============================================================================

summarise_states <- function(hmm_results) {
  
  cat("\n=== STATE SUMMARY ===\n")
  for (cc in names(hmm_results)) {
    r  <- hmm_results[[cc]]
    Mu <- r$model$Mu
    x  <- r$model$x
    N  <- nrow(Mu)
    u  <- r$unwind_state
    
    cat(sprintf("\n%s (unwind state = %d)\n", cc, u))
    cat(sprintf("  %-8s  %9s  %10s  %9s\n", "State", "d_y2y", "d_spread", "fx_ret"))
    for (s in seq_len(N)) {
      marker <- if (s == u) " <-- UNWIND" else ""
      cat(sprintf("  State %d   %9.5f  %10.5f  %9.5f%s\n",
                  s, Mu[s, 1], Mu[s, 2], Mu[s, 3], marker))
    }
    
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

plot_results <- function(signal, hmm_results) {
  
  countries <- names(hmm_results)
  dates     <- signal$dates
  
  old_par <- par(mfrow = c(ceiling((length(countries) + 1) / 2), 2),
                 mar   = c(2, 3, 2, 1))
  on.exit(par(old_par))
  
  for (cc in countries) {
    r   <- hmm_results[[cc]]
    idx <- match(as.character(dates), r$dates)
    up  <- r$unwind_prob[idx]
    plot(dates, up, type = "l", lwd = 1.2, col = "#D6604D",
         ylim = c(0, 1), main = cc, xlab = "", ylab = "P(unwind)", cex.main = 1)
    abline(h = 0.6, lty = 2, col = "grey50")
    grid()
  }
  
  dev.new()
  cum_ret <- cumprod(1 + ifelse(is.na(signal$total_pnl), 0, signal$total_pnl))
  plot(dates, cum_ret, type = "l", lwd = 2, col = "#2166AC",
       main = "Carry Strategy — Cumulative Return",
       xlab = "Date", ylab = "Growth of $1")
  abline(h = 1, lty = 2, col = "grey50")
  grid()
}


# =============================================================================
# SECTION 11: USAGE
# =============================================================================
#
# Yields CSV needs columns: date, X2Y, X10Y  (y3m no longer required)
# FX CSV needs columns:     date, rate
# TED CSV:                  date, ted
# VIX CSV:                  date, vix
#
#   countries <- c("AUD", "EUR", "GBP", "JPY", "CAD", "NZD", "SEK", "NOK")
#
#   country_data <- lapply(setNames(countries, countries), function(cc) {
#     prepare_country(
#       read.csv(paste0("data/yields_", cc, ".csv")),
#       read.csv(paste0("data/fx_",     cc, ".csv"))
#     )
#   })
#
#   ted <- read.csv("data/ted_spread.csv")
#   vix <- read.csv("data/vix.csv")
#
#   hmm_results <- fit_all_countries(country_data, ted, vix, N = 2, cyc = 100)
#   summarise_states(hmm_results)
#
#   # Override unwind state if misidentified:
#   # hmm_results[["JPY"]]$unwind_state <- 1
#   # hmm_results[["JPY"]]$unwind_prob  <- hmm_results[["JPY"]]$decoded$prb[1, ]
#
#   signal <- build_signal(hmm_results, fx_returns, rate_diffs,
#                          n_long = 3, n_short = 3, prob_threshold = 0.6)
#   evaluate_strategy(signal$total_pnl)
#   plot_results(signal, hmm_results)
