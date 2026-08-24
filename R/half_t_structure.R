# Structure, prior, and initial state for the coefficient-wise half-t BVAR.
#
# Uses scale-only standardization, Minnesota lag weights, and either one
# global shrinkage scale or separate self-diagonal and remaining scales.


# Structure -----

make_half_t_structure <- function(
    Y_list,
    p_lags,
    global_grouping = c(
      "all",
      "self_diagonal"
    ),
    standardize = TRUE) {
  
  global_grouping <- match.arg(
    global_grouping
  )
  
  prepared <- prepare_var_data(
    Y_list = Y_list,
    p_lags = p_lags,
    standardize = standardize
  )
  
  N <- prepared$n_series
  k <- prepared$n_coef
  
  lag_id <- rep(
    seq_len(prepared$p_lags),
    each = N
  )
  
  predictor_series <- rep(
    seq_len(N),
    times = prepared$p_lags
  )
  
  phi2 <- 1 / lag_id^2
  
  
  # Global shrinkage grouping -----
  
  if (global_grouping == "all") {
    
    grouping_code <- 0L
    n_tau <- 1L
    
    tau_group_counts <- as.integer(
      k * N
    )
    
  } else {
    
    grouping_code <- 1L
    n_tau <- 2L
    
    self_count <- N * prepared$p_lags
    
    tau_group_counts <- as.integer(
      c(
        self_count,
        k * N - self_count
      )
    )
  }
  
  
  # Structure -----
  
  list(
    data = list(
      Y = prepared$Y,
      X = prepared$X,
      XtX = crossprod(prepared$X),
      XtY = crossprod(
        prepared$X,
        prepared$Y
      ),
      YtY = colSums(
        prepared$Y^2
      ),
      T_p = prepared$T_eff,
      N = N,
      k = k
    ),
    
    maps = list(
      lag_id = as.integer(lag_id),
      predictor_series = as.integer(
        predictor_series
      ),
      phi2 = as.numeric(phi2),
      grouping_code = grouping_code,
      n_tau = n_tau,
      tau_group_counts = tau_group_counts
    ),
    
    dimensions = list(
      T_full = prepared$T_full,
      T_p = prepared$T_eff,
      n_units = prepared$n_units,
      m = prepared$m,
      N = N,
      p_lags = prepared$p_lags,
      k = k
    ),
    
    preprocessing = list(
      standardize = standardize,
      scale = prepared$scale,
      series_names = prepared$series_names
    ),
    
    global_grouping = global_grouping
  )
}


# Prior -----

make_half_t_prior <- function(
    lambda_df = 3,
    lambda_scale = 1,
    tau_df = 10,
    tau_scale = 1,
    sigma_a = 3,
    sigma_b = 2) {
  
  list(
    lambda_df = lambda_df,
    lambda_scale = lambda_scale,
    tau_df = tau_df,
    tau_scale = tau_scale,
    sigma_a = sigma_a,
    sigma_b = sigma_b
  )
}


# Initial state -----

make_half_t_state <- function(
    structure,
    beta = NULL,
    sigma2 = NULL,
    lambda2 = NULL,
    nu = NULL,
    tau2 = NULL,
    xi = NULL) {
  
  k <- structure$data$k
  N <- structure$data$N
  n_tau <- structure$maps$n_tau
  
  if (is.null(beta)) {
    beta <- matrix(
      0,
      nrow = k,
      ncol = N
    )
  }
  
  if (is.null(sigma2)) {
    sigma2 <- rep(
      1,
      N
    )
  }
  
  if (is.null(lambda2)) {
    lambda2 <- matrix(
      1,
      nrow = k,
      ncol = N
    )
  }
  
  if (is.null(nu)) {
    nu <- matrix(
      1,
      nrow = k,
      ncol = N
    )
  }
  
  if (is.null(tau2)) {
    tau2 <- rep(
      1,
      n_tau
    )
  }
  
  if (is.null(xi)) {
    xi <- rep(
      1,
      n_tau
    )
  }
  
  list(
    beta = beta,
    sigma2 = as.numeric(sigma2),
    lambda2 = lambda2,
    nu = nu,
    tau2 = as.numeric(tau2),
    xi = as.numeric(xi)
  )
}