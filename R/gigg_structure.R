# Structure, prior, and initial state for the GIGG block BVAR.
#
# Builds edge-lag block maps on the common standardized VAR design,
# defines the frozen GIGG prior, and constructs Gibbs initial states.


# Block map -----

make_gigg_block_id <- function(
    n_units,
    m,
    p_lags) {
  
  N <- n_units * m
  k <- N * p_lags
  
  predictor <- rep(
    seq_len(N),
    times = p_lags
  )
  
  lag <- rep(
    seq_len(p_lags),
    each = N
  )
  
  predictor_unit <- (
    predictor - 1L
  ) %/% m + 1L
  
  response_unit <- (
    seq_len(N) - 1L
  ) %/% m + 1L
  
  block_id <- matrix(
    rep(
      (
        lag - 1L
      ) *
        n_units^2 +
        predictor_unit,
      times = N
    ),
    nrow = k,
    ncol = N
  ) +
    matrix(
      rep(
        (
          response_unit - 1L
        ) *
          n_units,
        each = k
      ),
      nrow = k,
      ncol = N
    )
  
  storage.mode(
    block_id
  ) <- "integer"
  
  block_id
}


# Structure -----

make_gigg_structure <- function(
    Y_list,
    p_lags) {
  
  prepared <- prepare_var_data(
    Y_list = Y_list,
    p_lags = p_lags,
    standardize = TRUE
  )
  
  n_units <- prepared$n_units
  m <- prepared$m
  N <- prepared$n_series
  k <- prepared$n_coef
  p_lags <- prepared$p_lags
  
  block_id <- make_gigg_block_id(
    n_units = n_units,
    m = m,
    p_lags = p_lags
  )
  
  n_blocks <- as.integer(
    p_lags *
      n_units^2
  )
  
  block_counts <- tabulate(
    as.vector(
      block_id
    ),
    nbins = n_blocks
  )
  
  if (
    any(
      block_counts != m^2
    )
  ) {
    stop(
      "Every GIGG block must contain m^2 coefficients"
    )
  }
  
  list(
    data = list(
      Y =
        prepared$Y,
      
      X =
        prepared$X,
      
      T_p =
        prepared$T_eff,
      
      K =
        n_units,
      
      m =
        m,
      
      N =
        N,
      
      p_lags =
        p_lags,
      
      k =
        k
    ),
    
    maps = list(
      block_id =
        block_id,
      
      n_blocks =
        n_blocks
    ),
    
    dimensions = list(
      n_units =
        n_units,
      
      m =
        m,
      
      N =
        N,
      
      p_lags =
        p_lags,
      
      k =
        k,
      
      T_full =
        prepared$T_full,
      
      T_p =
        prepared$T_eff
    ),
    
    preprocessing = list(
      standardize =
        TRUE,
      
      centering =
        FALSE,
      
      scale =
        prepared$scale,
      
      series_names =
        prepared$series_names
    )
  )
}


# Prior -----

make_gigg_prior <- function(
    structure) {
  
  list(
    sigma_a =
      3,
    
    sigma_b =
      2,
    
    tau_df =
      5,
    
    tau_scale =
      3 /
      sqrt(
        5 *
          structure$data$T_p
      ),
    
    gamma_shape =
      0.5,
    
    gamma_rate =
      1,
    
    lambda_shape =
      2.5,
    
    lambda_scale =
      1
  )
}


# Initial state -----

make_gigg_state <- function(
    structure,
    beta = NULL,
    sigma2 = NULL,
    tau2 = NULL,
    psi_tau = NULL,
    gamma2 = NULL,
    lambda2 = NULL) {
  
  k <- structure$data$k
  N <- structure$data$N
  T_p <- structure$data$T_p
  n_blocks <- structure$maps$n_blocks
  
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
  
  if (is.null(tau2)) {
    
    tau2 <-
      1 / T_p
  }
  
  if (is.null(psi_tau)) {
    
    psi_tau <-
      1
  }
  
  if (is.null(gamma2)) {
    
    gamma2 <- rep(
      1,
      n_blocks
    )
  }
  
  if (is.null(lambda2)) {
    
    lambda2 <- matrix(
      1,
      nrow = k,
      ncol = N
    )
  }
  
  list(
    beta =
      unname(
        as.matrix(
          beta
        )
      ),
    
    sigma2 =
      as.numeric(
        sigma2
      ),
    
    tau2 =
      as.numeric(
        tau2
      ),
    
    psi_tau =
      as.numeric(
        psi_tau
      ),
    
    gamma2 =
      as.numeric(
        gamma2
      ),
    
    lambda2 =
      unname(
        as.matrix(
          lambda2
        )
      )
  )
}