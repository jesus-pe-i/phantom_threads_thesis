# R-to-C++ execution layer for the M3 block BVAR.
#
# Builds MCMC controls and C++ inputs, loads the native sampler,
# and runs one Gibbs chain with posterior and monitored draws.


# C++ loading -----

load_m3_cpp <- function(
    rebuild = FALSE) {
  
  Rcpp::sourceCpp(
    "src/m3_bvar_core.cpp",
    rebuild = rebuild
  )
}


# MCMC control -----

make_m3_control <- function(
    burnin,
    draws,
    thin = 1L,
    beta_algorithm = c(
      "auto",
      "chol",
      "bhattacharya"
    ),
    q_update = c(
      "gibbs_end",
      "gibbs_after_beta",
      "transport_after_beta"
    ),
    q_transport_global_probability = 0.25,
    use_c_asis = TRUE,
    c_asis_every = 1L,
    c_asis_slice_width = 1,
    c_asis_maximum_step_out = 20L,
    keep_all = FALSE,
    verbose = FALSE,
    monitor_beta = NULL,
    monitor_lambda = NULL) {
  
  beta_algorithm <- match.arg(
    beta_algorithm
  )
  
  q_update <- match.arg(
    q_update
  )
  
  if (is.null(monitor_beta)) {
    
    monitor_beta <- matrix(
      integer(0L),
      nrow = 0L,
      ncol = 2L
    )
  }
  
  monitor_beta <- as.matrix(
    monitor_beta
  )
  
  if (ncol(monitor_beta) != 2L) {
    stop(
      "monitor_beta must have two columns"
    )
  }
  
  storage.mode(monitor_beta) <-
    "integer"
  
  colnames(monitor_beta) <- c(
    "row",
    "equation"
  )
  
  if (is.null(monitor_lambda)) {
    monitor_lambda <- integer(0L)
  }
  
  monitor_lambda <- as.integer(
    monitor_lambda
  )
  
  list(
    burnin = as.integer(burnin),
    draws = as.integer(draws),
    thin = as.integer(thin),
    
    beta_algorithm =
      beta_algorithm,
    
    q_update =
      q_update,
    
    q_transport_global_probability =
      as.numeric(
        q_transport_global_probability
      ),
    
    use_c_asis =
      use_c_asis,
    
    c_asis_every =
      as.integer(c_asis_every),
    
    c_asis_slice_width =
      as.numeric(
        c_asis_slice_width
      ),
    
    c_asis_maximum_step_out =
      as.integer(
        c_asis_maximum_step_out
      ),
    
    scale_transport =
      "none",
    
    scale_transport_frequency =
      1L,
    
    scale_transport_global_log_sd =
      0.25,
    
    scale_transport_group_log_sd =
      0.25,
    
    keep_all = keep_all,
    verbose = verbose,
    
    monitor_beta =
      monitor_beta,
    
    monitor_lambda =
      monitor_lambda
  )
}

# C++ input -----

make_m3_input <- function(
    structure,
    prior,
    state,
    control) {
  
  list(
    data = structure$data,
    maps = structure$maps,
    prior = prior,
    state = state,
    control = control
  )
}


# Draw formatting -----

m3_draw_matrix <- function(
    x,
    draws,
    parameters) {
  
  x <- as.matrix(x)
  
  if (
    identical(
      dim(x),
      c(parameters, draws)
    )
  ) {
    x <- t(x)
  }
  
  if (
    !identical(
      dim(x),
      c(draws, parameters)
    )
  ) {
    stop(
      "Unexpected M3 draw dimensions"
    )
  }
  
  x
}


m3_monitor_names <- function(
    structure,
    beta_index,
    lambda_index) {
  
  N <- structure$data$N
  K <- structure$dimensions$n_units
  
  beta_names <- if (
    nrow(beta_index) == 0L
  ) {
    
    character(0L)
    
  } else {
    
    lag <- (
      beta_index[, 1L] - 1L
    ) %/% N + 1L
    
    sender <- (
      beta_index[, 1L] - 1L
    ) %% N + 1L
    
    paste0(
      "beta_l",
      lag,
      "_s",
      sender,
      "_r",
      beta_index[, 2L]
    )
  }
  
  if (length(lambda_index) == 0L) {
    
    return(
      list(
        beta = beta_names,
        lambda2 = character(0L),
        xi = character(0L)
      )
    )
  }
  
  zero <- lambda_index - 1L
  
  lag <- zero %/% K^2L + 1L
  
  within_lag <- zero %% K^2L
  
  receiver <- within_lag %/%
    K + 1L
  
  sender <- within_lag %%
    K + 1L
  
  suffix <- paste0(
    "_l",
    lag,
    "_s",
    sender,
    "_r",
    receiver
  )
  
  list(
    beta = beta_names,
    lambda2 = paste0(
      "lambda2",
      suffix
    ),
    xi = paste0(
      "xi",
      suffix
    )
  )
}


# One chain -----

run_m3_chain <- function(
    structure,
    prior,
    control,
    state = NULL,
    seed = 991L) {
  
  if (is.null(state)) {
    
    state <- make_m3_state(
      structure = structure,
      prior = prior
    )
  }
  
  input <- make_m3_input(
    structure = structure,
    prior = prior,
    state = state,
    control = control
  )
  
  seed <- as.integer(seed)
  
  set.seed(
    seed
  )
  
  elapsed <- system.time(
    raw <- m3_bvar_chain_cpp(
      input
    )
  )
  
  retained <- as.integer(
    raw$info$retained
  )
  
  
  # Monitors -----
  
  beta_index <- as.matrix(
    raw$monitor$beta_index
  )
  
  lambda_index <- as.integer(
    raw$monitor$lambda_index
  )
  
  draws <- list(
    beta = m3_draw_matrix(
      raw$monitor$beta,
      retained,
      nrow(beta_index)
    ),
    
    lambda2 = m3_draw_matrix(
      raw$monitor$lambda2,
      retained,
      length(lambda_index)
    ),
    
    xi = m3_draw_matrix(
      raw$monitor$xi,
      retained,
      length(lambda_index)
    ),
    
    sigma2 = m3_draw_matrix(
      raw$draws$sigma2,
      retained,
      structure$data$N
    ),
    
    tau2 = m3_draw_matrix(
      raw$draws$tau2,
      retained,
      1L
    ),
    
    psi_tau = m3_draw_matrix(
      raw$draws$psi_tau,
      retained,
      1L
    ),
    
    c2 = m3_draw_matrix(
      raw$draws$c2,
      retained,
      structure$maps$n_c
    ),
    
    psi_c = m3_draw_matrix(
      raw$draws$psi_c,
      retained,
      structure$maps$n_c
    ),
    
    q = m3_draw_matrix(
      raw$draws$q,
      retained,
      structure$maps$n_q
    ),
    
    q_index = m3_draw_matrix(
      raw$draws$q_index,
      retained,
      structure$maps$n_q
    )
  )
  
  monitor_names <- m3_monitor_names(
    structure = structure,
    beta_index = beta_index,
    lambda_index = lambda_index
  )
  
  colnames(draws$beta) <-
    monitor_names$beta
  
  colnames(draws$lambda2) <-
    monitor_names$lambda2
  
  colnames(draws$xi) <-
    monitor_names$xi
  
  
  # ASIS diagnostics -----
  
  c_asis <- raw$c_asis
  
  c_asis$enabled <- isTRUE(
    c_asis$enabled
  )
  
  c_asis$every <- as.integer(
    c_asis$every
  )
  
  c_asis$slice_width <- as.numeric(
    c_asis$slice_width
  )
  
  c_asis$maximum_step_out <-
    as.integer(
      c_asis$maximum_step_out
    )
  
  c_asis$updates <- as.integer(
    c_asis$updates
  )
  
  c_asis$moves <- as.integer(
    c_asis$moves
  )
  
  
  # Output -----
  
  resolved_algorithm <-
    raw$info$beta_algorithm
  
  if (is.null(resolved_algorithm)) {
    
    resolved_algorithm <-
      raw$info$
      beta_algorithm_resolved
  }
  
  list(
    seed = seed,
    retained = retained,
    
    beta_mean = as.matrix(
      raw$posterior$beta_mean
    ),
    
    beta_mean_square = as.matrix(
      raw$posterior$
        beta_mean_square
    ),
    
    sigma2_mean = colMeans(
      draws$sigma2
    ),
    
    tau2_mean = mean(
      draws$tau2
    ),
    
    c2_mean = colMeans(
      draws$c2
    ),
    
    q_mean = colMeans(
      draws$q
    ),
    
    draws = draws,
    
    full_draws = raw$full_draws,
    
    monitor = list(
      beta = beta_index,
      lambda2 = lambda_index
    ),
    
    final_state =
      raw$final_state,
    
    q_transport =
      raw$q_transport,
    
    c_asis = c_asis,
    timing = raw$timing,
    
    beta_algorithm =
      resolved_algorithm,
    
    q_update =
      raw$info$q_update,
    
    elapsed_seconds =
      unname(
        elapsed["elapsed"]
      )
  )
}