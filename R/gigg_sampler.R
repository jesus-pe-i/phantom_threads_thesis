# R-to-C++ execution layer for the GIGG block BVAR.
#
# Builds MCMC controls and C++ inputs, loads the native sampler,
# and runs one Gibbs chain with posterior and monitored draws.


# C++ loading -----

load_gigg_cpp <- function(
    rebuild = FALSE) {
  
  Rcpp::sourceCpp(
    "src/gigg_bvar_core.cpp",
    rebuild = rebuild
  )
}


# Monitoring -----

gigg_pair_monitor <- function(
    monitor) {
  
  if (is.null(monitor)) {
    
    return(
      matrix(
        integer(0L),
        nrow = 0L,
        ncol = 2L
      )
    )
  }
  
  monitor <- as.matrix(
    monitor
  )
  
  if (
    is.null(dim(monitor)) ||
    ncol(monitor) != 2L
  ) {
    stop(
      "Coefficient monitors must have two columns"
    )
  }
  
  storage.mode(monitor) <-
    "integer"
  
  colnames(monitor) <- c(
    "row",
    "equation"
  )
  
  monitor
}


gigg_block_monitor <- function(
    monitor) {
  
  if (is.null(monitor)) {
    return(
      integer(0L)
    )
  }
  
  as.integer(
    monitor
  )
}


# MCMC control -----

make_gigg_control <- function(
    burnin,
    draws,
    thin = 1L,
    beta_algorithm = c(
      "auto",
      "chol",
      "bhattacharya"
    ),
    use_asis = TRUE,
    asis_every = 1L,
    keep_all = FALSE,
    monitor_beta = NULL,
    monitor_gamma = NULL,
    monitor_lambda = NULL) {
  
  beta_algorithm <- match.arg(
    beta_algorithm
  )
  
  list(
    burnin =
      as.integer(
        burnin
      ),
    
    draws =
      as.integer(
        draws
      ),
    
    thin =
      as.integer(
        thin
      ),
    
    beta_algorithm =
      beta_algorithm,
    
    use_asis =
      use_asis,
    
    asis_every =
      as.integer(
        asis_every
      ),
    
    keep_all =
      keep_all,
    
    monitor_beta =
      gigg_pair_monitor(
        monitor_beta
      ),
    
    monitor_gamma =
      gigg_block_monitor(
        monitor_gamma
      ),
    
    monitor_lambda =
      gigg_pair_monitor(
        monitor_lambda
      )
  )
}


# C++ input -----

make_gigg_input <- function(
    structure,
    prior,
    state,
    control) {
  
  list(
    data =
      structure$data,
    
    maps =
      structure$maps,
    
    prior =
      prior,
    
    state = list(
      beta =
        unname(
          as.matrix(
            state$beta
          )
        ),
      
      sigma2 =
        as.numeric(
          state$sigma2
        ),
      
      tau2 =
        as.numeric(
          state$tau2
        ),
      
      psi_tau =
        as.numeric(
          state$psi_tau
        ),
      
      gamma2 =
        as.numeric(
          state$gamma2
        ),
      
      lambda2 =
        unname(
          as.matrix(
            state$lambda2
          )
        )
    ),
    
    control =
      control
  )
}


# Draw formatting -----

gigg_draw_matrix <- function(
    x,
    draws,
    parameters) {
  
  draws <- as.integer(
    draws
  )
  
  parameters <- as.integer(
    parameters
  )
  
  if (parameters == 0L) {
    
    return(
      matrix(
        numeric(0L),
        nrow = draws,
        ncol = 0L
      )
    )
  }
  
  x <- as.matrix(
    x
  )
  
  if (
    identical(
      dim(x),
      c(
        parameters,
        draws
      )
    )
  ) {
    x <- t(
      x
    )
  }
  
  if (
    !identical(
      dim(x),
      c(
        draws,
        parameters
      )
    )
  ) {
    stop(
      "Unexpected GIGG draw dimensions"
    )
  }
  
  x
}


gigg_monitor_names <- function(
    structure,
    beta_index,
    gamma_index,
    lambda_index) {
  
  N <- structure$data$N
  K <- structure$dimensions$n_units
  
  pair_suffix <- function(
    index) {
    
    if (nrow(index) == 0L) {
      return(
        character(0L)
      )
    }
    
    lag <- (
      index[, 1L] - 1L
    ) %/% N + 1L
    
    sender <- (
      index[, 1L] - 1L
    ) %% N + 1L
    
    receiver <-
      index[, 2L]
    
    paste0(
      "_l",
      lag,
      "_s",
      sender,
      "_r",
      receiver
    )
  }
  
  
  add_prefix <- function(
    prefix,
    suffix) {
    
    if (length(suffix) == 0L) {
      return(
        character(0L)
      )
    }
    
    paste0(
      prefix,
      suffix
    )
  }
  
  
  beta_suffix <- pair_suffix(
    beta_index
  )
  
  lambda_suffix <- pair_suffix(
    lambda_index
  )
  
  
  if (length(gamma_index) == 0L) {
    
    gamma_suffix <-
      character(0L)
    
  } else {
    
    zero <-
      gamma_index - 1L
    
    lag <- zero %/%
      K^2L + 1L
    
    within_lag <-
      zero %% K^2L
    
    receiver <-
      within_lag %/%
      K + 1L
    
    sender <-
      within_lag %%
      K + 1L
    
    gamma_suffix <- paste0(
      "_l",
      lag,
      "_su",
      sender,
      "_ru",
      receiver
    )
  }
  
  
  list(
    beta = add_prefix(
      "beta",
      beta_suffix
    ),
    
    gamma2 = add_prefix(
      "gamma2",
      gamma_suffix
    ),
    
    tau2_gamma2 = add_prefix(
      "tau2_gamma2",
      gamma_suffix
    ),
    
    lambda2 = add_prefix(
      "lambda2",
      lambda_suffix
    ),
    
    effective_variance = add_prefix(
      "effective_variance",
      lambda_suffix
    )
  )
}


# One chain -----

run_gigg_chain <- function(
    structure,
    prior,
    control,
    state = NULL,
    seed = 991L) {
  
  if (is.null(state)) {
    
    state <- make_gigg_state(
      structure
    )
  }
  
  input <- make_gigg_input(
    structure = structure,
    prior = prior,
    state = state,
    control = control
  )
  
  seed <- as.integer(
    seed
  )
  
  set.seed(
    seed
  )
  
  elapsed <- system.time(
    raw <- gigg_bvar_chain_cpp(
      input
    )
  )[["elapsed"]]
  
  retained <- as.integer(
    raw$info$retained
  )
  
  
  # Monitors -----
  
  beta_index <- as.matrix(
    raw$monitor$beta_index
  )
  
  gamma_index <- as.integer(
    raw$monitor$gamma_index
  )
  
  lambda_index <- as.matrix(
    raw$monitor$lambda_index
  )
  
  storage.mode(beta_index) <-
    "integer"
  
  storage.mode(lambda_index) <-
    "integer"
  
  draws <- list(
    beta =
      gigg_draw_matrix(
        raw$monitor$beta,
        retained,
        nrow(beta_index)
      ),
    
    gamma2 =
      gigg_draw_matrix(
        raw$monitor$gamma2,
        retained,
        length(gamma_index)
      ),
    
    tau2_gamma2 =
      gigg_draw_matrix(
        raw$monitor$tau2_gamma2,
        retained,
        length(gamma_index)
      ),
    
    lambda2 =
      gigg_draw_matrix(
        raw$monitor$lambda2,
        retained,
        nrow(lambda_index)
      ),
    
    effective_variance =
      gigg_draw_matrix(
        raw$monitor$effective_variance,
        retained,
        nrow(lambda_index)
      ),
    
    sigma2 =
      gigg_draw_matrix(
        raw$draws$sigma2,
        retained,
        structure$data$N
      ),
    
    tau2 =
      gigg_draw_matrix(
        raw$draws$tau2,
        retained,
        1L
      ),
    
    psi_tau =
      gigg_draw_matrix(
        raw$draws$psi_tau,
        retained,
        1L
      )
  )
  
  monitor_names <- gigg_monitor_names(
    structure = structure,
    beta_index = beta_index,
    gamma_index = gamma_index,
    lambda_index = lambda_index
  )
  
  colnames(draws$beta) <-
    monitor_names$beta
  
  colnames(draws$gamma2) <-
    monitor_names$gamma2
  
  colnames(draws$tau2_gamma2) <-
    monitor_names$tau2_gamma2
  
  colnames(draws$lambda2) <-
    monitor_names$lambda2
  
  colnames(draws$effective_variance) <-
    monitor_names$effective_variance
  
  colnames(draws$sigma2) <-
    paste0(
      "sigma2_",
      seq_len(
        ncol(
          draws$sigma2
        )
      )
    )
  
  colnames(draws$tau2) <-
    "tau2"
  
  colnames(draws$psi_tau) <-
    "psi_tau"
  
  
  # Output -----
  
  list(
    seed =
      seed,
    
    retained =
      retained,
    
    beta_mean =
      as.matrix(
        raw$posterior$beta_mean
      ),
    
    beta_mean_square =
      as.matrix(
        raw$posterior$beta_mean_square
      ),
    
    sigma2_mean =
      as.numeric(
        raw$posterior$sigma2_mean
      ),
    
    tau2_mean =
      as.numeric(
        raw$posterior$tau2_mean
      ),
    
    psi_tau_mean =
      as.numeric(
        raw$posterior$psi_tau_mean
      ),
    
    gamma2_mean =
      as.numeric(
        raw$posterior$gamma2_mean
      ),
    
    lambda2_mean =
      as.matrix(
        raw$posterior$lambda2_mean
      ),
    
    tau2_gamma2_mean =
      as.numeric(
        raw$posterior$tau2_gamma2_mean
      ),
    
    effective_variance_mean =
      as.matrix(
        raw$posterior$effective_variance_mean
      ),
    
    draws =
      draws,
    
    monitor = list(
      beta =
        beta_index,
      
      gamma2 =
        gamma_index,
      
      lambda2 =
        lambda_index
    ),
    
    full_draws =
      raw$full_draws,
    
    final_state =
      raw$final_state,
    
    gig =
      raw$gig,
    
    asis =
      raw$asis,
    
    timing =
      raw$timing,
    
    info =
      raw$info,
    
    beta_algorithm =
      raw$info$beta_algorithm_resolved,
    
    elapsed_seconds =
      unname(
        elapsed
      )
  )
}