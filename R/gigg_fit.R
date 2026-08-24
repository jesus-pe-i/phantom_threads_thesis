# Multi-chain fitting and network estimation for the GIGG block BVAR.
#
# Runs and pools Gibbs chains, restores coefficients to original units,
# and converts posterior coefficient means into block-network strengths.


# Chain setup -----

gigg_chain_seeds <- function(
    seed,
    chains) {
  
  as.integer(
    seed + seq_len(chains) - 1L
  )
}


make_gigg_initial_states <- function(
    structure,
    chains,
    overdispersed = TRUE) {
  
  if (!overdispersed) {
    
    return(
      replicate(
        chains,
        make_gigg_state(
          structure
        ),
        simplify = FALSE
      )
    )
  }
  
  profiles <- data.frame(
    sigma = c(
      1,
      0.5,
      2,
      1,
      1.5,
      0.75
    ),
    
    tau = c(
      1,
      0.25,
      4,
      2,
      0.5,
      1.5
    ),
    
    psi_tau = c(
      1,
      4,
      0.25,
      0.5,
      2,
      1.5
    ),
    
    gamma_odd = c(
      1,
      0.5,
      2,
      0.25,
      4,
      0.75
    ),
    
    gamma_even = c(
      1,
      2,
      0.5,
      4,
      0.25,
      1.5
    ),
    
    lambda_odd = c(
      1,
      0.5,
      2,
      4,
      0.25,
      1.5
    ),
    
    lambda_even = c(
      1,
      2,
      0.5,
      0.25,
      4,
      0.75
    )
  )
  
  base_tau2 <-
    1 /
    structure$data$T_p
  
  n_blocks <-
    structure$maps$n_blocks
  
  k <-
    structure$data$k
  
  N <-
    structure$data$N
  
  lapply(
    seq_len(chains),
    function(chain) {
      
      profile <- (
        chain - 1L
      ) %% nrow(profiles) + 1L
      
      gamma_odd <- seq_len(
        n_blocks
      ) %% 2L == 1L
      
      gamma2 <- rep(
        profiles$gamma_even[profile],
        n_blocks
      )
      
      gamma2[gamma_odd] <-
        profiles$gamma_odd[profile]
      
      lambda_odd <- seq_len(
        k * N
      ) %% 2L == 1L
      
      lambda2 <- rep(
        profiles$lambda_even[profile],
        k * N
      )
      
      lambda2[lambda_odd] <-
        profiles$lambda_odd[profile]
      
      lambda2 <- matrix(
        lambda2,
        nrow = k,
        ncol = N
      )
      
      make_gigg_state(
        structure = structure,
        
        sigma2 = rep(
          profiles$sigma[profile],
          N
        ),
        
        tau2 =
          base_tau2 *
          profiles$tau[profile],
        
        psi_tau =
          profiles$psi_tau[profile],
        
        gamma2 =
          gamma2,
        
        lambda2 =
          lambda2
      )
    }
  )
}

# Chain pooling -----

pool_gigg_chains <- function(
    chain_results) {
  
  retained <- vapply(
    chain_results,
    `[[`,
    integer(1L),
    "retained"
  )
  
  weights <- retained /
    sum(retained)
  
  weighted_mean <- function(
    name) {
    
    values <- lapply(
      chain_results,
      `[[`,
      name
    )
    
    Reduce(
      `+`,
      Map(
        `*`,
        values,
        weights
      )
    )
  }
  
  beta_mean <- weighted_mean(
    "beta_mean"
  )
  
  beta_mean_square <- weighted_mean(
    "beta_mean_square"
  )
  
  beta_variance <- pmax(
    beta_mean_square -
      beta_mean^2,
    0
  )
  
  list(
    beta_mean =
      beta_mean,
    
    beta_mean_square =
      beta_mean_square,
    
    beta_sd =
      sqrt(
        beta_variance
      ),
    
    sigma2_mean =
      weighted_mean(
        "sigma2_mean"
      ),
    
    tau2_mean =
      weighted_mean(
        "tau2_mean"
      ),
    
    psi_tau_mean =
      weighted_mean(
        "psi_tau_mean"
      ),
    
    gamma2_mean =
      weighted_mean(
        "gamma2_mean"
      ),
    
    lambda2_mean =
      weighted_mean(
        "lambda2_mean"
      ),
    
    tau2_gamma2_mean =
      weighted_mean(
        "tau2_gamma2_mean"
      ),
    
    effective_variance_mean =
      weighted_mean(
        "effective_variance_mean"
      ),
    
    retained =
      retained,
    
    retained_total =
      sum(retained)
  )
}


# Model fit -----

fit_gigg <- function(
    Y_list,
    p_lags,
    chains = 4L,
    burnin = 500,
    draws = 1000L,
    thin = 1L,
    seed = 991L,
    beta_algorithm = c(
      "auto",
      "chol",
      "bhattacharya"
    ),
    use_asis = TRUE,
    asis_every = 1L,
    selected_eps = 1e-12,
    monitor_beta = NULL,
    monitor_gamma = NULL,
    monitor_lambda = NULL,
    initial_states = NULL,
    overdispersed_initial_states = TRUE,
    keep_chain_results = FALSE) {
  
  beta_algorithm <- match.arg(
    beta_algorithm
  )
  
  
  # Model objects -----
  
  structure <- make_gigg_structure(
    Y_list = Y_list,
    p_lags = p_lags
  )
  
  prior <- make_gigg_prior(
    structure
  )
  
  control <- make_gigg_control(
    burnin = burnin,
    draws = draws,
    thin = thin,
    beta_algorithm =
      beta_algorithm,
    use_asis =
      use_asis,
    asis_every =
      asis_every,
    keep_all = FALSE,
    monitor_beta =
      monitor_beta,
    monitor_gamma =
      monitor_gamma,
    monitor_lambda =
      monitor_lambda
  )
  
  
  # Initial states -----
  
  if (is.null(initial_states)) {
    
    initial_states <-
      make_gigg_initial_states(
        structure = structure,
        chains = chains,
        overdispersed =
          overdispersed_initial_states
      )
  }
  
  if (length(initial_states) != chains) {
    stop(
      "initial_states must contain one state per chain"
    )
  }
  
  
  # Chains -----
  
  chain_seeds <- gigg_chain_seeds(
    seed = seed,
    chains = chains
  )
  
  chain_results <- lapply(
    seq_len(chains),
    function(chain) {
      
      run_gigg_chain(
        structure = structure,
        prior = prior,
        control = control,
        state =
          initial_states[[chain]],
        seed =
          chain_seeds[chain]
      )
    }
  )
  
  algorithms <- vapply(
    chain_results,
    `[[`,
    character(1L),
    "beta_algorithm"
  )
  
  if (
    length(
      unique(algorithms)
    ) != 1L
  ) {
    stop(
      "GIGG chains resolved different beta algorithms"
    )
  }
  
  
  # Posterior pooling -----
  
  pooled <- pool_gigg_chains(
    chain_results
  )
  
  beta_standardized <-
    pooled$beta_mean
  
  sigma2_standardized <-
    pooled$sigma2_mean
  
  beta_hat <- backtransform_var_beta(
    beta = beta_standardized,
    scale =
      structure$preprocessing$scale,
    p_lags =
      structure$dimensions$p_lags
  )
  
  sigma2_hat <-
    sigma2_standardized *
    structure$preprocessing$scale^2
  
  
  # Network estimate -----
  
  A_hat_lag <- beta_to_A_list(
    beta = beta_hat,
    p_lags =
      structure$dimensions$p_lags,
    series_names =
      structure$preprocessing$
      series_names
  )
  
  s_hat_lag <- A_list_to_s_lag(
    A_list = A_hat_lag,
    n_units =
      structure$dimensions$n_units,
    m =
      structure$dimensions$m
  )
  
  s_hat_unit_max <- s_lag_to_unit(
    s_lag = s_hat_lag,
    method = "max"
  )
  
  s_hat_unit_rms <- s_lag_to_unit(
    s_lag = s_hat_lag,
    method = "rms"
  )
  
  
  # Output -----
  
  fit <- list(
    model =
      "gigg",
    
    beta_hat =
      beta_hat,
    
    beta_standardized =
      beta_standardized,
    
    beta_sd_standardized =
      pooled$beta_sd,
    
    sigma2_hat =
      sigma2_hat,
    
    sigma2_standardized =
      sigma2_standardized,
    
    tau2_mean =
      pooled$tau2_mean,
    
    psi_tau_mean =
      pooled$psi_tau_mean,
    
    gamma2_mean =
      pooled$gamma2_mean,
    
    lambda2_mean =
      pooled$lambda2_mean,
    
    tau2_gamma2_mean =
      pooled$tau2_gamma2_mean,
    
    effective_variance_mean =
      pooled$effective_variance_mean,
    
    A_hat_lag =
      A_hat_lag,
    
    s_hat_lag =
      s_hat_lag,
    
    s_hat_unit_max =
      s_hat_unit_max,
    
    s_hat_unit_rms =
      s_hat_unit_rms,
    
    score_lag =
      s_hat_lag,
    
    selected_lag =
      s_hat_lag >
      selected_eps,
    
    retained_per_chain =
      pooled$retained,
    
    retained_total =
      pooled$retained_total,
    
    chain_seeds =
      chain_seeds,
    
    beta_algorithm_requested =
      beta_algorithm,
    
    beta_algorithm_resolved =
      algorithms[1L],
    
    dimensions =
      structure$dimensions,
    
    preprocessing =
      structure$preprocessing,
    
    prior =
      prior,
    
    control =
      control
  )
  
  if (keep_chain_results) {
    
    fit$chain_results <-
      chain_results
  }
  
  fit
}