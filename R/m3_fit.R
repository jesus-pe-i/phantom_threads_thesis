# Multi-chain fitting and network estimation for the M3 block BVAR.
#
# Runs and pools Gibbs chains, restores coefficients to original units,
# and converts posterior coefficient means into block-network strengths.


# Chain setup -----

m3_chain_seeds <- function(
    seed,
    chains) {
  
  as.integer(
    seed + seq_len(chains) - 1L
  )
}


make_m3_initial_states <- function(
    structure,
    prior,
    chains,
    overdispersed = TRUE) {
  
  if (!overdispersed) {
    
    return(
      replicate(
        chains,
        make_m3_state(
          structure,
          prior
        ),
        simplify = FALSE
      )
    )
  }
  
  profiles <- data.frame(
    tau = c(
      1,
      0.25,
      4,
      1,
      0.5,
      2
    ),
    foreign_c = c(
      1,
      1,
      1,
      0.25,
      4,
      0.5
    ),
    sigma = c(
      1,
      0.5,
      2,
      1,
      1.5,
      0.75
    ),
    lambda_odd = c(
      1,
      0.5,
      2,
      0.25,
      4,
      0.75
    ),
    lambda_even = c(
      1,
      2,
      0.5,
      4,
      0.25,
      1.5
    )
  )
  
  centre <- which.min(
    abs(
      prior$q_grid
    )
  )
  
  lower <- max(
    1L,
    centre - 4L
  )
  
  upper <- min(
    length(prior$q_grid),
    centre + 4L
  )
  
  q_profiles <- rbind(
    c(centre, centre),
    c(lower, upper),
    c(upper, lower),
    c(lower, lower),
    c(upper, upper),
    c(centre, upper)
  )
  
  lapply(
    seq_len(chains),
    function(chain) {
      
      profile <- (
        chain - 1L
      ) %% 6L + 1L
      
      state <- make_m3_state(
        structure,
        prior
      )
      
      state$sigma2[] <-
        profiles$sigma[profile]
      
      state$tau2 <-
        profiles$tau[profile]
      
      if (length(state$c2) > 1L) {
        
        state$c2[2L] <-
          profiles$foreign_c[
            profile
          ]
      }
      
      odd <- seq_along(
        state$lambda2
      ) %% 2L == 1L
      
      state$lambda2[odd] <-
        profiles$lambda_odd[
          profile
        ]
      
      state$lambda2[!odd] <-
        profiles$lambda_even[
          profile
        ]
      
      q_count <- min(
        2L,
        length(state$q_index)
      )
      
      state$q_index[
        seq_len(q_count)
      ] <- q_profiles[
        profile,
        seq_len(q_count)
      ]
      
      state
    }
  )
}


# Chain pooling -----

pool_m3_chains <- function(
    chain_results) {
  
  retained <- vapply(
    chain_results,
    `[[`,
    integer(1L),
    "retained"
  )
  
  weights <- retained /
    sum(retained)
  
  weighted_mean <- function(name) {
    
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
      sqrt(beta_variance),
    
    sigma2_mean =
      weighted_mean(
        "sigma2_mean"
      ),
    
    tau2_mean =
      weighted_mean(
        "tau2_mean"
      ),
    
    c2_mean =
      weighted_mean(
        "c2_mean"
      ),
    
    q_mean =
      weighted_mean(
        "q_mean"
      ),
    
    retained = retained,
    
    retained_total =
      sum(retained)
  )
}


# Model fit -----

fit_m3 <- function(
    Y_list,
    p_lags,
    chains = 4L,
    burnin = 500L,
    draws = 1000L,
    thin = 1L,
    seed = 991L,
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
    gc_id = NULL,
    gq_id = NULL,
    standardize = TRUE,
    selected_eps = 1e-12,
    q_transport_global_probability = 0.25,
    use_c_asis = TRUE,
    c_asis_every = 1L,
    c_asis_slice_width = 1,
    c_asis_maximum_step_out = 20L,
    monitor_beta = NULL,
    monitor_lambda = NULL,
    initial_states = NULL,
    overdispersed_initial_states = TRUE,
    tau_scale = NULL,
    p0 = NULL,
    q_grid = seq(
      -log(1000),
      log(1000),
      length.out = 49L
    ),
    q_prob = NULL,
    sigma_a = 3,
    sigma_b = 2,
    tau_df = 10,
    c_df = 5,
    lambda_df = 3,
    c_scale = 1,
    lambda_scale = 1,
    keep_chain_results = FALSE,
    keep_all = FALSE) {
  
  beta_algorithm <- match.arg(
    beta_algorithm
  )
  
  q_update <- match.arg(
    q_update
  )
  
  
  # Model objects -----
  
  structure <- make_m3_structure(
    Y_list = Y_list,
    p_lags = p_lags,
    gc_id = gc_id,
    gq_id = gq_id,
    standardize = standardize
  )
  
  prior <- make_m3_prior(
    structure = structure,
    tau_scale = tau_scale,
    p0 = p0,
    q_grid = q_grid,
    q_prob = q_prob,
    sigma_a = sigma_a,
    sigma_b = sigma_b,
    tau_df = tau_df,
    c_df = c_df,
    lambda_df = lambda_df,
    c_scale = c_scale,
    lambda_scale = lambda_scale
  )
  
  control <- make_m3_control(
    burnin = burnin,
    draws = draws,
    thin = thin,
    beta_algorithm =
      beta_algorithm,
    q_update =
      q_update,
    q_transport_global_probability =
      q_transport_global_probability,
    use_c_asis =
      use_c_asis,
    c_asis_every =
      c_asis_every,
    c_asis_slice_width =
      c_asis_slice_width,
    c_asis_maximum_step_out =
      c_asis_maximum_step_out,
    keep_all = keep_all,
    verbose = FALSE,
    monitor_beta =
      monitor_beta,
    monitor_lambda =
      monitor_lambda
  )
  
  
  # Initial states -----
  
  if (is.null(initial_states)) {
    
    initial_states <-
      make_m3_initial_states(
        structure = structure,
        prior = prior,
        chains = chains,
        overdispersed =
          overdispersed_initial_states
      )
  }
  
  
  # Chains -----
  
  chain_seeds <- m3_chain_seeds(
    seed = seed,
    chains = chains
  )
  
  chain_results <- lapply(
    seq_len(chains),
    function(chain) {
      
      run_m3_chain(
        structure = structure,
        prior = prior,
        control = control,
        state = initial_states[[chain]],
        seed = chain_seeds[chain]
      )
    }
  )
  
  full_draws <- if (keep_all) {
    lapply(chain_results, `[[`, "full_draws")
  } else {
    NULL
  }
  
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
      "M3 chains resolved different beta algorithms"
    )
  }
  
  
  # Posterior pooling -----
  
  pooled <- pool_m3_chains(
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
    model = "m3",
    
    c_group =
      structure$grouping$c_group,
    
    q_group =
      structure$grouping$q_group,
    
    beta_hat = beta_hat,
    
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
    
    c2_mean =
      pooled$c2_mean,
    
    q_mean =
      pooled$q_mean,
    
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
      s_hat_lag > selected_eps,
    
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
    
    q_update =
      q_update,
    
    dimensions =
      structure$dimensions,
    
    preprocessing =
      structure$preprocessing,
    
    prior = prior,
    control = control,
    
    full_draws =
      full_draws
  )
  
  if (keep_chain_results) {
    
    fit$chain_results <-
      chain_results
  }
  
  fit
}