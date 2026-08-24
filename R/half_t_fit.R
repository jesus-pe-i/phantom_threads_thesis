# Multi-chain fitting and network estimation for the coefficient-wise Half-t BVAR.
#
# Runs and pools Gibbs chains, restores coefficients to original units,
# and converts posterior coefficient means into block-network strengths.


# Chain setup -----

half_t_chain_seeds <- function(
    seed,
    chains) {

  as.integer(
    seed + seq_len(chains) - 1L
  )
}


make_half_t_initial_states <- function(
    structure,
    chains,
    overdispersed = TRUE) {

  if (!overdispersed) {

    return(
      replicate(
        chains,
        make_half_t_state(
          structure
        ),
        simplify = FALSE
      )
    )
  }

  sigma_profile <- c(
    0.25,
    4,
    1,
    2
  )

  lambda_profile <- c(
    0.125,
    8,
    1,
    0.5
  )

  tau_profile_two_groups <- rbind(
    c(0.10, 4),
    c(4, 0.10),
    c(1, 1),
    c(0.25, 2)
  )

  xi_profile_two_groups <- rbind(
    c(0.25, 4),
    c(4, 0.25),
    c(1, 1),
    c(2, 0.5)
  )

  tau_profile_one_group <- c(
    0.10,
    4,
    1,
    0.25
  )

  xi_profile_one_group <- c(
    0.25,
    4,
    1,
    2
  )

  lapply(
    seq_len(chains),
    function(chain) {

      profile <- (
        chain - 1L
      ) %% 4L + 1L

      if (structure$maps$n_tau == 2L) {

        tau2 <- tau_profile_two_groups[
          profile,
        ]

        xi <- xi_profile_two_groups[
          profile,
        ]

      } else {

        tau2 <- tau_profile_one_group[
          profile
        ]

        xi <- xi_profile_one_group[
          profile
        ]
      }

      make_half_t_state(
        structure = structure,

        sigma2 = rep(
          sigma_profile[profile],
          structure$data$N
        ),

        lambda2 = matrix(
          lambda_profile[profile],
          nrow = structure$data$k,
          ncol = structure$data$N
        ),

        tau2 = tau2,
        xi = xi
      )
    }
  )
}


# Chain pooling -----

pool_half_t_chains <- function(
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

    xi_mean =
      weighted_mean(
        "xi_mean"
      ),

    retained = retained,

    retained_total =
      sum(retained)
  )
}


# Model fit -----

fit_half_t <- function(
    Y_list,
    p_lags,
    global_grouping = c(
      "all",
      "self_diagonal"
    ),
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
    standardize = TRUE,
    selected_eps = 1e-12,
    monitor_beta = NULL,
    monitor_lambda = NULL,
    use_asis = NULL,
    asis_every = 1L,
    asis_slice_width = 1,
    asis_maximum_step_out = 100L,
    initial_states = NULL,
    overdispersed_initial_states = TRUE,
    lambda_df = 3,
    lambda_scale = 1,
    tau_df = 10,
    tau_scale = 1,
    sigma_a = 3,
    sigma_b = 2,
    keep_chain_results = FALSE) {

  global_grouping <- match.arg(
    global_grouping
  )

  beta_algorithm <- match.arg(
    beta_algorithm
  )

  if (is.null(use_asis)) {

    use_asis <-
      global_grouping ==
      "self_diagonal"
  }


  # Model objects -----

  structure <- make_half_t_structure(
    Y_list = Y_list,
    p_lags = p_lags,
    global_grouping =
      global_grouping,
    standardize =
      standardize
  )

  prior <- make_half_t_prior(
    lambda_df = lambda_df,
    lambda_scale = lambda_scale,
    tau_df = tau_df,
    tau_scale = tau_scale,
    sigma_a = sigma_a,
    sigma_b = sigma_b
  )

  control <- make_half_t_control(
    burnin = burnin,
    draws = draws,
    thin = thin,
    beta_algorithm =
      beta_algorithm,
    use_asis =
      use_asis,
    asis_every =
      asis_every,
    asis_slice_width =
      asis_slice_width,
    asis_maximum_step_out =
      asis_maximum_step_out,
    keep_all = FALSE,
    monitor_beta =
      monitor_beta,
    monitor_lambda =
      monitor_lambda
  )


  # Initial states -----

  if (is.null(initial_states)) {

    initial_states <-
      make_half_t_initial_states(
        structure = structure,
        chains = chains,
        overdispersed =
          overdispersed_initial_states
      )
  }


  # Chains -----

  chain_seeds <- half_t_chain_seeds(
    seed = seed,
    chains = chains
  )

  chain_results <- lapply(
    seq_len(chains),
    function(chain) {

      run_half_t_chain(
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
      "Half-t chains resolved different beta algorithms"
    )
  }


  # Posterior pooling -----

  pooled <- pool_half_t_chains(
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
    model = "half_t",

    global_grouping =
      global_grouping,

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

    xi_mean =
      pooled$xi_mean,

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

    dimensions =
      structure$dimensions,

    preprocessing =
      structure$preprocessing,

    prior = prior,
    control = control
  )

  if (keep_chain_results) {

    fit$chain_results <-
      chain_results
  }

  fit
}
