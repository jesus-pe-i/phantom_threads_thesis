# R-to-C++ execution layer for the coefficient-wise Half-t BVAR.
#
# Builds MCMC controls and C++ inputs, loads the native sampler,
# and runs one Gibbs chain with posterior and monitored draws.


# C++ loading -----

load_half_t_cpp <- function(
    rebuild = FALSE) {

  Rcpp::sourceCpp(
    "src/half_t_bvar_core.cpp",
    rebuild = rebuild
  )
}


# MCMC control -----

make_half_t_control <- function(
    burnin,
    draws,
    thin = 1L,
    beta_algorithm = c(
      "auto",
      "chol",
      "bhattacharya"
    ),
    use_asis = FALSE,
    asis_every = 1L,
    asis_slice_width = 1,
    asis_maximum_step_out = 100L,
    keep_all = FALSE,
    monitor_beta = NULL,
    monitor_lambda = NULL) {

  beta_algorithm <- match.arg(
    beta_algorithm
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

    monitor_lambda <- matrix(
      integer(0L),
      nrow = 0L,
      ncol = 2L
    )
  }

  monitor_lambda <- as.matrix(
    monitor_lambda
  )

  if (ncol(monitor_lambda) != 2L) {
    stop(
      "monitor_lambda must have two columns"
    )
  }

  storage.mode(monitor_lambda) <-
    "integer"

  colnames(monitor_lambda) <- c(
    "row",
    "equation"
  )

  list(
    burnin = as.integer(burnin),
    draws = as.integer(draws),
    thin = as.integer(thin),

    beta_algorithm =
      beta_algorithm,

    use_asis =
      use_asis,

    asis_every =
      as.integer(
        asis_every
      ),

    asis_slice_width =
      as.numeric(
        asis_slice_width
      ),

    asis_maximum_step_out =
      as.integer(
        asis_maximum_step_out
      ),

    keep_all =
      keep_all,

    monitor_beta =
      monitor_beta,

    monitor_lambda =
      monitor_lambda
  )
}


# C++ input -----

make_half_t_input <- function(
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

half_t_draw_matrix <- function(
    x,
    draws,
    parameters) {

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
      "Unexpected Half-t draw dimensions"
    )
  }

  x
}


half_t_monitor_names <- function(
    structure,
    beta_index,
    lambda_index) {

  N <- structure$data$N

  pair_names <- function(
      index,
      prefix) {

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

    paste0(
      prefix,
      "_l",
      lag,
      "_s",
      sender,
      "_r",
      index[, 2L]
    )
  }

  list(
    beta = pair_names(
      beta_index,
      "beta"
    ),

    lambda2 = pair_names(
      lambda_index,
      "lambda2"
    ),

    nu = pair_names(
      lambda_index,
      "nu"
    )
  )
}


# One chain -----

run_half_t_chain <- function(
    structure,
    prior,
    control,
    state = NULL,
    seed = 991L) {

  if (is.null(state)) {

    state <- make_half_t_state(
      structure = structure
    )
  }

  input <- make_half_t_input(
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
    raw <- half_t_bvar_chain_cpp(
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

  lambda_index <- as.matrix(
    raw$monitor$lambda_index
  )

  draws <- list(
    beta = half_t_draw_matrix(
      raw$monitor$beta,
      retained,
      nrow(beta_index)
    ),

    lambda2 = half_t_draw_matrix(
      raw$monitor$lambda2,
      retained,
      nrow(lambda_index)
    ),

    nu = half_t_draw_matrix(
      raw$monitor$nu,
      retained,
      nrow(lambda_index)
    ),

    sigma2 = half_t_draw_matrix(
      raw$draws$sigma2,
      retained,
      structure$data$N
    ),

    tau2 = half_t_draw_matrix(
      raw$draws$tau2,
      retained,
      structure$maps$n_tau
    ),

    xi = half_t_draw_matrix(
      raw$draws$xi,
      retained,
      structure$maps$n_tau
    )
  )

  monitor_names <- half_t_monitor_names(
    structure = structure,
    beta_index = beta_index,
    lambda_index = lambda_index
  )

  colnames(draws$beta) <-
    monitor_names$beta

  colnames(draws$lambda2) <-
    monitor_names$lambda2

  colnames(draws$nu) <-
    monitor_names$nu


  # ASIS diagnostics -----

  asis <- raw$asis

  asis$enabled <- isTRUE(
    asis$enabled
  )

  asis$every <- as.integer(
    asis$every
  )

  asis$slice_width <- as.numeric(
    asis$slice_width
  )

  asis$maximum_step_out <- as.integer(
    asis$maximum_step_out
  )

  asis$maximum_shrinks <- as.integer(
    asis$maximum_shrinks
  )

  asis$updates <- as.integer(
    asis$updates
  )

  asis$moves <- as.integer(
    asis$moves
  )


  # Output -----

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

    tau2_mean = colMeans(
      draws$tau2
    ),

    xi_mean = colMeans(
      draws$xi
    ),

    draws = draws,

    monitor = list(
      beta = beta_index,
      lambda2 = lambda_index
    ),

    final_state =
      raw$final_state,

    asis = asis,
    timing = raw$timing,

    beta_algorithm =
      raw$info$
        beta_algorithm_resolved,

    elapsed_seconds =
      unname(
        elapsed["elapsed"]
      )
  )
}
