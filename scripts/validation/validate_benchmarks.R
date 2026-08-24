# Permanent validation of the unified seven-model benchmark layer.
#
# Fits all benchmark models on one controlled asymmetric block VAR,
# checks the common interface and frozen model-specific settings,
# compares direct and dispatched fits, verifies deterministic coefficient
# and network reconstruction, and reports cross-model performance.


rm(
  list = ls()
)


# Setup -----

source(
  "R/benchmarks.R"
)

source_benchmark_models()


required_packages <- c(
  "Rcpp",
  "RcppArmadillo",
  "adelie",
  "BigVAR",
  "mclust"
)

package_available <- vapply(
  required_packages,
  requireNamespace,
  logical(1L),
  quietly = TRUE
)

cat(
  "\nPackage availability\n"
)

for (package in required_packages) {
  
  cat(
    sprintf(
      "  %-20s %s\n",
      package,
      if (package_available[package]) {
        "AVAILABLE"
      } else {
        "MISSING"
      }
    )
  )
}

if (!all(package_available)) {
  
  missing_packages <- required_packages[
    !package_available
  ]
  
  stop(
    "Missing packages: ",
    paste(
      missing_packages,
      collapse = ", "
    )
  )
}


RNGkind(
  kind = "Mersenne-Twister",
  normal.kind = "Inversion",
  sample.kind = "Rejection"
)


# Native engines -----

cat(
  "\nCompiling native Bayesian engines...\n"
)

Rcpp::sourceCpp(
  "src/half_t_bvar_core.cpp",
  rebuild = FALSE,
  verbose = FALSE
)

Rcpp::sourceCpp(
  "src/m3_bvar_core.cpp",
  rebuild = FALSE,
  verbose = FALSE
)

Rcpp::sourceCpp(
  "src/gigg_bvar_core.cpp",
  rebuild = FALSE,
  verbose = FALSE
)

cat(
  "Native engines loaded.\n"
)


# Validation helpers -----

set_validation_block <- function(
    A,
    receiver,
    sender,
    block,
    m) {
  
  rows <- (
    (receiver - 1L) *
      m +
      1L
  ):(
    receiver *
      m
  )
  
  columns <- (
    (sender - 1L) *
      m +
      1L
  ):(
    sender *
      m
  )
  
  A[
    rows,
    columns
  ] <- block
  
  A
}


var_companion_radius <- function(
    A_list) {
  
  p_lags <- length(
    A_list
  )
  
  N <- nrow(
    A_list[[1L]]
  )
  
  if (p_lags == 1L) {
    
    companion <- A_list[[1L]]
    
  } else {
    
    top <- do.call(
      cbind,
      A_list
    )
    
    lower_left <- diag(
      N * (
        p_lags -
          1L
      )
    )
    
    lower_right <- matrix(
      0,
      nrow =
        N *
        (
          p_lags -
            1L
        ),
      ncol =
        N
    )
    
    bottom <- cbind(
      lower_left,
      lower_right
    )
    
    companion <- rbind(
      top,
      bottom
    )
  }
  
  max(
    Mod(
      eigen(
        companion,
        only.values = TRUE
      )$values
    )
  )
}


maximum_A_error <- function(
    A_left,
    A_right) {
  
  max(
    vapply(
      seq_along(
        A_left
      ),
      function(lag) {
        
        max(
          abs(
            A_left[[lag]] -
              A_right[[lag]]
          )
        )
      },
      numeric(1L)
    )
  )
}


# Controlled DGP -----

set.seed(
  8601L
)

n_units <- 4L
m <- 2L
p_lags <- 2L
T_obs <- 160L
burn_in <- 250L

N <- n_units *
  m

A1 <- matrix(
  0,
  nrow = N,
  ncol = N
)

A2 <- matrix(
  0,
  nrow = N,
  ncol = N
)


## Self dynamics -----

diag(A1) <- 0.24
diag(A2) <- 0.05


## Lag 1 network -----

A1 <- set_validation_block(
  A = A1,
  receiver = 1L,
  sender = 2L,
  block = matrix(
    c(
      0.18,
      -0.05,
      0.06,
      0.14
    ),
    nrow = m,
    byrow = TRUE
  ),
  m = m
)

A1 <- set_validation_block(
  A = A1,
  receiver = 3L,
  sender = 1L,
  block = matrix(
    c(
      0.15,
      0.04,
      -0.05,
      0.12
    ),
    nrow = m,
    byrow = TRUE
  ),
  m = m
)

A1 <- set_validation_block(
  A = A1,
  receiver = 4L,
  sender = 3L,
  block = matrix(
    c(
      0.13,
      -0.03,
      0.04,
      0.11
    ),
    nrow = m,
    byrow = TRUE
  ),
  m = m
)


## Lag 2 network -----

A2 <- set_validation_block(
  A = A2,
  receiver = 2L,
  sender = 4L,
  block = matrix(
    c(
      0.08,
      0.02,
      -0.02,
      0.07
    ),
    nrow = m,
    byrow = TRUE
  ),
  m = m
)

A2 <- set_validation_block(
  A = A2,
  receiver = 3L,
  sender = 2L,
  block = matrix(
    c(
      0.07,
      -0.02,
      0.03,
      0.06
    ),
    nrow = m,
    byrow = TRUE
  ),
  m = m
)


truth_A_lag <- list(
  A1,
  A2
)

spectral_radius <- var_companion_radius(
  truth_A_lag
)

if (
  !is.finite(
    spectral_radius
  ) ||
  spectral_radius >= 1
) {
  stop(
    "Validation DGP is not stationary"
  )
}


## Simulation -----

Y_full <- matrix(
  0,
  nrow =
    T_obs +
    burn_in,
  ncol =
    N
)

for (t in 3:nrow(Y_full)) {
  
  Y_full[t, ] <- as.numeric(
    A1 %*%
      Y_full[
        t - 1L,
      ] +
      A2 %*%
      Y_full[
        t - 2L,
      ] +
      rnorm(
        N,
        sd = 0.35
      )
  )
}

Y <- tail(
  Y_full,
  T_obs
)

Y_list <- lapply(
  seq_len(
    n_units
  ),
  function(unit) {
    
    columns <- (
      (unit - 1L) *
        m +
        1L
    ):(
      unit *
        m
    )
    
    Y[
      ,
      columns,
      drop = FALSE
    ]
  }
)


# Truth objects -----

truth_beta <- A_list_to_beta(
  truth_A_lag
)

truth_s_lag <- A_list_to_s_lag(
  A_list = truth_A_lag,
  n_units = n_units,
  m = m
)

truth_s_unit_rms <- s_lag_to_unit(
  s_lag = truth_s_lag,
  method = "rms"
)

truth_active <- truth_s_lag > 0

n_truth_blocks <- sum(
  truth_active
)


# Fit settings -----

fit_seed <- 991L

bayesian_chains <- 2L
bayesian_burnin <- 100L
bayesian_draws <- 150L

model_args <- list(
  
  m3 = list(
    chains =
      bayesian_chains,
    
    burnin =
      bayesian_burnin,
    
    draws =
      bayesian_draws,
    
    thin =
      1L,
    
    beta_algorithm =
      "chol"
  ),
  
  gigg = list(
    chains =
      bayesian_chains,
    
    burnin =
      bayesian_burnin,
    
    draws =
      bayesian_draws,
    
    thin =
      1L,
    
    beta_algorithm =
      "chol"
  ),
  
  half_t = list(
    chains =
      bayesian_chains,
    
    burnin =
      bayesian_burnin,
    
    draws =
      bayesian_draws,
    
    thin =
      1L,
    
    beta_algorithm =
      "chol"
  ),
  
  adelie_gen = list(
    n_threads =
      1L
  ),
  
  bigvar_hlag =
    list(),
  
  mar = list(
    n_starts =
      3L,
    
    max_iter =
      100L
  ),
  
  nirvar = list(
    k =
      2L
  )
)


# Unified fits -----

cat(
  "\nRunning unified benchmark fits...\n\n"
)

fits <- fit_benchmarks(
  models =
    benchmark_models,
  Y_list =
    Y_list,
  p_lags =
    p_lags,
  seed =
    fit_seed,
  model_args =
    model_args,
  verbose =
    TRUE
)


# Direct fits -----

cat(
  "\nRunning direct fits for dispatcher comparison...\n"
)

set.seed(
  fit_seed
)

direct_m3 <- fit_m3(
  Y_list =
    Y_list,
  p_lags =
    p_lags,
  chains =
    bayesian_chains,
  burnin =
    bayesian_burnin,
  draws =
    bayesian_draws,
  thin =
    1L,
  seed =
    fit_seed,
  beta_algorithm =
    "chol"
)


set.seed(
  fit_seed
)

direct_gigg <- fit_gigg(
  Y_list =
    Y_list,
  p_lags =
    p_lags,
  chains =
    bayesian_chains,
  burnin =
    bayesian_burnin,
  draws =
    bayesian_draws,
  thin =
    1L,
  seed =
    fit_seed,
  beta_algorithm =
    "chol"
)


set.seed(
  fit_seed
)

direct_half_t <- fit_half_t(
  Y_list =
    Y_list,
  p_lags =
    p_lags,
  global_grouping =
    "self_diagonal",
  chains =
    bayesian_chains,
  burnin =
    bayesian_burnin,
  draws =
    bayesian_draws,
  thin =
    1L,
  seed =
    fit_seed,
  beta_algorithm =
    "chol"
)


set.seed(
  fit_seed
)

direct_adelie <- fit_adelie(
  Y_list =
    Y_list,
  p_lags =
    p_lags,
  n_threads =
    1L,
  seed =
    fit_seed
)


set.seed(
  fit_seed
)

direct_bigvar <- fit_bigvar(
  Y_list =
    Y_list,
  p_lags =
    p_lags
)


set.seed(
  fit_seed
)

direct_mar <- fit_mar(
  Y_list =
    Y_list,
  p_lags =
    p_lags,
  n_starts =
    3L,
  seed =
    fit_seed,
  max_iter =
    100L
)


set.seed(
  fit_seed
)

direct_nirvar <- fit_nirvar(
  Y_list =
    Y_list,
  p_lags =
    p_lags,
  k =
    2L,
  seed =
    fit_seed
)


direct_fits <- list(
  m3 =
    direct_m3,
  
  gigg =
    direct_gigg,
  
  half_t =
    direct_half_t,
  
  adelie_gen =
    direct_adelie,
  
  bigvar_hlag =
    direct_bigvar,
  
  mar =
    direct_mar,
  
  nirvar =
    direct_nirvar
)


# Registry validation -----

expected_models <- c(
  "m3",
  "gigg",
  "half_t",
  "adelie_gen",
  "bigvar_hlag",
  "mar",
  "nirvar"
)

expected_backends <- c(
  m3 =
    "RcppArmadillo",
  
  gigg =
    "RcppArmadillo",
  
  half_t =
    "RcppArmadillo",
  
  adelie_gen =
    "adelie",
  
  bigvar_hlag =
    "BigVAR",
  
  mar =
    "mar",
  
  nirvar =
    "nirvar"
)


registry_checks <- c(
  model_order =
    identical(
      benchmark_models,
      expected_models
    ),
  
  returned_models =
    identical(
      names(fits),
      expected_models
    ),
  
  registry_models =
    identical(
      names(
        benchmark_registry
      ),
      expected_models
    )
)

for (model in expected_models) {
  
  registry_checks[
    paste0(
      "backend_",
      model
    )
  ] <- identical(
    fits[[model]]$backend,
    expected_backends[[model]]
  )
}


# Common contract validation -----

contract_checks <- logical(
  length(
    expected_models
  )
)

names(
  contract_checks
) <- expected_models

for (model in expected_models) {
  
  result <- tryCatch(
    {
      
      validate_benchmark_fit(
        fit =
          fits[[model]],
        model =
          model,
        Y_list =
          Y_list,
        p_lags =
          p_lags
      )
      
      TRUE
    },
    error = function(error) {
      
      cat(
        sprintf(
          "Contract failure for %s: %s\n",
          model,
          conditionMessage(
            error
          )
        )
      )
      
      FALSE
    }
  )
  
  contract_checks[model] <- result
}


# Reconstruction validation -----

A_reconstruction_error <- numeric(
  length(
    expected_models
  )
)

network_reconstruction_error <- numeric(
  length(
    expected_models
  )
)

unit_max_reconstruction_error <- numeric(
  length(
    expected_models
  )
)

unit_rms_reconstruction_error <- numeric(
  length(
    expected_models
  )
)

score_error <- numeric(
  length(
    expected_models
  )
)

names(
  A_reconstruction_error
) <- expected_models

names(
  network_reconstruction_error
) <- expected_models

names(
  unit_max_reconstruction_error
) <- expected_models

names(
  unit_rms_reconstruction_error
) <- expected_models

names(
  score_error
) <- expected_models


for (model in expected_models) {
  
  fit <- fits[[model]]
  
  reconstructed_A <- beta_to_A_list(
    beta =
      fit$beta_hat,
    p_lags =
      p_lags
  )
  
  reconstructed_s <- A_list_to_s_lag(
    A_list =
      reconstructed_A,
    n_units =
      n_units,
    m =
      m
  )
  
  reconstructed_max <- s_lag_to_unit(
    s_lag =
      reconstructed_s,
    method =
      "max"
  )
  
  reconstructed_rms <- s_lag_to_unit(
    s_lag =
      reconstructed_s,
    method =
      "rms"
  )
  
  A_reconstruction_error[model] <-
    maximum_A_error(
      reconstructed_A,
      fit$A_hat_lag
    )
  
  network_reconstruction_error[model] <-
    max(
      abs(
        reconstructed_s -
          fit$s_hat_lag
      )
    )
  
  unit_max_reconstruction_error[model] <-
    max(
      abs(
        reconstructed_max -
          fit$s_hat_unit_max
      )
    )
  
  unit_rms_reconstruction_error[model] <-
    max(
      abs(
        reconstructed_rms -
          fit$s_hat_unit_rms
      )
    )
  
  score_error[model] <-
    max(
      abs(
        fit$score_lag -
          fit$s_hat_lag
      )
    )
}


reconstruction_tolerance <- 1e-10

reconstruction_checks <- c(
  A =
    all(
      A_reconstruction_error <
        reconstruction_tolerance
    ),
  
  lag_strength =
    all(
      network_reconstruction_error <
        reconstruction_tolerance
    ),
  
  unit_max =
    all(
      unit_max_reconstruction_error <
        reconstruction_tolerance
    ),
  
  unit_rms =
    all(
      unit_rms_reconstruction_error <
        reconstruction_tolerance
    ),
  
  score =
    all(
      score_error <
        reconstruction_tolerance
    )
)


# Direct-dispatch equivalence -----

direct_beta_error <- numeric(
  length(
    expected_models
  )
)

direct_A_error <- numeric(
  length(
    expected_models
  )
)

direct_network_error <- numeric(
  length(
    expected_models
  )
)

direct_selection_same <- logical(
  length(
    expected_models
  )
)

names(
  direct_beta_error
) <- expected_models

names(
  direct_A_error
) <- expected_models

names(
  direct_network_error
) <- expected_models

names(
  direct_selection_same
) <- expected_models


for (model in expected_models) {
  
  dispatched <- fits[[model]]
  direct <- direct_fits[[model]]
  
  direct_beta_error[model] <-
    max(
      abs(
        dispatched$beta_hat -
          direct$beta_hat
      )
    )
  
  direct_A_error[model] <-
    maximum_A_error(
      dispatched$A_hat_lag,
      direct$A_hat_lag
    )
  
  direct_network_error[model] <-
    max(
      abs(
        dispatched$s_hat_lag -
          direct$s_hat_lag
      )
    )
  
  direct_selection_same[model] <-
    identical(
      dispatched$selected_lag,
      direct$selected_lag
    )
}


direct_tolerance <- 1e-10

direct_checks <- c(
  beta =
    all(
      direct_beta_error <
        direct_tolerance
    ),
  
  A =
    all(
      direct_A_error <
        direct_tolerance
    ),
  
  network =
    all(
      direct_network_error <
        direct_tolerance
    ),
  
  selection =
    all(
      direct_selection_same
    )
)


# Bayesian invariants -----

expected_retained <-
  bayesian_chains *
  bayesian_draws


## Half-t -----

half_t_checks <- c(
  grouping =
    identical(
      fits$half_t$global_grouping,
      "self_diagonal"
    ),
  
  tau_df =
    identical(
      as.numeric(
        fits$half_t$prior$tau_df
      ),
      10
    ),
  
  lambda_df =
    identical(
      as.numeric(
        fits$half_t$prior$lambda_df
      ),
      3
    ),
  
  tau_scale =
    identical(
      as.numeric(
        fits$half_t$prior$tau_scale
      ),
      1
    ),
  
  standardize =
    identical(
      fits$half_t$preprocessing$standardize,
      TRUE
    ),
  
  asis =
    identical(
      fits$half_t$control$use_asis,
      TRUE
    ),
  
  beta_algorithm =
    identical(
      fits$half_t$beta_algorithm_resolved,
      "chol"
    ),
  
  retained =
    identical(
      as.integer(
        fits$half_t$retained_total
      ),
      as.integer(
        expected_retained
      )
    )
)


## M3 -----

m3_checks <- c(
  c_group =
    identical(
      fits$m3$c_group,
      "self_foreign"
    ),
  
  q_group =
    identical(
      fits$m3$q_group,
      "self_foreign"
    ),
  
  tau_df =
    identical(
      as.numeric(
        fits$m3$prior$tau_df
      ),
      10
    ),
  
  c_df =
    identical(
      as.numeric(
        fits$m3$prior$c_df
      ),
      5
    ),
  
  lambda_df =
    identical(
      as.numeric(
        fits$m3$prior$lambda_df
      ),
      3
    ),
  
  standardize =
    identical(
      fits$m3$preprocessing$standardize,
      TRUE
    ),
  
  c_asis =
    identical(
      fits$m3$control$use_c_asis,
      TRUE
    ),
  
  beta_algorithm =
    identical(
      fits$m3$beta_algorithm_resolved,
      "chol"
    ),
  
  retained =
    identical(
      as.integer(
        fits$m3$retained_total
      ),
      as.integer(
        expected_retained
      )
    )
)


## GIGG -----

gigg_checks <- c(
  tau_df =
    identical(
      as.numeric(
        fits$gigg$prior$tau_df
      ),
      5
    ),
  
  gamma_shape =
    identical(
      as.numeric(
        fits$gigg$prior$gamma_shape
      ),
      0.5
    ),
  
  gamma_rate =
    identical(
      as.numeric(
        fits$gigg$prior$gamma_rate
      ),
      1
    ),
  
  lambda_shape =
    identical(
      as.numeric(
        fits$gigg$prior$lambda_shape
      ),
      2.5
    ),
  
  lambda_scale =
    identical(
      as.numeric(
        fits$gigg$prior$lambda_scale
      ),
      1
    ),
  
  standardize =
    identical(
      fits$gigg$preprocessing$standardize,
      TRUE
    ),
  
  global_tau =
    identical(
      fits$gigg$control$use_global_tau,
      TRUE
    ),
  
  asis =
    identical(
      fits$gigg$control$use_asis,
      TRUE
    ),
  
  beta_algorithm =
    identical(
      fits$gigg$beta_algorithm_resolved,
      "chol"
    ),
  
  retained =
    identical(
      as.integer(
        fits$gigg$retained_total
      ),
      as.integer(
        expected_retained
      )
    )
)


# External-model invariants -----

## Adelie -----

adelie_active_total <- sum(
  lengths(
    fits$adelie_gen$active_groups
  )
)

adelie_selected_total <- sum(
  fits$adelie_gen$selected_lag
)

adelie_checks <- c(
  alpha =
    abs(
      fits$adelie_gen$alpha -
        0.5
    ) <
    1e-12,
  
  lambda_rule =
    identical(
      fits$adelie_gen$lambda_rule,
      "lambda.min"
    ),
  
  folds =
    identical(
      fits$adelie_gen$n_folds,
      5L
    ),
  
  blocked_folds =
    all(
      diff(
        fits$adelie_gen$foldid
      ) >= 0
    ),
  
  fold_levels =
    identical(
      as.integer(
        sort(
          unique(
            fits$adelie_gen$foldid
          )
        )
      ),
      seq_len(
        fits$adelie_gen$n_folds
      )
    ),
  
  standardize =
    identical(
      fits$adelie_gen$preprocessing$adelie_standardize,
      TRUE
    ),
  
  no_intercept =
    identical(
      fits$adelie_gen$control$intercept,
      FALSE
    ),
  
  active_selection =
    identical(
      as.integer(
        adelie_active_total
      ),
      as.integer(
        adelie_selected_total
      )
    )
)


## BigVAR -----

bigvar_checks <- c(
  structure =
    identical(
      fits$bigvar_hlag$struct,
      "HLAGOO"
    ),
  
  rolling_cv =
    identical(
      fits$bigvar_hlag$cv,
      "Rolling"
    ),
  
  no_IC =
    identical(
      fits$bigvar_hlag$IC,
      FALSE
    ),
  
  T1 =
    identical(
      fits$bigvar_hlag$T1,
      as.integer(
        floor(
          T_obs /
            3
        )
      )
    ),
  
  T2 =
    identical(
      fits$bigvar_hlag$T2,
      as.integer(
        floor(
          2 *
            T_obs /
            3
        )
      )
    ),
  
  lambda_finite =
    is.finite(
      fits$bigvar_hlag$optimal_lambda
    ),
  
  lambda_index =
    fits$bigvar_hlag$lambda_index >=
    1L &&
    fits$bigvar_hlag$lambda_index <=
    length(
      fits$bigvar_hlag$lambda_grid
    )
)


## MAR -----

mar_B_norm <- sqrt(
  sum(
    fits$mar$mar_B^2
  )
)

mar_checks <- c(
  B_normalization =
    abs(
      mar_B_norm -
        sqrt(m)
    ) <
    1e-10,
  
  multistarts =
    identical(
      fits$mar$n_starts,
      3L
    ),
  
  best_start =
    fits$mar$best_start >=
    1L &&
    fits$mar$best_start <=
    fits$mar$n_starts,
  
  finite_rss =
    is.finite(
      fits$mar$final_rss
    ),
  
  improved_rss =
    fits$mar$final_rss <=
    fits$mar$initial_rss +
    1e-10,
  
  converged =
    isTRUE(
      fits$mar$converged
    )
)


## NIRVAR -----

nirvar_checks <- c(
  k_requested =
    identical(
      as.integer(
        fits$nirvar$k_requested
      ),
      2L
    ),
  
  k_used =
    fits$nirvar$k_used >=
    1L &&
    fits$nirvar$k_used <=
    fits$nirvar$k_requested,
  
  dimension =
    fits$nirvar$d >=
    1L &&
    fits$nirvar$d <=
    n_units -
    1L,
  
  restriction_dimensions =
    identical(
      dim(
        fits$nirvar$restriction
      ),
      as.integer(
        c(
          m,
          m,
          n_units,
          n_units
        )
      )
    ),
  
  labels_dimensions =
    identical(
      dim(
        fits$nirvar$labels
      ),
      as.integer(
        c(
          m,
          n_units
        )
      )
    ),
  
  equation_count =
    nrow(
      fits$nirvar$equation_diagnostics
    ) ==
    N,
  
  finite_covariance =
    all(
      is.finite(
        fits$nirvar$Sigma_hat
      )
    ),
  
  cluster_method =
    fits$nirvar$cluster_method %in%
    c(
      "gmm",
      "kmeans_fallback",
      "single_cluster",
      "single_cluster_distinct_rows",
      "single_cluster_fallback"
    )
)


# Seed and runtime invariants -----

seed_checks <- vapply(
  fits,
  function(fit) {
    
    identical(
      fit$fit_seed,
      fit_seed
    )
  },
  logical(1L)
)

runtime_checks <- vapply(
  fits,
  function(fit) {
    
    length(
      fit$runtime_seconds
    ) ==
      1L &&
      is.finite(
        fit$runtime_seconds
      ) &&
      fit$runtime_seconds >=
      0
  },
  logical(1L)
)


# Truth comparison -----

comparison <- data.frame(
  model =
    expected_models,
  
  backend =
    unname(
      expected_backends[
        expected_models
      ]
    ),
  
  runtime_seconds =
    NA_real_,
  
  selected_blocks =
    NA_integer_,
  
  beta_rmse =
    NA_real_,
  
  lag_strength_rmse =
    NA_real_,
  
  unit_rms_rmse =
    NA_real_,
  
  active_mean_score =
    NA_real_,
  
  inactive_mean_score =
    NA_real_,
  
  score_separation =
    NA_real_,
  
  spearman_truth =
    NA_real_,
  
  stringsAsFactors =
    FALSE
)


for (row in seq_len(
  nrow(
    comparison
  )
)) {
  
  model <- comparison$model[row]
  fit <- fits[[model]]
  
  comparison$runtime_seconds[row] <-
    fit$runtime_seconds
  
  comparison$selected_blocks[row] <-
    sum(
      fit$selected_lag
    )
  
  comparison$beta_rmse[row] <-
    sqrt(
      mean(
        (
          fit$beta_hat -
            truth_beta
        )^2
      )
    )
  
  comparison$lag_strength_rmse[row] <-
    sqrt(
      mean(
        (
          fit$s_hat_lag -
            truth_s_lag
        )^2
      )
    )
  
  comparison$unit_rms_rmse[row] <-
    sqrt(
      mean(
        (
          fit$s_hat_unit_rms -
            truth_s_unit_rms
        )^2
      )
    )
  
  comparison$active_mean_score[row] <-
    mean(
      fit$score_lag[
        truth_active
      ]
    )
  
  comparison$inactive_mean_score[row] <-
    mean(
      fit$score_lag[
        !truth_active
      ]
    )
  
  comparison$score_separation[row] <-
    comparison$active_mean_score[row] -
    comparison$inactive_mean_score[row]
  
  comparison$spearman_truth[row] <-
    suppressWarnings(
      stats::cor(
        as.vector(
          fit$score_lag
        ),
        as.vector(
          truth_s_lag
        ),
        method =
          "spearman"
      )
    )
}


comparison_finite <- all(
  is.finite(
    comparison$runtime_seconds
  )
) &&
  all(
    is.finite(
      comparison$beta_rmse
    )
  ) &&
  all(
    is.finite(
      comparison$lag_strength_rmse
    )
  ) &&
  all(
    is.finite(
      comparison$unit_rms_rmse
    )
  ) &&
  all(
    is.finite(
      comparison$active_mean_score
    )
  ) &&
  all(
    is.finite(
      comparison$inactive_mean_score
    )
  ) &&
  all(
    is.finite(
      comparison$score_separation
    )
  ) &&
  all(
    is.finite(
      comparison$spearman_truth
    )
  )


# Combined checks -----

checks <- c(
  
  registry_checks,
  
  setNames(
    contract_checks,
    paste0(
      "contract_",
      names(
        contract_checks
      )
    )
  ),
  
  setNames(
    reconstruction_checks,
    paste0(
      "reconstruction_",
      names(
        reconstruction_checks
      )
    )
  ),
  
  setNames(
    direct_checks,
    paste0(
      "direct_",
      names(
        direct_checks
      )
    )
  ),
  
  setNames(
    half_t_checks,
    paste0(
      "half_t_",
      names(
        half_t_checks
      )
    )
  ),
  
  setNames(
    m3_checks,
    paste0(
      "m3_",
      names(
        m3_checks
      )
    )
  ),
  
  setNames(
    gigg_checks,
    paste0(
      "gigg_",
      names(
        gigg_checks
      )
    )
  ),
  
  setNames(
    adelie_checks,
    paste0(
      "adelie_",
      names(
        adelie_checks
      )
    )
  ),
  
  setNames(
    bigvar_checks,
    paste0(
      "bigvar_",
      names(
        bigvar_checks
      )
    )
  ),
  
  setNames(
    mar_checks,
    paste0(
      "mar_",
      names(
        mar_checks
      )
    )
  ),
  
  setNames(
    nirvar_checks,
    paste0(
      "nirvar_",
      names(
        nirvar_checks
      )
    )
  ),
  
  seeds =
    all(
      seed_checks
    ),
  
  runtimes =
    all(
      runtime_checks
    ),
  
  comparison_finite =
    comparison_finite
)


# Console report -----

separator <- paste0(
  rep(
    "=",
    76L
  ),
  collapse = ""
)


cat(
  "\n",
  separator,
  "\nFULL BENCHMARK VALIDATION\n",
  separator,
  "\n\n",
  sep = ""
)


cat(
  "Controlled DGP\n",
  sprintf(
    "  Units / variables / lags ............ %d / %d / %d\n",
    n_units,
    m,
    p_lags
  ),
  sprintf(
    "  Observations ........................ %d\n",
    T_obs
  ),
  sprintf(
    "  Total series ........................ %d\n",
    N
  ),
  sprintf(
    "  True active lag blocks .............. %d / %d\n",
    n_truth_blocks,
    length(
      truth_active
    )
  ),
  sprintf(
    "  Companion spectral radius ........... %.6f\n",
    spectral_radius
  ),
  "\n",
  sep = ""
)


cat(
  "Benchmark registry\n"
)

for (model in expected_models) {
  
  cat(
    sprintf(
      "  %-14s -> %s\n",
      model,
      fits[[model]]$backend
    )
  )
}


cat(
  "\nDirect versus dispatched\n"
)

for (model in expected_models) {
  
  cat(
    sprintf(
      paste0(
        "  %-14s",
        " beta %.3e",
        " | A %.3e",
        " | network %.3e",
        " | selection %s\n"
      ),
      model,
      direct_beta_error[model],
      direct_A_error[model],
      direct_network_error[model],
      if (direct_selection_same[model]) {
        "MATCH"
      } else {
        "DIFF"
      }
    )
  )
}


cat(
  "\nReconstruction errors\n"
)

for (model in expected_models) {
  
  cat(
    sprintf(
      paste0(
        "  %-14s",
        " A %.3e",
        " | lag %.3e",
        " | max %.3e",
        " | rms %.3e\n"
      ),
      model,
      A_reconstruction_error[model],
      network_reconstruction_error[model],
      unit_max_reconstruction_error[model],
      unit_rms_reconstruction_error[model]
    )
  )
}


cat(
  "\nBayesian defaults\n",
  sprintf(
    "  Half-t df (tau / lambda) ............ %.0f / %.0f\n",
    fits$half_t$prior$tau_df,
    fits$half_t$prior$lambda_df
  ),
  sprintf(
    "  Half-t grouping ..................... %s\n",
    fits$half_t$global_grouping
  ),
  sprintf(
    "  M3 df (tau / c / lambda) ............ %.0f / %.0f / %.0f\n",
    fits$m3$prior$tau_df,
    fits$m3$prior$c_df,
    fits$m3$prior$lambda_df
  ),
  sprintf(
    "  M3 groups (c / q) ................... %s / %s\n",
    fits$m3$c_group,
    fits$m3$q_group
  ),
  sprintf(
    "  GIGG tau df ........................ %.0f\n",
    fits$gigg$prior$tau_df
  ),
  sprintf(
    "  GIGG gamma shape / rate ............. %.1f / %.1f\n",
    fits$gigg$prior$gamma_shape,
    fits$gigg$prior$gamma_rate
  ),
  sprintf(
    "  GIGG lambda shape / scale ........... %.1f / %.1f\n",
    fits$gigg$prior$lambda_shape,
    fits$gigg$prior$lambda_scale
  ),
  sprintf(
    "  Retained draws per Bayesian fit ..... %d\n",
    expected_retained
  ),
  "\n",
  sep = ""
)


cat(
  "External-model diagnostics\n",
  sprintf(
    "  Adelie active / selected blocks ..... %d / %d\n",
    adelie_active_total,
    adelie_selected_total
  ),
  sprintf(
    "  BigVAR lambda index ................. %d / %d\n",
    fits$bigvar_hlag$lambda_index,
    length(
      fits$bigvar_hlag$lambda_grid
    )
  ),
  sprintf(
    "  BigVAR lambda boundary .............. %s\n",
    fits$bigvar_hlag$lambda_boundary
  ),
  sprintf(
    "  MAR best start ...................... %d / %d\n",
    fits$mar$best_start,
    fits$mar$n_starts
  ),
  sprintf(
    "  MAR converged ....................... %s\n",
    fits$mar$converged
  ),
  sprintf(
    "  MAR ||B||_F / target ................ %.6f / %.6f\n",
    mar_B_norm,
    sqrt(m)
  ),
  sprintf(
    "  NIRVAR d ............................ %d\n",
    fits$nirvar$d
  ),
  sprintf(
    "  NIRVAR clusters ..................... %d / %d requested\n",
    fits$nirvar$k_used,
    fits$nirvar$k_requested
  ),
  sprintf(
    "  NIRVAR clustering method ............ %s\n",
    fits$nirvar$cluster_method
  ),
  "\n",
  sep = ""
)


cat(
  "Cross-model comparison\n\n"
)

comparison_print <- comparison

comparison_print$runtime_seconds <-
  round(
    comparison_print$runtime_seconds,
    3L
  )

comparison_print$beta_rmse <-
  round(
    comparison_print$beta_rmse,
    5L
  )

comparison_print$lag_strength_rmse <-
  round(
    comparison_print$lag_strength_rmse,
    5L
  )

comparison_print$unit_rms_rmse <-
  round(
    comparison_print$unit_rms_rmse,
    5L
  )

comparison_print$active_mean_score <-
  round(
    comparison_print$active_mean_score,
    5L
  )

comparison_print$inactive_mean_score <-
  round(
    comparison_print$inactive_mean_score,
    5L
  )

comparison_print$score_separation <-
  round(
    comparison_print$score_separation,
    5L
  )

comparison_print$spearman_truth <-
  round(
    comparison_print$spearman_truth,
    4L
  )

print(
  comparison_print,
  row.names = FALSE
)


cat(
  "\nChecks\n"
)

for (name in names(checks)) {
  
  cat(
    sprintf(
      "  %-48s %s\n",
      name,
      if (isTRUE(
        checks[name]
      )) {
        "PASSED"
      } else {
        "FAILED"
      }
    )
  )
}


failed_checks <- names(
  checks
)[
  !checks
]


cat(
  "\n",
  separator,
  "\n",
  if (length(
    failed_checks
  ) == 0L) {
    "ALL BENCHMARK VALIDATION CHECKS PASSED"
  } else {
    "BENCHMARK VALIDATION FAILED"
  },
  "\n",
  separator,
  "\n",
  sep = ""
)


if (length(
  failed_checks
) > 0L) {
  
  cat(
    "\nFailed checks:\n"
  )
  
  for (name in failed_checks) {
    
    cat(
      "  - ",
      name,
      "\n",
      sep = ""
    )
  }
  
  stop(
    "Benchmark validation failed after completing the full report"
  )
}