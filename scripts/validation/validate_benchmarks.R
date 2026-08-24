# Permanent integration validation for the seven-model benchmark layer.
#
# Checks model loading and dispatch, the shared fit contract, deterministic
# reconstruction, nested benchmark records, benchmark metrics, frozen Bayesian
# settings, direct-dispatch equivalence, seeds, runtimes, and model performance.


rm(
  list = ls()
)


# Helpers -----

checks <- logical(0L)


check <- function(
    label,
    condition) {
  
  checks[label] <<- isTRUE(
    condition
  )
  
  invisible(
    checks[label]
  )
}


same <- function(
    x,
    y,
    tolerance = 1e-12) {
  
  isTRUE(
    all.equal(
      x,
      y,
      tolerance = tolerance,
      check.attributes = FALSE
    )
  )
}


finite <- function(
    x) {
  
  all(
    is.finite(
      as.numeric(x)
    )
  )
}


maximum_A_error <- function(
    A_left,
    A_right) {
  
  max(
    vapply(
      seq_along(A_left),
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


set_validation_block <- function(
    A,
    receiver,
    sender,
    block,
    m) {
  
  rows <- ((receiver - 1L) * m + 1L):(receiver * m)
  columns <- ((sender - 1L) * m + 1L):(sender * m)
  
  A[
    rows,
    columns
  ] <- block
  
  A
}


# Setup -----

source(
  "R/benchmarks.R"
)

source_benchmark_models()

source(
  "R/simulation.R"
)

source(
  "R/benchmark_records.R"
)

source(
  "R/benchmark_metrics.R"
)


required_packages <- c(
  "Rcpp",
  "RcppArmadillo",
  "adelie",
  "BigVAR"
)


package_available <- vapply(
  required_packages,
  requireNamespace,
  logical(1L),
  quietly = TRUE
)


if (!all(package_available)) {
  
  stop(
    "Missing packages: ",
    paste(
      required_packages[!package_available],
      collapse = ", "
    ),
    call. = FALSE
  )
}


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
  m3 = "RcppArmadillo",
  gigg = "RcppArmadillo",
  half_t = "RcppArmadillo",
  adelie_gen = "adelie",
  bigvar_hlag = "BigVAR",
  mar = "mar",
  nirvar = "nirvar"
)


expected_fit_functions <- vapply(
  benchmark_registry,
  function(specification) {
    specification$fit
  },
  character(1L)
)


check(
  "Benchmark registry",
  identical(
    benchmark_models,
    expected_models
  ) &&
    identical(
      names(benchmark_registry),
      expected_models
    )
)


check(
  "Fit functions loaded",
  all(
    vapply(
      expected_fit_functions,
      exists,
      logical(1L),
      mode = "function",
      inherits = TRUE
    )
  )
)


# Native engines -----

load_half_t_cpp(
  rebuild = FALSE
)

load_m3_cpp(
  rebuild = FALSE
)

load_gigg_cpp(
  rebuild = FALSE
)


check(
  "Native engines loaded",
  exists(
    "half_t_bvar_chain_cpp",
    mode = "function"
  ) &&
    exists(
      "m3_bvar_chain_cpp",
      mode = "function"
    ) &&
    exists(
      "gigg_bvar_chain_cpp",
      mode = "function"
    )
)


# Controlled DGP -----

n_units <- 4L
m <- 2L
p_lags <- 2L
T_obs <- 120L
burn_in <- 150L

N <- n_units * m
k <- N * p_lags


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


diag(A1) <- 0.24
diag(A2) <- 0.05


A1 <- set_validation_block(
  A = A1,
  receiver = 1L,
  sender = 2L,
  block = matrix(
    c(
      0.18, -0.05,
      0.06, 0.14
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
      0.15, 0.04,
      -0.05, 0.12
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
      0.13, -0.03,
      0.04, 0.11
    ),
    nrow = m,
    byrow = TRUE
  ),
  m = m
)


A2 <- set_validation_block(
  A = A2,
  receiver = 2L,
  sender = 4L,
  block = matrix(
    c(
      0.08, 0.02,
      -0.02, 0.07
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
      0.07, -0.02,
      0.03, 0.06
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


stability <- check_A_radius(
  truth_A_lag
)


check(
  "Controlled DGP stable",
  finite(
    stability$radius
  ) &&
    stability$radius < 1
)


simulation <- simulate_var_from_A(
  A_list = truth_A_lag,
  Sigma = make_sigma_homoscedastic(
    N = N,
    sigma = 0.35
  ),
  T_obs = T_obs,
  burn_in = burn_in,
  seed = 8601L,
  n_units = n_units,
  m = m
)


Y_list <- simulation$Y_list


check(
  "Controlled data dimensions",
  length(Y_list) == n_units &&
    all(
      vapply(
        Y_list,
        function(Y) {
          
          identical(
            dim(Y),
            c(
              T_obs,
              m
            )
          ) &&
            finite(Y)
        },
        logical(1L)
      )
    )
)


truth_beta <- A_list_to_beta(
  truth_A_lag
)


truth <- make_benchmark_truth(
  A_list = truth_A_lag,
  n_units = n_units,
  m = m
)


# Benchmark records -----

validation_root <- tempfile(
  "benchmark_validation_"
)

dgp_id <- "validation_n4_m2_p2"

dgp_path <- file.path(
  validation_root,
  dgp_id
)


dir.create(
  dgp_path,
  recursive = TRUE
)


saveRDS(
  truth_A_lag,
  file.path(
    dgp_path,
    "A_list.rds"
  )
)


loaded_dgp <- load_benchmark_dgp(
  dgp_id = dgp_id,
  dgp_root = validation_root
)


check(
  "DGP loading",
  loaded_dgp$n_units == n_units &&
    loaded_dgp$m == m &&
    loaded_dgp$N == N &&
    loaded_dgp$p_lags == p_lags &&
    maximum_A_error(
      loaded_dgp$A_list,
      truth_A_lag
    ) < 1e-12
)


record_T_grid <- c(
  60L,
  90L
)

record_seeds <- c(
  7001L,
  7002L
)


records <- make_benchmark_records(
  dgp_id = dgp_id,
  T_grid = record_T_grid,
  seeds = record_seeds,
  burn_in = 100L,
  sigma = 0.35,
  dgp_root = validation_root
)


expected_record_names <- c(
  paste0(
    dgp_id,
    "_seed_7001_T_60"
  ),
  paste0(
    dgp_id,
    "_seed_7001_T_90"
  ),
  paste0(
    dgp_id,
    "_seed_7002_T_60"
  ),
  paste0(
    dgp_id,
    "_seed_7002_T_90"
  )
)


check(
  "Benchmark record identities",
  length(records) == 4L &&
    identical(
      names(records),
      expected_record_names
    )
)


record_truth_same <- all(
  vapply(
    records,
    function(record) {
      
      identical(
        record$truth,
        records[[1L]]$truth
      )
    },
    logical(1L)
  )
)


check(
  "Frozen record truth",
  record_truth_same &&
    same(
      records[[1L]]$truth$s_lag,
      truth$s_lag
    ) &&
    identical(
      records[[1L]]$truth$G_lag,
      truth$G_lag
    )
)


nested_checks <- vapply(
  record_seeds,
  function(seed) {
    
    short_id <- paste0(
      dgp_id,
      "_seed_",
      seed,
      "_T_60"
    )
    
    long_id <- paste0(
      dgp_id,
      "_seed_",
      seed,
      "_T_90"
    )
    
    short_record <- records[[short_id]]
    long_record <- records[[long_id]]
    
    all(
      vapply(
        seq_len(n_units),
        function(unit) {
          
          identical(
            short_record$data$Y_list[[unit]],
            long_record$data$Y_list[[unit]][
              seq_len(60L),
              ,
              drop = FALSE
            ]
          )
        },
        logical(1L)
      )
    )
  },
  logical(1L)
)


check(
  "Nested sample paths",
  all(
    nested_checks
  )
)


different_seed_paths <- !identical(
  records[[expected_record_names[2L]]]$data$Y_list,
  records[[expected_record_names[4L]]]$data$Y_list
)


check(
  "Independent record seeds",
  different_seed_paths
)


# Metric helpers -----

lag_mask <- benchmark_network_mask(
  n_units = n_units,
  p_lags = p_lags
)

unit_mask <- benchmark_network_mask(
  n_units = n_units
)


check(
  "Network masks",
  identical(
    dim(lag_mask),
    c(
      n_units,
      n_units,
      p_lags
    )
  ) &&
    identical(
      dim(unit_mask),
      c(
        n_units,
        n_units
      )
    ) &&
    all(
      vapply(
        seq_len(p_lags),
        function(lag) {
          !any(
            diag(
              lag_mask[, , lag]
            )
          )
        },
        logical(1L)
      )
    ) &&
    !any(
      diag(unit_mask)
    )
)


check(
  "PR-AUC helper",
  same(
    benchmark_pr_auc(
      y_true = c(
        1L,
        0L,
        1L,
        0L
      ),
      score = c(
        4,
        2,
        3,
        1
      )
    ),
    1
  )
)


check(
  "RMSE helper",
  same(
    benchmark_rmse(
      truth = c(
        0,
        1,
        2
      ),
      estimate = c(
        0,
        1,
        2
      )
    ),
    0
  )
)


# Fit settings -----

fit_seed <- 991L

bayesian_chains <- 2L
bayesian_burnin <- 50L
bayesian_draws <- 80L


model_args <- list(
  m3 = list(
    chains = bayesian_chains,
    burnin = bayesian_burnin,
    draws = bayesian_draws,
    thin = 1L,
    beta_algorithm = "chol"
  ),
  
  gigg = list(
    chains = bayesian_chains,
    burnin = bayesian_burnin,
    draws = bayesian_draws,
    thin = 1L,
    beta_algorithm = "chol"
  ),
  
  half_t = list(
    chains = bayesian_chains,
    burnin = bayesian_burnin,
    draws = bayesian_draws,
    thin = 1L,
    beta_algorithm = "chol"
  ),
  
  adelie_gen = list(
    n_threads = 1L
  ),
  
  bigvar_hlag = list(),
  
  mar = list(
    n_starts = 3L,
    max_iter = 100L
  ),
  
  nirvar = list(
    k = 2L
  )
)


# Dispatched fits -----

cat(
  "\nRunning benchmark dispatcher...\n\n"
)


fits <- fit_benchmarks(
  models = benchmark_models,
  Y_list = Y_list,
  p_lags = p_lags,
  seed = fit_seed,
  model_args = model_args,
  verbose = TRUE
)


check(
  "Returned model order",
  identical(
    names(fits),
    expected_models
  )
)


backend_checks <- vapply(
  expected_models,
  function(model) {
    
    identical(
      fits[[model]]$backend,
      expected_backends[[model]]
    )
  },
  logical(1L)
)


check(
  "Benchmark backends",
  all(
    backend_checks
  )
)


# Common contract -----

contract_checks <- vapply(
  expected_models,
  function(model) {
    
    tryCatch(
      {
        validate_benchmark_fit(
          fit = fits[[model]],
          model = model,
          Y_list = Y_list,
          p_lags = p_lags
        )
        
        TRUE
      },
      error = function(error) {
        
        cat(
          sprintf(
            "Contract failure for %s: %s\n",
            model,
            conditionMessage(error)
          )
        )
        
        FALSE
      }
    )
  },
  logical(1L)
)


check(
  "Common fit contract",
  all(
    contract_checks
  )
)


# Reconstruction -----

A_reconstruction_error <- setNames(
  numeric(
    length(expected_models)
  ),
  expected_models
)

lag_reconstruction_error <- setNames(
  numeric(
    length(expected_models)
  ),
  expected_models
)

max_reconstruction_error <- setNames(
  numeric(
    length(expected_models)
  ),
  expected_models
)

rms_reconstruction_error <- setNames(
  numeric(
    length(expected_models)
  ),
  expected_models
)


for (model in expected_models) {
  
  fit <- fits[[model]]
  
  reconstructed_A <- beta_to_A_list(
    beta = fit$beta_hat,
    p_lags = p_lags
  )
  
  reconstructed_s <- A_list_to_s_lag(
    A_list = reconstructed_A,
    n_units = n_units,
    m = m
  )
  
  reconstructed_max <- s_lag_to_unit(
    s_lag = reconstructed_s,
    method = "max"
  )
  
  reconstructed_rms <- s_lag_to_unit(
    s_lag = reconstructed_s,
    method = "rms"
  )
  
  A_reconstruction_error[model] <- maximum_A_error(
    reconstructed_A,
    fit$A_hat_lag
  )
  
  lag_reconstruction_error[model] <- max(
    abs(
      reconstructed_s -
        fit$s_hat_lag
    )
  )
  
  max_reconstruction_error[model] <- max(
    abs(
      reconstructed_max -
        fit$s_hat_unit_max
    )
  )
  
  rms_reconstruction_error[model] <- max(
    abs(
      reconstructed_rms -
        fit$s_hat_unit_rms
    )
  )
}


reconstruction_tolerance <- 1e-10


check(
  "Coefficient reconstruction",
  all(
    A_reconstruction_error <
      reconstruction_tolerance
  )
)


check(
  "Lag network reconstruction",
  all(
    lag_reconstruction_error <
      reconstruction_tolerance
  )
)


check(
  "Unit max reconstruction",
  all(
    max_reconstruction_error <
      reconstruction_tolerance
  )
)


check(
  "Unit RMS reconstruction",
  all(
    rms_reconstruction_error <
      reconstruction_tolerance
  )
)


# Direct fits -----

cat(
  "\nRunning direct fits...\n\n"
)


set.seed(
  fit_seed
)

direct_m3 <- fit_m3(
  Y_list = Y_list,
  p_lags = p_lags,
  chains = bayesian_chains,
  burnin = bayesian_burnin,
  draws = bayesian_draws,
  thin = 1L,
  seed = fit_seed,
  beta_algorithm = "chol"
)


set.seed(
  fit_seed
)

direct_gigg <- fit_gigg(
  Y_list = Y_list,
  p_lags = p_lags,
  chains = bayesian_chains,
  burnin = bayesian_burnin,
  draws = bayesian_draws,
  thin = 1L,
  seed = fit_seed,
  beta_algorithm = "chol"
)


set.seed(
  fit_seed
)

direct_half_t <- fit_half_t(
  Y_list = Y_list,
  p_lags = p_lags,
  global_grouping = "self_diagonal",
  chains = bayesian_chains,
  burnin = bayesian_burnin,
  draws = bayesian_draws,
  thin = 1L,
  seed = fit_seed,
  beta_algorithm = "chol"
)


set.seed(
  fit_seed
)

direct_adelie <- fit_adelie(
  Y_list = Y_list,
  p_lags = p_lags,
  n_threads = 1L,
  seed = fit_seed
)


set.seed(
  fit_seed
)

direct_bigvar <- fit_bigvar(
  Y_list = Y_list,
  p_lags = p_lags
)


set.seed(
  fit_seed
)

direct_mar <- fit_mar(
  Y_list = Y_list,
  p_lags = p_lags,
  n_starts = 3L,
  seed = fit_seed,
  max_iter = 100L
)


set.seed(
  fit_seed
)

direct_nirvar <- fit_nirvar(
  Y_list = Y_list,
  p_lags = p_lags,
  k = 2L,
  seed = fit_seed
)


direct_fits <- list(
  m3 = direct_m3,
  gigg = direct_gigg,
  half_t = direct_half_t,
  adelie_gen = direct_adelie,
  bigvar_hlag = direct_bigvar,
  mar = direct_mar,
  nirvar = direct_nirvar
)


# Direct-dispatch equivalence -----

direct_beta_error <- setNames(
  numeric(
    length(expected_models)
  ),
  expected_models
)

direct_A_error <- setNames(
  numeric(
    length(expected_models)
  ),
  expected_models
)

direct_network_error <- setNames(
  numeric(
    length(expected_models)
  ),
  expected_models
)

direct_score_error <- setNames(
  numeric(
    length(expected_models)
  ),
  expected_models
)

direct_selection_same <- setNames(
  logical(
    length(expected_models)
  ),
  expected_models
)


for (model in expected_models) {
  
  dispatched <- fits[[model]]
  direct <- direct_fits[[model]]
  
  direct_beta_error[model] <- max(
    abs(
      dispatched$beta_hat -
        direct$beta_hat
    )
  )
  
  direct_A_error[model] <- maximum_A_error(
    dispatched$A_hat_lag,
    direct$A_hat_lag
  )
  
  direct_network_error[model] <- max(
    abs(
      dispatched$s_hat_lag -
        direct$s_hat_lag
    )
  )
  
  direct_score_error[model] <- max(
    abs(
      dispatched$score_lag -
        direct$score_lag
    )
  )
  
  direct_selection_same[model] <- identical(
    dispatched$selected_lag,
    direct$selected_lag
  )
}


direct_tolerance <- 1e-10


check(
  "Direct-dispatch coefficients",
  all(
    direct_beta_error <
      direct_tolerance
  )
)


check(
  "Direct-dispatch lag matrices",
  all(
    direct_A_error <
      direct_tolerance
  )
)


check(
  "Direct-dispatch networks",
  all(
    direct_network_error <
      direct_tolerance
  )
)


check(
  "Direct-dispatch scores",
  all(
    direct_score_error <
      direct_tolerance
  )
)


check(
  "Direct-dispatch selections",
  all(
    direct_selection_same
  )
)


# Bayesian settings -----

expected_retained <-
  bayesian_chains *
  bayesian_draws


half_t_checks <- c(
  grouping =
    identical(
      fits$half_t$global_grouping,
      "self_diagonal"
    ),
  
  tau_df =
    same(
      fits$half_t$prior$tau_df,
      10
    ),
  
  lambda_df =
    same(
      fits$half_t$prior$lambda_df,
      3
    ),
  
  tau_scale =
    same(
      fits$half_t$prior$tau_scale,
      1
    ),
  
  standardize =
    isTRUE(
      fits$half_t$preprocessing$standardize
    ),
  
  asis =
    isTRUE(
      fits$half_t$control$use_asis
    ),
  
  beta_algorithm =
    identical(
      fits$half_t$beta_algorithm_resolved,
      "chol"
    ),
  
  retained =
    fits$half_t$retained_total ==
    expected_retained
)


check(
  "Half-t benchmark settings",
  all(
    half_t_checks
  )
)


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
    same(
      fits$m3$prior$tau_df,
      10
    ),
  
  c_df =
    same(
      fits$m3$prior$c_df,
      5
    ),
  
  lambda_df =
    same(
      fits$m3$prior$lambda_df,
      3
    ),
  
  standardize =
    isTRUE(
      fits$m3$preprocessing$standardize
    ),
  
  c_asis =
    isTRUE(
      fits$m3$control$use_c_asis
    ),
  
  q_update =
    identical(
      fits$m3$q_update,
      "gibbs_end"
    ),
  
  beta_algorithm =
    identical(
      fits$m3$beta_algorithm_resolved,
      "chol"
    ),
  
  retained =
    fits$m3$retained_total ==
    expected_retained
)


check(
  "M3 benchmark settings",
  all(
    m3_checks
  )
)


stale_gigg_controls <- c(
  "use_global_tau",
  "numerical_epsilon",
  "invariant_tolerance",
  "gig_maximum_attempts"
)


gigg_checks <- c(
  tau_df =
    same(
      fits$gigg$prior$tau_df,
      5
    ),
  
  gamma_shape =
    same(
      fits$gigg$prior$gamma_shape,
      0.5
    ),
  
  gamma_rate =
    same(
      fits$gigg$prior$gamma_rate,
      1
    ),
  
  lambda_shape =
    same(
      fits$gigg$prior$lambda_shape,
      2.5
    ),
  
  lambda_scale =
    same(
      fits$gigg$prior$lambda_scale,
      1
    ),
  
  standardize =
    isTRUE(
      fits$gigg$preprocessing$standardize
    ),
  
  asis =
    isTRUE(
      fits$gigg$control$use_asis
    ),
  
  lean_control =
    !any(
      stale_gigg_controls %in%
        names(
          fits$gigg$control
        )
    ),
  
  beta_algorithm =
    identical(
      fits$gigg$beta_algorithm_resolved,
      "chol"
    ),
  
  retained =
    fits$gigg$retained_total ==
    expected_retained
)


check(
  "GIGG benchmark settings",
  all(
    gigg_checks
  )
)


# Seeds and runtimes -----

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
    ) == 1L &&
      is.finite(
        fit$runtime_seconds
      ) &&
      fit$runtime_seconds >= 0
  },
  logical(1L)
)


check(
  "Fit seeds",
  all(
    seed_checks
  )
)


check(
  "Fit runtimes",
  all(
    runtime_checks
  )
)


# Benchmark metrics -----

metric_rows <- lapply(
  expected_models,
  function(model) {
    
    metrics <- evaluate_benchmark_fit(
      fit = fits[[model]],
      truth = truth
    )
    
    data.frame(
      model = model,
      backend = fits[[model]]$backend,
      runtime_seconds = fits[[model]]$runtime_seconds,
      selected_blocks = sum(
        fits[[model]]$selected_lag
      ),
      metrics,
      stringsAsFactors = FALSE
    )
  }
)


comparison <- do.call(
  rbind,
  metric_rows
)

rownames(
  comparison
) <- NULL


expected_lag_positives <- sum(
  truth$G_lag[
    lag_mask
  ]
)

expected_lag_negatives <- sum(
  !truth$G_lag[
    lag_mask
  ]
)

expected_unit_positives <- sum(
  truth$G_unit[
    unit_mask
  ]
)

expected_unit_negatives <- sum(
  !truth$G_unit[
    unit_mask
  ]
)


check(
  "Metric truth counts",
  all(
    comparison$lag_positives ==
      expected_lag_positives
  ) &&
    all(
      comparison$lag_negatives ==
        expected_lag_negatives
    ) &&
    all(
      comparison$unit_positives ==
        expected_unit_positives
    ) &&
    all(
      comparison$unit_negatives ==
        expected_unit_negatives
    )
)


metric_numeric_fields <- c(
  "lag_pr_auc",
  "unit_max_pr_auc",
  "unit_rms_rmse",
  "unit_rms_rmse_active",
  "unit_rms_inactive_leakage",
  "mean_false",
  "maximum_false",
  "self_block_rmse"
)


check(
  "Benchmark metrics finite",
  all(
    vapply(
      metric_numeric_fields,
      function(field) {
        
        all(
          is.finite(
            comparison[[field]]
          )
        )
      },
      logical(1L)
    )
  )
)


check(
  "PR-AUC bounds",
  all(
    comparison$lag_pr_auc >= 0 &
      comparison$lag_pr_auc <= 1
  ) &&
    all(
      comparison$unit_max_pr_auc >= 0 &
        comparison$unit_max_pr_auc <= 1
    )
)


check(
  "Magnitude metrics non-negative",
  all(
    comparison$unit_rms_rmse >= 0
  ) &&
    all(
      comparison$unit_rms_rmse_active >= 0
    ) &&
    all(
      comparison$unit_rms_inactive_leakage >= 0
    ) &&
    all(
      comparison$mean_false >= 0
    ) &&
    all(
      comparison$maximum_false >= 0
    ) &&
    all(
      comparison$self_block_rmse >= 0
    )
)


# Auxiliary coefficient comparison -----

beta_rmse <- vapply(
  expected_models,
  function(model) {
    
    sqrt(
      mean(
        (
          fits[[model]]$beta_hat -
            truth_beta
        )^2
      )
    )
  },
  numeric(1L)
)


check(
  "Coefficient RMSE finite",
  finite(
    beta_rmse
  )
)


comparison$beta_rmse <- unname(
  beta_rmse[
    comparison$model
  ]
)


# Report -----

comparison_print <- comparison


round_fields <- c(
  "runtime_seconds",
  "beta_rmse",
  "lag_pr_auc",
  "unit_max_pr_auc",
  "unit_rms_rmse",
  "unit_rms_rmse_active",
  "unit_rms_inactive_leakage",
  "mean_false",
  "maximum_false",
  "self_block_rmse"
)


for (field in round_fields) {
  
  comparison_print[[field]] <- round(
    comparison_print[[field]],
    5L
  )
}


cat(
  "\nBENCHMARK VALIDATION\n\n"
)


cat(
  "Controlled DGP\n"
)

print(
  data.frame(
    n_units = n_units,
    m = m,
    N = N,
    p_lags = p_lags,
    T_obs = T_obs,
    spectral_radius = round(
      stability$radius,
      6L
    )
  ),
  row.names = FALSE
)


cat(
  "\nModel comparison\n"
)

print(
  comparison_print,
  row.names = FALSE
)


cat(
  "\nChecks\n"
)

print(
  data.frame(
    check = names(checks),
    passed = unname(checks)
  ),
  row.names = FALSE
)


cat(
  "\nOVERALL:",
  if (all(checks)) {
    "PASS\n"
  } else {
    "FAIL\n"
  }
)