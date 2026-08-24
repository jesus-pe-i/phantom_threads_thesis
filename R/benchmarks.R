# Common fitting interface for all benchmark estimators.
#
# Registers model-specific fit functions, loads their dependencies,
# dispatches fits through a common interface, records runtime and seeds,
# and validates the shared network-output contract.


# Model registry -----

benchmark_registry <- list(
  m3 = list(
    fit = "fit_m3",
    backend = "RcppArmadillo",
    uses_seed = TRUE,
    defaults = list()
  ),
  
  gigg = list(
    fit = "fit_gigg",
    backend = "RcppArmadillo",
    uses_seed = TRUE,
    defaults = list()
  ),
  
  half_t = list(
    fit = "fit_half_t",
    backend = "RcppArmadillo",
    uses_seed = TRUE,
    defaults = list(
      global_grouping = "self_diagonal"
    )
  ),
  
  adelie_gen = list(
    fit = "fit_adelie",
    backend = "adelie",
    uses_seed = TRUE,
    defaults = list()
  ),
  
  bigvar_hlag = list(
    fit = "fit_bigvar",
    backend = "BigVAR",
    uses_seed = FALSE,
    defaults = list()
  ),
  
  mar = list(
    fit = "fit_mar",
    backend = "mar",
    uses_seed = TRUE,
    defaults = list()
  ),
  
  nirvar = list(
    fit = "fit_nirvar",
    backend = "nirvar",
    uses_seed = TRUE,
    defaults = list()
  )
)

benchmark_models <- names(
  benchmark_registry
)


# Source loading -----

benchmark_source_files <- c(
  "R/var_data.R",
  "R/var_coefficients.R",
  "R/network_blocks.R",
  
  "R/half_t_structure.R",
  "R/half_t_sampler.R",
  "R/half_t_fit.R",
  
  "R/m3_structure.R",
  "R/m3_sampler.R",
  "R/m3_fit.R",
  
  "R/gigg_structure.R",
  "R/gigg_sampler.R",
  "R/gigg_fit.R",
  
  "R/fit_adelie.R",
  "R/fit_bigvar.R",
  "R/fit_mar.R",
  "R/fit_nirvar.R"
)


source_benchmark_models <- function(
    root = ".") {
  
  target_environment <- parent.frame()
  
  for (file in benchmark_source_files) {
    
    path <- file.path(
      root,
      file
    )
    
    if (!file.exists(path)) {
      stop(
        "Missing benchmark source file: ",
        path
      )
    }
    
    source(
      path,
      local = target_environment
    )
  }
  
  invisible(
    benchmark_source_files
  )
}


# Common contract -----

validate_benchmark_fit <- function(
    fit,
    model,
    Y_list,
    p_lags) {
  
  required_fields <- c(
    "model",
    "backend",
    "beta_hat",
    "A_hat_lag",
    "s_hat_lag",
    "s_hat_unit_max",
    "s_hat_unit_rms",
    "score_lag",
    "selected_lag",
    "dimensions",
    "preprocessing",
    "fit_seed",
    "runtime_seconds"
  )
  
  missing_fields <- setdiff(
    required_fields,
    names(fit)
  )
  
  if (length(missing_fields) > 0L) {
    stop(
      "Benchmark fit is missing: ",
      paste(
        missing_fields,
        collapse = ", "
      )
    )
  }
  
  n_units <- length(
    Y_list
  )
  
  m <- ncol(
    Y_list[[1L]]
  )
  
  N <- n_units * m
  k <- N * p_lags
  
  expected_network_dim <- as.integer(
    c(
      n_units,
      n_units,
      p_lags
    )
  )
  
  expected_unit_dim <- as.integer(
    c(
      n_units,
      n_units
    )
  )
  
  if (!identical(
    fit$model,
    model
  )) {
    stop(
      "Returned model identifier does not match requested model"
    )
  }
  
  if (!identical(
    dim(fit$beta_hat),
    as.integer(
      c(
        k,
        N
      )
    )
  )) {
    stop(
      "beta_hat does not follow the common k x N convention"
    )
  }
  
  if (
    length(fit$A_hat_lag) != p_lags ||
    any(
      vapply(
        fit$A_hat_lag,
        function(A) {
          !identical(
            dim(A),
            as.integer(
              c(
                N,
                N
              )
            )
          )
        },
        logical(1L)
      )
    )
  ) {
    stop(
      "A_hat_lag does not contain p_lags N x N matrices"
    )
  }
  
  numeric_networks <- c(
    "s_hat_lag",
    "score_lag"
  )
  
  for (field in numeric_networks) {
    
    if (
      !is.numeric(fit[[field]]) ||
      !identical(
        dim(fit[[field]]),
        expected_network_dim
      ) ||
      any(
        !is.finite(
          fit[[field]]
        )
      )
    ) {
      stop(
        field,
        " does not satisfy the common network contract"
      )
    }
  }
  
  if (
    !is.logical(fit$selected_lag) ||
    !identical(
      dim(fit$selected_lag),
      expected_network_dim
    ) ||
    anyNA(
      fit$selected_lag
    )
  ) {
    stop(
      "selected_lag does not satisfy the common network contract"
    )
  }
  
  if (
    !identical(
      dim(fit$s_hat_unit_max),
      expected_unit_dim
    ) ||
    !identical(
      dim(fit$s_hat_unit_rms),
      expected_unit_dim
    )
  ) {
    stop(
      "Unit-level network strengths have incorrect dimensions"
    )
  }
  
  expected_dimensions <- c(
    n_units = n_units,
    m = m,
    N = N,
    p_lags = p_lags,
    k = k
  )
  
  observed_dimensions <- unlist(
    fit$dimensions[
      names(
        expected_dimensions
      )
    ],
    use.names = TRUE
  )
  
  if (
    length(observed_dimensions) !=
    length(expected_dimensions) ||
    any(
      observed_dimensions !=
      expected_dimensions
    )
  ) {
    stop(
      "Benchmark dimensions record is inconsistent"
    )
  }
  
  invisible(
    TRUE
  )
}


# Single-model dispatch -----

fit_benchmark <- function(
    model,
    Y_list,
    p_lags,
    seed = 991L,
    args = list()) {
  
  if (!model %in% benchmark_models) {
    stop(
      "Unknown benchmark model: ",
      model
    )
  }
  
  if (!is.list(args)) {
    stop(
      "args must be a list"
    )
  }
  
  reserved_arguments <- intersect(
    names(args),
    c(
      "Y_list",
      "p_lags",
      "seed"
    )
  )
  
  if (length(reserved_arguments) > 0L) {
    stop(
      "Model-specific args cannot override: ",
      paste(
        reserved_arguments,
        collapse = ", "
      )
    )
  }
  
  seed <- as.integer(
    seed
  )
  
  specification <-
    benchmark_registry[[model]]
  
  fit_function_name <-
    specification$fit
  
  if (!exists(
    fit_function_name,
    mode = "function",
    inherits = TRUE
  )) {
    stop(
      "Fit function '",
      fit_function_name,
      "' is not loaded. ",
      "Run source_benchmark_models() first."
    )
  }
  
  fit_function <- get(
    fit_function_name,
    mode = "function",
    inherits = TRUE
  )
  
  model_arguments <- utils::modifyList(
    specification$defaults,
    args
  )
  
  call_arguments <- c(
    list(
      Y_list = Y_list,
      p_lags = p_lags
    ),
    model_arguments
  )
  
  if (specification$uses_seed) {
    call_arguments$seed <- seed
  }
  
  set.seed(
    seed
  )
  
  timing <- system.time(
    fit <- do.call(
      fit_function,
      call_arguments
    )
  )
  
  if (
    is.null(fit$model) ||
    !identical(
      fit$model,
      model
    )
  ) {
    stop(
      "Model '",
      model,
      "' returned an inconsistent model identifier"
    )
  }
  
  if (is.null(fit$backend)) {
    fit$backend <-
      specification$backend
  }
  
  fit$fit_seed <-
    seed
  
  fit$runtime_seconds <-
    as.numeric(
      timing[["elapsed"]]
    )
  
  validate_benchmark_fit(
    fit = fit,
    model = model,
    Y_list = Y_list,
    p_lags = p_lags
  )
  
  fit
}


# Multi-model dispatch -----

fit_benchmarks <- function(
    models = benchmark_models,
    Y_list,
    p_lags,
    seed = 991L,
    model_args = list(),
    verbose = TRUE) {
  
  if (
    length(models) == 0L ||
    any(
      !models %in%
      benchmark_models
    )
  ) {
    stop(
      "models must contain valid benchmark model names"
    )
  }
  
  if (anyDuplicated(models)) {
    stop(
      "models must not contain duplicates"
    )
  }
  
  if (!is.list(model_args)) {
    stop(
      "model_args must be a named list"
    )
  }
  
  unknown_args <- setdiff(
    names(model_args),
    benchmark_models
  )
  
  if (length(unknown_args) > 0L) {
    stop(
      "Unknown model_args entries: ",
      paste(
        unknown_args,
        collapse = ", "
      )
    )
  }
  
  fits <- setNames(
    vector(
      "list",
      length(models)
    ),
    models
  )
  
  for (index in seq_along(models)) {
    
    model <- models[index]
    
    args <- model_args[[model]]
    
    if (is.null(args)) {
      args <- list()
    }
    
    if (verbose) {
      cat(
        sprintf(
          "[%d/%d] Fitting %s...\n",
          index,
          length(models),
          model
        )
      )
      
      flush.console()
    }
    
    fits[[model]] <- fit_benchmark(
      model = model,
      Y_list = Y_list,
      p_lags = p_lags,
      seed = seed,
      args = args
    )
    
    if (verbose) {
      cat(
        sprintf(
          "      done in %.3f seconds\n",
          fits[[model]]$runtime_seconds
        )
      )
      
      flush.console()
    }
  }
  
  fits
}