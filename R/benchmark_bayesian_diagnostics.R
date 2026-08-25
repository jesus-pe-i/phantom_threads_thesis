# Builds shared Bayesian monitoring plans and reduces retained MCMC traces
# to compact convergence and mixing diagnostics for benchmark campaigns.


# Diagnostic seed -----

make_bayesian_diagnostic_seed <- function(
    record_id,
    seed) {
  
  if (
    length(record_id) != 1L ||
    is.na(record_id) ||
    !nzchar(record_id)
  ) {
    stop(
      "record_id must be one non-empty string"
    )
  }
  
  if (
    length(seed) != 1L ||
    !is.finite(seed)
  ) {
    stop(
      "seed must be one finite value"
    )
  }
  
  modulus <- 2147483646
  
  hash <-
    abs(
      as.numeric(seed)
    ) %%
    modulus
  
  bytes <- utf8ToInt(
    enc2utf8(
      record_id
    )
  )
  
  for (byte in bytes) {
    
    hash <- (
      hash * 131 +
        byte
    ) %%
      modulus
  }
  
  as.integer(
    hash + 1
  )
}


# Monitor plan -----

make_bayesian_monitor_plan <- function(
    record,
    seed,
    max_blocks = 100L,
    extra_beta = 200L) {
  
  if (
    is.null(record$id) ||
    is.null(record$dimensions)
  ) {
    stop(
      "record must contain id and dimensions"
    )
  }
  
  required_dimensions <- c(
    "n_units",
    "m",
    "N",
    "p_lags"
  )
  
  if (
    !all(
      required_dimensions %in%
      names(
        record$dimensions
      )
    )
  ) {
    stop(
      "record dimensions are incomplete"
    )
  }
  
  n_units <- as.integer(
    record$dimensions$n_units
  )
  
  m <- as.integer(
    record$dimensions$m
  )
  
  N <- as.integer(
    record$dimensions$N
  )
  
  p_lags <- as.integer(
    record$dimensions$p_lags
  )
  
  max_blocks <- as.integer(
    max_blocks
  )
  
  extra_beta <- as.integer(
    extra_beta
  )
  
  if (
    n_units < 1L ||
    m < 1L ||
    N != n_units * m ||
    p_lags < 1L
  ) {
    stop(
      "Invalid benchmark dimensions"
    )
  }
  
  if (
    max_blocks < 0L ||
    extra_beta < 0L
  ) {
    stop(
      "Monitor counts must be non-negative"
    )
  }
  
  diagnostic_seed <- make_bayesian_diagnostic_seed(
    record_id = record$id,
    seed = seed
  )
  
  
  ## Preserve RNG state -----
  
  had_random_seed <- exists(
    ".Random.seed",
    envir = .GlobalEnv,
    inherits = FALSE
  )
  
  if (had_random_seed) {
    
    previous_random_seed <- get(
      ".Random.seed",
      envir = .GlobalEnv,
      inherits = FALSE
    )
  }
  
  on.exit(
    {
      
      if (had_random_seed) {
        
        assign(
          ".Random.seed",
          previous_random_seed,
          envir = .GlobalEnv
        )
        
      } else if (
        exists(
          ".Random.seed",
          envir = .GlobalEnv,
          inherits = FALSE
        )
      ) {
        
        rm(
          list = ".Random.seed",
          envir = .GlobalEnv
        )
      }
    },
    add = TRUE
  )
  
  set.seed(
    diagnostic_seed
  )
  
  
  ## Blocks -----
  
  n_blocks <-
    n_units^2L *
    p_lags
  
  n_monitor_blocks <- min(
    max_blocks,
    n_blocks
  )
  
  if (n_monitor_blocks > 0L) {
    
    blocks <- sort(
      sample.int(
        n = n_blocks,
        size = n_monitor_blocks,
        replace = FALSE
      )
    )
    
  } else {
    
    blocks <- integer(0L)
  }
  
  
  ## One beta per block -----
  
  block_beta <- matrix(
    integer(0L),
    nrow = n_monitor_blocks,
    ncol = 2L
  )
  
  if (n_monitor_blocks > 0L) {
    
    for (index in seq_len(
      n_monitor_blocks
    )) {
      
      block <- blocks[index]
      
      zero <- block - 1L
      
      lag <-
        zero %/%
        n_units^2L +
        1L
      
      within_lag <-
        zero %%
        n_units^2L
      
      receiver_unit <-
        within_lag %/%
        n_units +
        1L
      
      sender_unit <-
        within_lag %%
        n_units +
        1L
      
      sender_feature <- sample.int(
        m,
        1L
      )
      
      receiver_feature <- sample.int(
        m,
        1L
      )
      
      sender_series <- (
        sender_unit - 1L
      ) * m +
        sender_feature
      
      receiver_series <- (
        receiver_unit - 1L
      ) * m +
        receiver_feature
      
      row <- (
        lag - 1L
      ) * N +
        sender_series
      
      equation <-
        receiver_series
      
      block_beta[
        index,
      ] <- c(
        row,
        equation
      )
    }
  }
  
  
  ## Extra beta coordinates -----
  
  k <-
    N *
    p_lags
  
  n_coefficients <-
    k *
    N
  
  block_beta_ids <- if (
    nrow(block_beta) > 0L
  ) {
    
    block_beta[, 1L] +
      (
        block_beta[, 2L] -
          1L
      ) *
      k
    
  } else {
    
    integer(0L)
  }
  
  available_beta <- setdiff(
    seq_len(
      n_coefficients
    ),
    block_beta_ids
  )
  
  n_extra_beta <- min(
    extra_beta,
    length(
      available_beta
    )
  )
  
  if (n_extra_beta > 0L) {
    
    extra_ids <- sample(
      available_beta,
      size = n_extra_beta,
      replace = FALSE
    )
    
    extra_beta_coordinates <- cbind(
      row = (
        (
          extra_ids -
            1L
        ) %%
          k
      ) +
        1L,
      
      equation = (
        (
          extra_ids -
            1L
        ) %/%
          k
      ) +
        1L
    )
    
  } else {
    
    extra_beta_coordinates <- matrix(
      integer(0L),
      nrow = 0L,
      ncol = 2L
    )
  }
  
  
  ## Combined beta monitor -----
  
  beta <- rbind(
    block_beta,
    extra_beta_coordinates
  )
  
  storage.mode(beta) <-
    "integer"
  
  colnames(beta) <- c(
    "row",
    "equation"
  )
  
  
  # Output -----
  
  list(
    seed =
      diagnostic_seed,
    
    blocks =
      as.integer(
        blocks
      ),
    
    beta =
      beta
  )
}


# Model mapping -----

make_bayesian_monitor_args <- function(
    model,
    plan) {
  
  model <- match.arg(
    model,
    c(
      "m3",
      "half_t",
      "gigg"
    )
  )
  
  if (
    is.null(plan$beta) ||
    is.null(plan$blocks)
  ) {
    stop(
      "Invalid Bayesian monitor plan"
    )
  }
  
  beta <- as.matrix(
    plan$beta
  )
  
  if (
    ncol(beta) != 2L
  ) {
    stop(
      "plan$beta must have two columns"
    )
  }
  
  storage.mode(beta) <-
    "integer"
  
  colnames(beta) <- c(
    "row",
    "equation"
  )
  
  blocks <- as.integer(
    plan$blocks
  )
  
  switch(
    model,
    
    m3 = list(
      monitor_beta =
        beta,
      
      monitor_lambda =
        blocks
    ),
    
    half_t = list(
      monitor_beta =
        beta,
      
      monitor_lambda =
        beta
    ),
    
    gigg = list(
      monitor_beta =
        beta,
      
      monitor_gamma =
        blocks,
      
      monitor_lambda =
        beta
    )
  )
}


# Parameter names -----

make_bayesian_diagnostic_parameter_names <- function(
    model,
    family,
    draw_name,
    draws) {
  
  parameter <- colnames(
    draws
  )
  
  if (
    is.null(parameter) ||
    length(parameter) !=
    ncol(draws) ||
    any(
      !nzchar(parameter)
    )
  ) {
    
    parameter <- paste0(
      draw_name,
      "_",
      seq_len(
        ncol(draws)
      )
    )
  }
  
  if (family == "beta") {
    
    parameter <- sub(
      "^beta_",
      "",
      parameter
    )
  }
  
  if (family == "lambda") {
    
    parameter <- sub(
      "^lambda2_",
      "",
      parameter
    )
  }
  
  if (family == "gamma") {
    
    parameter <- sub(
      "^gamma2_",
      "",
      parameter
    )
  }
  
  if (
    model == "m3" &&
    family == "lambda"
  ) {
    
    parameter <- sub(
      "_s([0-9]+)_r([0-9]+)$",
      "_su\\1_ru\\2",
      parameter
    )
  }
  
  parameter
}


# Diagnostic extraction -----

extract_bayesian_diagnostic_draws <- function(
    fit) {
  
  if (
    is.null(fit$model) ||
    is.null(fit$chain_results)
  ) {
    stop(
      "Bayesian fit must retain chain_results"
    )
  }
  
  model <- fit$model
  
  if (
    !model %in%
    c(
      "m3",
      "half_t",
      "gigg"
    )
  ) {
    stop(
      "Unsupported Bayesian model: ",
      model
    )
  }
  
  chain_results <- fit$chain_results
  
  if (length(chain_results) < 2L) {
    stop(
      "At least two chains are required for diagnostics"
    )
  }
  
  draw_spec <- switch(
    model,
    
    m3 = c(
      beta = "beta",
      lambda = "lambda2",
      tau = "tau2",
      c = "c2",
      q = "q"
    ),
    
    half_t = c(
      beta = "beta",
      lambda = "lambda2",
      tau = "tau2"
    ),
    
    gigg = c(
      beta = "beta",
      lambda = "lambda2",
      gamma = "gamma2",
      tau = "tau2"
    )
  )
  
  output <- list()
  
  for (family in names(
    draw_spec
  )) {
    
    draw_name <- unname(
      draw_spec[family]
    )
    
    chain_draws <- lapply(
      chain_results,
      function(chain) {
        
        if (
          is.null(
            chain$draws[[draw_name]]
          )
        ) {
          stop(
            "Missing ",
            draw_name,
            " draws for ",
            model
          )
        }
        
        as.matrix(
          chain$draws[[draw_name]]
        )
      }
    )
    
    if (
      model == "m3" &&
      family == "c"
    ) {
      
      chain_draws <- lapply(
        chain_draws,
        function(draws) {
          
          parameter <- colnames(
            draws
          )
          
          if (is.null(parameter)) {
            
            parameter <- paste0(
              draw_name,
              "_",
              seq_len(
                ncol(draws)
              )
            )
          }
          
          draws <- draws[
            ,
            -1L,
            drop = FALSE
          ]
          
          colnames(draws) <-
            parameter[-1L]
          
          draws
        }
      )
    }
    
    draws_per_chain <- vapply(
      chain_draws,
      nrow,
      integer(1L)
    )
    
    parameters_per_chain <- vapply(
      chain_draws,
      ncol,
      integer(1L)
    )
    
    if (
      length(
        unique(
          draws_per_chain
        )
      ) != 1L
    ) {
      stop(
        "Chains retain different numbers of draws"
      )
    }
    
    if (
      length(
        unique(
          parameters_per_chain
        )
      ) != 1L
    ) {
      stop(
        "Chains contain different monitor dimensions"
      )
    }
    
    if (parameters_per_chain[1L] == 0L) {
      next
    }
    
    parameter <- make_bayesian_diagnostic_parameter_names(
      model = model,
      family = family,
      draw_name = draw_name,
      draws = chain_draws[[1L]]
    )
    
    for (chain in seq_along(
      chain_draws
    )) {
      
      chain_parameter <-
        make_bayesian_diagnostic_parameter_names(
          model = model,
          family = family,
          draw_name = draw_name,
          draws = chain_draws[[chain]]
        )
      
      if (
        !identical(
          parameter,
          chain_parameter
        )
      ) {
        stop(
          "Diagnostic monitor order differs across chains"
        )
      }
    }
    
    output[[family]] <- list(
      parameter =
        parameter,
      
      draws =
        chain_draws
    )
  }
  
  output
}


# Diagnostic summary -----

summarise_bayesian_diagnostics <- function(
    fit) {
  
  if (
    !requireNamespace(
      "posterior",
      quietly = TRUE
    )
  ) {
    stop(
      "Package 'posterior' is required for Bayesian diagnostics"
    )
  }
  
  extracted <- extract_bayesian_diagnostic_draws(
    fit
  )
  
  rows <- list()
  row_index <- 0L
  
  for (family in names(
    extracted
  )) {
    
    family_draws <- extracted[[family]]
    
    n_chains <- length(
      family_draws$draws
    )
    
    draws_per_chain <- nrow(
      family_draws$draws[[1L]]
    )
    
    n_parameters <- length(
      family_draws$parameter
    )
    
    total_draws <-
      n_chains *
      draws_per_chain
    
    draw_array <- array(
      NA_real_,
      dim = c(
        draws_per_chain,
        n_chains,
        n_parameters
      ),
      dimnames = list(
        iteration = NULL,
        chain = NULL,
        variable =
          family_draws$parameter
      )
    )
    
    for (chain in seq_len(
      n_chains
    )) {
      
      draw_array[
        ,
        chain,
      ] <- family_draws$draws[[chain]]
    }
    
    posterior_draws <-
      posterior::as_draws_array(
        draw_array
      )
    
    diagnostic_summary <- as.data.frame(
      posterior::summarise_draws(
        posterior_draws,
        "rhat",
        "ess_bulk",
        "ess_tail"
      )
    )
    
    if (
      nrow(diagnostic_summary) !=
      n_parameters
    ) {
      stop(
        "Unexpected posterior diagnostic dimensions"
      )
    }
    
    row_index <- row_index + 1L
    
    rows[[row_index]] <- data.frame(
      family =
        family,
      
      parameter =
        family_draws$parameter,
      
      rhat =
        as.numeric(
          diagnostic_summary[["rhat"]]
        ),
      
      ess_bulk =
        as.numeric(
          diagnostic_summary[["ess_bulk"]]
        ),
      
      ess_tail =
        as.numeric(
          diagnostic_summary[["ess_tail"]]
        ),
      
      inefficiency_factor =
        total_draws /
        as.numeric(
          diagnostic_summary[["ess_bulk"]]
        ),
      
      n_chains =
        n_chains,
      
      draws_per_chain =
        draws_per_chain,
      
      total_draws =
        total_draws,
      
      stringsAsFactors =
        FALSE
    )
  }
  
  diagnostics <- do.call(
    rbind,
    rows
  )
  
  rownames(
    diagnostics
  ) <- NULL
  
  diagnostics
}