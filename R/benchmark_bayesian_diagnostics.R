# Builds shared Bayesian monitoring plans and reduces retained MCMC traces
# to compact convergence and mixing diagnostics for benchmark campaigns.
# All Bayesian models also report effective coefficient shrinkage multipliers.


# Diagnostic seed -----

make_bayesian_diagnostic_seed <- function(
    record_id,
    seed) {
  
  if (
    length(record_id) != 1L ||
    is.na(record_id) ||
    !nzchar(record_id)
  ) {
    stop("record_id must be one non-empty string")
  }
  
  if (
    length(seed) != 1L ||
    !is.finite(seed)
  ) {
    stop("seed must be one finite value")
  }
  
  modulus <- 2147483646
  hash <- abs(as.numeric(seed)) %% modulus
  bytes <- utf8ToInt(enc2utf8(record_id))
  
  for (byte in bytes) {
    hash <- (
      hash * 131 +
        byte
    ) %% modulus
  }
  
  as.integer(hash + 1)
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
    stop("record must contain id and dimensions")
  }
  
  dimensions <- record$dimensions
  
  required <- c(
    "n_units",
    "m",
    "N",
    "p_lags"
  )
  
  if (!all(required %in% names(dimensions))) {
    stop("record dimensions are incomplete")
  }
  
  n_units <- as.integer(dimensions$n_units)
  m <- as.integer(dimensions$m)
  N <- as.integer(dimensions$N)
  p_lags <- as.integer(dimensions$p_lags)
  max_blocks <- as.integer(max_blocks)
  extra_beta <- as.integer(extra_beta)
  
  if (
    n_units < 1L ||
    m < 1L ||
    N != n_units * m ||
    p_lags < 1L
  ) {
    stop("Invalid benchmark dimensions")
  }
  
  if (
    max_blocks < 0L ||
    extra_beta < 0L
  ) {
    stop("Monitor counts must be non-negative")
  }
  
  diagnostic_seed <- make_bayesian_diagnostic_seed(
    record_id = record$id,
    seed = seed
  )
  
  
  ## Preserve RNG state -----
  
  had_seed <- exists(
    ".Random.seed",
    envir = .GlobalEnv,
    inherits = FALSE
  )
  
  if (had_seed) {
    old_seed <- get(
      ".Random.seed",
      envir = .GlobalEnv,
      inherits = FALSE
    )
  }
  
  on.exit(
    {
      if (had_seed) {
        assign(
          ".Random.seed",
          old_seed,
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
  
  set.seed(diagnostic_seed)
  
  
  ## Blocks -----
  
  n_blocks <- n_units^2L * p_lags
  n_monitor_blocks <- min(
    max_blocks,
    n_blocks
  )
  
  blocks <- if (n_monitor_blocks > 0L) {
    sort(
      sample.int(
        n_blocks,
        n_monitor_blocks,
        replace = FALSE
      )
    )
  } else {
    integer(0L)
  }
  
  
  ## One beta per block -----
  
  block_beta <- matrix(
    integer(0L),
    nrow = n_monitor_blocks,
    ncol = 2L
  )
  
  if (n_monitor_blocks > 0L) {
    
    for (index in seq_len(n_monitor_blocks)) {
      
      zero <- blocks[index] - 1L
      lag <- zero %/% n_units^2L + 1L
      within_lag <- zero %% n_units^2L
      
      receiver_unit <-
        within_lag %/% n_units +
        1L
      
      sender_unit <-
        within_lag %% n_units +
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
      
      block_beta[index, ] <- c(
        (
          lag - 1L
        ) * N +
          sender_series,
        receiver_series
      )
    }
  }
  
  
  ## Extra beta coordinates -----
  
  k <- N * p_lags
  n_coefficients <- k * N
  
  block_beta_ids <- if (
    nrow(block_beta) > 0L
  ) {
    block_beta[, 1L] +
      (
        block_beta[, 2L] - 1L
      ) * k
  } else {
    integer(0L)
  }
  
  available_beta <- setdiff(
    seq_len(n_coefficients),
    block_beta_ids
  )
  
  n_extra <- min(
    extra_beta,
    length(available_beta)
  )
  
  extra <- if (n_extra > 0L) {
    
    ids <- sample(
      available_beta,
      n_extra,
      replace = FALSE
    )
    
    cbind(
      row = (
        ids - 1L
      ) %% k + 1L,
      
      equation = (
        ids - 1L
      ) %/% k + 1L
    )
    
  } else {
    
    matrix(
      integer(0L),
      nrow = 0L,
      ncol = 2L
    )
  }
  
  
  ## Combined beta monitor -----
  
  beta <- rbind(
    block_beta,
    extra
  )
  
  storage.mode(beta) <- "integer"
  
  colnames(beta) <- c(
    "row",
    "equation"
  )
  
  
  ## Blocks containing monitored betas -----
  
  if (nrow(beta) > 0L) {
    
    beta_lag <- (
      beta[, 1L] - 1L
    ) %/% N + 1L
    
    beta_sender_series <- (
      beta[, 1L] - 1L
    ) %% N + 1L
    
    beta_receiver_series <-
      beta[, 2L]
    
    beta_sender_unit <- (
      beta_sender_series - 1L
    ) %/% m + 1L
    
    beta_receiver_unit <- (
      beta_receiver_series - 1L
    ) %/% m + 1L
    
    beta_blocks <- sort(
      unique(
        as.integer(
          (
            beta_lag - 1L
          ) * n_units^2L +
            (
              beta_receiver_unit - 1L
            ) * n_units +
            beta_sender_unit
        )
      )
    )
    
  } else {
    
    beta_blocks <- integer(0L)
  }
  
  
  # Output -----
  
  list(
    seed = diagnostic_seed,
    blocks = as.integer(blocks),
    beta_blocks = as.integer(beta_blocks),
    beta = beta
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
    is.null(plan$blocks) ||
    is.null(plan$beta_blocks)
  ) {
    stop("Invalid Bayesian monitor plan")
  }
  
  beta <- as.matrix(plan$beta)
  blocks <- as.integer(plan$blocks)
  beta_blocks <- as.integer(
    plan$beta_blocks
  )
  
  if (ncol(beta) != 2L) {
    stop("plan$beta must have two columns")
  }
  
  storage.mode(beta) <- "integer"
  
  colnames(beta) <- c(
    "row",
    "equation"
  )
  
  switch(
    model,
    
    m3 = list(
      monitor_beta = beta,
      monitor_lambda = beta_blocks
    ),
    
    half_t = list(
      monitor_beta = beta,
      monitor_lambda = beta
    ),
    
    gigg = list(
      monitor_beta = beta,
      monitor_gamma = blocks,
      monitor_lambda = beta
    )
  )
}


# Parameter names -----

make_bayesian_diagnostic_parameter_names <- function(
    model,
    family,
    draw_name,
    draws) {
  
  parameter <- colnames(draws)
  
  if (
    is.null(parameter) ||
    length(parameter) != ncol(draws) ||
    any(!nzchar(parameter))
  ) {
    parameter <- paste0(
      draw_name,
      "_",
      seq_len(ncol(draws))
    )
  }
  
  prefix <- switch(
    family,
    beta = "^beta_",
    lambda = "^lambda2_",
    gamma = "^gamma2_",
    omega = "^(effective_variance|omega)_",
    NULL
  )
  
  if (!is.null(prefix)) {
    parameter <- sub(
      prefix,
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


# M3 effective scales -----

make_m3_omega_draws <- function(
    fit) {
  
  if (!identical(fit$model, "m3")) {
    stop("M3 omega draws require an M3 fit")
  }
  
  if (
    is.null(fit$chain_results) ||
    is.null(fit$dimensions) ||
    is.null(fit$c_group) ||
    is.null(fit$q_group)
  ) {
    stop("M3 fit is incomplete for omega diagnostics")
  }
  
  if (
    fit$c_group != "self_foreign" ||
    fit$q_group != "self_foreign"
  ) {
    stop(
      "M3 omega diagnostics require self_foreign c and q groups"
    )
  }
  
  chains <- fit$chain_results
  
  n_units <- as.integer(
    fit$dimensions$n_units
  )
  
  m <- as.integer(
    fit$dimensions$m
  )
  
  N <- as.integer(
    fit$dimensions$N
  )
  
  p_lags <- as.integer(
    fit$dimensions$p_lags
  )
  
  k <- N * p_lags
  
  
  ## Monitors -----
  
  beta_monitors <- lapply(
    chains,
    function(chain) {
      
      monitor <- as.matrix(
        chain$monitor$beta
      )
      
      if (
        is.null(dim(monitor)) ||
        ncol(monitor) != 2L
      ) {
        stop("Invalid M3 beta monitor")
      }
      
      storage.mode(monitor) <- "integer"
      monitor
    }
  )
  
  lambda_monitors <- lapply(
    chains,
    function(chain) {
      as.integer(
        chain$monitor$lambda2
      )
    }
  )
  
  beta_monitor <-
    beta_monitors[[1L]]
  
  lambda_monitor <-
    lambda_monitors[[1L]]
  
  if (
    !all(
      vapply(
        beta_monitors,
        identical,
        logical(1L),
        beta_monitor
      )
    )
  ) {
    stop(
      "M3 beta monitor order differs across chains"
    )
  }
  
  if (
    !all(
      vapply(
        lambda_monitors,
        identical,
        logical(1L),
        lambda_monitor
      )
    )
  ) {
    stop(
      "M3 lambda monitor order differs across chains"
    )
  }
  
  if (
    any(
      beta_monitor[, 1L] < 1L |
      beta_monitor[, 1L] > k |
      beta_monitor[, 2L] < 1L |
      beta_monitor[, 2L] > N
    )
  ) {
    stop(
      "M3 beta monitor contains invalid coordinates"
    )
  }
  
  
  ## Coefficient map -----
  
  beta_reference <- as.matrix(
    chains[[1L]]$draws$beta
  )
  
  parameter <-
    make_bayesian_diagnostic_parameter_names(
      model = "m3",
      family = "beta",
      draw_name = "beta",
      draws = beta_reference
    )
  
  if (
    ncol(beta_reference) !=
    nrow(beta_monitor) ||
    length(parameter) !=
    nrow(beta_monitor)
  ) {
    stop(
      "M3 beta monitor and retained draws do not align"
    )
  }
  
  if (nrow(beta_monitor) == 0L) {
    
    return(
      list(
        parameter = character(0L),
        
        draws = lapply(
          chains,
          function(chain) {
            matrix(
              numeric(0L),
              nrow = nrow(
                as.matrix(
                  chain$draws$beta
                )
              ),
              ncol = 0L
            )
          }
        )
      )
    )
  }
  
  row <- beta_monitor[, 1L]
  equation <- beta_monitor[, 2L]
  
  lag <- (
    row - 1L
  ) %/% N + 1L
  
  sender_series <- (
    row - 1L
  ) %% N + 1L
  
  receiver_series <- equation
  
  sender_unit <- (
    sender_series - 1L
  ) %/% m + 1L
  
  receiver_unit <- (
    receiver_series - 1L
  ) %/% m + 1L
  
  sender_variable <- (
    sender_series - 1L
  ) %% m + 1L
  
  receiver_variable <- (
    receiver_series - 1L
  ) %% m + 1L
  
  block <- as.integer(
    (
      lag - 1L
    ) * n_units^2L +
      (
        receiver_unit - 1L
      ) * n_units +
      sender_unit
  )
  
  lambda_position <- match(
    block,
    lambda_monitor
  )
  
  if (anyNA(lambda_position)) {
    stop(
      "M3 lambda monitor does not cover every monitored beta block"
    )
  }
  
  self_block <-
    sender_unit ==
    receiver_unit
  
  same_variable <-
    sender_variable ==
    receiver_variable
  
  c_group <- ifelse(
    self_block,
    1L,
    2L
  )
  
  q_group <- c_group
  
  q_weight <- ifelse(
    same_variable,
    (
      m - 1
    ) / m,
    -1 / m
  )
  
  
  ## Effective draws -----
  
  omega_draws <- lapply(
    chains,
    function(chain) {
      
      beta_draws <- as.matrix(
        chain$draws$beta
      )
      
      lambda_draws <- as.matrix(
        chain$draws$lambda2
      )
      
      tau_draws <- as.matrix(
        chain$draws$tau2
      )
      
      c_draws <- as.matrix(
        chain$draws$c2
      )
      
      q_draws <- as.matrix(
        chain$draws$q
      )
      
      retained <- nrow(
        beta_draws
      )
      
      if (
        ncol(beta_draws) !=
        nrow(beta_monitor) ||
        ncol(lambda_draws) !=
        length(lambda_monitor)
      ) {
        stop(
          "Unexpected M3 monitor dimensions"
        )
      }
      
      if (
        any(
          c(
            nrow(lambda_draws),
            nrow(tau_draws),
            nrow(c_draws),
            nrow(q_draws)
          ) != retained
        )
      ) {
        stop(
          "M3 hierarchy draws retain different iterations"
        )
      }
      
      if (
        ncol(tau_draws) < 1L ||
        ncol(c_draws) < max(c_group) ||
        ncol(q_draws) < max(q_group)
      ) {
        stop(
          "M3 hierarchy draws do not contain the required groups"
        )
      }
      
      chain_parameter <-
        make_bayesian_diagnostic_parameter_names(
          model = "m3",
          family = "beta",
          draw_name = "beta",
          draws = beta_draws
        )
      
      if (!identical(parameter, chain_parameter)) {
        stop(
          "M3 beta parameter order differs across chains"
        )
      }
      
      lambda_selected <- lambda_draws[
        ,
        lambda_position,
        drop = FALSE
      ]
      
      c_selected <- c_draws[
        ,
        c_group,
        drop = FALSE
      ]
      
      q_selected <- q_draws[
        ,
        q_group,
        drop = FALSE
      ]
      
      anatomy <- exp(
        sweep(
          q_selected,
          2L,
          q_weight,
          FUN = "*"
        )
      )
      
      omega <-
        lambda_selected *
        c_selected *
        anatomy
      
      omega <- sweep(
        omega,
        1L,
        tau_draws[, 1L],
        FUN = "*"
      )
      
      colnames(omega) <- paste0(
        "omega_",
        parameter
      )
      
      omega
    }
  )
  
  list(
    parameter = parameter,
    draws = omega_draws
  )
}


# Half-t effective scales -----

make_half_t_omega_draws <- function(
    fit) {
  
  if (!identical(fit$model, "half_t")) {
    stop(
      "Half-t omega draws require a Half-t fit"
    )
  }
  
  if (
    is.null(fit$chain_results) ||
    is.null(fit$dimensions) ||
    is.null(fit$global_grouping)
  ) {
    stop(
      "Half-t fit is incomplete for omega diagnostics"
    )
  }
  
  chains <- fit$chain_results
  
  N <- as.integer(
    fit$dimensions$N
  )
  
  p_lags <- as.integer(
    fit$dimensions$p_lags
  )
  
  k <- N * p_lags
  
  if (
    !fit$global_grouping %in%
    c(
      "all",
      "self_diagonal"
    )
  ) {
    stop(
      "Unknown Half-t global grouping"
    )
  }
  
  
  ## Monitor -----
  
  monitors <- lapply(
    chains,
    function(chain) {
      
      monitor <- as.matrix(
        chain$monitor$lambda2
      )
      
      if (
        is.null(dim(monitor)) ||
        ncol(monitor) != 2L
      ) {
        stop(
          "Invalid Half-t lambda monitor"
        )
      }
      
      storage.mode(monitor) <- "integer"
      monitor
    }
  )
  
  monitor <- monitors[[1L]]
  
  if (
    !all(
      vapply(
        monitors,
        identical,
        logical(1L),
        monitor
      )
    )
  ) {
    stop(
      "Half-t lambda monitor order differs across chains"
    )
  }
  
  if (
    any(
      monitor[, 1L] < 1L |
      monitor[, 1L] > k |
      monitor[, 2L] < 1L |
      monitor[, 2L] > N
    )
  ) {
    stop(
      "Half-t lambda monitor contains invalid coordinates"
    )
  }
  
  
  ## Coefficient map -----
  
  lambda_reference <- as.matrix(
    chains[[1L]]$draws$lambda2
  )
  
  parameter <-
    make_bayesian_diagnostic_parameter_names(
      model = "half_t",
      family = "lambda",
      draw_name = "lambda2",
      draws = lambda_reference
    )
  
  if (
    ncol(lambda_reference) !=
    nrow(monitor) ||
    length(parameter) !=
    nrow(monitor)
  ) {
    stop(
      "Half-t lambda monitor and retained draws do not align"
    )
  }
  
  if (nrow(monitor) == 0L) {
    
    return(
      list(
        parameter = character(0L),
        
        draws = lapply(
          chains,
          function(chain) {
            as.matrix(
              chain$draws$lambda2
            )
          }
        )
      )
    )
  }
  
  row <- monitor[, 1L]
  equation <- monitor[, 2L]
  
  lag <- (
    row - 1L
  ) %/% N + 1L
  
  sender <- (
    row - 1L
  ) %% N + 1L
  
  phi2 <- 1 / lag^2
  
  tau_group <- if (
    fit$global_grouping == "all"
  ) {
    rep(
      1L,
      nrow(monitor)
    )
  } else {
    ifelse(
      sender == equation,
      1L,
      2L
    )
  }
  
  
  ## Effective draws -----
  
  omega_draws <- lapply(
    chains,
    function(chain) {
      
      lambda_draws <- as.matrix(
        chain$draws$lambda2
      )
      
      tau_draws <- as.matrix(
        chain$draws$tau2
      )
      
      if (
        ncol(lambda_draws) !=
        nrow(monitor) ||
        nrow(lambda_draws) !=
        nrow(tau_draws)
      ) {
        stop(
          "Unexpected Half-t hierarchy draw dimensions"
        )
      }
      
      if (
        ncol(tau_draws) <
        max(tau_group)
      ) {
        stop(
          "Half-t tau draws do not contain the required groups"
        )
      }
      
      chain_parameter <-
        make_bayesian_diagnostic_parameter_names(
          model = "half_t",
          family = "lambda",
          draw_name = "lambda2",
          draws = lambda_draws
        )
      
      if (!identical(parameter, chain_parameter)) {
        stop(
          "Half-t lambda parameter order differs across chains"
        )
      }
      
      omega <-
        lambda_draws *
        tau_draws[
          ,
          tau_group,
          drop = FALSE
        ]
      
      omega <- sweep(
        omega,
        2L,
        phi2,
        FUN = "*"
      )
      
      colnames(omega) <- paste0(
        "omega_",
        parameter
      )
      
      omega
    }
  )
  
  list(
    parameter = parameter,
    draws = omega_draws
  )
}


# Stored diagnostics -----

extract_bayesian_stored_draws <- function(
    fit,
    family,
    draw_name) {
  
  model <- fit$model
  chains <- fit$chain_results
  
  chain_draws <- lapply(
    chains,
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
            seq_len(ncol(draws))
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
  
  draw_counts <- vapply(
    chain_draws,
    nrow,
    integer(1L)
  )
  
  parameter_counts <- vapply(
    chain_draws,
    ncol,
    integer(1L)
  )
  
  if (
    length(
      unique(draw_counts)
    ) != 1L
  ) {
    stop(
      "Chains retain different numbers of draws"
    )
  }
  
  if (
    length(
      unique(parameter_counts)
    ) != 1L
  ) {
    stop(
      "Chains contain different monitor dimensions"
    )
  }
  
  if (parameter_counts[1L] == 0L) {
    return(NULL)
  }
  
  parameter <-
    make_bayesian_diagnostic_parameter_names(
      model = model,
      family = family,
      draw_name = draw_name,
      draws = chain_draws[[1L]]
    )
  
  for (chain in seq_along(chain_draws)) {
    
    chain_parameter <-
      make_bayesian_diagnostic_parameter_names(
        model = model,
        family = family,
        draw_name = draw_name,
        draws = chain_draws[[chain]]
      )
    
    if (!identical(parameter, chain_parameter)) {
      stop(
        "Diagnostic monitor order differs across chains"
      )
    }
  }
  
  list(
    parameter = parameter,
    draws = chain_draws
  )
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
  
  if (
    length(fit$chain_results) <
    2L
  ) {
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
      omega = "effective_variance",
      gamma = "gamma2",
      tau = "tau2"
    )
  )
  
  output <- list()
  
  for (family in names(draw_spec)) {
    
    stored <- extract_bayesian_stored_draws(
      fit = fit,
      family = family,
      draw_name = unname(
        draw_spec[family]
      )
    )
    
    if (!is.null(stored)) {
      output[[family]] <- stored
    }
  }
  
  
  ## Derived omega -----
  
  if (model == "m3") {
    
    omega <- make_m3_omega_draws(
      fit
    )
    
    if (
      length(omega$parameter) >
      0L
    ) {
      output$omega <- omega
    }
  }
  
  if (model == "half_t") {
    
    omega <- make_half_t_omega_draws(
      fit
    )
    
    if (
      length(omega$parameter) >
      0L
    ) {
      output$omega <- omega
    }
  }
  
  
  ## Pairing checks -----
  
  if (
    "omega" %in%
    names(output)
  ) {
    
    if (
      !identical(
        output$beta$parameter,
        output$omega$parameter
      )
    ) {
      stop(
        "Beta and omega diagnostic coordinates do not align"
      )
    }
  }
  
  if (
    model %in%
    c(
      "half_t",
      "gigg"
    ) &&
    "omega" %in%
    names(output)
  ) {
    
    if (
      !identical(
        output$lambda$parameter,
        output$omega$parameter
      )
    ) {
      stop(
        "Lambda and omega diagnostic coordinates do not align"
      )
    }
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
  
  rows <- lapply(
    names(extracted),
    function(family) {
      
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
      
      for (chain in seq_len(n_chains)) {
        draw_array[
          ,
          chain,
        ] <- family_draws$draws[[chain]]
      }
      
      diagnostic_summary <- as.data.frame(
        posterior::summarise_draws(
          posterior::as_draws_array(
            draw_array
          ),
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
      
      data.frame(
        family = family,
        parameter =
          family_draws$parameter,
        
        rhat = as.numeric(
          diagnostic_summary[["rhat"]]
        ),
        
        ess_bulk = as.numeric(
          diagnostic_summary[["ess_bulk"]]
        ),
        
        ess_tail = as.numeric(
          diagnostic_summary[["ess_tail"]]
        ),
        
        inefficiency_factor =
          total_draws /
          as.numeric(
            diagnostic_summary[["ess_bulk"]]
          ),
        
        n_chains = n_chains,
        draws_per_chain =
          draws_per_chain,
        total_draws = total_draws,
        stringsAsFactors = FALSE
      )
    }
  )
  
  diagnostics <- do.call(
    rbind,
    rows
  )
  
  rownames(diagnostics) <- NULL
  
  diagnostics
}