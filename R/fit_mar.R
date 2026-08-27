# Matrix autoregressive benchmark for block VAR network recovery.
#
# Fits a separable MAR(p) coefficient structure by multistart alternating
# least squares and converts the estimate to the common block-network form.


# Data construction -----

make_mar_array <- function(
    Y,
    n_units,
    m) {
  
  T_obs <- nrow(Y)
  
  Y_array <- array(
    NA_real_,
    dim = c(
      T_obs,
      n_units,
      m
    )
  )
  
  for (unit in seq_len(n_units)) {
    
    columns <- (
      (unit - 1L) * m + 1L
    ):(
      unit * m
    )
    
    Y_array[, unit, ] <-
      Y[
        ,
        columns,
        drop = FALSE
      ]
  }
  
  Y_array
}


make_mar_lags <- function(
    Y_array,
    p_lags) {
  
  T_obs <- dim(Y_array)[1L]
  n_units <- dim(Y_array)[2L]
  m <- dim(Y_array)[3L]
  
  T_eff <- T_obs - p_lags
  
  current <- Y_array[
    (p_lags + 1L):T_obs,
    ,
    ,
    drop = FALSE
  ]
  
  lags <- array(
    NA_real_,
    dim = c(
      T_eff,
      n_units * p_lags,
      m
    )
  )
  
  for (t in seq_len(T_eff)) {
    
    for (lag in seq_len(p_lags)) {
      
      columns <- (
        (lag - 1L) *
          n_units +
          1L
      ):(
        lag *
          n_units
      )
      
      lags[t, columns, ] <-
        Y_array[
          t + p_lags - lag,
          ,
        ]
    }
  }
  
  list(
    current = current,
    lags = lags,
    T_eff = T_eff,
    n_units = n_units,
    m = m,
    p_lags = p_lags
  )
}


# Least squares -----

mar_lstsq <- function(
    X,
    Y,
    tol = 1e-10) {
  
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  rank_max <- min(
    dim(X)
  )
  
  decomposition <- svd(
    X,
    nu = rank_max,
    nv = rank_max
  )
  
  cutoff <-
    tol *
    max(
      dim(X)
    ) *
    max(
      decomposition$d
    )
  
  keep <- which(
    decomposition$d >
      cutoff
  )
  
  coefficients <- matrix(
    0,
    nrow = ncol(X),
    ncol = ncol(Y)
  )
  
  if (length(keep) > 0L) {
    
    projection <- crossprod(
      decomposition$u[
        ,
        keep,
        drop = FALSE
      ],
      Y
    )
    
    projection <- sweep(
      projection,
      1L,
      decomposition$d[keep],
      "/"
    )
    
    coefficients <-
      decomposition$v[
        ,
        keep,
        drop = FALSE
      ] %*%
      projection
  }
  
  list(
    coefficients = coefficients,
    
    rank = as.integer(
      length(keep)
    ),
    
    condition = if (
      length(keep) > 0L
    ) {
      
      max(
        decomposition$d
      ) /
        min(
          decomposition$d[keep]
        )
      
    } else {
      
      NA_real_
    }
  )
}


# Objective -----

evaluate_mar <- function(
    A,
    B,
    current,
    lags) {
  
  T_eff <- dim(current)[1L]
  n_units <- dim(current)[2L]
  m <- dim(current)[3L]
  
  n_predictors <- dim(lags)[2L]
  
  fitted <- array(
    0,
    dim = dim(current)
  )
  
  for (t in seq_len(T_eff)) {
    
    X_t <- matrix(
      lags[t, , ],
      nrow = n_predictors,
      ncol = m
    )
    
    fitted[t, , ] <-
      A %*%
      X_t %*%
      t(B)
  }
  
  residuals <-
    current -
    fitted
  
  list(
    fitted = fitted,
    residuals = residuals,
    rss = sum(
      residuals^2
    )
  )
}


# Alternating least squares -----

fit_mar_core <- function(
    Y_array,
    p_lags,
    max_iter = 200L,
    tol = 1e-8,
    solver_tol = 1e-10,
    init_B = NULL) {
  
  lag_data <- make_mar_lags(
    Y_array = Y_array,
    p_lags = p_lags
  )
  
  current <- lag_data$current
  lags <- lag_data$lags
  
  T_eff <- lag_data$T_eff
  n_units <- lag_data$n_units
  m <- lag_data$m
  
  n_predictors <-
    n_units *
    p_lags
  
  
  # Initial state -----
  
  if (is.null(init_B)) {
    
    B <- diag(m)
    
  } else {
    
    B <- as.matrix(
      init_B
    )
  }
  
  B <-
    B *
    sqrt(m) /
    sqrt(
      sum(
        B^2
      )
    )
  
  A <- matrix(
    0,
    nrow = n_units,
    ncol = n_predictors
  )
  
  A[
    ,
    seq_len(n_units)
  ] <- diag(
    n_units
  )
  
  initial <- evaluate_mar(
    A = A,
    B = B,
    current = current,
    lags = lags
  )
  
  objective <- rep(
    NA_real_,
    max_iter + 1L
  )
  
  objective[1L] <-
    initial$rss
  
  converged <- FALSE
  
  
  # Alternating updates -----
  
  for (iter in seq_len(max_iter)) {
    
    
    ## Unit-lag matrix -----
    
    X_A <- matrix(
      NA_real_,
      nrow = T_eff * m,
      ncol = n_predictors
    )
    
    Y_A <- matrix(
      NA_real_,
      nrow = T_eff * m,
      ncol = n_units
    )
    
    for (t in seq_len(T_eff)) {
      
      X_t <- matrix(
        lags[t, , ],
        nrow = n_predictors,
        ncol = m
      )
      
      Y_t <- matrix(
        current[t, , ],
        nrow = n_units,
        ncol = m
      )
      
      rows <- (
        (t - 1L) * m + 1L
      ):(
        t * m
      )
      
      X_A[rows, ] <-
        t(
          X_t %*%
            t(B)
        )
      
      Y_A[rows, ] <-
        t(Y_t)
    }
    
    A_fit <- mar_lstsq(
      X = X_A,
      Y = Y_A,
      tol = solver_tol
    )
    
    A <- t(
      A_fit$coefficients
    )
    
    
    ## Within-unit matrix -----
    
    X_B <- matrix(
      NA_real_,
      nrow =
        T_eff *
        n_units,
      ncol = m
    )
    
    Y_B <- matrix(
      NA_real_,
      nrow =
        T_eff *
        n_units,
      ncol = m
    )
    
    for (t in seq_len(T_eff)) {
      
      X_t <- matrix(
        lags[t, , ],
        nrow = n_predictors,
        ncol = m
      )
      
      Y_t <- matrix(
        current[t, , ],
        nrow = n_units,
        ncol = m
      )
      
      rows <- (
        (t - 1L) *
          n_units +
          1L
      ):(
        t *
          n_units
      )
      
      X_B[rows, ] <-
        A %*%
        X_t
      
      Y_B[rows, ] <-
        Y_t
    }
    
    B_fit <- mar_lstsq(
      X = X_B,
      Y = Y_B,
      tol = solver_tol
    )
    
    B <- t(
      B_fit$coefficients
    )
    
    
    ## Scale normalization -----
    
    B_norm <- sqrt(
      sum(
        B^2
      )
    )
    
    if (B_norm == 0) {
      stop(
        "MAR B estimate collapsed to zero"
      )
    }
    
    scale <- sqrt(m) /
      B_norm
    
    B <- B *
      scale
    
    A <- A /
      scale
    
    
    ## Convergence -----
    
    current_fit <- evaluate_mar(
      A = A,
      B = B,
      current = current,
      lags = lags
    )
    
    objective[iter + 1L] <-
      current_fit$rss
    
    relative_change <- abs(
      objective[iter + 1L] -
        objective[iter]
    ) /
      max(
        1,
        objective[iter]
      )
    
    if (
      relative_change <
      tol
    ) {
      
      converged <- TRUE
      break
    }
  }
  
  
  # Final state -----
  
  objective <- objective[
    seq_len(
      iter + 1L
    )
  ]
  
  final_fit <- evaluate_mar(
    A = A,
    B = B,
    current = current,
    lags = lags
  )
  
  final_relative_change <- abs(
    objective[
      length(objective)
    ] -
      objective[
        length(objective) -
          1L
      ]
  ) /
    max(
      1,
      objective[
        length(objective) -
          1L
      ]
    )
  
  list(
    A = A,
    B = B,
    
    fitted =
      final_fit$fitted,
    
    residuals =
      final_fit$residuals,
    
    objective =
      objective,
    
    initial_rss =
      objective[1L],
    
    final_rss =
      final_fit$rss,
    
    final_relative_change =
      final_relative_change,
    
    iter =
      as.integer(iter),
    
    converged =
      converged,
    
    A_rank =
      A_fit$rank,
    
    B_rank =
      B_fit$rank,
    
    A_condition =
      A_fit$condition,
    
    B_condition =
      B_fit$condition
  )
}


# Multistart fitting -----

fit_mar_multistart <- function(
    Y_array,
    p_lags,
    n_starts = 3L,
    seed = 991L,
    max_iter = 200L,
    tol = 1e-8,
    solver_tol = 1e-10) {
  
  n_starts <- as.integer(
    n_starts
  )
  
  if (n_starts < 1L) {
    stop(
      "n_starts must be positive"
    )
  }
  
  m <- dim(Y_array)[3L]
  
  set.seed(
    seed
  )
  
  initial_B <- vector(
    "list",
    n_starts
  )
  
  initial_B[[1L]] <-
    diag(m)
  
  if (n_starts > 1L) {
    
    for (
      start in 2:n_starts
    ) {
      
      B <- matrix(
        rnorm(
          m * m
        ),
        nrow = m,
        ncol = m
      )
      
      initial_B[[start]] <-
        B *
        sqrt(m) /
        sqrt(
          sum(
            B^2
          )
        )
    }
  }
  
  fits <- lapply(
    initial_B,
    function(B) {
      
      fit_mar_core(
        Y_array = Y_array,
        p_lags = p_lags,
        max_iter = max_iter,
        tol = tol,
        solver_tol = solver_tol,
        init_B = B
      )
    }
  )
  
  start_rss <- vapply(
    fits,
    `[[`,
    numeric(1L),
    "final_rss"
  )
  
  best_start <- which.min(
    start_rss
  )
  
  fit <- fits[[
      best_start
    ]]
  
  fit$n_starts <-
    n_starts
  
  fit$best_start <-
    as.integer(
      best_start
    )
  
  fit$start_rss <-
    start_rss
  
  fit$start_iter <- vapply(
    fits,
    `[[`,
    integer(1L),
    "iter"
  )
  
  fit$start_converged <- vapply(
    fits,
    `[[`,
    logical(1L),
    "converged"
  )
  
  fit
}


# VAR reconstruction -----

mar_to_A_list <- function(
    A,
    B,
    n_units,
    m,
    p_lags) {
  
  A_list <- vector(
    "list",
    p_lags
  )
  
  N <- n_units *
    m
  
  for (lag in seq_len(p_lags)) {
    
    A_unit <- A[
      ,
      (
        (lag - 1L) *
          n_units +
          1L
      ):(
        lag *
          n_units
      ),
      drop = FALSE
    ]
    
    A_lag <- matrix(
      0,
      nrow = N,
      ncol = N
    )
    
    for (
      receiver in seq_len(n_units)
    ) {
      
      rows <- (
        (receiver - 1L) *
          m +
          1L
      ):(
        receiver *
          m
      )
      
      for (
        sender in seq_len(n_units)
      ) {
        
        columns <- (
          (sender - 1L) *
            m +
            1L
        ):(
          sender *
            m
        )
        
        A_lag[
          rows,
          columns
        ] <-
          A_unit[
            receiver,
            sender
          ] *
          B
      }
    }
    
    A_list[[lag]] <-
      A_lag
  }
  
  A_list
}


# Model fit -----

fit_mar <- function(
    Y_list,
    p_lags,
    n_starts = 3L,
    seed = 991L,
    max_iter = 200L,
    tol = 1e-8,
    solver_tol = 1e-10,
    selected_eps = 1e-12) {
  
  
  # Data -----
  
  prepared <- prepare_var_data(
    Y_list = Y_list,
    p_lags = p_lags,
    standardize = FALSE
  )
  
  Y_array <- make_mar_array(
    Y = prepared$Y_full,
    n_units =
      prepared$n_units,
    m =
      prepared$m
  )
  
  
  # Estimation -----
  
  mar_fit <- fit_mar_multistart(
    Y_array = Y_array,
    p_lags =
      prepared$p_lags,
    n_starts = n_starts,
    seed = seed,
    max_iter = max_iter,
    tol = tol,
    solver_tol = solver_tol
  )
  
  
  # VAR coefficients -----
  
  A_hat_lag <- mar_to_A_list(
    A = mar_fit$A,
    B = mar_fit$B,
    n_units =
      prepared$n_units,
    m =
      prepared$m,
    p_lags =
      prepared$p_lags
  )
  
  beta_hat <- A_list_to_beta(
    A_hat_lag
  )
  
  
  # Network estimate -----
  
  s_hat_lag <- A_list_to_s_lag(
    A_list = A_hat_lag,
    n_units =
      prepared$n_units,
    m =
      prepared$m
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
  
  list(
    model = "mar",
    backend = "mar",
    
    beta_hat = beta_hat,
    
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
    
    mar_A =
      mar_fit$A,
    
    mar_B =
      mar_fit$B,
    
    fitted =
      mar_fit$fitted,
    
    residuals =
      mar_fit$residuals,
    
    objective =
      mar_fit$objective,
    
    initial_rss =
      mar_fit$initial_rss,
    
    final_rss =
      mar_fit$final_rss,
    
    iter =
      mar_fit$iter,
    
    converged =
      mar_fit$converged,
    
    final_relative_change =
      mar_fit$final_relative_change,
    
    n_starts =
      mar_fit$n_starts,
    
    best_start =
      mar_fit$best_start,
    
    start_rss =
      mar_fit$start_rss,
    
    start_iter =
      mar_fit$start_iter,
    
    start_converged =
      mar_fit$start_converged,
    
    A_rank =
      mar_fit$A_rank,
    
    B_rank =
      mar_fit$B_rank,
    
    A_condition =
      mar_fit$A_condition,
    
    B_condition =
      mar_fit$B_condition,
    
    fit_seed =
      as.integer(seed),
    
    dimensions = list(
      n_units =
        prepared$n_units,
      
      m =
        prepared$m,
      
      N =
        prepared$n_series,
      
      p_lags =
        prepared$p_lags,
      
      k =
        prepared$n_coef,
      
      T_full =
        prepared$T_full,
      
      T_p =
        prepared$T_eff
    ),
    
    preprocessing = list(
      standardize = FALSE,
      centering = FALSE,
      scale =
        prepared$scale,
      series_names =
        prepared$series_names
    ),
    
    control = list(
      n_starts =
        as.integer(n_starts),
      
      max_iter =
        as.integer(max_iter),
      
      tol =
        as.numeric(tol),
      
      solver_tol =
        as.numeric(
          solver_tol
        ),
      
      selected_eps =
        as.numeric(
          selected_eps
        )
    )
  )
}