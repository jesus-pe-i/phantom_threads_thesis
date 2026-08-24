# Adelie group elastic-net benchmark for block VAR network recovery.
#
# Fits one multivariate group elastic-net model per receiver unit using
# sender-unit by lag predictor groups and converts the result to the
# common block-network representation.


# Cross-validation folds -----

make_adelie_block_folds <- function(
    n_obs,
    n_folds = 5L) {
  
  floor(
    (seq_len(n_obs) - 1L) *
      n_folds /
      n_obs
  ) + 1L
}


# Sender-lag groups -----

make_adelie_sender_lag_groups <- function(
    n_units,
    m,
    n_series,
    p_lags) {
  
  n_groups <-
    n_units *
    p_lags
  
  groups <- integer(
    n_groups
  )
  
  sender_unit <- integer(
    n_groups
  )
  
  lag <- integer(
    n_groups
  )
  
  group <- 0L
  
  for (ell in seq_len(p_lags)) {
    
    for (sender in seq_len(n_units)) {
      
      group <- group + 1L
      
      groups[group] <- (
        ell - 1L
      ) * n_series +
        (
          sender - 1L
        ) * m +
        1L
      
      sender_unit[group] <-
        sender
      
      lag[group] <-
        ell
    }
  }
  
  list(
    groups =
      groups,
    
    sender_unit =
      sender_unit,
    
    lag =
      lag,
    
    group_size =
      rep(
        m,
        n_groups
      )
  )
}


# Receiver-unit fit -----

fit_adelie_receiver_unit <- function(
    X,
    Y,
    receiver_unit,
    n_units,
    m,
    n_series,
    p_lags,
    foldid,
    groups,
    alpha = 0.5,
    lambda_rule = c(
      "lambda.min",
      "lambda.1se"
    ),
    standardize = TRUE,
    intercept = FALSE,
    min_ratio = 0.01,
    lmda_path_size = 100L,
    n_threads = 1L,
    ...) {
  
  if (!requireNamespace(
    "adelie",
    quietly = TRUE
  )) {
    stop(
      "Package 'adelie' is required"
    )
  }
  
  lambda_rule <- match.arg(
    lambda_rule
  )
  
  n_coef <-
    n_series *
    p_lags
  
  response_columns <- (
    (receiver_unit - 1L) *
      m +
      1L
  ):(
    receiver_unit *
      m
  )
  
  Y_receiver <- Y[
    ,
    response_columns,
    drop = FALSE
  ]
  
  fit <- adelie::cv.grpnet(
    X = X,
    glm = adelie::glm.multigaussian(
      Y_receiver
    ),
    groups = groups$groups,
    alpha = alpha,
    foldid = foldid,
    n_folds = max(foldid),
    standardize = standardize,
    intercept = intercept,
    min_ratio = min_ratio,
    lmda_path_size = lmda_path_size,
    n_threads = n_threads,
    progress_bar = FALSE,
    ...
  )
  
  coefficients <- stats::coef(
    fit,
    lambda = lambda_rule
  )
  
  beta_vector <- as.numeric(
    coefficients$betas
  )
  
  beta_receiver <- matrix(
    beta_vector,
    nrow = n_coef,
    ncol = m,
    byrow = TRUE
  )
  
  intercept_receiver <- if (
    intercept
  ) {
    
    tail(
      as.numeric(
        coefficients$intercepts
      ),
      m
    )
    
  } else {
    
    rep(
      0,
      m
    )
  }
  
  active_groups <- stats::predict(
    fit,
    lambda = lambda_rule,
    type = "nonzero"
  )
  
  if (is.list(active_groups)) {
    active_groups <- active_groups[[1L]]
  }
  
  lambda_used <- fit[[lambda_rule]]
  
  list(
    beta =
      beta_receiver,
    
    intercept =
      intercept_receiver,
    
    active_groups =
      as.integer(
        active_groups
      ),
    
    lambda_used =
      lambda_used,
    
    lambda_min =
      fit$lambda.min,
    
    lambda_1se =
      fit$lambda.1se,
    
    lambda_index =
      which.min(
        abs(
          fit$lambda -
            lambda_used
        )
      ),
    
    lambda_path =
      fit$lambda,
    
    cvm =
      fit$cvm,
    
    cvsd =
      fit$cvsd,
    
    fit =
      fit
  )
}


# Coefficient fitting -----

fit_adelie_coefficients <- function(
    prepared,
    alpha = 0.5,
    lambda_rule = c(
      "lambda.min",
      "lambda.1se"
    ),
    n_folds = 5L,
    standardize = TRUE,
    intercept = FALSE,
    min_ratio = 0.01,
    lmda_path_size = 100L,
    n_threads = 1L,
    seed = NULL,
    keep_fits = FALSE,
    ...) {
  
  lambda_rule <- match.arg(
    lambda_rule
  )
  
  groups <- make_adelie_sender_lag_groups(
    n_units =
      prepared$n_units,
    m =
      prepared$m,
    n_series =
      prepared$n_series,
    p_lags =
      prepared$p_lags
  )
  
  foldid <- make_adelie_block_folds(
    n_obs =
      prepared$T_eff,
    n_folds =
      n_folds
  )
  
  beta_hat <- matrix(
    0,
    nrow =
      prepared$n_coef,
    ncol =
      prepared$n_series
  )
  
  selected_beta <- matrix(
    0,
    nrow =
      prepared$n_coef,
    ncol =
      prepared$n_series
  )
  
  intercept_hat <- matrix(
    0,
    nrow =
      prepared$n_units,
    ncol =
      prepared$m
  )
  
  active_groups <- vector(
    "list",
    prepared$n_units
  )
  
  lambda_used <- numeric(
    prepared$n_units
  )
  
  lambda_min <- numeric(
    prepared$n_units
  )
  
  lambda_1se <- numeric(
    prepared$n_units
  )
  
  lambda_index <- integer(
    prepared$n_units
  )
  
  path_length <- integer(
    prepared$n_units
  )
  
  fits <- if (
    keep_fits
  ) {
    
    vector(
      "list",
      prepared$n_units
    )
    
  } else {
    
    NULL
  }
  
  
  # Receiver-unit models -----
  
  for (
    receiver in seq_len(
      prepared$n_units
    )
  ) {
    
    fit_receiver <- fit_adelie_receiver_unit(
      X =
        prepared$X,
      Y =
        prepared$Y,
      receiver_unit =
        receiver,
      n_units =
        prepared$n_units,
      m =
        prepared$m,
      n_series =
        prepared$n_series,
      p_lags =
        prepared$p_lags,
      foldid =
        foldid,
      groups =
        groups,
      alpha =
        alpha,
      lambda_rule =
        lambda_rule,
      standardize =
        standardize,
      intercept =
        intercept,
      min_ratio =
        min_ratio,
      lmda_path_size =
        lmda_path_size,
      n_threads =
        n_threads,
      ...
    )
    
    response_columns <- (
      (receiver - 1L) *
        prepared$m +
        1L
    ):(
      receiver *
        prepared$m
    )
    
    beta_hat[
      ,
      response_columns
    ] <- fit_receiver$beta
    
    intercept_hat[
      receiver,
    ] <- fit_receiver$intercept
    
    active_groups[[receiver]] <-
      fit_receiver$active_groups
    
    for (
      group in fit_receiver$active_groups
    ) {
      
      predictor_rows <-
        groups$groups[group] +
        seq_len(
          prepared$m
        ) -
        1L
      
      selected_beta[
        predictor_rows,
        response_columns
      ] <- 1
    }
    
    lambda_used[receiver] <-
      fit_receiver$lambda_used
    
    lambda_min[receiver] <-
      fit_receiver$lambda_min
    
    lambda_1se[receiver] <-
      fit_receiver$lambda_1se
    
    lambda_index[receiver] <-
      fit_receiver$lambda_index
    
    path_length[receiver] <-
      length(
        fit_receiver$lambda_path
      )
    
    if (keep_fits) {
      fits[[receiver]] <- fit_receiver
    }
  }
  
  out <- list(
    beta_hat =
      beta_hat,
    
    selected_beta =
      selected_beta,
    
    intercept =
      intercept_hat,
    
    active_groups =
      active_groups,
    
    groups =
      groups,
    
    lambda_used =
      lambda_used,
    
    lambda_min =
      lambda_min,
    
    lambda_1se =
      lambda_1se,
    
    lambda_index =
      lambda_index,
    
    path_length =
      path_length,
    
    alpha =
      alpha,
    
    lambda_rule =
      lambda_rule,
    
    foldid =
      foldid,
    
    n_folds =
      as.integer(
        n_folds
      ),
    
    standardize =
      standardize,
    
    intercept_enabled =
      intercept,
    
    min_ratio =
      min_ratio,
    
    lmda_path_size =
      as.integer(
        lmda_path_size
      ),
    
    n_threads =
      as.integer(
        n_threads
      ),
    
    fit_seed =
      seed
  )
  
  if (keep_fits) {
    out$fits <- fits
  }
  
  out
}


# Model fit -----

fit_adelie <- function(
    Y_list,
    p_lags,
    alpha = 0.5,
    lambda_rule = c(
      "lambda.min",
      "lambda.1se"
    ),
    n_folds = 5L,
    standardize = TRUE,
    intercept = FALSE,
    min_ratio = 0.01,
    lmda_path_size = 100L,
    n_threads = 1L,
    seed = NULL,
    keep_fits = FALSE,
    ...) {
  
  lambda_rule <- match.arg(
    lambda_rule
  )
  
  
  # Data -----
  
  prepared <- prepare_var_data(
    Y_list =
      Y_list,
    p_lags =
      p_lags,
    standardize =
      FALSE
  )
  
  
  # Adelie fits -----
  
  coefficient_fit <- fit_adelie_coefficients(
    prepared =
      prepared,
    alpha =
      alpha,
    lambda_rule =
      lambda_rule,
    n_folds =
      n_folds,
    standardize =
      standardize,
    intercept =
      intercept,
    min_ratio =
      min_ratio,
    lmda_path_size =
      lmda_path_size,
    n_threads =
      n_threads,
    seed =
      seed,
    keep_fits =
      keep_fits,
    ...
  )
  
  
  # VAR coefficients -----
  
  A_hat_lag <- beta_to_A_list(
    beta =
      coefficient_fit$beta_hat,
    p_lags =
      prepared$p_lags,
    series_names =
      prepared$series_names
  )
  
  selected_A_lag <- beta_to_A_list(
    beta =
      coefficient_fit$selected_beta,
    p_lags =
      prepared$p_lags
  )
  
  
  # Network estimate -----
  
  s_hat_lag <- A_list_to_s_lag(
    A_list =
      A_hat_lag,
    n_units =
      prepared$n_units,
    m =
      prepared$m
  )
  
  selected_lag <- A_list_to_s_lag(
    A_list =
      selected_A_lag,
    n_units =
      prepared$n_units,
    m =
      prepared$m
  ) > 0
  
  s_hat_unit_max <- s_lag_to_unit(
    s_lag =
      s_hat_lag,
    method =
      "max"
  )
  
  s_hat_unit_rms <- s_lag_to_unit(
    s_lag =
      s_hat_lag,
    method =
      "rms"
  )
  
  
  # Output -----
  
  out <- list(
    model =
      "adelie_gen",
    
    backend =
      "adelie",
    
    beta_hat =
      coefficient_fit$beta_hat,
    
    selected_beta =
      coefficient_fit$selected_beta,
    
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
      selected_lag,
    
    intercept =
      coefficient_fit$intercept,
    
    active_groups =
      coefficient_fit$active_groups,
    
    groups =
      coefficient_fit$groups,
    
    lambda_used =
      coefficient_fit$lambda_used,
    
    lambda_min =
      coefficient_fit$lambda_min,
    
    lambda_1se =
      coefficient_fit$lambda_1se,
    
    lambda_index =
      coefficient_fit$lambda_index,
    
    path_length =
      coefficient_fit$path_length,
    
    alpha =
      coefficient_fit$alpha,
    
    lambda_rule =
      coefficient_fit$lambda_rule,
    
    foldid =
      coefficient_fit$foldid,
    
    n_folds =
      coefficient_fit$n_folds,
    
    fit_seed =
      coefficient_fit$fit_seed,
    
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
      input_standardize =
        FALSE,
      
      adelie_standardize =
        standardize,
      
      centering =
        FALSE,
      
      input_scale =
        prepared$scale,
      
      series_names =
        prepared$series_names
    ),
    
    control = list(
      alpha =
        as.numeric(
          alpha
        ),
      
      lambda_rule =
        lambda_rule,
      
      n_folds =
        as.integer(
          n_folds
        ),
      
      standardize =
        standardize,
      
      intercept =
        intercept,
      
      min_ratio =
        as.numeric(
          min_ratio
        ),
      
      lmda_path_size =
        as.integer(
          lmda_path_size
        ),
      
      n_threads =
        as.integer(
          n_threads
        )
    )
  )
  
  if (keep_fits) {
    out$fits <- coefficient_fit$fits
  }
  
  out
}