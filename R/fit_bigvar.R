# BigVAR HLAG-OO benchmark for block VAR network recovery.
#
# Selects the HLAG-OO penalty by rolling cross-validation, refits at the
# selected lambda, and converts the estimate to the common network form.


# Coefficient extraction -----

extract_bigvar_coefficients <- function(
    beta,
    n_series,
    p_lags) {
  
  beta <- beta[, , 1L, drop = TRUE]
  
  expected_columns <-
    n_series *
    p_lags +
    1L
  
  if (!identical(
    dim(beta),
    c(
      n_series,
      expected_columns
    )
  )) {
    stop(
      "Unexpected BigVAR coefficient dimensions: ",
      paste(
        dim(beta),
        collapse = " x "
      )
    )
  }
  
  bigvar_beta <- beta[
    ,
    -1L,
    drop = FALSE
  ]
  
  A_hat_lag <- vector(
    "list",
    p_lags
  )
  
  for (lag in seq_len(p_lags)) {
    
    columns <- (
      (lag - 1L) *
        n_series +
        1L
    ):(
      lag *
        n_series
    )
    
    A_hat_lag[[lag]] <- bigvar_beta[
      ,
      columns,
      drop = FALSE
    ]
  }
  
  list(
    A_hat_lag =
      A_hat_lag,
    
    bigvar_beta =
      bigvar_beta
  )
}


# Model fit -----

fit_bigvar <- function(
    Y_list,
    p_lags,
    struct = "HLAGOO",
    gran = c(200, 20),
    h = 1L,
    cv = "Rolling",
    IC = FALSE,
    T1 = NULL,
    T2 = NULL,
    tol = 1e-4,
    verbose = FALSE,
    selected_eps = 1e-12,
    keep_cv_fit = FALSE,
    ...) {
  
  if (!requireNamespace(
    "BigVAR",
    quietly = TRUE
  )) {
    stop(
      "Package 'BigVAR' is required"
    )
  }
  
  
  # Data -----
  
  prepared <- prepare_var_data(
    Y_list = Y_list,
    p_lags = p_lags,
    standardize = FALSE
  )
  
  Y <- prepared$Y_full
  
  colnames(Y) <- paste0(
    "Y",
    seq_len(
      prepared$n_series
    )
  )
  
  if (is.null(T1)) {
    T1 <- floor(
      nrow(Y) /
        3
    )
  }
  
  if (is.null(T2)) {
    T2 <- floor(
      2 *
        nrow(Y) /
        3
    )
  }
  
  
  # Cross-validation -----
  
  model <- BigVAR::constructModel(
    Y = Y,
    p = prepared$p_lags,
    struct = struct,
    gran = gran,
    h = h,
    cv = cv,
    IC = IC,
    T1 = T1,
    T2 = T2,
    verbose = verbose,
    model.controls = list(
      intercept = FALSE,
      tol = tol
    ),
    ...
  )
  
  cv_fit <- BigVAR::cv.BigVAR(
    model
  )
  
  optimal_lambda <- as.numeric(
    cv_fit@OptimalLambda
  )
  
  lambda_grid <- as.numeric(
    cv_fit@LambdaGrid
  )
  
  lambda_index <- which.min(
    abs(
      lambda_grid -
        optimal_lambda
    )
  )
  
  
  # Final fit -----
  
  final_beta <- BigVAR::BigVAR.fit(
    Y = Y,
    p = prepared$p_lags,
    struct = struct,
    lambda = optimal_lambda,
    intercept = FALSE,
    tol = tol
  )
  
  extracted <- extract_bigvar_coefficients(
    beta = final_beta,
    n_series = prepared$n_series,
    p_lags = prepared$p_lags
  )
  
  A_hat_lag <-
    extracted$A_hat_lag
  
  beta_hat <- A_list_to_beta(
    A_hat_lag
  )
  
  
  # Network estimate -----
  
  s_hat_lag <- A_list_to_s_lag(
    A_list = A_hat_lag,
    n_units = prepared$n_units,
    m = prepared$m
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
  
  out <- list(
    model =
      "bigvar_hlag",
    
    backend =
      "BigVAR",
    
    beta_hat =
      beta_hat,
    
    bigvar_beta =
      extracted$bigvar_beta,
    
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
    
    struct =
      struct,
    
    gran =
      gran,
    
    h =
      as.integer(h),
    
    cv =
      cv,
    
    IC =
      IC,
    
    T1 =
      as.integer(T1),
    
    T2 =
      as.integer(T2),
    
    tol =
      as.numeric(tol),
    
    optimal_lambda =
      optimal_lambda,
    
    lambda_grid =
      lambda_grid,
    
    lambda_index =
      as.integer(
        lambda_index
      ),
    
    lambda_boundary =
      lambda_index %in%
      c(
        1L,
        length(
          lambda_grid
        )
      ),
    
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
      standardize =
        FALSE,
      
      centering =
        FALSE,
      
      scale =
        prepared$scale,
      
      series_names =
        prepared$series_names
    ),
    
    control = list(
      struct =
        struct,
      
      gran =
        gran,
      
      h =
        as.integer(h),
      
      cv =
        cv,
      
      IC =
        IC,
      
      T1 =
        as.integer(T1),
      
      T2 =
        as.integer(T2),
      
      tol =
        as.numeric(tol),
      
      selected_eps =
        as.numeric(
          selected_eps
        )
    )
  )
  
  if (keep_cv_fit) {
    out$cv_fit <- cv_fit
  }
  
  out
}