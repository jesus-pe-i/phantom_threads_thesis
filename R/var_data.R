# Core data preparation for block VAR models.
#
# Stacks unit-level series, applies optional scale-only standardization,
# and constructs the response and lagged predictor matrices.


# Unit data -----

Y_list_to_matrix <- function(Y_list) {
  
  if (!is.list(Y_list) || length(Y_list) == 0L) {
    stop("Y_list must be a non-empty list")
  }
  
  Y_list <- lapply(
    Y_list,
    function(Y) {
      Y <- as.matrix(Y)
      storage.mode(Y) <- "double"
      Y
    }
  )
  
  n_rows <- vapply(Y_list, nrow, integer(1L))
  n_cols <- vapply(Y_list, ncol, integer(1L))
  
  if (
    any(n_rows != n_rows[1L]) ||
    any(n_cols != n_cols[1L])
  ) {
    stop("All unit matrices must have the same dimensions")
  }
  
  Y <- do.call(cbind, Y_list)
  
  if (any(!is.finite(Y))) {
    stop("Y_list contains non-finite observations")
  }
  
  Y
}


# Scaling -----

scale_var_data <- function(Y, standardize = TRUE) {
  
  if (!standardize) {
    return(
      list(
        Y = Y,
        scale = rep(1, ncol(Y))
      )
    )
  }
  
  series_scale <- apply(
    Y,
    2L,
    stats::sd
  )
  
  if (any(!is.finite(series_scale)) || any(series_scale <= 0)) {
    stop("Every series must have a finite positive scale")
  }
  
  list(
    Y = sweep(
      Y,
      2L,
      series_scale,
      "/"
    ),
    scale = series_scale
  )
}


# VAR design -----

make_var_design_matrix <- function(Y, p_lags) {
  
  if (
    length(p_lags) != 1L ||
    !is.finite(p_lags) ||
    p_lags < 1L ||
    p_lags != as.integer(p_lags)
  ) {
    stop("p_lags must be a positive integer")
  }
  
  p_lags <- as.integer(p_lags)
  T_full <- nrow(Y)
  n_series <- ncol(Y)
  
  if (T_full <= p_lags) {
    stop("The number of observations must exceed p_lags")
  }
  
  response_rows <- (p_lags + 1L):T_full
  
  Y_response <- Y[
    response_rows,
    ,
    drop = FALSE
  ]
  
  X <- do.call(
    cbind,
    lapply(
      seq_len(p_lags),
      function(ell) {
        Y[
          response_rows - ell,
          ,
          drop = FALSE
        ]
      }
    )
  )
  
  list(
    Y = unname(Y_response),
    X = unname(X),
    T_full = T_full,
    T_eff = T_full - p_lags,
    n_series = n_series,
    p_lags = p_lags,
    n_coef = n_series * p_lags
  )
}


# Combined preparation -----

prepare_var_data <- function(
    Y_list,
    p_lags,
    standardize = TRUE) {
  
  Y_full <- Y_list_to_matrix(Y_list)
  
  n_units <- length(Y_list)
  m <- ncol(Y_list[[1L]])
  n_series <- ncol(Y_full)
  
  series_names <- colnames(Y_full)
  
  if (
    is.null(series_names) ||
    anyNA(series_names) ||
    any(series_names == "") ||
    anyDuplicated(series_names)
  ) {
    series_names <- paste0(
      "series_",
      seq_len(n_series)
    )
  }
  
  scaled <- scale_var_data(
    Y = Y_full,
    standardize = standardize
  )
  
  design <- make_var_design_matrix(
    Y = scaled$Y,
    p_lags = p_lags
  )
  
  list(
    Y = design$Y,
    X = design$X,
    Y_full = unname(Y_full),
    scale = as.numeric(scaled$scale),
    series_names = as.character(series_names),
    n_units = n_units,
    m = m,
    n_series = n_series,
    p_lags = design$p_lags,
    n_coef = design$n_coef,
    T_full = design$T_full,
    T_eff = design$T_eff,
    standardize = standardize
  )
}