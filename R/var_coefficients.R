# Coefficient transformations for VAR models.
#
# Converts between regression-oriented beta matrices and VAR lag matrices,
# and restores coefficients from scale-standardized to original units.


# Scale back-transformation -----

backtransform_var_beta <- function(
    beta,
    scale,
    p_lags) {
  
  beta <- as.matrix(beta)
  scale <- as.numeric(scale)
  p_lags <- as.integer(p_lags)
  
  n_series <- ncol(beta)
  
  if (
    length(scale) != n_series ||
    nrow(beta) != n_series * p_lags
  ) {
    stop("Inconsistent beta, scale, and lag dimensions")
  }
  
  if (
    any(!is.finite(scale)) ||
    any(scale <= 0)
  ) {
    stop("scale must contain finite positive values")
  }
  
  predictor_scale <- rep(
    scale,
    times = p_lags
  )
  
  beta <- sweep(
    beta,
    1L,
    predictor_scale,
    "/"
  )
  
  sweep(
    beta,
    2L,
    scale,
    "*"
  )
}


# Beta to VAR matrices -----

beta_to_A_list <- function(
    beta,
    p_lags,
    series_names = NULL) {
  
  beta <- as.matrix(beta)
  p_lags <- as.integer(p_lags)
  
  n_series <- ncol(beta)
  
  if (nrow(beta) != n_series * p_lags) {
    stop("beta must have n_series * p_lags rows")
  }
  
  A_list <- lapply(
    seq_len(p_lags),
    function(ell) {
      
      rows <- (
        (ell - 1L) * n_series + 1L
      ):(
        ell * n_series
      )
      
      t(
        beta[
          rows,
          ,
          drop = FALSE
        ]
      )
    }
  )
  
  names(A_list) <- paste0(
    "lag_",
    seq_len(p_lags)
  )
  
  if (!is.null(series_names)) {
    
    if (length(series_names) != n_series) {
      stop("series_names must have one name per series")
    }
    
    for (ell in seq_len(p_lags)) {
      dimnames(A_list[[ell]]) <- list(
        response = series_names,
        predictor = series_names
      )
    }
  }
  
  A_list
}


# VAR matrices to beta -----

A_list_to_beta <- function(A_list) {
  
  if (!is.list(A_list) || length(A_list) == 0L) {
    stop("A_list must be a non-empty list")
  }
  
  A_list <- lapply(
    A_list,
    as.matrix
  )
  
  dimensions <- lapply(
    A_list,
    dim
  )
  
  if (
    any(vapply(
      dimensions,
      function(x) length(x) != 2L || x[1L] != x[2L],
      logical(1L)
    )) ||
    any(vapply(
      dimensions,
      function(x) !identical(x, dimensions[[1L]]),
      logical(1L)
    ))
  ) {
    stop("All lag matrices must be square and have the same dimensions")
  }
  
  do.call(
    rbind,
    lapply(
      A_list,
      t
    )
  )
}