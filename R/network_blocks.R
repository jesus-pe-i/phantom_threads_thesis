# Block-level network transformations for VAR coefficient matrices.
#
# Converts VAR lag matrices into directed unit-level block strengths
# and aggregates lag-specific networks to the unit level.


# Block strength -----

block_strength <- function(B) {
  
  B <- as.matrix(B)
  
  sqrt(
    sum(B^2)
  ) / nrow(B)
}


# Lag-level strengths -----

A_list_to_s_lag <- function(
    A_list,
    n_units,
    m) {
  
  if (!is.list(A_list) || length(A_list) == 0L) {
    stop("A_list must be a non-empty list")
  }
  
  n_series <- n_units * m
  p_lags <- length(A_list)
  
  A_list <- lapply(
    A_list,
    as.matrix
  )
  
  if (
    any(vapply(
      A_list,
      function(A) {
        !identical(
          dim(A),
          c(n_series, n_series)
        )
      },
      logical(1L)
    ))
  ) {
    stop("All lag matrices must have dimensions n_units * m")
  }
  
  s_lag <- array(
    0,
    dim = c(
      n_units,
      n_units,
      p_lags
    )
  )
  
  for (ell in seq_len(p_lags)) {
    
    A <- A_list[[ell]]
    
    for (receiver in seq_len(n_units)) {
      
      rows <- (
        (receiver - 1L) * m + 1L
      ):(
        receiver * m
      )
      
      for (sender in seq_len(n_units)) {
        
        cols <- (
          (sender - 1L) * m + 1L
        ):(
          sender * m
        )
        
        s_lag[
          receiver,
          sender,
          ell
        ] <- block_strength(
          A[
            rows,
            cols,
            drop = FALSE
          ]
        )
      }
    }
  }
  
  s_lag
}


# Lag aggregation -----

s_lag_to_unit <- function(
    s_lag,
    method = c("max", "rms")) {
  
  method <- match.arg(method)
  
  if (
    length(dim(s_lag)) != 3L ||
    dim(s_lag)[1L] != dim(s_lag)[2L]
  ) {
    stop("s_lag must be a square unit x unit x lag array")
  }
  
  if (method == "max") {
    
    apply(
      s_lag,
      c(1L, 2L),
      max
    )
    
  } else {
    
    apply(
      s_lag,
      c(1L, 2L),
      function(x) {
        sqrt(
          mean(x^2)
        )
      }
    )
  }
}


# Graph aggregation -----

G_lag_to_G_unit <- function(G_lag) {
  
  if (
    length(dim(G_lag)) != 3L ||
    dim(G_lag)[1L] != dim(G_lag)[2L]
  ) {
    stop("G_lag must be a square unit x unit x lag array")
  }
  
  apply(
    G_lag,
    c(1L, 2L),
    any
  )
}