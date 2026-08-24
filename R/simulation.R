# Generic Gaussian VAR simulation from fixed coefficient matrices.
#
# Validates frozen VAR coefficients, checks stability, constructs exact
# coefficient and network truth, and simulates homoscedastic Gaussian paths.


# Coefficient matrices -----

validate_A_list <- function(
    A_list) {
  
  if (
    !is.list(A_list) ||
    length(A_list) == 0L
  ) {
    stop(
      "A_list must be a non-empty list"
    )
  }
  
  dimensions <- lapply(
    A_list,
    dim
  )
  
  if (
    any(
      vapply(
        A_list,
        function(A) {
          !is.matrix(A)
        },
        logical(1L)
      )
    )
  ) {
    stop(
      "All elements of A_list must be matrices"
    )
  }
  
  N <- dimensions[[1L]][1L]
  
  valid_dimensions <- vapply(
    dimensions,
    function(dimension) {
      
      length(dimension) == 2L &&
        all(
          dimension ==
            c(
              N,
              N
            )
        )
    },
    logical(1L)
  )
  
  if (!all(valid_dimensions)) {
    stop(
      "All lag matrices must be square with common dimensions"
    )
  }
  
  if (
    any(
      vapply(
        A_list,
        function(A) {
          any(
            !is.finite(A)
          )
        },
        logical(1L)
      )
    )
  ) {
    stop(
      "A_list contains non-finite coefficients"
    )
  }
  
  invisible(
    TRUE
  )
}


# Stability -----

companion_matrix <- function(
    A_list) {
  
  validate_A_list(
    A_list
  )
  
  N <- nrow(
    A_list[[1L]]
  )
  
  p_lags <- length(
    A_list
  )
  
  top <- do.call(
    cbind,
    A_list
  )
  
  if (p_lags == 1L) {
    return(
      top
    )
  }
  
  bottom <- cbind(
    diag(
      N *
        (
          p_lags -
            1L
        )
    ),
    
    matrix(
      0,
      nrow =
        N *
        (
          p_lags -
            1L
        ),
      ncol =
        N
    )
  )
  
  rbind(
    top,
    bottom
  )
}


check_A_radius <- function(
    A_list) {
  
  companion <- companion_matrix(
    A_list
  )
  
  list(
    radius =
      max(
        Mod(
          eigen(
            companion,
            only.values = TRUE
          )$values
        )
      ),
    
    N =
      nrow(
        A_list[[1L]]
      ),
    
    p_lags =
      length(
        A_list
      )
  )
}


# Data layout -----

matrix_to_Y_list <- function(
    Y,
    n_units,
    m) {
  
  Y <- as.matrix(
    Y
  )
  
  if (
    ncol(Y) !=
    n_units *
    m
  ) {
    stop(
      "ncol(Y) must equal n_units * m"
    )
  }
  
  lapply(
    seq_len(
      n_units
    ),
    function(unit) {
      
      columns <- (
        (unit - 1L) *
          m +
          1L
      ):(
        unit *
          m
      )
      
      Y[
        ,
        columns,
        drop = FALSE
      ]
    }
  )
}


# Innovation covariance -----

make_sigma_homoscedastic <- function(
    N,
    sigma = 1) {
  
  N <- as.integer(
    N
  )
  
  sigma <- as.numeric(
    sigma
  )
  
  if (length(sigma) == 1L) {
    
    sigma <- rep(
      sigma,
      N
    )
  }
  
  if (
    length(sigma) != N ||
    any(
      !is.finite(sigma)
    ) ||
    any(
      sigma < 0
    )
  ) {
    stop(
      "sigma must contain one finite non-negative value per series"
    )
  }
  
  diag(
    sigma^2,
    nrow = N,
    ncol = N
  )
}


make_innovation_factor <- function(
    Sigma,
    N,
    tolerance = 1e-10) {
  
  Sigma <- as.matrix(
    Sigma
  )
  
  if (
    !all(
      dim(Sigma) ==
      c(
        N,
        N
      )
    ) ||
    any(
      !is.finite(Sigma)
    )
  ) {
    stop(
      "Sigma must be one finite N x N matrix"
    )
  }
  
  Sigma <- (
    Sigma +
      t(Sigma)
  ) / 2
  
  eigen_decomposition <- eigen(
    Sigma,
    symmetric = TRUE
  )
  
  if (
    min(
      eigen_decomposition$values
    ) <
    -tolerance
  ) {
    stop(
      "Sigma must be positive semidefinite"
    )
  }
  
  sweep(
    eigen_decomposition$vectors,
    2L,
    sqrt(
      pmax(
        eigen_decomposition$values,
        0
      )
    ),
    "*"
  )
}


# Truth -----

make_truth_from_A <- function(
    A_list,
    n_units,
    m,
    truth_eps = 1e-12) {
  
  validate_A_list(
    A_list
  )
  
  N <- nrow(
    A_list[[1L]]
  )
  
  p_lags <- length(
    A_list
  )
  
  if (
    N !=
    n_units *
    m
  ) {
    stop(
      "A matrix dimension must equal n_units * m"
    )
  }
  
  beta_true <- A_list_to_beta(
    A_list
  )
  
  s_true_lag <- A_list_to_s_lag(
    A_list =
      A_list,
    n_units =
      n_units,
    m =
      m
  )
  
  G_true_lag <-
    s_true_lag >
    truth_eps
  
  list(
    beta_true =
      beta_true,
    
    s_true_lag =
      s_true_lag,
    
    s_true_unit_max =
      s_lag_to_unit(
        s_lag =
          s_true_lag,
        method =
          "max"
      ),
    
    s_true_unit_rms =
      s_lag_to_unit(
        s_lag =
          s_true_lag,
        method =
          "rms"
      ),
    
    G_true_lag =
      G_true_lag,
    
    G_true_unit =
      G_lag_to_G_unit(
        G_true_lag
      ),
    
    dimensions =
      list(
        n_units =
          as.integer(
            n_units
          ),
        
        m =
          as.integer(
            m
          ),
        
        N =
          as.integer(
            N
          ),
        
        p_lags =
          as.integer(
            p_lags
          ),
        
        k =
          as.integer(
            N *
              p_lags
          )
      )
  )
}


# VAR simulation -----

simulate_var_from_A <- function(
    A_list,
    Sigma,
    T_obs,
    burn_in = 300L,
    seed = NULL,
    n_units,
    m,
    truth_eps = 1e-12,
    allow_unstable = FALSE) {
  
  validate_A_list(
    A_list
  )
  
  N <- nrow(
    A_list[[1L]]
  )
  
  p_lags <- length(
    A_list
  )
  
  T_obs <- as.integer(
    T_obs
  )
  
  burn_in <- as.integer(
    burn_in
  )
  
  n_units <- as.integer(
    n_units
  )
  
  m <- as.integer(
    m
  )
  
  if (
    T_obs < 1L ||
    burn_in < 0L
  ) {
    stop(
      "T_obs must be positive and burn_in non-negative"
    )
  }
  
  if (
    N !=
    n_units *
    m
  ) {
    stop(
      "A matrix dimension must equal n_units * m"
    )
  }
  
  stability <- check_A_radius(
    A_list
  )
  
  if (
    !allow_unstable &&
    stability$radius >= 1
  ) {
    stop(
      "A_list is unstable: radius = ",
      stability$radius
    )
  }
  
  if (!is.null(seed)) {
    set.seed(
      seed
    )
  }
  
  T_total <-
    T_obs +
    burn_in
  
  if (
    T_total <=
    p_lags
  ) {
    stop(
      "T_obs + burn_in must exceed p_lags"
    )
  }
  
  innovation_factor <-
    make_innovation_factor(
      Sigma =
        Sigma,
      N =
        N
    )
  
  innovations <- matrix(
    rnorm(
      T_total *
        N
    ),
    nrow =
      T_total,
    ncol =
      N,
    byrow =
      TRUE
  ) %*%
    t(
      innovation_factor
    )
  
  Y_full <- matrix(
    0,
    nrow =
      T_total,
    ncol =
      N
  )
  
  Y_full[
    seq_len(
      p_lags
    ),
  ] <- innovations[
    seq_len(
      p_lags
    ),
    ,
    drop = FALSE
  ]
  
  for (
    time in seq.int(
      p_lags + 1L,
      T_total
    )
  ) {
    
    y_time <- innovations[
      time,
    ]
    
    for (lag in seq_len(
      p_lags
    )) {
      
      y_time <-
        y_time +
        as.numeric(
          A_list[[lag]] %*%
            Y_full[
              time - lag,
            ]
        )
    }
    
    Y_full[
      time,
    ] <- y_time
  }
  
  rows_keep <- seq.int(
    burn_in + 1L,
    T_total
  )
  
  Y <- Y_full[
    rows_keep,
    ,
    drop = FALSE
  ]
  
  Y_list <- matrix_to_Y_list(
    Y =
      Y,
    n_units =
      n_units,
    m =
      m
  )
  
  truth <- make_truth_from_A(
    A_list =
      A_list,
    n_units =
      n_units,
    m =
      m,
    truth_eps =
      truth_eps
  )
  
  c(
    list(
      Y =
        Y,
      
      Y_list =
        Y_list,
      
      A_true =
        A_list,
      
      stability =
        stability,
      
      simulation =
        list(
          T_obs =
            T_obs,
          
          burn_in =
            burn_in,
          
          seed =
            seed
        )
    ),
    
    truth
  )
}