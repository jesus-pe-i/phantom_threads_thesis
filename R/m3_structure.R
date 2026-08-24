# Structure, prior, and initial state for the M3 block BVAR.
#
# Builds block/group maps on the common VAR design, calibrates the
# global shrinkage prior, and constructs the initial Gibbs state.


# Group maps -----

normalize_m3_group_map <- function(
    map,
    k,
    N,
    self_block,
    anchor_self = FALSE) {
  
  map <- as.matrix(map)
  
  if (
    !identical(
      dim(map),
      c(k, N)
    ) ||
    any(!is.finite(map)) ||
    any(map != round(map))
  ) {
    stop(
      "M3 group maps must be finite integer k x N matrices"
    )
  }
  
  storage.mode(map) <- "integer"
  
  if (
    anchor_self &&
    length(
      unique(
        map[self_block]
      )
    ) != 1L
  ) {
    stop(
      "All self blocks must belong to one reference c group"
    )
  }
  
  levels <- sort(
    unique(
      as.vector(map)
    )
  )
  
  map <- matrix(
    match(
      as.vector(map),
      levels
    ),
    nrow = k,
    ncol = N
  )
  
  storage.mode(map) <- "integer"
  
  if (anchor_self) {
    
    self_group <- unique(
      map[self_block]
    )
    
    if (self_group != 1L) {
      
      map[map == 1L] <- -1L
      map[map == self_group] <- 1L
      map[map == -1L] <- self_group
    }
  }
  
  list(
    id = map,
    n = as.integer(
      length(
        unique(
          as.vector(map)
        )
      )
    )
  )
}


# Structure -----

make_m3_structure <- function(
    Y_list,
    p_lags,
    gc_id = NULL,
    gq_id = NULL,
    standardize = TRUE) {
  
  prepared <- prepare_var_data(
    Y_list = Y_list,
    p_lags = p_lags,
    standardize = standardize
  )
  
  n_units <- prepared$n_units
  m <- prepared$m
  N <- prepared$n_series
  k <- prepared$n_coef
  p_lags <- prepared$p_lags
  
  
  # Coefficient indices -----
  
  predictor <- rep(
    seq_len(N),
    times = p_lags
  )
  
  lag_vector <- rep(
    seq_len(p_lags),
    each = N
  )
  
  response <- seq_len(N)
  
  predictor_unit <- (
    (predictor - 1L) %/% m
  ) + 1L
  
  predictor_variable <- (
    (predictor - 1L) %% m
  ) + 1L
  
  response_unit <- (
    (response - 1L) %/% m
  ) + 1L
  
  response_variable <- (
    (response - 1L) %% m
  ) + 1L
  
  lag_id <- matrix(
    rep(
      lag_vector,
      times = N
    ),
    nrow = k,
    ncol = N
  )
  
  pred_unit <- matrix(
    rep(
      predictor_unit,
      times = N
    ),
    nrow = k,
    ncol = N
  )
  
  pred_var <- matrix(
    rep(
      predictor_variable,
      times = N
    ),
    nrow = k,
    ncol = N
  )
  
  resp_unit <- matrix(
    rep(
      response_unit,
      each = k
    ),
    nrow = k,
    ncol = N
  )
  
  resp_var <- matrix(
    rep(
      response_variable,
      each = k
    ),
    nrow = k,
    ncol = N
  )
  
  edge_id <- matrix(
    as.integer(
      (resp_unit - 1L) *
        n_units +
        pred_unit
    ),
    nrow = k,
    ncol = N
  )
  
  block_id <- matrix(
    as.integer(
      (lag_id - 1L) *
        n_units^2 +
        edge_id
    ),
    nrow = k,
    ncol = N
  )
  
  same_var <- pred_var ==
    resp_var
  
  self_block <- pred_unit ==
    resp_unit
  
  storage.mode(lag_id) <- "integer"
  storage.mode(block_id) <- "integer"
  storage.mode(same_var) <- "integer"
  
  
  # Shrinkage groups -----
  
  if (is.null(gc_id)) {
    
    gc_id <- ifelse(
      self_block,
      1L,
      2L
    )
    
    c_group <- "self_foreign"
    
  } else {
    
    c_group <- "custom"
  }
  
  if (is.null(gq_id)) {
    
    gq_id <- ifelse(
      self_block,
      1L,
      2L
    )
    
    q_group <- "self_foreign"
    
  } else {
    
    q_group <- "custom"
  }
  
  c_map <- normalize_m3_group_map(
    map = gc_id,
    k = k,
    N = N,
    self_block = self_block,
    anchor_self = TRUE
  )
  
  q_map <- normalize_m3_group_map(
    map = gq_id,
    k = k,
    N = N,
    self_block = self_block,
    anchor_self = FALSE
  )
  
  n_blocks <- as.integer(
    n_units^2 * p_lags
  )
  
  for (block in seq_len(n_blocks)) {
    
    location <- block_id ==
      block
    
    if (
      length(
        unique(
          c_map$id[location]
        )
      ) != 1L ||
      length(
        unique(
          q_map$id[location]
        )
      ) != 1L
    ) {
      stop(
        "Each M3 block must belong to one c group and one q group"
      )
    }
  }
  
  
  # Output -----
  
  list(
    data = list(
      Y = prepared$Y,
      X = prepared$X,
      T_p = prepared$T_eff,
      K = n_units,
      m = m,
      N = N,
      p_lags = p_lags,
      k = k
    ),
    
    maps = list(
      lag_id = lag_id,
      gc_id = c_map$id,
      gq_id = q_map$id,
      block_id = block_id,
      same_var = same_var,
      phi2 = rep(
        1,
        p_lags
      ),
      n_c = c_map$n,
      n_q = q_map$n,
      n_blocks = n_blocks
    ),
    
    dimensions = list(
      n_units = n_units,
      m = m,
      N = N,
      p_lags = p_lags,
      k = k,
      T_full = prepared$T_full,
      T_p = prepared$T_eff
    ),
    
    preprocessing = list(
      standardize = standardize,
      centering = FALSE,
      scale = prepared$scale,
      series_names =
        prepared$series_names
    ),
    
    grouping = list(
      c_group = c_group,
      q_group = q_group
    )
  )
}


# Prior -----

calibrate_m3_tau_scale <- function(
    n_blocks,
    T_p,
    p0 = NULL) {
  
  if (is.null(p0)) {
    p0 <- n_blocks / 2
  }
  
  if (
    length(p0) != 1L ||
    !is.finite(p0) ||
    p0 <= 0 ||
    p0 >= n_blocks
  ) {
    stop(
      "p0 must lie strictly between zero and n_blocks"
    )
  }
  
  as.numeric(
    (
      p0 /
        (n_blocks - p0)
    ) /
      sqrt(T_p)
  )
}


make_m3_prior <- function(
    structure,
    tau_scale = NULL,
    p0 = NULL,
    q_grid = seq(
      -log(10),
      log(10),
      length.out = 17L
    ),
    q_prob = NULL,
    sigma_a = 3,
    sigma_b = 2,
    tau_df = 10,
    c_df = 5,
    lambda_df = 3,
    c_scale = 1,
    lambda_scale = 1) {
  
  positive_values <- c(
    sigma_a = sigma_a,
    sigma_b = sigma_b,
    tau_df = tau_df,
    c_df = c_df,
    lambda_df = lambda_df,
    c_scale = c_scale,
    lambda_scale = lambda_scale
  )
  
  if (
    any(!is.finite(positive_values)) ||
    any(positive_values <= 0)
  ) {
    stop(
      "M3 prior scales and degrees of freedom must be positive"
    )
  }
  
  n_blocks <- structure$maps$n_blocks
  T_p <- structure$data$T_p
  
  if (is.null(p0)) {
    p0 <- n_blocks / 2
  }
  
  if (
    length(p0) != 1L ||
    !is.finite(p0) ||
    p0 <= 0 ||
    p0 >= n_blocks
  ) {
    stop(
      "p0 must lie strictly between zero and n_blocks"
    )
  }
  
  if (is.null(tau_scale)) {
    
    tau_scale <- calibrate_m3_tau_scale(
      n_blocks = n_blocks,
      T_p = T_p,
      p0 = p0
    )
  }
  
  if (
    length(tau_scale) != 1L ||
    !is.finite(tau_scale) ||
    tau_scale <= 0
  ) {
    stop(
      "tau_scale must be finite and positive"
    )
  }
  
  q_grid <- as.numeric(
    q_grid
  )
  
  if (
    length(q_grid) < 2L ||
    any(!is.finite(q_grid)) ||
    is.unsorted(
      q_grid,
      strictly = TRUE
    )
  ) {
    stop(
      "q_grid must be strictly increasing and finite"
    )
  }
  
  if (is.null(q_prob)) {
    
    q_prob <- rep(
      1 / length(q_grid),
      length(q_grid)
    )
  }
  
  q_prob <- as.numeric(
    q_prob
  )
  
  if (
    length(q_prob) !=
    length(q_grid) ||
    any(!is.finite(q_prob)) ||
    any(q_prob < 0) ||
    sum(q_prob) <= 0
  ) {
    stop(
      "q_prob must contain one non-negative weight per q-grid value"
    )
  }
  
  q_prob <- q_prob /
    sum(q_prob)
  
  list(
    sigma_a = as.numeric(sigma_a),
    sigma_b = as.numeric(sigma_b),
    tau_df = as.numeric(tau_df),
    tau_scale = as.numeric(tau_scale),
    c_df = as.numeric(c_df),
    c_scale = as.numeric(c_scale),
    lambda_df = as.numeric(lambda_df),
    lambda_scale =
      as.numeric(lambda_scale),
    p0 = as.numeric(p0),
    q_grid = q_grid,
    q_prob = q_prob
  )
}


# Initial state -----

make_m3_state <- function(
    structure,
    prior,
    beta = NULL,
    sigma2 = NULL,
    tau2 = 1,
    psi_tau = 1,
    c2 = NULL,
    psi_c = NULL,
    lambda2 = NULL,
    xi = NULL,
    q = 0) {
  
  k <- structure$data$k
  N <- structure$data$N
  n_c <- structure$maps$n_c
  n_q <- structure$maps$n_q
  n_blocks <- structure$maps$n_blocks
  
  if (is.null(beta)) {
    
    beta <- matrix(
      0,
      nrow = k,
      ncol = N
    )
  }
  
  if (is.null(sigma2)) {
    
    sigma2 <- rep(
      1,
      N
    )
  }
  
  if (is.null(c2)) {
    
    c2 <- rep(
      1,
      n_c
    )
  }
  
  if (is.null(psi_c)) {
    
    psi_c <- rep(
      1,
      n_c
    )
  }
  
  if (is.null(lambda2)) {
    
    lambda2 <- rep(
      1,
      n_blocks
    )
  }
  
  if (is.null(xi)) {
    
    xi <- rep(
      1,
      n_blocks
    )
  }
  
  q <- rep_len(
    as.numeric(q),
    n_q
  )
  
  if (any(!is.finite(q))) {
    stop(
      "Initial q values must be finite"
    )
  }
  
  q_index <- vapply(
    q,
    function(value) {
      which.min(
        abs(
          prior$q_grid -
            value
        )
      )
    },
    integer(1L)
  )
  
  state <- list(
    beta = unname(
      as.matrix(beta)
    ),
    sigma2 = as.numeric(sigma2),
    tau2 = as.numeric(tau2),
    psi_tau =
      as.numeric(psi_tau),
    c2 = as.numeric(c2),
    psi_c = as.numeric(psi_c),
    lambda2 =
      as.numeric(lambda2),
    xi = as.numeric(xi),
    q_index =
      as.integer(q_index)
  )
  
  state$c2[1L] <- 1
  state$psi_c[1L] <- 1
  
  state
}