# Builds nested benchmark datasets from one frozen DGP.
#
# For each seed, simulates one path at the largest requested sample size and
# returns smaller sample sizes as prefixes with the same frozen network truth.


# DGP loading -----

parse_benchmark_dgp_id <- function(dgp_id) {
  match <- regexec("^(.*)_n([0-9]+)_m([0-9]+)_p([0-9]+)$", dgp_id)
  parts <- regmatches(dgp_id, match)[[1L]]

  if (length(parts) != 5L) {
    stop("Unable to parse DGP ID: ", dgp_id)
  }

  list(
    n_units = as.integer(parts[3L]),
    m = as.integer(parts[4L]),
    p_lags = as.integer(parts[5L])
  )
}


load_benchmark_dgp <- function(dgp_id, dgp_root = file.path("data", "dgp_bank")) {
  if (!is.character(dgp_id) || length(dgp_id) != 1L || !nzchar(dgp_id)) {
    stop("dgp_id must be one non-empty string")
  }

  dimensions <- parse_benchmark_dgp_id(dgp_id)
  A_path <- file.path(dgp_root, dgp_id, "A_list.rds")

  if (!file.exists(A_path)) {
    stop("Missing DGP coefficient file: ", A_path)
  }

  A_list <- readRDS(A_path)
  validate_A_list(A_list)

  N <- dimensions$n_units * dimensions$m
  expected_dim <- c(N, N)

  if (length(A_list) != dimensions$p_lags ||
      any(vapply(A_list, function(A) !identical(dim(A), expected_dim), logical(1L)))) {
    stop("DGP dimensions do not match its ID: ", dgp_id)
  }

  list(
    dgp_id = dgp_id,
    n_units = dimensions$n_units,
    m = dimensions$m,
    N = N,
    p_lags = dimensions$p_lags,
    A_list = A_list
  )
}


# Frozen truth -----

make_benchmark_truth <- function(A_list, n_units, m, truth_eps = 1e-12) {
  s_lag <- A_list_to_s_lag(A_list, n_units = n_units, m = m)
  G_lag <- s_lag > truth_eps

  list(
    G_lag = G_lag,
    G_unit = G_lag_to_G_unit(G_lag),
    s_lag = s_lag,
    s_unit_max = s_lag_to_unit(s_lag, method = "max"),
    s_unit_rms = s_lag_to_unit(s_lag, method = "rms")
  )
}


# Nested benchmark records -----

make_benchmark_records <- function(
    dgp_id,
    T_grid,
    seeds,
    burn_in = 300L,
    sigma = 1,
    dgp_root = file.path("data", "dgp_bank"),
    truth_eps = 1e-12) {

  T_grid <- sort(unique(as.integer(T_grid)))
  seeds <- sort(unique(as.integer(seeds)))
  burn_in <- as.integer(burn_in)

  if (length(T_grid) == 0L || anyNA(T_grid) || any(T_grid < 1L)) {
    stop("T_grid must contain positive integers")
  }

  if (length(seeds) == 0L || anyNA(seeds)) {
    stop("seeds must contain valid integers")
  }

  if (length(burn_in) != 1L || is.na(burn_in) || burn_in < 0L) {
    stop("burn_in must be a non-negative integer")
  }

  if (length(sigma) != 1L || !is.finite(sigma) || sigma <= 0) {
    stop("sigma must be a positive scalar")
  }

  if (length(truth_eps) != 1L || !is.finite(truth_eps) || truth_eps < 0) {
    stop("truth_eps must be a non-negative scalar")
  }

  dgp <- load_benchmark_dgp(dgp_id, dgp_root = dgp_root)
  truth <- make_benchmark_truth(dgp$A_list, dgp$n_units, dgp$m, truth_eps = truth_eps)

  T_max <- max(T_grid)
  Sigma <- diag(sigma^2, dgp$N)
  records <- vector("list", length(seeds) * length(T_grid))
  record_index <- 0L

  for (seed in seeds) {
    simulation <- simulate_var_from_A(
      A_list = dgp$A_list,
      Sigma = Sigma,
      T_obs = T_max,
      burn_in = burn_in,
      seed = seed,
      n_units = dgp$n_units,
      m = dgp$m,
      truth_eps = truth_eps
    )

    for (T_obs in T_grid) {
      record_index <- record_index + 1L
      record_id <- paste0(dgp_id, "_seed_", seed, "_T_", T_obs)

      Y_list <- lapply(
        simulation$Y_list,
        function(Y) Y[seq_len(T_obs), , drop = FALSE]
      )

      records[[record_index]] <- list(
        id = record_id,
        design = list(dgp_id = dgp_id, seed = seed, T_obs = T_obs),
        dimensions = list(n_units = dgp$n_units, m = dgp$m, N = dgp$N, p_lags = dgp$p_lags),
        data = list(Y_list = Y_list),
        truth = truth
      )
    }
  }

  names(records) <- vapply(records, `[[`, character(1L), "id")
  records
}
