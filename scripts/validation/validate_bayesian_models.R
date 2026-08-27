# Permanent integration validation for the three Bayesian BVAR engines.
#
# Checks model structure, frozen prior defaults, principal execution paths,
# posterior outputs, selected alternate paths, and fixed-seed reproducibility.


rm(
  list = ls()
)


# Helpers -----

checks <- logical(0L)


check <- function(
    label,
    condition) {
  
  checks[label] <<-
    isTRUE(
      condition
    )
  
  cat(
    sprintf(
      "%-44s %s\n",
      label,
      if (checks[label]) {
        "PASS"
      } else {
        "FAIL"
      }
    )
  )
  
  invisible(
    checks[label]
  )
}


same <- function(
    x,
    y,
    tolerance = 1e-12) {
  
  isTRUE(
    all.equal(
      x,
      y,
      tolerance = tolerance,
      check.attributes = FALSE
    )
  )
}


finite <- function(
    x) {
  
  all(
    is.finite(
      as.numeric(
        x
      )
    )
  )
}


positive_finite <- function(
    x) {
  
  finite(
    x
  ) &&
    all(
      as.numeric(
        x
      ) > 0
    )
}


check_fit_output <- function(
    fit,
    model,
    chains,
    draws,
    n_units,
    m,
    p_lags) {
  
  N <- n_units * m
  k <- N * p_lags
  
  s_check <- A_list_to_s_lag(
    A_list = fit$A_hat_lag,
    n_units = n_units,
    m = m
  )
  
  all(
    identical(
      fit$model,
      model
    ),
    
    identical(
      dim(
        fit$beta_hat
      ),
      c(
        k,
        N
      )
    ),
    
    identical(
      dim(
        fit$beta_standardized
      ),
      c(
        k,
        N
      )
    ),
    
    finite(
      fit$beta_hat
    ),
    
    finite(
      fit$beta_standardized
    ),
    
    length(
      fit$sigma2_hat
    ) == N,
    
    length(
      fit$sigma2_standardized
    ) == N,
    
    positive_finite(
      fit$sigma2_hat
    ),
    
    positive_finite(
      fit$sigma2_standardized
    ),
    
    length(
      fit$A_hat_lag
    ) == p_lags,
    
    all(
      vapply(
        fit$A_hat_lag,
        function(A) {
          
          identical(
            dim(A),
            c(
              N,
              N
            )
          ) &&
            finite(
              A
            )
        },
        logical(1L)
      )
    ),
    
    identical(
      dim(
        fit$s_hat_lag
      ),
      c(
        n_units,
        n_units,
        p_lags
      )
    ),
    
    identical(
      dim(
        fit$s_hat_unit_max
      ),
      c(
        n_units,
        n_units
      )
    ),
    
    identical(
      dim(
        fit$s_hat_unit_rms
      ),
      c(
        n_units,
        n_units
      )
    ),
    
    finite(
      fit$s_hat_lag
    ),
    
    all(
      fit$s_hat_lag >= 0
    ),
    
    same(
      fit$s_hat_lag,
      s_check
    ),
    
    same(
      fit$s_hat_unit_max,
      s_lag_to_unit(
        fit$s_hat_lag,
        method = "max"
      )
    ),
    
    same(
      fit$s_hat_unit_rms,
      s_lag_to_unit(
        fit$s_hat_lag,
        method = "rms"
      )
    ),
    
    same(
      fit$score_lag,
      fit$s_hat_lag
    ),
    
    is.logical(
      fit$selected_lag
    ),
    
    identical(
      dim(
        fit$selected_lag
      ),
      dim(
        fit$s_hat_lag
      )
    ),
    
    identical(
      as.integer(
        fit$retained_per_chain
      ),
      rep(
        as.integer(draws),
        chains
      )
    ),
    
    fit$retained_total ==
      chains * draws,
    
    identical(
      as.integer(
        fit$chain_seeds
      ),
      as.integer(
        fit$chain_seeds[1L] +
          seq_len(chains) -
          1L
      )
    ),
    
    fit$dimensions$n_units ==
      n_units,
    
    fit$dimensions$m ==
      m,
    
    fit$dimensions$N ==
      N,
    
    fit$dimensions$p_lags ==
      p_lags,
    
    fit$dimensions$k ==
      k
  )
}


# Source files -----

source_paths <- c(
  "R/var_data.R",
  "R/var_coefficients.R",
  "R/network_blocks.R",
  "R/simulation.R",
  "R/half_t_structure.R",
  "R/half_t_sampler.R",
  "R/half_t_fit.R",
  "R/m3_structure.R",
  "R/m3_sampler.R",
  "R/m3_fit.R",
  "R/gigg_structure.R",
  "R/gigg_sampler.R",
  "R/gigg_fit.R"
)


cpp_paths <- c(
  "src/half_t_bvar_core.cpp",
  "src/m3_bvar_core.cpp",
  "src/gigg_bvar_core.cpp"
)


required_paths <- c(
  source_paths,
  cpp_paths
)


missing_paths <- required_paths[
  !file.exists(
    required_paths
  )
]


if (length(missing_paths) > 0L) {
  
  stop(
    paste0(
      "Run from the project root. Missing: ",
      paste(
        missing_paths,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}


for (path in source_paths) {
  source(
    path
  )
}


# Native engines -----

if (
  !requireNamespace(
    "Rcpp",
    quietly = TRUE
  ) ||
  !requireNamespace(
    "RcppArmadillo",
    quietly = TRUE
  )
) {
  
  stop(
    "Rcpp and RcppArmadillo are required",
    call. = FALSE
  )
}


cat(
  "\nBAYESIAN MODEL VALIDATION\n\n"
)


cat(
  "Compiling native engines...\n\n"
)


load_half_t_cpp(
  rebuild = FALSE
)

load_m3_cpp(
  rebuild = FALSE
)

load_gigg_cpp(
  rebuild = FALSE
)


check(
  "Native engines loaded",
  exists(
    "half_t_bvar_chain_cpp",
    mode = "function"
  ) &&
    exists(
      "m3_bvar_chain_cpp",
      mode = "function"
    ) &&
    exists(
      "gigg_bvar_chain_cpp",
      mode = "function"
    )
)


# Synthetic VAR -----

n_units <- 2L
m <- 2L
p_lags <- 2L

N <- n_units * m
k <- N * p_lags

T_obs <- 80L


A1 <- matrix(
  0,
  nrow = N,
  ncol = N
)

diag(A1) <- c(
  0.32,
  0.28,
  0.30,
  0.26
)

A1[1L, 3L] <- 0.12
A1[2L, 4L] <- -0.10
A1[3L, 1L] <- 0.08
A1[4L, 2L] <- 0.10


A2 <- matrix(
  0,
  nrow = N,
  ncol = N
)

diag(A2) <- c(
  0.08,
  0.06,
  0.07,
  0.05
)

A2[1L, 4L] <- 0.05
A2[3L, 2L] <- -0.04


A_list <- list(
  A1,
  A2
)


stability <- check_A_radius(
  A_list
)


simulation <- simulate_var_from_A(
  A_list = A_list,
  Sigma = make_sigma_homoscedastic(
    N = N,
    sigma = 0.7
  ),
  T_obs = T_obs,
  burn_in = 120L,
  seed = 1201L,
  n_units = n_units,
  m = m
)


Y_list <- simulation$Y_list


check(
  "Synthetic VAR",
  stability$radius < 1 &&
    length(Y_list) ==
    n_units &&
    all(
      vapply(
        Y_list,
        function(Y) {
          
          identical(
            dim(Y),
            c(
              T_obs,
              m
            )
          ) &&
            finite(
              Y
            )
        },
        logical(1L)
      )
    )
)


prepared <- prepare_var_data(
  Y_list = Y_list,
  p_lags = p_lags,
  standardize = TRUE
)


check(
  "Common VAR structure",
  prepared$n_units ==
    n_units &&
    prepared$m ==
    m &&
    prepared$n_series ==
    N &&
    prepared$n_coef ==
    k &&
    prepared$p_lags ==
    p_lags &&
    prepared$T_eff ==
    T_obs - p_lags
)


# Half-t -----

## Structure -----

half_structure_all <- make_half_t_structure(
  Y_list = Y_list,
  p_lags = p_lags,
  global_grouping = "all"
)


half_structure_sd <- make_half_t_structure(
  Y_list = Y_list,
  p_lags = p_lags,
  global_grouping = "self_diagonal"
)


expected_lag_id <- rep(
  seq_len(p_lags),
  each = N
)


expected_phi2 <-
  1 /
  expected_lag_id^2


check(
  "Half-t structure",
  identical(
    half_structure_sd$maps$lag_id,
    as.integer(
      expected_lag_id
    )
  ) &&
    same(
      half_structure_sd$maps$phi2,
      expected_phi2
    ) &&
    half_structure_all$maps$n_tau ==
    1L &&
    half_structure_sd$maps$n_tau ==
    2L
)


## Prior and state -----

half_prior <- make_half_t_prior()


check(
  "Half-t canonical prior",
  half_prior$tau_df ==
    10 &&
    half_prior$tau_scale ==
    1 &&
    half_prior$lambda_df ==
    3 &&
    half_prior$lambda_scale ==
    1 &&
    half_prior$sigma_a ==
    3 &&
    half_prior$sigma_b ==
    2
)


half_state <- make_half_t_state(
  half_structure_sd
)


check(
  "Half-t initial state",
  identical(
    dim(
      half_state$beta
    ),
    c(
      k,
      N
    )
  ) &&
    identical(
      dim(
        half_state$lambda2
      ),
      c(
        k,
        N
      )
    ) &&
    length(
      half_state$tau2
    ) ==
    2L &&
    length(
      half_state$xi
    ) ==
    2L &&
    positive_finite(
      half_state$sigma2
    ) &&
    positive_finite(
      half_state$tau2
    ) &&
    positive_finite(
      half_state$xi
    ) &&
    positive_finite(
      half_state$lambda2
    )
)


## Main fit -----

half_fit <- fit_half_t(
  Y_list = Y_list,
  p_lags = p_lags,
  global_grouping = "self_diagonal",
  chains = 2L,
  burnin = 40L,
  draws = 80L,
  thin = 1L,
  seed = 2101L,
  keep_chain_results = TRUE
)


check(
  "Half-t fit output",
  check_fit_output(
    fit = half_fit,
    model = "half_t",
    chains = 2L,
    draws = 80L,
    n_units = n_units,
    m = m,
    p_lags = p_lags
  )
)


check(
  "Half-t posterior scales",
  positive_finite(
    half_fit$tau2_mean
  ) &&
    positive_finite(
      half_fit$xi_mean
    )
)


check(
  "Half-t default path",
  identical(
    half_fit$global_grouping,
    "self_diagonal"
  ) &&
    isTRUE(
      half_fit$control$use_asis
    ) &&
    half_fit$prior$tau_df ==
    10 &&
    half_fit$prior$lambda_df ==
    3
)


## Alternate fit -----

half_fit_alt <- fit_half_t(
  Y_list = Y_list,
  p_lags = p_lags,
  global_grouping = "all",
  chains = 1L,
  burnin = 20L,
  draws = 30L,
  thin = 1L,
  seed = 2201L,
  beta_algorithm = "bhattacharya",
  tau_df = 8,
  lambda_df = 4
)


check(
  "Half-t alternate path",
  identical(
    half_fit_alt$global_grouping,
    "all"
  ) &&
    !isTRUE(
      half_fit_alt$control$use_asis
    ) &&
    identical(
      half_fit_alt$beta_algorithm_resolved,
      "bhattacharya"
    ) &&
    half_fit_alt$prior$tau_df ==
    8 &&
    half_fit_alt$prior$lambda_df ==
    4
)


## Reproducibility -----

half_fit_repeat <- fit_half_t(
  Y_list = Y_list,
  p_lags = p_lags,
  global_grouping = "self_diagonal",
  chains = 2L,
  burnin = 40L,
  draws = 80L,
  thin = 1L,
  seed = 2101L
)


check(
  "Half-t reproducibility",
  same(
    half_fit$beta_hat,
    half_fit_repeat$beta_hat
  ) &&
    same(
      half_fit$s_hat_lag,
      half_fit_repeat$s_hat_lag
    ) &&
    same(
      half_fit$tau2_mean,
      half_fit_repeat$tau2_mean
    ) &&
    same(
      half_fit$xi_mean,
      half_fit_repeat$xi_mean
    )
)


# M3 -----

## Structure -----

m3_structure <- make_m3_structure(
  Y_list = Y_list,
  p_lags = p_lags
)


expected_blocks <-
  n_units^2 *
  p_lags


m3_block_counts <- tabulate(
  as.vector(
    m3_structure$maps$block_id
  ),
  nbins = expected_blocks
)


check(
  "M3 structure",
  identical(
    dim(
      m3_structure$maps$lag_id
    ),
    c(
      k,
      N
    )
  ) &&
    identical(
      dim(
        m3_structure$maps$gc_id
      ),
      c(
        k,
        N
      )
    ) &&
    identical(
      dim(
        m3_structure$maps$gq_id
      ),
      c(
        k,
        N
      )
    ) &&
    identical(
      dim(
        m3_structure$maps$block_id
      ),
      c(
        k,
        N
      )
    ) &&
    m3_structure$maps$n_blocks ==
    expected_blocks &&
    all(
      m3_block_counts ==
        m^2
    ) &&
    same(
      m3_structure$maps$phi2,
      rep(
        1,
        p_lags
      )
    )
)


check(
  "M3 default grouping",
  identical(
    m3_structure$grouping$c_group,
    "self_foreign"
  ) &&
    identical(
      m3_structure$grouping$q_group,
      "self_foreign"
    ) &&
    m3_structure$maps$n_c ==
    2L &&
    m3_structure$maps$n_q ==
    2L
)


## Prior and state -----

m3_prior <- make_m3_prior(
  structure = m3_structure
)


expected_m3_tau_scale <-
  1 /
  sqrt(
    m3_structure$data$T_p
  )


expected_q_grid <- seq(
  -log(1000),
  log(1000),
  length.out = 49L
)


check(
  "M3 canonical prior",
  m3_prior$tau_df ==
    10 &&
    m3_prior$c_df ==
    5 &&
    m3_prior$lambda_df ==
    3 &&
    m3_prior$c_scale ==
    1 &&
    m3_prior$lambda_scale ==
    1 &&
    m3_prior$sigma_a ==
    3 &&
    m3_prior$sigma_b ==
    2 &&
    m3_prior$p0 ==
    expected_blocks / 2 &&
    same(
      m3_prior$tau_scale,
      expected_m3_tau_scale
    ) &&
    same(
      m3_prior$q_grid,
      expected_q_grid
    ) &&
    same(
      m3_prior$q_prob,
      rep(
        1 / 49,
        49L
      )
    )
)


m3_state <- make_m3_state(
  structure = m3_structure,
  prior = m3_prior
)


check(
  "M3 initial state",
  identical(
    dim(
      m3_state$beta
    ),
    c(
      k,
      N
    )
  ) &&
    length(
      m3_state$sigma2
    ) ==
    N &&
    length(
      m3_state$tau2
    ) ==
    1L &&
    length(
      m3_state$psi_tau
    ) ==
    1L &&
    length(
      m3_state$c2
    ) ==
    m3_structure$maps$n_c &&
    length(
      m3_state$psi_c
    ) ==
    m3_structure$maps$n_c &&
    length(
      m3_state$lambda2
    ) ==
    expected_blocks &&
    length(
      m3_state$xi
    ) ==
    expected_blocks &&
    length(
      m3_state$q_index
    ) ==
    m3_structure$maps$n_q
)


## Control -----

m3_control <- make_m3_control(
  burnin = 40L,
  draws = 80L
)



## Main fit -----

m3_fit <- fit_m3(
  Y_list = Y_list,
  p_lags = p_lags,
  chains = 2L,
  burnin = 40L,
  draws = 80L,
  thin = 1L,
  seed = 3101L,
  keep_chain_results = TRUE
)


check(
  "M3 fit output",
  check_fit_output(
    fit = m3_fit,
    model = "m3",
    chains = 2L,
    draws = 80L,
    n_units = n_units,
    m = m,
    p_lags = p_lags
  )
)


check(
  "M3 posterior scales",
  positive_finite(
    m3_fit$tau2_mean
  ) &&
    positive_finite(
      m3_fit$c2_mean
    ) &&
    finite(
      m3_fit$q_mean
    )
)


check(
  "M3 default path",
  identical(
    m3_fit$q_update,
    "gibbs_end"
  ) &&
    identical(
      m3_fit$control$q_update,
      "gibbs_end"
    ) &&
    identical(
      m3_fit$control$scale_transport,
      "none"
    ) &&
    isTRUE(
      m3_fit$control$use_c_asis
    ) &&
    m3_fit$control$c_asis_every ==
    1L
)


m3_chain_contract <- vapply(
  m3_fit$chain_results,
  function(chain) {
    
    all(
      identical(
        chain$q_update,
        "gibbs_end"
      ),
      
      identical(
        dim(
          chain$draws$c2
        ),
        c(
          chain$retained,
          m3_structure$maps$n_c
        )
      ),
      
      identical(
        dim(
          chain$draws$q
        ),
        c(
          chain$retained,
          m3_structure$maps$n_q
        )
      ),
      
      identical(
        dim(
          chain$draws$q_index
        ),
        c(
          chain$retained,
          m3_structure$maps$n_q
        )
      ),
      
      positive_finite(
        chain$draws$c2
      ),
      
      finite(
        chain$draws$q
      ),
      
      isTRUE(
        chain$c_asis$enabled
      ),
      
      chain$c_asis$every ==
        1L,
      
      length(
        chain$c_asis$updates
      ) ==
        m3_structure$maps$n_c,
      
      length(
        chain$c_asis$moves
      ) ==
        m3_structure$maps$n_c,
      
      chain$c_asis$updates[1L] ==
        0L,
      
      chain$c_asis$moves[1L] ==
        0L,
      
      all(
        chain$c_asis$updates[-1L] >
          0L
      ),
      
      all(
        chain$c_asis$moves[-1L] >=
          0L
      )
    )
  },
  logical(1L)
)


check(
  "M3 chain diagnostics",
  all(
    m3_chain_contract
  )
)


## Alternate fit -----

m3_fit_alt <- fit_m3(
  Y_list = Y_list,
  p_lags = p_lags,
  chains = 1L,
  burnin = 20L,
  draws = 30L,
  thin = 1L,
  seed = 3201L,
  q_update = "gibbs_end",
  use_c_asis = FALSE,
  tau_df = 8,
  c_df = 6,
  lambda_df = 4,
  keep_chain_results = TRUE
)


m3_alt_chain <-
  m3_fit_alt$chain_results[[1L]]


check(
  "M3 alternate path",
  identical(
    m3_fit_alt$q_update,
    "gibbs_end"
  ) &&
    identical(
      m3_alt_chain$q_update,
      "gibbs_end"
    ) &&
    !isTRUE(
      m3_fit_alt$control$use_c_asis
    ) &&
    !isTRUE(
      m3_alt_chain$c_asis$enabled
    ) &&
    all(
      m3_alt_chain$c_asis$updates ==
        0L
    ) &&
    all(
      m3_alt_chain$c_asis$moves ==
        0L
    ) &&
    m3_fit_alt$prior$tau_df ==
    8 &&
    m3_fit_alt$prior$c_df ==
    6 &&
    m3_fit_alt$prior$lambda_df ==
    4
)

## Reproducibility -----

m3_fit_repeat <- fit_m3(
  Y_list = Y_list,
  p_lags = p_lags,
  chains = 2L,
  burnin = 40L,
  draws = 80L,
  thin = 1L,
  seed = 3101L
)


check(
  "M3 reproducibility",
  same(
    m3_fit$beta_hat,
    m3_fit_repeat$beta_hat
  ) &&
    same(
      m3_fit$s_hat_lag,
      m3_fit_repeat$s_hat_lag
    ) &&
    same(
      m3_fit$tau2_mean,
      m3_fit_repeat$tau2_mean
    ) &&
    same(
      m3_fit$c2_mean,
      m3_fit_repeat$c2_mean
    ) &&
    same(
      m3_fit$q_mean,
      m3_fit_repeat$q_mean
    )
)


# GIGG -----

## Structure -----

gigg_structure <- make_gigg_structure(
  Y_list = Y_list,
  p_lags = p_lags
)


gigg_block_counts <- tabulate(
  as.vector(
    gigg_structure$maps$block_id
  ),
  nbins = expected_blocks
)


check(
  "GIGG structure",
  gigg_structure$maps$n_blocks ==
    expected_blocks &&
    all(
      gigg_block_counts ==
        m^2
    ) &&
    isTRUE(
      gigg_structure$preprocessing$standardize
    ) &&
    !isTRUE(
      gigg_structure$preprocessing$centering
    )
)


## Prior and state -----

gigg_prior <- make_gigg_prior(
  gigg_structure
)


expected_gigg_tau_scale <-
  3 /
  sqrt(
    5 *
      gigg_structure$data$T_p
  )


check(
  "GIGG canonical prior",
  gigg_prior$sigma_a ==
    3 &&
    gigg_prior$sigma_b ==
    2 &&
    gigg_prior$tau_df ==
    5 &&
    same(
      gigg_prior$tau_scale,
      expected_gigg_tau_scale
    ) &&
    gigg_prior$gamma_shape ==
    0.5 &&
    gigg_prior$gamma_rate ==
    1 &&
    gigg_prior$lambda_shape ==
    2.5 &&
    gigg_prior$lambda_scale ==
    1
)


gigg_state <- make_gigg_state(
  gigg_structure
)


check(
  "GIGG initial state",
  identical(
    dim(
      gigg_state$beta
    ),
    c(
      k,
      N
    )
  ) &&
    length(
      gigg_state$sigma2
    ) ==
    N &&
    same(
      gigg_state$tau2,
      1 /
        gigg_structure$data$T_p
    ) &&
    length(
      gigg_state$psi_tau
    ) ==
    1L &&
    length(
      gigg_state$gamma2
    ) ==
    expected_blocks &&
    identical(
      dim(
        gigg_state$lambda2
      ),
      c(
        k,
        N
      )
    ) &&
    positive_finite(
      gigg_state$sigma2
    ) &&
    positive_finite(
      gigg_state$tau2
    ) &&
    positive_finite(
      gigg_state$psi_tau
    ) &&
    positive_finite(
      gigg_state$gamma2
    ) &&
    positive_finite(
      gigg_state$lambda2
    )
)


## Main fit -----

monitor_beta <- matrix(
  c(
    1L, 1L,
    k, N
  ),
  ncol = 2L,
  byrow = TRUE
)


monitor_gamma <- c(
  1L,
  expected_blocks
)


monitor_lambda <- matrix(
  c(
    1L, 1L,
    k, N
  ),
  ncol = 2L,
  byrow = TRUE
)


gigg_fit <- fit_gigg(
  Y_list = Y_list,
  p_lags = p_lags,
  chains = 2L,
  burnin = 40L,
  draws = 80L,
  thin = 1L,
  seed = 4101L,
  monitor_beta = monitor_beta,
  monitor_gamma = monitor_gamma,
  monitor_lambda = monitor_lambda,
  keep_chain_results = TRUE
)


check(
  "GIGG fit output",
  check_fit_output(
    fit = gigg_fit,
    model = "gigg",
    chains = 2L,
    draws = 80L,
    n_units = n_units,
    m = m,
    p_lags = p_lags
  )
)


stale_gigg_controls <- c(
  "use_global_tau",
  "numerical_epsilon",
  "invariant_tolerance",
  "gig_maximum_attempts"
)


check(
  "GIGG control",
  !any(
    stale_gigg_controls %in%
      names(
        gigg_fit$control
      )
  ) &&
    isTRUE(
      gigg_fit$control$use_asis
    ) &&
    gigg_fit$control$asis_every ==
    1L
)


check(
  "GIGG posterior scales",
  positive_finite(
    gigg_fit$tau2_mean
  ) &&
    positive_finite(
      gigg_fit$psi_tau_mean
    ) &&
    length(
      gigg_fit$gamma2_mean
    ) ==
    expected_blocks &&
    positive_finite(
      gigg_fit$gamma2_mean
    ) &&
    identical(
      dim(
        gigg_fit$lambda2_mean
      ),
      c(
        k,
        N
      )
    ) &&
    positive_finite(
      gigg_fit$lambda2_mean
    ) &&
    length(
      gigg_fit$tau2_gamma2_mean
    ) ==
    expected_blocks &&
    positive_finite(
      gigg_fit$tau2_gamma2_mean
    ) &&
    identical(
      dim(
        gigg_fit$effective_variance_mean
      ),
      c(
        k,
        N
      )
    ) &&
    positive_finite(
      gigg_fit$effective_variance_mean
    )
)


gigg_chain_contract <- vapply(
  gigg_fit$chain_results,
  function(chain) {
    
    all(
      chain$gig$draws >
        0,
      
      chain$gig$total_attempts >
        0,
      
      isTRUE(
        chain$asis$enabled
      ),
      
      chain$asis$
        block_coefficient$
        invariant_failures ==
        0,
      
      chain$asis$
        global_block$
        invariant_failures ==
        0,
      
      chain$asis$
        block_coefficient$
        maximum_relative_invariant_error <
        1e-10,
      
      chain$asis$
        global_block$
        maximum_relative_invariant_error <
        1e-10,
      
      identical(
        dim(
          chain$draws$beta
        ),
        c(
          chain$retained,
          2L
        )
      ),
      
      identical(
        dim(
          chain$draws$gamma2
        ),
        c(
          chain$retained,
          2L
        )
      ),
      
      identical(
        dim(
          chain$draws$lambda2
        ),
        c(
          chain$retained,
          2L
        )
      )
    )
  },
  logical(1L)
)


check(
  "GIGG chain diagnostics",
  all(
    gigg_chain_contract
  )
)


monitor_names <- c(
  colnames(
    gigg_fit$
      chain_results[[1L]]$
      draws$gamma2
  ),
  colnames(
    gigg_fit$
      chain_results[[1L]]$
      draws$tau2_gamma2
  ),
  colnames(
    gigg_fit$
      chain_results[[1L]]$
      draws$lambda2
  ),
  colnames(
    gigg_fit$
      chain_results[[1L]]$
      draws$effective_variance
  )
)


check(
  "GIGG monitor naming",
  !any(
    grepl(
      "^log_",
      monitor_names
    )
  )
)


## Alternate fit -----

gigg_fit_alt <- fit_gigg(
  Y_list = Y_list,
  p_lags = p_lags,
  chains = 1L,
  burnin = 20L,
  draws = 30L,
  thin = 1L,
  seed = 4201L,
  use_asis = FALSE,
  keep_chain_results = TRUE
)


check(
  "GIGG no-ASIS path",
  !isTRUE(
    gigg_fit_alt$control$use_asis
  ) &&
    !isTRUE(
      gigg_fit_alt$
        chain_results[[1L]]$
        asis$enabled
    ) &&
    positive_finite(
      gigg_fit_alt$tau2_mean
    ) &&
    positive_finite(
      gigg_fit_alt$gamma2_mean
    ) &&
    positive_finite(
      gigg_fit_alt$lambda2_mean
    )
)


## Reproducibility -----

gigg_fit_repeat <- fit_gigg(
  Y_list = Y_list,
  p_lags = p_lags,
  chains = 2L,
  burnin = 40L,
  draws = 80L,
  thin = 1L,
  seed = 4101L,
  monitor_beta = monitor_beta,
  monitor_gamma = monitor_gamma,
  monitor_lambda = monitor_lambda
)


check(
  "GIGG reproducibility",
  same(
    gigg_fit$beta_hat,
    gigg_fit_repeat$beta_hat
  ) &&
    same(
      gigg_fit$s_hat_lag,
      gigg_fit_repeat$s_hat_lag
    ) &&
    same(
      gigg_fit$tau2_mean,
      gigg_fit_repeat$tau2_mean
    ) &&
    same(
      gigg_fit$psi_tau_mean,
      gigg_fit_repeat$psi_tau_mean
    ) &&
    same(
      gigg_fit$gamma2_mean,
      gigg_fit_repeat$gamma2_mean
    ) &&
    same(
      gigg_fit$lambda2_mean,
      gigg_fit_repeat$lambda2_mean
    ) &&
    same(
      gigg_fit$tau2_gamma2_mean,
      gigg_fit_repeat$tau2_gamma2_mean
    ) &&
    same(
      gigg_fit$effective_variance_mean,
      gigg_fit_repeat$
        effective_variance_mean
    )
)


# Report -----

cat(
  "\nSUMMARY\n",
  "=======\n\n",
  sep = ""
)


print(
  data.frame(
    check = names(checks),
    passed = unname(
      checks
    )
  ),
  row.names = FALSE
)


cat(
  "\nOVERALL:",
  if (all(checks)) {
    "PASS\n"
  } else {
    "FAIL\n"
  }
)


if (!all(checks)) {
  stop(
    "Bayesian model validation failed",
    call. = FALSE
  )
}