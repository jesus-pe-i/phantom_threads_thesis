# Permanent integration validation for the three Bayesian BVAR engines.
#
# Checks structural contracts, frozen prior defaults, principal execution
# paths, posterior output integrity, selected alternate paths, and
# fixed-seed reproducibility on one small synthetic block VAR.


rm(
  list = ls()
)


# Helpers -----

check <- function(
    label,
    condition) {
  
  if (
    length(condition) != 1L ||
    is.na(condition) ||
    !isTRUE(condition)
  ) {
    stop(
      paste0(
        "FAILED: ",
        label
      ),
      call. = FALSE
    )
  }
  
  cat(
    sprintf(
      "%-42s PASSED\n",
      label
    )
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


positive_finite <- function(x) {
  
  all(
    is.finite(x)
  ) &&
    all(
      x > 0
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
  
  conditions <- c(
    identical(
      fit$model,
      model
    ),
    
    identical(
      dim(fit$beta_hat),
      c(k, N)
    ),
    
    identical(
      dim(fit$beta_standardized),
      c(k, N)
    ),
    
    length(fit$sigma2_hat) ==
      N,
    
    positive_finite(
      fit$sigma2_hat
    ),
    
    length(fit$A_hat_lag) ==
      p_lags,
    
    all(
      vapply(
        fit$A_hat_lag,
        function(A) {
          
          identical(
            dim(A),
            c(N, N)
          ) &&
            all(
              is.finite(A)
            )
        },
        logical(1L)
      )
    ),
    
    identical(
      dim(fit$s_hat_lag),
      c(
        n_units,
        n_units,
        p_lags
      )
    ),
    
    identical(
      dim(fit$s_hat_unit_max),
      c(
        n_units,
        n_units
      )
    ),
    
    identical(
      dim(fit$s_hat_unit_rms),
      c(
        n_units,
        n_units
      )
    ),
    
    all(
      is.finite(
        fit$s_hat_lag
      )
    ),
    
    all(
      fit$s_hat_lag >= 0
    ),
    
    same(
      fit$score_lag,
      fit$s_hat_lag
    ),
    
    identical(
      dim(fit$selected_lag),
      dim(fit$s_hat_lag)
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
  
  all(
    conditions
  )
}


# Source files -----

source_paths <- c(
  "R/var_data.R",
  "R/var_coefficients.R",
  "R/network_blocks.R",
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
  source(path)
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
  "Compiling native engines...\n"
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

n_units <- 3L
m <- 2L
p_lags <- 2L

N <- n_units * m
k <- N * p_lags

T_obs <- 80L
simulation_burnin <- 100L


A1 <- diag(
  0.25,
  N
)

A2 <- diag(
  0.08,
  N
)


A1[1L, 3L] <- 0.14
A1[4L, 6L] <- -0.12
A1[5L, 1L] <- 0.10

A2[2L, 5L] <- 0.06
A2[6L, 3L] <- -0.05


set.seed(
  1201L
)


Y <- matrix(
  0,
  nrow =
    T_obs +
    simulation_burnin,
  ncol = N
)


for (
  t in 3:nrow(Y)
) {
  
  Y[t, ] <- as.numeric(
    A1 %*%
      Y[t - 1L, ] +
      A2 %*%
      Y[t - 2L, ] +
      stats::rnorm(
        N,
        sd = 0.5
      )
  )
}


Y <- Y[
  (
    simulation_burnin +
      1L
  ):nrow(Y),
  ,
  drop = FALSE
]


colnames(Y) <- paste0(
  "series_",
  seq_len(N)
)


Y_list <- lapply(
  seq_len(n_units),
  function(unit) {
    
    columns <- (
      (unit - 1L) *
        m +
        1L
    ):(
      unit * m
    )
    
    Y[
      ,
      columns,
      drop = FALSE
    ]
  }
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

cat(
  "\nHALF-T\n\n"
)


half_structure_all <-
  make_half_t_structure(
    Y_list = Y_list,
    p_lags = p_lags,
    global_grouping = "all"
  )


half_structure_sd <-
  make_half_t_structure(
    Y_list = Y_list,
    p_lags = p_lags,
    global_grouping =
      "self_diagonal"
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
    half_structure_sd$
      maps$
      lag_id,
    as.integer(
      expected_lag_id
    )
  ) &&
    same(
      half_structure_sd$
        maps$
        phi2,
      expected_phi2
    ) &&
    half_structure_all$
    maps$
    n_tau ==
    1L &&
    half_structure_sd$
    maps$
    n_tau ==
    2L
)


half_prior <-
  make_half_t_prior()


check(
  "Half-t canonical defaults",
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


half_prior_override <-
  make_half_t_prior(
    tau_df = 8,
    lambda_df = 4
  )


check(
  "Half-t prior override",
  half_prior_override$tau_df ==
    8 &&
    half_prior_override$lambda_df ==
    4
)


half_fit <- fit_half_t(
  Y_list = Y_list,
  p_lags = p_lags,
  global_grouping =
    "self_diagonal",
  chains = 2L,
  burnin = 50L,
  draws = 100L,
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
    draws = 100L,
    n_units = n_units,
    m = m,
    p_lags = p_lags
  )
)


check(
  "Half-t fit defaults",
  half_fit$prior$tau_df ==
    10 &&
    half_fit$prior$lambda_df ==
    3 &&
    half_fit$prior$tau_scale ==
    1 &&
    half_fit$prior$lambda_scale ==
    1 &&
    isTRUE(
      half_fit$
        control$
        use_asis
    )
)


check(
  "Half-t posterior scales",
  positive_finite(
    half_fit$tau2_mean
  ) &&
    positive_finite(
      half_fit$
        sigma2_standardized
    )
)


half_fit_alt <- fit_half_t(
  Y_list = Y_list,
  p_lags = p_lags,
  global_grouping = "all",
  chains = 1L,
  burnin = 20L,
  draws = 30L,
  thin = 1L,
  seed = 2201L,
  beta_algorithm =
    "bhattacharya",
  tau_df = 8,
  lambda_df = 4
)


check(
  "Half-t alternate path",
  identical(
    half_fit_alt$
      global_grouping,
    "all"
  ) &&
    !isTRUE(
      half_fit_alt$
        control$
        use_asis
    ) &&
    identical(
      half_fit_alt$
        beta_algorithm_resolved,
      "bhattacharya"
    ) &&
    half_fit_alt$
    prior$
    tau_df ==
    8 &&
    half_fit_alt$
    prior$
    lambda_df ==
    4
)


half_fit_repeat <- fit_half_t(
  Y_list = Y_list,
  p_lags = p_lags,
  global_grouping =
    "self_diagonal",
  chains = 2L,
  burnin = 50L,
  draws = 100L,
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
    )
)


# M3 -----

cat(
  "\nM3\n\n"
)


m3_structure <- make_m3_structure(
  Y_list = Y_list,
  p_lags = p_lags
)


expected_blocks <-
  n_units^2 *
  p_lags


m3_block_counts <- tabulate(
  as.vector(
    m3_structure$
      maps$
      block_id
  ),
  nbins = expected_blocks
)


check(
  "M3 structure",
  identical(
    dim(
      m3_structure$
        maps$
        lag_id
    ),
    c(k, N)
  ) &&
    identical(
      dim(
        m3_structure$
          maps$
          gc_id
      ),
      c(k, N)
    ) &&
    identical(
      dim(
        m3_structure$
          maps$
          gq_id
      ),
      c(k, N)
    ) &&
    identical(
      dim(
        m3_structure$
          maps$
          block_id
      ),
      c(k, N)
    ) &&
    m3_structure$
    maps$
    n_blocks ==
    expected_blocks &&
    all(
      m3_block_counts ==
        m^2
    ) &&
    identical(
      as.numeric(
        m3_structure$
          maps$
          phi2
      ),
      rep(
        1,
        p_lags
      )
    )
)


check(
  "M3 default grouping",
  identical(
    m3_structure$
      grouping$
      c_group,
    "self_foreign"
  ) &&
    identical(
      m3_structure$
        grouping$
        q_group,
      "self_foreign"
    ) &&
    m3_structure$
    maps$
    n_c ==
    2L &&
    m3_structure$
    maps$
    n_q ==
    2L
)


m3_state <- make_m3_state(
  structure = m3_structure,
  prior = make_m3_prior(
    structure = m3_structure
  )
)


check(
  "M3 state dimensions",
  identical(
    dim(
      m3_state$beta
    ),
    c(k, N)
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
    m3_structure$
    maps$
    n_c &&
    length(
      m3_state$psi_c
    ) ==
    m3_structure$
    maps$
    n_c &&
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
    m3_structure$
    maps$
    n_q &&
    m3_state$c2[1L] ==
    1 &&
    m3_state$psi_c[1L] ==
    1
)


m3_prior <- make_m3_prior(
  structure = m3_structure
)


expected_m3_tau_scale <-
  1 /
  sqrt(
    m3_structure$
      data$
      T_p
  )


expected_q_grid <- seq(
  -log(10),
  log(10),
  length.out = 17L
)


check(
  "M3 canonical defaults",
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
    )
)


check(
  "M3 q-grid defaults",
  same(
    m3_prior$q_grid,
    expected_q_grid
  ) &&
    same(
      m3_prior$q_prob,
      rep(
        1 / 17,
        17L
      )
    ) &&
    same(
      sum(
        m3_prior$q_prob
      ),
      1
    )
)


m3_prior_override <- make_m3_prior(
  structure = m3_structure,
  tau_df = 8,
  c_df = 6,
  lambda_df = 4
)


check(
  "M3 prior override",
  m3_prior_override$tau_df ==
    8 &&
    m3_prior_override$c_df ==
    6 &&
    m3_prior_override$
    lambda_df ==
    4
)


m3_control <- make_m3_control(
  burnin = 50L,
  draws = 100L
)


check(
  "M3 control defaults",
  identical(
    m3_control$q_update,
    "transport_mh"
  ) &&
    identical(
      m3_control$
        scale_transport,
      "none"
    ) &&
    same(
      m3_control$
        q_transport_global_probability,
      0.25
    ) &&
    isTRUE(
      m3_control$
        use_c_asis
    ) &&
    m3_control$
    c_asis_every ==
    1L &&
    m3_control$
    c_asis_slice_width ==
    1 &&
    m3_control$
    c_asis_maximum_step_out ==
    20L
)


m3_fit <- fit_m3(
  Y_list = Y_list,
  p_lags = p_lags,
  chains = 2L,
  burnin = 50L,
  draws = 100L,
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
    draws = 100L,
    n_units = n_units,
    m = m,
    p_lags = p_lags
  )
)


check(
  "M3 fit defaults",
  m3_fit$prior$tau_df ==
    10 &&
    m3_fit$prior$c_df ==
    5 &&
    m3_fit$prior$lambda_df ==
    3 &&
    same(
      m3_fit$
        prior$
        tau_scale,
      expected_m3_tau_scale
    ) &&
    identical(
      m3_fit$c_group,
      "self_foreign"
    ) &&
    identical(
      m3_fit$q_group,
      "self_foreign"
    ) &&
    identical(
      m3_fit$
        control$
        q_update,
      "transport_mh"
    ) &&
    identical(
      m3_fit$
        control$
        scale_transport,
      "none"
    ) &&
    isTRUE(
      m3_fit$
        control$
        use_c_asis
    )
)


check(
  "M3 pooled dimensions",
  length(
    m3_fit$tau2_mean
  ) ==
    1L &&
    length(
      m3_fit$c2_mean
    ) ==
    m3_structure$
    maps$
    n_c &&
    length(
      m3_fit$q_mean
    ) ==
    m3_structure$
    maps$
    n_q &&
    positive_finite(
      m3_fit$tau2_mean
    ) &&
    positive_finite(
      m3_fit$c2_mean
    ) &&
    all(
      is.finite(
        m3_fit$q_mean
      )
    )
)


m3_chain_contract <- vapply(
  m3_fit$chain_results,
  function(chain) {
    
    c_updates <-
      chain$c_asis$updates
    
    c_moves <-
      chain$c_asis$moves
    
    q_proposals <-
      chain$q_transport$proposals
    
    q_acceptance <-
      chain$q_transport$
      acceptance_rate
    
    conditions <- c(
      length(
        chain$tau2_mean
      ) ==
        1L,
      
      length(
        chain$c2_mean
      ) ==
        m3_structure$
        maps$
        n_c,
      
      length(
        chain$q_mean
      ) ==
        m3_structure$
        maps$
        n_q,
      
      identical(
        dim(
          chain$draws$c2
        ),
        c(
          chain$retained,
          m3_structure$
            maps$
            n_c
        )
      ),
      
      identical(
        dim(
          chain$draws$q
        ),
        c(
          chain$retained,
          m3_structure$
            maps$
            n_q
        )
      ),
      
      isTRUE(
        chain$c_asis$enabled
      ),
      
      length(c_updates) ==
        m3_structure$
        maps$
        n_c,
      
      length(c_moves) ==
        m3_structure$
        maps$
        n_c,
      
      c_updates[1L] ==
        0L,
      
      c_moves[1L] ==
        0L,
      
      all(
        c_updates[-1L] >
          0L
      ),
      
      all(
        c_moves >=
          0L
      ),
      
      length(q_proposals) ==
        m3_structure$
        maps$
        n_q,
      
      length(q_acceptance) ==
        m3_structure$
        maps$
        n_q,
      
      all(
        q_proposals >
          0L
      ),
      
      all(
        is.finite(
          q_acceptance
        )
      ),
      
      all(
        q_acceptance >=
          0 &
          q_acceptance <=
          1
      )
    )
    
    all(
      conditions
    )
  },
  logical(1L)
)


check(
  "M3 chain-level contracts",
  all(
    m3_chain_contract
  )
)


m3_fit_alt <- fit_m3(
  Y_list = Y_list,
  p_lags = p_lags,
  chains = 1L,
  burnin = 20L,
  draws = 30L,
  thin = 1L,
  seed = 3201L,
  use_c_asis = FALSE,
  tau_df = 8,
  c_df = 6,
  lambda_df = 4,
  keep_chain_results = TRUE
)


m3_alt_chain <-
  m3_fit_alt$
  chain_results[[1L]]


check(
  "M3 alternate path",
  !isTRUE(
    m3_fit_alt$
      control$
      use_c_asis
  ) &&
    !isTRUE(
      m3_alt_chain$
        c_asis$
        enabled
    ) &&
    length(
      m3_alt_chain$
        c_asis$
        updates
    ) ==
    m3_structure$
    maps$
    n_c &&
    all(
      m3_alt_chain$
        c_asis$
        updates ==
        0L
    ) &&
    all(
      m3_alt_chain$
        c_asis$
        moves ==
        0L
    ) &&
    m3_fit_alt$
    prior$
    tau_df ==
    8 &&
    m3_fit_alt$
    prior$
    c_df ==
    6 &&
    m3_fit_alt$
    prior$
    lambda_df ==
    4
)


m3_fit_repeat <- fit_m3(
  Y_list = Y_list,
  p_lags = p_lags,
  chains = 2L,
  burnin = 50L,
  draws = 100L,
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

cat(
  "\nGIGG\n\n"
)


gigg_structure <- make_gigg_structure(
  Y_list = Y_list,
  p_lags = p_lags
)


gigg_block_counts <- tabulate(
  as.vector(
    gigg_structure$
      maps$
      block_id
  ),
  nbins = expected_blocks
)


check(
  "GIGG structure",
  gigg_structure$
    maps$
    n_blocks ==
    expected_blocks &&
    all(
      gigg_block_counts ==
        m^2
    ) &&
    isTRUE(
      gigg_structure$
        preprocessing$
        standardize
    ) &&
    !isTRUE(
      gigg_structure$
        preprocessing$
        centering
    )
)


gigg_prior <- make_gigg_prior(
  gigg_structure
)


expected_gigg_tau_scale <-
  3 /
  sqrt(
    5 *
      gigg_structure$
      data$
      T_p
  )


check(
  "GIGG canonical defaults",
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
  same(
    gigg_state$tau2,
    1 /
      gigg_structure$
      data$
      T_p
  ) &&
    length(
      gigg_state$gamma2
    ) ==
    expected_blocks &&
    all(
      gigg_state$gamma2 ==
        1
    ) &&
    identical(
      dim(
        gigg_state$lambda2
      ),
      c(k, N)
    )
)


gigg_fit <- fit_gigg(
  Y_list = Y_list,
  p_lags = p_lags,
  chains = 2L,
  burnin = 50L,
  draws = 100L,
  thin = 1L,
  seed = 4101L,
  keep_chain_results = TRUE
)


check(
  "GIGG fit output",
  check_fit_output(
    fit = gigg_fit,
    model = "gigg",
    chains = 2L,
    draws = 100L,
    n_units = n_units,
    m = m,
    p_lags = p_lags
  )
)


check(
  "GIGG fit defaults",
  gigg_fit$prior$tau_df ==
    5 &&
    same(
      gigg_fit$
        prior$
        tau_scale,
      expected_gigg_tau_scale
    ) &&
    gigg_fit$prior$
    gamma_shape ==
    0.5 &&
    gigg_fit$prior$
    gamma_rate ==
    1 &&
    gigg_fit$prior$
    lambda_shape ==
    2.5 &&
    gigg_fit$prior$
    lambda_scale ==
    1 &&
    isTRUE(
      gigg_fit$
        control$
        use_global_tau
    ) &&
    isTRUE(
      gigg_fit$
        control$
        use_asis
    ) &&
    gigg_fit$
    control$
    asis_every ==
    1L
)


check(
  "GIGG posterior scales",
  positive_finite(
    gigg_fit$tau2_mean
  ) &&
    positive_finite(
      gigg_fit$gamma2_mean
    ) &&
    positive_finite(
      gigg_fit$
        tau2_gamma2_mean
    )
)


gigg_asis_contract <- vapply(
  gigg_fit$chain_results,
  function(chain) {
    
    isTRUE(
      chain$asis$enabled
    ) &&
      isTRUE(
        chain$asis$
          global_block$
          applicable
      )
  },
  logical(1L)
)


check(
  "GIGG ASIS path",
  all(
    gigg_asis_contract
  )
)


gigg_fit_alt <- fit_gigg(
  Y_list = Y_list,
  p_lags = p_lags,
  chains = 1L,
  burnin = 20L,
  draws = 30L,
  thin = 1L,
  seed = 4201L,
  use_global_tau = FALSE,
  use_asis = FALSE,
  keep_chain_results = TRUE
)


gigg_alt_chain <-
  gigg_fit_alt$
  chain_results[[1L]]


check(
  "GIGG alternate path",
  !isTRUE(
    gigg_fit_alt$
      control$
      use_global_tau
  ) &&
    !isTRUE(
      gigg_fit_alt$
        control$
        use_asis
    ) &&
    !isTRUE(
      gigg_alt_chain$
        asis$
        enabled
    ) &&
    !isTRUE(
      gigg_alt_chain$
        asis$
        global_block$
        applicable
    )
)


gigg_fit_repeat <- fit_gigg(
  Y_list = Y_list,
  p_lags = p_lags,
  chains = 2L,
  burnin = 50L,
  draws = 100L,
  thin = 1L,
  seed = 4101L
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
      gigg_fit$gamma2_mean,
      gigg_fit_repeat$gamma2_mean
    ) &&
    same(
      gigg_fit$
        tau2_gamma2_mean,
      gigg_fit_repeat$
        tau2_gamma2_mean
    )
)


# Final report -----

separator <- paste0(
  rep(
    "=",
    72L
  ),
  collapse = ""
)


cat(
  "\n",
  separator,
  "\nBAYESIAN MODEL VALIDATION SUMMARY\n",
  separator,
  "\n\n",
  
  "Common VAR\n",
  sprintf(
    "  Units / variables / lags ............ %d / %d / %d\n",
    n_units,
    m,
    p_lags
  ),
  sprintf(
    "  Effective observations .............. %d\n",
    prepared$T_eff
  ),
  "  Native engines ....................... PASSED\n",
  
  "\nHalf-t\n",
  sprintf(
    "  Prior df (tau / lambda) .............. %g / %g\n",
    half_fit$prior$tau_df,
    half_fit$prior$lambda_df
  ),
  "  Lag variance decay ................... 1 / lag^2\n",
  "  Groupings tested ..................... all, self_diagonal\n",
  "  Default ASIS path .................... PASSED\n",
  "  Bhattacharya path .................... PASSED\n",
  "  Fixed-seed reproducibility ........... PASSED\n",
  
  "\nM3\n",
  sprintf(
    "  Prior df (tau / c / lambda) .......... %g / %g / %g\n",
    m3_fit$prior$tau_df,
    m3_fit$prior$c_df,
    m3_fit$prior$lambda_df
  ),
  sprintf(
    "  Default groups (c / q) ............... %s / %s\n",
    m3_fit$c_group,
    m3_fit$q_group
  ),
  sprintf(
    "  c groups / q groups / blocks ........ %d / %d / %d\n",
    m3_structure$maps$n_c,
    m3_structure$maps$n_q,
    m3_structure$maps$n_blocks
  ),
  sprintf(
    "  Global tau scale ..................... %.6f\n",
    m3_fit$prior$tau_scale
  ),
  "  Lag variance weights ................. all 1\n",
  sprintf(
    "  q transport .......................... %s\n",
    m3_fit$control$q_update
  ),
  "  c-ASIS default path .................. PASSED\n",
  "  c-ASIS disabled path ................. PASSED\n",
  "  Fixed-seed reproducibility ........... PASSED\n",
  
  "\nGIGG\n",
  sprintf(
    "  tau df ............................... %g\n",
    gigg_fit$prior$tau_df
  ),
  sprintf(
    "  gamma shape / rate ................... %g / %g\n",
    gigg_fit$prior$gamma_shape,
    gigg_fit$prior$gamma_rate
  ),
  sprintf(
    "  lambda shape / scale ................. %g / %g\n",
    gigg_fit$prior$lambda_shape,
    gigg_fit$prior$lambda_scale
  ),
  sprintf(
    "  Global tau scale ..................... %.6f\n",
    gigg_fit$prior$tau_scale
  ),
  "  Global tau + ASIS path ............... PASSED\n",
  "  Disabled path ........................ PASSED\n",
  "  Fixed-seed reproducibility ........... PASSED\n",
  
  "\n",
  separator,
  "\nALL BAYESIAN MODEL CHECKS PASSED\n",
  separator,
  "\n",
  sep = ""
)