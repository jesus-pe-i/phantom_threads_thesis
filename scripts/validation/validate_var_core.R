# Permanent validation of the core block-VAR transformations.
#
# Checks data stacking, lag construction, coefficient orientation,
# scale back-transformation, block strengths, and lag aggregation.


rm(list = ls())

source("R/var_data.R")
source("R/var_coefficients.R")
source("R/network_blocks.R")


# Data preparation -----

Y_list <- list(
  unit_1 = cbind(
    x1 = 1:8,
    x2 = 11:18
  ),
  unit_2 = cbind(
    x1 = 21:28,
    x2 = 31:38
  )
)

prepared <- prepare_var_data(
  Y_list = Y_list,
  p_lags = 2L,
  standardize = FALSE
)

expected_Y_full <- cbind(
  1:8,
  11:18,
  21:28,
  31:38
)

expected_Y <- expected_Y_full[
  3:8,
  ,
  drop = FALSE
]

expected_X <- cbind(
  expected_Y_full[
    2:7,
    ,
    drop = FALSE
  ],
  expected_Y_full[
    1:6,
    ,
    drop = FALSE
  ]
)

data_ok <- all(
  isTRUE(
    all.equal(
      prepared$Y_full,
      unname(expected_Y_full),
      tolerance = 0
    )
  ),
  isTRUE(
    all.equal(
      prepared$Y,
      unname(expected_Y),
      tolerance = 0
    )
  ),
  isTRUE(
    all.equal(
      prepared$X,
      unname(expected_X),
      tolerance = 0
    )
  )
)


# Coefficient orientation -----

A1 <- matrix(
  c(
    1,  2,  3,  4,
    5,  6,  7,  8,
    9, 10, 11, 12,
    13, 14, 15, 16
  ),
  nrow = 4L,
  byrow = TRUE
)

A2 <- matrix(
  c(
    17, 18, 19, 20,
    21, 22, 23, 24,
    25, 26, 27, 28,
    29, 30, 31, 32
  ),
  nrow = 4L,
  byrow = TRUE
)

A_list <- list(
  lag_1 = A1,
  lag_2 = A2
)

beta <- A_list_to_beta(
  A_list
)

A_roundtrip <- beta_to_A_list(
  beta = beta,
  p_lags = 2L
)

orientation_ok <- all(
  beta[1L, 2L] == A1[2L, 1L],
  beta[4L, 1L] == A1[1L, 4L],
  beta[5L, 3L] == A2[3L, 1L]
)

coefficient_roundtrip_ok <- all(
  isTRUE(
    all.equal(
      unname(A_roundtrip[[1L]]),
      A1,
      tolerance = 0
    )
  ),
  isTRUE(
    all.equal(
      unname(A_roundtrip[[2L]]),
      A2,
      tolerance = 0
    )
  )
)


# Scale back-transformation -----

scale <- c(
  2,
  4,
  5,
  10
)

predictor_scale <- rep(
  scale,
  times = 2L
)

beta_standardized <- beta *
  outer(
    predictor_scale,
    1 / scale,
    "*"
  )

beta_back <- backtransform_var_beta(
  beta = beta_standardized,
  scale = scale,
  p_lags = 2L
)

backtransform_ok <- isTRUE(
  all.equal(
    beta_back,
    beta,
    tolerance = 1e-14
  )
)


# Block network -----

B12_lag1 <- matrix(
  c(
    3, 4,
    0, 0
  ),
  nrow = 2L,
  byrow = TRUE
)

B21_lag2 <- matrix(
  c(
    0, 6,
    8, 0
  ),
  nrow = 2L,
  byrow = TRUE
)

A_network_1 <- matrix(
  0,
  nrow = 4L,
  ncol = 4L
)

A_network_2 <- A_network_1

# Unit 1 -> Unit 2
A_network_1[
  3:4,
  1:2
] <- B12_lag1

# Unit 2 -> Unit 1
A_network_2[
  1:2,
  3:4
] <- B21_lag2

s_lag <- A_list_to_s_lag(
  A_list = list(
    lag_1 = A_network_1,
    lag_2 = A_network_2
  ),
  n_units = 2L,
  m = 2L
)

expected_1_to_2 <- 5 / 2
expected_2_to_1 <- 10 / 2

strength_ok <- all(
  isTRUE(
    all.equal(
      s_lag[2L, 1L, 1L],
      expected_1_to_2,
      tolerance = 1e-14
    )
  ),
  isTRUE(
    all.equal(
      s_lag[1L, 2L, 2L],
      expected_2_to_1,
      tolerance = 1e-14
    )
  ),
  s_lag[1L, 2L, 1L] == 0,
  s_lag[2L, 1L, 2L] == 0
)


# Lag aggregation -----

s_unit_max <- s_lag_to_unit(
  s_lag,
  method = "max"
)

s_unit_rms <- s_lag_to_unit(
  s_lag,
  method = "rms"
)

aggregation_ok <- all(
  isTRUE(
    all.equal(
      s_unit_max[2L, 1L],
      expected_1_to_2,
      tolerance = 1e-14
    )
  ),
  isTRUE(
    all.equal(
      s_unit_rms[2L, 1L],
      expected_1_to_2 / sqrt(2),
      tolerance = 1e-14
    )
  ),
  isTRUE(
    all.equal(
      s_unit_max[1L, 2L],
      expected_2_to_1,
      tolerance = 1e-14
    )
  )
)

G_unit <- G_lag_to_G_unit(
  s_lag > 0
)

graph_ok <- all(
  G_unit[2L, 1L],
  G_unit[1L, 2L],
  !G_unit[1L, 1L],
  !G_unit[2L, 2L]
)


# Validation -----

checks <- c(
  data_preparation = data_ok,
  coefficient_orientation = orientation_ok,
  coefficient_roundtrip = coefficient_roundtrip_ok,
  scale_backtransform = backtransform_ok,
  block_strength = strength_ok,
  lag_aggregation = aggregation_ok,
  graph_aggregation = graph_ok
)

stopifnot(
  all(checks)
)


# Console report -----

cat("\nVAR CORE VALIDATION\n\n")

print(checks)

cat("\nCORE DIMENSIONS\n")

print(
  c(
    n_units = prepared$n_units,
    m = prepared$m,
    n_series = prepared$n_series,
    p_lags = prepared$p_lags,
    T_full = prepared$T_full,
    T_eff = prepared$T_eff
  )
)

cat("\nNETWORK LAG 1\n")
print(s_lag[, , 1L])

cat("\nNETWORK LAG 2\n")
print(s_lag[, , 2L])

cat("\nUNIT MAX NETWORK\n")
print(s_unit_max)

cat("\nUNIT RMS NETWORK\n")
print(round(s_unit_rms, 6))

cat("\nVAR CORE VALIDATION PASSED\n")