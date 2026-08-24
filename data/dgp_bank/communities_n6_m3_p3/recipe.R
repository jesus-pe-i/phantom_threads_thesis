source("R/network_diagnostics.R")
source("R/utils_network.R")
source("R/simulate_var_v3.R")

# DGP: communities_n6_m3_p3
dgp <- eval({
  {
      n_units <- 6L
      m <- 3L
      p_lags <- 3L
      make_A_community_block(n_units = n_units, m = m, p_lags = p_lags, 
          communities = list(c(1L, 2L, 3L), c(4L, 5L, 6L)), self_diag = c(0.43, 
              0.04, 0.02), self_off = c(0.3, 0.1, 0.02)/6, self_diag_sign = 1, 
          self_off_sign = 1, self_off_pattern = "alternating", 
          within_diag = c(0.04, 0.025, 0.01), between_diag = c(0, 
              0, 0), within_diag_sign = 1, within_off_sign = 1, 
          between_diag_sign = 1, between_off_sign = 1, community_off_pattern = "positive", 
          off_ratio = 0.9, coef_jitter_log_sd = 0.2, jitter_seed = 11L)
  }
})

A_list <- dgp$A_list
meta <- dgp$meta
