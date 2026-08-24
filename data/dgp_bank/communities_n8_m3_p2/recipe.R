source("R/network_diagnostics.R")
source("R/utils_network.R")
source("R/simulate_var_v3.R")

# DGP: communities_n8_m3_p2
dgp <- eval({
  {
      n_units <- 8L
      m <- 3L
      p_lags <- 2L
      make_A_community_block(n_units = n_units, m = m, p_lags = p_lags, 
          communities = list(c(1L, 2L), c(3L, 4L, 5L, 6L), c(7L, 
              8L)), self_diag = c(0.4, 0.037), self_off = c(0.37, 
              0.037)/6, self_diag_sign = 1, self_off_sign = 1, 
          self_off_pattern = "alternating", within_diag = c(0.04, 
              0.02), between_diag = c(0, 0), within_diag_sign = 1, 
          within_off_sign = 1, between_diag_sign = 1, between_off_sign = 1, 
          community_off_pattern = "positive", off_ratio = 0.9, 
          coef_jitter_log_sd = 0.2, jitter_seed = 11L)
  }
})

A_list <- dgp$A_list
meta <- dgp$meta
