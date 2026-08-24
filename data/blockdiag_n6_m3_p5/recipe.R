source("R/network_diagnostics.R")
source("R/utils_network.R")
source("R/simulate_var_v3.R")

# DGP: blockdiag_n6_m3_p5
dgp <- eval({
  {
      n_units <- 6L
      m <- 3L
      p_lags <- 5L
      make_A_block_diagonal(n_units = n_units, m = m, p_lags = p_lags, 
          self_diag = c(0.47, 0.08, 0.05, 0.025, 0.0125), self_off = c(0.47, 
              0.08, 0.05, 0.025, 0.0125)/6, self_diag_sign = 1, 
          self_off_sign = 1, self_off_pattern = "alternating", 
          coef_jitter_log_sd = 0.2, jitter_seed = 11L)
  }
})

A_list <- dgp$A_list
meta <- dgp$meta
