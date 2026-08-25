source("R/network_diagnostics.R")
source("R/utils_network.R")
source("R/simulate_var_v3.R")

# DGP: coreperiphery_n6_m3_p5
dgp <- eval({
  {
      n_units <- 6L
      m <- 3L
      p_lags <- 5L
      make_A_core_periphery(n_units = n_units, m = m, p_lags = p_lags, 
          core_units = c(1L, 2L), self_diag = c(0.37, 0.04, 0.035, 
              0.03, 0.02), self_off = c(0.35, 0.04, 0.035, 0.03, 
              0.02)/6, self_diag_sign = 1, self_off_sign = 1, self_off_pattern = "alternating", 
          core_core_diag = c(0.055, 0.03, 0.023, 0.02, 0.015), 
          core_to_periphery_diag = c(0.045, 0.022, 0.022, 0.018, 
              0.01), periphery_to_core_diag = c(0, 0, 0, 0, 0), 
          periphery_to_periphery_diag = c(0, 0, 0, 0, 0), core_core_diag_sign = 1, 
          core_core_off_sign = 1, core_to_periphery_diag_sign = 1, 
          core_to_periphery_off_sign = 1, periphery_to_core_diag_sign = 1, 
          periphery_to_core_off_sign = 1, periphery_to_periphery_diag_sign = 1, 
          periphery_to_periphery_off_sign = 1, core_core_off_pattern = "positive", 
          core_to_periphery_off_pattern = "positive", periphery_to_core_off_pattern = "positive", 
          periphery_to_periphery_off_pattern = "positive", off_ratio = 0.9, 
          coef_jitter_log_sd = 0.2, jitter_seed = 11L)
  }
})

A_list <- dgp$A_list
meta <- dgp$meta
