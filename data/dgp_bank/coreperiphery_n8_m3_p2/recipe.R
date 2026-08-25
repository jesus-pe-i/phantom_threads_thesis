source("R/network_diagnostics.R")
source("R/utils_network.R")
source("R/simulate_var_v3.R")

# DGP: coreperiphery_n8_m3_p2
dgp <- eval({
  {
      n_units <- 8L
      m <- 3L
      p_lags <- 2L
      make_A_core_periphery(n_units = n_units, m = m, p_lags = p_lags, 
          core_units = c(1L, 2L, 3L), self_diag = c(0.41, 0.05), 
          self_off = c(0.4, 0.05)/6, self_diag_sign = 1, self_off_sign = 1, 
          self_off_pattern = "alternating", core_core_diag = c(0.05, 
              0.029), core_to_periphery_diag = c(0.045, 0.025), 
          periphery_to_core_diag = c(0, 0), periphery_to_periphery_diag = c(0, 
              0), core_core_diag_sign = 1, core_core_off_sign = 1, 
          core_to_periphery_diag_sign = 1, core_to_periphery_off_sign = 1, 
          periphery_to_core_diag_sign = 1, periphery_to_core_off_sign = 1, 
          periphery_to_periphery_diag_sign = 1, periphery_to_periphery_off_sign = 1, 
          core_core_off_pattern = "positive", core_to_periphery_off_pattern = "positive", 
          periphery_to_core_off_pattern = "positive", periphery_to_periphery_off_pattern = "positive", 
          off_ratio = 0.9, coef_jitter_log_sd = 0.2, jitter_seed = 11L)
  }
})

A_list <- dgp$A_list
meta <- dgp$meta
