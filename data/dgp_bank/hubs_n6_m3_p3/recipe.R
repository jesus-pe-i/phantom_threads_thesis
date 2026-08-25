source("R/network_diagnostics.R")
source("R/utils_network.R")
source("R/simulate_var_v3.R")

# DGP: hubs_n6_m3_p3
dgp <- eval({
  {
      n_units <- 6L
      m <- 3L
      p_lags <- 3L
      make_A_hub(n_units = n_units, m = m, p_lags = p_lags, sender_hub_units = c(1L), 
          receiver_hub_units = c(6L), self_diag = c(0.47, 0.06, 
              0.04), self_off = c(0.47, 0.06, 0.04)/6, self_diag_sign = 1, 
          self_off_sign = 1, self_off_pattern = "alternating", 
          sender_hub_diag = c(0.05, 0.035, 0.025), receiver_hub_diag = c(0.055, 
              0.03, 0.022), hub_to_hub_diag = c(0.074, 0.045, 0.03), 
          background_diag = c(0, 0, 0), sender_hub_diag_sign = 1, 
          sender_hub_off_sign = 1, receiver_hub_diag_sign = 1, 
          receiver_hub_off_sign = 1, hub_to_hub_diag_sign = 1, 
          hub_to_hub_off_sign = 1, background_diag_sign = 1, background_off_sign = 1, 
          sender_hub_off_pattern = "positive", receiver_hub_off_pattern = "positive", 
          hub_to_hub_off_pattern = "positive", background_off_pattern = "positive", 
          off_ratio = 0.9, coef_jitter_log_sd = 0.2, jitter_seed = 11L)
  }
})

A_list <- dgp$A_list
meta <- dgp$meta
