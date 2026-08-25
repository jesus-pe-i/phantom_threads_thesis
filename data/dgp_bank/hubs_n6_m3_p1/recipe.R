source("R/network_diagnostics.R")
source("R/utils_network.R")
source("R/simulate_var_v3.R")

# DGP: hubs_n6_m3_p1
dgp <- eval({
  {
      n_units <- 6L
      m <- 3L
      p_lags <- 1L
      make_A_hub(n_units = n_units, m = m, p_lags = p_lags, sender_hub_units = c(1L), 
          receiver_hub_units = c(6L), self_diag = c(0.56), self_off = c(0.56)/6, 
          self_diag_sign = 1, self_off_sign = 1, self_off_pattern = "alternating", 
          sender_hub_diag = c(0.07), receiver_hub_diag = c(0.055), 
          hub_to_hub_diag = c(0.09), background_diag = c(0), sender_hub_diag_sign = 1, 
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
