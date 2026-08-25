source("R/network_diagnostics.R")
source("R/utils_network.R")
source("R/simulate_var_v3.R")

# DGP: complex_n40_m4_p2
dgp <- eval({
  {
      n_units <- 40L
      m <- 4L
      p_lags <- 2L
      phi_scale <- 1
      phi_by_lag <- lapply(phi_by_lag, function(phi) {
          phi_scale * phi
      })
      phi_by_lag <- list(c(-0.05, 0.2, 0.45, 0.7), c(0, 0.02, 0.04, 
          0.06))
      phi_scale <- 0.915
      phi_by_lag <- lapply(phi_by_lag, function(phi) {
          phi_scale * phi
      })
      within_ratio <- c(0.1, 0.1)
      foreign_ratio <- c(0.445, 0.8)
      role_weights <- list(c(community_a = 0.6, community_b = 0.6, 
          community_bridge = 0.3, disconnected = 0.5, receiver_hub = 0.8, 
          core_out = 1), c(community_core = 1, community_a = 1, 
          community_b = 1, lag_2_bridge = 0.5, disconnected = 0.8))
      make_A_complex_scaling(n_units = n_units, phi_by_lag = phi_by_lag, 
          within_ratio = within_ratio, foreign_ratio = foreign_ratio, 
          role_weights = role_weights, foreign_off_ratio = 0.8, 
          foreign_negative_share = 0.4, self_off_pattern = "positive", 
          self_phi_assignment = "balanced_irregular")
  }
})

A_list <- dgp$A_list
meta <- dgp$meta
