# Computes the network-recovery metrics used in benchmark experiments.
#
# Detection is evaluated on cross-unit relationships only. Magnitude metrics
# use the same block-strength definitions as the fitted benchmark models.


# Helpers -----

benchmark_network_mask <- function(n_units, p_lags = NULL) {
  if (is.null(p_lags)) {
    mask <- matrix(TRUE, n_units, n_units)
    diag(mask) <- FALSE
    return(mask)
  }

  mask <- array(TRUE, c(n_units, n_units, p_lags))

  for (lag in seq_len(p_lags)) {
    diag(mask[, , lag]) <- FALSE
  }

  mask
}


benchmark_pr_auc <- function(y_true, score) {
  y_true <- as.integer(y_true)
  score <- as.numeric(score)

  if (length(y_true) != length(score) || anyNA(y_true) || any(!is.finite(score))) {
    stop("y_true and score must be finite vectors of equal length")
  }

  if (any(!y_true %in% c(0L, 1L))) {
    stop("y_true must be binary")
  }

  positives <- sum(y_true)

  if (positives == 0L) {
    return(NA_real_)
  }

  order_index <- order(score, decreasing = TRUE)
  y_true <- y_true[order_index]
  score <- score[order_index]

  tied_end <- c(which(diff(score) != 0), length(score))
  true_positive <- cumsum(y_true)[tied_end]
  false_positive <- tied_end - true_positive

  recall <- true_positive / positives
  precision <- true_positive / (true_positive + false_positive)

  sum(diff(c(0, recall)) * precision)
}


benchmark_rmse <- function(truth, estimate) {
  sqrt(mean((estimate - truth)^2))
}


# Benchmark evaluation -----

evaluate_benchmark_fit <- function(fit, truth) {
  required_truth <- c("G_lag", "G_unit", "s_lag", "s_unit_rms")

  if (!all(required_truth %in% names(truth))) {
    stop("truth must contain G_lag, G_unit, s_lag, and s_unit_rms")
  }

  if (!identical(dim(fit$s_hat_lag), dim(truth$s_lag)) ||
      !identical(dim(fit$score_lag), dim(truth$G_lag)) ||
      !identical(dim(fit$s_hat_unit_rms), dim(truth$s_unit_rms)) ||
      !identical(dim(truth$G_unit), dim(truth$s_unit_rms))) {
    stop("Benchmark fit and truth dimensions are inconsistent")
  }

  n_units <- dim(truth$G_lag)[1L]
  p_lags <- dim(truth$G_lag)[3L]

  lag_mask <- benchmark_network_mask(n_units, p_lags)
  unit_mask <- benchmark_network_mask(n_units)

  lag_truth <- truth$G_lag[lag_mask]
  lag_score <- fit$score_lag[lag_mask]
  lag_true_strength <- truth$s_lag[lag_mask]
  lag_hat_strength <- fit$s_hat_lag[lag_mask]

  unit_truth <- truth$G_unit[unit_mask]
  unit_score <- apply(fit$score_lag, c(1L, 2L), max)[unit_mask]
  unit_true_strength <- truth$s_unit_rms[unit_mask]
  unit_hat_strength <- fit$s_hat_unit_rms[unit_mask]

  lag_active <- as.logical(lag_truth)
  lag_inactive <- !lag_active
  unit_active <- as.logical(unit_truth)
  unit_inactive <- !unit_active

  self_true <- unlist(lapply(seq_len(p_lags), function(lag) diag(truth$s_lag[, , lag])), use.names = FALSE)
  self_hat <- unlist(lapply(seq_len(p_lags), function(lag) diag(fit$s_hat_lag[, , lag])), use.names = FALSE)
  self_active <- unlist(lapply(seq_len(p_lags), function(lag) diag(truth$G_lag[, , lag])), use.names = FALSE)

  data.frame(
    lag_positives = sum(lag_active),
    lag_negatives = sum(lag_inactive),
    unit_positives = sum(unit_active),
    unit_negatives = sum(unit_inactive),

    lag_pr_auc = benchmark_pr_auc(lag_active, lag_score),
    unit_max_pr_auc = benchmark_pr_auc(unit_active, unit_score),

    unit_rms_rmse = benchmark_rmse(unit_true_strength, unit_hat_strength),
    unit_rms_rmse_active = if (any(unit_active)) benchmark_rmse(unit_true_strength[unit_active], unit_hat_strength[unit_active]) else NA_real_,
    unit_rms_inactive_leakage = if (any(unit_inactive)) mean(abs(unit_hat_strength[unit_inactive])) else NA_real_,

    mean_false = if (any(lag_inactive)) mean(abs(lag_hat_strength[lag_inactive])) else NA_real_,
    maximum_false = if (any(lag_inactive)) max(abs(lag_hat_strength[lag_inactive])) else NA_real_,
    self_block_rmse = if (any(self_active)) benchmark_rmse(self_true[self_active], self_hat[self_active]) else NA_real_,

    row.names = NULL
  )
}
