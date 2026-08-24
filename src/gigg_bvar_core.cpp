// Fast GIGG BVAR one-chain engine.
//
// Lag-receiver-sender block GIG shrinkage with a global half-t
// scale, exact native GIG block updates, coefficient-local
// inverse-gamma scales, and the two GIGG ASIS interweaving moves.

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>

#include "gigg_bvar_updates.h"


namespace {

using Clock = std::chrono::steady_clock;


inline std::uint64_t checked_total_iterations(
    const int burnin,
    const int draws,
    const int thin) {

  if (burnin < 0) {
    Rcpp::stop(
      "burnin must be non-negative"
    );
  }

  if (draws < 1) {
    Rcpp::stop(
      "draws must be positive"
    );
  }

  if (thin < 1) {
    Rcpp::stop(
      "thin must be positive"
    );
  }

  const std::uint64_t burnin_u =
    static_cast<std::uint64_t>(
      burnin
    );

  const std::uint64_t draws_u =
    static_cast<std::uint64_t>(
      draws
    );

  const std::uint64_t thin_u =
    static_cast<std::uint64_t>(
      thin
    );

  if (
    draws_u >
    (
      std::numeric_limits<std::uint64_t>::max() -
      burnin_u
    ) /
    thin_u
  ) {
    Rcpp::stop(
      "Requested chain length overflows the iteration counter"
    );
  }

  return burnin_u +
    draws_u * thin_u;
}


inline double elapsed_seconds(
    const Clock::time_point& start,
    const Clock::time_point& end) {

  return std::chrono::duration<double>(
    end - start
  ).count();
}


inline bool list_has(
    const Rcpp::List& object,
    const char* name) {

  return object.containsElementNamed(name);
}


struct PairMonitorIndex {

  arma::uvec row;
  arma::uvec equation;
};


inline PairMonitorIndex validate_pair_monitor(
    const arma::imat& raw,
    const arma::uword k,
    const arma::uword N,
    const std::string& name) {

  if (raw.n_rows == 0u) {

    if (
      raw.n_cols != 0u &&
      raw.n_cols != 2u
    ) {
      Rcpp::stop(
        "Empty " + name + " must have zero or two columns"
      );
    }

    PairMonitorIndex result;

    result.row.set_size(0u);
    result.equation.set_size(0u);

    return result;
  }

  if (raw.n_cols != 2u) {
    Rcpp::stop(
      name + " must have columns row and equation"
    );
  }

  PairMonitorIndex result;

  result.row.set_size(raw.n_rows);
  result.equation.set_size(raw.n_rows);

  for (
    arma::uword j = 0u;
    j < raw.n_rows;
    ++j
  ) {

    const int row = raw(j, 0u);
    const int equation = raw(j, 1u);

    if (
      row < 1 ||
      row > static_cast<int>(k) ||
      equation < 1 ||
      equation > static_cast<int>(N)
    ) {
      Rcpp::stop(
        name +
        " contains an index outside the coefficient matrix"
      );
    }

    for (
      arma::uword previous = 0u;
      previous < j;
      ++previous
    ) {

      if (
        raw(previous, 0u) == row &&
        raw(previous, 1u) == equation
      ) {
        Rcpp::stop(
          name + " contains duplicated coefficient indices"
        );
      }
    }

    result.row[j] =
      static_cast<arma::uword>(row - 1);

    result.equation[j] =
      static_cast<arma::uword>(equation - 1);
  }

  return result;
}


inline arma::uvec validate_block_monitor(
    const arma::ivec& raw,
    const arma::uword n_blocks) {

  arma::uvec result(raw.n_elem);

  for (
    arma::uword j = 0u;
    j < raw.n_elem;
    ++j
  ) {

    if (
      raw[j] < 1 ||
      raw[j] > static_cast<int>(n_blocks)
    ) {
      Rcpp::stop(
        "monitor_gamma contains an index outside 1:n_blocks"
      );
    }

    for (
      arma::uword previous = 0u;
      previous < j;
      ++previous
    ) {

      if (raw[previous] == raw[j]) {
        Rcpp::stop(
          "monitor_gamma contains duplicated block indices"
        );
      }
    }

    result[j] =
      static_cast<arma::uword>(raw[j] - 1);
  }

  return result;
}


struct ChainTiming {

  double omega_seconds = 0.0;
  double beta_seconds = 0.0;
  double sigma2_seconds = 0.0;
  double tau_seconds = 0.0;
  double gamma_seconds = 0.0;
  double lambda_seconds = 0.0;
  double block_asis_seconds = 0.0;
  double global_asis_seconds = 0.0;
  double storage_seconds = 0.0;
};


inline arma::vec tau_gamma_product(
    const double tau2,
    const arma::vec& gamma2) {

  return tau2 * gamma2;
}


inline arma::mat effective_product(
    const double tau2,
    const arma::vec& gamma2,
    const arma::mat& lambda2,
    const arma::imat& block_id_zero) {

  arma::mat result(
    lambda2.n_rows,
    lambda2.n_cols
  );

  for (
    arma::uword equation = 0u;
    equation < lambda2.n_cols;
    ++equation
  ) {
    for (
      arma::uword row = 0u;
      row < lambda2.n_rows;
      ++row
    ) {

      const arma::uword block =
        static_cast<arma::uword>(
          block_id_zero(row, equation)
        );

      result(row, equation) =
        tau2 *
        gamma2[block] *
        lambda2(row, equation);
    }
  }

  return result;
}


}  // anonymous namespace


// ============================================================
// Low-level exports used only by disposable smoke tests
// ============================================================

// [[Rcpp::export]]
arma::vec gigg_bvar_test_rinvgamma_cpp(
    const int n,
    const double shape,
    const double scale) {

  Rcpp::RNGScope rng_scope;

  if (n < 1) {
    Rcpp::stop("n must be positive");
  }

  arma::vec draws(
    static_cast<arma::uword>(n)
  );

  for (int j = 0; j < n; ++j) {

    draws[static_cast<arma::uword>(j)] =
      gigg_bvar::rinvgamma(
        shape,
        scale
      );
  }

  return draws;
}


// [[Rcpp::export]]
arma::vec gigg_bvar_test_gamma_rate_cpp(
    const int n,
    const double shape,
    const double rate) {

  Rcpp::RNGScope rng_scope;

  if (n < 1) {
    Rcpp::stop("n must be positive");
  }

  arma::vec draws(
    static_cast<arma::uword>(n)
  );

  for (int j = 0; j < n; ++j) {

    draws[static_cast<arma::uword>(j)] =
      gigg_bvar::rgamma_rate(
        shape,
        rate
      );
  }

  return draws;
}


// [[Rcpp::export]]
Rcpp::List gigg_bvar_test_gig_cpp(
    const int n,
    const double p,
    const double chi,
    const double psi,
    const int maximum_attempts = 100000) {

  Rcpp::RNGScope rng_scope;

  if (n < 1) {
    Rcpp::stop("n must be positive");
  }

  arma::vec draws(
    static_cast<arma::uword>(n)
  );

  Rcpp::IntegerVector attempts(n);

  for (int j = 0; j < n; ++j) {

    const gigg_bvar::GIGDrawResult result =
      gigg_bvar::draw_gig(
        p,
        chi,
        psi,
        maximum_attempts
      );

    draws[static_cast<arma::uword>(j)] =
      result.value;

    attempts[j] =
      result.attempts;
  }

  return Rcpp::List::create(
    Rcpp::Named("draws") = draws,
    Rcpp::Named("attempts") = attempts
  );
}


// [[Rcpp::export]]
arma::mat gigg_bvar_test_beta_cpp(
    const int n,
    const arma::vec& y,
    const arma::mat& X,
    const arma::vec& omega,
    const double sigma2,
    const std::string beta_algorithm = "auto") {

  Rcpp::RNGScope rng_scope;

  if (
    n < 1 ||
    y.n_elem != X.n_rows ||
    omega.n_elem != X.n_cols ||
    !y.is_finite() ||
    !X.is_finite() ||
    !omega.is_finite() ||
    arma::any(omega <= 0.0) ||
    !std::isfinite(sigma2) ||
    sigma2 <= 0.0
  ) {
    Rcpp::stop("Invalid beta-test inputs");
  }

  const arma::mat Y =
    arma::reshape(
      y,
      X.n_rows,
      1u
    );

  const arma::mat XtX =
    X.t() * X;

  const arma::mat XtY =
    X.t() * Y;

  const int resolved =
    gigg_bvar::resolve_beta_algorithm(
      gigg_bvar::parse_beta_algorithm(
        beta_algorithm
      ),
      X.n_cols,
      X.n_rows
    );

  gigg_bvar::CholWorkspace chol_workspace(
    X.n_cols
  );

  gigg_bvar::BhattacharyaWorkspace
    bhattacharya_workspace(
      X.n_rows,
      X.n_cols
    );

  arma::mat draws(
    static_cast<arma::uword>(n),
    X.n_cols
  );

  arma::vec beta(X.n_cols);

  for (int j = 0; j < n; ++j) {

    gigg_bvar::draw_beta(
      beta,
      X,
      Y,
      XtX,
      XtY,
      0u,
      omega,
      sigma2,
      resolved,
      chol_workspace,
      bhattacharya_workspace
    );

    draws.row(
      static_cast<arma::uword>(j)
    ) = beta.t();
  }

  return draws;
}


// [[Rcpp::export]]
arma::vec gigg_bvar_test_sigma2_cpp(
    const int n,
    const arma::vec& residual,
    const arma::vec& beta,
    const arma::vec& omega,
    const double prior_shape,
    const double prior_scale) {

  Rcpp::RNGScope rng_scope;

  if (n < 1) {
    Rcpp::stop("n must be positive");
  }

  arma::vec draws(
    static_cast<arma::uword>(n)
  );

  for (int j = 0; j < n; ++j) {

    draws[static_cast<arma::uword>(j)] =
      gigg_bvar::draw_sigma2_from_residual(
        residual,
        beta,
        omega,
        prior_shape,
        prior_scale
      );
  }

  return draws;
}


// [[Rcpp::export]]
arma::vec gigg_bvar_test_tau2_cpp(
    const int n,
    const arma::mat& beta,
    const arma::vec& sigma2,
    const arma::vec& gamma2,
    const arma::mat& lambda2,
    const arma::imat& block_id,
    const double tau_df,
    const double psi_tau) {

  Rcpp::RNGScope rng_scope;

  if (n < 1) {
    Rcpp::stop("n must be positive");
  }

  const arma::imat block_id_zero =
    gigg_bvar::zero_based_block_map(
      block_id,
      static_cast<int>(gamma2.n_elem)
    );

  arma::vec draws(
    static_cast<arma::uword>(n)
  );

  for (int j = 0; j < n; ++j) {

    draws[static_cast<arma::uword>(j)] =
      gigg_bvar::draw_tau2(
        beta,
        sigma2,
        gamma2,
        lambda2,
        block_id_zero,
        tau_df,
        psi_tau
      );
  }

  return draws;
}


// [[Rcpp::export]]
arma::vec gigg_bvar_test_half_t_auxiliary_cpp(
    const int n,
    const double tau2,
    const double tau_df,
    const double tau_scale) {

  Rcpp::RNGScope rng_scope;

  if (n < 1) {
    Rcpp::stop("n must be positive");
  }

  arma::vec draws(
    static_cast<arma::uword>(n)
  );

  for (int j = 0; j < n; ++j) {

    draws[static_cast<arma::uword>(j)] =
      gigg_bvar::draw_half_t_auxiliary(
        tau2,
        tau_df,
        tau_scale
      );
  }

  return draws;
}


// [[Rcpp::export]]
arma::vec gigg_bvar_test_local_scale_cpp(
    const int n,
    const double coefficient,
    const double sigma2,
    const double tau2,
    const double gamma2,
    const double lambda_shape = 2.5,
    const double lambda_prior_scale = 1.0) {

  Rcpp::RNGScope rng_scope;

  if (n < 1) {
    Rcpp::stop("n must be positive");
  }

  const double conditional_shape =
    lambda_shape + 0.5;

  const double conditional_scale =
    lambda_prior_scale +
    coefficient *
    coefficient /
    (
      2.0 *
      sigma2 *
      tau2 *
      gamma2
    );

  arma::vec draws(
    static_cast<arma::uword>(n)
  );

  for (int j = 0; j < n; ++j) {

    draws[static_cast<arma::uword>(j)] =
      gigg_bvar::rinvgamma(
        conditional_shape,
        conditional_scale
      );
  }

  return draws;
}


// [[Rcpp::export]]
Rcpp::List gigg_bvar_test_block_asis_cpp(
    arma::vec gamma2,
    arma::mat lambda2,
    const arma::imat& block_id,
    const double gamma_shape = 0.5,
    const double gamma_rate = 1.0,
    const double lambda_shape = 2.5,
    const double lambda_prior_scale = 1.0,
    const double invariant_tolerance = 1e-10) {

  Rcpp::RNGScope rng_scope;

  const gigg_bvar::BlockMembershipCache cache =
    gigg_bvar::build_block_membership_cache(
      block_id,
      static_cast<int>(gamma2.n_elem)
    );

  arma::mat before(
    lambda2.n_rows,
    lambda2.n_cols
  );

  for (
    arma::uword equation = 0u;
    equation < lambda2.n_cols;
    ++equation
  ) {
    for (
      arma::uword row = 0u;
      row < lambda2.n_rows;
      ++row
    ) {

      before(row, equation) =
        gamma2[
          static_cast<arma::uword>(
            block_id(row, equation) - 1
          )
        ] *
        lambda2(row, equation);
    }
  }

  const arma::uword maximum_block_size =
    cache.coefficient_count.max();

  gigg_bvar::BlockCoefficientAsisWorkspace
    workspace(maximum_block_size);

  const gigg_bvar::BlockCoefficientAsisResult result =
    gigg_bvar::update_block_coefficient_asis(
      gamma2,
      lambda2,
      cache,
      gamma_shape,
      gamma_rate,
      lambda_shape,
      lambda_prior_scale,
      workspace,
      invariant_tolerance
    );

  arma::mat after(
    lambda2.n_rows,
    lambda2.n_cols
  );

  for (
    arma::uword equation = 0u;
    equation < lambda2.n_cols;
    ++equation
  ) {
    for (
      arma::uword row = 0u;
      row < lambda2.n_rows;
      ++row
    ) {

      after(row, equation) =
        gamma2[
          static_cast<arma::uword>(
            block_id(row, equation) - 1
          )
        ] *
        lambda2(row, equation);
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("gamma2") =
      gamma2,

    Rcpp::Named("lambda2") =
      lambda2,

    Rcpp::Named("invariant_before") =
      before,

    Rcpp::Named("invariant_after") =
      after,

    Rcpp::Named("block_draws_attempted") =
      static_cast<double>(
        result.block_draws_attempted
      ),

    Rcpp::Named("successful_block_draws") =
      static_cast<double>(
        result.successful_block_draws
      ),

    Rcpp::Named("nonfinite_failures") =
      static_cast<double>(
        result.nonfinite_failures
      ),

    Rcpp::Named("invariant_failures") =
      static_cast<double>(
        result.invariant_failures
      ),

    Rcpp::Named("absolute_log_movement_sum") =
      result.absolute_log_movement_sum,

    Rcpp::Named("maximum_absolute_log_movement") =
      result.maximum_absolute_log_movement,

    Rcpp::Named("maximum_relative_invariant_error") =
      result.maximum_relative_invariant_error
  );
}


// [[Rcpp::export]]
Rcpp::List gigg_bvar_test_global_asis_cpp(
    double tau2,
    arma::vec gamma2,
    const double psi_tau,
    const double tau_df = 5.0,
    const double gamma_shape = 0.5,
    const double gamma_rate = 1.0,
    const double invariant_tolerance = 1e-10) {

  Rcpp::RNGScope rng_scope;

  const arma::vec before =
    tau2 * gamma2;

  gigg_bvar::GlobalBlockAsisWorkspace
    workspace(gamma2.n_elem);

  const gigg_bvar::GlobalBlockAsisResult result =
    gigg_bvar::update_global_block_asis(
      tau2,
      gamma2,
      psi_tau,
      tau_df,
      gamma_shape,
      gamma_rate,
      workspace,
      invariant_tolerance
    );

  return Rcpp::List::create(
    Rcpp::Named("tau2") =
      tau2,

    Rcpp::Named("gamma2") =
      gamma2,

    Rcpp::Named("invariant_before") =
      before,

    Rcpp::Named("invariant_after") =
      tau2 * gamma2,

    Rcpp::Named("attempted") =
      result.attempted,

    Rcpp::Named("successful") =
      result.successful,

    Rcpp::Named("nonfinite_failures") =
      static_cast<double>(
        result.nonfinite_failures
      ),

    Rcpp::Named("invariant_failures") =
      static_cast<double>(
        result.invariant_failures
      ),

    Rcpp::Named("absolute_log_movement") =
      result.absolute_log_movement,

    Rcpp::Named("maximum_relative_invariant_error") =
      result.maximum_relative_invariant_error
  );
}


// ============================================================
// One-chain kernel
// ============================================================

// [[Rcpp::export]]
Rcpp::List gigg_bvar_chain_cpp(
    const Rcpp::List& input) {

  Rcpp::RNGScope rng_scope;

  const Clock::time_point total_start =
    Clock::now();

  const Rcpp::List data =
    input["data"];

  const Rcpp::List maps =
    input["maps"];

  const Rcpp::List prior =
    input["prior"];

  const Rcpp::List initial_state =
    input["state"];

  const Rcpp::List control =
    input["control"];

  const arma::mat Y =
    Rcpp::as<arma::mat>(
      data["Y"]
    );

  const arma::mat X =
    Rcpp::as<arma::mat>(
      data["X"]
    );

  const int T_p =
    Rcpp::as<int>(
      data["T_p"]
    );

  const int m =
    Rcpp::as<int>(
      data["m"]
    );

  const int N =
    Rcpp::as<int>(
      data["N"]
    );

  const int k =
    Rcpp::as<int>(
      data["k"]
    );

  const arma::imat block_id =
    Rcpp::as<arma::imat>(
      maps["block_id"]
    );

  const int n_blocks =
    Rcpp::as<int>(
      maps["n_blocks"]
    );

  const double sigma_prior_shape =
    Rcpp::as<double>(
      prior["sigma_a"]
    );

  const double sigma_prior_scale =
    Rcpp::as<double>(
      prior["sigma_b"]
    );

  const double tau_df =
    Rcpp::as<double>(
      prior["tau_df"]
    );

  const double tau_scale =
    Rcpp::as<double>(
      prior["tau_scale"]
    );

  const double gamma_shape =
    Rcpp::as<double>(
      prior["gamma_shape"]
    );

  const double gamma_rate =
    Rcpp::as<double>(
      prior["gamma_rate"]
    );

  const double lambda_shape =
    Rcpp::as<double>(
      prior["lambda_shape"]
    );

  const double lambda_prior_scale =
    Rcpp::as<double>(
      prior["lambda_scale"]
    );

  arma::mat beta =
    Rcpp::as<arma::mat>(
      initial_state["beta"]
    );

  arma::vec sigma2 =
    Rcpp::as<arma::vec>(
      initial_state["sigma2"]
    );

  double tau2 =
    Rcpp::as<double>(
      initial_state["tau2"]
    );

  double psi_tau =
    Rcpp::as<double>(
      initial_state["psi_tau"]
    );

  arma::vec gamma2 =
    Rcpp::as<arma::vec>(
      initial_state["gamma2"]
    );

  arma::mat lambda2 =
    Rcpp::as<arma::mat>(
      initial_state["lambda2"]
    );

  const int burnin =
    Rcpp::as<int>(
      control["burnin"]
    );

  const int n_draws =
    Rcpp::as<int>(
      control["draws"]
    );

  const int thin =
    Rcpp::as<int>(
      control["thin"]
    );

  const std::string requested_algorithm_name =
    Rcpp::as<std::string>(
      control["beta_algorithm"]
    );

  const bool use_asis =
    list_has(
      control,
      "use_asis"
    )
      ? Rcpp::as<bool>(
          control["use_asis"]
        )
      : true;

  const int asis_every =
    list_has(
      control,
      "asis_every"
    )
      ? Rcpp::as<int>(
          control["asis_every"]
        )
      : 1;

  const bool keep_all =
    list_has(
      control,
      "keep_all"
    )
      ? Rcpp::as<bool>(
          control["keep_all"]
        )
      : false;

  arma::imat monitor_beta_raw(
    0u,
    2u
  );

  arma::ivec monitor_gamma_raw(0u);

  arma::imat monitor_lambda_raw(
    0u,
    2u
  );

  if (
    list_has(
      control,
      "monitor_beta"
    )
  ) {
    monitor_beta_raw =
      Rcpp::as<arma::imat>(
        control["monitor_beta"]
      );
  }

  if (
    list_has(
      control,
      "monitor_gamma"
    )
  ) {
    monitor_gamma_raw =
      Rcpp::as<arma::ivec>(
        control["monitor_gamma"]
      );
  }

  if (
    list_has(
      control,
      "monitor_lambda"
    )
  ) {
    monitor_lambda_raw =
      Rcpp::as<arma::imat>(
        control["monitor_lambda"]
      );
  }

  if (
    T_p < 1 ||
    m < 1 ||
    N < 1 ||
    k < 1 ||
    n_blocks < 1 ||
    burnin < 0 ||
    n_draws < 1 ||
    thin < 1 ||
    asis_every < 1
  ) {
    Rcpp::stop(
      "Invalid dimensions or MCMC controls"
    );
  }

  const arma::uword T_p_u =
    static_cast<arma::uword>(T_p);

  const arma::uword N_u =
    static_cast<arma::uword>(N);

  const arma::uword k_u =
    static_cast<arma::uword>(k);

  const arma::uword n_blocks_u =
    static_cast<arma::uword>(n_blocks);

  const arma::uword draws_u =
    static_cast<arma::uword>(n_draws);

  if (
    Y.n_rows != T_p_u ||
    Y.n_cols != N_u ||
    X.n_rows != T_p_u ||
    X.n_cols != k_u ||
    beta.n_rows != k_u ||
    beta.n_cols != N_u ||
    sigma2.n_elem != N_u ||
    gamma2.n_elem != n_blocks_u ||
    lambda2.n_rows != k_u ||
    lambda2.n_cols != N_u ||
    block_id.n_rows != k_u ||
    block_id.n_cols != N_u ||
    N % m != 0 ||
    k % N != 0
  ) {
    Rcpp::stop(
      "Inconsistent GIGG data, state or map dimensions"
    );
  }

  if (
    !Y.is_finite() ||
    !X.is_finite() ||
    !beta.is_finite() ||
    !sigma2.is_finite() ||
    arma::any(sigma2 <= 0.0) ||
    !gamma2.is_finite() ||
    arma::any(gamma2 <= 0.0) ||
    !lambda2.is_finite() ||
    arma::any(
      arma::vectorise(lambda2) <= 0.0
    )
  ) {
    Rcpp::stop(
      "Initial GIGG state and data must be finite and valid"
    );
  }

  gigg_bvar::require_positive_finite(
    sigma_prior_shape,
    "sigma_a"
  );

  gigg_bvar::require_positive_finite(
    sigma_prior_scale,
    "sigma_b"
  );

  gigg_bvar::require_positive_finite(
    tau_df,
    "tau_df"
  );

  gigg_bvar::require_positive_finite(
    tau_scale,
    "tau_scale"
  );

  gigg_bvar::require_positive_finite(
    gamma_shape,
    "gamma_shape"
  );

  gigg_bvar::require_positive_finite(
    gamma_rate,
    "gamma_rate"
  );

  gigg_bvar::require_positive_finite(
    lambda_shape,
    "lambda_shape"
  );

  gigg_bvar::require_positive_finite(
    lambda_prior_scale,
    "lambda_scale"
  );

  gigg_bvar::require_positive_finite(
    tau2,
    "initial tau2"
  );

  gigg_bvar::require_positive_finite(
    psi_tau,
    "initial psi_tau"
  );

  const gigg_bvar::BlockMembershipCache block_cache =
    gigg_bvar::build_block_membership_cache(
      block_id,
      n_blocks,
      static_cast<arma::uword>(m * m)
    );

  const arma::imat block_id_zero =
    gigg_bvar::zero_based_block_map(
      block_id,
      n_blocks
    );

  const PairMonitorIndex monitor_beta =
    validate_pair_monitor(
      monitor_beta_raw,
      k_u,
      N_u,
      "monitor_beta"
    );

  const PairMonitorIndex monitor_lambda =
    validate_pair_monitor(
      monitor_lambda_raw,
      k_u,
      N_u,
      "monitor_lambda"
    );

  const arma::uvec monitor_gamma =
    validate_block_monitor(
      monitor_gamma_raw,
      n_blocks_u
    );

  const int requested_algorithm =
    gigg_bvar::parse_beta_algorithm(
      requested_algorithm_name
    );

  const int resolved_algorithm =
    gigg_bvar::resolve_beta_algorithm(
      requested_algorithm,
      k_u,
      T_p_u
    );

  const std::string resolved_algorithm_name =
    gigg_bvar::beta_algorithm_name(
      resolved_algorithm
    );

  const arma::mat XtX =
    X.t() * X;

  const arma::mat XtY =
    X.t() * Y;

  gigg_bvar::CholWorkspace chol_workspace(
    k_u
  );

  gigg_bvar::BhattacharyaWorkspace
    bhattacharya_workspace(
      T_p_u,
      k_u
    );

  gigg_bvar::BlockCoefficientAsisWorkspace
    block_asis_workspace(
      block_cache.coefficient_count.max()
    );

  gigg_bvar::GlobalBlockAsisWorkspace
    global_asis_workspace(
      n_blocks_u
    );

  arma::mat omega(
    k_u,
    N_u
  );

  arma::vec omega_equation(k_u);

  arma::mat residual(
    T_p_u,
    N_u
  );

  arma::mat beta_sum(
    k_u,
    N_u,
    arma::fill::zeros
  );

  arma::mat beta_sum_square(
    k_u,
    N_u,
    arma::fill::zeros
  );

  arma::vec sigma2_sum(
    N_u,
    arma::fill::zeros
  );

  double tau2_sum = 0.0;
  double psi_tau_sum = 0.0;

  arma::vec gamma2_sum(
    n_blocks_u,
    arma::fill::zeros
  );

  arma::mat lambda2_sum(
    k_u,
    N_u,
    arma::fill::zeros
  );

  arma::vec tau2_gamma2_sum(
    n_blocks_u,
    arma::fill::zeros
  );

  arma::mat effective_sum(
    k_u,
    N_u,
    arma::fill::zeros
  );

  arma::mat sigma2_draws(
    draws_u,
    N_u,
    arma::fill::zeros
  );

  arma::mat tau2_draws(
    draws_u,
    1u,
    arma::fill::zeros
  );

  arma::mat psi_tau_draws(
    draws_u,
    1u,
    arma::fill::zeros
  );

  arma::mat monitored_beta_draws(
    draws_u,
    monitor_beta.row.n_elem,
    arma::fill::zeros
  );

  arma::mat monitored_gamma2_draws(
    draws_u,
    monitor_gamma.n_elem,
    arma::fill::zeros
  );

  arma::mat monitored_tau2_gamma2_draws(
    draws_u,
    monitor_gamma.n_elem,
    arma::fill::zeros
  );

  arma::mat monitored_lambda2_draws(
    draws_u,
    monitor_lambda.row.n_elem,
    arma::fill::zeros
  );

  arma::mat monitored_effective_draws(
    draws_u,
    monitor_lambda.row.n_elem,
    arma::fill::zeros
  );

  arma::cube full_beta_draws;
  arma::mat full_gamma2_draws;
  arma::cube full_lambda2_draws;

  if (keep_all) {

    full_beta_draws.set_size(
      k_u,
      N_u,
      draws_u
    );

    full_gamma2_draws.set_size(
      n_blocks_u,
      draws_u
    );

    full_lambda2_draws.set_size(
      k_u,
      N_u,
      draws_u
    );
  }

  std::uint64_t gig_draws = 0u;
  std::uint64_t gig_attempts = 0u;
  std::uint64_t gig_chi_floor_count = 0u;

  int gig_max_attempts_one_draw = 0;

  double gig_min_raw_chi =
    std::numeric_limits<double>::infinity();

  double gig_max_raw_chi = 0.0;

  std::uint64_t block_asis_updates = 0u;
  std::uint64_t block_asis_successful_updates = 0u;
  std::uint64_t block_asis_draws_attempted = 0u;
  std::uint64_t block_asis_successful_draws = 0u;
  std::uint64_t block_asis_nonfinite_failures = 0u;
  std::uint64_t block_asis_invariant_failures = 0u;

  double block_asis_log_movement_sum = 0.0;
  double block_asis_max_log_movement = 0.0;
  double block_asis_max_invariant_error = 0.0;

  std::uint64_t global_asis_updates = 0u;
  std::uint64_t global_asis_successful_updates = 0u;
  std::uint64_t global_asis_nonfinite_failures = 0u;
  std::uint64_t global_asis_invariant_failures = 0u;

  double global_asis_log_movement_sum = 0.0;
  double global_asis_max_log_movement = 0.0;
  double global_asis_max_invariant_error = 0.0;

  const std::uint64_t burnin_u =
    static_cast<std::uint64_t>(burnin);

  const std::uint64_t draws_count_u =
    static_cast<std::uint64_t>(n_draws);

  const std::uint64_t thin_u =
    static_cast<std::uint64_t>(thin);

  const std::uint64_t total_iterations =
    checked_total_iterations(
      burnin,
      n_draws,
      thin
    );

  std::uint64_t retained = 0u;

  ChainTiming timing;

  for (
    std::uint64_t iteration = 1u;
    iteration <= total_iterations;
    ++iteration
  ) {

    Clock::time_point step_start =
      Clock::now();

    for (
      arma::uword equation = 0u;
      equation < N_u;
      ++equation
    ) {

      gigg_bvar::fill_omega(
        omega_equation,
        equation,
        tau2,
        gamma2,
        lambda2,
        block_id_zero
      );

      omega.col(equation) =
        omega_equation;
    }

    timing.omega_seconds +=
      elapsed_seconds(
        step_start,
        Clock::now()
      );

    // Step 1: beta.

    step_start =
      Clock::now();

    for (
      arma::uword equation = 0u;
      equation < N_u;
      ++equation
    ) {

      arma::vec beta_equation(
        beta.col(equation)
      );

      omega_equation =
        omega.col(equation);

      gigg_bvar::draw_beta(
        beta_equation,
        X,
        Y,
        XtX,
        XtY,
        equation,
        omega_equation,
        sigma2[equation],
        resolved_algorithm,
        chol_workspace,
        bhattacharya_workspace
      );

      beta.col(equation) =
        beta_equation;
    }

    timing.beta_seconds +=
      elapsed_seconds(
        step_start,
        Clock::now()
      );

    // Step 2: innovation variances.

    step_start =
      Clock::now();

    residual =
      Y - X * beta;

    for (
      arma::uword equation = 0u;
      equation < N_u;
      ++equation
    ) {

      omega_equation =
        omega.col(equation);

      sigma2[equation] =
        gigg_bvar::draw_sigma2_from_residual(
          residual.col(equation),
          beta.col(equation),
          omega_equation,
          sigma_prior_shape,
          sigma_prior_scale
        );
    }

    timing.sigma2_seconds +=
      elapsed_seconds(
        step_start,
        Clock::now()
      );

    // Steps 3-4: global half-t hierarchy.

    step_start =
      Clock::now();

    tau2 =
      gigg_bvar::draw_tau2(
        beta,
        sigma2,
        gamma2,
        lambda2,
        block_id_zero,
        tau_df,
        psi_tau
      );

    psi_tau =
      gigg_bvar::draw_half_t_auxiliary(
        tau2,
        tau_df,
        tau_scale
      );

    timing.tau_seconds +=
      elapsed_seconds(
        step_start,
        Clock::now()
      );

    // Step 5: ordinary block GIG update.

    step_start =
      Clock::now();

    const gigg_bvar::GIGUpdateSummary gig_summary =
      gigg_bvar::update_block_scales_gig(
        gamma2,
        beta,
        sigma2,
        lambda2,
        block_cache,
        tau2,
        gamma_shape,
        gamma_rate
      );

    gig_draws +=
      gig_summary.draws;

    gig_attempts +=
      gig_summary.total_attempts;

    gig_chi_floor_count +=
      gig_summary.chi_floor_count;

    gig_max_attempts_one_draw =
      std::max(
        gig_max_attempts_one_draw,
        gig_summary.maximum_attempts_for_one_draw
      );

    gig_min_raw_chi =
      std::min(
        gig_min_raw_chi,
        gig_summary.minimum_raw_chi
      );

    gig_max_raw_chi =
      std::max(
        gig_max_raw_chi,
        gig_summary.maximum_raw_chi
      );

    timing.gamma_seconds +=
      elapsed_seconds(
        step_start,
        Clock::now()
      );

    // Step 6: coefficient-local inverse-gamma update.

    step_start =
      Clock::now();

    gigg_bvar::update_local_scales(
      lambda2,
      beta,
      sigma2,
      gamma2,
      block_id_zero,
      tau2,
      lambda_shape,
      lambda_prior_scale
    );

    timing.lambda_seconds +=
      elapsed_seconds(
        step_start,
        Clock::now()
      );

    // Step 7: block-coefficient ASIS.

    if (
      use_asis &&
      iteration %
        static_cast<std::uint64_t>(
          asis_every
        ) ==
        0u
    ) {

      step_start =
        Clock::now();

      ++block_asis_updates;

      const gigg_bvar::BlockCoefficientAsisResult result =
        gigg_bvar::update_block_coefficient_asis(
          gamma2,
          lambda2,
          block_cache,
          gamma_shape,
          gamma_rate,
          lambda_shape,
          lambda_prior_scale,
          block_asis_workspace
        );

      block_asis_draws_attempted +=
        result.block_draws_attempted;

      block_asis_successful_draws +=
        result.successful_block_draws;

      block_asis_nonfinite_failures +=
        result.nonfinite_failures;

      block_asis_invariant_failures +=
        result.invariant_failures;

      block_asis_log_movement_sum +=
        result.absolute_log_movement_sum;

      block_asis_max_log_movement =
        std::max(
          block_asis_max_log_movement,
          result.maximum_absolute_log_movement
        );

      block_asis_max_invariant_error =
        std::max(
          block_asis_max_invariant_error,
          result.maximum_relative_invariant_error
        );

      if (
        result.successful_block_draws ==
          result.block_draws_attempted &&
        result.nonfinite_failures == 0u &&
        result.invariant_failures == 0u
      ) {
        ++block_asis_successful_updates;
      }

      timing.block_asis_seconds +=
        elapsed_seconds(
          step_start,
          Clock::now()
        );
    }

    // Step 8: global-block ASIS.

    if (
      use_asis &&
      iteration %
        static_cast<std::uint64_t>(
          asis_every
        ) ==
        0u
    ) {

      step_start =
        Clock::now();

      ++global_asis_updates;

      const gigg_bvar::GlobalBlockAsisResult result =
        gigg_bvar::update_global_block_asis(
          tau2,
          gamma2,
          psi_tau,
          tau_df,
          gamma_shape,
          gamma_rate,
          global_asis_workspace
        );

      global_asis_nonfinite_failures +=
        result.nonfinite_failures;

      global_asis_invariant_failures +=
        result.invariant_failures;

      global_asis_max_invariant_error =
        std::max(
          global_asis_max_invariant_error,
          result.maximum_relative_invariant_error
        );

      if (result.successful) {

        ++global_asis_successful_updates;

        global_asis_log_movement_sum +=
          result.absolute_log_movement;

        global_asis_max_log_movement =
          std::max(
            global_asis_max_log_movement,
            result.absolute_log_movement
          );
      }

      timing.global_asis_seconds +=
        elapsed_seconds(
          step_start,
          Clock::now()
        );
    }

    // Step 9: storage.

    if (
      iteration > burnin_u &&
      (
        iteration - burnin_u
      ) %
        thin_u ==
        0u
    ) {

      step_start =
        Clock::now();

      const arma::uword draw =
        static_cast<arma::uword>(retained);

      const arma::vec tau2_gamma2 =
        tau_gamma_product(
          tau2,
          gamma2
        );

      const arma::mat effective =
        effective_product(
          tau2,
          gamma2,
          lambda2,
          block_id_zero
        );

      beta_sum +=
        beta;

      beta_sum_square +=
        beta % beta;

      sigma2_sum +=
        sigma2;

      tau2_sum +=
        tau2;

      psi_tau_sum +=
        psi_tau;

      gamma2_sum +=
        gamma2;

      lambda2_sum +=
        lambda2;

      tau2_gamma2_sum +=
        tau2_gamma2;

      effective_sum +=
        effective;

      sigma2_draws.row(draw) =
        sigma2.t();

      tau2_draws(draw, 0u) =
        tau2;

      psi_tau_draws(draw, 0u) =
        psi_tau;

      for (
        arma::uword j = 0u;
        j < monitor_beta.row.n_elem;
        ++j
      ) {

        monitored_beta_draws(draw, j) =
          beta(
            monitor_beta.row[j],
            monitor_beta.equation[j]
          );
      }

      for (
        arma::uword j = 0u;
        j < monitor_gamma.n_elem;
        ++j
      ) {

        const arma::uword block =
          monitor_gamma[j];

        monitored_gamma2_draws(draw, j) =
          gamma2[block];

        monitored_tau2_gamma2_draws(draw, j) =
          tau2_gamma2[block];
      }

      for (
        arma::uword j = 0u;
        j < monitor_lambda.row.n_elem;
        ++j
      ) {

        const arma::uword row =
          monitor_lambda.row[j];

        const arma::uword equation =
          monitor_lambda.equation[j];

        monitored_lambda2_draws(draw, j) =
          lambda2(row, equation);

        monitored_effective_draws(draw, j) =
          effective(row, equation);
      }

      if (keep_all) {

        full_beta_draws.slice(draw) =
          beta;

        full_gamma2_draws.col(draw) =
          gamma2;

        full_lambda2_draws.slice(draw) =
          lambda2;
      }

      ++retained;

      timing.storage_seconds +=
        elapsed_seconds(
          step_start,
          Clock::now()
        );
    }

    if (iteration % 1000u == 0u) {
      Rcpp::checkUserInterrupt();
    }
  }

  if (retained != draws_count_u) {
    Rcpp::stop(
      "Retained draw count does not match control$draws"
    );
  }

  const double retained_double =
    static_cast<double>(retained);

  const arma::mat beta_mean =
    beta_sum / retained_double;

  const arma::mat beta_mean_square =
    beta_sum_square / retained_double;

  const arma::vec sigma2_mean =
    sigma2_sum / retained_double;

  const double tau2_mean =
    tau2_sum / retained_double;

  const double psi_tau_mean =
    psi_tau_sum / retained_double;

  const arma::vec gamma2_mean =
    gamma2_sum / retained_double;

  const arma::mat lambda2_mean =
    lambda2_sum / retained_double;

  const arma::vec tau2_gamma2_mean =
    tau2_gamma2_sum / retained_double;

  const arma::mat effective_mean =
    effective_sum / retained_double;

  const double block_asis_mean_log_movement =
    block_asis_successful_draws > 0u
      ? block_asis_log_movement_sum /
          static_cast<double>(
            block_asis_successful_draws
          )
      : NA_REAL;

  const double global_asis_mean_log_movement =
    global_asis_successful_updates > 0u
      ? global_asis_log_movement_sum /
          static_cast<double>(
            global_asis_successful_updates
          )
      : NA_REAL;

  Rcpp::List final_state =
    Rcpp::clone(initial_state);

  final_state["beta"] =
    beta;

  final_state["sigma2"] =
    sigma2;

  final_state["tau2"] =
    tau2;

  final_state["psi_tau"] =
    psi_tau;

  final_state["gamma2"] =
    gamma2;

  final_state["lambda2"] =
    lambda2;

  Rcpp::List full_draws =
    Rcpp::List::create(
      Rcpp::Named("beta") =
        R_NilValue,

      Rcpp::Named("gamma2") =
        R_NilValue,

      Rcpp::Named("lambda2") =
        R_NilValue
    );

  if (keep_all) {

    full_draws["beta"] =
      full_beta_draws;

    full_draws["gamma2"] =
      full_gamma2_draws;

    full_draws["lambda2"] =
      full_lambda2_draws;
  }

  const double total_seconds =
    elapsed_seconds(
      total_start,
      Clock::now()
    );

  const double accounted_seconds =
    timing.omega_seconds +
    timing.beta_seconds +
    timing.sigma2_seconds +
    timing.tau_seconds +
    timing.gamma_seconds +
    timing.lambda_seconds +
    timing.block_asis_seconds +
    timing.global_asis_seconds +
    timing.storage_seconds;

  std::string stage =
    "beta_sigma2_tau2_psi_tau_gamma2_lambda2";

  if (use_asis) {
    stage += "_block_asis_global_asis";
  }

  return Rcpp::List::create(

    Rcpp::Named("posterior") =
      Rcpp::List::create(

        Rcpp::Named("beta_mean") =
          beta_mean,

        Rcpp::Named("beta_mean_square") =
          beta_mean_square,

        Rcpp::Named("sigma2_mean") =
          sigma2_mean,

        Rcpp::Named("tau2_mean") =
          tau2_mean,

        Rcpp::Named("psi_tau_mean") =
          psi_tau_mean,

        Rcpp::Named("gamma2_mean") =
          gamma2_mean,

        Rcpp::Named("lambda2_mean") =
          lambda2_mean,

        Rcpp::Named("tau2_gamma2_mean") =
          tau2_gamma2_mean,

        Rcpp::Named("effective_variance_mean") =
          effective_mean
      ),

    Rcpp::Named("draws") =
      Rcpp::List::create(

        Rcpp::Named("sigma2") =
          sigma2_draws,

        Rcpp::Named("tau2") =
          tau2_draws,

        Rcpp::Named("psi_tau") =
          psi_tau_draws
      ),

    Rcpp::Named("monitor") =
      Rcpp::List::create(

        Rcpp::Named("beta_index") =
          monitor_beta_raw,

        Rcpp::Named("beta") =
          monitored_beta_draws,

        Rcpp::Named("gamma_index") =
          monitor_gamma_raw,

        Rcpp::Named("gamma2") =
          monitored_gamma2_draws,

        Rcpp::Named("tau2_gamma2") =
          monitored_tau2_gamma2_draws,

        Rcpp::Named("lambda_index") =
          monitor_lambda_raw,

        Rcpp::Named("lambda2") =
          monitored_lambda2_draws,

        Rcpp::Named("effective_variance") =
          monitored_effective_draws
      ),

    Rcpp::Named("full_draws") =
      full_draws,

    Rcpp::Named("final_state") =
      final_state,

    Rcpp::Named("gig") =
      Rcpp::List::create(

        Rcpp::Named("parameterization") =
          "x^(p-1) exp(-(psi*x + chi/x)/2)",

        Rcpp::Named("sampler") =
          "log_scale_three_piece_tangent_rejection",

        Rcpp::Named("chi_floor") =
          gigg_bvar::gig_chi_floor,

        Rcpp::Named("maximum_attempts") =
          gigg_bvar::gig_maximum_attempts,

        Rcpp::Named("draws") =
          static_cast<double>(gig_draws),

        Rcpp::Named("total_attempts") =
          static_cast<double>(gig_attempts),

        Rcpp::Named("mean_attempts") =
          gig_draws > 0u
            ? static_cast<double>(gig_attempts) /
                static_cast<double>(gig_draws)
            : NA_REAL,

        Rcpp::Named("maximum_attempts_one_draw") =
          gig_max_attempts_one_draw,

        Rcpp::Named("chi_floor_count") =
          static_cast<double>(
            gig_chi_floor_count
          ),

        Rcpp::Named("minimum_raw_chi") =
          gig_min_raw_chi,

        Rcpp::Named("maximum_raw_chi") =
          gig_max_raw_chi
      ),

    Rcpp::Named("asis") =
      Rcpp::List::create(

        Rcpp::Named("enabled") =
          use_asis,

        Rcpp::Named("every") =
          asis_every,

        Rcpp::Named("invariant_tolerance") =
          gigg_bvar::asis_invariant_tolerance,

        Rcpp::Named("block_coefficient") =
          Rcpp::List::create(

            Rcpp::Named("updates") =
              static_cast<double>(
                block_asis_updates
              ),

            Rcpp::Named("successful_updates") =
              static_cast<double>(
                block_asis_successful_updates
              ),

            Rcpp::Named("block_draws_attempted") =
              static_cast<double>(
                block_asis_draws_attempted
              ),

            Rcpp::Named("successful_block_draws") =
              static_cast<double>(
                block_asis_successful_draws
              ),

            Rcpp::Named("mean_absolute_log_movement_gamma2") =
              block_asis_mean_log_movement,

            Rcpp::Named("maximum_absolute_log_movement_gamma2") =
              block_asis_max_log_movement,

            Rcpp::Named("seconds") =
              timing.block_asis_seconds,

            Rcpp::Named("nonfinite_failures") =
              static_cast<double>(
                block_asis_nonfinite_failures
              ),

            Rcpp::Named("invariant_failures") =
              static_cast<double>(
                block_asis_invariant_failures
              ),

            Rcpp::Named("maximum_relative_invariant_error") =
              block_asis_max_invariant_error
          ),

        Rcpp::Named("global_block") =
          Rcpp::List::create(

            Rcpp::Named("updates") =
              static_cast<double>(
                global_asis_updates
              ),

            Rcpp::Named("successful_updates") =
              static_cast<double>(
                global_asis_successful_updates
              ),

            Rcpp::Named("mean_absolute_log_movement_tau2") =
              global_asis_mean_log_movement,

            Rcpp::Named("maximum_absolute_log_movement_tau2") =
              global_asis_max_log_movement,

            Rcpp::Named("seconds") =
              timing.global_asis_seconds,

            Rcpp::Named("nonfinite_failures") =
              static_cast<double>(
                global_asis_nonfinite_failures
              ),

            Rcpp::Named("invariant_failures") =
              static_cast<double>(
                global_asis_invariant_failures
              ),

            Rcpp::Named("maximum_relative_invariant_error") =
              global_asis_max_invariant_error
          )
      ),

    Rcpp::Named("timing") =
      Rcpp::List::create(

        Rcpp::Named("total_seconds") =
          total_seconds,

        Rcpp::Named("omega_seconds") =
          timing.omega_seconds,

        Rcpp::Named("beta_seconds") =
          timing.beta_seconds,

        Rcpp::Named("sigma2_seconds") =
          timing.sigma2_seconds,

        Rcpp::Named("tau_seconds") =
          timing.tau_seconds,

        Rcpp::Named("gamma_seconds") =
          timing.gamma_seconds,

        Rcpp::Named("lambda_seconds") =
          timing.lambda_seconds,

        Rcpp::Named("block_asis_seconds") =
          timing.block_asis_seconds,

        Rcpp::Named("global_asis_seconds") =
          timing.global_asis_seconds,

        Rcpp::Named("storage_seconds") =
          timing.storage_seconds,

        Rcpp::Named("other_seconds") =
          std::max(
            0.0,
            total_seconds - accounted_seconds
          )
      ),

    Rcpp::Named("info") =
      Rcpp::List::create(

        Rcpp::Named("model") =
          "gigg_bvar",

        Rcpp::Named("model_set") =
          "gigg_v1",

        Rcpp::Named("backend") =
          "RcppArmadillo",

        Rcpp::Named("engine") =
          "fast",

        Rcpp::Named("stage") =
          stage,

        Rcpp::Named("beta_algorithm_requested") =
          requested_algorithm_name,

        Rcpp::Named("beta_algorithm_resolved") =
          resolved_algorithm_name,

        Rcpp::Named("use_asis") =
          use_asis,

        Rcpp::Named("asis_every") =
          asis_every,

        Rcpp::Named("iterations") =
          static_cast<double>(
            total_iterations
          ),

        Rcpp::Named("burnin") =
          burnin,

        Rcpp::Named("draws") =
          n_draws,

        Rcpp::Named("thin") =
          thin,

        Rcpp::Named("retained") =
          static_cast<double>(
            retained
          ),

        Rcpp::Named("keep_all") =
          keep_all,

        Rcpp::Named("storage_mode") =
          keep_all
            ? "full_and_monitored"
            : "monitored",

        Rcpp::Named("n_monitor_beta") =
          static_cast<int>(
            monitor_beta.row.n_elem
          ),

        Rcpp::Named("n_monitor_gamma") =
          static_cast<int>(
            monitor_gamma.n_elem
          ),

        Rcpp::Named("n_monitor_lambda") =
          static_cast<int>(
            monitor_lambda.row.n_elem
          ),

        Rcpp::Named("n_blocks") =
          n_blocks,

        Rcpp::Named("coefficients") =
          static_cast<double>(
            beta.n_elem
          )
      )
  );
}