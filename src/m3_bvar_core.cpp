// Fast M3 BVAR one-chain engine.
//
// Supports exact q Gibbs at either scan position, joint q-beta
// transport-MH after beta, and optional group-c ASIS interweaving.

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>
#include "m3_bvar_updates.h"


// Internal helpers -----

namespace {
  using Clock = std::chrono::steady_clock;

  inline double elapsed_seconds(
      const Clock::time_point& start,
      const Clock::time_point& end) {

    return std::chrono::duration<double>(
      end - start
    ).count();
  }

  inline void require_positive_finite(
      const double value,
      const std::string& name) {

    if (
      !std::isfinite(value) ||
      value <= 0.0
    ) {
      Rcpp::stop(
        name + " must be finite and positive"
      );
    }
  }

  inline void require_probability(
      const double value,
      const std::string& name) {

    if (
      !std::isfinite(value) ||
      value < 0.0 ||
      value > 1.0
    ) {
      Rcpp::stop(
        name + " must lie between zero and one"
      );
    }
  }

  inline void require_finite_matrix(
      const arma::mat& value,
      const std::string& name) {

    if (!value.is_finite()) {
      Rcpp::stop(
        name + " must contain only finite values"
      );
    }
  }

  inline void require_positive_matrix(
      const arma::mat& value,
      const std::string& name) {

    if (
      !value.is_finite() ||
      arma::any(
        arma::vectorise(value) <= 0.0
      )
    ) {
      Rcpp::stop(
        name + " must contain only finite positive values"
      );
    }
  }

  inline void require_positive_vector(
      const arma::vec& value,
      const std::string& name) {

    if (
      !value.is_finite() ||
      arma::any(
        value <= 0.0
      )
    ) {
      Rcpp::stop(
        name + " must contain only finite positive values"
      );
    }
  }

  inline bool list_has(
      const Rcpp::List& object,
      const char* name) {

    return object.containsElementNamed(
      name
    );
  }

  struct BetaMonitorIndex {

    arma::uvec row;
    arma::uvec equation;
  };

  inline BetaMonitorIndex validate_beta_monitor(
      const arma::imat& raw,
      const int k,
      const int N) {

    if (raw.n_rows == 0u) {

      if (
        raw.n_cols != 0u &&
        raw.n_cols != 2u
      ) {
        Rcpp::stop(
          "Empty monitor_beta must have zero or two columns"
        );
      }

      BetaMonitorIndex empty;

      empty.row.set_size(
        0u
      );

      empty.equation.set_size(
        0u
      );

      return empty;
    }

    if (raw.n_cols != 2u) {
      Rcpp::stop(
        "monitor_beta must have columns row and equation"
      );
    }

    BetaMonitorIndex result;

    result.row.set_size(
      raw.n_rows
    );

    result.equation.set_size(
      raw.n_rows
    );

    for (
      arma::uword j = 0u;
      j < raw.n_rows;
      ++j
    ) {

      const int row =
        raw(j, 0u);

      const int equation =
        raw(j, 1u);

      if (
        row < 1 ||
        row > k ||
        equation < 1 ||
        equation > N
      ) {
        Rcpp::stop(
          "monitor_beta contains an index outside the coefficient matrix"
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
            "monitor_beta contains duplicated coefficient indices"
          );
        }
      }

      result.row[j] =
        static_cast<arma::uword>(
          row - 1
        );

      result.equation[j] =
        static_cast<arma::uword>(
          equation - 1
        );
    }

    return result;
  }

  inline arma::uvec validate_lambda_monitor(
      const arma::ivec& raw,
      const int n_blocks) {

    arma::uvec result(
      raw.n_elem
    );

    for (
      arma::uword j = 0u;
      j < raw.n_elem;
      ++j
    ) {

      const int block =
        raw[j];

      if (
        block < 1 ||
        block > n_blocks
      ) {
        Rcpp::stop(
          "monitor_lambda contains an index outside 1:n_blocks"
        );
      }

      for (
        arma::uword previous = 0u;
        previous < j;
        ++previous
      ) {

        if (
          raw[previous] ==
          block
        ) {
          Rcpp::stop(
            "monitor_lambda contains duplicated block indices"
          );
        }
      }

      result[j] =
        static_cast<arma::uword>(
          block - 1
        );
    }

    return result;
  }

  struct ChainTiming {

    double variance_refresh_seconds =
      0.0;

    double beta_seconds =
      0.0;

    double q_transport_seconds =
      0.0;

    double sigma2_seconds =
      0.0;

    double tau_seconds =
      0.0;

    double c_seconds =
      0.0;

    double c_asis_seconds =
      0.0;

    double lambda_seconds =
      0.0;

    double q_gibbs_seconds =
      0.0;

    double storage_seconds =
      0.0;
  };
}


// Conditional test wrappers -----

// [[Rcpp::export]]
arma::vec m3_bvar_test_rinvgamma(
    const int n,
    const double shape,
    const double scale) {

  Rcpp::RNGScope rng_scope;

  if (n < 1) {
    Rcpp::stop(
      "n must be positive"
    );
  }

  require_positive_finite(
    shape,
    "shape"
  );

  require_positive_finite(
    scale,
    "scale"
  );

  arma::vec draws(
    static_cast<arma::uword>(
      n
    )
  );

  for (
    int draw = 0;
    draw < n;
    ++draw
  ) {

    draws[
      static_cast<arma::uword>(
        draw
      )
    ] =
      m3_bvar::rinvgamma(
        shape,
        scale
      );
  }

  return draws;
}


// [[Rcpp::export]]
Rcpp::IntegerVector m3_bvar_test_log_categorical(
    const arma::vec& log_weights,
    const int n) {

  Rcpp::RNGScope rng_scope;

  if (n < 1) {
    Rcpp::stop(
      "n must be positive"
    );
  }

  Rcpp::IntegerVector draws(
    n
  );

  for (
    int draw = 0;
    draw < n;
    ++draw
  ) {

    draws[draw] =
      static_cast<int>(
        m3_bvar::sample_log_weights(
          log_weights
        )
      ) +
      1;
  }

  return draws;
}


// [[Rcpp::export]]
arma::mat m3_bvar_test_beta_chol(
    const arma::vec& y,
    const arma::mat& X,
    const arma::vec& omega,
    const double sigma2,
    const int n_draws = 1) {

  Rcpp::RNGScope rng_scope;

  if (n_draws < 1) {
    Rcpp::stop(
      "n_draws must be positive"
    );
  }

  arma::mat draws(
    X.n_cols,
    static_cast<arma::uword>(
      n_draws
    )
  );

  for (
    int draw = 0;
    draw < n_draws;
    ++draw
  ) {

    draws.col(
      static_cast<arma::uword>(
        draw
      )
    ) =
      m3_bvar::draw_beta_chol(
        y,
        X,
        omega,
        sigma2
      );
  }

  return draws;
}


// [[Rcpp::export]]
arma::mat m3_bvar_test_beta_bhattacharya(
    const arma::vec& y,
    const arma::mat& X,
    const arma::vec& omega,
    const double sigma2,
    const int n_draws = 1) {

  Rcpp::RNGScope rng_scope;

  if (n_draws < 1) {
    Rcpp::stop(
      "n_draws must be positive"
    );
  }

  arma::mat draws(
    X.n_cols,
    static_cast<arma::uword>(
      n_draws
    )
  );

  for (
    int draw = 0;
    draw < n_draws;
    ++draw
  ) {

    draws.col(
      static_cast<arma::uword>(
        draw
      )
    ) =
      m3_bvar::draw_beta_bhattacharya(
        y,
        X,
        omega,
        sigma2
      );
  }

  return draws;
}


// [[Rcpp::export]]
arma::vec m3_bvar_test_beta_dispatch(
    const arma::vec& y,
    const arma::mat& X,
    const arma::vec& omega,
    const double sigma2,
    const std::string algorithm = "auto") {

  Rcpp::RNGScope rng_scope;

  return m3_bvar::draw_beta(
    y,
    X,
    omega,
    sigma2,
    algorithm
  );
}


// [[Rcpp::export]]
arma::vec m3_bvar_test_sigma2(
    const arma::vec& y,
    const arma::mat& X,
    const arma::vec& beta,
    const arma::vec& omega,
    const double prior_shape,
    const double prior_scale,
    const int n_draws = 1) {

  Rcpp::RNGScope rng_scope;

  if (n_draws < 1) {
    Rcpp::stop(
      "n_draws must be positive"
    );
  }

  arma::vec draws(
    static_cast<arma::uword>(
      n_draws
    )
  );

  for (
    int draw = 0;
    draw < n_draws;
    ++draw
  ) {

    draws[
      static_cast<arma::uword>(
        draw
      )
    ] =
      m3_bvar::draw_sigma2(
        y,
        X,
        beta,
        omega,
        prior_shape,
        prior_scale
      );
  }

  return draws;
}


// [[Rcpp::export]]
arma::vec m3_bvar_test_tau2(
    const arma::mat& beta,
    const arma::vec& sigma2,
    const arma::mat& base_omega,
    const double tau_df,
    const double psi_tau,
    const int n_draws = 1) {

  Rcpp::RNGScope rng_scope;

  if (n_draws < 1) {
    Rcpp::stop(
      "n_draws must be positive"
    );
  }

  arma::vec draws(
    static_cast<arma::uword>(
      n_draws
    )
  );

  for (
    int draw = 0;
    draw < n_draws;
    ++draw
  ) {

    draws[
      static_cast<arma::uword>(
        draw
      )
    ] =
      m3_bvar::draw_tau2(
        beta,
        sigma2,
        base_omega,
        tau_df,
        psi_tau
      );
  }

  return draws;
}


// [[Rcpp::export]]
arma::vec m3_bvar_test_half_t_auxiliary(
    const double scale2,
    const double degrees_freedom,
    const double prior_scale,
    const int n_draws = 1) {

  Rcpp::RNGScope rng_scope;

  if (n_draws < 1) {
    Rcpp::stop(
      "n_draws must be positive"
    );
  }

  arma::vec draws(
    static_cast<arma::uword>(
      n_draws
    )
  );

  for (
    int draw = 0;
    draw < n_draws;
    ++draw
  ) {

    draws[
      static_cast<arma::uword>(
        draw
      )
    ] =
      m3_bvar::draw_half_t_auxiliary(
        scale2,
        degrees_freedom,
        prior_scale
      );
  }

  return draws;
}


// [[Rcpp::export]]
arma::vec m3_bvar_test_group_scale2(
    const arma::mat& beta,
    const arma::vec& sigma2,
    const arma::mat& base_omega,
    const arma::imat& group_id,
    const int target_group,
    const double degrees_freedom,
    const double auxiliary,
    const int n_draws = 1) {

  Rcpp::RNGScope rng_scope;

  if (n_draws < 1) {
    Rcpp::stop(
      "n_draws must be positive"
    );
  }

  arma::vec draws(
    static_cast<arma::uword>(
      n_draws
    )
  );

  for (
    int draw = 0;
    draw < n_draws;
    ++draw
  ) {

    draws[
      static_cast<arma::uword>(
        draw
      )
    ] =
      m3_bvar::draw_group_scale2(
        beta,
        sigma2,
        base_omega,
        group_id,
        target_group,
        degrees_freedom,
        auxiliary
      );
  }

  return draws;
}


// [[Rcpp::export]]
Rcpp::List m3_bvar_test_q_log_weights(
    const arma::mat& beta,
    const arma::vec& sigma2,
    const arma::mat& base_omega,
    const arma::imat& group_id,
    const arma::umat& same_var,
    const int target_group,
    const int m,
    const arma::vec& q_grid,
    const arma::vec& q_prob) {

  const m3_bvar::QLogWeightResult result =
    m3_bvar::compute_q_log_weights(
      beta,
      sigma2,
      base_omega,
      group_id,
      same_var,
      target_group,
      m,
      q_grid,
      q_prob
    );

  return Rcpp::List::create(

    Rcpp::Named("log_weights") =
      result.log_weights,

    Rcpp::Named("n_diagonal") =
      result.n_diagonal,

    Rcpp::Named("n_off_diagonal") =
      result.n_off_diagonal,

    Rcpp::Named("diagonal_quadratic") =
      result.diagonal_quadratic,

    Rcpp::Named("off_diagonal_quadratic") =
      result.off_diagonal_quadratic
  );
}


// [[Rcpp::export]]
Rcpp::IntegerVector m3_bvar_test_q_draw(
    const arma::mat& beta,
    const arma::vec& sigma2,
    const arma::mat& base_omega,
    const arma::imat& group_id,
    const arma::umat& same_var,
    const int target_group,
    const int m,
    const arma::vec& q_grid,
    const arma::vec& q_prob,
    const int n_draws = 1) {

  Rcpp::RNGScope rng_scope;

  if (n_draws < 1) {
    Rcpp::stop(
      "n_draws must be positive"
    );
  }

  Rcpp::IntegerVector draws(
    n_draws
  );

  for (
    int draw = 0;
    draw < n_draws;
    ++draw
  ) {

    draws[draw] =
      static_cast<int>(
        m3_bvar::draw_q_index(
          beta,
          sigma2,
          base_omega,
          group_id,
          same_var,
          target_group,
          m,
          q_grid,
          q_prob
        )
      ) +
      1;
  }

  return draws;
}


// One-chain kernel -----

// [[Rcpp::export]]
Rcpp::List m3_bvar_chain_cpp(
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

  const Rcpp::List state =
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

  const arma::imat lag_id =
    Rcpp::as<arma::imat>(
      maps["lag_id"]
    );

  const arma::imat gc_id =
    Rcpp::as<arma::imat>(
      maps["gc_id"]
    );

  const arma::imat gq_id =
    Rcpp::as<arma::imat>(
      maps["gq_id"]
    );

  const arma::imat block_id =
    Rcpp::as<arma::imat>(
      maps["block_id"]
    );

  const arma::umat same_var =
    Rcpp::as<arma::umat>(
      maps["same_var"]
    );

  const arma::vec phi2 =
    Rcpp::as<arma::vec>(
      maps["phi2"]
    );

  const int n_c =
    Rcpp::as<int>(
      maps["n_c"]
    );

  const int n_q =
    Rcpp::as<int>(
      maps["n_q"]
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

  const double c_df =
    Rcpp::as<double>(
      prior["c_df"]
    );

  const double c_scale =
    Rcpp::as<double>(
      prior["c_scale"]
    );

  const double lambda_df =
    Rcpp::as<double>(
      prior["lambda_df"]
    );

  const double lambda_scale =
    Rcpp::as<double>(
      prior["lambda_scale"]
    );

  const arma::vec q_grid =
    Rcpp::as<arma::vec>(
      prior["q_grid"]
    );

  const arma::vec q_prob =
    Rcpp::as<arma::vec>(
      prior["q_prob"]
    );

  arma::mat beta =
    Rcpp::as<arma::mat>(
      state["beta"]
    );

  arma::vec sigma2 =
    Rcpp::as<arma::vec>(
      state["sigma2"]
    );

  double tau2 =
    Rcpp::as<double>(
      state["tau2"]
    );

  double psi_tau =
    Rcpp::as<double>(
      state["psi_tau"]
    );

  arma::vec c2 =
    Rcpp::as<arma::vec>(
      state["c2"]
    );

  arma::vec psi_c =
    Rcpp::as<arma::vec>(
      state["psi_c"]
    );

  arma::vec lambda2 =
    Rcpp::as<arma::vec>(
      state["lambda2"]
    );

  arma::vec xi =
    Rcpp::as<arma::vec>(
      state["xi"]
    );

  arma::ivec q_index =
    Rcpp::as<arma::ivec>(
      state["q_index"]
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

  const std::string q_update =
    Rcpp::as<std::string>(
      control["q_update"]
    );

  const double q_transport_global_probability =
    Rcpp::as<double>(
      control["q_transport_global_probability"]
    );

  const bool use_c_asis =
    list_has(
      control,
      "use_c_asis"
    )
      ? Rcpp::as<bool>(
          control["use_c_asis"]
        )
      : false;

  const int c_asis_every =
    list_has(
      control,
      "c_asis_every"
    )
      ? Rcpp::as<int>(
          control["c_asis_every"]
        )
      : 1;

  const double c_asis_slice_width =
    list_has(
      control,
      "c_asis_slice_width"
    )
      ? Rcpp::as<double>(
          control["c_asis_slice_width"]
        )
      : 1.0;

  const int c_asis_maximum_step_out =
    list_has(
      control,
      "c_asis_maximum_step_out"
    )
      ? Rcpp::as<int>(
          control["c_asis_maximum_step_out"]
        )
      : 20;

  const bool keep_all =
    list_has(
      control,
      "keep_all"
    )
      ? Rcpp::as<bool>(
          control["keep_all"]
        )
      : false;

  if (
    list_has(
      control,
      "scale_transport"
    )
  ) {

    const std::string scale_transport =
      Rcpp::as<std::string>(
        control["scale_transport"]
      );

    if (
      scale_transport !=
      "none"
    ) {
      Rcpp::stop(
        "The fast M3 engine does not implement scale-ridge transport"
      );
    }
  }

  arma::imat monitor_beta_raw(
    0u,
    2u
  );

  arma::ivec monitor_lambda_raw(
    0u
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
      "monitor_lambda"
    )
  ) {

    monitor_lambda_raw =
      Rcpp::as<arma::ivec>(
        control["monitor_lambda"]
      );
  }

  if (
    T_p < 1 ||
    m < 1 ||
    N < 1 ||
    k < 1
  ) {
    Rcpp::stop(
      "T_p, m, N and k must be positive"
    );
  }

  if (
    Y.n_rows !=
      static_cast<arma::uword>(
        T_p
      ) ||
    Y.n_cols !=
      static_cast<arma::uword>(
        N
      )
  ) {
    Rcpp::stop(
      "Y dimensions do not match T_p x N"
    );
  }

  if (
    X.n_rows !=
      static_cast<arma::uword>(
        T_p
      ) ||
    X.n_cols !=
      static_cast<arma::uword>(
        k
      )
  ) {
    Rcpp::stop(
      "X dimensions do not match T_p x k"
    );
  }

  const arma::uword k_u =
    static_cast<arma::uword>(
      k
    );

  const arma::uword N_u =
    static_cast<arma::uword>(
      N
    );

  const arma::uword T_p_u =
    static_cast<arma::uword>(
      T_p
    );

  if (
    lag_id.n_rows != k_u ||
    lag_id.n_cols != N_u ||
    gc_id.n_rows != k_u ||
    gc_id.n_cols != N_u ||
    gq_id.n_rows != k_u ||
    gq_id.n_cols != N_u ||
    block_id.n_rows != k_u ||
    block_id.n_cols != N_u ||
    same_var.n_rows != k_u ||
    same_var.n_cols != N_u
  ) {
    Rcpp::stop(
      "Map dimensions do not match k x N"
    );
  }

  if (
    beta.n_rows != k_u ||
    beta.n_cols != N_u
  ) {
    Rcpp::stop(
      "Initial beta dimensions do not match k x N"
    );
  }

  require_finite_matrix(
    Y,
    "Y"
  );

  require_finite_matrix(
    X,
    "X"
  );

  require_finite_matrix(
    beta,
    "beta"
  );

  if (
    sigma2.n_elem != N_u
  ) {
    Rcpp::stop(
      "Initial sigma2 length does not match N"
    );
  }

  require_positive_vector(
    sigma2,
    "sigma2"
  );

  require_positive_finite(
    sigma_prior_shape,
    "sigma_a"
  );

  require_positive_finite(
    sigma_prior_scale,
    "sigma_b"
  );

  require_positive_finite(
    tau2,
    "tau2"
  );

  require_positive_finite(
    psi_tau,
    "psi_tau"
  );

  require_positive_finite(
    tau_df,
    "tau_df"
  );

  require_positive_finite(
    tau_scale,
    "tau_scale"
  );

  if (
    n_c < 1 ||
    c2.n_elem !=
      static_cast<arma::uword>(
        n_c
      ) ||
    psi_c.n_elem !=
      static_cast<arma::uword>(
        n_c
      )
  ) {
    Rcpp::stop(
      "Initial c hierarchy dimensions are invalid"
    );
  }

  require_positive_vector(
    c2,
    "c2"
  );

  require_positive_vector(
    psi_c,
    "psi_c"
  );

  if (
    std::abs(
      c2[0u] -
      1.0
    ) >
      1e-12 ||
    std::abs(
      psi_c[0u] -
      1.0
    ) >
      1e-12
  ) {
    Rcpp::stop(
      "Anchored c2[1] and psi_c[1] must equal one"
    );
  }

  require_positive_finite(
    c_df,
    "c_df"
  );

  require_positive_finite(
    c_scale,
    "c_scale"
  );

  if (
    n_blocks < 1 ||
    lambda2.n_elem !=
      static_cast<arma::uword>(
        n_blocks
      ) ||
    xi.n_elem !=
      static_cast<arma::uword>(
        n_blocks
      )
  ) {
    Rcpp::stop(
      "Initial lambda hierarchy dimensions are invalid"
    );
  }

  require_positive_vector(
    lambda2,
    "lambda2"
  );

  require_positive_vector(
    xi,
    "xi"
  );

  require_positive_finite(
    lambda_df,
    "lambda_df"
  );

  require_positive_finite(
    lambda_scale,
    "lambda_scale"
  );

  require_positive_vector(
    phi2,
    "phi2"
  );

  if (
    n_q < 1 ||
    q_index.n_elem !=
      static_cast<arma::uword>(
        n_q
      )
  ) {
    Rcpp::stop(
      "Initial q_index length does not match n_q"
    );
  }

  if (
    q_grid.n_elem < 1u ||
    q_prob.n_elem !=
      q_grid.n_elem ||
    !q_grid.is_finite() ||
    !q_prob.is_finite() ||
    arma::any(
      q_prob < 0.0
    ) ||
    arma::accu(
      q_prob
    ) <=
      0.0
  ) {
    Rcpp::stop(
      "Invalid q_grid or q_prob"
    );
  }

  for (
    int group = 0;
    group < n_q;
    ++group
  ) {

    if (
      q_index[
        static_cast<arma::uword>(
          group
        )
      ] <
        1 ||
      q_index[
        static_cast<arma::uword>(
          group
        )
      ] >
        static_cast<int>(
          q_grid.n_elem
        )
    ) {
      Rcpp::stop(
        "Initial q_index lies outside q_grid"
      );
    }

    if (
      q_prob[
        static_cast<arma::uword>(
          q_index[
            static_cast<arma::uword>(
              group
            )
          ] -
          1
        )
      ] <=
        0.0
    ) {
      Rcpp::stop(
        "Initial q state has zero prior probability"
      );
    }
  }

  if (
    burnin < 0 ||
    n_draws < 1 ||
    thin < 1
  ) {
    Rcpp::stop(
      "Invalid MCMC controls"
    );
  }

  if (
    q_update !=
      "transport_after_beta" &&
    q_update !=
      "gibbs_after_beta" &&
    q_update !=
      "gibbs_end"
  ) {
    Rcpp::stop(
      "q_update must be 'transport_after_beta', 'gibbs_after_beta' or 'gibbs_end'"
    );
  }

  require_probability(
    q_transport_global_probability,
    "q_transport_global_probability"
  );

  if (
    c_asis_every < 1
  ) {
    Rcpp::stop(
      "c_asis_every must be positive"
    );
  }

  require_positive_finite(
    c_asis_slice_width,
    "c_asis_slice_width"
  );

  if (
    c_asis_maximum_step_out < 1
  ) {
    Rcpp::stop(
      "c_asis_maximum_step_out must be positive"
    );
  }

  if (
    q_update ==
      "transport_after_beta" &&
    q_grid.n_elem <
      2u
  ) {
    Rcpp::stop(
      "transport_after_beta requires at least two q-grid values"
    );
  }

  arma::imat gc_id_zero(
    k_u,
    N_u
  );

  arma::imat gq_id_zero(
    k_u,
    N_u
  );

  arma::imat block_id_zero(
    k_u,
    N_u
  );

  arma::mat phi_omega(
    k_u,
    N_u
  );

  for (
    arma::uword equation = 0u;
    equation < N_u;
    ++equation
  ) {

    for (
      arma::uword row = 0u;
      row < k_u;
      ++row
    ) {

      const int lag =
        lag_id(
          row,
          equation
        ) -
        1;

      const int c_group =
        gc_id(
          row,
          equation
        ) -
        1;

      const int q_group =
        gq_id(
          row,
          equation
        ) -
        1;

      const int block =
        block_id(
          row,
          equation
        ) -
        1;

      if (
        lag < 0 ||
        lag >=
          static_cast<int>(
            phi2.n_elem
          )
      ) {
        Rcpp::stop(
          "lag_id lies outside phi2"
        );
      }

      if (
        c_group < 0 ||
        c_group >= n_c
      ) {
        Rcpp::stop(
          "gc_id lies outside 1:n_c"
        );
      }

      if (
        q_group < 0 ||
        q_group >= n_q
      ) {
        Rcpp::stop(
          "gq_id lies outside 1:n_q"
        );
      }

      if (
        block < 0 ||
        block >= n_blocks
      ) {
        Rcpp::stop(
          "block_id lies outside 1:n_blocks"
        );
      }

      if (
        same_var(
          row,
          equation
        ) !=
          0u &&
        same_var(
          row,
          equation
        ) !=
          1u
      ) {
        Rcpp::stop(
          "same_var must contain only zero or one"
        );
      }

      gc_id_zero(
        row,
        equation
      ) =
        c_group;

      gq_id_zero(
        row,
        equation
      ) =
        q_group;

      block_id_zero(
        row,
        equation
      ) =
        block;

      phi_omega(
        row,
        equation
      ) =
        phi2[
          static_cast<arma::uword>(
            lag
          )
        ];
    }
  }

  const m3_bvar::GroupMembershipCache c_cache =
    m3_bvar::build_group_membership_cache(
      gc_id,
      n_c,
      "gc_id"
    );

  std::vector<m3_bvar::CAsisGroupCache>
    c_asis_cache;

  if (
    use_c_asis &&
    n_c > 1
  ) {

    c_asis_cache =
      m3_bvar::build_c_asis_cache(
        gc_id,
        n_c
      );
  }

  const m3_bvar::GroupMembershipCache block_cache =
    m3_bvar::build_group_membership_cache(
      block_id,
      n_blocks,
      "block_id"
    );

  const std::vector<m3_bvar::QGibbsGroupCache>
    q_gibbs_cache =
      m3_bvar::build_q_gibbs_cache(
        gq_id,
        same_var,
        n_q
      );

  std::vector<m3_bvar::QGroupTransportCache>
    q_transport_cache;

  if (
    q_update ==
      "transport_after_beta"
  ) {

    q_transport_cache =
      m3_bvar::build_q_transport_cache(
        gq_id,
        same_var,
        n_q
      );
  }

  const BetaMonitorIndex monitor_beta =
    validate_beta_monitor(
      monitor_beta_raw,
      k,
      N
    );

  const arma::uvec monitor_lambda =
    validate_lambda_monitor(
      monitor_lambda_raw,
      n_blocks
    );

  const int requested_algorithm =
    m3_bvar::parse_beta_algorithm(
      requested_algorithm_name
    );

  const int resolved_algorithm =
    m3_bvar::resolve_beta_algorithm(
      requested_algorithm,
      k_u,
      T_p_u
    );

  const std::string resolved_algorithm_name =
    m3_bvar::beta_algorithm_name(
      resolved_algorithm
    );

  const arma::mat XtX =
    X.t() *
    X;

  const arma::mat XtY =
    X.t() *
    Y;

  m3_bvar::CholWorkspace chol_workspace(
    k_u
  );

  m3_bvar::BhattacharyaWorkspace
    bhattacharya_workspace(
      T_p_u,
      k_u
    );

  std::vector<m3_bvar::QTransportWorkspace>
    q_transport_workspace;

  if (
    q_update ==
      "transport_after_beta"
  ) {

    q_transport_workspace.reserve(
      static_cast<std::size_t>(
        n_q
      )
    );

    for (
      int group = 0;
      group < n_q;
      ++group
    ) {

      q_transport_workspace.emplace_back(
        T_p_u
      );
    }
  }

  m3_bvar::CAsisWorkspace c_asis_workspace;

  arma::vec c_quadratics(
    static_cast<arma::uword>(
      n_c
    ),
    arma::fill::zeros
  );

  arma::vec lambda_quadratics(
    static_cast<arma::uword>(
      n_blocks
    ),
    arma::fill::zeros
  );

  arma::mat q_omega(
    k_u,
    N_u
  );

  arma::mat tau_base_omega(
    k_u,
    N_u
  );

  arma::mat c_base_omega(
    k_u,
    N_u
  );

  arma::mat lambda_base_omega(
    k_u,
    N_u
  );

  arma::mat q_base_omega(
    k_u,
    N_u
  );

  arma::mat omega(
    k_u,
    N_u
  );

  arma::mat residual(
    T_p_u,
    N_u,
    arma::fill::zeros
  );

  auto refresh_q_and_tau_omega =
    [&]() {

      for (
        arma::uword equation = 0u;
        equation < N_u;
        ++equation
      ) {

        for (
          arma::uword row = 0u;
          row < k_u;
          ++row
        ) {

          const arma::uword c_group =
            static_cast<arma::uword>(
              gc_id_zero(
                row,
                equation
              )
            );

          const arma::uword q_group =
            static_cast<arma::uword>(
              gq_id_zero(
                row,
                equation
              )
            );

          const arma::uword block =
            static_cast<arma::uword>(
              block_id_zero(
                row,
                equation
              )
            );

          const arma::uword q_grid_index =
            static_cast<arma::uword>(
              q_index[
                q_group
              ] -
              1
            );

          const double q =
            q_grid[
              q_grid_index
            ];

          const double log_q_multiplier =
            same_var(
              row,
              equation
            ) ==
              1u
              ? (
                  static_cast<double>(
                    m -
                    1
                  ) /
                  static_cast<double>(
                    m
                  )
                ) *
                q
              : -q /
                static_cast<double>(
                  m
                );

          const double q_multiplier =
            std::exp(
              log_q_multiplier
            );

          q_omega(
            row,
            equation
          ) =
            q_multiplier;

          const double tau_base =
            c2[
              c_group
            ] *
            phi_omega(
              row,
              equation
            ) *
            lambda2[
              block
            ] *
            q_multiplier;

          tau_base_omega(
            row,
            equation
          ) =
            tau_base;

          omega(
            row,
            equation
          ) =
            tau2 *
            tau_base;
        }
      }
    };

  auto construct_c_base_omega =
    [&]() {

      for (
        arma::uword equation = 0u;
        equation < N_u;
        ++equation
      ) {

        for (
          arma::uword row = 0u;
          row < k_u;
          ++row
        ) {

          const arma::uword block =
            static_cast<arma::uword>(
              block_id_zero(
                row,
                equation
              )
            );

          c_base_omega(
            row,
            equation
          ) =
            tau2 *
            phi_omega(
              row,
              equation
            ) *
            lambda2[
              block
            ] *
            q_omega(
              row,
              equation
            );
        }
      }
    };

  auto construct_lambda_base_omega =
    [&]() {

      for (
        arma::uword equation = 0u;
        equation < N_u;
        ++equation
      ) {

        for (
          arma::uword row = 0u;
          row < k_u;
          ++row
        ) {

          const arma::uword c_group =
            static_cast<arma::uword>(
              gc_id_zero(
                row,
                equation
              )
            );

          lambda_base_omega(
            row,
            equation
          ) =
            tau2 *
            c2[
              c_group
            ] *
            phi_omega(
              row,
              equation
            ) *
            q_omega(
              row,
              equation
            );
        }
      }
    };

  auto construct_q_base_omega =
    [&]() {

      for (
        arma::uword equation = 0u;
        equation < N_u;
        ++equation
      ) {

        for (
          arma::uword row = 0u;
          row < k_u;
          ++row
        ) {

          const arma::uword c_group =
            static_cast<arma::uword>(
              gc_id_zero(
                row,
                equation
              )
            );

          const arma::uword block =
            static_cast<arma::uword>(
              block_id_zero(
                row,
                equation
              )
            );

          q_base_omega(
            row,
            equation
          ) =
            tau2 *
            c2[
              c_group
            ] *
            phi_omega(
              row,
              equation
            ) *
            lambda2[
              block
            ];
        }
      }
    };

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

  double tau2_sum =
    0.0;

  double psi_tau_sum =
    0.0;

  arma::vec c2_sum(
    static_cast<arma::uword>(
      n_c
    ),
    arma::fill::zeros
  );

  arma::vec psi_c_sum(
    static_cast<arma::uword>(
      n_c
    ),
    arma::fill::zeros
  );

  arma::vec lambda2_sum(
    static_cast<arma::uword>(
      n_blocks
    ),
    arma::fill::zeros
  );

  arma::vec xi_sum(
    static_cast<arma::uword>(
      n_blocks
    ),
    arma::fill::zeros
  );

  arma::vec q_sum(
    static_cast<arma::uword>(
      n_q
    ),
    arma::fill::zeros
  );

  arma::vec q_index_sum(
    static_cast<arma::uword>(
      n_q
    ),
    arma::fill::zeros
  );

  const arma::uword draws_u =
    static_cast<arma::uword>(
      n_draws
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

  arma::mat c2_draws(
    draws_u,
    static_cast<arma::uword>(
      n_c
    ),
    arma::fill::zeros
  );

  arma::mat psi_c_draws(
    draws_u,
    static_cast<arma::uword>(
      n_c
    ),
    arma::fill::zeros
  );

  arma::mat q_draws(
    draws_u,
    static_cast<arma::uword>(
      n_q
    ),
    arma::fill::zeros
  );

  arma::imat q_index_draws(
    draws_u,
    static_cast<arma::uword>(
      n_q
    ),
    arma::fill::zeros
  );

  arma::mat monitored_beta_draws(
    draws_u,
    monitor_beta.row.n_elem,
    arma::fill::zeros
  );

  arma::mat monitored_lambda2_draws(
    draws_u,
    monitor_lambda.n_elem,
    arma::fill::zeros
  );

  arma::mat monitored_xi_draws(
    draws_u,
    monitor_lambda.n_elem,
    arma::fill::zeros
  );

  arma::cube full_beta_draws;

  arma::mat full_lambda2_draws;

  arma::mat full_xi_draws;

  if (keep_all) {

    full_beta_draws.set_size(
      k_u,
      N_u,
      draws_u
    );

    full_lambda2_draws.set_size(
      static_cast<arma::uword>(
        n_blocks
      ),
      draws_u
    );

    full_xi_draws.set_size(
      static_cast<arma::uword>(
        n_blocks
      ),
      draws_u
    );
  }

  Rcpp::IntegerVector q_transport_proposals(
    n_q
  );

  Rcpp::IntegerVector q_transport_acceptances(
    n_q
  );

  Rcpp::IntegerVector q_transport_local_proposals(
    n_q
  );

  Rcpp::IntegerVector q_transport_local_acceptances(
    n_q
  );

  Rcpp::IntegerVector q_transport_global_proposals(
    n_q
  );

  Rcpp::IntegerVector q_transport_global_acceptances(
    n_q
  );

  Rcpp::NumericVector q_transport_proposed_jump_sum(
    n_q
  );

  Rcpp::NumericVector q_transport_accepted_jump_sum(
    n_q
  );

  Rcpp::NumericVector q_transport_last_log_likelihood_ratio(
    n_q,
    NA_REAL
  );

  Rcpp::NumericVector q_transport_last_log_q_prior_ratio(
    n_q,
    NA_REAL
  );

  Rcpp::NumericVector q_transport_last_log_proposal_ratio(
    n_q,
    NA_REAL
  );

  Rcpp::NumericVector q_transport_last_log_acceptance_ratio(
    n_q,
    NA_REAL
  );

  Rcpp::NumericVector q_transport_last_acceptance_probability(
    n_q,
    NA_REAL
  );

  Rcpp::IntegerVector c_asis_updates(
    n_c
  );

  Rcpp::IntegerVector c_asis_moves(
    n_c
  );

  Rcpp::NumericVector c_asis_absolute_log_move_sum(
    n_c
  );

  Rcpp::NumericVector c_asis_evaluation_sum(
    n_c
  );

  Rcpp::NumericVector c_asis_left_expansion_sum(
    n_c
  );

  Rcpp::NumericVector c_asis_right_expansion_sum(
    n_c
  );

  Rcpp::NumericVector c_asis_shrink_step_sum(
    n_c
  );

  Rcpp::LogicalVector c_asis_last_moved(
    n_c,
    NA_LOGICAL
  );

  Rcpp::NumericVector c_asis_last_current_scale2(
    n_c,
    NA_REAL
  );

  Rcpp::NumericVector c_asis_last_new_scale2(
    n_c,
    NA_REAL
  );

  Rcpp::NumericVector c_asis_last_transport_scale(
    n_c,
    NA_REAL
  );

  Rcpp::NumericVector c_asis_last_log_density(
    n_c,
    NA_REAL
  );

  Rcpp::NumericVector c_asis_last_log_slice_height(
    n_c,
    NA_REAL
  );

  const std::uint64_t burnin_u =
    static_cast<std::uint64_t>(
      burnin
    );

  const std::uint64_t draws_count_u =
    static_cast<std::uint64_t>(
      n_draws
    );

  const std::uint64_t thin_u =
    static_cast<std::uint64_t>(
      thin
    );

  if (
    draws_count_u >
    (
      std::numeric_limits<std::uint64_t>::max() -
      burnin_u
    ) /
    thin_u
  ) {
    Rcpp::stop(
      "Requested MCMC iteration count overflows uint64"
    );
  }

  const std::uint64_t total_iterations =
    burnin_u +
    draws_count_u *
    thin_u;

  std::uint64_t retained =
    0u;

  ChainTiming timing;

  for (
    std::uint64_t iteration = 1u;
    iteration <= total_iterations;
    ++iteration
  ) {

    Clock::time_point step_start =
      Clock::now();

    refresh_q_and_tau_omega();

    timing.variance_refresh_seconds +=
      elapsed_seconds(
        step_start,
        Clock::now()
      );


    // Beta -----

    step_start =
      Clock::now();

    for (
      arma::uword equation = 0u;
      equation < N_u;
      ++equation
    ) {

      arma::vec beta_equation(
        beta.col(
          equation
        )
      );

      m3_bvar::draw_beta(
        beta_equation,
        X,
        Y,
        XtX,
        XtY,
        equation,
        omega.col(
          equation
        ),
        sigma2[
          equation
        ],
        resolved_algorithm,
        chol_workspace,
        bhattacharya_workspace
      );

      beta.col(
        equation
      ) =
        beta_equation;
    }

    timing.beta_seconds +=
      elapsed_seconds(
        step_start,
        Clock::now()
      );


    // q transport after beta -----

    if (
      q_update ==
        "transport_after_beta"
    ) {

      step_start =
        Clock::now();

      residual =
        Y -
        X *
        beta;

      for (
        int group = 0;
        group < n_q;
        ++group
      ) {

        const m3_bvar::QTransportMHResult result =
          m3_bvar::update_q_transport_mh(
            beta,
            residual,
            X,
            sigma2,
            q_transport_cache[
              static_cast<std::size_t>(
                group
              )
            ],
            static_cast<arma::uword>(
              q_index[
                static_cast<arma::uword>(
                  group
                )
              ] -
              1
            ),
            m,
            q_grid,
            q_prob,
            q_transport_global_probability,
            q_transport_workspace[
              static_cast<std::size_t>(
                group
              )
            ]
          );

        q_transport_proposals[
          group
        ] +=
          1;

        q_transport_proposed_jump_sum[
          group
        ] +=
          static_cast<double>(
            result.absolute_jump
          );

        if (result.used_global) {

          q_transport_global_proposals[
            group
          ] +=
            1;

        } else {

          q_transport_local_proposals[
            group
          ] +=
            1;
        }

        if (result.accepted) {

          q_transport_acceptances[
            group
          ] +=
            1;

          q_transport_accepted_jump_sum[
            group
          ] +=
            static_cast<double>(
              result.absolute_jump
            );

          if (result.used_global) {

            q_transport_global_acceptances[
              group
            ] +=
              1;

          } else {

            q_transport_local_acceptances[
              group
            ] +=
              1;
          }
        }

        q_transport_last_log_likelihood_ratio[
          group
        ] =
          result.log_likelihood_ratio;

        q_transport_last_log_q_prior_ratio[
          group
        ] =
          result.log_q_prior_ratio;

        q_transport_last_log_proposal_ratio[
          group
        ] =
          result.log_proposal_ratio;

        q_transport_last_log_acceptance_ratio[
          group
        ] =
          result.log_acceptance_ratio;

        q_transport_last_acceptance_probability[
          group
        ] =
          result.acceptance_probability;

        q_index[
          static_cast<arma::uword>(
            group
          )
        ] =
          static_cast<int>(
            result.new_index
          ) +
          1;
      }

      timing.q_transport_seconds +=
        elapsed_seconds(
          step_start,
          Clock::now()
        );

      step_start =
        Clock::now();

      refresh_q_and_tau_omega();

      timing.variance_refresh_seconds +=
        elapsed_seconds(
          step_start,
          Clock::now()
        );
    }


    // Exact q Gibbs after beta -----

    if (
      q_update ==
        "gibbs_after_beta"
    ) {

      step_start =
        Clock::now();

      construct_q_base_omega();

      for (
        int group = 0;
        group < n_q;
        ++group
      ) {

        q_index[
          static_cast<arma::uword>(
            group
          )
        ] =
          static_cast<int>(
            m3_bvar::draw_q_index_cached(
              beta,
              sigma2,
              q_base_omega,
              q_gibbs_cache[
                static_cast<std::size_t>(
                  group
                )
              ],
              m,
              q_grid,
              q_prob
            )
          ) +
          1;
      }

      timing.q_gibbs_seconds +=
        elapsed_seconds(
          step_start,
          Clock::now()
        );

      step_start =
        Clock::now();

      refresh_q_and_tau_omega();

      timing.variance_refresh_seconds +=
        elapsed_seconds(
          step_start,
          Clock::now()
        );
    }


    // Sigma2 -----

    step_start =
      Clock::now();

    if (
      q_update !=
        "transport_after_beta"
    ) {

      residual =
        Y -
        X *
        beta;
    }

    for (
      arma::uword equation = 0u;
      equation < N_u;
      ++equation
    ) {

      sigma2[
        equation
      ] =
        m3_bvar::draw_sigma2_from_residual(
          residual.col(
            equation
          ),
          beta.col(
            equation
          ),
          omega.col(
            equation
          ),
          sigma_prior_shape,
          sigma_prior_scale
        );
    }

    timing.sigma2_seconds +=
      elapsed_seconds(
        step_start,
        Clock::now()
      );


    // Tau hierarchy -----

    step_start =
      Clock::now();

    tau2 =
      m3_bvar::draw_tau2(
        beta,
        sigma2,
        tau_base_omega,
        tau_df,
        psi_tau
      );

    psi_tau =
      m3_bvar::draw_half_t_auxiliary(
        tau2,
        tau_df,
        tau_scale
      );

    timing.tau_seconds +=
      elapsed_seconds(
        step_start,
        Clock::now()
      );


    // c hierarchy -----

    step_start =
      Clock::now();

    construct_c_base_omega();

    if (n_c > 1) {

      m3_bvar::draw_all_group_scales2(
        c2,
        c_quadratics,
        beta,
        sigma2,
        c_base_omega,
        gc_id_zero,
        c_cache.coefficient_count,
        c_df,
        psi_c,
        1u
      );

      for (
        arma::uword group = 1u;
        group <
          static_cast<arma::uword>(
            n_c
          );
        ++group
      ) {

        psi_c[
          group
        ] =
          m3_bvar::draw_half_t_auxiliary(
            c2[
              group
            ],
            c_df,
            c_scale
          );
      }
    }

    c2[0u] =
      1.0;

    psi_c[0u] =
      1.0;

    timing.c_seconds +=
      elapsed_seconds(
        step_start,
        Clock::now()
      );


    // c-ASIS -----

    if (
      use_c_asis &&
      n_c > 1 &&
      iteration %
        static_cast<std::uint64_t>(
          c_asis_every
        ) ==
        0u
    ) {

      step_start =
        Clock::now();

      for (
        arma::uword group = 1u;
        group <
          static_cast<arma::uword>(
            n_c
          );
        ++group
      ) {

        const m3_bvar::CAsisSliceResult result =
          m3_bvar::update_c_group_asis(
            beta,
            residual,
            X,
            sigma2,
            c2[
              group
            ],
            psi_c[
              group
            ],
            c_df,
            c_asis_cache[
              static_cast<std::size_t>(
                group
              )
            ],
            c_asis_slice_width,
            c_asis_maximum_step_out,
            c_asis_workspace
          );

        const int output_group =
          static_cast<int>(
            group
          );

        c_asis_updates[
          output_group
        ] +=
          1;

        if (
          result.moved
        ) {

          c_asis_moves[
            output_group
          ] +=
            1;
        }

        c_asis_absolute_log_move_sum[
          output_group
        ] +=
          std::abs(
            result.new_log_scale2 -
            result.current_log_scale2
          );

        c_asis_evaluation_sum[
          output_group
        ] +=
          static_cast<double>(
            result.evaluations
          );

        c_asis_left_expansion_sum[
          output_group
        ] +=
          static_cast<double>(
            result.left_expansions
          );

        c_asis_right_expansion_sum[
          output_group
        ] +=
          static_cast<double>(
            result.right_expansions
          );

        c_asis_shrink_step_sum[
          output_group
        ] +=
          static_cast<double>(
            result.shrink_steps
          );

        c_asis_last_moved[
          output_group
        ] =
          result.moved;

        c_asis_last_current_scale2[
          output_group
        ] =
          result.current_scale2;

        c_asis_last_new_scale2[
          output_group
        ] =
          result.new_scale2;

        c_asis_last_transport_scale[
          output_group
        ] =
          result.transport_scale;

        c_asis_last_log_density[
          output_group
        ] =
          result.new_log_density;

        c_asis_last_log_slice_height[
          output_group
        ] =
          result.log_slice_height;
      }

      c2[0u] =
        1.0;

      psi_c[0u] =
        1.0;

      timing.c_asis_seconds +=
        elapsed_seconds(
          step_start,
          Clock::now()
        );
    }


    // Lambda hierarchy -----

    step_start =
      Clock::now();

    construct_lambda_base_omega();

    m3_bvar::draw_all_group_scales2(
      lambda2,
      lambda_quadratics,
      beta,
      sigma2,
      lambda_base_omega,
      block_id_zero,
      block_cache.coefficient_count,
      lambda_df,
      xi,
      0u
    );

    for (
      arma::uword block = 0u;
      block <
        static_cast<arma::uword>(
          n_blocks
        );
      ++block
    ) {

      xi[
        block
      ] =
        m3_bvar::draw_half_t_auxiliary(
          lambda2[
            block
          ],
          lambda_df,
          lambda_scale
        );
    }

    timing.lambda_seconds +=
      elapsed_seconds(
        step_start,
        Clock::now()
      );


    // Exact q Gibbs at end -----

    if (
      q_update ==
        "gibbs_end"
    ) {

      step_start =
        Clock::now();

      construct_q_base_omega();

      for (
        int group = 0;
        group < n_q;
        ++group
      ) {

        q_index[
          static_cast<arma::uword>(
            group
          )
        ] =
          static_cast<int>(
            m3_bvar::draw_q_index_cached(
              beta,
              sigma2,
              q_base_omega,
              q_gibbs_cache[
                static_cast<std::size_t>(
                  group
                )
              ],
              m,
              q_grid,
              q_prob
            )
          ) +
          1;
      }

      timing.q_gibbs_seconds +=
        elapsed_seconds(
          step_start,
          Clock::now()
        );
    }


    // Storage -----

    if (
      iteration >
        burnin_u &&
      (
        iteration -
        burnin_u
      ) %
        thin_u ==
        0u
    ) {

      step_start =
        Clock::now();

      if (
        retained >=
        draws_count_u
      ) {
        Rcpp::stop(
          "Internal retained-draw index exceeds control$draws"
        );
      }

      const arma::uword retained_index =
        static_cast<arma::uword>(
          retained
        );

      beta_sum +=
        beta;

      beta_sum_square +=
        beta %
        beta;

      sigma2_sum +=
        sigma2;

      tau2_sum +=
        tau2;

      psi_tau_sum +=
        psi_tau;

      c2_sum +=
        c2;

      psi_c_sum +=
        psi_c;

      lambda2_sum +=
        lambda2;

      xi_sum +=
        xi;

      sigma2_draws.row(
        retained_index
      ) =
        sigma2.t();

      tau2_draws(
        retained_index,
        0u
      ) =
        tau2;

      psi_tau_draws(
        retained_index,
        0u
      ) =
        psi_tau;

      c2_draws.row(
        retained_index
      ) =
        c2.t();

      psi_c_draws.row(
        retained_index
      ) =
        psi_c.t();

      for (
        arma::uword group = 0u;
        group <
          static_cast<arma::uword>(
            n_q
          );
        ++group
      ) {

        const int index =
          q_index[
            group
          ];

        const double q =
          q_grid[
            static_cast<arma::uword>(
              index -
              1
            )
          ];

        q_index_draws(
          retained_index,
          group
        ) =
          index;

        q_draws(
          retained_index,
          group
        ) =
          q;

        q_sum[
          group
        ] +=
          q;

        q_index_sum[
          group
        ] +=
          static_cast<double>(
            index
          );
      }

      for (
        arma::uword monitor = 0u;
        monitor <
          monitor_beta.row.n_elem;
        ++monitor
      ) {

        monitored_beta_draws(
          retained_index,
          monitor
        ) =
          beta(
            monitor_beta.row[
              monitor
            ],
            monitor_beta.equation[
              monitor
            ]
          );
      }

      for (
        arma::uword monitor = 0u;
        monitor <
          monitor_lambda.n_elem;
        ++monitor
      ) {

        const arma::uword block =
          monitor_lambda[
            monitor
          ];

        monitored_lambda2_draws(
          retained_index,
          monitor
        ) =
          lambda2[
            block
          ];

        monitored_xi_draws(
          retained_index,
          monitor
        ) =
          xi[
            block
          ];
      }

      if (keep_all) {

        full_beta_draws.slice(
          retained_index
        ) =
          beta;

        full_lambda2_draws.col(
          retained_index
        ) =
          lambda2;

        full_xi_draws.col(
          retained_index
        ) =
          xi;
      }

      retained +=
        1u;

      timing.storage_seconds +=
        elapsed_seconds(
          step_start,
          Clock::now()
        );
    }

    if (
      iteration %
        1000u ==
        0u
    ) {

      Rcpp::checkUserInterrupt();

      if (
        !beta.is_finite() ||
        !sigma2.is_finite() ||
        !c2.is_finite() ||
        !psi_c.is_finite() ||
        !lambda2.is_finite() ||
        !xi.is_finite() ||
        arma::any(
          sigma2 <= 0.0
        ) ||
        arma::any(
          c2 <= 0.0
        ) ||
        arma::any(
          psi_c <= 0.0
        ) ||
        arma::any(
          lambda2 <= 0.0
        ) ||
        arma::any(
          xi <= 0.0
        ) ||
        !std::isfinite(
          tau2
        ) ||
        tau2 <= 0.0 ||
        !std::isfinite(
          psi_tau
        ) ||
        psi_tau <= 0.0
      ) {
        Rcpp::stop(
          "Non-finite or non-positive state detected during MCMC"
        );
      }
    }
  }

  if (
    retained !=
      draws_count_u
  ) {
    Rcpp::stop(
      "The number of retained draws does not match control$draws"
    );
  }

  require_finite_matrix(
    beta_sum,
    "beta_sum"
  );

  require_finite_matrix(
    beta_sum_square,
    "beta_sum_square"
  );

  require_finite_matrix(
    sigma2_draws,
    "sigma2_draws"
  );

  require_finite_matrix(
    tau2_draws,
    "tau2_draws"
  );

  require_finite_matrix(
    psi_tau_draws,
    "psi_tau_draws"
  );

  require_finite_matrix(
    c2_draws,
    "c2_draws"
  );

  require_finite_matrix(
    psi_c_draws,
    "psi_c_draws"
  );

  require_finite_matrix(
    q_draws,
    "q_draws"
  );

  require_finite_matrix(
    monitored_beta_draws,
    "monitored_beta_draws"
  );

  require_positive_matrix(
    monitored_lambda2_draws,
    "monitored_lambda2_draws"
  );

  require_positive_matrix(
    monitored_xi_draws,
    "monitored_xi_draws"
  );

  const double retained_double =
    static_cast<double>(
      retained
    );

  const arma::mat beta_mean =
    beta_sum /
    retained_double;

  const arma::mat beta_mean_square =
    beta_sum_square /
    retained_double;

  const arma::vec sigma2_mean =
    sigma2_sum /
    retained_double;

  const double tau2_mean =
    tau2_sum /
    retained_double;

  const double psi_tau_mean =
    psi_tau_sum /
    retained_double;

  const arma::vec c2_mean =
    c2_sum /
    retained_double;

  const arma::vec psi_c_mean =
    psi_c_sum /
    retained_double;

  const arma::vec lambda2_mean =
    lambda2_sum /
    retained_double;

  const arma::vec xi_mean =
    xi_sum /
    retained_double;

  const arma::vec q_mean =
    q_sum /
    retained_double;

  const arma::vec q_index_mean =
    q_index_sum /
    retained_double;

  Rcpp::NumericVector q_transport_acceptance_rate(
    n_q,
    NA_REAL
  );

  Rcpp::NumericVector q_transport_local_acceptance_rate(
    n_q,
    NA_REAL
  );

  Rcpp::NumericVector q_transport_global_acceptance_rate(
    n_q,
    NA_REAL
  );

  Rcpp::NumericVector q_transport_mean_proposed_jump(
    n_q,
    NA_REAL
  );

  Rcpp::NumericVector q_transport_mean_accepted_jump(
    n_q,
    NA_REAL
  );

  for (
    int group = 0;
    group < n_q;
    ++group
  ) {

    if (
      q_transport_proposals[
        group
      ] >
        0
    ) {

      q_transport_acceptance_rate[
        group
      ] =
        static_cast<double>(
          q_transport_acceptances[
            group
          ]
        ) /
        static_cast<double>(
          q_transport_proposals[
            group
          ]
        );

      q_transport_mean_proposed_jump[
        group
      ] =
        q_transport_proposed_jump_sum[
          group
        ] /
        static_cast<double>(
          q_transport_proposals[
            group
          ]
        );
    }

    if (
      q_transport_acceptances[
        group
      ] >
        0
    ) {

      q_transport_mean_accepted_jump[
        group
      ] =
        q_transport_accepted_jump_sum[
          group
        ] /
        static_cast<double>(
          q_transport_acceptances[
            group
          ]
        );
    }

    if (
      q_transport_local_proposals[
        group
      ] >
        0
    ) {

      q_transport_local_acceptance_rate[
        group
      ] =
        static_cast<double>(
          q_transport_local_acceptances[
            group
          ]
        ) /
        static_cast<double>(
          q_transport_local_proposals[
            group
          ]
        );
    }

    if (
      q_transport_global_proposals[
        group
      ] >
        0
    ) {

      q_transport_global_acceptance_rate[
        group
      ] =
        static_cast<double>(
          q_transport_global_acceptances[
            group
          ]
        ) /
        static_cast<double>(
          q_transport_global_proposals[
            group
          ]
        );
    }
  }

  Rcpp::NumericVector c_asis_move_rate(
    n_c,
    NA_REAL
  );

  Rcpp::NumericVector c_asis_mean_absolute_log_move(
    n_c,
    NA_REAL
  );

  Rcpp::NumericVector c_asis_mean_evaluations(
    n_c,
    NA_REAL
  );

  Rcpp::NumericVector c_asis_mean_left_expansions(
    n_c,
    NA_REAL
  );

  Rcpp::NumericVector c_asis_mean_right_expansions(
    n_c,
    NA_REAL
  );

  Rcpp::NumericVector c_asis_mean_shrink_steps(
    n_c,
    NA_REAL
  );

  for (
    int group = 1;
    group < n_c;
    ++group
  ) {

    const int updates =
      c_asis_updates[
        group
      ];

    if (
      updates >
        0
    ) {

      const double updates_double =
        static_cast<double>(
          updates
        );

      c_asis_move_rate[
        group
      ] =
        static_cast<double>(
          c_asis_moves[
            group
          ]
        ) /
        updates_double;

      c_asis_mean_absolute_log_move[
        group
      ] =
        c_asis_absolute_log_move_sum[
          group
        ] /
        updates_double;

      c_asis_mean_evaluations[
        group
      ] =
        c_asis_evaluation_sum[
          group
        ] /
        updates_double;

      c_asis_mean_left_expansions[
        group
      ] =
        c_asis_left_expansion_sum[
          group
        ] /
        updates_double;

      c_asis_mean_right_expansions[
        group
      ] =
        c_asis_right_expansion_sum[
          group
        ] /
        updates_double;

      c_asis_mean_shrink_steps[
        group
      ] =
        c_asis_shrink_step_sum[
          group
        ] /
        updates_double;
    }
  }

  Rcpp::List final_state =
    Rcpp::clone(
      state
    );

  final_state["beta"] =
    beta;

  final_state["sigma2"] =
    sigma2;

  final_state["tau2"] =
    tau2;

  final_state["psi_tau"] =
    psi_tau;

  final_state["c2"] =
    c2;

  final_state["psi_c"] =
    psi_c;

  final_state["lambda2"] =
    lambda2;

  final_state["xi"] =
    xi;

  final_state["q_index"] =
    q_index;

  const double total_seconds =
    elapsed_seconds(
      total_start,
      Clock::now()
    );

  const double accounted_seconds =
    timing.variance_refresh_seconds +
    timing.beta_seconds +
    timing.q_transport_seconds +
    timing.sigma2_seconds +
    timing.tau_seconds +
    timing.c_seconds +
    timing.c_asis_seconds +
    timing.lambda_seconds +
    timing.q_gibbs_seconds +
    timing.storage_seconds;

  const double other_seconds =
    std::max(
      0.0,
      total_seconds -
      accounted_seconds
    );

  const Rcpp::List timing_output =
    Rcpp::List::create(

      Rcpp::Named("total_seconds") =
        total_seconds,

      Rcpp::Named("variance_refresh_seconds") =
        timing.variance_refresh_seconds,

      Rcpp::Named("beta_seconds") =
        timing.beta_seconds,

      Rcpp::Named("q_transport_seconds") =
        timing.q_transport_seconds,

      Rcpp::Named("sigma2_seconds") =
        timing.sigma2_seconds,

      Rcpp::Named("tau_seconds") =
        timing.tau_seconds,

      Rcpp::Named("c_seconds") =
        timing.c_seconds,

      Rcpp::Named("c_asis_seconds") =
        timing.c_asis_seconds,

      Rcpp::Named("lambda_seconds") =
        timing.lambda_seconds,

      Rcpp::Named("q_gibbs_seconds") =
        timing.q_gibbs_seconds,

      Rcpp::Named("storage_seconds") =
        timing.storage_seconds,

      Rcpp::Named("other_seconds") =
        other_seconds
    );

  Rcpp::List full_draws =
    Rcpp::List::create(

      Rcpp::Named("beta") =
        R_NilValue,

      Rcpp::Named("lambda2") =
        R_NilValue,

      Rcpp::Named("xi") =
        R_NilValue
    );

  if (keep_all) {

    full_draws["beta"] =
      full_beta_draws;

    full_draws["lambda2"] =
      full_lambda2_draws;

    full_draws["xi"] =
      full_xi_draws;
  }

  std::string stage;

  if (
    q_update ==
      "transport_after_beta"
  ) {

    stage =
      "beta_q_transport_sigma2_tau2_psi_tau_c2_psi_c";

  } else if (
    q_update ==
      "gibbs_after_beta"
  ) {

    stage =
      "beta_q_gibbs_sigma2_tau2_psi_tau_c2_psi_c";

  } else {

    stage =
      "beta_sigma2_tau2_psi_tau_c2_psi_c";
  }

  if (
    use_c_asis
  ) {

    stage +=
      "_c_asis";
  }

  stage +=
    q_update ==
      "gibbs_end"
      ? "_lambda2_xi_q_gibbs"
      : "_lambda2_xi";

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

        Rcpp::Named("c2_mean") =
          c2_mean,

        Rcpp::Named("psi_c_mean") =
          psi_c_mean,

        Rcpp::Named("lambda2_mean") =
          lambda2_mean,

        Rcpp::Named("xi_mean") =
          xi_mean,

        Rcpp::Named("q_mean") =
          q_mean,

        Rcpp::Named("q_index_mean") =
          q_index_mean
      ),

    Rcpp::Named("draws") =
      Rcpp::List::create(

        Rcpp::Named("sigma2") =
          sigma2_draws,

        Rcpp::Named("tau2") =
          tau2_draws,

        Rcpp::Named("psi_tau") =
          psi_tau_draws,

        Rcpp::Named("c2") =
          c2_draws,

        Rcpp::Named("psi_c") =
          psi_c_draws,

        Rcpp::Named("q") =
          q_draws,

        Rcpp::Named("q_index") =
          q_index_draws
      ),

    Rcpp::Named("monitor") =
      Rcpp::List::create(

        Rcpp::Named("beta_index") =
          monitor_beta_raw,

        Rcpp::Named("beta") =
          monitored_beta_draws,

        Rcpp::Named("lambda_index") =
          monitor_lambda_raw,

        Rcpp::Named("lambda2") =
          monitored_lambda2_draws,

        Rcpp::Named("xi") =
          monitored_xi_draws
      ),

    Rcpp::Named("full_draws") =
      full_draws,

    Rcpp::Named("final_state") =
      final_state,

    Rcpp::Named("q_transport") =
      Rcpp::List::create(

        Rcpp::Named("mode") =
          q_update,

        Rcpp::Named("global_probability") =
          q_transport_global_probability,

        Rcpp::Named("proposals") =
          q_transport_proposals,

        Rcpp::Named("acceptances") =
          q_transport_acceptances,

        Rcpp::Named("acceptance_rate") =
          q_transport_acceptance_rate,

        Rcpp::Named("local_proposals") =
          q_transport_local_proposals,

        Rcpp::Named("local_acceptances") =
          q_transport_local_acceptances,

        Rcpp::Named("local_acceptance_rate") =
          q_transport_local_acceptance_rate,

        Rcpp::Named("global_proposals") =
          q_transport_global_proposals,

        Rcpp::Named("global_acceptances") =
          q_transport_global_acceptances,

        Rcpp::Named("global_acceptance_rate") =
          q_transport_global_acceptance_rate,

        Rcpp::Named("mean_proposed_jump") =
          q_transport_mean_proposed_jump,

        Rcpp::Named("mean_accepted_jump") =
          q_transport_mean_accepted_jump,

        Rcpp::Named("last_log_likelihood_ratio") =
          q_transport_last_log_likelihood_ratio,

        Rcpp::Named("last_log_q_prior_ratio") =
          q_transport_last_log_q_prior_ratio,

        Rcpp::Named("last_log_proposal_ratio") =
          q_transport_last_log_proposal_ratio,

        Rcpp::Named("last_log_acceptance_ratio") =
          q_transport_last_log_acceptance_ratio,

        Rcpp::Named("last_acceptance_probability") =
          q_transport_last_acceptance_probability
      ),

    Rcpp::Named("c_asis") =
      Rcpp::List::create(

        Rcpp::Named("enabled") =
          use_c_asis,

        Rcpp::Named("every") =
          c_asis_every,

        Rcpp::Named("slice_width") =
          c_asis_slice_width,

        Rcpp::Named("maximum_step_out") =
          c_asis_maximum_step_out,

        Rcpp::Named("updates") =
          c_asis_updates,

        Rcpp::Named("moves") =
          c_asis_moves,

        Rcpp::Named("move_rate") =
          c_asis_move_rate,

        Rcpp::Named("mean_absolute_log_move") =
          c_asis_mean_absolute_log_move,

        Rcpp::Named("mean_evaluations") =
          c_asis_mean_evaluations,

        Rcpp::Named("mean_left_expansions") =
          c_asis_mean_left_expansions,

        Rcpp::Named("mean_right_expansions") =
          c_asis_mean_right_expansions,

        Rcpp::Named("mean_shrink_steps") =
          c_asis_mean_shrink_steps,

        Rcpp::Named("last_moved") =
          c_asis_last_moved,

        Rcpp::Named("last_current_scale2") =
          c_asis_last_current_scale2,

        Rcpp::Named("last_new_scale2") =
          c_asis_last_new_scale2,

        Rcpp::Named("last_transport_scale") =
          c_asis_last_transport_scale,

        Rcpp::Named("last_log_density") =
          c_asis_last_log_density,

        Rcpp::Named("last_log_slice_height") =
          c_asis_last_log_slice_height
      ),

    Rcpp::Named("timing") =
      timing_output,

    Rcpp::Named("info") =
      Rcpp::List::create(

        Rcpp::Named("engine") =
          "fast",

        Rcpp::Named("stage") =
          stage,

        Rcpp::Named("beta_algorithm_requested") =
          requested_algorithm_name,

        Rcpp::Named("beta_algorithm_resolved") =
          resolved_algorithm_name,

        Rcpp::Named("q_update") =
          q_update,

        Rcpp::Named("use_c_asis") =
          use_c_asis,

        Rcpp::Named("c_asis_every") =
          c_asis_every,

        Rcpp::Named("c_asis_slice_width") =
          c_asis_slice_width,

        Rcpp::Named("c_asis_maximum_step_out") =
          c_asis_maximum_step_out,

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

        Rcpp::Named("n_monitor_lambda") =
          static_cast<int>(
            monitor_lambda.n_elem
          ),

        Rcpp::Named("n_c") =
          n_c,

        Rcpp::Named("n_q") =
          n_q,

        Rcpp::Named("n_blocks") =
          n_blocks
      )
  );
}