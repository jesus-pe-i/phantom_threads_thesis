// Fast half-t BVAR one-chain engine.
//
// Coefficient-wise half-t shrinkage with Minnesota lag weights,
// optional self-diagonal global grouping, and non-self tau ASIS.

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include "half_t_bvar_updates.h"


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


  inline bool list_has(
      const Rcpp::List& object,
      const char* name) {

    return object.containsElementNamed(
      name
    );
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
      draws_u *
      thin_u;
  }


  struct PairMonitorIndex {

    arma::uvec row;
    arma::uvec equation;
  };


  inline PairMonitorIndex validate_pair_monitor(
      const arma::imat& raw,
      const int k,
      const int N,
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

      PairMonitorIndex empty;

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
        name + " must have columns row and equation"
      );
    }

    PairMonitorIndex result;

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
          name + " contains an index outside the coefficient matrix"
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


  inline void validate_chain_input(
      const arma::mat& Y,
      const arma::mat& X,
      const arma::mat& XtX,
      const arma::mat& XtY,
      const arma::vec& YtY,
      const arma::vec& phi2,
      const arma::ivec& predictor_series,
      const int grouping_code,
      const int n_tau,
      const arma::uvec& tau_group_counts,
      const arma::mat& beta,
      const arma::vec& sigma2,
      const arma::mat& lambda2,
      const arma::mat& nu,
      const arma::vec& tau2,
      const arma::vec& xi,
      const int T_p,
      const int N,
      const int k) {

    if (
      T_p < 1 ||
      N < 1 ||
      k < 1
    ) {
      Rcpp::stop(
        "T_p, N and k must be positive"
      );
    }

    const arma::uword T_p_u =
      static_cast<arma::uword>(
        T_p
      );

    const arma::uword N_u =
      static_cast<arma::uword>(
        N
      );

    const arma::uword k_u =
      static_cast<arma::uword>(
        k
      );

    if (
      Y.n_rows != T_p_u ||
      Y.n_cols != N_u ||
      X.n_rows != T_p_u ||
      X.n_cols != k_u ||
      XtX.n_rows != k_u ||
      XtX.n_cols != k_u ||
      XtY.n_rows != k_u ||
      XtY.n_cols != N_u ||
      YtY.n_elem != N_u
    ) {
      Rcpp::stop(
        "Data dimensions are internally inconsistent"
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
      XtX,
      "XtX"
    );

    require_finite_matrix(
      XtY,
      "XtY"
    );

    if (
      !YtY.is_finite() ||
      arma::any(
        YtY < 0.0
      )
    ) {
      Rcpp::stop(
        "YtY must contain finite non-negative values"
      );
    }

    if (
      phi2.n_elem != k_u ||
      !phi2.is_finite() ||
      arma::any(
        phi2 <= 0.0
      )
    ) {
      Rcpp::stop(
        "phi2 must be a finite positive vector of length k"
      );
    }

    if (
      predictor_series.n_elem != k_u ||
      arma::any(
        predictor_series < 1
      ) ||
      arma::any(
        predictor_series > N
      )
    ) {
      Rcpp::stop(
        "predictor_series must contain 1-based series indices"
      );
    }

    if (
      grouping_code !=
        half_t_bvar::grouping_all &&
      grouping_code !=
        half_t_bvar::grouping_self_diagonal
    ) {
      Rcpp::stop(
        "Unknown global grouping code"
      );
    }

    const int expected_n_tau =
      grouping_code ==
        half_t_bvar::grouping_all ?
        1 :
        2;

    if (
      n_tau != expected_n_tau ||
      tau_group_counts.n_elem !=
        static_cast<arma::uword>(
          n_tau
        )
    ) {
      Rcpp::stop(
        "Global-scale dimensions do not match grouping_code"
      );
    }

    std::uint64_t counted_coefficients =
      0u;

    for (
      arma::uword group = 0u;
      group < tau_group_counts.n_elem;
      ++group
    ) {

      if (
        tau_group_counts[group] <
        1u
      ) {
        Rcpp::stop(
          "Every global shrinkage group must contain coefficients"
        );
      }

      counted_coefficients +=
        static_cast<std::uint64_t>(
          tau_group_counts[group]
        );
    }

    const std::uint64_t expected_coefficients =
      static_cast<std::uint64_t>(
        k_u
      ) *
      static_cast<std::uint64_t>(
        N_u
      );

    if (
      counted_coefficients !=
      expected_coefficients
    ) {
      Rcpp::stop(
        "tau_group_counts do not cover all coefficients"
      );
    }

    if (
      beta.n_rows != k_u ||
      beta.n_cols != N_u ||
      lambda2.n_rows != k_u ||
      lambda2.n_cols != N_u ||
      nu.n_rows != k_u ||
      nu.n_cols != N_u ||
      sigma2.n_elem != N_u ||
      tau2.n_elem !=
        static_cast<arma::uword>(
          n_tau
        ) ||
      xi.n_elem !=
        static_cast<arma::uword>(
          n_tau
        )
    ) {
      Rcpp::stop(
        "Initial-state dimensions do not match the data"
      );
    }

    require_finite_matrix(
      beta,
      "beta"
    );

    require_positive_vector(
      sigma2,
      "sigma2"
    );

    require_positive_matrix(
      lambda2,
      "lambda2"
    );

    require_positive_matrix(
      nu,
      "nu"
    );

    require_positive_vector(
      tau2,
      "tau2"
    );

    require_positive_vector(
      xi,
      "xi"
    );

    const double symmetry_error =
      arma::abs(
        XtX -
        XtX.t()
      ).max();

    const double XtX_scale =
      std::max(
        1.0,
        arma::abs(
          XtX
        ).max()
      );

    if (
      !std::isfinite(
        symmetry_error
      ) ||
      symmetry_error >
        1e-10 *
        XtX_scale
    ) {
      Rcpp::stop(
        "XtX must be symmetric"
      );
    }
  }


  struct ChainTiming {

    double beta_seconds =
      0.0;

    double sigma2_seconds =
      0.0;

    double local_seconds =
      0.0;

    double global_seconds =
      0.0;

    double asis_seconds =
      0.0;

    double storage_seconds =
      0.0;
  };
}


// Conditional test wrappers -----

// [[Rcpp::export]]
arma::vec half_t_bvar_test_rinvgamma(
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

  arma::vec result(
    static_cast<arma::uword>(
      n
    )
  );

  for (
    int draw = 0;
    draw < n;
    ++draw
  ) {
    result[
      static_cast<arma::uword>(
        draw
      )
    ] =
      half_t_bvar::rinvgamma(
        shape,
        scale
      );
  }

  return result;
}


// [[Rcpp::export]]
Rcpp::List half_t_bvar_draw_beta_cpp(
    const arma::mat& X,
    const arma::mat& Y,
    const int equation,
    const double sigma2,
    const arma::vec& omega,
    const std::string beta_algorithm = "auto") {

  Rcpp::RNGScope rng_scope;

  if (
    X.n_rows != Y.n_rows ||
    equation < 1 ||
    equation >
      static_cast<int>(
        Y.n_cols
      ) ||
    omega.n_elem != X.n_cols
  ) {
    Rcpp::stop(
      "Invalid dimensions in beta smoke wrapper"
    );
  }

  require_finite_matrix(
    X,
    "X"
  );

  require_finite_matrix(
    Y,
    "Y"
  );

  require_positive_vector(
    omega,
    "omega"
  );

  require_positive_finite(
    sigma2,
    "sigma2"
  );

  const int requested_algorithm =
    half_t_bvar::parse_beta_algorithm(
      beta_algorithm
    );

  const int resolved_algorithm =
    half_t_bvar::resolve_beta_algorithm(
      requested_algorithm,
      X.n_cols,
      X.n_rows
    );

  const arma::mat XtX =
    X.t() *
    X;

  const arma::mat XtY =
    X.t() *
    Y;

  half_t_bvar::CholWorkspace chol_workspace(
    X.n_cols
  );

  half_t_bvar::BhattacharyaWorkspace
    bhattacharya_workspace(
      X.n_rows,
      X.n_cols
    );

  arma::vec beta(
    X.n_cols
  );

  half_t_bvar::draw_beta(
    beta,
    X,
    Y,
    XtX,
    XtY,
    static_cast<arma::uword>(
      equation - 1
    ),
    omega,
    sigma2,
    resolved_algorithm,
    chol_workspace,
    bhattacharya_workspace
  );

  return Rcpp::List::create(

    Rcpp::Named("beta") =
      beta,

    Rcpp::Named("algorithm") =
      half_t_bvar::beta_algorithm_name(
        resolved_algorithm
      )
  );
}


// [[Rcpp::export]]
double half_t_bvar_draw_sigma2_cpp(
    const arma::mat& X,
    const arma::vec& y,
    const arma::vec& beta,
    const arma::vec& omega,
    const double sigma_a = 3.0,
    const double sigma_b = 2.0) {

  Rcpp::RNGScope rng_scope;

  return half_t_bvar::draw_sigma2(
    y,
    X,
    beta,
    omega,
    sigma_a,
    sigma_b
  );
}


// [[Rcpp::export]]
Rcpp::List half_t_bvar_update_hierarchy_cpp(
    arma::mat lambda2,
    arma::mat nu,
    arma::vec tau2,
    arma::vec xi,
    const arma::mat& beta,
    const arma::vec& sigma2,
    const arma::vec& phi2,
    const arma::ivec& predictor_series,
    const int grouping_code,
    const arma::uvec& tau_group_counts,
    const double lambda_df = 3.0,
    const double lambda_scale = 1.0,
    const double tau_df = 10.0,
    const double tau_scale = 1.0) {

  Rcpp::RNGScope rng_scope;

  const int expected_n_tau =
    grouping_code ==
      half_t_bvar::grouping_all ?
      1 :
      grouping_code ==
        half_t_bvar::grouping_self_diagonal ?
        2 :
        0;

  if (
    expected_n_tau == 0 ||
    lambda2.n_rows != beta.n_rows ||
    lambda2.n_cols != beta.n_cols ||
    nu.n_rows != beta.n_rows ||
    nu.n_cols != beta.n_cols ||
    sigma2.n_elem != beta.n_cols ||
    phi2.n_elem != beta.n_rows ||
    predictor_series.n_elem != beta.n_rows ||
    tau2.n_elem !=
      static_cast<arma::uword>(
        expected_n_tau
      ) ||
    xi.n_elem !=
      static_cast<arma::uword>(
        expected_n_tau
      ) ||
    tau_group_counts.n_elem !=
      static_cast<arma::uword>(
        expected_n_tau
      )
  ) {
    Rcpp::stop(
      "Invalid dimensions or grouping in hierarchy smoke wrapper"
    );
  }

  if (
    arma::any(
      predictor_series < 1
    ) ||
    arma::any(
      predictor_series >
        static_cast<int>(
          beta.n_cols
        )
    )
  ) {
    Rcpp::stop(
      "predictor_series must be 1-based"
    );
  }

  const arma::uvec predictor_series_zero =
    arma::conv_to<arma::uvec>::from(
      predictor_series -
      1
    );

  arma::vec tau_quadratic(
    tau2.n_elem,
    arma::fill::zeros
  );

  half_t_bvar::update_local_hierarchy(
    lambda2,
    nu,
    tau_quadratic,
    beta,
    sigma2,
    phi2,
    predictor_series_zero,
    tau2,
    grouping_code,
    lambda_df,
    lambda_scale
  );

  half_t_bvar::update_global_hierarchy(
    tau2,
    xi,
    tau_quadratic,
    tau_group_counts,
    tau_df,
    tau_scale
  );

  return Rcpp::List::create(

    Rcpp::Named("lambda2") =
      lambda2,

    Rcpp::Named("nu") =
      nu,

    Rcpp::Named("tau2") =
      tau2,

    Rcpp::Named("xi") =
      xi,

    Rcpp::Named("tau_quadratic") =
      tau_quadratic
  );
}


// One-chain kernel -----

// [[Rcpp::export]]
Rcpp::List half_t_bvar_chain_cpp(
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


  // Inputs -----

  const arma::mat Y =
    Rcpp::as<arma::mat>(
      data["Y"]
    );

  const arma::mat X =
    Rcpp::as<arma::mat>(
      data["X"]
    );

  const arma::mat XtX =
    Rcpp::as<arma::mat>(
      data["XtX"]
    );

  const arma::mat XtY =
    Rcpp::as<arma::mat>(
      data["XtY"]
    );

  const arma::vec YtY =
    Rcpp::as<arma::vec>(
      data["YtY"]
    );

  const int T_p =
    Rcpp::as<int>(
      data["T_p"]
    );

  const int N =
    Rcpp::as<int>(
      data["N"]
    );

  const int k =
    Rcpp::as<int>(
      data["k"]
    );

  const arma::vec phi2 =
    Rcpp::as<arma::vec>(
      maps["phi2"]
    );

  const arma::ivec predictor_series_raw =
    Rcpp::as<arma::ivec>(
      maps["predictor_series"]
    );

  const int grouping_code =
    Rcpp::as<int>(
      maps["grouping_code"]
    );

  const int n_tau =
    Rcpp::as<int>(
      maps["n_tau"]
    );

  const arma::uvec tau_group_counts =
    Rcpp::as<arma::uvec>(
      maps["tau_group_counts"]
    );

  const double lambda_df =
    Rcpp::as<double>(
      prior["lambda_df"]
    );

  const double lambda_scale =
    Rcpp::as<double>(
      prior["lambda_scale"]
    );

  const double tau_df =
    Rcpp::as<double>(
      prior["tau_df"]
    );

  const double tau_scale =
    Rcpp::as<double>(
      prior["tau_scale"]
    );

  const double sigma_prior_shape =
    Rcpp::as<double>(
      prior["sigma_a"]
    );

  const double sigma_prior_scale =
    Rcpp::as<double>(
      prior["sigma_b"]
    );

  arma::mat beta =
    Rcpp::as<arma::mat>(
      state["beta"]
    );

  arma::vec sigma2 =
    Rcpp::as<arma::vec>(
      state["sigma2"]
    );

  arma::mat lambda2 =
    Rcpp::as<arma::mat>(
      state["lambda2"]
    );

  arma::mat nu =
    Rcpp::as<arma::mat>(
      state["nu"]
    );

  arma::vec tau2 =
    Rcpp::as<arma::vec>(
      state["tau2"]
    );

  arma::vec xi =
    Rcpp::as<arma::vec>(
      state["xi"]
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
      : false;

  const int asis_every =
    list_has(
      control,
      "asis_every"
    )
      ? Rcpp::as<int>(
          control["asis_every"]
        )
      : 1;

  const double asis_slice_width =
    list_has(
      control,
      "asis_slice_width"
    )
      ? Rcpp::as<double>(
          control["asis_slice_width"]
        )
      : 1.0;

  const int asis_maximum_step_out =
    list_has(
      control,
      "asis_maximum_step_out"
    )
      ? Rcpp::as<int>(
          control["asis_maximum_step_out"]
        )
      : 100;

  const int asis_maximum_shrinks =
    1000;

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
      "monitor_lambda"
    )
  ) {
    monitor_lambda_raw =
      Rcpp::as<arma::imat>(
        control["monitor_lambda"]
      );
  }


  // Validation and fixed choices -----

  validate_chain_input(
    Y,
    X,
    XtX,
    XtY,
    YtY,
    phi2,
    predictor_series_raw,
    grouping_code,
    n_tau,
    tau_group_counts,
    beta,
    sigma2,
    lambda2,
    nu,
    tau2,
    xi,
    T_p,
    N,
    k
  );

  require_positive_finite(
    lambda_df,
    "lambda_df"
  );

  require_positive_finite(
    lambda_scale,
    "lambda_scale"
  );

  require_positive_finite(
    tau_df,
    "tau_df"
  );

  require_positive_finite(
    tau_scale,
    "tau_scale"
  );

  require_positive_finite(
    sigma_prior_shape,
    "sigma_a"
  );

  require_positive_finite(
    sigma_prior_scale,
    "sigma_b"
  );

  if (asis_every < 1) {
    Rcpp::stop(
      "asis_every must be positive"
    );
  }

  require_positive_finite(
    asis_slice_width,
    "asis_slice_width"
  );

  if (
    asis_maximum_step_out <
    1
  ) {
    Rcpp::stop(
      "asis_maximum_step_out must be positive"
    );
  }

  if (
    use_asis &&
    (
      grouping_code !=
        half_t_bvar::grouping_self_diagonal ||
      n_tau != 2
    )
  ) {
    Rcpp::stop(
      "ASIS requires self-diagonal grouping with two global scales"
    );
  }

  const PairMonitorIndex monitor_beta =
    validate_pair_monitor(
      monitor_beta_raw,
      k,
      N,
      "monitor_beta"
    );

  const PairMonitorIndex monitor_lambda =
    validate_pair_monitor(
      monitor_lambda_raw,
      k,
      N,
      "monitor_lambda"
    );

  const int requested_algorithm =
    half_t_bvar::parse_beta_algorithm(
      requested_algorithm_name
    );

  const int resolved_algorithm =
    half_t_bvar::resolve_beta_algorithm(
      requested_algorithm,
      static_cast<arma::uword>(
        k
      ),
      static_cast<arma::uword>(
        T_p
      )
    );

  const std::string resolved_algorithm_name =
    half_t_bvar::beta_algorithm_name(
      resolved_algorithm
    );

  const std::uint64_t total_iterations =
    checked_total_iterations(
      burnin,
      n_draws,
      thin
    );

  const std::uint64_t burnin_u =
    static_cast<std::uint64_t>(
      burnin
    );

  const std::uint64_t thin_u =
    static_cast<std::uint64_t>(
      thin
    );

  const std::uint64_t draws_count_u =
    static_cast<std::uint64_t>(
      n_draws
    );

  const arma::uword T_p_u =
    static_cast<arma::uword>(
      T_p
    );

  const arma::uword N_u =
    static_cast<arma::uword>(
      N
    );

  const arma::uword k_u =
    static_cast<arma::uword>(
      k
    );

  const arma::uword n_tau_u =
    static_cast<arma::uword>(
      n_tau
    );

  const arma::uvec predictor_series =
    arma::conv_to<arma::uvec>::from(
      predictor_series_raw -
      1
    );


  // Workspaces -----

  half_t_bvar::CholWorkspace chol_workspace(
    k_u
  );

  half_t_bvar::BhattacharyaWorkspace
    bhattacharya_workspace(
      T_p_u,
      k_u
    );

  half_t_bvar::AsisWorkspace asis_workspace(
    T_p_u
  );

  arma::vec omega(
    k_u,
    arma::fill::zeros
  );

  arma::vec beta_equation(
    k_u,
    arma::fill::zeros
  );

  arma::vec sigma2_workspace(
    k_u,
    arma::fill::zeros
  );

  arma::vec tau_quadratic(
    n_tau_u,
    arma::fill::zeros
  );


  // Streaming summaries -----

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

  arma::mat lambda2_sum(
    k_u,
    N_u,
    arma::fill::zeros
  );

  arma::mat nu_sum(
    k_u,
    N_u,
    arma::fill::zeros
  );

  arma::vec tau2_sum(
    n_tau_u,
    arma::fill::zeros
  );

  arma::vec xi_sum(
    n_tau_u,
    arma::fill::zeros
  );


  // Retained draws -----

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
    n_tau_u,
    arma::fill::zeros
  );

  arma::mat xi_draws(
    draws_u,
    n_tau_u,
    arma::fill::zeros
  );

  arma::mat monitored_beta_draws(
    draws_u,
    monitor_beta.row.n_elem,
    arma::fill::zeros
  );

  arma::mat monitored_lambda2_draws(
    draws_u,
    monitor_lambda.row.n_elem,
    arma::fill::zeros
  );

  arma::mat monitored_nu_draws(
    draws_u,
    monitor_lambda.row.n_elem,
    arma::fill::zeros
  );

  arma::cube full_beta_draws;

  arma::cube full_lambda2_draws;

  arma::cube full_nu_draws;

  if (keep_all) {

    full_beta_draws.set_size(
      k_u,
      N_u,
      draws_u
    );

    full_lambda2_draws.set_size(
      k_u,
      N_u,
      draws_u
    );

    full_nu_draws.set_size(
      k_u,
      N_u,
      draws_u
    );
  }


  // ASIS diagnostics -----

  std::uint64_t asis_updates =
    0u;

  std::uint64_t asis_moves =
    0u;

  std::uint64_t asis_total_evaluations =
    0u;

  std::uint64_t asis_total_left_expansions =
    0u;

  std::uint64_t asis_total_right_expansions =
    0u;

  std::uint64_t asis_total_shrink_steps =
    0u;

  double asis_absolute_log_move_sum =
    0.0;

  double asis_maximum_absolute_log_move =
    0.0;

  double asis_maximum_invariant_error =
    0.0;

  double asis_last_tau2_before =
    NA_REAL;

  double asis_last_tau2_after =
    NA_REAL;

  double asis_last_xi_before =
    NA_REAL;

  double asis_last_xi_after =
    NA_REAL;

  double asis_last_likelihood_quadratic =
    NA_REAL;

  double asis_last_likelihood_linear =
    NA_REAL;

  double asis_last_log_tau_move =
    NA_REAL;

  int asis_last_evaluations =
    0;

  int asis_last_left_expansions =
    0;

  int asis_last_right_expansions =
    0;

  int asis_last_shrink_steps =
    0;


  // Gibbs loop -----

  std::uint64_t retained =
    0u;

  ChainTiming timing;

  for (
    std::uint64_t iteration = 1u;
    iteration <= total_iterations;
    ++iteration
  ) {


    // Beta and sigma2 -----

    for (
      arma::uword equation = 0u;
      equation < N_u;
      ++equation
    ) {

      half_t_bvar::fill_omega(
        omega,
        phi2,
        lambda2,
        tau2,
        predictor_series,
        equation,
        grouping_code
      );

      Clock::time_point step_start =
        Clock::now();

      half_t_bvar::draw_beta(
        beta_equation,
        X,
        Y,
        XtX,
        XtY,
        equation,
        omega,
        sigma2[equation],
        resolved_algorithm,
        chol_workspace,
        bhattacharya_workspace
      );

      beta.col(
        equation
      ) =
        beta_equation;

      timing.beta_seconds +=
        elapsed_seconds(
          step_start,
          Clock::now()
        );

      step_start =
        Clock::now();

      sigma2[equation] =
        half_t_bvar::draw_sigma2_from_sufficient_statistics(
          YtY[equation],
          XtY.col(
            equation
          ),
          XtX,
          beta_equation,
          omega,
          sigma_prior_shape,
          sigma_prior_scale,
          T_p_u,
          sigma2_workspace
        );

      timing.sigma2_seconds +=
        elapsed_seconds(
          step_start,
          Clock::now()
        );
    }


    // Local hierarchy -----

    Clock::time_point step_start =
      Clock::now();

    half_t_bvar::update_local_hierarchy(
      lambda2,
      nu,
      tau_quadratic,
      beta,
      sigma2,
      phi2,
      predictor_series,
      tau2,
      grouping_code,
      lambda_df,
      lambda_scale
    );

    timing.local_seconds +=
      elapsed_seconds(
        step_start,
        Clock::now()
      );


    // Global hierarchy -----

    step_start =
      Clock::now();

    half_t_bvar::update_global_hierarchy(
      tau2,
      xi,
      tau_quadratic,
      tau_group_counts,
      tau_df,
      tau_scale
    );

    timing.global_seconds +=
      elapsed_seconds(
        step_start,
        Clock::now()
      );


    // ASIS -----

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

      const half_t_bvar::AsisResult result =
        half_t_bvar::update_nonself_asis(
          beta,
          tau2,
          xi,
          Y,
          X,
          sigma2,
          predictor_series,
          grouping_code,
          tau_df,
          tau_scale,
          asis_slice_width,
          asis_maximum_step_out,
          asis_maximum_shrinks,
          asis_workspace
        );

      ++asis_updates;

      const double absolute_log_move =
        std::abs(
          result.log_tau_move
        );

      if (
        absolute_log_move >
        1e-12
      ) {
        ++asis_moves;
      }

      asis_absolute_log_move_sum +=
        absolute_log_move;

      asis_maximum_absolute_log_move =
        std::max(
          asis_maximum_absolute_log_move,
          absolute_log_move
        );

      asis_maximum_invariant_error =
        std::max(
          asis_maximum_invariant_error,
          result.maximum_relative_invariant_error
        );

      asis_total_evaluations +=
        static_cast<std::uint64_t>(
          result.evaluations
        );

      asis_total_left_expansions +=
        static_cast<std::uint64_t>(
          result.left_expansions
        );

      asis_total_right_expansions +=
        static_cast<std::uint64_t>(
          result.right_expansions
        );

      asis_total_shrink_steps +=
        static_cast<std::uint64_t>(
          result.shrink_steps
        );

      asis_last_tau2_before =
        result.tau2_before;

      asis_last_tau2_after =
        result.tau2_after;

      asis_last_xi_before =
        result.xi_before;

      asis_last_xi_after =
        result.xi_after;

      asis_last_likelihood_quadratic =
        result.likelihood_quadratic;

      asis_last_likelihood_linear =
        result.likelihood_linear;

      asis_last_log_tau_move =
        result.log_tau_move;

      asis_last_evaluations =
        result.evaluations;

      asis_last_left_expansions =
        result.left_expansions;

      asis_last_right_expansions =
        result.right_expansions;

      asis_last_shrink_steps =
        result.shrink_steps;

      timing.asis_seconds +=
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

      lambda2_sum +=
        lambda2;

      nu_sum +=
        nu;

      tau2_sum +=
        tau2;

      xi_sum +=
        xi;

      sigma2_draws.row(
        retained_index
      ) =
        sigma2.t();

      tau2_draws.row(
        retained_index
      ) =
        tau2.t();

      xi_draws.row(
        retained_index
      ) =
        xi.t();

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
          monitor_lambda.row.n_elem;
        ++monitor
      ) {

        const arma::uword row =
          monitor_lambda.row[
            monitor
          ];

        const arma::uword equation =
          monitor_lambda.equation[
            monitor
          ];

        monitored_lambda2_draws(
          retained_index,
          monitor
        ) =
          lambda2(
            row,
            equation
          );

        monitored_nu_draws(
          retained_index,
          monitor
        ) =
          nu(
            row,
            equation
          );
      }

      if (keep_all) {

        full_beta_draws.slice(
          retained_index
        ) =
          beta;

        full_lambda2_draws.slice(
          retained_index
        ) =
          lambda2;

        full_nu_draws.slice(
          retained_index
        ) =
          nu;
      }

      ++retained;

      timing.storage_seconds +=
        elapsed_seconds(
          step_start,
          Clock::now()
        );
    }


    // Periodic state validation -----

    if (
      iteration %
        1000u ==
        0u
    ) {

      Rcpp::checkUserInterrupt();

      if (
        !beta.is_finite() ||
        !sigma2.is_finite() ||
        !lambda2.is_finite() ||
        !nu.is_finite() ||
        !tau2.is_finite() ||
        !xi.is_finite() ||
        arma::any(
          sigma2 <= 0.0
        ) ||
        arma::any(
          arma::vectorise(
            lambda2
          ) <= 0.0
        ) ||
        arma::any(
          arma::vectorise(
            nu
          ) <= 0.0
        ) ||
        arma::any(
          tau2 <= 0.0
        ) ||
        arma::any(
          xi <= 0.0
        )
      ) {
        Rcpp::stop(
          "Non-finite or non-positive state detected during MCMC"
        );
      }
    }
  }


  // Completion validation -----

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
    xi_draws,
    "xi_draws"
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
    monitored_nu_draws,
    "monitored_nu_draws"
  );

  if (
    arma::any(
      arma::vectorise(
        sigma2_draws
      ) <= 0.0
    ) ||
    arma::any(
      arma::vectorise(
        tau2_draws
      ) <= 0.0
    ) ||
    arma::any(
      arma::vectorise(
        xi_draws
      ) <= 0.0
    )
  ) {
    Rcpp::stop(
      "Invalid retained scale draws"
    );
  }


  // Posterior summaries -----

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

  const arma::mat lambda2_mean =
    lambda2_sum /
    retained_double;

  const arma::mat nu_mean =
    nu_sum /
    retained_double;

  const arma::vec tau2_mean =
    tau2_sum /
    retained_double;

  const arma::vec xi_mean =
    xi_sum /
    retained_double;


  // ASIS summaries -----

  const double asis_updates_double =
    static_cast<double>(
      asis_updates
    );

  const double asis_move_rate =
    asis_updates > 0u
      ? static_cast<double>(
          asis_moves
        ) /
        asis_updates_double
      : NA_REAL;

  const double asis_mean_absolute_log_move =
    asis_updates > 0u
      ? asis_absolute_log_move_sum /
        asis_updates_double
      : NA_REAL;

  const double asis_mean_evaluations =
    asis_updates > 0u
      ? static_cast<double>(
          asis_total_evaluations
        ) /
        asis_updates_double
      : NA_REAL;

  const double asis_mean_left_expansions =
    asis_updates > 0u
      ? static_cast<double>(
          asis_total_left_expansions
        ) /
        asis_updates_double
      : NA_REAL;

  const double asis_mean_right_expansions =
    asis_updates > 0u
      ? static_cast<double>(
          asis_total_right_expansions
        ) /
        asis_updates_double
      : NA_REAL;

  const double asis_mean_shrink_steps =
    asis_updates > 0u
      ? static_cast<double>(
          asis_total_shrink_steps
        ) /
        asis_updates_double
      : NA_REAL;


  // Final state and timing -----

  Rcpp::List final_state =
    Rcpp::clone(
      state
    );

  final_state["beta"] =
    beta;

  final_state["sigma2"] =
    sigma2;

  final_state["lambda2"] =
    lambda2;

  final_state["nu"] =
    nu;

  final_state["tau2"] =
    tau2;

  final_state["xi"] =
    xi;

  const double total_seconds =
    elapsed_seconds(
      total_start,
      Clock::now()
    );

  const double accounted_seconds =
    timing.beta_seconds +
    timing.sigma2_seconds +
    timing.local_seconds +
    timing.global_seconds +
    timing.asis_seconds +
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

      Rcpp::Named("beta_seconds") =
        timing.beta_seconds,

      Rcpp::Named("sigma2_seconds") =
        timing.sigma2_seconds,

      Rcpp::Named("local_seconds") =
        timing.local_seconds,

      Rcpp::Named("global_seconds") =
        timing.global_seconds,

      Rcpp::Named("asis_seconds") =
        timing.asis_seconds,

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

      Rcpp::Named("nu") =
        R_NilValue
    );

  if (keep_all) {

    full_draws["beta"] =
      full_beta_draws;

    full_draws["lambda2"] =
      full_lambda2_draws;

    full_draws["nu"] =
      full_nu_draws;
  }

  std::string stage =
    "beta_sigma2_lambda2_nu_tau2_xi";

  if (use_asis) {
    stage +=
      "_asis";
  }


  // Output -----

  return Rcpp::List::create(

    Rcpp::Named("posterior") =
      Rcpp::List::create(

        Rcpp::Named("beta_mean") =
          beta_mean,

        Rcpp::Named("beta_mean_square") =
          beta_mean_square,

        Rcpp::Named("sigma2_mean") =
          sigma2_mean,

        Rcpp::Named("lambda2_mean") =
          lambda2_mean,

        Rcpp::Named("nu_mean") =
          nu_mean,

        Rcpp::Named("tau2_mean") =
          tau2_mean,

        Rcpp::Named("xi_mean") =
          xi_mean
      ),

    Rcpp::Named("draws") =
      Rcpp::List::create(

        Rcpp::Named("sigma2") =
          sigma2_draws,

        Rcpp::Named("tau2") =
          tau2_draws,

        Rcpp::Named("xi") =
          xi_draws
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

        Rcpp::Named("nu") =
          monitored_nu_draws
      ),

    Rcpp::Named("full_draws") =
      full_draws,

    Rcpp::Named("final_state") =
      final_state,

    Rcpp::Named("asis") =
      Rcpp::List::create(

        Rcpp::Named("enabled") =
          use_asis,

        Rcpp::Named("every") =
          asis_every,

        Rcpp::Named("slice_width") =
          asis_slice_width,

        Rcpp::Named("maximum_step_out") =
          asis_maximum_step_out,

        Rcpp::Named("maximum_shrinks") =
          asis_maximum_shrinks,

        Rcpp::Named("updates") =
          static_cast<double>(
            asis_updates
          ),

        Rcpp::Named("moves") =
          static_cast<double>(
            asis_moves
          ),

        Rcpp::Named("move_rate") =
          asis_move_rate,

        Rcpp::Named("mean_absolute_log_move") =
          asis_mean_absolute_log_move,

        Rcpp::Named("maximum_absolute_log_move") =
          asis_updates > 0u
            ? asis_maximum_absolute_log_move
            : NA_REAL,

        Rcpp::Named("mean_evaluations") =
          asis_mean_evaluations,

        Rcpp::Named("mean_left_expansions") =
          asis_mean_left_expansions,

        Rcpp::Named("mean_right_expansions") =
          asis_mean_right_expansions,

        Rcpp::Named("mean_shrink_steps") =
          asis_mean_shrink_steps,

        Rcpp::Named("maximum_relative_invariant_error") =
          asis_updates > 0u
            ? asis_maximum_invariant_error
            : NA_REAL,

        Rcpp::Named("last_tau2_before") =
          asis_last_tau2_before,

        Rcpp::Named("last_tau2_after") =
          asis_last_tau2_after,

        Rcpp::Named("last_xi_before") =
          asis_last_xi_before,

        Rcpp::Named("last_xi_after") =
          asis_last_xi_after,

        Rcpp::Named("last_likelihood_quadratic") =
          asis_last_likelihood_quadratic,

        Rcpp::Named("last_likelihood_linear") =
          asis_last_likelihood_linear,

        Rcpp::Named("last_log_tau_move") =
          asis_last_log_tau_move,

        Rcpp::Named("last_evaluations") =
          asis_last_evaluations,

        Rcpp::Named("last_left_expansions") =
          asis_last_left_expansions,

        Rcpp::Named("last_right_expansions") =
          asis_last_right_expansions,

        Rcpp::Named("last_shrink_steps") =
          asis_last_shrink_steps
      ),

    Rcpp::Named("timing") =
      timing_output,

    Rcpp::Named("info") =
      Rcpp::List::create(

        Rcpp::Named("engine") =
          "fast",

        Rcpp::Named("stage") =
          stage,

        Rcpp::Named("T_p") =
          T_p,

        Rcpp::Named("N") =
          N,

        Rcpp::Named("k") =
          k,

        Rcpp::Named("n_tau") =
          n_tau,

        Rcpp::Named("grouping_code") =
          grouping_code,

        Rcpp::Named("beta_algorithm_requested") =
          requested_algorithm_name,

        Rcpp::Named("beta_algorithm_resolved") =
          resolved_algorithm_name,

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

        Rcpp::Named("use_asis") =
          use_asis,

        Rcpp::Named("asis_every") =
          asis_every,

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
            monitor_lambda.row.n_elem
          ),

        Rcpp::Named("draw_orientation") =
          "draw_by_parameter"
      )
  );
}
