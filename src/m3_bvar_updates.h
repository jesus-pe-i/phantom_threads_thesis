#ifndef M3_BVAR_UPDATES_H
#define M3_BVAR_UPDATES_H

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>


// M3 BVAR conditional-update layer.
// Production kernels, compatibility wrappers, c-ASIS, and q transport.

namespace m3_bvar {


  // Helpers -----

  constexpr int beta_auto = 0;
  constexpr int beta_chol = 1;
  constexpr int beta_bhattacharya = 2;


  inline double rinvgamma(
      const double shape,
      const double scale) {

    const double gamma_draw =
      R::rgamma(
        shape,
        1.0 / scale
      );

    if (
      !std::isfinite(gamma_draw) ||
      gamma_draw <= 0.0
    ) {
      Rcpp::stop(
        "Invalid inverse-gamma draw"
      );
    }

    const double result =
      1.0 / gamma_draw;

    if (
      !std::isfinite(result) ||
      result <= 0.0
    ) {
      Rcpp::stop(
        "Invalid inverse-gamma result"
      );
    }

    return result;
  }


  inline arma::uword sample_log_weights(
      const arma::vec& log_weights) {

    if (
      log_weights.n_elem < 1u ||
      log_weights.has_nan()
    ) {
      Rcpp::stop(
        "log_weights must contain at least one non-NaN value"
      );
    }

    const double maximum =
      log_weights.max();

    if (!std::isfinite(maximum)) {

      if (
        maximum ==
        -std::numeric_limits<double>::infinity()
      ) {
        Rcpp::stop(
          "All categorical log weights are minus infinity"
        );
      }

      Rcpp::stop(
        "Categorical log weights contain an invalid value"
      );
    }

    arma::vec probabilities =
      arma::exp(
        log_weights -
        maximum
      );

    const double total =
      arma::accu(
        probabilities
      );

    if (
      !std::isfinite(total) ||
      total <= 0.0
    ) {
      Rcpp::stop(
        "Categorical probabilities have invalid total mass"
      );
    }

    probabilities /= total;

    const double uniform =
      R::runif(
        0.0,
        1.0
      );

    double cumulative = 0.0;

    for (
      arma::uword j = 0u;
      j < probabilities.n_elem;
      ++j
    ) {

      cumulative += probabilities[j];

      if (
        uniform <= cumulative ||
        j + 1u ==
        probabilities.n_elem
      ) {
        return j;
      }
    }

    return probabilities.n_elem - 1u;
  }


  inline int parse_beta_algorithm(
      const std::string& algorithm) {

    if (algorithm == "auto") {
      return beta_auto;
    }

    if (algorithm == "chol") {
      return beta_chol;
    }

    if (algorithm == "bhattacharya") {
      return beta_bhattacharya;
    }

    Rcpp::stop(
      "beta_algorithm must be 'auto', 'chol' or 'bhattacharya'"
    );

    return beta_auto;
  }


  inline std::string beta_algorithm_name(
      const int algorithm) {

    if (algorithm == beta_chol) {
      return "chol";
    }

    if (algorithm == beta_bhattacharya) {
      return "bhattacharya";
    }

    if (algorithm == beta_auto) {
      return "auto";
    }

    Rcpp::stop(
      "Unknown beta algorithm code"
    );

    return "unknown";
  }


  inline int resolve_beta_algorithm(
      const int requested_algorithm,
      const arma::uword k,
      const arma::uword T_p) {

    if (
      requested_algorithm == beta_chol ||
      requested_algorithm ==
        beta_bhattacharya
    ) {
      return requested_algorithm;
    }

    if (
      requested_algorithm !=
      beta_auto
    ) {
      Rcpp::stop(
        "Unknown beta algorithm code"
      );
    }

    return (
      k <=
      2ULL *
      T_p
    ) ?
      beta_chol :
      beta_bhattacharya;
  }


  // Gaussian workspaces -----

  struct CholWorkspace {

    arma::mat precision;
    arma::mat U;

    arma::vec rhs;
    arma::vec intermediate;
    arma::vec posterior_mean;
    arma::vec z;
    arma::vec perturbation;

    explicit CholWorkspace(
        const arma::uword k) :

      precision(
        k,
        k,
        arma::fill::zeros
      ),

      U(
        k,
        k,
        arma::fill::zeros
      ),

      rhs(
        k,
        arma::fill::zeros
      ),

      intermediate(
        k,
        arma::fill::zeros
      ),

      posterior_mean(
        k,
        arma::fill::zeros
      ),

      z(
        k,
        arma::fill::zeros
      ),

      perturbation(
        k,
        arma::fill::zeros
      ) {}
  };


  struct BhattacharyaWorkspace {

    arma::vec u;
    arma::vec delta;
    arma::vec v;
    arma::vec rhs;
    arma::vec intermediate;
    arma::vec w;
    arma::vec Xt_w;

    arma::mat XD;
    arma::mat system;
    arma::mat U;

    BhattacharyaWorkspace(
        const arma::uword T_p,
        const arma::uword k) :

      u(
        k,
        arma::fill::zeros
      ),

      delta(
        T_p,
        arma::fill::zeros
      ),

      v(
        T_p,
        arma::fill::zeros
      ),

      rhs(
        T_p,
        arma::fill::zeros
      ),

      intermediate(
        T_p,
        arma::fill::zeros
      ),

      w(
        T_p,
        arma::fill::zeros
      ),

      Xt_w(
        k,
        arma::fill::zeros
      ),

      XD(
        T_p,
        k,
        arma::fill::zeros
      ),

      system(
        T_p,
        T_p,
        arma::fill::zeros
      ),

      U(
        T_p,
        T_p,
        arma::fill::zeros
      ) {}
  };


  // Beta updates -----

  inline void draw_beta_chol(
      arma::vec& beta,
      const arma::mat& XtX,
      const arma::mat& XtY,
      const arma::uword equation,
      const arma::vec& omega,
      const double sigma2,
      CholWorkspace& workspace) {

    const arma::uword k =
      XtX.n_rows;

    workspace.precision =
      XtX;

    for (
      arma::uword r = 0u;
      r < k;
      ++r
    ) {
      workspace.precision(r, r) +=
        1.0 /
        omega[r];
    }

    if (
      !arma::chol(
        workspace.U,
        workspace.precision
      )
    ) {
      Rcpp::stop(
        "Beta posterior precision is not positive definite"
      );
    }

    workspace.rhs =
      XtY.col(
        equation
      );

    workspace.intermediate =
      arma::solve(
        arma::trimatl(
          workspace.U.t()
        ),
        workspace.rhs,
        arma::solve_opts::fast
      );

    workspace.posterior_mean =
      arma::solve(
        arma::trimatu(
          workspace.U
        ),
        workspace.intermediate,
        arma::solve_opts::fast
      );

    for (
      arma::uword r = 0u;
      r < k;
      ++r
    ) {

      workspace.z[r] =
        R::rnorm(
          0.0,
          1.0
        );
    }

    workspace.perturbation =
      arma::solve(
        arma::trimatu(
          workspace.U
        ),
        workspace.z,
        arma::solve_opts::fast
      );

    beta =
      workspace.posterior_mean +
      std::sqrt(
        sigma2
      ) *
      workspace.perturbation;
  }


  inline void draw_beta_bhattacharya(
      arma::vec& beta,
      const arma::mat& X,
      const arma::mat& Y,
      const arma::uword equation,
      const arma::vec& omega,
      const double sigma2,
      BhattacharyaWorkspace& workspace) {

    const arma::uword T_p =
      X.n_rows;

    const arma::uword k =
      X.n_cols;

    const double sigma =
      std::sqrt(
        sigma2
      );

    for (
      arma::uword r = 0u;
      r < k;
      ++r
    ) {

      workspace.u[r] =
        sigma *
        std::sqrt(
          omega[r]
        ) *
        R::rnorm(
          0.0,
          1.0
        );
    }

    for (
      arma::uword t = 0u;
      t < T_p;
      ++t
    ) {

      workspace.delta[t] =
        sigma *
        R::rnorm(
          0.0,
          1.0
        );
    }

    workspace.v =
      X *
      workspace.u +
      workspace.delta;

    workspace.XD =
      X;

    workspace.XD.each_row() %=
      omega.t();

    workspace.system =
      workspace.XD *
      X.t();

    workspace.system.diag() +=
      1.0;

    if (
      !arma::chol(
        workspace.U,
        workspace.system
      )
    ) {
      Rcpp::stop(
        "Bhattacharya system is not positive definite"
      );
    }

    workspace.rhs =
      Y.col(
        equation
      ) -
      workspace.v;

    workspace.intermediate =
      arma::solve(
        arma::trimatl(
          workspace.U.t()
        ),
        workspace.rhs,
        arma::solve_opts::fast
      );

    workspace.w =
      arma::solve(
        arma::trimatu(
          workspace.U
        ),
        workspace.intermediate,
        arma::solve_opts::fast
      );

    workspace.Xt_w =
      X.t() *
      workspace.w;

    beta =
      workspace.u +
      omega %
      workspace.Xt_w;
  }


  inline void draw_beta(
      arma::vec& beta,
      const arma::mat& X,
      const arma::mat& Y,
      const arma::mat& XtX,
      const arma::mat& XtY,
      const arma::uword equation,
      const arma::vec& omega,
      const double sigma2,
      const int resolved_algorithm,
      CholWorkspace& chol_workspace,
      BhattacharyaWorkspace& bhattacharya_workspace) {

    if (
      resolved_algorithm ==
      beta_chol
    ) {

      draw_beta_chol(
        beta,
        XtX,
        XtY,
        equation,
        omega,
        sigma2,
        chol_workspace
      );

      return;
    }

    if (
      resolved_algorithm ==
      beta_bhattacharya
    ) {

      draw_beta_bhattacharya(
        beta,
        X,
        Y,
        equation,
        omega,
        sigma2,
        bhattacharya_workspace
      );

      return;
    }

    Rcpp::stop(
      "Unknown resolved beta algorithm"
    );
  }


  inline arma::vec draw_beta_chol(
      const arma::vec& y,
      const arma::mat& X,
      const arma::vec& omega,
      const double sigma2) {

    const arma::uword T_p =
      X.n_rows;

    const arma::uword k =
      X.n_cols;

    if (
      y.n_elem != T_p ||
      omega.n_elem != k
    ) {
      Rcpp::stop(
        "Invalid beta Cholesky dimensions"
      );
    }

    if (
      !y.is_finite() ||
      !X.is_finite() ||
      !omega.is_finite() ||
      arma::any(
        omega <= 0.0
      ) ||
      !std::isfinite(
        sigma2
      ) ||
      sigma2 <= 0.0
    ) {
      Rcpp::stop(
        "Invalid beta Cholesky inputs"
      );
    }

    const arma::mat Y =
      arma::reshape(
        y,
        T_p,
        1u
      );

    const arma::mat XtX =
      X.t() *
      X;

    const arma::mat XtY =
      X.t() *
      Y;

    arma::vec beta(
      k
    );

    CholWorkspace workspace(
      k
    );

    draw_beta_chol(
      beta,
      XtX,
      XtY,
      0u,
      omega,
      sigma2,
      workspace
    );

    return beta;
  }


  inline arma::vec draw_beta_bhattacharya(
      const arma::vec& y,
      const arma::mat& X,
      const arma::vec& omega,
      const double sigma2) {

    const arma::uword T_p =
      X.n_rows;

    const arma::uword k =
      X.n_cols;

    if (
      y.n_elem != T_p ||
      omega.n_elem != k
    ) {
      Rcpp::stop(
        "Invalid Bhattacharya dimensions"
      );
    }

    if (
      !y.is_finite() ||
      !X.is_finite() ||
      !omega.is_finite() ||
      arma::any(
        omega <= 0.0
      ) ||
      !std::isfinite(
        sigma2
      ) ||
      sigma2 <= 0.0
    ) {
      Rcpp::stop(
        "Invalid Bhattacharya inputs"
      );
    }

    const arma::mat Y =
      arma::reshape(
        y,
        T_p,
        1u
      );

    arma::vec beta(
      k
    );

    BhattacharyaWorkspace workspace(
      T_p,
      k
    );

    draw_beta_bhattacharya(
      beta,
      X,
      Y,
      0u,
      omega,
      sigma2,
      workspace
    );

    return beta;
  }


  inline arma::vec draw_beta(
      const arma::vec& y,
      const arma::mat& X,
      const arma::vec& omega,
      const double sigma2,
      const std::string& algorithm) {

    const int resolved_algorithm =
      resolve_beta_algorithm(
        parse_beta_algorithm(
          algorithm
        ),
        X.n_cols,
        X.n_rows
      );

    if (
      resolved_algorithm ==
      beta_chol
    ) {
      return draw_beta_chol(
        y,
        X,
        omega,
        sigma2
      );
    }

    return draw_beta_bhattacharya(
      y,
      X,
      omega,
      sigma2
    );
  }


  // Sigma2 -----

  inline double draw_sigma2_from_residual(
      const arma::vec& residual,
      const arma::vec& beta,
      const arma::vec& omega,
      const double prior_shape,
      const double prior_scale) {

    const double rss =
      arma::dot(
        residual,
        residual
      );

    const double beta_quadratic =
      arma::dot(
        beta % beta,
        1.0 /
        omega
      );

    const double posterior_shape =
      prior_shape +
      0.5 *
      static_cast<double>(
        residual.n_elem +
        beta.n_elem
      );

    const double posterior_scale =
      prior_scale +
      0.5 *
      (
        rss +
        beta_quadratic
      );

    return rinvgamma(
      posterior_shape,
      posterior_scale
    );
  }


  inline double draw_sigma2(
      const arma::vec& y,
      const arma::mat& X,
      const arma::vec& beta,
      const arma::vec& omega,
      const double prior_shape,
      const double prior_scale) {

    if (
      X.n_rows != y.n_elem ||
      X.n_cols != beta.n_elem ||
      omega.n_elem != beta.n_elem
    ) {
      Rcpp::stop(
        "Invalid sigma2 dimensions"
      );
    }

    if (
      !y.is_finite() ||
      !X.is_finite() ||
      !beta.is_finite() ||
      !omega.is_finite() ||
      arma::any(
        omega <= 0.0
      ) ||
      !std::isfinite(
        prior_shape
      ) ||
      prior_shape <= 0.0 ||
      !std::isfinite(
        prior_scale
      ) ||
      prior_scale <= 0.0
    ) {
      Rcpp::stop(
        "Invalid sigma2 inputs"
      );
    }

    const arma::vec residual =
      y -
      X *
      beta;

    return draw_sigma2_from_residual(
      residual,
      beta,
      omega,
      prior_shape,
      prior_scale
    );
  }


  // Half-t hierarchy -----

  inline double draw_tau2(
      const arma::mat& beta,
      const arma::vec& sigma2,
      const arma::mat& base_omega,
      const double tau_df,
      const double psi_tau) {

    const arma::uword k =
      beta.n_rows;

    const arma::uword N =
      beta.n_cols;

    double quadratic =
      0.0;

    for (
      arma::uword i = 0u;
      i < N;
      ++i
    ) {

      const double inverse_sigma2 =
        1.0 /
        sigma2[i];

      for (
        arma::uword r = 0u;
        r < k;
        ++r
      ) {

        const double coefficient =
          beta(r, i);

        quadratic +=
          coefficient *
          coefficient *
          inverse_sigma2 /
          base_omega(r, i);
      }
    }

    const double posterior_shape =
      0.5 *
      (
        tau_df +
        static_cast<double>(
          beta.n_elem
        )
      );

    const double posterior_scale =
      tau_df /
      psi_tau +
      0.5 *
      quadratic;

    return rinvgamma(
      posterior_shape,
      posterior_scale
    );
  }


  inline double draw_half_t_auxiliary(
      const double scale2,
      const double degrees_freedom,
      const double prior_scale) {

    const double posterior_shape =
      0.5 *
      (
        degrees_freedom +
        1.0
      );

    const double posterior_scale =
      degrees_freedom /
      scale2 +
      1.0 /
      (
        prior_scale *
        prior_scale
      );

    return rinvgamma(
      posterior_shape,
      posterior_scale
    );
  }


  // Group caches and scale updates -----

  struct GroupMembershipCache {

    arma::uword n_rows;
    arma::uword n_equations;

    std::vector<arma::uvec> linear_indices;
    arma::uvec coefficient_count;
  };


  inline arma::uvec uvec_from_std(
      const std::vector<arma::uword>& values) {

    arma::uvec output(
      values.size()
    );

    for (
      arma::uword j = 0u;
      j < output.n_elem;
      ++j
    ) {

      output[j] =
        values[
          static_cast<std::size_t>(
            j
          )
        ];
    }

    return output;
  }


  inline GroupMembershipCache build_group_membership_cache(
      const arma::imat& group_id,
      const int n_groups,
      const std::string& map_name =
        "group_id") {

    if (n_groups < 1) {
      Rcpp::stop(
        "n_groups must be positive"
      );
    }

    const arma::uword k =
      group_id.n_rows;

    const arma::uword N =
      group_id.n_cols;

    if (
      k < 1u ||
      N < 1u
    ) {
      Rcpp::stop(
        map_name +
        " must have positive dimensions"
      );
    }

    std::vector<
      std::vector<arma::uword>
    > membership(
      static_cast<std::size_t>(
        n_groups
      )
    );

    for (
      arma::uword i = 0u;
      i < N;
      ++i
    ) {

      for (
        arma::uword r = 0u;
        r < k;
        ++r
      ) {

        const int group =
          group_id(r, i);

        if (
          group < 1 ||
          group > n_groups
        ) {
          Rcpp::stop(
            map_name +
            " contains a value outside 1:n_groups"
          );
        }

        membership[
          static_cast<std::size_t>(
            group -
            1
          )
        ].push_back(
          r +
          i *
          k
        );
      }
    }

    GroupMembershipCache cache;

    cache.n_rows =
      k;

    cache.n_equations =
      N;

    cache.linear_indices.resize(
      static_cast<std::size_t>(
        n_groups
      )
    );

    cache.coefficient_count.set_size(
      static_cast<arma::uword>(
        n_groups
      )
    );

    for (
      int group = 0;
      group < n_groups;
      ++group
    ) {

      const std::vector<arma::uword>&
        group_membership =
          membership[
            static_cast<std::size_t>(
              group
            )
          ];

      if (
        group_membership.empty()
      ) {
        Rcpp::stop(
          map_name +
          " contains an empty group"
        );
      }

      cache.linear_indices[
        static_cast<std::size_t>(
          group
        )
      ] =
        uvec_from_std(
          group_membership
        );

      cache.coefficient_count[
        static_cast<arma::uword>(
          group
        )
      ] =
        static_cast<arma::uword>(
          group_membership.size()
        );
    }

    return cache;
  }


  inline double group_scale_quadratic(
      const arma::mat& beta,
      const arma::vec& sigma2,
      const arma::mat& base_omega,
      const GroupMembershipCache& cache,
      const arma::uword group) {

    const arma::uvec& indices =
      cache.linear_indices[
        static_cast<std::size_t>(
          group
        )
      ];

    const arma::uword k =
      cache.n_rows;

    double quadratic =
      0.0;

    for (
      arma::uword position = 0u;
      position < indices.n_elem;
      ++position
    ) {

      const arma::uword linear =
        indices[position];

      const arma::uword equation =
        linear /
        k;

      const double coefficient =
        beta[linear];

      quadratic +=
        coefficient *
        coefficient /
        (
          sigma2[equation] *
          base_omega[linear]
        );
    }

    return quadratic;
  }


  inline double draw_group_scale2_cached(
      const arma::mat& beta,
      const arma::vec& sigma2,
      const arma::mat& base_omega,
      const GroupMembershipCache& cache,
      const arma::uword group,
      const double degrees_freedom,
      const double auxiliary) {

    const double quadratic =
      group_scale_quadratic(
        beta,
        sigma2,
        base_omega,
        cache,
        group
      );

    const double posterior_shape =
      0.5 *
      (
        degrees_freedom +
        static_cast<double>(
          cache.coefficient_count[
            group
          ]
        )
      );

    const double posterior_scale =
      degrees_freedom /
      auxiliary +
      0.5 *
      quadratic;

    return rinvgamma(
      posterior_shape,
      posterior_scale
    );
  }


  inline void accumulate_group_quadratics(
      arma::vec& quadratics,
      const arma::mat& beta,
      const arma::vec& sigma2,
      const arma::mat& base_omega,
      const arma::imat& group_id_zero,
      const arma::uword n_groups) {

    quadratics.zeros(
      n_groups
    );

    const arma::uword k =
      beta.n_rows;

    const arma::uword N =
      beta.n_cols;

    for (
      arma::uword i = 0u;
      i < N;
      ++i
    ) {

      const double inverse_sigma2 =
        1.0 /
        sigma2[i];

      for (
        arma::uword r = 0u;
        r < k;
        ++r
      ) {

        const arma::uword group =
          static_cast<arma::uword>(
            group_id_zero(r, i)
          );

        const double coefficient =
          beta(r, i);

        quadratics[group] +=
          coefficient *
          coefficient *
          inverse_sigma2 /
          base_omega(r, i);
      }
    }
  }


  inline void draw_all_group_scales2(
      arma::vec& scale2,
      arma::vec& quadratics,
      const arma::mat& beta,
      const arma::vec& sigma2,
      const arma::mat& base_omega,
      const arma::imat& group_id_zero,
      const arma::uvec& coefficient_count,
      const double degrees_freedom,
      const arma::vec& auxiliary,
      const arma::uword first_updated_group =
        0u) {

    const arma::uword n_groups =
      scale2.n_elem;

    accumulate_group_quadratics(
      quadratics,
      beta,
      sigma2,
      base_omega,
      group_id_zero,
      n_groups
    );

    for (
      arma::uword group =
        first_updated_group;
      group < n_groups;
      ++group
    ) {

      const double posterior_shape =
        0.5 *
        (
          degrees_freedom +
          static_cast<double>(
            coefficient_count[
              group
            ]
          )
        );

      const double posterior_scale =
        degrees_freedom /
        auxiliary[group] +
        0.5 *
        quadratics[group];

      scale2[group] =
        rinvgamma(
          posterior_shape,
          posterior_scale
        );
    }
  }


  inline double draw_group_scale2(
      const arma::mat& beta,
      const arma::vec& sigma2,
      const arma::mat& base_omega,
      const arma::imat& group_id,
      const int target_group,
      const double degrees_freedom,
      const double auxiliary) {

    if (
      beta.n_rows !=
        base_omega.n_rows ||
      beta.n_cols !=
        base_omega.n_cols ||
      beta.n_rows !=
        group_id.n_rows ||
      beta.n_cols !=
        group_id.n_cols ||
      sigma2.n_elem !=
        beta.n_cols
    ) {
      Rcpp::stop(
        "Group-scale dimensions are inconsistent"
      );
    }

    if (
      target_group < 1 ||
      !std::isfinite(
        degrees_freedom
      ) ||
      degrees_freedom <= 0.0 ||
      !std::isfinite(
        auxiliary
      ) ||
      auxiliary <= 0.0
    ) {
      Rcpp::stop(
        "Invalid group-scale parameters"
      );
    }

    const GroupMembershipCache cache =
      build_group_membership_cache(
        group_id,
        group_id.max(),
        "group_id"
      );

    const arma::uword group =
      static_cast<arma::uword>(
        target_group -
        1
      );

    if (
      group >=
      cache.linear_indices.size()
    ) {
      Rcpp::stop(
        "target_group lies outside group_id"
      );
    }

    return draw_group_scale2_cached(
      beta,
      sigma2,
      base_omega,
      cache,
      group,
      degrees_freedom,
      auxiliary
    );
  }


  // Group-c ASIS -----

  struct CAsisEquationCache {

    arma::uword response_column;
    arma::uvec rows;
  };


  struct CAsisGroupCache {

    std::vector<CAsisEquationCache>
      equations;

    arma::uword coefficient_count;
  };


  inline std::vector<CAsisGroupCache>
  build_c_asis_cache(
      const arma::imat& group_id,
      const int n_groups) {

    if (n_groups < 1) {
      Rcpp::stop(
        "n_groups must be positive"
      );
    }

    const arma::uword k =
      group_id.n_rows;

    const arma::uword N =
      group_id.n_cols;

    if (
      k < 1u ||
      N < 1u
    ) {
      Rcpp::stop(
        "c ASIS group map must have positive dimensions"
      );
    }

    const std::size_t cell_count =
      static_cast<std::size_t>(
        n_groups
      ) *
      static_cast<std::size_t>(
        N
      );

    std::vector<
      std::vector<arma::uword>
    > rows_by_cell(
      cell_count
    );

    for (
      arma::uword equation = 0u;
      equation < N;
      ++equation
    ) {

      for (
        arma::uword row = 0u;
        row < k;
        ++row
      ) {

        const int group =
          group_id(
            row,
            equation
          );

        if (
          group < 1 ||
          group > n_groups
        ) {
          Rcpp::stop(
            "c ASIS group map contains a value outside 1:n_groups"
          );
        }

        const std::size_t cell =
          static_cast<std::size_t>(
            group -
            1
          ) *
          static_cast<std::size_t>(
            N
          ) +
          static_cast<std::size_t>(
            equation
          );

        rows_by_cell[
          cell
        ].push_back(
          row
        );
      }
    }

    std::vector<CAsisGroupCache> cache(
      static_cast<std::size_t>(
        n_groups
      )
    );

    for (
      int group = 0;
      group < n_groups;
      ++group
    ) {

      CAsisGroupCache& group_cache =
        cache[
          static_cast<std::size_t>(
            group
          )
        ];

      group_cache.coefficient_count =
        0u;

      for (
        arma::uword equation = 0u;
        equation < N;
        ++equation
      ) {

        const std::size_t cell =
          static_cast<std::size_t>(
            group
          ) *
          static_cast<std::size_t>(
            N
          ) +
          static_cast<std::size_t>(
            equation
          );

        if (
          rows_by_cell[
            cell
          ].empty()
        ) {
          continue;
        }

        CAsisEquationCache equation_cache;

        equation_cache.response_column =
          equation;

        equation_cache.rows =
          uvec_from_std(
            rows_by_cell[
              cell
            ]
          );

        group_cache.coefficient_count +=
          equation_cache.rows.n_elem;

        group_cache.equations.push_back(
          equation_cache
        );
      }

      if (
        group_cache.coefficient_count ==
        0u
      ) {
        Rcpp::stop(
          "Every c ASIS group must contain coefficients"
        );
      }
    }

    return cache;
  }


  struct CAsisWorkspace {

    std::vector<arma::vec>
      fitted_values;

    double likelihood_quadratic;
    double likelihood_linear;

    CAsisWorkspace() :

      fitted_values(),

      likelihood_quadratic(
        0.0
      ),

      likelihood_linear(
        0.0
      ) {}
  };


  inline void prepare_c_asis_workspace(
      const arma::mat& beta,
      const arma::mat& residual,
      const arma::mat& X,
      const arma::vec& sigma2,
      const CAsisGroupCache& cache,
      CAsisWorkspace& workspace) {

    workspace.fitted_values.resize(
      cache.equations.size()
    );

    workspace.likelihood_quadratic =
      0.0;

    workspace.likelihood_linear =
      0.0;

    for (
      std::size_t cell = 0u;
      cell < cache.equations.size();
      ++cell
    ) {

      const CAsisEquationCache&
        equation_cache =
          cache.equations[
            cell
          ];

      const arma::uword equation =
        equation_cache.response_column;

      const arma::vec beta_column =
        beta.col(
          equation
        );

      arma::vec fitted =
        X.cols(
          equation_cache.rows
        ) *
        beta_column.elem(
          equation_cache.rows
        );

      const double inverse_sigma2 =
        1.0 /
        sigma2[
          equation
        ];

      workspace.likelihood_quadratic +=
        arma::dot(
          fitted,
          fitted
        ) *
        inverse_sigma2;

      workspace.likelihood_linear +=
        arma::dot(
          residual.col(
            equation
          ),
          fitted
        ) *
        inverse_sigma2;

      workspace.fitted_values[
        cell
      ] =
        fitted;
    }

    if (
      !std::isfinite(
        workspace.likelihood_quadratic
      ) ||
      workspace.likelihood_quadratic <
        0.0 ||
      !std::isfinite(
        workspace.likelihood_linear
      )
    ) {
      Rcpp::stop(
        "Invalid c ASIS likelihood workspace"
      );
    }
  }


  inline double c_asis_log_density(
      const double candidate_log_scale2,
      const double current_log_scale2,
      const double degrees_freedom,
      const double auxiliary,
      const CAsisWorkspace& workspace) {

    if (
      !std::isfinite(
        candidate_log_scale2
      ) ||
      !std::isfinite(
        current_log_scale2
      ) ||
      !std::isfinite(
        degrees_freedom
      ) ||
      degrees_freedom <= 0.0 ||
      !std::isfinite(
        auxiliary
      ) ||
      auxiliary <= 0.0
    ) {
      return
        -std::numeric_limits<
          double
        >::infinity();
    }

    const double log_transport =
      0.5 *
      (
        candidate_log_scale2 -
        current_log_scale2
      );

    if (
      log_transport >
        700.0 ||
      log_transport <
        -745.0
    ) {
      return
        -std::numeric_limits<
          double
        >::infinity();
    }

    const double transport =
      std::exp(
        log_transport
      );

    if (
      !std::isfinite(
        transport
      ) ||
      transport <= 0.0
    ) {
      return
        -std::numeric_limits<
          double
        >::infinity();
    }

    const double transport_change =
      transport -
      1.0;

    const double log_likelihood_change =
      -0.5 *
      (
        transport_change *
        transport_change *
        workspace.likelihood_quadratic -
        2.0 *
        transport_change *
        workspace.likelihood_linear
      );

    if (
      !std::isfinite(
        log_likelihood_change
      )
    ) {
      return
        -std::numeric_limits<
          double
        >::infinity();
    }

    const double inverse_log_scale =
      -candidate_log_scale2;

    if (
      inverse_log_scale >
      700.0
    ) {
      return
        -std::numeric_limits<
          double
        >::infinity();
    }

    const double inverse_scale2 =
      std::exp(
        inverse_log_scale
      );

    if (
      !std::isfinite(
        inverse_scale2
      )
    ) {
      return
        -std::numeric_limits<
          double
        >::infinity();
    }

    const double log_prior =
      -0.5 *
      degrees_freedom *
      candidate_log_scale2 -
      (
        degrees_freedom /
        auxiliary
      ) *
      inverse_scale2;

    const double result =
      log_likelihood_change +
      log_prior;

    return std::isfinite(
      result
    )
      ? result
      : -std::numeric_limits<
          double
        >::infinity();
  }


  struct CAsisSliceResult {

    bool moved;

    double current_log_scale2;
    double new_log_scale2;

    double current_scale2;
    double new_scale2;

    double transport_scale;

    double current_log_density;
    double new_log_density;
    double log_slice_height;

    int evaluations;
    int left_expansions;
    int right_expansions;
    int shrink_steps;
  };


  inline CAsisSliceResult draw_c_asis_log_scale2(
      const double current_scale2,
      const double degrees_freedom,
      const double auxiliary,
      const double slice_width,
      const int maximum_step_out,
      const CAsisWorkspace& workspace) {

    if (
      !std::isfinite(
        current_scale2
      ) ||
      current_scale2 <= 0.0 ||
      !std::isfinite(
        degrees_freedom
      ) ||
      degrees_freedom <= 0.0 ||
      !std::isfinite(
        auxiliary
      ) ||
      auxiliary <= 0.0 ||
      !std::isfinite(
        slice_width
      ) ||
      slice_width <= 0.0 ||
      maximum_step_out < 1
    ) {
      Rcpp::stop(
        "Invalid c ASIS slice inputs"
      );
    }

    const double current_log_scale2 =
      std::log(
        current_scale2
      );

    const double current_log_density =
      c_asis_log_density(
        current_log_scale2,
        current_log_scale2,
        degrees_freedom,
        auxiliary,
        workspace
      );

    if (
      !std::isfinite(
        current_log_density
      )
    ) {
      Rcpp::stop(
        "Current c ASIS state has non-finite log density"
      );
    }

    double slice_uniform =
      R::runif(
        0.0,
        1.0
      );

    if (
      slice_uniform <= 0.0
    ) {
      slice_uniform =
        std::numeric_limits<
          double
        >::min();
    }

    const double log_slice_height =
      current_log_density +
      std::log(
        slice_uniform
      );

    const double initial_offset =
      R::runif(
        0.0,
        slice_width
      );

    double left =
      current_log_scale2 -
      initial_offset;

    double right =
      left +
      slice_width;

    int left_budget =
      static_cast<int>(
        std::floor(
          R::runif(
            0.0,
            1.0
          ) *
          static_cast<double>(
            maximum_step_out
          )
        )
      );

    if (
      left_budget >=
      maximum_step_out
    ) {
      left_budget =
        maximum_step_out -
        1;
    }

    int right_budget =
      maximum_step_out -
      1 -
      left_budget;

    int evaluations =
      1;

    int left_expansions =
      0;

    int right_expansions =
      0;

    while (
      left_budget >
      0
    ) {

      const double left_density =
        c_asis_log_density(
          left,
          current_log_scale2,
          degrees_freedom,
          auxiliary,
          workspace
        );

      ++evaluations;

      if (
        !std::isfinite(
          left_density
        ) ||
        left_density <=
          log_slice_height
      ) {
        break;
      }

      left -=
        slice_width;

      --left_budget;

      ++left_expansions;
    }

    while (
      right_budget >
      0
    ) {

      const double right_density =
        c_asis_log_density(
          right,
          current_log_scale2,
          degrees_freedom,
          auxiliary,
          workspace
        );

      ++evaluations;

      if (
        !std::isfinite(
          right_density
        ) ||
        right_density <=
          log_slice_height
      ) {
        break;
      }

      right +=
        slice_width;

      --right_budget;

      ++right_expansions;
    }

    constexpr int maximum_shrink_steps =
      10000;

    double accepted_log_scale2 =
      current_log_scale2;

    double accepted_log_density =
      current_log_density;

    int shrink_steps =
      0;

    bool accepted =
      false;

    for (
      int attempt = 0;
      attempt < maximum_shrink_steps;
      ++attempt
    ) {

      const double candidate =
        R::runif(
          left,
          right
        );

      const double candidate_density =
        c_asis_log_density(
          candidate,
          current_log_scale2,
          degrees_freedom,
          auxiliary,
          workspace
        );

      ++evaluations;

      if (
        std::isfinite(
          candidate_density
        ) &&
        candidate_density >=
          log_slice_height
      ) {

        accepted_log_scale2 =
          candidate;

        accepted_log_density =
          candidate_density;

        accepted =
          true;

        break;
      }

      ++shrink_steps;

      if (
        candidate <
        current_log_scale2
      ) {

        left =
          candidate;

      } else {

        right =
          candidate;
      }
    }

    if (!accepted) {
      Rcpp::stop(
        "c ASIS slice sampler exceeded the shrink-step limit"
      );
    }

    const double new_scale2 =
      std::exp(
        accepted_log_scale2
      );

    const double transport_scale =
      std::exp(
        0.5 *
        (
          accepted_log_scale2 -
          current_log_scale2
        )
      );

    if (
      !std::isfinite(
        new_scale2
      ) ||
      new_scale2 <= 0.0 ||
      !std::isfinite(
        transport_scale
      ) ||
      transport_scale <= 0.0
    ) {
      Rcpp::stop(
        "c ASIS slice sampler produced an invalid scale"
      );
    }

    CAsisSliceResult result;

    result.moved =
      std::abs(
        accepted_log_scale2 -
        current_log_scale2
      ) >
      1e-14;

    result.current_log_scale2 =
      current_log_scale2;

    result.new_log_scale2 =
      accepted_log_scale2;

    result.current_scale2 =
      current_scale2;

    result.new_scale2 =
      new_scale2;

    result.transport_scale =
      transport_scale;

    result.current_log_density =
      current_log_density;

    result.new_log_density =
      accepted_log_density;

    result.log_slice_height =
      log_slice_height;

    result.evaluations =
      evaluations;

    result.left_expansions =
      left_expansions;

    result.right_expansions =
      right_expansions;

    result.shrink_steps =
      shrink_steps;

    return result;
  }


  inline void apply_c_asis_transport(
      arma::mat& beta,
      arma::mat& residual,
      const CAsisGroupCache& cache,
      const CAsisWorkspace& workspace,
      const double transport_scale) {

    if (
      !std::isfinite(
        transport_scale
      ) ||
      transport_scale <= 0.0
    ) {
      Rcpp::stop(
        "Invalid c ASIS transport scale"
      );
    }

    const double transport_change =
      transport_scale -
      1.0;

    for (
      std::size_t cell = 0u;
      cell < cache.equations.size();
      ++cell
    ) {

      const CAsisEquationCache&
        equation_cache =
          cache.equations[
            cell
          ];

      const arma::uword equation =
        equation_cache.response_column;

      for (
        arma::uword position = 0u;
        position < equation_cache.rows.n_elem;
        ++position
      ) {

        beta(
          equation_cache.rows[
            position
          ],
          equation
        ) *=
          transport_scale;
      }

      residual.col(
        equation
      ) -=
        transport_change *
        workspace.fitted_values[
          cell
        ];
    }
  }


  inline CAsisSliceResult update_c_group_asis(
      arma::mat& beta,
      arma::mat& residual,
      const arma::mat& X,
      const arma::vec& sigma2,
      double& scale2,
      const double auxiliary,
      const double degrees_freedom,
      const CAsisGroupCache& cache,
      const double slice_width,
      const int maximum_step_out,
      CAsisWorkspace& workspace) {

    if (
      beta.n_rows != X.n_cols ||
      beta.n_cols != residual.n_cols ||
      residual.n_rows != X.n_rows ||
      sigma2.n_elem != beta.n_cols
    ) {
      Rcpp::stop(
        "c ASIS dimensions are inconsistent"
      );
    }

    if (
      !beta.is_finite() ||
      !residual.is_finite() ||
      !X.is_finite() ||
      !sigma2.is_finite() ||
      arma::any(
        sigma2 <= 0.0
      )
    ) {
      Rcpp::stop(
        "c ASIS state contains invalid values"
      );
    }

    prepare_c_asis_workspace(
      beta,
      residual,
      X,
      sigma2,
      cache,
      workspace
    );

    const CAsisSliceResult result =
      draw_c_asis_log_scale2(
        scale2,
        degrees_freedom,
        auxiliary,
        slice_width,
        maximum_step_out,
        workspace
      );

    apply_c_asis_transport(
      beta,
      residual,
      cache,
      workspace,
      result.transport_scale
    );

    scale2 =
      result.new_scale2;

    return result;
  }


  // Discrete q Gibbs -----

  struct QLogWeightResult {

    arma::vec log_weights;

    double n_diagonal;
    double n_off_diagonal;

    double diagonal_quadratic;
    double off_diagonal_quadratic;
  };


  struct QGibbsGroupCache {

    arma::uvec diagonal_linear_indices;
    arma::uvec off_diagonal_linear_indices;

    arma::uword n_rows;
  };


  inline std::vector<QGibbsGroupCache>
  build_q_gibbs_cache(
      const arma::imat& group_id,
      const arma::umat& same_var,
      const int n_groups) {

    if (
      group_id.n_rows !=
        same_var.n_rows ||
      group_id.n_cols !=
        same_var.n_cols
    ) {
      Rcpp::stop(
        "group_id and same_var must have equal dimensions"
      );
    }

    if (n_groups < 1) {
      Rcpp::stop(
        "n_groups must be positive"
      );
    }

    const arma::uword k =
      group_id.n_rows;

    const arma::uword N =
      group_id.n_cols;

    std::vector<
      std::vector<arma::uword>
    > diagonal(
      static_cast<std::size_t>(
        n_groups
      )
    );

    std::vector<
      std::vector<arma::uword>
    > off_diagonal(
      static_cast<std::size_t>(
        n_groups
      )
    );

    for (
      arma::uword i = 0u;
      i < N;
      ++i
    ) {

      for (
        arma::uword r = 0u;
        r < k;
        ++r
      ) {

        const int group =
          group_id(r, i);

        const arma::uword same =
          same_var(r, i);

        if (
          group < 1 ||
          group > n_groups
        ) {
          Rcpp::stop(
            "q group_id contains a value outside 1:n_groups"
          );
        }

        if (
          same != 0u &&
          same != 1u
        ) {
          Rcpp::stop(
            "same_var must contain only zero or one"
          );
        }

        const arma::uword linear =
          r +
          i *
          k;

        if (same == 1u) {

          diagonal[
            static_cast<std::size_t>(
              group -
              1
            )
          ].push_back(
            linear
          );

        } else {

          off_diagonal[
            static_cast<std::size_t>(
              group -
              1
            )
          ].push_back(
            linear
          );
        }
      }
    }

    std::vector<QGibbsGroupCache> cache(
      static_cast<std::size_t>(
        n_groups
      )
    );

    for (
      int group = 0;
      group < n_groups;
      ++group
    ) {

      QGibbsGroupCache& group_cache =
        cache[
          static_cast<std::size_t>(
            group
          )
        ];

      group_cache.n_rows =
        k;

      group_cache.diagonal_linear_indices =
        uvec_from_std(
          diagonal[
            static_cast<std::size_t>(
              group
            )
          ]
        );

      group_cache.off_diagonal_linear_indices =
        uvec_from_std(
          off_diagonal[
            static_cast<std::size_t>(
              group
            )
          ]
        );

      if (
        group_cache.diagonal_linear_indices.n_elem +
        group_cache.off_diagonal_linear_indices.n_elem ==
        0u
      ) {
        Rcpp::stop(
          "Every q group must contain coefficients"
        );
      }
    }

    return cache;
  }


  inline QLogWeightResult compute_q_log_weights_cached(
      const arma::mat& beta,
      const arma::vec& sigma2,
      const arma::mat& base_omega,
      const QGibbsGroupCache& cache,
      const int m,
      const arma::vec& q_grid,
      const arma::vec& q_prob) {

    double diagonal_quadratic =
      0.0;

    double off_diagonal_quadratic =
      0.0;

    const arma::uword k =
      cache.n_rows;

    for (
      arma::uword position = 0u;
      position <
        cache.diagonal_linear_indices.n_elem;
      ++position
    ) {

      const arma::uword linear =
        cache.diagonal_linear_indices[
          position
        ];

      const arma::uword equation =
        linear /
        k;

      const double coefficient =
        beta[
          linear
        ];

      diagonal_quadratic +=
        coefficient *
        coefficient /
        (
          sigma2[equation] *
          base_omega[linear]
        );
    }

    for (
      arma::uword position = 0u;
      position <
        cache.off_diagonal_linear_indices.n_elem;
      ++position
    ) {

      const arma::uword linear =
        cache.off_diagonal_linear_indices[
          position
        ];

      const arma::uword equation =
        linear /
        k;

      const double coefficient =
        beta[
          linear
        ];

      off_diagonal_quadratic +=
        coefficient *
        coefficient /
        (
          sigma2[equation] *
          base_omega[linear]
        );
    }

    const double n_diagonal =
      static_cast<double>(
        cache.diagonal_linear_indices.n_elem
      );

    const double n_off_diagonal =
      static_cast<double>(
        cache.off_diagonal_linear_indices.n_elem
      );

    const double m_double =
      static_cast<double>(
        m
      );

    const double q_probability_sum =
      arma::accu(
        q_prob
      );

    arma::vec log_weights(
      q_grid.n_elem
    );

    for (
      arma::uword candidate = 0u;
      candidate < q_grid.n_elem;
      ++candidate
    ) {

      const double q =
        q_grid[
          candidate
        ];

      const double log_r_diagonal =
        (
          static_cast<double>(
            m -
            1
          ) /
          m_double
        ) *
        q;

      const double log_r_off_diagonal =
        -q /
        m_double;

      const double determinant_term =
        n_diagonal *
        log_r_diagonal +
        n_off_diagonal *
        log_r_off_diagonal;

      const double quadratic_term =
        diagonal_quadratic *
        std::exp(
          -log_r_diagonal
        ) +
        off_diagonal_quadratic *
        std::exp(
          -log_r_off_diagonal
        );

      const double prior_probability =
        q_prob[
          candidate
        ] /
        q_probability_sum;

      const double log_prior =
        prior_probability > 0.0
          ? std::log(
              prior_probability
            )
          : -std::numeric_limits<
              double
            >::infinity();

      log_weights[
        candidate
      ] =
        log_prior -
        0.5 *
        (
          determinant_term +
          quadratic_term
        );
    }

    QLogWeightResult result;

    result.log_weights =
      log_weights;

    result.n_diagonal =
      n_diagonal;

    result.n_off_diagonal =
      n_off_diagonal;

    result.diagonal_quadratic =
      diagonal_quadratic;

    result.off_diagonal_quadratic =
      off_diagonal_quadratic;

    return result;
  }


  inline arma::uword draw_q_index_cached(
      const arma::mat& beta,
      const arma::vec& sigma2,
      const arma::mat& base_omega,
      const QGibbsGroupCache& cache,
      const int m,
      const arma::vec& q_grid,
      const arma::vec& q_prob) {

    const QLogWeightResult result =
      compute_q_log_weights_cached(
        beta,
        sigma2,
        base_omega,
        cache,
        m,
        q_grid,
        q_prob
      );

    return sample_log_weights(
      result.log_weights
    );
  }


  inline QLogWeightResult compute_q_log_weights(
      const arma::mat& beta,
      const arma::vec& sigma2,
      const arma::mat& base_omega,
      const arma::imat& group_id,
      const arma::umat& same_var,
      const int target_group,
      const int m,
      const arma::vec& q_grid,
      const arma::vec& q_prob) {

    if (
      target_group < 1
    ) {
      Rcpp::stop(
        "target_group must be positive"
      );
    }

    const int n_groups =
      group_id.max();

    const std::vector<QGibbsGroupCache> cache =
      build_q_gibbs_cache(
        group_id,
        same_var,
        n_groups
      );

    if (
      target_group >
      n_groups
    ) {
      Rcpp::stop(
        "target_group lies outside group_id"
      );
    }

    return compute_q_log_weights_cached(
      beta,
      sigma2,
      base_omega,
      cache[
        static_cast<std::size_t>(
          target_group -
          1
        )
      ],
      m,
      q_grid,
      q_prob
    );
  }


  inline arma::uword draw_q_index(
      const arma::mat& beta,
      const arma::vec& sigma2,
      const arma::mat& base_omega,
      const arma::imat& group_id,
      const arma::umat& same_var,
      const int target_group,
      const int m,
      const arma::vec& q_grid,
      const arma::vec& q_prob) {

    const QLogWeightResult result =
      compute_q_log_weights(
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

    return sample_log_weights(
      result.log_weights
    );
  }


  // Joint q + beta transport MH -----

  struct QEquationTransportCache {

    arma::uword response_column;

    arma::uvec diagonal_rows;
    arma::uvec off_diagonal_rows;
  };


  struct QGroupTransportCache {

    std::vector<QEquationTransportCache>
      equations;

    arma::uword n_diagonal;
    arma::uword n_off_diagonal;
  };


  struct QIndexProposal {

    arma::uword proposed_index;

    bool used_global;

    double log_proposal_ratio;

    arma::uword absolute_jump;
  };


  struct QTransportWorkspace {

    arma::vec fitted_change;
    arma::vec diagonal_coefficient_change;
    arma::vec off_diagonal_coefficient_change;

    std::vector<arma::vec>
      fitted_changes;

    explicit QTransportWorkspace(
        const arma::uword T_p) :

      fitted_change(
        T_p,
        arma::fill::zeros
      ),

      diagonal_coefficient_change(),

      off_diagonal_coefficient_change(),

      fitted_changes() {}
  };


  struct QTransportEvaluation {

    double diagonal_scale;
    double off_diagonal_scale;

    double log_likelihood_ratio;
  };


  struct QTransportMHResult {

    bool accepted;
    bool used_global;

    arma::uword current_index;
    arma::uword proposed_index;
    arma::uword new_index;
    arma::uword absolute_jump;

    double diagonal_scale;
    double off_diagonal_scale;

    double log_likelihood_ratio;
    double log_q_prior_ratio;
    double log_proposal_ratio;
    double log_acceptance_ratio;
    double acceptance_probability;
  };


  inline std::vector<QGroupTransportCache>
  build_q_transport_cache(
      const arma::imat& group_id,
      const arma::umat& same_var,
      const int n_groups) {

    if (
      group_id.n_rows !=
        same_var.n_rows ||
      group_id.n_cols !=
        same_var.n_cols
    ) {
      Rcpp::stop(
        "group_id and same_var must have equal dimensions"
      );
    }

    if (n_groups < 1) {
      Rcpp::stop(
        "n_groups must be positive"
      );
    }

    const arma::uword k =
      group_id.n_rows;

    const arma::uword N =
      group_id.n_cols;

    const std::size_t cell_count =
      static_cast<std::size_t>(
        n_groups
      ) *
      static_cast<std::size_t>(
        N
      );

    std::vector<
      std::vector<arma::uword>
    > diagonal_rows(
      cell_count
    );

    std::vector<
      std::vector<arma::uword>
    > off_diagonal_rows(
      cell_count
    );

    for (
      arma::uword i = 0u;
      i < N;
      ++i
    ) {

      for (
        arma::uword r = 0u;
        r < k;
        ++r
      ) {

        const int group =
          group_id(r, i);

        const arma::uword same =
          same_var(r, i);

        if (
          group < 1 ||
          group > n_groups
        ) {
          Rcpp::stop(
            "q transport group_id lies outside 1:n_groups"
          );
        }

        if (
          same != 0u &&
          same != 1u
        ) {
          Rcpp::stop(
            "same_var must contain only zero or one"
          );
        }

        const std::size_t cell =
          static_cast<std::size_t>(
            group -
            1
          ) *
          static_cast<std::size_t>(
            N
          ) +
          static_cast<std::size_t>(
            i
          );

        if (same == 1u) {

          diagonal_rows[
            cell
          ].push_back(
            r
          );

        } else {

          off_diagonal_rows[
            cell
          ].push_back(
            r
          );
        }
      }
    }

    std::vector<QGroupTransportCache> cache(
      static_cast<std::size_t>(
        n_groups
      )
    );

    for (
      int group = 0;
      group < n_groups;
      ++group
    ) {

      QGroupTransportCache& group_cache =
        cache[
          static_cast<std::size_t>(
            group
          )
        ];

      group_cache.n_diagonal =
        0u;

      group_cache.n_off_diagonal =
        0u;

      for (
        arma::uword i = 0u;
        i < N;
        ++i
      ) {

        const std::size_t cell =
          static_cast<std::size_t>(
            group
          ) *
          static_cast<std::size_t>(
            N
          ) +
          static_cast<std::size_t>(
            i
          );

        if (
          diagonal_rows[
            cell
          ].empty() &&
          off_diagonal_rows[
            cell
          ].empty()
        ) {
          continue;
        }

        QEquationTransportCache
          equation_cache;

        equation_cache.response_column =
          i;

        equation_cache.diagonal_rows =
          uvec_from_std(
            diagonal_rows[
              cell
            ]
          );

        equation_cache.off_diagonal_rows =
          uvec_from_std(
            off_diagonal_rows[
              cell
            ]
          );

        group_cache.n_diagonal +=
          equation_cache
            .diagonal_rows
            .n_elem;

        group_cache.n_off_diagonal +=
          equation_cache
            .off_diagonal_rows
            .n_elem;

        group_cache.equations.push_back(
          equation_cache
        );
      }

      if (
        group_cache.n_diagonal +
        group_cache.n_off_diagonal ==
        0u
      ) {
        Rcpp::stop(
          "Every q transport group must contain coefficients"
        );
      }
    }

    return cache;
  }


  inline arma::uword q_transport_local_neighbor_count(
      const arma::uword index,
      const arma::uword n_grid) {

    return (
      index == 0u ||
      index + 1u ==
        n_grid
    ) ?
      1u :
      2u;
  }


  inline QIndexProposal propose_q_transport_index(
      const arma::uword current_index,
      const arma::uword n_grid,
      const double global_probability) {

    QIndexProposal proposal;

    proposal.used_global =
      R::runif(
        0.0,
        1.0
      ) <
      global_probability;

    if (proposal.used_global) {

      arma::uword draw =
        static_cast<arma::uword>(
          std::floor(
            R::runif(
              0.0,
              static_cast<double>(
                n_grid -
                1u
              )
            )
          )
        );

      if (
        draw >=
        n_grid -
        1u
      ) {
        draw =
          n_grid -
          2u;
      }

      proposal.proposed_index =
        draw >=
          current_index
          ? draw +
            1u
          : draw;

    } else {

      if (current_index == 0u) {

        proposal.proposed_index =
          1u;

      } else if (
        current_index +
        1u ==
        n_grid
      ) {

        proposal.proposed_index =
          n_grid -
          2u;

      } else {

        proposal.proposed_index =
          R::runif(
            0.0,
            1.0
          ) <
          0.5
            ? current_index -
              1u
            : current_index +
              1u;
      }
    }

    const arma::uword forward_neighbors =
      q_transport_local_neighbor_count(
        current_index,
        n_grid
      );

    const arma::uword reverse_neighbors =
      q_transport_local_neighbor_count(
        proposal.proposed_index,
        n_grid
      );

    const bool adjacent =
      proposal.proposed_index +
        1u ==
        current_index ||
      current_index +
        1u ==
        proposal.proposed_index;

    const double global_component =
      global_probability /
      static_cast<double>(
        n_grid -
        1u
      );

    const double forward_probability =
      global_component +
      (
        adjacent
          ? (
              1.0 -
              global_probability
            ) /
            static_cast<double>(
              forward_neighbors
            )
          : 0.0
      );

    const double reverse_probability =
      global_component +
      (
        adjacent
          ? (
              1.0 -
              global_probability
            ) /
            static_cast<double>(
              reverse_neighbors
            )
          : 0.0
      );

    proposal.log_proposal_ratio =
      std::log(
        reverse_probability
      ) -
      std::log(
        forward_probability
      );

    proposal.absolute_jump =
      proposal.proposed_index >
        current_index
        ? proposal.proposed_index -
          current_index
        : current_index -
          proposal.proposed_index;

    return proposal;
  }


  inline QTransportEvaluation evaluate_q_transport(
      const arma::mat& beta,
      const arma::mat& residual,
      const arma::mat& X,
      const arma::vec& sigma2,
      const QGroupTransportCache& cache,
      const double current_q,
      const double proposed_q,
      const int m,
      QTransportWorkspace& workspace) {

    const double m_double =
      static_cast<double>(
        m
      );

    const double delta_q =
      proposed_q -
      current_q;

    QTransportEvaluation evaluation;

    evaluation.diagonal_scale =
      std::exp(
        0.5 *
        (
          static_cast<double>(
            m -
            1
          ) /
          m_double
        ) *
        delta_q
      );

    evaluation.off_diagonal_scale =
      std::exp(
        -0.5 *
        delta_q /
        m_double
      );

    evaluation.log_likelihood_ratio =
      0.0;

    workspace.fitted_changes.resize(
      cache.equations.size()
    );

    for (
      std::size_t cell = 0u;
      cell < cache.equations.size();
      ++cell
    ) {

      const QEquationTransportCache&
        equation_cache =
          cache.equations[
            cell
          ];

      const arma::uword equation =
        equation_cache
          .response_column;

      const arma::vec beta_column =
        beta.col(
          equation
        );

      workspace.fitted_change.zeros();

      if (
        equation_cache
          .diagonal_rows
          .n_elem >
        0u
      ) {

        workspace
          .diagonal_coefficient_change =
            (
              evaluation.diagonal_scale -
              1.0
            ) *
            beta_column.elem(
              equation_cache
                .diagonal_rows
            );

        workspace.fitted_change +=
          X.cols(
            equation_cache
              .diagonal_rows
          ) *
          workspace
            .diagonal_coefficient_change;
      }

      if (
        equation_cache
          .off_diagonal_rows
          .n_elem >
        0u
      ) {

        workspace
          .off_diagonal_coefficient_change =
            (
              evaluation.off_diagonal_scale -
              1.0
            ) *
            beta_column.elem(
              equation_cache
                .off_diagonal_rows
            );

        workspace.fitted_change +=
          X.cols(
            equation_cache
              .off_diagonal_rows
          ) *
          workspace
            .off_diagonal_coefficient_change;
      }

      const double rss_change =
        arma::dot(
          workspace.fitted_change,
          workspace.fitted_change
        ) -
        2.0 *
        arma::dot(
          residual.col(
            equation
          ),
          workspace.fitted_change
        );

      evaluation.log_likelihood_ratio +=
        -0.5 *
        rss_change /
        sigma2[
          equation
        ];

      workspace.fitted_changes[
        cell
      ] =
        workspace.fitted_change;
    }

    return evaluation;
  }


  inline void apply_q_transport(
      arma::mat& beta,
      arma::mat& residual,
      const QGroupTransportCache& cache,
      const QTransportEvaluation& evaluation,
      const QTransportWorkspace& workspace) {

    for (
      std::size_t cell = 0u;
      cell < cache.equations.size();
      ++cell
    ) {

      const QEquationTransportCache&
        equation_cache =
          cache.equations[
            cell
          ];

      const arma::uword equation =
        equation_cache
          .response_column;

      for (
        arma::uword position = 0u;
        position <
          equation_cache
            .diagonal_rows
            .n_elem;
        ++position
      ) {

        beta(
          equation_cache
            .diagonal_rows[
              position
            ],
          equation
        ) *=
          evaluation
            .diagonal_scale;
      }

      for (
        arma::uword position = 0u;
        position <
          equation_cache
            .off_diagonal_rows
            .n_elem;
        ++position
      ) {

        beta(
          equation_cache
            .off_diagonal_rows[
              position
            ],
          equation
        ) *=
          evaluation
            .off_diagonal_scale;
      }

      residual.col(
        equation
      ) -=
        workspace
          .fitted_changes[
            cell
          ];
    }
  }


  inline QTransportMHResult update_q_transport_mh(
      arma::mat& beta,
      arma::mat& residual,
      const arma::mat& X,
      const arma::vec& sigma2,
      const QGroupTransportCache& cache,
      const arma::uword current_index,
      const int m,
      const arma::vec& q_grid,
      const arma::vec& q_prob,
      const double global_probability,
      QTransportWorkspace& workspace) {

    const QIndexProposal proposal =
      propose_q_transport_index(
        current_index,
        q_grid.n_elem,
        global_probability
      );

    const QTransportEvaluation evaluation =
      evaluate_q_transport(
        beta,
        residual,
        X,
        sigma2,
        cache,
        q_grid[
          current_index
        ],
        q_grid[
          proposal.proposed_index
        ],
        m,
        workspace
      );

    const double log_q_prior_ratio =
      q_prob[
        proposal.proposed_index
      ] >
      0.0
        ? std::log(
            q_prob[
              proposal.proposed_index
            ]
          ) -
          std::log(
            q_prob[
              current_index
            ]
          )
        : -std::numeric_limits<
            double
          >::infinity();

    const double log_acceptance_ratio =
      evaluation
        .log_likelihood_ratio +
      log_q_prior_ratio +
      proposal
        .log_proposal_ratio;

    if (
      std::isnan(
        log_acceptance_ratio
      )
    ) {
      Rcpp::stop(
        "NaN acceptance ratio in q transport"
      );
    }

    const double acceptance_probability =
      log_acceptance_ratio >=
        0.0
        ? 1.0
        : (
            std::isfinite(
              log_acceptance_ratio
            )
              ? std::exp(
                  log_acceptance_ratio
                )
              : 0.0
          );

    const bool accepted =
      R::runif(
        0.0,
        1.0
      ) <=
      acceptance_probability;

    if (accepted) {

      apply_q_transport(
        beta,
        residual,
        cache,
        evaluation,
        workspace
      );
    }

    QTransportMHResult result;

    result.accepted =
      accepted;

    result.used_global =
      proposal.used_global;

    result.current_index =
      current_index;

    result.proposed_index =
      proposal.proposed_index;

    result.new_index =
      accepted
        ? proposal.proposed_index
        : current_index;

    result.absolute_jump =
      proposal.absolute_jump;

    result.diagonal_scale =
      evaluation.diagonal_scale;

    result.off_diagonal_scale =
      evaluation.off_diagonal_scale;

    result.log_likelihood_ratio =
      evaluation.log_likelihood_ratio;

    result.log_q_prior_ratio =
      log_q_prior_ratio;

    result.log_proposal_ratio =
      proposal.log_proposal_ratio;

    result.log_acceptance_ratio =
      log_acceptance_ratio;

    result.acceptance_probability =
      acceptance_probability;

    return result;
  }


}  // namespace m3_bvar

#endif