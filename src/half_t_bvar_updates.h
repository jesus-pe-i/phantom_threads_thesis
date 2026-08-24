#ifndef HALF_T_BVAR_UPDATES_H
#define HALF_T_BVAR_UPDATES_H

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>


// Half-t BVAR conditional-update layer.
// Production kernels, compatibility wrappers, and non-self tau ASIS.

namespace half_t_bvar {


  // Helpers -----

  constexpr int grouping_all = 0;
  constexpr int grouping_self_diagonal = 1;

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


  // Prior variance construction -----

  inline void fill_omega(
      arma::vec& omega,
      const arma::vec& phi2,
      const arma::mat& lambda2,
      const arma::vec& tau2,
      const arma::uvec& predictor_series,
      const arma::uword equation,
      const int grouping_code) {

    const arma::uword k =
      phi2.n_elem;

    if (
      grouping_code ==
      grouping_all
    ) {

      const double tau2_all =
        tau2[0u];

      for (
        arma::uword r = 0u;
        r < k;
        ++r
      ) {

        const double omega_value =
          phi2[r] *
          lambda2(r, equation) *
          tau2_all;

        if (
          !std::isfinite(omega_value) ||
          omega_value <= 0.0
        ) {
          Rcpp::stop(
            "Invalid coefficient prior variance multiplier"
          );
        }

        omega[r] =
          omega_value;
      }

      return;
    }

    if (
      grouping_code ==
      grouping_self_diagonal
    ) {

      for (
        arma::uword r = 0u;
        r < k;
        ++r
      ) {

        const arma::uword group =
          predictor_series[r] ==
          equation ?
            0u :
            1u;

        const double omega_value =
          phi2[r] *
          lambda2(r, equation) *
          tau2[group];

        if (
          !std::isfinite(omega_value) ||
          omega_value <= 0.0
        ) {
          Rcpp::stop(
            "Invalid coefficient prior variance multiplier"
          );
        }

        omega[r] =
          omega_value;
      }

      return;
    }

    Rcpp::stop(
      "Unknown global grouping code"
    );
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

    if (!beta.is_finite()) {
      Rcpp::stop(
        "Invalid Cholesky beta draw"
      );
    }
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

    if (!beta.is_finite()) {
      Rcpp::stop(
        "Invalid Bhattacharya beta draw"
      );
    }
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
      !std::isfinite(sigma2) ||
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
      !std::isfinite(sigma2) ||
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

  inline double draw_sigma2_from_sufficient_statistics(
      const double Yty,
      const arma::vec& Xty,
      const arma::mat& XtX,
      const arma::vec& beta,
      const arma::vec& omega,
      const double prior_shape,
      const double prior_scale,
      const arma::uword T_p,
      arma::vec& XtX_beta) {

    XtX_beta =
      XtX *
      beta;

    double rss =
      Yty -
      2.0 *
      arma::dot(
        beta,
        Xty
      ) +
      arma::dot(
        beta,
        XtX_beta
      );

    const double rss_tolerance =
      1e-10 *
      std::max(
        1.0,
        std::abs(
          Yty
        )
      );

    if (
      rss < 0.0 &&
      rss >=
        -rss_tolerance
    ) {
      rss = 0.0;
    }

    if (
      !std::isfinite(rss) ||
      rss < 0.0
    ) {
      Rcpp::stop(
        "Invalid residual sum of squares"
      );
    }

    const double beta_quadratic =
      arma::dot(
        beta % beta,
        1.0 /
        omega
      );

    if (
      !std::isfinite(beta_quadratic) ||
      beta_quadratic < 0.0
    ) {
      Rcpp::stop(
        "Invalid coefficient-prior quadratic form"
      );
    }

    const double posterior_shape =
      prior_shape +
      0.5 *
      static_cast<double>(
        T_p +
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
      !std::isfinite(prior_shape) ||
      prior_shape <= 0.0 ||
      !std::isfinite(prior_scale) ||
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
        y.n_elem +
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


  // Local and global half-t hierarchy -----

  inline void update_local_hierarchy(
      arma::mat& lambda2,
      arma::mat& nu,
      arma::vec& tau_quadratic,
      const arma::mat& beta,
      const arma::vec& sigma2,
      const arma::vec& phi2,
      const arma::uvec& predictor_series,
      const arma::vec& tau2,
      const int grouping_code,
      const double lambda_df,
      const double lambda_scale) {

    const arma::uword k =
      beta.n_rows;

    const arma::uword N =
      beta.n_cols;

    const double conditional_shape =
      0.5 *
      (
        lambda_df +
        1.0
      );

    tau_quadratic.zeros();

    if (
      grouping_code ==
      grouping_all
    ) {

      const double tau2_all =
        tau2[0u];

      for (
        arma::uword equation = 0u;
        equation < N;
        ++equation
      ) {

        const double sigma2_i =
          sigma2[equation];

        for (
          arma::uword r = 0u;
          r < k;
          ++r
        ) {

          const double beta_value =
            beta(r, equation);

          const double beta2 =
            beta_value *
            beta_value;

          const double lambda2_new =
            rinvgamma(
              conditional_shape,
              lambda_df /
              nu(r, equation) +
              beta2 /
              (
                2.0 *
                sigma2_i *
                phi2[r] *
                tau2_all
              )
            );

          lambda2(r, equation) =
            lambda2_new;

          nu(r, equation) =
            draw_half_t_auxiliary(
              lambda2_new,
              lambda_df,
              lambda_scale
            );

          tau_quadratic[0u] +=
            beta2 /
            (
              sigma2_i *
              phi2[r] *
              lambda2_new
            );
        }
      }

    } else if (
      grouping_code ==
      grouping_self_diagonal
    ) {

      for (
        arma::uword equation = 0u;
        equation < N;
        ++equation
      ) {

        const double sigma2_i =
          sigma2[equation];

        for (
          arma::uword r = 0u;
          r < k;
          ++r
        ) {

          const arma::uword group =
            predictor_series[r] ==
            equation ?
              0u :
              1u;

          const double beta_value =
            beta(r, equation);

          const double beta2 =
            beta_value *
            beta_value;

          const double lambda2_new =
            rinvgamma(
              conditional_shape,
              lambda_df /
              nu(r, equation) +
              beta2 /
              (
                2.0 *
                sigma2_i *
                phi2[r] *
                tau2[group]
              )
            );

          lambda2(r, equation) =
            lambda2_new;

          nu(r, equation) =
            draw_half_t_auxiliary(
              lambda2_new,
              lambda_df,
              lambda_scale
            );

          tau_quadratic[group] +=
            beta2 /
            (
              sigma2_i *
              phi2[r] *
              lambda2_new
            );
        }
      }

    } else {
      Rcpp::stop(
        "Unknown global grouping code"
      );
    }

    if (
      !tau_quadratic.is_finite() ||
      arma::any(
        tau_quadratic < 0.0
      )
    ) {
      Rcpp::stop(
        "Invalid global-scale quadratic terms"
      );
    }
  }


  inline void update_global_hierarchy(
      arma::vec& tau2,
      arma::vec& xi,
      const arma::vec& tau_quadratic,
      const arma::uvec& tau_group_counts,
      const double tau_df,
      const double tau_scale) {

    const arma::uword n_tau =
      tau2.n_elem;

    for (
      arma::uword group = 0u;
      group < n_tau;
      ++group
    ) {

      const double posterior_shape =
        0.5 *
        (
          tau_df +
          static_cast<double>(
            tau_group_counts[group]
          )
        );

      const double posterior_scale =
        tau_df /
        xi[group] +
        0.5 *
        tau_quadratic[group];

      tau2[group] =
        rinvgamma(
          posterior_shape,
          posterior_scale
        );

      xi[group] =
        draw_half_t_auxiliary(
          tau2[group],
          tau_df,
          tau_scale
        );
    }
  }


  // Non-self tau ASIS -----

  struct AsisWorkspace {

    arma::vec foreign_signal;
    arma::vec residual_without_foreign;

    explicit AsisWorkspace(
        const arma::uword T_p) :

      foreign_signal(
        T_p,
        arma::fill::zeros
      ),

      residual_without_foreign(
        T_p,
        arma::fill::zeros
      ) {}
  };


  struct AsisSliceResult {

    double eta;

    int evaluations;
    int left_expansions;
    int right_expansions;
    int shrink_steps;
  };


  struct AsisResult {

    double tau2_before;
    double tau2_after;

    double xi_before;
    double xi_after;

    double likelihood_quadratic;
    double likelihood_linear;

    double log_tau_move;
    double maximum_relative_invariant_error;

    int evaluations;
    int left_expansions;
    int right_expansions;
    int shrink_steps;

    AsisResult() :

      tau2_before(
        NA_REAL
      ),

      tau2_after(
        NA_REAL
      ),

      xi_before(
        NA_REAL
      ),

      xi_after(
        NA_REAL
      ),

      likelihood_quadratic(
        NA_REAL
      ),

      likelihood_linear(
        NA_REAL
      ),

      log_tau_move(
        NA_REAL
      ),

      maximum_relative_invariant_error(
        NA_REAL
      ),

      evaluations(
        0
      ),

      left_expansions(
        0
      ),

      right_expansions(
        0
      ),

      shrink_steps(
        0
      ) {}
  };


  inline double nonself_asis_log_density(
      const double eta,
      const double likelihood_quadratic,
      const double likelihood_linear,
      const double xi_foreign,
      const double tau_df) {

    if (
      !std::isfinite(eta) ||
      !std::isfinite(
        likelihood_quadratic
      ) ||
      likelihood_quadratic < 0.0 ||
      !std::isfinite(
        likelihood_linear
      ) ||
      !std::isfinite(
        xi_foreign
      ) ||
      xi_foreign <= 0.0 ||
      !std::isfinite(
        tau_df
      ) ||
      tau_df <= 0.0
    ) {
      return
        -std::numeric_limits<
          double
        >::infinity();
    }

    if (
      eta < -300.0 ||
      eta > 300.0
    ) {
      return
        -std::numeric_limits<
          double
        >::infinity();
    }

    const double tau =
      std::exp(
        eta
      );

    const double tau2 =
      tau *
      tau;

    const double inverse_tau2 =
      std::exp(
        -2.0 *
        eta
      );

    const double log_density =
      -0.5 *
      likelihood_quadratic *
      tau2 +
      likelihood_linear *
      tau -
      tau_df *
      eta -
      (
        tau_df /
        xi_foreign
      ) *
      inverse_tau2;

    if (
      !std::isfinite(
        log_density
      )
    ) {
      return
        -std::numeric_limits<
          double
        >::infinity();
    }

    return log_density;
  }


  inline AsisSliceResult draw_nonself_asis_eta(
      const double current_eta,
      const double likelihood_quadratic,
      const double likelihood_linear,
      const double xi_foreign,
      const double tau_df,
      const double slice_width,
      const int maximum_step_out,
      const int maximum_shrinks) {

    if (
      !std::isfinite(
        slice_width
      ) ||
      slice_width <= 0.0 ||
      maximum_step_out < 1 ||
      maximum_shrinks < 1
    ) {
      Rcpp::stop(
        "Invalid non-self ASIS slice controls"
      );
    }

    const double current_log_density =
      nonself_asis_log_density(
        current_eta,
        likelihood_quadratic,
        likelihood_linear,
        xi_foreign,
        tau_df
      );

    if (
      !std::isfinite(
        current_log_density
      )
    ) {
      Rcpp::stop(
        "Current non-self ASIS state has non-finite log density"
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
      current_eta -
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
      left_budget > 0
    ) {

      const double left_density =
        nonself_asis_log_density(
          left,
          likelihood_quadratic,
          likelihood_linear,
          xi_foreign,
          tau_df
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
      right_budget > 0
    ) {

      const double right_density =
        nonself_asis_log_density(
          right,
          likelihood_quadratic,
          likelihood_linear,
          xi_foreign,
          tau_df
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

    for (
      int shrink = 0;
      shrink < maximum_shrinks;
      ++shrink
    ) {

      const double proposal =
        R::runif(
          left,
          right
        );

      const double proposal_density =
        nonself_asis_log_density(
          proposal,
          likelihood_quadratic,
          likelihood_linear,
          xi_foreign,
          tau_df
        );

      ++evaluations;

      if (
        std::isfinite(
          proposal_density
        ) &&
        proposal_density >=
          log_slice_height
      ) {

        return AsisSliceResult{
          proposal,
          evaluations,
          left_expansions,
          right_expansions,
          shrink
        };
      }

      if (
        proposal <
        current_eta
      ) {
        left =
          proposal;

      } else {
        right =
          proposal;
      }
    }

    Rcpp::stop(
      "Non-self ASIS slice sampler exceeded the shrink limit"
    );

    return AsisSliceResult{
      current_eta,
      evaluations,
      left_expansions,
      right_expansions,
      maximum_shrinks
    };
  }


  inline void compute_nonself_asis_likelihood(
      double& likelihood_quadratic,
      double& likelihood_linear,
      const arma::mat& beta,
      const arma::mat& Y,
      const arma::mat& X,
      const arma::vec& sigma2,
      const arma::uvec& predictor_series,
      const double tau_foreign,
      AsisWorkspace& workspace) {

    const arma::uword T_p =
      X.n_rows;

    const arma::uword k =
      X.n_cols;

    const arma::uword N =
      Y.n_cols;

    if (
      Y.n_rows != T_p ||
      beta.n_rows != k ||
      beta.n_cols != N ||
      sigma2.n_elem != N ||
      predictor_series.n_elem != k ||
      workspace.foreign_signal.n_elem !=
        T_p ||
      workspace.residual_without_foreign.n_elem !=
        T_p
    ) {
      Rcpp::stop(
        "Invalid dimensions in the non-self ASIS likelihood"
      );
    }

    if (
      !std::isfinite(
        tau_foreign
      ) ||
      tau_foreign <= 0.0
    ) {
      Rcpp::stop(
        "Non-self ASIS tau must be positive"
      );
    }

    likelihood_quadratic =
      0.0;

    likelihood_linear =
      0.0;

    for (
      arma::uword equation = 0u;
      equation < N;
      ++equation
    ) {

      const double sigma2_i =
        sigma2[equation];

      if (
        !std::isfinite(
          sigma2_i
        ) ||
        sigma2_i <= 0.0
      ) {
        Rcpp::stop(
          "Invalid innovation variance in the non-self ASIS update"
        );
      }

      workspace.foreign_signal.zeros();

      workspace.residual_without_foreign =
        Y.col(
          equation
        );

      for (
        arma::uword r = 0u;
        r < k;
        ++r
      ) {

        const double beta_value =
          beta(
            r,
            equation
          );

        if (
          predictor_series[r] ==
          equation
        ) {

          workspace.residual_without_foreign -=
            X.col(
              r
            ) *
            beta_value;

        } else {

          workspace.foreign_signal +=
            X.col(
              r
            ) *
            (
              beta_value /
              tau_foreign
            );
        }
      }

      const double inverse_sigma2 =
        1.0 /
        sigma2_i;

      likelihood_quadratic +=
        arma::dot(
          workspace.foreign_signal,
          workspace.foreign_signal
        ) *
        inverse_sigma2;

      likelihood_linear +=
        arma::dot(
          workspace.residual_without_foreign,
          workspace.foreign_signal
        ) *
        inverse_sigma2;
    }

    if (
      !std::isfinite(
        likelihood_quadratic
      ) ||
      likelihood_quadratic < 0.0 ||
      !std::isfinite(
        likelihood_linear
      )
    ) {
      Rcpp::stop(
        "Invalid non-self ASIS likelihood coefficients"
      );
    }
  }


  inline AsisResult update_nonself_asis(
      arma::mat& beta,
      arma::vec& tau2,
      arma::vec& xi,
      const arma::mat& Y,
      const arma::mat& X,
      const arma::vec& sigma2,
      const arma::uvec& predictor_series,
      const int grouping_code,
      const double tau_df,
      const double tau_scale,
      const double slice_width,
      const int maximum_step_out,
      const int maximum_shrinks,
      AsisWorkspace& workspace) {

    if (
      grouping_code !=
      grouping_self_diagonal
    ) {
      Rcpp::stop(
        "Non-self ASIS requires self-diagonal grouping"
      );
    }

    if (
      tau2.n_elem != 2u ||
      xi.n_elem != 2u
    ) {
      Rcpp::stop(
        "Non-self ASIS requires exactly two global groups"
      );
    }

    if (
      !std::isfinite(
        tau_df
      ) ||
      tau_df <= 0.0 ||
      !std::isfinite(
        tau_scale
      ) ||
      tau_scale <= 0.0
    ) {
      Rcpp::stop(
        "Invalid non-self ASIS tau parameters"
      );
    }

    const arma::uword foreign_group =
      1u;

    AsisResult result;

    result.tau2_before =
      tau2[
        foreign_group
      ];

    result.xi_before =
      xi[
        foreign_group
      ];

    if (
      !std::isfinite(
        result.tau2_before
      ) ||
      result.tau2_before <= 0.0 ||
      !std::isfinite(
        result.xi_before
      ) ||
      result.xi_before <= 0.0
    ) {
      Rcpp::stop(
        "Invalid non-self global state before ASIS"
      );
    }

    const double tau_before =
      std::sqrt(
        result.tau2_before
      );

    compute_nonself_asis_likelihood(
      result.likelihood_quadratic,
      result.likelihood_linear,
      beta,
      Y,
      X,
      sigma2,
      predictor_series,
      tau_before,
      workspace
    );

    const double eta_before =
      std::log(
        tau_before
      );

    const AsisSliceResult slice =
      draw_nonself_asis_eta(
        eta_before,
        result.likelihood_quadratic,
        result.likelihood_linear,
        result.xi_before,
        tau_df,
        slice_width,
        maximum_step_out,
        maximum_shrinks
      );

    const double tau_after =
      std::exp(
        slice.eta
      );

    result.tau2_after =
      tau_after *
      tau_after;

    result.log_tau_move =
      slice.eta -
      eta_before;

    result.evaluations =
      slice.evaluations;

    result.left_expansions =
      slice.left_expansions;

    result.right_expansions =
      slice.right_expansions;

    result.shrink_steps =
      slice.shrink_steps;

    if (
      !std::isfinite(
        result.tau2_after
      ) ||
      result.tau2_after <= 0.0 ||
      !std::isfinite(
        result.log_tau_move
      )
    ) {
      Rcpp::stop(
        "Invalid non-self global scale after ASIS"
      );
    }

    const double beta_rescaling =
      tau_after /
      tau_before;

    const arma::uword k =
      beta.n_rows;

    const arma::uword N =
      beta.n_cols;

    double maximum_relative_error =
      0.0;

    for (
      arma::uword equation = 0u;
      equation < N;
      ++equation
    ) {

      for (
        arma::uword r = 0u;
        r < k;
        ++r
      ) {

        if (
          predictor_series[r] !=
          equation
        ) {

          const double beta_before =
            beta(
              r,
              equation
            );

          const double alpha_before =
            beta_before /
            tau_before;

          const double beta_after =
            beta_before *
            beta_rescaling;

          if (
            !std::isfinite(
              beta_after
            )
          ) {
            Rcpp::stop(
              "Non-finite coefficient produced by non-self ASIS"
            );
          }

          beta(
            r,
            equation
          ) =
            beta_after;

          const double alpha_after =
            beta_after /
            tau_after;

          const double denominator =
            std::max(
              1.0,
              std::abs(
                alpha_before
              )
            );

          const double relative_error =
            std::abs(
              alpha_after -
              alpha_before
            ) /
            denominator;

          maximum_relative_error =
            std::max(
              maximum_relative_error,
              relative_error
            );
        }
      }
    }

    result.maximum_relative_invariant_error =
      maximum_relative_error;

    tau2[
      foreign_group
    ] =
      result.tau2_after;

    result.xi_after =
      draw_half_t_auxiliary(
        result.tau2_after,
        tau_df,
        tau_scale
      );

    xi[
      foreign_group
    ] =
      result.xi_after;

    return result;
  }


}  // namespace half_t_bvar


#endif
