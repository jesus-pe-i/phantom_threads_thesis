#ifndef GIGG_BVAR_UPDATES_H
#define GIGG_BVAR_UPDATES_H

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>


// GIGG-BVAR conditional-update layer.
// Gaussian coefficient updates, exact GIG block updates,
// coefficient-local scales, and the two GIGG ASIS moves.


namespace gigg_bvar {


// ============================================================
// 1. Codes, constants and validation helpers
// ============================================================

constexpr int beta_auto = 0;
constexpr int beta_chol = 1;
constexpr int beta_bhattacharya = 2;

constexpr double gig_chi_floor = 1e-12;
constexpr double asis_invariant_tolerance = 1e-10;
constexpr int gig_maximum_attempts = 100000;


inline void require_positive_finite(
    const double value,
    const std::string& name) {

  if (!std::isfinite(value) || value <= 0.0) {
    Rcpp::stop(name + " must be finite and positive");
  }
}


inline int parse_beta_algorithm(const std::string& algorithm) {

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


inline std::string beta_algorithm_name(const int algorithm) {

  if (algorithm == beta_auto) {
    return "auto";
  }

  if (algorithm == beta_chol) {
    return "chol";
  }

  if (algorithm == beta_bhattacharya) {
    return "bhattacharya";
  }

  Rcpp::stop("Unknown beta algorithm code");
  return "unknown";
}


// Frozen empirical crossover:
//
//   k <= 2 T_p  -> Cholesky
//   k >  2 T_p  -> Bhattacharya

inline int resolve_beta_algorithm(
    const int requested_algorithm,
    const arma::uword k,
    const arma::uword T_p) {

  if (
    requested_algorithm == beta_chol ||
    requested_algorithm == beta_bhattacharya
  ) {
    return requested_algorithm;
  }

  if (requested_algorithm != beta_auto) {
    Rcpp::stop("Unknown beta algorithm code");
  }

  return k <= 2ULL * T_p ? beta_chol : beta_bhattacharya;
}


// ============================================================
// 2. Scalar RNG helpers
// ============================================================

// IG(shape, scale) uses density
//
//   f(x) proportional to x^(-shape - 1) exp(-scale / x).

inline double rinvgamma(
    const double shape,
    const double scale) {

  require_positive_finite(shape, "inverse-gamma shape");
  require_positive_finite(scale, "inverse-gamma scale");

  const double gamma_draw = R::rgamma(shape, 1.0 / scale);

  if (!std::isfinite(gamma_draw) || gamma_draw <= 0.0) {
    Rcpp::stop("Invalid inverse-gamma denominator draw");
  }

  const double result = 1.0 / gamma_draw;

  if (!std::isfinite(result) || result <= 0.0) {
    Rcpp::stop("Invalid inverse-gamma result");
  }

  return result;
}


// Gamma(shape, rate) helper. R::rgamma uses shape-scale.

inline double rgamma_rate(
    const double shape,
    const double rate) {

  require_positive_finite(shape, "gamma shape");
  require_positive_finite(rate, "gamma rate");

  const double result = R::rgamma(shape, 1.0 / rate);

  if (!std::isfinite(result) || result <= 0.0) {
    Rcpp::stop("Invalid gamma result");
  }

  return result;
}


// Non-throwing scalar draws used by ASIS. A failed ASIS draw is
// handled as an identity move, which preserves the target.

inline bool try_rinvgamma(
    const double shape,
    const double scale,
    double& result) {

  result = NA_REAL;

  if (
    !std::isfinite(shape) || shape <= 0.0 ||
    !std::isfinite(scale) || scale <= 0.0
  ) {
    return false;
  }

  const double gamma_draw = R::rgamma(shape, 1.0 / scale);

  if (!std::isfinite(gamma_draw) || gamma_draw <= 0.0) {
    return false;
  }

  result = 1.0 / gamma_draw;

  return std::isfinite(result) && result > 0.0;
}


inline bool try_rgamma_rate(
    const double shape,
    const double rate,
    double& result) {

  result = NA_REAL;

  if (
    !std::isfinite(shape) || shape <= 0.0 ||
    !std::isfinite(rate) || rate <= 0.0
  ) {
    return false;
  }

  result = R::rgamma(shape, 1.0 / rate);

  return std::isfinite(result) && result > 0.0;
}


// ============================================================
// 3. Native GIG acceptance-rejection sampler
// ============================================================

// Parameter convention:
//
//   X ~ GIG(p, chi, psi)
//
// has density
//
//   f(x) proportional to
//     x^(p - 1) exp{-(psi x + chi / x) / 2}, x > 0.
//
// The implementation is an exact rejection sampler on Y = log(X).
// Its log density is strictly concave:
//
//   ell(y) = p y - (psi exp(y) + chi exp(-y)) / 2 + constant.
//
// A three-piece envelope is used: the tangent at a left point,
// the horizontal tangent at the mode, and the tangent at a right
// point. The left and right points are where the log density has
// fallen by one unit. Concavity proves that all three pieces are
// valid upper bounds. Their integrals and inverse CDFs are analytic.
//
// Centering at the log-density mode gives the stable identity
//
//   ell(y_mode + d) - ell(y_mode)
//     = -p {sinh(d) - d}
//       - kappa {cosh(d) - 1},
//
//   kappa = sqrt(p^2 + chi psi).
//
// This supports negative p directly and remains well behaved for
// small chi and concentrated laws. No normalising constant, Bessel
// function, R callback, or reciprocal approximation is required.

struct GIGDrawResult {

  double value;
  int attempts;

  GIGDrawResult() :
    value(NA_REAL),
    attempts(0) {}
};


inline double safe_exp(const double log_value) {

  static const double log_max = std::log(
    std::numeric_limits<double>::max()
  );

  static const double log_min = std::log(
    std::numeric_limits<double>::denorm_min()
  );

  if (log_value > log_max) {
    return std::numeric_limits<double>::infinity();
  }

  if (log_value < log_min) {
    return 0.0;
  }

  return std::exp(log_value);
}


// Positive root of
//
//   psi x^2 - 2 c x - chi = 0,
//
// returned as log(x). When c < 0 the conjugate expression avoids
// subtractive cancellation.

inline double log_positive_quadratic_root(
    const double c,
    const double chi,
    const double psi) {

  const double root_product =
    std::sqrt(chi) * std::sqrt(psi);

  const double discriminant = std::hypot(c, root_product);

  if (!std::isfinite(discriminant) || discriminant <= 0.0) {
    Rcpp::stop("Invalid GIG quadratic-root discriminant");
  }

  double log_root;

  if (c >= 0.0) {

    const double numerator = c + discriminant;

    if (!std::isfinite(numerator) || numerator <= 0.0) {
      Rcpp::stop("Invalid direct GIG quadratic root");
    }

    log_root = std::log(numerator) - std::log(psi);

  } else {

    const double denominator = discriminant - c;

    if (!std::isfinite(denominator) || denominator <= 0.0) {
      Rcpp::stop("Invalid conjugate GIG quadratic root");
    }

    log_root = std::log(chi) - std::log(denominator);
  }

  if (!std::isfinite(log_root)) {
    Rcpp::stop("Non-finite log GIG quadratic root");
  }

  return log_root;
}


inline double cosh_minus_one(const double value) {

  if (!std::isfinite(value) || std::abs(value) > 700.0) {
    return std::numeric_limits<double>::infinity();
  }

  const double half_sinh = std::sinh(0.5 * value);

  return 2.0 * half_sinh * half_sinh;
}


inline double sinh_minus_identity(const double value) {

  if (!std::isfinite(value) || std::abs(value) > 700.0) {
    return std::copysign(
      std::numeric_limits<double>::infinity(),
      value
    );
  }

  const double absolute_value = std::abs(value);

  if (absolute_value < 1e-3) {

    const double value2 = value * value;

    return value * value2 *
      (
        1.0 / 6.0 +
        value2 *
        (
          1.0 / 120.0 +
          value2 / 5040.0
        )
      );
  }

  return std::sinh(value) - value;
}


inline double gig_shifted_log_kernel(
    const double displacement,
    const double p,
    const double kappa) {

  if (!std::isfinite(displacement)) {
    return -std::numeric_limits<double>::infinity();
  }

  const double sinh_remainder =
    sinh_minus_identity(displacement);

  const double cosh_remainder =
    cosh_minus_one(displacement);

  const double result =
    -p * sinh_remainder -
    kappa * cosh_remainder;

  return std::isfinite(result)
    ? std::min(0.0, result)
    : -std::numeric_limits<double>::infinity();
}


inline double gig_shifted_log_derivative(
    const double displacement,
    const double p,
    const double kappa) {

  const double cosh_remainder =
    cosh_minus_one(displacement);

  if (
    !std::isfinite(cosh_remainder) ||
    !std::isfinite(displacement) ||
    std::abs(displacement) > 700.0
  ) {
    return displacement < 0.0
      ? std::numeric_limits<double>::infinity()
      : -std::numeric_limits<double>::infinity();
  }

  return
    -p * cosh_remainder -
    kappa * std::sinh(displacement);
}


inline double find_gig_drop_point(
    const int direction,
    const double p,
    const double kappa,
    const double log_drop = 1.0) {

  if (direction != -1 && direction != 1) {
    Rcpp::stop("GIG drop-point direction must be -1 or 1");
  }

  require_positive_finite(kappa, "GIG kappa");
  require_positive_finite(log_drop, "GIG envelope log drop");

  double step = std::sqrt(2.0 * log_drop / kappa);

  if (!std::isfinite(step) || step <= 0.0) {
    Rcpp::stop("Invalid initial GIG envelope step");
  }

  step = std::min(1.0, step);

  double outside =
    static_cast<double>(direction) * step;

  bool bracketed = false;

  for (int expansion = 0; expansion < 1024; ++expansion) {

    if (
      gig_shifted_log_kernel(outside, p, kappa) <=
      -log_drop
    ) {
      bracketed = true;
      break;
    }

    outside *= 2.0;

    if (!std::isfinite(outside)) {
      break;
    }
  }

  if (!bracketed) {
    Rcpp::stop("Could not bracket a GIG envelope drop point");
  }

  double inside = 0.0;

  for (int iteration = 0; iteration < 80; ++iteration) {

    const double midpoint =
      0.5 * (outside + inside);

    if (
      gig_shifted_log_kernel(midpoint, p, kappa) <=
      -log_drop
    ) {
      outside = midpoint;
    } else {
      inside = midpoint;
    }
  }

  const double result =
    0.5 * (outside + inside);

  if (!std::isfinite(result) || result == 0.0) {
    Rcpp::stop("Invalid GIG envelope drop point");
  }

  return result;
}


inline GIGDrawResult draw_gig(
    const double p,
    const double chi,
    const double psi,
    const int maximum_attempts = 100000) {

  if (
    !std::isfinite(p) ||
    !std::isfinite(chi) || chi <= 0.0 ||
    !std::isfinite(psi) || psi <= 0.0 ||
    maximum_attempts < 1
  ) {
    Rcpp::stop("Invalid GIG parameters");
  }

  GIGDrawResult result;

  const double root_product =
    std::sqrt(chi) * std::sqrt(psi);

  const double kappa =
    std::hypot(p, root_product);

  require_positive_finite(kappa, "GIG kappa");

  // Mode of the log-X density, not the mode of X.
  const double log_mode =
    log_positive_quadratic_root(p, chi, psi);

  const double left =
    find_gig_drop_point(-1, p, kappa);

  const double right =
    find_gig_drop_point(1, p, kappa);

  const double left_log_kernel =
    gig_shifted_log_kernel(left, p, kappa);

  const double right_log_kernel =
    gig_shifted_log_kernel(right, p, kappa);

  const double left_slope =
    gig_shifted_log_derivative(left, p, kappa);

  const double right_slope =
    gig_shifted_log_derivative(right, p, kappa);

  if (
    !std::isfinite(left_slope) || left_slope <= 0.0 ||
    !std::isfinite(right_slope) || right_slope >= 0.0
  ) {
    Rcpp::stop("Invalid GIG envelope slopes");
  }

  const double left_weight =
    std::exp(left_log_kernel) / left_slope;

  const double center_weight =
    right - left;

  const double right_weight =
    std::exp(right_log_kernel) / (-right_slope);

  const double total_weight =
    left_weight + center_weight + right_weight;

  if (
    !std::isfinite(left_weight) || left_weight <= 0.0 ||
    !std::isfinite(center_weight) || center_weight <= 0.0 ||
    !std::isfinite(right_weight) || right_weight <= 0.0 ||
    !std::isfinite(total_weight) || total_weight <= 0.0
  ) {
    Rcpp::stop("Invalid GIG envelope weights");
  }

  double accepted_displacement = NA_REAL;

  for (
    int attempt = 1;
    attempt <= maximum_attempts;
    ++attempt
  ) {

    ++result.attempts;

    const double component_draw =
      R::runif(0.0, total_weight);

    double position_draw =
      R::runif(0.0, 1.0);

    double acceptance_draw =
      R::runif(0.0, 1.0);

    if (
      !std::isfinite(component_draw) ||
      !std::isfinite(position_draw) ||
      !std::isfinite(acceptance_draw)
    ) {
      continue;
    }

    position_draw = std::max(
      position_draw,
      std::numeric_limits<double>::min()
    );

    acceptance_draw = std::max(
      acceptance_draw,
      std::numeric_limits<double>::min()
    );

    double displacement;
    double envelope_log_kernel;

    if (component_draw < left_weight) {

      displacement =
        left +
        std::log(position_draw) / left_slope;

      envelope_log_kernel =
        left_log_kernel +
        left_slope * (displacement - left);

    } else if (
      component_draw <
      left_weight + center_weight
    ) {

      displacement =
        left +
        (right - left) * position_draw;

      envelope_log_kernel = 0.0;

    } else {

      displacement =
        right -
        std::log(position_draw) / (-right_slope);

      envelope_log_kernel =
        right_log_kernel +
        right_slope * (displacement - right);
    }

    const double target_log_kernel =
      gig_shifted_log_kernel(
        displacement,
        p,
        kappa
      );

    if (!std::isfinite(target_log_kernel)) {
      continue;
    }

    const double log_acceptance = std::min(
      0.0,
      target_log_kernel - envelope_log_kernel
    );

    if (
      std::log(acceptance_draw) <=
      log_acceptance
    ) {
      accepted_displacement = displacement;
      break;
    }
  }

  if (!std::isfinite(accepted_displacement)) {
    Rcpp::stop(
      "GIG rejection sampler exceeded maximum_attempts"
    );
  }

  result.value = safe_exp(
    log_mode + accepted_displacement
  );

  if (!std::isfinite(result.value) || result.value <= 0.0) {
    Rcpp::stop(
      "GIG sampler produced a non-finite or non-positive draw"
    );
  }

  return result;
}


// ============================================================
// 4. Gaussian coefficient updates
// ============================================================

struct CholWorkspace {

  arma::mat precision;
  arma::mat U;

  arma::vec rhs;
  arma::vec intermediate;
  arma::vec posterior_mean;
  arma::vec z;
  arma::vec perturbation;

  explicit CholWorkspace(const arma::uword k) :
    precision(k, k, arma::fill::zeros),
    U(k, k, arma::fill::zeros),
    rhs(k, arma::fill::zeros),
    intermediate(k, arma::fill::zeros),
    posterior_mean(k, arma::fill::zeros),
    z(k, arma::fill::zeros),
    perturbation(k, arma::fill::zeros) {}
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
    u(k, arma::fill::zeros),
    delta(T_p, arma::fill::zeros),
    v(T_p, arma::fill::zeros),
    rhs(T_p, arma::fill::zeros),
    intermediate(T_p, arma::fill::zeros),
    w(T_p, arma::fill::zeros),
    Xt_w(k, arma::fill::zeros),
    XD(T_p, k, arma::fill::zeros),
    system(T_p, T_p, arma::fill::zeros),
    U(T_p, T_p, arma::fill::zeros) {}
};


inline void draw_beta_chol(
    arma::vec& beta,
    const arma::mat& XtX,
    const arma::mat& XtY,
    const arma::uword equation,
    const arma::vec& omega,
    const double sigma2,
    CholWorkspace& workspace) {

  const arma::uword k = XtX.n_rows;

  workspace.precision = XtX;

  for (arma::uword row = 0u; row < k; ++row) {
    workspace.precision(row, row) += 1.0 / omega[row];
  }

  if (!arma::chol(workspace.U, workspace.precision)) {
    Rcpp::stop(
      "Beta posterior precision is not positive definite"
    );
  }

  workspace.rhs = XtY.col(equation);

  bool solve_ok = arma::solve(
    workspace.intermediate,
    arma::trimatl(workspace.U.t()),
    workspace.rhs,
    arma::solve_opts::fast
  );

  if (!solve_ok) {
    Rcpp::stop("Lower-triangular beta-mean solve failed");
  }

  solve_ok = arma::solve(
    workspace.posterior_mean,
    arma::trimatu(workspace.U),
    workspace.intermediate,
    arma::solve_opts::fast
  );

  if (!solve_ok) {
    Rcpp::stop("Upper-triangular beta-mean solve failed");
  }

  for (arma::uword row = 0u; row < k; ++row) {
    workspace.z[row] = R::rnorm(0.0, 1.0);
  }

  solve_ok = arma::solve(
    workspace.perturbation,
    arma::trimatu(workspace.U),
    workspace.z,
    arma::solve_opts::fast
  );

  if (!solve_ok) {
    Rcpp::stop("Upper-triangular beta-noise solve failed");
  }

  beta =
    workspace.posterior_mean +
    std::sqrt(sigma2) * workspace.perturbation;

  if (!beta.is_finite()) {
    Rcpp::stop(
      "Cholesky beta update produced non-finite values"
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

  const arma::uword T_p = X.n_rows;
  const arma::uword k = X.n_cols;
  const double sigma = std::sqrt(sigma2);

  for (arma::uword row = 0u; row < k; ++row) {
    workspace.u[row] =
      sigma *
      std::sqrt(omega[row]) *
      R::rnorm(0.0, 1.0);
  }

  for (arma::uword time = 0u; time < T_p; ++time) {
    workspace.delta[time] =
      sigma * R::rnorm(0.0, 1.0);
  }

  workspace.v =
    X * workspace.u + workspace.delta;

  workspace.XD = X;
  workspace.XD.each_row() %= omega.t();

  workspace.system =
    workspace.XD * X.t();

  workspace.system.diag() += 1.0;

  if (!arma::chol(workspace.U, workspace.system)) {
    Rcpp::stop(
      "Bhattacharya system is not positive definite"
    );
  }

  workspace.rhs =
    Y.col(equation) - workspace.v;

  bool solve_ok = arma::solve(
    workspace.intermediate,
    arma::trimatl(workspace.U.t()),
    workspace.rhs,
    arma::solve_opts::fast
  );

  if (!solve_ok) {
    Rcpp::stop(
      "Lower-triangular Bhattacharya solve failed"
    );
  }

  solve_ok = arma::solve(
    workspace.w,
    arma::trimatu(workspace.U),
    workspace.intermediate,
    arma::solve_opts::fast
  );

  if (!solve_ok) {
    Rcpp::stop(
      "Upper-triangular Bhattacharya solve failed"
    );
  }

  workspace.Xt_w =
    X.t() * workspace.w;

  beta =
    workspace.u +
    omega % workspace.Xt_w;

  if (!beta.is_finite()) {
    Rcpp::stop(
      "Bhattacharya beta update produced non-finite values"
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

  if (resolved_algorithm == beta_chol) {

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

  if (resolved_algorithm == beta_bhattacharya) {

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

  Rcpp::stop("Unknown resolved beta algorithm");
}


// Checked convenience wrappers for low-level conditional tests.

inline arma::vec draw_beta_chol(
    const arma::vec& y,
    const arma::mat& X,
    const arma::vec& omega,
    const double sigma2) {

  if (
    y.n_elem != X.n_rows ||
    omega.n_elem != X.n_cols ||
    !y.is_finite() ||
    !X.is_finite() ||
    !omega.is_finite() ||
    arma::any(omega <= 0.0) ||
    !std::isfinite(sigma2) ||
    sigma2 <= 0.0
  ) {
    Rcpp::stop(
      "Invalid checked Cholesky beta inputs"
    );
  }

  const arma::mat Y =
    arma::reshape(y, X.n_rows, 1u);

  const arma::mat XtX =
    X.t() * X;

  const arma::mat XtY =
    X.t() * Y;

  arma::vec beta(X.n_cols);
  CholWorkspace workspace(X.n_cols);

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

  if (
    y.n_elem != X.n_rows ||
    omega.n_elem != X.n_cols ||
    !y.is_finite() ||
    !X.is_finite() ||
    !omega.is_finite() ||
    arma::any(omega <= 0.0) ||
    !std::isfinite(sigma2) ||
    sigma2 <= 0.0
  ) {
    Rcpp::stop(
      "Invalid checked Bhattacharya beta inputs"
    );
  }

  const arma::mat Y =
    arma::reshape(y, X.n_rows, 1u);

  arma::vec beta(X.n_cols);

  BhattacharyaWorkspace workspace(
    X.n_rows,
    X.n_cols
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


// ============================================================
// 5. Exact lag-receiver-sender block cache
// ============================================================

struct BlockMembershipCache {

  arma::uword n_rows;
  arma::uword n_equations;

  std::vector<arma::uvec> linear_indices;
  arma::uvec coefficient_count;
};


inline arma::uvec uvec_from_std(
    const std::vector<arma::uword>& values) {

  arma::uvec output(values.size());

  for (
    arma::uword index = 0u;
    index < output.n_elem;
    ++index
  ) {
    output[index] =
      values[static_cast<std::size_t>(index)];
  }

  return output;
}


inline BlockMembershipCache build_block_membership_cache(
    const arma::imat& block_id,
    const int n_blocks,
    const arma::uword expected_block_size = 0u) {

  if (n_blocks < 1) {
    Rcpp::stop("n_blocks must be positive");
  }

  const arma::uword k = block_id.n_rows;
  const arma::uword N = block_id.n_cols;

  if (k < 1u || N < 1u) {
    Rcpp::stop(
      "block_id must have positive dimensions"
    );
  }

  std::vector<std::vector<arma::uword> > membership(
    static_cast<std::size_t>(n_blocks)
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

      const int block =
        block_id(row, equation);

      if (block < 1 || block > n_blocks) {
        Rcpp::stop(
          "block_id contains a value outside 1:n_blocks"
        );
      }

      membership[
        static_cast<std::size_t>(block - 1)
      ].push_back(
        row + equation * k
      );
    }
  }

  BlockMembershipCache cache;

  cache.n_rows = k;
  cache.n_equations = N;

  cache.linear_indices.resize(
    static_cast<std::size_t>(n_blocks)
  );

  cache.coefficient_count.set_size(
    static_cast<arma::uword>(n_blocks)
  );

  for (int block = 0; block < n_blocks; ++block) {

    const std::vector<arma::uword>& block_membership =
      membership[static_cast<std::size_t>(block)];

    if (block_membership.empty()) {
      Rcpp::stop(
        "block_id contains an empty block"
      );
    }

    if (
      expected_block_size > 0u &&
      block_membership.size() !=
        static_cast<std::size_t>(expected_block_size)
    ) {
      Rcpp::stop(
        "A GIGG block does not contain the expected m^2 coefficients"
      );
    }

    cache.linear_indices[
      static_cast<std::size_t>(block)
    ] = uvec_from_std(block_membership);

    cache.coefficient_count[
      static_cast<arma::uword>(block)
    ] = static_cast<arma::uword>(
      block_membership.size()
    );
  }

  return cache;
}


inline arma::imat zero_based_block_map(
    const arma::imat& block_id,
    const int n_blocks) {

  arma::imat output(
    block_id.n_rows,
    block_id.n_cols
  );

  for (
    arma::uword equation = 0u;
    equation < block_id.n_cols;
    ++equation
  ) {
    for (
      arma::uword row = 0u;
      row < block_id.n_rows;
      ++row
    ) {

      const int block =
        block_id(row, equation);

      if (block < 1 || block > n_blocks) {
        Rcpp::stop(
          "block_id contains a value outside 1:n_blocks"
        );
      }

      output(row, equation) = block - 1;
    }
  }

  return output;
}


// ============================================================
// 6. Prior variance multipliers and sigma2 update
// ============================================================

inline void fill_omega(
    arma::vec& omega,
    const arma::uword equation,
    const double tau2,
    const arma::vec& gamma2,
    const arma::mat& lambda2,
    const arma::imat& block_id_zero) {

  const arma::uword k = lambda2.n_rows;

  for (arma::uword row = 0u; row < k; ++row) {

    const arma::uword block =
      static_cast<arma::uword>(
        block_id_zero(row, equation)
      );

    if (block >= gamma2.n_elem) {
      Rcpp::stop(
        "Invalid block index while constructing omega"
      );
    }

    const double value =
      tau2 *
      gamma2[block] *
      lambda2(row, equation);

    if (!std::isfinite(value) || value <= 0.0) {
      Rcpp::stop(
        "Invalid GIGG coefficient prior variance multiplier"
      );
    }

    omega[row] = value;
  }
}


inline double draw_sigma2_from_residual(
    const arma::vec& residual,
    const arma::vec& beta,
    const arma::vec& omega,
    const double prior_shape,
    const double prior_scale) {

  if (
    residual.n_elem < 1u ||
    beta.n_elem != omega.n_elem ||
    !residual.is_finite() ||
    !beta.is_finite() ||
    !omega.is_finite() ||
    arma::any(omega <= 0.0) ||
    !std::isfinite(prior_shape) ||
    prior_shape <= 0.0 ||
    !std::isfinite(prior_scale) ||
    prior_scale <= 0.0
  ) {
    Rcpp::stop(
      "Invalid sigma2 conditional inputs"
    );
  }

  const double rss =
    arma::dot(residual, residual);

  const double beta_quadratic =
    arma::dot(
      beta % beta,
      1.0 / omega
    );

  if (
    !std::isfinite(rss) || rss < 0.0 ||
    !std::isfinite(beta_quadratic) ||
    beta_quadratic < 0.0
  ) {
    Rcpp::stop(
      "Invalid sigma2 conditional quadratic"
    );
  }

  const double posterior_shape =
    prior_shape +
    0.5 * static_cast<double>(
      residual.n_elem + beta.n_elem
    );

  const double posterior_scale =
    prior_scale +
    0.5 * (rss + beta_quadratic);

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
      "Invalid checked sigma2 dimensions"
    );
  }

  const arma::vec residual =
    y - X * beta;

  return draw_sigma2_from_residual(
    residual,
    beta,
    omega,
    prior_shape,
    prior_scale
  );
}


// ============================================================
// 7. Global half-t hierarchy
// ============================================================

inline double global_tau_quadratic(
    const arma::mat& beta,
    const arma::vec& sigma2,
    const arma::vec& gamma2,
    const arma::mat& lambda2,
    const arma::imat& block_id_zero) {

  const arma::uword k = beta.n_rows;
  const arma::uword N = beta.n_cols;

  if (
    sigma2.n_elem != N ||
    lambda2.n_rows != k ||
    lambda2.n_cols != N ||
    block_id_zero.n_rows != k ||
    block_id_zero.n_cols != N ||
    gamma2.n_elem < 1u ||
    !beta.is_finite() ||
    !sigma2.is_finite() ||
    arma::any(sigma2 <= 0.0) ||
    !gamma2.is_finite() ||
    arma::any(gamma2 <= 0.0) ||
    !lambda2.is_finite() ||
    arma::any(arma::vectorise(lambda2) <= 0.0)
  ) {
    Rcpp::stop(
      "Invalid global tau conditional inputs"
    );
  }

  double quadratic = 0.0;

  for (
    arma::uword equation = 0u;
    equation < N;
    ++equation
  ) {

    const double inverse_sigma2 =
      1.0 / sigma2[equation];

    for (
      arma::uword row = 0u;
      row < k;
      ++row
    ) {

      const arma::uword block =
        static_cast<arma::uword>(
          block_id_zero(row, equation)
        );

      if (block >= gamma2.n_elem) {
        Rcpp::stop(
          "Invalid block index in global tau conditional"
        );
      }

      const double coefficient =
        beta(row, equation);

      quadratic +=
        coefficient *
        coefficient *
        inverse_sigma2 /
        (
          gamma2[block] *
          lambda2(row, equation)
        );
    }
  }

  if (!std::isfinite(quadratic) || quadratic < 0.0) {
    Rcpp::stop("Invalid global tau quadratic");
  }

  return quadratic;
}


inline double draw_tau2(
    const arma::mat& beta,
    const arma::vec& sigma2,
    const arma::vec& gamma2,
    const arma::mat& lambda2,
    const arma::imat& block_id_zero,
    const double tau_df,
    const double psi_tau) {

  require_positive_finite(tau_df, "tau_df");
  require_positive_finite(psi_tau, "psi_tau");

  const double quadratic =
    global_tau_quadratic(
      beta,
      sigma2,
      gamma2,
      lambda2,
      block_id_zero
    );

  const double posterior_shape =
    0.5 *
    (
      tau_df +
      static_cast<double>(beta.n_elem)
    );

  const double posterior_scale =
    tau_df / psi_tau +
    0.5 * quadratic;

  return rinvgamma(
    posterior_shape,
    posterior_scale
  );
}


// For tau ~ half-t_nu(0, A):
//
//   tau^2 | psi_tau ~ IG(nu/2, nu/psi_tau)
//   psi_tau         ~ IG(1/2, 1/A^2).

inline double draw_half_t_auxiliary(
    const double tau2,
    const double tau_df,
    const double tau_scale) {

  require_positive_finite(tau2, "tau2");
  require_positive_finite(tau_df, "tau_df");
  require_positive_finite(tau_scale, "tau_scale");

  const double posterior_shape =
    0.5 * (tau_df + 1.0);

  const double posterior_scale =
    tau_df / tau2 +
    1.0 / (tau_scale * tau_scale);

  return rinvgamma(
    posterior_shape,
    posterior_scale
  );
}


// ============================================================
// 8. Ordinary GIG block-scale update
// ============================================================

inline double block_gamma_quadratic(
    const arma::mat& beta,
    const arma::vec& sigma2,
    const arma::mat& lambda2,
    const BlockMembershipCache& cache,
    const arma::uword block) {

  const arma::uvec& indices =
    cache.linear_indices[
      static_cast<std::size_t>(block)
    ];

  const arma::uword k = cache.n_rows;
  double quadratic = 0.0;

  for (
    arma::uword position = 0u;
    position < indices.n_elem;
    ++position
  ) {

    const arma::uword linear =
      indices[position];

    const arma::uword equation =
      linear / k;

    const double coefficient =
      beta[linear];

    quadratic +=
      coefficient *
      coefficient /
      (
        sigma2[equation] *
        lambda2[linear]
      );
  }

  if (!std::isfinite(quadratic) || quadratic < 0.0) {
    Rcpp::stop(
      "Invalid block-scale quadratic"
    );
  }

  return quadratic;
}


struct GIGUpdateSummary {

  std::uint64_t draws;
  std::uint64_t total_attempts;
  std::uint64_t chi_floor_count;

  int maximum_attempts_for_one_draw;

  double minimum_raw_chi;
  double maximum_raw_chi;

  GIGUpdateSummary() :
    draws(0u),
    total_attempts(0u),
    chi_floor_count(0u),
    maximum_attempts_for_one_draw(0),
    minimum_raw_chi(
      std::numeric_limits<double>::infinity()
    ),
    maximum_raw_chi(0.0) {}
};


inline GIGUpdateSummary update_block_scales_gig(
    arma::vec& gamma2,
    const arma::mat& beta,
    const arma::vec& sigma2,
    const arma::mat& lambda2,
    const BlockMembershipCache& cache,
    const double tau2,
    const double gamma_shape,
    const double gamma_rate,
    const double numerical_epsilon = gig_chi_floor,
    const int maximum_attempts = gig_maximum_attempts) {

  if (
    gamma2.n_elem !=
      cache.coefficient_count.n_elem ||
    beta.n_rows != cache.n_rows ||
    beta.n_cols != cache.n_equations ||
    lambda2.n_rows != cache.n_rows ||
    lambda2.n_cols != cache.n_equations ||
    sigma2.n_elem != cache.n_equations ||
    !beta.is_finite() ||
    !sigma2.is_finite() ||
    arma::any(sigma2 <= 0.0) ||
    !lambda2.is_finite() ||
    arma::any(arma::vectorise(lambda2) <= 0.0)
  ) {
    Rcpp::stop(
      "Inconsistent dimensions in the block GIG update"
    );
  }

  require_positive_finite(tau2, "tau2");

  require_positive_finite(
    gamma_shape,
    "gamma_shape"
  );

  require_positive_finite(
    gamma_rate,
    "gamma_rate"
  );

  require_positive_finite(
    numerical_epsilon,
    "numerical_epsilon"
  );

  GIGUpdateSummary summary;

  const double gig_psi =
    2.0 * gamma_rate;

  for (
    arma::uword block = 0u;
    block < gamma2.n_elem;
    ++block
  ) {

    const double S_b =
      block_gamma_quadratic(
        beta,
        sigma2,
        lambda2,
        cache,
        block
      );

    const double raw_chi =
      S_b / tau2;

    if (
      !std::isfinite(raw_chi) ||
      raw_chi < 0.0
    ) {
      Rcpp::stop(
        "Invalid raw chi in the block GIG update"
      );
    }

    summary.minimum_raw_chi = std::min(
      summary.minimum_raw_chi,
      raw_chi
    );

    summary.maximum_raw_chi = std::max(
      summary.maximum_raw_chi,
      raw_chi
    );

    double gig_chi = raw_chi;

    if (gig_chi < numerical_epsilon) {
      gig_chi = numerical_epsilon;
      ++summary.chi_floor_count;
    }

    const double gig_p =
      gamma_shape -
      0.5 *
      static_cast<double>(
        cache.coefficient_count[block]
      );

    const GIGDrawResult draw =
      draw_gig(
        gig_p,
        gig_chi,
        gig_psi,
        maximum_attempts
      );

    gamma2[block] = draw.value;

    ++summary.draws;

    summary.total_attempts +=
      static_cast<std::uint64_t>(
        draw.attempts
      );

    summary.maximum_attempts_for_one_draw =
      std::max(
        summary.maximum_attempts_for_one_draw,
        draw.attempts
      );
  }

  if (summary.draws == 0u) {
    summary.minimum_raw_chi = NA_REAL;
    summary.maximum_raw_chi = NA_REAL;
  }

  return summary;
}


// ============================================================
// 9. Ordinary coefficient-local inverse-gamma update
// ============================================================

inline void update_local_scales(
    arma::mat& lambda2,
    const arma::mat& beta,
    const arma::vec& sigma2,
    const arma::vec& gamma2,
    const arma::imat& block_id_zero,
    const double tau2,
    const double lambda_shape,
    const double lambda_prior_scale) {

  if (
    lambda2.n_rows != beta.n_rows ||
    lambda2.n_cols != beta.n_cols ||
    block_id_zero.n_rows != beta.n_rows ||
    block_id_zero.n_cols != beta.n_cols ||
    sigma2.n_elem != beta.n_cols ||
    gamma2.n_elem < 1u ||
    !beta.is_finite() ||
    !sigma2.is_finite() ||
    arma::any(sigma2 <= 0.0) ||
    !gamma2.is_finite() ||
    arma::any(gamma2 <= 0.0)
  ) {
    Rcpp::stop(
      "Inconsistent dimensions in the local-scale update"
    );
  }

  require_positive_finite(
    tau2,
    "tau2"
  );

  require_positive_finite(
    lambda_shape,
    "lambda_shape"
  );

  require_positive_finite(
    lambda_prior_scale,
    "lambda_prior_scale"
  );

  const double conditional_shape =
    lambda_shape + 0.5;

  for (
    arma::uword equation = 0u;
    equation < beta.n_cols;
    ++equation
  ) {
    for (
      arma::uword row = 0u;
      row < beta.n_rows;
      ++row
    ) {

      const arma::uword block =
        static_cast<arma::uword>(
          block_id_zero(
            row,
            equation
          )
        );

      if (block >= gamma2.n_elem) {
        Rcpp::stop(
          "Invalid block index in the local-scale update"
        );
      }

      const double coefficient =
        beta(
          row,
          equation
        );

      const double conditional_scale =
        lambda_prior_scale +
        coefficient *
        coefficient /
        (
          2.0 *
          sigma2[equation] *
          tau2 *
          gamma2[block]
        );

      lambda2(
        row,
        equation
      ) =
        rinvgamma(
          conditional_shape,
          conditional_scale
        );
    }
  }
}


// ============================================================
// 10. Block-coefficient ASIS
// ============================================================

// Define h_bj = gamma_b^2 lambda_bj^2. Conditional on h in block b,
//
//   gamma_b^2 | h
//     ~ Gamma(
//         a + d_b b,
//         rate = gamma_rate + lambda_prior_scale sum_j 1/h_bj
//       ).
//
// After the gamma draw, lambda_bj^2 = h_bj / gamma_b^2. Thus every
// product gamma_b^2 lambda_bj^2 is unchanged up to floating-point
// roundoff.

struct BlockCoefficientAsisWorkspace {

  arma::vec h;

  explicit BlockCoefficientAsisWorkspace(
      const arma::uword maximum_block_size) :
    h(maximum_block_size, arma::fill::zeros) {}
};


struct BlockCoefficientAsisResult {

  std::uint64_t block_draws_attempted;
  std::uint64_t successful_block_draws;
  std::uint64_t nonfinite_failures;
  std::uint64_t invariant_failures;

  double absolute_log_movement_sum;
  double maximum_absolute_log_movement;
  double maximum_relative_invariant_error;

  BlockCoefficientAsisResult() :
    block_draws_attempted(0u),
    successful_block_draws(0u),
    nonfinite_failures(0u),
    invariant_failures(0u),
    absolute_log_movement_sum(0.0),
    maximum_absolute_log_movement(0.0),
    maximum_relative_invariant_error(0.0) {}
};


inline BlockCoefficientAsisResult
update_block_coefficient_asis(
    arma::vec& gamma2,
    arma::mat& lambda2,
    const BlockMembershipCache& cache,
    const double gamma_shape,
    const double gamma_rate,
    const double lambda_shape,
    const double lambda_prior_scale,
    BlockCoefficientAsisWorkspace& workspace,
    const double invariant_tolerance =
      asis_invariant_tolerance) {

  if (
    gamma2.n_elem !=
      cache.coefficient_count.n_elem ||
    lambda2.n_rows != cache.n_rows ||
    lambda2.n_cols != cache.n_equations
  ) {
    Rcpp::stop(
      "Inconsistent dimensions in block-coefficient ASIS"
    );
  }

  require_positive_finite(
    gamma_shape,
    "gamma_shape"
  );

  require_positive_finite(
    gamma_rate,
    "gamma_rate"
  );

  require_positive_finite(
    lambda_shape,
    "lambda_shape"
  );

  require_positive_finite(
    lambda_prior_scale,
    "lambda_prior_scale"
  );

  require_positive_finite(
    invariant_tolerance,
    "invariant_tolerance"
  );

  BlockCoefficientAsisResult result;

  for (
    arma::uword block = 0u;
    block < gamma2.n_elem;
    ++block
  ) {

    ++result.block_draws_attempted;

    const arma::uvec& indices =
      cache.linear_indices[
        static_cast<std::size_t>(block)
      ];

    if (workspace.h.n_elem < indices.n_elem) {
      workspace.h.set_size(indices.n_elem);
    }

    const double gamma_before =
      gamma2[block];

    if (
      !std::isfinite(gamma_before) ||
      gamma_before <= 0.0
    ) {
      Rcpp::stop(
        "Invalid gamma state before block-coefficient ASIS"
      );
    }

    double inverse_h_sum = 0.0;

    for (
      arma::uword position = 0u;
      position < indices.n_elem;
      ++position
    ) {

      const arma::uword linear =
        indices[position];

      const double lambda_before =
        lambda2[linear];

      if (
        !std::isfinite(lambda_before) ||
        lambda_before <= 0.0
      ) {
        Rcpp::stop(
          "Invalid lambda state before block-coefficient ASIS"
        );
      }

      const double h =
        gamma_before * lambda_before;

      if (!std::isfinite(h) || h <= 0.0) {
        Rcpp::stop(
          "Invalid h invariant before block-coefficient ASIS"
        );
      }

      workspace.h[position] = h;
      inverse_h_sum += 1.0 / h;
    }

    const double conditional_shape =
      gamma_shape +
      static_cast<double>(indices.n_elem) *
      lambda_shape;

    const double conditional_rate =
      gamma_rate +
      lambda_prior_scale * inverse_h_sum;

    double gamma_after;

    if (
      !try_rgamma_rate(
        conditional_shape,
        conditional_rate,
        gamma_after
      )
    ) {
      ++result.nonfinite_failures;
      continue;
    }

    bool proposed_state_valid = true;
    double maximum_block_error = 0.0;

    for (
      arma::uword position = 0u;
      position < indices.n_elem;
      ++position
    ) {

      const double lambda_after =
        workspace.h[position] / gamma_after;

      if (
        !std::isfinite(lambda_after) ||
        lambda_after <= 0.0
      ) {
        proposed_state_valid = false;
        break;
      }

      const double reconstructed_h =
        gamma_after * lambda_after;

      const double relative_error =
        std::abs(
          reconstructed_h -
          workspace.h[position]
        ) /
        std::max(
          std::numeric_limits<double>::min(),
          std::abs(workspace.h[position])
        );

      if (!std::isfinite(relative_error)) {
        proposed_state_valid = false;
        break;
      }

      maximum_block_error =
        std::max(
          maximum_block_error,
          relative_error
        );
    }

    if (!proposed_state_valid) {
      ++result.nonfinite_failures;
      continue;
    }

    result.maximum_relative_invariant_error =
      std::max(
        result.maximum_relative_invariant_error,
        maximum_block_error
      );

    if (
      maximum_block_error >
      invariant_tolerance
    ) {
      ++result.invariant_failures;
      continue;
    }

    gamma2[block] = gamma_after;

    for (
      arma::uword position = 0u;
      position < indices.n_elem;
      ++position
    ) {
      lambda2[indices[position]] =
        workspace.h[position] / gamma_after;
    }

    const double absolute_log_movement =
      std::abs(
        std::log(gamma_after) -
        std::log(gamma_before)
      );

    result.absolute_log_movement_sum +=
      absolute_log_movement;

    result.maximum_absolute_log_movement =
      std::max(
        result.maximum_absolute_log_movement,
        absolute_log_movement
      );

    ++result.successful_block_draws;
  }

  return result;
}


// ============================================================
// 11. Global-block ASIS
// ============================================================

// Define u_b = tau^2 gamma_b^2. Conditional on psi_tau and all u_b,
//
//   tau^2 | u, psi_tau
//     ~ IG(
//         nu/2 + B a,
//         nu/psi_tau + gamma_rate sum_b u_b
//       ).
//
// After the tau draw, gamma_b^2 = u_b / tau^2. Thus every product
// tau^2 gamma_b^2 is unchanged up to floating-point roundoff.

struct GlobalBlockAsisWorkspace {

  arma::vec u;

  explicit GlobalBlockAsisWorkspace(
      const arma::uword n_blocks) :
    u(n_blocks, arma::fill::zeros) {}
};


struct GlobalBlockAsisResult {

  bool attempted;
  bool successful;

  std::uint64_t nonfinite_failures;
  std::uint64_t invariant_failures;

  double tau2_before;
  double tau2_after;
  double absolute_log_movement;
  double maximum_relative_invariant_error;

  GlobalBlockAsisResult() :
    attempted(false),
    successful(false),
    nonfinite_failures(0u),
    invariant_failures(0u),
    tau2_before(NA_REAL),
    tau2_after(NA_REAL),
    absolute_log_movement(NA_REAL),
    maximum_relative_invariant_error(0.0) {}
};


inline GlobalBlockAsisResult update_global_block_asis(
    double& tau2,
    arma::vec& gamma2,
    const double psi_tau,
    const double tau_df,
    const double gamma_shape,
    const double gamma_rate,
    GlobalBlockAsisWorkspace& workspace,
    const double invariant_tolerance =
      asis_invariant_tolerance) {

  if (workspace.u.n_elem != gamma2.n_elem) {
    workspace.u.set_size(gamma2.n_elem);
  }

  require_positive_finite(tau2, "tau2");

  require_positive_finite(
    psi_tau,
    "psi_tau"
  );

  require_positive_finite(
    tau_df,
    "tau_df"
  );

  require_positive_finite(
    gamma_shape,
    "gamma_shape"
  );

  require_positive_finite(
    gamma_rate,
    "gamma_rate"
  );

  require_positive_finite(
    invariant_tolerance,
    "invariant_tolerance"
  );

  GlobalBlockAsisResult result;

  result.attempted = true;
  result.tau2_before = tau2;

  double u_sum = 0.0;

  for (
    arma::uword block = 0u;
    block < gamma2.n_elem;
    ++block
  ) {

    if (
      !std::isfinite(gamma2[block]) ||
      gamma2[block] <= 0.0
    ) {
      Rcpp::stop(
        "Invalid gamma state before global-block ASIS"
      );
    }

    const double u =
      tau2 * gamma2[block];

    if (!std::isfinite(u) || u <= 0.0) {
      Rcpp::stop(
        "Invalid u invariant before global-block ASIS"
      );
    }

    workspace.u[block] = u;
    u_sum += u;
  }

  const double conditional_shape =
    0.5 * tau_df +
    static_cast<double>(gamma2.n_elem) *
    gamma_shape;

  const double conditional_scale =
    tau_df / psi_tau +
    gamma_rate * u_sum;

  double tau2_after;

  if (
    !try_rinvgamma(
      conditional_shape,
      conditional_scale,
      tau2_after
    )
  ) {
    ++result.nonfinite_failures;
    result.tau2_after = tau2;
    return result;
  }

  arma::vec proposed_gamma(gamma2.n_elem);

  bool proposed_state_valid = true;
  double maximum_error = 0.0;

  for (
    arma::uword block = 0u;
    block < gamma2.n_elem;
    ++block
  ) {

    proposed_gamma[block] =
      workspace.u[block] / tau2_after;

    if (
      !std::isfinite(proposed_gamma[block]) ||
      proposed_gamma[block] <= 0.0
    ) {
      proposed_state_valid = false;
      break;
    }

    const double reconstructed_u =
      tau2_after * proposed_gamma[block];

    const double relative_error =
      std::abs(
        reconstructed_u -
        workspace.u[block]
      ) /
      std::max(
        std::numeric_limits<double>::min(),
        std::abs(workspace.u[block])
      );

    if (!std::isfinite(relative_error)) {
      proposed_state_valid = false;
      break;
    }

    maximum_error =
      std::max(
        maximum_error,
        relative_error
      );
  }

  result.maximum_relative_invariant_error =
    maximum_error;

  if (!proposed_state_valid) {
    ++result.nonfinite_failures;
    result.tau2_after = tau2;
    return result;
  }

  if (maximum_error > invariant_tolerance) {
    ++result.invariant_failures;
    result.tau2_after = tau2;
    return result;
  }

  result.tau2_after = tau2_after;

  result.absolute_log_movement =
    std::abs(
      std::log(tau2_after) -
      std::log(tau2)
    );

  tau2 = tau2_after;
  gamma2 = proposed_gamma;

  result.successful = true;

  return result;
}


// ============================================================
// 12. Complete effective-scale helpers
// ============================================================

inline double global_block_scale(
    const double tau2,
    const double gamma2) {

  const double value =
    tau2 * gamma2;

  if (!std::isfinite(value) || value <= 0.0) {
    Rcpp::stop(
      "Invalid global-block scale"
    );
  }

  return value;
}


inline double effective_coefficient_variance(
    const double tau2,
    const double gamma2,
    const double lambda2) {

  const double value =
    tau2 * gamma2 * lambda2;

  if (!std::isfinite(value) || value <= 0.0) {
    Rcpp::stop(
      "Invalid complete effective coefficient variance"
    );
  }

  return value;
}


}  // namespace gigg_bvar


#endif