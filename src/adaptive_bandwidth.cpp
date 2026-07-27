#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(openmp)]]

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
SEXP ComputeAdaptiveBandwidths(SEXP xin) {
  try {
    const arma::mat& x = as<arma::mat>(xin);
    const int n = x.n_rows;
    if (n < 10) Rcpp::stop("At least 10 spatial locations are required.");
    if (x.n_cols < 2) Rcpp::stop("Spatial locations must have at least two columns.");

    // Preserve the original global-distance construction for the middle scales.
    double lmin_val = 1e10;
    double lmax_val = 0.0;
    for (int i = 0; i < n - 1; ++i) {
      const double xi = x(i, 0), yi = x(i, 1);
      for (int j = i + 1; j < n; ++j) {
        const double d_sq = std::pow(xi - x(j, 0), 2.0) +
          std::pow(yi - x(j, 1), 2.0);
        if (d_sq > 1e-16) {
          const double d = std::sqrt(d_sq);
          lmin_val = std::min(lmin_val, d);
          lmax_val = std::max(lmax_val, d);
        }
      }
    }
    if (!std::isfinite(lmin_val) || lmax_val <= 0.0) {
      Rcpp::stop("Spatial locations must contain at least two distinct points.");
    }

    const double old_lmin = lmin_val / 2.0;
    const double old_lmax = lmax_val * 2.0;
    const double log_low_old = std::log10(old_lmin);
    const double log_high_old = std::log10(old_lmax);
    const double sigma_3_orig = std::pow(10.0,
      log_low_old + (2.0 / 9.0) * (log_high_old - log_low_old));
    const double sigma_7_orig = std::pow(10.0,
      log_low_old + (6.0 / 9.0) * (log_high_old - log_low_old));

    // Estimate local sampling density with the median 8-nearest-neighbour radius.
    const int n_sample = std::min(n, 1000);
    arma::uvec sample_idx = arma::randperm(n, n_sample);
    arma::vec d8_sq(n_sample);

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n_sample; ++i) {
      const int idx_i = sample_idx(i);
      const double xi = x(idx_i, 0), yi = x(idx_i, 1);
      arma::vec dists_sq(n);
      for (int j = 0; j < n; ++j) {
        const double dx = xi - x(j, 0), dy = yi - x(j, 1);
        dists_sq[j] = dx * dx + dy * dy;
      }
      std::nth_element(dists_sq.begin(), dists_sq.begin() + 8, dists_sq.end());
      d8_sq[i] = dists_sq[8];
    }

    const double sigma_k8 = std::sqrt(arma::median(d8_sq)) / 2.0;
    const double final_low = std::max(sigma_3_orig, sigma_k8);
    double final_high = sigma_7_orig;
    if (final_low >= final_high) final_high = final_low * 2.0;

    const arma::vec kparam = arma::exp10(arma::linspace(
      std::log10(final_low), std::log10(final_high), 5));
    return List::create(Named("kparam") = kparam);
  } catch (std::exception& ex) {
    forward_exception_to_r(ex);
  }
  return R_NilValue;
}
