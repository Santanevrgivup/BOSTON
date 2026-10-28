#include <RcppArmadillo.h>

// Apply P through its diagonal and rank-q correction, without allocating P.
static arma::vec project_working(const arma::vec& v, const arma::vec& h,
                                const arma::mat& u, const arma::mat& a) {
  return h % v - u * (a * (u.t() * v));
}

// [[Rcpp::export(rng = false)]]
Rcpp::List BOSTONCovariateAI(const arma::vec& y, const arma::mat& x,
                            const arma::vec& d2, arma::vec tau,
                            bool fix_nugget, double tol) {
  if (x.n_rows != y.n_elem || d2.n_elem != y.n_elem || tau.n_elem != 2 ||
      !y.is_finite() || !x.is_finite() || !d2.is_finite() || arma::any(d2 <= 0))
    Rcpp::stop("Invalid working quantities in covariate null-model update.");
  arma::vec h = 1.0 / (tau(0) / (d2 + 1e-5) + tau(1) + 1e-5);
  arma::mat u = x.each_col() % h;
  arma::mat a;
  if (!arma::inv_sympd(a, x.t() * u))
    Rcpp::stop("Weighted covariate design is singular.");
  arma::vec alpha = a * (u.t() * y);
  arma::vec py = project_working(y, h, u, a);
  arma::vec eta = y - tau(0) * py / d2;
  if (!fix_nugget) {
    const arma::vec papy = project_working(py, h, u, a);
    const double trp = arma::sum(h) - arma::trace(a * (u.t() * u));
    const double ai = arma::dot(py, papy);
    if (!std::isfinite(ai) || ai <= 0)
      Rcpp::stop("Degenerate average information in covariate null model.");
    const double delta = (arma::dot(y, papy) - trp) / ai;
    double step = 1.0;
    double candidate = tau(1) + delta;
    for (int backtrack = 0; candidate < 0 && backtrack < 64; ++backtrack) {
      step *= 0.5;
      candidate = tau(1) + step * delta;
      if (tau(1) < tol && candidate < tol) candidate = 0;
    }
    if (!std::isfinite(candidate) || candidate < 0)
      Rcpp::stop("Covariate variance update failed to find a finite nonnegative step.");
    tau(1) = candidate < tol ? 0 : candidate;
  }
  return Rcpp::List::create(Rcpp::Named("tau") = tau,
    Rcpp::Named("alpha") = alpha, Rcpp::Named("eta") = eta);
}

// Row-block Gaussian/custom-kernel score moments for an arbitrary design X.
// The two moments intentionally match the unadjusted quadratic statistic in
// the published working-model derivation, not a nuisance-adjusted variance.
// [[Rcpp::export(rng = false)]]
Rcpp::List BOSTONCovariateScore(const arma::mat& y, const arma::mat& location,
                               double sigma, const arma::mat& d2,
                               const arma::mat& tau, const arma::mat& x,
                               SEXP kernel = R_NilValue, int block_size = 64) {
  const arma::uword n = y.n_rows, m = y.n_cols, q = x.n_cols;
  const bool custom = !Rf_isNull(kernel);
  if (d2.n_rows != n || d2.n_cols != m || tau.n_rows != 2 || tau.n_cols != m ||
      x.n_rows != n || location.n_rows != n || location.n_cols != 2 ||
      block_size < 1 || q < 1 || !y.is_finite() || !d2.is_finite() ||
      !tau.is_finite() || !x.is_finite() || arma::any(arma::vectorise(d2) <= 0))
    Rcpp::stop("Invalid dimensions or working quantities for kernel score test.");
  arma::mat supplied;
  if (custom) {
    supplied = Rcpp::as<arma::mat>(kernel);
    if (supplied.n_rows != n || supplied.n_cols != n || !supplied.is_finite())
      Rcpp::stop("Custom kernel has invalid dimensions or values.");
  } else if (!std::isfinite(sigma) || sigma <= 0) {
    Rcpp::stop("Gaussian bandwidth must be positive and finite.");
  }

  arma::mat h(n, m), py(n, m);
  std::vector<arma::mat> us(m), as(m), t(m), v(m);
  arma::vec s0(m, arma::fill::zeros), base(m, arma::fill::zeros);
  for (arma::uword g = 0; g < m; ++g) {
    h.col(g) = 1.0 / (tau(0, g) / (d2.col(g) + 1e-5) + tau(1, g) + 1e-5);
    us[g] = x.each_col() % h.col(g);
    if (!arma::inv_sympd(as[g], x.t() * us[g]))
      Rcpp::stop("Weighted covariate design is singular in score test.");
    py.col(g) = project_working(y.col(g), h.col(g), us[g], as[g]);
    t[g].zeros(q, q);
    v[g].zeros(q, q);
  }
  arma::vec diagonal(n, arma::fill::ones);
  if (custom) diagonal = supplied.diag();
  for (arma::uword start = 0; start < n; start += block_size) {
    Rcpp::checkUserInterrupt();
    const arma::uword end = std::min(n - 1, start + block_size - 1);
    arma::mat k(end - start + 1, n);
    if (custom) k = supplied.rows(start, end);
    else {
      for (arma::uword i = start; i <= end; ++i)
        for (arma::uword j = 0; j < n; ++j) {
          const double dx = location(i, 0) - location(j, 0);
          const double dy = location(i, 1) - location(j, 1);
          k(i - start, j) = std::exp(-(dx * dx + dy * dy) / (2 * sigma * sigma));
        }
    }
    const arma::mat squared = arma::square(k);
    for (arma::uword g = 0; g < m; ++g) {
      const arma::vec hb = h.col(g).rows(start, end);
      const arma::mat ku = k * us[g];
      t[g] += us[g].rows(start, end).t() * ku;
      v[g] += ku.t() * (ku.each_col() % hb);
      base(g) += arma::dot(hb, squared * h.col(g));
      s0(g) += 0.5 * arma::dot(py.col(g).rows(start, end), k * py.col(g));
    }
  }
  arma::vec mean(m), variance(m);
  for (arma::uword g = 0; g < m; ++g) {
    const double diagonal_trace = arma::dot(h.col(g), diagonal);
    mean(g) = 0.5 * (diagonal_trace - arma::trace(as[g] * t[g]));
    variance(g) = 0.5 * (base(g) - 2 * arma::trace(as[g] * v[g]) +
      arma::trace(as[g] * t[g] * as[g] * t[g]));
    if (mean(g) <= 1e-10 * std::abs(diagonal_trace)) {
      mean(g) = 0;
      variance(g) = 0;
    }
  }
  return Rcpp::List::create(Rcpp::Named("S0") = s0,
    Rcpp::Named("ee") = mean, Rcpp::Named("infoMp1") = variance);
}
