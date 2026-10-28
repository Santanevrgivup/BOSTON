#include <Rcpp.h>
#include "bandwidth_knn_core.hpp"

// [[Rcpp::plugins(cpp11)]]

// kparam is sigma in exp(-distance^2/(2*sigma^2)).
// [[Rcpp::export(rng = false)]]
Rcpp::List ComputeAdaptiveBandwidthsKNN(
    Rcpp::NumericMatrix xin,
    int random_seed = 123,
    int k_local = 16,
    int n_query = 1024,
    int n_pairs = 50000,
    double local_quantile = 0.50,
    double global_quantile = 0.50) {
  if (random_seed < 0)
    Rcpp::stop("random_seed must be a nonnegative integer.");
  std::vector<boston_bw_knn::Point> x(
      xin.nrow(), boston_bw_knn::Point(xin.ncol()));
  for (int i = 0; i < xin.nrow(); ++i)
    for (int j = 0; j < xin.ncol(); ++j)
      x[i][j] = xin(i, j);

  const auto result = boston_bw_knn::compute(
      x, static_cast<std::uint32_t>(random_seed), k_local, n_query, n_pairs,
      local_quantile, global_quantile);

  Rcpp::IntegerVector unique_index(result.unique_indices.size());
  for (std::size_t i = 0; i < result.unique_indices.size(); ++i)
    unique_index[i] = static_cast<int>(result.unique_indices[i] + 1);

  return Rcpp::List::create(
      Rcpp::Named("kparam") = result.bandwidths,
      Rcpp::Named("n_locations") = xin.nrow(),
      Rcpp::Named("n_distinct_coordinates") =
          static_cast<double>(result.unique_indices.size()),
      Rcpp::Named("unique_coordinate_index") = unique_index,
      Rcpp::Named("k_local") = k_local,
      Rcpp::Named("local_quantile") = local_quantile,
      Rcpp::Named("global_quantile") = global_quantile,
      Rcpp::Named("local_candidate") = result.local_candidate,
      Rcpp::Named("global_candidate") = result.global_candidate,
      Rcpp::Named("range_expanded") = result.range_expanded,
      Rcpp::Named("local_radii") = result.local_radii,
      Rcpp::Named("pair_distances") = result.pair_distances);
}
