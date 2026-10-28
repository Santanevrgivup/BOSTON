#ifndef BOSTON_BANDWIDTH_KNN_CORE_HPP
#define BOSTON_BANDWIDTH_KNN_CORE_HPP

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>

namespace boston_bw_knn {

using Point = std::vector<double>;

struct Result {
  std::vector<double> bandwidths;
  std::vector<double> local_radii;
  std::vector<double> pair_distances;
  std::vector<std::size_t> unique_indices;
  double local_candidate;
  double global_candidate;
  bool range_expanded;
};

inline double distance2(const Point& a, const Point& b) {
  double ans = 0.0;
  for (std::size_t d = 0; d < a.size(); ++d) {
    const double delta = a[d] - b[d];
    ans += delta * delta;
  }
  return ans;
}

inline double quantile(std::vector<double> values, double probability) {
  if (values.empty()) throw std::runtime_error("The distance sample is empty.");
  if (!(probability > 0.0 && probability < 1.0))
    throw std::runtime_error("The distance quantile must be between 0 and 1.");
  std::sort(values.begin(), values.end());
  const double position = (values.size() - 1) * probability;
  const std::size_t lower = static_cast<std::size_t>(position);
  const std::size_t upper = std::min(lower + 1, values.size() - 1);
  return values[lower] + (position - lower) * (values[upper] - values[lower]);
}

inline Result compute(const std::vector<Point>& input,
                      std::uint32_t seed = 123,
                      int k_local = 16,
                      int n_query = 1024,
                      int n_pairs = 50000,
                      double local_quantile = 0.50,
                      double global_quantile = 0.50) {
  if (input.size() < 10)
    throw std::runtime_error("At least 10 spatial locations are required.");
  const std::size_t dimensions = input[0].size();
  if (dimensions != 2 && dimensions != 3)
    throw std::runtime_error("Use exactly two or three coordinate columns.");
  if (k_local < 1 || n_query < 1 || n_pairs < 1 ||
      !(local_quantile > 0.0 && local_quantile < 1.0) ||
      !(global_quantile > 0.0 && global_quantile < 1.0))
    throw std::runtime_error("Invalid kNN bandwidth parameters.");

  for (const auto& row : input) {
    if (row.size() != dimensions)
      throw std::runtime_error("Inconsistent coordinate dimensions.");
    for (double value : row)
      if (!std::isfinite(value))
        throw std::runtime_error("Coordinates must be finite.");
  }

  // Exact duplicates are collapsed only for pilot bandwidth estimation. The
  // original observations remain available for the downstream kernel matrix.
  std::map<Point, std::size_t> first_index;
  for (std::size_t i = 0; i < input.size(); ++i)
    if (first_index.find(input[i]) == first_index.end())
      first_index.emplace(input[i], i);

  std::vector<Point> points;
  Result result;
  points.reserve(first_index.size());
  result.unique_indices.reserve(first_index.size());
  for (const auto& entry : first_index) {
    points.push_back(entry.first);
    result.unique_indices.push_back(entry.second);
  }
  const std::size_t n = points.size();
  if (n <= static_cast<std::size_t>(k_local))
    throw std::runtime_error(
      "Too few distinct coordinates for k_local; reduce k_local or check the coordinate data.");

  std::mt19937 generator(seed);
  std::vector<std::size_t> query_order(n);
  std::iota(query_order.begin(), query_order.end(), 0);
  const std::size_t n_queries = std::min(n, static_cast<std::size_t>(n_query));
  if (n_queries < n) {
    for (std::size_t i = 0; i < n_queries; ++i) {
      std::uniform_int_distribution<std::size_t> pick(i, n - 1);
      std::swap(query_order[i], query_order[pick(generator)]);
    }
  }

  std::vector<double> squared_distances;
  squared_distances.reserve(n - 1);
  for (std::size_t q = 0; q < n_queries; ++q) {
    const std::size_t i = query_order[q];
    squared_distances.clear();
    for (std::size_t j = 0; j < n; ++j)
      if (j != i) squared_distances.push_back(distance2(points[i], points[j]));
    std::nth_element(squared_distances.begin(),
                     squared_distances.begin() + (k_local - 1),
                     squared_distances.end());
    result.local_radii.push_back(std::sqrt(squared_distances[k_local - 1]));
  }

  // Robust global upper scale. Pairs are sampled uniformly from distinct
  // coordinates. Exact coordinate multiplicities do not add pair weight;
  // spatially dense clusters of different coordinates still affect the pilot.
  std::uniform_int_distribution<std::size_t> first_pick(0, n - 1);
  std::uniform_int_distribution<std::size_t> second_pick(0, n - 2);
  result.pair_distances.reserve(n_pairs);
  for (int b = 0; b < n_pairs; ++b) {
    const std::size_t i = first_pick(generator);
    std::size_t j = second_pick(generator);
    if (j >= i) ++j;
    result.pair_distances.push_back(std::sqrt(distance2(points[i], points[j])));
  }

  // The factor two follows the convention K(d)=exp[-d^2/(2*sigma^2)]:
  // a distance of 2*sigma has kernel weight exp(-2).
  const double lower = quantile(result.local_radii, local_quantile) / 2.0;
  double upper = quantile(result.pair_distances, global_quantile) / 2.0;
  if (!(lower > 0.0) || !std::isfinite(lower) ||
      !(upper > 0.0) || !std::isfinite(upper))
    throw std::runtime_error("The coordinate geometry did not yield positive finite scales.");

  result.local_candidate = lower;
  result.global_candidate = upper;
  result.range_expanded = upper <= lower;
  if (result.range_expanded) upper = 2.0 * lower;

  for (int j = 0; j < 5; ++j) {
    const double fraction = j / 4.0;
    const double sigma = std::exp(std::log(lower) +
                                  fraction * (std::log(upper) - std::log(lower)));
    if (!std::isfinite(sigma) || sigma <= 0.0 ||
        (!result.bandwidths.empty() && sigma <= result.bandwidths.back()))
      throw std::runtime_error("Five distinct positive bandwidths cannot be represented.");
    result.bandwidths.push_back(sigma);
  }
  return result;
}

} // namespace boston_bw_knn

#endif
