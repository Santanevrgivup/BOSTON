test_that("adaptive bandwidths are valid", {
  set.seed(1)
  location <- cbind(runif(100), runif(100))
  bandwidths <- adaptive_bandwidths(location)

  expect_length(bandwidths, 5)
  expect_true(all(is.finite(bandwidths)))
  expect_true(all(bandwidths > 0))
  expect_true(!is.unsorted(bandwidths, strictly = TRUE))
})

test_that("duplicate locations are handled when locations are not all duplicated", {
  location <- cbind(x = c(0, 0, seq_len(18)), y = c(0, 0, rep(1, 18)))
  expect_length(suppressWarnings(adaptive_bandwidths(location)), 5)
})

test_that("pilot is reproducible and does not change R's RNG state", {
  set.seed(42)
  location <- cbind(runif(1100), runif(1100))
  before <- .Random.seed
  bw <- adaptive_bandwidths(location)
  expect_identical(.Random.seed, before)
  expect_identical(adaptive_bandwidths(location), bw)
  expect_equal(adaptive_bandwidths(location[1100:1, ]), bw)
  expect_equal(adaptive_bandwidths(rbind(location, location[1:20, ])), bw)
  expect_equal(adaptive_bandwidths(location * 100), bw * 100, tolerance = 1e-10)
  expect_equal(adaptive_bandwidths(sweep(location, 2, c(20, -5), "+")), bw,
    tolerance = 1e-10)
})

test_that("local radii match brute-force distinct-coordinate neighbours", {
  set.seed(21)
  location <- cbind(runif(100), runif(100))
  pilot <- adaptive_bandwidths(location, diagnostics = TRUE)
  u <- location[pilot$unique_coordinate_index, ]
  distances <- as.matrix(dist(u))
  diag(distances) <- Inf
  radii <- apply(distances, 1, function(x) sort(x)[16])
  expect_equal(as.numeric(pilot$local_radii), radii, ignore_attr = TRUE)
  expect_equal(pilot$local_candidate, median(radii) / 2)
  expect_equal(pilot$global_candidate, median(pilot$pair_distances) / 2)
  expect_equal(unname(pilot$kparam), exp(seq(log(pilot$local_candidate),
    log(pilot$global_candidate), length.out = 5)))
})

test_that("small inputs and invalid pilot parameters have explicit errors", {
  location <- cbind(seq_len(12), rep(0, 12))
  expect_error(adaptive_bandwidths(location), "Too few distinct")
  expect_length(suppressWarnings(adaptive_bandwidths(location, k_local = 3)), 5)
  expect_error(adaptive_bandwidths(location, local_quantile = c(.2, .5)), "scalars")
  expect_error(adaptive_bandwidths(location, random_seed = NA_real_), "integer")
  expect_error(adaptive_bandwidths(location, n_pairs = 0), "integer")
  expect_error(adaptive_bandwidths(matrix(0, 20, 2)), "Too few distinct")
})
