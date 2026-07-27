test_that("adaptive bandwidths are valid", {
  set.seed(1)
  location <- cbind(runif(30), runif(30))
  bandwidths <- adaptive_bandwidths(location)

  expect_length(bandwidths, 5)
  expect_true(all(is.finite(bandwidths)))
  expect_true(all(bandwidths > 0))
  expect_true(!is.unsorted(bandwidths, strictly = TRUE))
})

test_that("duplicate locations are handled when locations are not all duplicated", {
  location <- cbind(x = c(0, 0, seq_len(18)), y = c(0, 0, rep(1, 18)))
  expect_length(adaptive_bandwidths(location), 5)
})
