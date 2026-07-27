test_that("sparkb returns a complete result table", {
  set.seed(2)
  location <- cbind(x = runif(30), y = runif(30))
  counts <- rbind(
    signal = rbinom(30, 1, plogis(5 * (location[, 1] - 0.5))),
    null = rbinom(30, 1, 0.35),
    absent = rep(0, 30)
  )
  result <- sparkb(counts, location, fit_maxiter = 50)

  expect_s3_class(result, "sparkb_result")
  expect_equal(nrow(result), 3)
  expect_true(all(c("pvalue", "padj_bh", "status", "detected_locations") %in% names(result)))
  expect_true(all(paste0("pvalue_gaussian_", 1:5) %in% names(result)))
  expect_equal(result["absent", "status"], "constant")
  expect_length(attr(result, "bandwidths"), 5)
})

test_that("sparkb validates the input shape", {
  expect_error(sparkb(matrix(0, 2, 10), matrix(1, 10, 1)), "exactly two columns")
  expect_error(sparkb(matrix(0, 2, 10), matrix(1, 9, 2)), "must equal")
})

test_that("sparkb accepts dgCMatrix input", {
  set.seed(3)
  location <- cbind(x = runif(12), y = runif(12))
  counts <- Matrix::Matrix(rbind(rbinom(12, 1, .4), rbinom(12, 1, .5)), sparse = TRUE)
  expect_s3_class(sparkb(counts, location, fit_maxiter = 30), "sparkb_result")
})
