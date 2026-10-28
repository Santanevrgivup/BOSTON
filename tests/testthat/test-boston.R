test_that("boston returns a complete result table", {
  set.seed(2)
  location <- cbind(x = runif(30), y = runif(30))
  counts <- rbind(
    signal = rbinom(30, 1, plogis(5 * (location[, 1] - 0.5))),
    null = rbinom(30, 1, 0.35),
    absent = rep(0, 30)
  )
  ids <- paste0("spot_", seq_len(ncol(counts)))
  colnames(counts) <- ids
  rownames(location) <- ids
  result <- boston(counts, location, bandwidths = c(.1, .2, .4, .8, 1.6),
    fit_maxiter = 50)

  expect_s3_class(result, "boston_result")
  expect_equal(nrow(result), 3)
  expect_true(all(c("pvalue", "padj_BY", "status", "detected_locations") %in% names(result)))
  expect_true(all(paste0("pvalue_gaussian_", 1:5) %in% names(result)))
  expect_equal(result["absent", "status"], "constant")
  expect_true(is.na(result["absent", "pvalue"]))
  expect_true(all(is.finite(result[c("signal", "null"), "pvalue"])))
  expect_length(attr(result, "bandwidths"), 5)
})

test_that("default pilot and explicitly supplied pilot give identical results", {
  set.seed(91)
  location <- cbind(runif(80), runif(80))
  counts <- rbind(a = rbinom(80, 1, .3), b = rbinom(80, 1, .6))
  colnames(counts) <- rownames(location) <- paste0("spot_", 1:80)
  bw <- adaptive_bandwidths(location, random_seed = 91L)
  automatic <- boston(counts, location, seed = 91L)
  explicit <- boston(counts, location, bandwidths = bw)
  expect_equal(automatic, explicit)
  bh <- boston(counts, location, bandwidths = bw, padj_method = "BH")
  expect_equal(bh$padj_BH, p.adjust(bh$pvalue, "BH"))
})

test_that("boston validates the input shape", {
  expect_error(boston(matrix(0, 2, 10), matrix(1, 10, 1)), "exactly two columns")
  expect_error(boston(matrix(0, 2, 10), matrix(1, 9, 2)), "must equal")
})

test_that("boston accepts dgCMatrix input", {
  set.seed(3)
  location <- cbind(x = runif(12), y = runif(12))
  counts <- Matrix::Matrix(rbind(rbinom(12, 1, .4), rbinom(12, 1, .5)), sparse = TRUE)
  ids <- paste0("spot_", seq_len(ncol(counts)))
  colnames(counts) <- ids
  rownames(location) <- ids
  bw <- suppressWarnings(adaptive_bandwidths(location, k_local = 8L))
  expect_s3_class(boston(counts, location, bandwidths = bw, fit_maxiter = 30), "boston_result")
})

test_that("boston rejects unverified or misaligned spatial IDs", {
  set.seed(4)
  counts <- rbind(feature_1 = rbinom(12, 1, .4), feature_2 = rbinom(12, 1, .5))
  location <- cbind(x = runif(12), y = runif(12))
  ids <- paste0("spot_", seq_len(ncol(counts)))
  colnames(counts) <- ids
  rownames(location) <- ids

  expect_error(boston(counts, location[rev(ids), , drop = FALSE]), "different orders")
  rownames(location)[1] <- "not_a_spot"
  expect_error(boston(counts, location), "different location IDs")
  rownames(location) <- NULL
  expect_error(boston(counts, location), "supply both")
})

test_that("duplicate feature IDs are made unique without altering punctuation", {
  set.seed(5)
  counts <- rbind(rbinom(12, 1, .4), rbinom(12, 1, .5))
  rownames(counts) <- c("chr1-1-100", "chr1-1-100")
  location <- cbind(x = runif(12), y = runif(12))
  ids <- paste0("spot_", seq_len(ncol(counts)))
  colnames(counts) <- ids
  rownames(location) <- ids

  bw <- suppressWarnings(adaptive_bandwidths(location, k_local = 8L))
  expect_warning(result <- boston(counts, location, bandwidths = bw, fit_maxiter = 30), "made unique")
  expect_setequal(rownames(result), c("chr1-1-100", "chr1-1-100_1"))
})
