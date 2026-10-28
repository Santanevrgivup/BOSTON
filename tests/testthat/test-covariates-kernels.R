extension_fixture <- function(n = 80L) {
  set.seed(732)
  location <- cbind(x = runif(n), y = runif(n))
  counts <- rbind(a = rbinom(n, 1, .35), b = rbinom(n, 1, .55))
  colnames(counts) <- rownames(location) <- paste0("spot_", seq_len(n))
  covariates <- data.frame(depth = seq(-1, 1, length.out = n),
    batch = factor(rep(c("A", "B"), length.out = n)), row.names = colnames(counts))
  list(counts = counts, location = location, covariates = covariates)
}

test_that("blocked covariate score agrees with explicit independent P algebra", {
  f <- extension_fixture(35L)
  x <- model.matrix(~ ., f$covariates)
  set.seed(51)
  y <- matrix(rnorm(70), 35, 2)
  d2 <- matrix(runif(70, .08, .24), 35, 2)
  tau <- rbind(c(1, 1), c(.1, .3))
  k <- exp(-as.matrix(dist(f$location))^2 / (2 * .3^2))
  built <- BOSTONCovariateScore(y, f$location, .3, d2, tau, x, block_size = 7L)
  custom <- BOSTONCovariateScore(y, f$location, 1, d2, tau, x, k, 11L)
  for (g in 1:2) {
    h <- 1 / (tau[1, g] / (d2[, g] + 1e-5) + tau[2, g] + 1e-5)
    hmat <- diag(h)
    p <- hmat - hmat %*% x %*% solve(t(x) %*% hmat %*% x) %*% t(x) %*% hmat
    pk <- p %*% k
    score <- drop(.5 * t(y[, g]) %*% p %*% k %*% p %*% y[, g])
    mean <- .5 * sum(diag(pk))
    variance <- .5 * sum(diag(pk %*% pk))
    expect_equal(as.numeric(built$S0)[g], score, tolerance = 1e-10)
    expect_equal(as.numeric(built$ee)[g], mean, tolerance = 1e-10)
    expect_equal(as.numeric(built$infoMp1)[g], variance, tolerance = 1e-10)
  }
  expect_equal(built, custom, tolerance = 1e-10)
})

test_that("low-rank null-model update matches inherited dense covariate update", {
  f <- extension_fixture(35L)
  x <- model.matrix(~ ., f$covariates)
  set.seed(12)
  y <- rnorm(35)
  d2 <- runif(35, .1, .24)
  for (fixed in c(FALSE, TRUE)) {
    dense <- CovariatesAI(y, x, d2, c(1, .4), c(1, as.integer(fixed)), 1e-5)
    lowrank <- BOSTONCovariateAI(y, x, d2, c(1, .4), fixed, 1e-5)
    for (field in c("tau", "alpha", "eta")) {
      expect_equal(as.numeric(lowrank[[field]]), as.numeric(dense[[field]]), tolerance = 1e-10)
    }
  }
})

test_that("covariates run with Gaussian and custom kernels", {
  f <- extension_fixture()
  bw <- c(.05, .1, .2, .4, .8)
  gaussian <- boston(f$counts, f$location, covariates = f$covariates, bandwidths = bw)
  kernels <- lapply(bw, function(s) exp(-as.matrix(dist(f$location))^2 / (2 * s^2)))
  names(kernels) <- paste0("gaussian_", 1:5)
  custom <- boston(f$counts, f$location, covariates = f$covariates, kernels = kernels)
  fields <- c(paste0("pvalue_gaussian_", 1:5), "pvalue", "padj_BY")
  expect_true(all(is.finite(gaussian$pvalue)))
  expect_equal(as.matrix(gaussian[, fields]), as.matrix(custom[, fields]), tolerance = 1e-9)
  expect_equal(attr(gaussian, "covariate_columns"), colnames(model.matrix(~ ., f$covariates)))
})

test_that("custom Gaussian kernels preserve the intercept score path", {
  f <- extension_fixture()
  bw <- c(.05, .1, .2, .4, .8)
  kernels <- lapply(bw, function(s) exp(-as.matrix(dist(f$location))^2 / (2 * s^2)))
  names(kernels) <- paste0("gaussian_", 1:5)
  original <- boston(f$counts, f$location, bandwidths = bw)
  custom <- boston(f$counts, f$location, kernels = kernels)
  fields <- c(paste0("pvalue_gaussian_", 1:5), "pvalue", "padj_BY")
  expect_equal(as.matrix(original[, fields]), as.matrix(custom[, fields]), tolerance = 1e-9)
})

test_that("covariate IDs, rank and missing values are validated", {
  f <- extension_fixture()
  expect_error(boston(f$counts, f$location, covariates = f$covariates[80:1, ]), "IDs")
  bad <- f$covariates
  bad$depth[1] <- NA_real_
  expect_error(boston(f$counts, f$location, covariates = bad), "missing")
  bad <- cbind(f$covariates, duplicate = f$covariates$depth)
  expect_error(boston(f$counts, f$location, covariates = bad), "full rank")
  bad <- matrix(1, 80, 1, dimnames = list(colnames(f$counts), "intercept"))
  expect_error(boston(f$counts, f$location, covariates = bad), "intercept")
})

test_that("custom kernel IDs, symmetry and PSD are checked before fitting", {
  f <- extension_fixture(20L)
  k <- diag(20)
  dimnames(k) <- list(colnames(f$counts), colnames(f$counts))
  expect_error(boston(f$counts, f$location, kernels = k[20:1, ]), "IDs")
  bad <- k
  bad[1, 2] <- 1
  expect_error(boston(f$counts, f$location, kernels = bad), "symmetric")
  bad <- k
  bad[1, 1] <- -1
  expect_error(boston(f$counts, f$location, kernels = bad), "semidefinite")
  expect_error(boston(f$counts, f$location, kernels = k, bandwidths = 1:5), "not both")
  expect_error(boston(f$counts, f$location, kernels = k, weights = 0), "weights")
  expect_error(boston(f$counts, f$location, kernels = k, weights = NA_real_), "weights")
  constant <- k * 0 + 1
  out <- boston(f$counts, f$location, kernels = list(constant = constant))
  expect_true(all(is.na(out$pvalue)))
  expect_true(all(out$status == "test_failed"))
})
