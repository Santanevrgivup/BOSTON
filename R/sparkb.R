#' Select five adaptive Gaussian-kernel bandwidths
#'
#' The bandwidth range starts at the larger of the original third Gaussian scale
#' and half the median distance to the eighth nearest spatial neighbour. It ends
#' at the original seventh Gaussian scale; when necessary, the upper endpoint is
#' expanded to preserve an increasing range. Five values equally spaced on the
#' log10 scale are returned.
#'
#' @param location Numeric matrix or data frame with one row per spatial
#'   location and exactly two coordinate columns.
#'
#' @return A named numeric vector of five strictly increasing Gaussian-kernel
#'   bandwidths.
#' @export
#'
#' @examples
#' coords <- cbind(x = seq_len(20), y = rep(0, 20))
#' adaptive_bandwidths(coords)
adaptive_bandwidths <- function(location) {
  location <- .validate_location(location)
  bandwidths <- as.numeric(ComputeAdaptiveBandwidths(location)$kparam)
  names(bandwidths) <- paste0("gaussian_", seq_along(bandwidths))
  bandwidths
}

#' Identify spatially variable features in sparse spatial omics data
#'
#' `sparkb()` binarizes a feature-by-location count matrix, fits a binomial null
#' model to each non-constant feature, tests five Gaussian spatial kernels, and
#' combines their p-values using the aggregated Cauchy association test (ACAT).
#' The default bandwidths are selected adaptively from the supplied coordinates.
#'
#' The scalable score-test implementation currently supports an intercept-only
#' null model. `covariates` is retained in the interface for forward compatibility
#' and must be `NULL` in this release.
#'
#' @param counts A numeric matrix or a `Matrix::dgCMatrix` with features in rows
#'   and spatial locations in columns. Values must be finite and non-negative.
#' @param location Numeric matrix or data frame with exactly two coordinate
#'   columns and one row for every column of `counts`.
#' @param covariates Reserved for a future covariate-adjusted scalable test; must
#'   be `NULL`.
#' @param features Optional feature names or integer indices to test. By default,
#'   all features are tested.
#' @param bandwidths Optional numeric vector of exactly five positive Gaussian
#'   bandwidths. Defaults to [adaptive_bandwidths()].
#' @param detection_threshold A count is treated as detected when it is strictly
#'   greater than this value.
#' @param fit_maxiter Maximum iterations of the binomial null-model fit.
#' @param fit_tol Convergence tolerance of the binomial null-model fit.
#' @param num_cores Number of cores for fitting independent feature-wise null
#'   models. On Windows, fitting is performed serially.
#' @param verbose Logical; report progress messages.
#'
#' @return A `sparkb_result` data frame ordered by combined p-value. It includes
#'   five kernel-specific p-values (`pvalue_gaussian_*`), the ACAT p-value
#'   (`pvalue`), Benjamini-Hochberg FDR (`padj_bh`), and model status. The
#'   bandwidths used are stored in `attr(result, "bandwidths")`.
#' @export
#'
#' @examples
#' set.seed(1)
#' location <- cbind(x = runif(40), y = runif(40))
#' counts <- rbind(
#'   spatial = rbinom(40, 1, plogis(5 * (location[, 1] - 0.5))),
#'   null = rbinom(40, 1, 0.3)
#' )
#' sparkb(counts, location, fit_maxiter = 50)
sparkb <- function(counts,
                   location,
                   covariates = NULL,
                   features = NULL,
                   bandwidths = NULL,
                   detection_threshold = 0,
                   fit_maxiter = 500,
                   fit_tol = 1e-5,
                   num_cores = 1,
                   verbose = FALSE) {
  if (!is.null(covariates)) {
    stop("`covariates` is not supported by the scalable SPARKb score test in this release.",
      call. = FALSE)
  }
  if (!is.numeric(detection_threshold) || length(detection_threshold) != 1L ||
      !is.finite(detection_threshold)) {
    stop("`detection_threshold` must be one finite numeric value.", call. = FALSE)
  }
  if (!is.numeric(fit_maxiter) || length(fit_maxiter) != 1L || fit_maxiter < 1L ||
      !is.numeric(fit_tol) || length(fit_tol) != 1L || fit_tol <= 0 ||
      !is.numeric(num_cores) || length(num_cores) != 1L || num_cores < 1L) {
    stop("`fit_maxiter`, `fit_tol`, and `num_cores` must be positive scalar values.",
      call. = FALSE)
  }

  counts <- .validate_counts(counts)
  location <- .validate_location(location)
  if (ncol(counts) != nrow(location)) {
    stop("`ncol(counts)` must equal `nrow(location)`.", call. = FALSE)
  }
  if (ncol(counts) < 10L) {
    stop("At least 10 spatial locations are required.", call. = FALSE)
  }

  selected <- .select_features(counts, features)
  counts <- counts[selected, , drop = FALSE]
  feature_names <- rownames(counts)

  if (is.null(bandwidths)) {
    bandwidths <- adaptive_bandwidths(location)
  } else {
    bandwidths <- as.numeric(bandwidths)
    if (length(bandwidths) != 5L || any(!is.finite(bandwidths)) ||
        any(bandwidths <= 0) || is.unsorted(bandwidths, strictly = TRUE)) {
      stop("`bandwidths` must contain five finite, strictly increasing positive values.",
        call. = FALSE)
    }
    names(bandwidths) <- paste0("gaussian_", seq_along(bandwidths))
  }

  binary_counts <- Matrix::Matrix(0 + (counts > detection_threshold), sparse = TRUE)
  detected <- Matrix::rowSums(binary_counts)
  testable <- detected > 0 & detected < ncol(binary_counts)
  status <- rep("constant", nrow(binary_counts))
  status[testable] <- "fit_failed"

  pvalues <- matrix(NA_real_, nrow(binary_counts), length(bandwidths),
    dimnames = list(feature_names, paste0("pvalue_", names(bandwidths))))
  if (any(testable)) {
    valid_idx <- which(testable)
    if (verbose) message("Fitting ", length(valid_idx), " binomial null models.")
    fitted <- .map_features(valid_idx, function(i) {
      tryCatch(.sparkb_null_model(as.numeric(binary_counts[i, ]),
        maxiter = as.integer(fit_maxiter), tol = fit_tol),
      error = function(e) structure(list(message = conditionMessage(e)), class = "sparkb_fit_error"))
    }, num_cores)
    good_fit <- !vapply(fitted, inherits, logical(1), what = "sparkb_fit_error")
    if (any(good_fit)) {
      good_idx <- valid_idx[good_fit]
      status[good_idx] <- "ok"
      y_working <- do.call(cbind, lapply(fitted[good_fit], `[[`, "Y"))
      d_squared <- do.call(cbind, lapply(fitted[good_fit], function(x) x$D^2))
      theta <- do.call(cbind, lapply(fitted[good_fit], `[[`, "theta"))

      if (verbose) message("Testing five Gaussian kernels.")
      for (k in seq_along(bandwidths)) {
        kernel_result <- tryCatch(
          fastComputeQuantitiesWithoutCVT(y_working, location, bandwidths[k],
            d_squared, theta, "gaussian"),
          error = function(e) e
        )
        if (inherits(kernel_result, "error")) {
          status[good_idx] <- "test_failed"
          break
        }
        scale <- kernel_result$infoMp1 / (2 * kernel_result$ee)
        df <- 2 * kernel_result$ee^2 / kernel_result$infoMp1
        pvalues[good_idx, k] <- stats::pchisq(kernel_result$S0 / scale, df,
          lower.tail = FALSE)
      }
    }
  }

  combined <- apply(pvalues, 1L, .acat)
  result <- data.frame(
    feature = feature_names,
    pvalues,
    pvalue = combined,
    padj_bh = stats::p.adjust(combined, method = "BH"),
    detected_locations = as.integer(detected),
    status = status,
    row.names = feature_names,
    check.names = FALSE
  )
  result <- result[order(result$pvalue, na.last = TRUE), , drop = FALSE]
  class(result) <- c("sparkb_result", "data.frame")
  attr(result, "bandwidths") <- bandwidths
  attr(result, "detection_threshold") <- detection_threshold
  result
}

.validate_location <- function(location) {
  if (is.data.frame(location)) location <- as.matrix(location)
  if (!is.matrix(location) || !is.numeric(location) || ncol(location) != 2L ||
      nrow(location) < 1L || any(!is.finite(location))) {
    stop("`location` must be a finite numeric matrix or data frame with exactly two columns.",
      call. = FALSE)
  }
  storage.mode(location) <- "double"
  location
}

.validate_counts <- function(counts) {
  if (!is.matrix(counts) && !inherits(counts, "Matrix")) {
    stop("`counts` must be a numeric matrix or Matrix object.", call. = FALSE)
  }
  values <- if (inherits(counts, "Matrix")) counts@x else counts
  if (!is.numeric(values) || any(!is.finite(values)) || any(values < 0)) {
    stop("`counts` must contain finite, non-negative numeric values.", call. = FALSE)
  }
  if (is.null(rownames(counts))) rownames(counts) <- paste0("feature_", seq_len(nrow(counts)))
  if (anyDuplicated(rownames(counts))) {
    stop("`counts` row names must be unique.", call. = FALSE)
  }
  counts
}

.select_features <- function(counts, features) {
  if (is.null(features)) return(seq_len(nrow(counts)))
  if (is.character(features)) {
    missing <- setdiff(features, rownames(counts))
    if (length(missing)) stop("Unknown feature(s): ", paste(missing, collapse = ", "), call. = FALSE)
    return(match(features, rownames(counts)))
  }
  if (is.numeric(features) && all(features == as.integer(features)) &&
      all(features >= 1L) && all(features <= nrow(counts))) return(as.integer(features))
  stop("`features` must be feature names or valid positive integer indices.", call. = FALSE)
}

.map_features <- function(indices, fun, num_cores) {
  if (.Platform$OS.type == "windows" || num_cores <= 1L) return(lapply(indices, fun))
  parallel::mclapply(indices, fun, mc.cores = as.integer(num_cores))
}

.acat <- function(pvalues) {
  pvalues <- pvalues[is.finite(pvalues)]
  if (!length(pvalues)) return(NA_real_)
  pvalues[pvalues <= 0] <- 5.55e-17
  pvalues[pvalues >= 1] <- 0.99
  statistic <- mean(tan((0.5 - pvalues) * pi))
  if (statistic > 1e15) return(1 / (pi * statistic))
  max(5.55e-17, 1 - stats::pcauchy(statistic))
}

.sparkb_null_model <- function(y, maxiter, tol) {
  model0 <- stats::glm(y ~ 1, family = stats::binomial(link = "logit"))
  .sparkb_fit(model0, maxiter = maxiter, tol = tol)
}

.sparkb_fit <- function(model0, maxiter, tol) {
  fit <- .sparkb_ai(model0, num_vc = 1L, maxiter = maxiter, tol = tol)
  fixed_old <- rep(0, 2L)
  fixed_new <- 1 * (fit$theta < 1.01 * tol)
  count <- 1L
  while (any(fixed_new != fixed_old) && count < maxiter) {
    count <- count + 1L
    fixed_old <- fixed_new
    fit <- .sparkb_ai(model0, num_vc = 1L, fixtau = fixed_old,
      maxiter = maxiter, tol = tol)
    fixed_new <- 1 * (fit$theta < 1.01 * tol)
  }
  fit
}

.sparkb_ai <- function(model0, num_vc, tau = rep(0.1, num_vc + 1L),
                       fixtau = rep(0, num_vc + 1L), maxiter, tol) {
  y <- model0$y
  n <- length(y)
  family <- model0$family
  eta <- model0$linear.predictors
  mu <- model0$fitted.values
  derivative <- family$mu.eta(eta)
  d <- derivative / sqrt(family$variance(mu))
  y_working <- eta + (y - mu) / derivative
  x <- stats::model.matrix(model0)
  alpha <- model0$coefficients

  tau[1] <- 1
  fixtau[1] <- 1
  free <- which(fixtau == 0)
  if (length(free)) {
    tau[free] <- min(0.9, stats::var(y_working) / (num_vc + 1L))
    h_inverse <- 1 / (tau[1] / d^2 + tau[2])
    xhix <- sum(h_inverse)
    py_initial <- h_inverse * y_working - h_inverse * sum(h_inverse * y_working) / xhix
    diag_p <- h_inverse - h_inverse^2 / xhix
    papy <- h_inverse * py_initial - h_inverse * sum(h_inverse * py_initial) / xhix
    tau[free] <- max(0, tau[free] + tau[free]^2 *
      (sum(y_working * papy) - sum(diag_p)) / n)
  }
  initial_tau <- tau
  tau[tau > 10 | tau < 0] <- 0.5
  model1 <- NULL
  for (iteration in seq_len(maxiter)) {
    alpha_old <- alpha
    tau_old <- tau
    model1 <- noCovariatesAI(y_working, x, d^2, tau, fixtau, tol)
    tau <- as.numeric(model1$tau)
    alpha <- as.numeric(model1$alpha)
    eta <- as.numeric(model1$eta)
    mu <- family$linkinv(eta)
    derivative <- family$mu.eta(eta)
    d <- derivative / sqrt(family$variance(mu))
    y_working <- eta + (y - mu) / derivative
    y_working[abs(y_working) > 1e3] <- stats::median(y_working)
    change <- 2 * max(abs(alpha - alpha_old) / (abs(alpha) + abs(alpha_old) + tol),
      abs(tau - tau_old) / (abs(tau) + abs(tau_old) + tol))
    if (is.finite(change) && change < tol) break
    if (max(tau) > tol^(-2) || any(!is.finite(c(d, mu, eta)))) {
      tau <- initial_tau
      break
    }
  }
  list(theta = tau, Y = y_working, D = d, converged = iteration < maxiter)
}
