.validate_covariates <- function(covariates, counts) {
  if (is.null(covariates)) return(NULL)
  if ((!is.matrix(covariates) && !is.data.frame(covariates)) ||
      nrow(covariates) != ncol(counts) || ncol(covariates) < 1L) {
    stop("`covariates` must be a matrix or data.frame with one row per location and at least one column.",
      call. = FALSE)
  }
  ids <- rownames(covariates)
  if (is.data.frame(covariates) && .row_names_info(covariates, 1L) < 0) ids <- NULL
  .validate_external_ids(ids, colnames(counts), "covariate rows")
  if (is.matrix(covariates)) {
    if (!is.numeric(covariates) || any(!is.finite(covariates))) {
      stop("Matrix `covariates` must contain finite numeric values.", call. = FALSE)
    }
    x <- cbind(`(Intercept)` = 1, covariates)
  } else {
    valid <- vapply(covariates, function(z) is.numeric(z) || is.factor(z) || is.logical(z),
      logical(1))
    if (!all(valid) || anyNA(covariates)) {
      stop("Covariate data.frame columns must be numeric, factor or logical, without missing values.",
        call. = FALSE)
    }
    x <- stats::model.matrix(~ ., data = covariates, na.action = stats::na.fail)
  }
  if (any(!is.finite(x)) || ncol(x) >= nrow(x) || qr(x)$rank != ncol(x)) {
    stop("Covariates plus the automatic intercept must be finite and full rank, with fewer columns than locations. Remove constant, redundant or intercept columns.",
      call. = FALSE)
  }
  storage.mode(x) <- "double"
  x
}

.validate_external_ids <- function(ids, expected, label) {
  if (is.null(expected) && is.null(ids)) return(invisible(NULL))
  if (is.null(expected) || is.null(ids) || anyNA(ids) || any(!nzchar(ids)) ||
      anyDuplicated(ids) || !identical(ids, expected)) {
    stop("The IDs of ", label,
      " must be unique and match `colnames(counts)` in exactly the same order; supply IDs on both inputs or neither.",
      call. = FALSE)
  }
  invisible(NULL)
}

.validate_kernels <- function(kernels, counts) {
  if (is.null(kernels)) return(NULL)
  if (is.matrix(kernels)) kernels <- list(custom_1 = kernels)
  if (!is.list(kernels) || !length(kernels)) {
    stop("`kernels` must be a numeric matrix or a non-empty named list of matrices.", call. = FALSE)
  }
  if (is.null(names(kernels))) names(kernels) <- paste0("custom_", seq_along(kernels))
  if (anyNA(names(kernels)) || anyDuplicated(names(kernels)) ||
      any(!grepl("^[A-Za-z][A-Za-z0-9_]*$", names(kernels)))) {
    stop("Kernel names must be unique and use letters, digits or underscores, starting with a letter.",
      call. = FALSE)
  }
  n <- ncol(counts)
  if (n > 2000L) {
    stop("Dense custom kernels are currently limited to 2000 locations because exact PSD validation requires an eigendecomposition; use built-in Gaussian kernels for larger inputs.",
      call. = FALSE)
  }
  for (j in seq_along(kernels)) {
    k <- kernels[[j]]
    if (!is.matrix(k) || !is.numeric(k) || !identical(dim(k), c(n, n)) ||
        any(!is.finite(k))) {
      stop("Kernel '", names(kernels)[j], "' must be a finite numeric location-by-location matrix.",
        call. = FALSE)
    }
    .validate_external_ids(rownames(k), colnames(counts), "kernel rows")
    .validate_external_ids(colnames(k), colnames(counts), "kernel columns")
    tolerance <- 1e-8 * max(1, max(abs(k)))
    if (max(abs(k - t(k))) > tolerance) {
      stop("Kernel '", names(kernels)[j], "' must be symmetric.", call. = FALSE)
    }
    k <- (k + t(k)) / 2
    values <- eigen(k, symmetric = TRUE, only.values = TRUE)$values
    if (min(values) < -1e-8 * max(1, max(abs(values))) || max(values) <= 0) {
      stop("Kernel '", names(kernels)[j], "' must be nonzero and positive semidefinite.", call. = FALSE)
    }
    storage.mode(k) <- "double"
    kernels[[j]] <- k
  }
  kernels
}
