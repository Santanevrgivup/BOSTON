# BOSTON 0.3.0

* Support full-rank numeric/factor covariates with an automatic intercept.
* Add low-rank covariate null updates and blocked Gaussian score moments.
* Accept named dense PSD kernels (up to 2000 locations), with alignment,
  symmetry, PSD and kernel-weight checks. Degenerate projected kernels fail
  explicitly. The ACAT function is unchanged.
* Include an executable real-data QC, adjustment, plotting and saving tutorial.
* Guard the OpenMP header for toolchains that compile without OpenMP.

# BOSTON 0.2.0

* Rename the package and primary function to BOSTON and `boston()`.
* Integrate the coordinate-only kNN bandwidth pilot: K=16, up to 1024 queries,
  50000 sampled distinct-coordinate pairs, median endpoints divided by two,
  and five log-spaced Gaussian bandwidths.
* Preserve all original observations for testing, and preserve R's RNG state
  during bandwidth sampling. Expose pilot diagnostics and tuning parameters.
* Document the existing BY adjustment default and kernel weights. Retain the
  original ACAT implementation and the fixed-bandwidth score-test calculations.

# SPARKb 0.1.1 (predecessor package)

* Validate unique spatial location IDs and require exact count-column to
  coordinate-row alignment when IDs are supplied.
* Make duplicate feature IDs unique with `make.unique()` while preserving
  biologically meaningful punctuation in feature names.
* Report null-model non-convergence and feature-level fitting errors in results.

# SPARKb 0.1.0 (predecessor package)

* First independent package release of the binary spatial feature test.
* Uses five adaptive Gaussian kernels and ACAT p-value combination.
* Provides input validation, unit tests, documentation, and GitHub Actions checks.
