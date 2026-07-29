# SPARKb 0.1.1

* Validate unique spatial location IDs and require exact count-column to
  coordinate-row alignment when IDs are supplied.
* Make duplicate feature IDs unique with `make.unique()` while preserving
  biologically meaningful punctuation in feature names.
* Report null-model non-convergence and feature-level fitting errors in results.

# SPARKb 0.1.0

* First independent package release of the binary spatial feature test.
* Uses five adaptive Gaussian kernels and ACAT p-value combination.
* Provides input validation, unit tests, documentation, and GitHub Actions checks.
