# BOSTON

**Binary-based inference fOr large-Scale sparse spaTial Omics data aNalysis**
detects spatially variable features from binary detection patterns in spatial
omics data. It fits an intercept-only binomial null model, tests five Gaussian
spatial kernels, combines their p-values with ACAT, and applies BY adjustment
by default. BOSTON succeeds the independent SPARKb package.

## Installation

```r
install.packages("remotes")
remotes::install_github("Santanevrgivup/BOSTON")
```

Compilation requires an R-compatible C++ toolchain (Rtools on Windows,
Xcode Command Line Tools on macOS). OpenMP is optional when supported by R's
toolchain. R >= 4.1.0, Matrix, Rcpp and RcppArmadillo are required.

## Quick start

```r
library(BOSTON)
set.seed(1)
location <- cbind(x = runif(100), y = runif(100))
counts <- rbind(
  spatial_feature = rbinom(100, 1, plogis(6 * (location[, 1] - 0.5))),
  null_feature = rbinom(100, 1, 0.3)
)
ids <- paste0("spot_", seq_len(nrow(location)))
colnames(counts) <- rownames(location) <- ids

result <- boston(counts, location, seed = 123L)
head(result)
attr(result, "bandwidths")
```

`counts` is a non-negative numeric feature-by-location matrix or a sparse
`Matrix::dgCMatrix`. `location` is a finite numeric location-by-2 matrix or
data frame. When supplied, count column IDs and coordinate row IDs must be
unique and in exactly the same order. Both may be omitted, with a warning
that alignment is assumed. Coordinate axes must use commensurate units;
automatic axis-wise standardization is not performed.

`covariates` is reserved and must be `NULL` in this release. At least 10
observations are required for testing. The default bandwidth pilot additionally
requires more than 16 distinct coordinates.

## Adaptive bandwidths

```r
pilot <- adaptive_bandwidths(location, random_seed = 123L, diagnostics = TRUE)
bandwidths <- pilot$kparam
result <- boston(counts, location, bandwidths = bandwidths)

# For smaller datasets or a coordinate-only sensitivity analysis:
bw_k8 <- adaptive_bandwidths(location, k_local = 8L)
```

For bandwidth estimation only, exact duplicate coordinates are collapsed.
The default lower endpoint is half the median distance to the 16th nearest
distinct coordinate, using at most 1024 query points. The upper endpoint is
half the median distance of 50000 pairs sampled with replacement from distinct
coordinates, excluding self-pairs. If the upper candidate is no larger than
the lower candidate, the upper endpoint becomes twice the lower, with a warning.
Five bandwidths are log-spaced between the endpoints. They are sigma in
`exp(-distance^2 / (2 * sigma^2))`, in the original coordinate units.

The pilot uses a local C++ RNG and does not change R's random-number state.
Repeated calls with the same seed and coordinates are reproducible within a
toolchain; bit-for-bit agreement across C++ standard libraries is not promised.
All original observations remain in the downstream analysis. Diagnostic output
contains the distinct-coordinate indices, local radii, sampled pair distances,
candidate bounds and expansion flag. KNN scales reflect sampling geometry;
they are heuristic reference scales, not estimated biological correlation lengths.
No automatic grid or removal of near-duplicate coordinates is performed.

## Results and multiple testing

`boston()` returns a `boston_result` data frame ordered by combined p-value,
containing `feature`, `pvalue_gaussian_1` through `pvalue_gaussian_5`, `pvalue`,
`detected_locations`, `status`, `message`, and `padj_BY` by default.

```r
result_bh <- boston(counts, location, bandwidths = bandwidths, padj_method = "BH")
# Adjusted column: padj_BH
```

`weights = NULL` gives equal ACAT weights. Supplied weights correspond to the
five bandwidths and should be finite, non-negative, and have a positive sum.
The inherited ACAT implementation is retained unchanged, including its handling
of boundary p-values. All five kernel p-values must be finite before combination.
Constant or failed features retain an NA combined p-value; inspect `status` and
`message`. `not_converged` rows may still contain p-values and require review.
Adjustment is over the non-missing p-values supplied in a call. If you split
features across calls, collect all raw combined p-values and adjust once over
the intended full feature universe.

The bandwidths and detection threshold are stored as result attributes. Save
the result with `saveRDS()` to preserve these attributes, and separately record
the pilot settings, feature universe and package version for reproducibility.

## Large-data execution

The bandwidth pilot avoids a dense all-pairs matrix. Exact Gaussian score tests
still have substantial quadratic-in-location computation and dense intermediate
location-by-feature matrices; this release does not add feature batching.

`num_cores` controls R's feature-wise null-model workers, not OpenMP or BLAS
threads. To avoid nested thread oversubscription, configure the R process before
starting it, for example on Linux/macOS:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript analysis.R
```

Then select `num_cores` within the analysis script according to allocated CPUs.
Windows uses serial null-model fitting. Benchmark memory on a subset before a
large run; if processing feature chunks, reuse the same bandwidth vector and
apply one final adjustment across all chunks.

## Migration and provenance

Use `library(BOSTON)` and `boston()` in place of `library(SPARKb)` and `sparkb()`.
The new automatic bandwidth rule differs from SPARKb 0.1.1, so automatic results
can change. With identical explicit bandwidths and analysis inputs, the inherited
null-model, Gaussian score-test and ACAT calculations are preserved.

BOSTON retains GPL-3 licensing and the SpatialX-derived statistical core.
Authors: Yuhan Zhao and Shiquan Sun. Cite the software repository; manuscript
citation details will be added when available.
