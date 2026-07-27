# SPARKb

`SPARKb` identifies spatially variable features in sparse spatial omics data. It
models feature detection rather than raw counts, tests five Gaussian spatial
scales, and combines their evidence with ACAT.

## Installation

```r
install.packages("remotes")
remotes::install_github("Santanevrgivup/SPARKb")
```

## Quick start

```r
library(SPARKb)

set.seed(1)
location <- cbind(x = runif(100), y = runif(100))
counts <- rbind(
  spatial_feature = rbinom(100, 1, plogis(6 * (location[, 1] - 0.5))),
  null_feature = rbinom(100, 1, 0.3)
)

result <- sparkb(counts, location)
head(result)
attr(result, "bandwidths")
```

`counts` must contain features in rows and locations in columns; `location`
must be an `n_locations` by 2 coordinate matrix. The present release implements
the scalable intercept-only score test, so `covariates` must be `NULL`.

## Method defaults

SPARKb first binarizes counts (`count > 0` by default). It uses five Gaussian
kernels. Their lower scale is adjusted using the median eighth-nearest-neighbour
distance, then the five bandwidths are log-spaced to the original seventh
global-distance scale. Kernel p-values are combined per feature with ACAT, and
the resulting p-values are adjusted using Benjamini-Hochberg FDR.

## Citation

If you use SPARKb, please cite the accompanying manuscript (citation details
will be added upon publication) and the software repository.
