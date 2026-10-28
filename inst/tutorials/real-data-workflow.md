# A reproducible real-data workflow

This tutorial covers input preparation, location alignment, detection QC,
covariate adjustment, model diagnostics, BY adjustment, plotting and saving.
The script is installed with BOSTON:

```r
script <- system.file("examples", "run_real_data.R", package = "BOSTON")
```

## Prepare independent inputs

Save a plain R list containing the complete count matrix and coordinates:

```r
# counts: features x locations; location: locations x 2.
# IDs must be unique. The script explicitly reorders matching coordinate IDs.
saveRDS(list(counts = counts, location = location), "my_spatial_data.rds")
```

For a legacy SpatialX object, one-time export requires its class definitions:

```r
library(SpatialX)
spe <- readRDS("legacy_spatialx.rds")
assay <- spe@assays$Spatial
saveRDS(list(counts = assay@counts, location = assay@location), "my_spatial_data.rds")
```

SpatialX is not a BOSTON dependency. Keep separate slices in separate analyses
unless you have explicitly constructed a suitable block kernel. Overlapping
local coordinates from different slices do not define shared physical space.
Use the same units for both coordinate axes; independent x/y standardization
changes the spatial metric. Inputs are not silently normalized.

## Run a small, explicitly labelled demonstration

In a shell, replace SCRIPT with the installed script path printed above:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript SCRIPT \
  my_spatial_data.rds demo_output --max-features=20 --cores=2
```

The script chooses the first 20 input features detected in at least 10 and fewer
than all locations. This is a demonstration, not a full discovery analysis.
The BY universe is precisely these 20 features. Set `--max-features=0` to test
every feature passing the filter, after benchmarking available memory. Use a new,
empty output directory for every run. `--max-locations=500` optionally draws a
reproducible subset for a small illustration; this changes the scientific sample.

For the author's Slide-tags human melanoma ATAC legacy file, the tested command
is equivalent to:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript SCRIPT \
  LEGACY_SLIDETAGS.rds slidetags_demo --legacy --depth --max-features=2 --cores=2
```

The legacy input has 44,957 features and 2,535 locations. This repository does
not redistribute that processed object. Obtain source data through the original
study and its access terms; the public study entry is
https://singlecell.broadinstitute.org/single_cell/study/SCP2171/slide-tags-snrna-seq-on-human-melanoma .
That link describes the related RNA study, not a download of the processed ATAC
RDS used here. Users may run the same workflow on their own legitimately obtained
count and coordinate matrices. No institutional filesystem is required by the script.

## Adjust for covariates when the scientific question calls for it

`--depth` calculates standardized log1p library size from the COMPLETE input
feature matrix before demonstration subsetting. This is an illustrative predictor,
not a universal remedy for technical confounding. Counts-derived depth includes
the tested features; an independently measured QC variable or predefined robust
depth estimate may be preferable. Consider the assay and study design.

Alternatively, put named covariates into the input list:

```r
covariates <- data.frame(log_depth = log1p(library_size),
                        batch = factor(batch), row.names = colnames(counts))
saveRDS(list(counts = counts, location = location, covariates = covariates),
        "adjusted_inputs.rds")
result <- BOSTON::boston(counts, location, covariates = covariates)
```

An intercept is added automatically. Avoid constant predictors, duplicate columns,
an explicit intercept, missing values, or a design with as many columns as locations.
Specify character categories as factors and preserve their reference levels.
Covariate-adjusted results test residual spatial structure conditional on the
design. Adjusting for spatial domains or cell-type composition can remove the
biological signal you intended to detect. Choose predictors before inspecting
favourable p-values, and document the unadjusted/adjusted scientific questions.

## Inspect and save results

The output directory contains:

| File | Purpose |
| --- | --- |
| `location_qc.csv` | Library size and number of detected input features |
| `feature_qc.csv` | Detection counts for the selected feature universe |
| `results.csv` | Per-kernel/combined p-values, BY values and statuses |
| `model_status.csv` | Counts of model states |
| `analysis.rds` | Result attributes, full pilot diagnostics, input checksum, exact features/locations, covariates, runtime and session information |
| `sessionInfo.txt` | Human-readable environment record |
| `spatial_detection.pdf` | QC depth map and a selected feature's detection map |

Inspect `fit_failed`, `test_failed` and `not_converged` rows before interpretation.
The inherited intercept-only convergence flag is not a complete numerical diagnostic;
an `ok` label is not a statistical calibration guarantee. Significant results require
biological validation and assessment of residual technical structure.

BY is computed once over non-missing combined p-values in the call. If manually
chunking features, reuse the same locations, covariates and bandwidths, collect
all raw combined p-values, and run `p.adjust(all_pvalues, "BY")` once. Report both
the intended feature universe and the number of failed/missing tests.

## Custom positive-semidefinite kernels

For a small dataset, an explicit Gaussian matrix is a useful reference:

```r
bw <- BOSTON::adaptive_bandwidths(location)
d2 <- as.matrix(dist(location))^2
k <- exp(-d2 / (2 * bw[3]^2))
dimnames(k) <- list(colnames(counts), colnames(counts))
custom_result <- BOSTON::boston(counts, location, covariates = covariates,
                               kernels = list(gaussian_reference = k))
```

Custom kernels replace, rather than supplement, the default five Gaussian
kernels. `weights` must have one value per supplied kernel. Each matrix is checked
for dimensions, matching IDs on both axes, finite entries, symmetry and PSD.
The exact eigendecomposition makes this path limited to 2,000 locations.
The 2,535-location Slide-tags example therefore uses built-in Gaussian kernels;
do not materialize its matrices merely to follow this custom-kernel illustration.
A constant kernel is removed by the intercept projection and yields `test_failed`.
A kernel constructed from the tested feature's own measurements requires separate
inference validation; PSD alone does not make such a choice statistically valid.

## Evidence and scope

The implementation tests compare blocked score moments against explicit dense
P-matrix algebra, and the covariate null update against the inherited dense
covariate implementation. Real-data runs establish execution and numerical
agreement, not type-I-error control or universal superiority of a covariate model.
Non-spatial variance identifiability and PQL/Satterthwaite calibration remain
methodological questions requiring separate theoretical and simulation evidence.
