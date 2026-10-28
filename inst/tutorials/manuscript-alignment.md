# Manuscript / implementation correspondence (BOSTON 0.3.0)

| Manuscript element | Implemented route | Boundary to state |
| --- | --- | --- |
| Binary detection indicators | `counts > detection_threshold`, default zero | Count magnitude is discarded after detection |
| General fixed-effect design X | `covariates=` numeric matrix or numeric/factor/logical data.frame; automatic intercept | Full rank, finite predictors and exact location alignment required; no formula or offset interface |
| PQL nuisance fitting | Intercept path inherited; covariate path uses low-rank weighted projection / AI updates | Numerical tests do not establish identifiability or inferential calibration |
| Shared null fit across kernels | One fitted null model per nonconstant feature | Refit if predictors or the observation set change |
| Five Gaussian kernels | Built-in coordinate-only kNN pilot and block kernel calculations | Two coordinate dimensions, common coordinate units; no nonstationary covariance |
| Alternative PSD kernels | `kernels=` matrix or named list | Dense interface limited to 2000 locations; PSD and matching IDs required; data-dependent selection needs separate inference validation |
| Score statistic | S = working_y' P K P working_y / 2 | Low-rank/block calculations checked against explicit dense P algebra |
| Satterthwaite moments | e = tr(PK)/2; v = tr(PKPK)/2 | Working-model moment approximation; not a claim of exact finite-sample calibration |
| Multi-kernel combination | Inherited `.acat()` function, unchanged | All supplied kernel p-values must be finite; inherited boundary handling retained |
| BY multiple testing | Default `padj_method="BY"`, output `padj_BY` | One full intended feature universe; do not combine already adjusted chunk values |

Suggested software-scope wording:

> BOSTON accepts two-dimensional spatial coordinates and optional observation-level
> covariates supplied as a full-rank numeric matrix or data frame, with an
> automatically included intercept. A feature-specific null model is fitted
> once and shared across the spatial kernels. The default analysis evaluates
> five stationary Gaussian kernels with coordinate-derived bandwidths. Users
> may alternatively supply a named set of symmetric positive-semidefinite
> kernel matrices for analyses with at most 2,000 observations. Gaussian
> covariate-adjusted score calculations use row-block kernel operations and
> low-rank projection updates. Kernel-specific p-values are combined using
> the inherited ACAT implementation and adjusted using the BY procedure.

This is implementation scope, not evidence of type-I-error control. Report
the exact released version/commit and separate the software's numerical
verification from statistical calibration experiments. The original notebook
was not overwritten by this package update; authors should incorporate the
scope wording into the manuscript where appropriate.
