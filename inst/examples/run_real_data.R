# Run with: Rscript run_real_data.R INPUT.rds OUTPUT_DIR [options]
# Input: list(counts = feature-by-location matrix, location = location-by-2 matrix,
#             covariates = optional location-by-predictor matrix/data.frame)
main <- function(args = commandArgs(TRUE)) {
  if (length(args) < 2L) stop(paste(
    "Usage: Rscript run_real_data.R INPUT.rds OUTPUT_DIR",
    "[--legacy] [--depth] [--max-features=20] [--max-locations=0] [--cores=1]"))
  input <- normalizePath(args[1], mustWork = TRUE)
  output <- args[2]
  option_integer <- function(name, default) {
    found <- args[startsWith(args, paste0(name, "="))]
    if (!length(found)) return(default)
    if (length(found) != 1L) stop("Repeated option: ", name)
    value <- suppressWarnings(as.numeric(sub(paste0("^", name, "="), "", found)))
    if (length(value) != 1L || !is.finite(value) || value < 0 ||
        value != floor(value) || value > .Machine$integer.max) stop("Invalid option: ", name)
    as.integer(value)
  }
  max_features <- option_integer("--max-features", 20L)
  max_locations <- option_integer("--max-locations", 0L)
  cores <- option_integer("--cores", 1L)
  if (cores < 1L) stop("--cores must be positive.")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  if (length(list.files(output, all.files = TRUE, no.. = TRUE))) {
    stop("Use an empty output directory to preserve previous analyses.")
  }
  if ("--legacy" %in% args) {
    if (!requireNamespace("SpatialX", quietly = TRUE)) {
      stop("Reading a legacy S4 RDS requires SpatialX. Export a plain list first.")
    }
    legacy <- readRDS(input)
    assay <- legacy@assays$Spatial
    data <- list(counts = assay@counts, location = assay@location)
    rm(legacy, assay)
  } else data <- readRDS(input)
  if (!is.list(data) || is.null(data$counts) || is.null(data$location)) {
    stop("Input must contain counts and location.")
  }
  counts <- data$counts
  location <- as.matrix(data$location)
  covariates <- data$covariates
  ids <- colnames(counts)
  if (is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("Tutorial inputs require unique count-column location IDs.")
  }
  if (is.null(rownames(location)) || anyDuplicated(rownames(location)) ||
      !setequal(ids, rownames(location))) stop("Coordinate IDs do not match counts.")
  location <- location[ids, , drop = FALSE]
  if (ncol(location) != 2L || any(!is.finite(location))) stop("Use finite 2D coordinates.")
  if (!is.null(covariates)) {
    if (is.null(rownames(covariates)) || anyDuplicated(rownames(covariates)) ||
        !setequal(ids, rownames(covariates))) stop("Covariate IDs do not match counts.")
    covariates <- covariates[ids, , drop = FALSE]
  }

  # Derive QC/depth from the complete input feature matrix, before subsetting.
  totals <- as.numeric(Matrix::colSums(counts))
  detections <- as.numeric(Matrix::colSums(counts > 0))
  if ("--depth" %in% args) {
    if (!is.null(covariates)) stop("Supply input covariates OR --depth, not both.")
    log_depth <- log1p(totals)
    if (!is.finite(stats::sd(log_depth)) || stats::sd(log_depth) == 0) {
      stop("Depth is constant; an adjusted model cannot include it.")
    }
    covariates <- data.frame(log_depth = as.numeric(scale(log_depth)), row.names = ids)
  }
  selected_locations <- seq_along(ids)
  if (max_locations > 0L && max_locations < length(ids)) {
    set.seed(123L)
    selected_locations <- sort(sample(selected_locations, max_locations))
  }
  counts <- counts[, selected_locations, drop = FALSE]
  location <- location[selected_locations, , drop = FALSE]
  if (!is.null(covariates)) covariates <- covariates[selected_locations, , drop = FALSE]
  detected <- as.numeric(Matrix::rowSums(counts > 0))
  # A prespecified demonstration filter, independent of spatial p-values.
  selected <- which(detected >= 10L & detected < ncol(counts))
  if (max_features > 0L) selected <- head(selected, max_features)
  if (!length(selected)) stop("No feature passed the detection filter.")
  if (is.null(rownames(counts))) rownames(counts) <- paste0("feature_", seq_len(nrow(counts)))
  tested_features <- rownames(counts)[selected]
  counts <- counts[selected, , drop = FALSE]
  qc <- data.frame(location_id = ids[selected_locations],
    total_counts = totals[selected_locations], detected_features = detections[selected_locations])
  utils::write.csv(qc, file.path(output, "location_qc.csv"), row.names = FALSE)
  utils::write.csv(data.frame(feature = tested_features, detected_locations = detected[selected]),
    file.path(output, "feature_qc.csv"), row.names = FALSE)
  pilot <- BOSTON::adaptive_bandwidths(location, random_seed = 123L, diagnostics = TRUE)
  elapsed <- system.time(result <- BOSTON::boston(counts, location,
    covariates = covariates, bandwidths = pilot$kparam, num_cores = cores, verbose = TRUE))
  utils::write.csv(as.data.frame(result), file.path(output, "results.csv"), row.names = FALSE)
  utils::write.csv(as.data.frame(table(result$status)), file.path(output, "model_status.csv"),
    row.names = FALSE)
  saveRDS(list(results = result, bandwidth_pilot = pilot,
    metadata = list(input_file = input, input_md5 = unname(tools::md5sum(input)),
      package_version = as.character(utils::packageVersion("BOSTON")),
      feature_filter = "detected in >=10 and <n selected locations",
      tested_features = tested_features, selected_location_ids = rownames(location),
      covariates = covariates, adjustment_universe = tested_features,
      elapsed_seconds = unname(elapsed[["elapsed"]]),
      thread_environment = Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS")),
      session_info = utils::sessionInfo())), file.path(output, "analysis.rds"))
  capture.output(utils::sessionInfo(), file = file.path(output, "sessionInfo.txt"))

  grDevices::pdf(file.path(output, "spatial_detection.pdf"), width = 9, height = 4)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(1, 2))
  graphics::plot(location, pch = 16, cex = .45,
    col = grDevices::hcl.colors(100, "YlOrRd")[cut(log1p(qc$total_counts), 100,
      labels = FALSE)], main = "Total counts (log scale)", asp = 1)
  feature <- result$feature[which(is.finite(result$pvalue))[1L]]
  if (!length(feature) || is.na(feature)) feature <- tested_features[1]
  binary <- as.numeric(counts[feature, ] > 0)
  graphics::plot(location, pch = 16, cex = .5,
    col = ifelse(binary > 0, "#C53929", "#D9D9D9"),
    main = paste("Detection:", feature), asp = 1)
  message("Saved analysis to ", normalizePath(output))
  invisible(result)
}
main()
