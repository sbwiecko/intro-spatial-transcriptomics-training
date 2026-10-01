# =========================================================================
# data_prep_helpers.R
# 
# This script provides utility functions for the Spatial Transcriptomics 
# data ingestion pipeline. It handles robust file downloading (with resume 
# capabilities), archive extraction, path resolution, and basic spatial 
# subsetting.
# =========================================================================

# -------------------------------------------------------------------------
# Path Resolution & Environment Setup
# -------------------------------------------------------------------------

# Increase the default download timeout from 60s to 1 hour (3600s) to 
# prevent massive 10x Genomics dataset downloads from failing prematurely.
options(timeout = 3600)

# Finds the project root directory by searching iteratively upwards for .Rprofile
find_project_root <- function(path = getwd()) {
  path <- normalizePath(path, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, ".Rprofile"))) {
      return(path)
    }
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not find project root containing .Rprofile")
    }
    path <- parent
  }
}

# Initialize core paths globally
project_root <- find_project_root()

# Convenience function to quickly generate paths relative to the data directory
data_path <- function(...) file.path(project_root, "data", ...)

# Ensure required data directories exist
dir.create(data_path("raw"), showWarnings = FALSE, recursive = TRUE)
dir.create(data_path("extracted"), showWarnings = FALSE, recursive = TRUE)
dir.create(data_path("prep_reports"), showWarnings = FALSE, recursive = TRUE)
dir.create(data_path("Human_Colon_Cancer_P2"), showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------------------
# Safe Downloading Functions
# -------------------------------------------------------------------------

# Basic download function that skips the download if the file already exists
download_if_missing <- function(url, destfile) {
  dir.create(dirname(destfile), showWarnings = FALSE, recursive = TRUE)
  size <- suppressWarnings(file.info(destfile)$size)
  if (!file.exists(destfile) || is.na(size) || size == 0) {
    message("Downloading: ", basename(destfile))
    download.file(url, destfile = destfile, mode = "wb", method = "curl")
  } else {
    message("Already present: ", destfile)
  }
}

# Advanced resilient download function capable of resuming broken downloads 
# using curl's --continue-at argument and checking expected file sizes
download_if_missing_or_incomplete <- function(url, destfile, expected_size = NA_real_) {
  dir.create(dirname(destfile), showWarnings = FALSE, recursive = TRUE)
  current_size <- if (file.exists(destfile)) file.info(destfile)$size else 0

  # Skip if perfectly matched
  if (!is.na(expected_size) && current_size == expected_size) {
    message("Already present with expected size: ", destfile)
    return(invisible(destfile))
  }

  # Skip if file exists and we don't have an expected size constraint
  if (is.na(expected_size) && current_size > 0) {
    message("Already present: ", destfile)
    return(invisible(destfile))
  }

  message("Downloading/resuming: ", basename(destfile))
  status <- system2(
    "curl",
    args = c(
      "--location",
      "--fail",
      "--continue-at", "-",
      "--retry", "10",
      "--retry-delay", "20",
      "--output", destfile,
      url
    )
  )

  if (!identical(status, 0L)) {
    stop("Download failed for: ", url)
  }

  if (!is.na(expected_size)) {
    final_size <- file.info(destfile)$size
    if (!identical(as.numeric(final_size), as.numeric(expected_size))) {
      stop(
        "Downloaded file has unexpected size: ", destfile,
        "\nExpected: ", expected_size,
        "\nObserved: ", final_size,
        "\nRerun this chunk to resume."
      )
    }
  }

  invisible(destfile)
}

# -------------------------------------------------------------------------
# Extraction Utilities
# -------------------------------------------------------------------------

# Safely untars an archive only if a specific marker file/directory is missing
extract_tar_if_missing <- function(tarfile, exdir, marker) {
  if (!file.exists(marker) && !dir.exists(marker)) {
    message("Extracting: ", basename(tarfile))
    dir.create(exdir, showWarnings = FALSE, recursive = TRUE)
    untar(tarfile, exdir = exdir)
  } else {
    message("Already extracted: ", marker)
  }
}

# Safely unzips an archive only if a specific marker file/directory is missing
unzip_if_missing <- function(zipfile, exdir, marker) {
  if (!file.exists(marker) && !dir.exists(marker)) {
    message("Unzipping: ", basename(zipfile))
    dir.create(exdir, showWarnings = FALSE, recursive = TRUE)
    unzip(zipfile, exdir = exdir)
  } else {
    message("Already unzipped: ", marker)
  }
}

# -------------------------------------------------------------------------
# Spatial Utilities
# -------------------------------------------------------------------------

# Recursively hunts for the experiment.xenium file to locate the correct 
# Xenium output folder
find_xenium_outs <- function(root) {
  candidates <- c(
    file.path(root, "outs"),
    list.dirs(root, recursive = TRUE, full.names = TRUE)
  )
  candidates <- unique(candidates)
  hits <- candidates[file.exists(file.path(candidates, "experiment.xenium"))]
  if (length(hits) == 0) {
    stop(
      "Could not find a Xenium outs directory containing experiment.xenium under: ", 
      root)
  }
  hits[[1]]
}

# Helper to easily subset a SpatialExperiment object based on a coordinate bounding box
subset_by_roi <- function(x, roi) {
  xy <- spatialCoords(x)
  keep <- xy[, 1] >= roi[["xmin"]] & xy[, 1] <= roi[["xmax"]] &
    xy[, 2] >= roi[["ymin"]] & xy[, 2] <= roi[["ymax"]]
  x[, keep]
}
