#!/usr/bin/env Rscript

# Reproduce the coalescent experiments used in the manuscript. Both worker
# pipelines checkpoint at the dataset level, so this driver can safely be
# restarted after interruption.

run_script <- function(script) {
  message("Running ", script)
  status <- system2("Rscript", c("--vanilla", script))
  if (!identical(status, 0L)) stop("Failed: ", script)
}

scripts <- c(
  file.path("experiments", "run-eta0-calibration-coalescent-L.R"),
  file.path("experiments", "summarize-eta0-calibration-coalescent-L.R"),
  file.path("experiments", "run-coalescent-susie-iw-benchmark.R"),
  file.path("experiments", "summarize-coalescent-susie-iw-benchmark.R")
)
invisible(lapply(scripts, run_script))
