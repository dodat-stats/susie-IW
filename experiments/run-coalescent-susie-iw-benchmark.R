#!/usr/bin/env Rscript

library(data.table)

numeric_grid <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.numeric(default))
  as.numeric(strsplit(value, ",", fixed = TRUE)[[1L]])
}

J_grid <- as.integer(numeric_grid(
  "COAL_BENCH_J", c(500, 1000, 2000)
))
number_of_replications <- as.integer(Sys.getenv("COAL_BENCH_N_REP", "20"))
panel_seeds <- as.integer(numeric_grid(
  "COAL_BENCH_PANEL_SEEDS", seq.int(61L, length.out = number_of_replications)
))
if (length(panel_seeds) != number_of_replications) {
  stop("COAL_BENCH_PANEL_SEEDS must contain COAL_BENCH_N_REP values")
}
cores <- as.integer(Sys.getenv("COAL_BENCH_CORES", "3"))
output_directory <- Sys.getenv(
  "COAL_BENCH_OUTPUT_DIRECTORY",
  file.path("experiments", "coalescent-susie-iw-benchmark-20rep-workers")
)
log_directory <- file.path(output_directory, "logs")
dir.create(log_directory, recursive = TRUE, showWarnings = FALSE)

tasks <- data.table(
  replication = seq_len(number_of_replications),
  J = rep(J_grid, length.out = number_of_replications),
  panel_seed = panel_seeds
)
run_task <- function(index) {
  task <- tasks[index]
  log_file <- file.path(
    log_directory, sprintf("J%d-panel%d.log", task$J, task$panel_seed)
  )
  status <- system2(
    "Rscript",
    c("--vanilla", file.path(
      "analysis", "coalescent-susie-iw-benchmark-worker.R"
    ), task$J, task$panel_seed),
    stdout = log_file, stderr = log_file
  )
  data.table(
    replication = task$replication, J = task$J,
    panel_seed = task$panel_seed,
    status = status, log_file = log_file
  )
}

results <- if (.Platform$OS.type == "unix" && cores > 1L) {
  rbindlist(parallel::mclapply(
    seq_len(nrow(tasks)), run_task, mc.cores = cores,
    mc.preschedule = FALSE
  ))
} else {
  rbindlist(lapply(seq_len(nrow(tasks)), run_task))
}
print(results)
if (any(results$status != 0L)) quit(save = "no", status = 1L)
