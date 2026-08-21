#!/usr/bin/env Rscript

library(data.table)

numeric_grid <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.numeric(default))
  as.numeric(strsplit(value, ",", fixed = TRUE)[[1L]])
}

J_grid <- as.integer(numeric_grid(
  "ETA0_L_J_GRID", c(100, 250, 500, 1000, 2000)
))
panel_seeds <- as.integer(numeric_grid(
  "ETA0_L_PANEL_SEEDS", c(51, 52)
))
cores <- as.integer(Sys.getenv("ETA0_L_CORES", "3"))
output_directory <- Sys.getenv(
  "ETA0_L_OUTPUT_DIRECTORY",
  file.path("experiments", "eta0-calibration-coalescent-L-workers")
)
log_directory <- file.path(output_directory, "logs")
dir.create(log_directory, recursive = TRUE, showWarnings = FALSE)

tasks <- CJ(J = J_grid, panel_seed = panel_seeds)
run_task <- function(index) {
  task <- tasks[index]
  log_file <- file.path(
    log_directory,
    sprintf("J%d-panel%d.log", task$J, task$panel_seed)
  )
  status <- system2(
    "Rscript",
    c(
      "--vanilla",
      file.path("experiments", "eta0-calibration-coalescent-L-worker.R"),
      task$J, task$panel_seed
    ),
    stdout = log_file, stderr = log_file
  )
  data.table(
    J = task$J, panel_seed = task$panel_seed,
    status = status, log_file = log_file
  )
}

task_results <- if (.Platform$OS.type == "unix" && cores > 1L) {
  rbindlist(parallel::mclapply(
    seq_len(nrow(tasks)), run_task,
    mc.cores = cores, mc.preschedule = FALSE
  ))
} else {
  rbindlist(lapply(seq_len(nrow(tasks)), run_task))
}
print(task_results)
if (any(task_results$status != 0L)) quit(save = "no", status = 1L)
