#!/usr/bin/env Rscript

source(file.path("experiments", "ancestry-mismatch-coalescent-functions.R"))
dir.create(ancestry_worker_directory, recursive = TRUE, showWarnings = FALSE)

expected_methods <- c(
  "SuSiE-RSS: in-sample R", "SuSiE-RSS: reference R0",
  "collapsed-R", "AIW-N"
)
valid_worker <- function(T_generations, panel_seed) {
  filename <- file.path(
    ancestry_worker_directory,
    sprintf("results-T%d-seed%d.csv", T_generations, panel_seed)
  )
  if (!file.exists(filename)) return(FALSE)
  result <- tryCatch(utils::read.csv(filename), error = function(e) NULL)
  if (is.null(result)) return(FALSE)
  nrow(result) == length(ancestry_N0_grid) * length(expected_methods) &&
    identical(sort(unique(result$N0)), sort(ancestry_N0_grid)) &&
    identical(sort(unique(result$method)), sort(expected_methods)) &&
    all(result$T_generations == T_generations) &&
    all(result$panel_seed == panel_seed) && all(result$J == ancestry_J)
}

tasks <- expand.grid(
  T_generations = ancestry_T_grid,
  panel_seed = ancestry_panel_seeds,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
tasks$complete <- mapply(
  valid_worker, tasks$T_generations, tasks$panel_seed
)
remaining <- tasks[!tasks$complete, , drop = FALSE]
message(sprintf(
  "%d/%d workers already complete; %d remain",
  sum(tasks$complete), nrow(tasks), nrow(remaining)
))

if (nrow(remaining)) {
  number_of_cores <- ancestry_get_integer(
    "ANCESTRY_CORES", min(3L, parallel::detectCores(logical = FALSE))
  )
  number_of_cores <- max(1L, min(number_of_cores, nrow(remaining)))
  run_one <- function(index) {
    task <- remaining[index, ]
    log_file <- file.path(
      ancestry_worker_directory,
      sprintf("worker-T%d-seed%d.log", task$T_generations, task$panel_seed)
    )
    status <- system2(
      file.path(R.home("bin"), "Rscript"),
      c(
        "--vanilla", "experiments/ancestry-mismatch-coalescent-worker.R",
        as.character(task$T_generations), as.character(task$panel_seed)
      ),
      stdout = log_file, stderr = log_file
    )
    data.frame(
      T_generations = task$T_generations,
      panel_seed = task$panel_seed,
      status = status,
      log_file = log_file,
      stringsAsFactors = FALSE
    )
  }
  completed <- parallel::mclapply(
    seq_len(nrow(remaining)), run_one, mc.cores = number_of_cores,
    mc.preschedule = FALSE
  )
  completed <- do.call(rbind, completed)
  failed <- completed$status != 0L | !mapply(
    valid_worker, completed$T_generations, completed$panel_seed
  )
  if (any(failed)) {
    print(completed[failed, ], row.names = FALSE)
    stop("One or more ancestry-mismatch workers failed; inspect their logs")
  }
}

message("All ancestry-mismatch workers are complete")
