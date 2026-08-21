#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) {
  stop("Usage: ancestry-mismatch-coalescent-worker.R T_generations panel_seed")
}
T_generations <- as.integer(arguments[1L])
panel_seed <- as.integer(arguments[2L])

source(file.path("experiments", "ancestry-mismatch-coalescent-functions.R"))
source(file.path("code", "01-collapsed-R.R"))
source(file.path("code", "02-AIW-and-AIW-N.R"))
source(file.path("code", "03-SER-fallback-and-plot.R"))

if (!T_generations %in% ancestry_T_grid) stop("T is not in ancestry_T_grid")
if (!panel_seed %in% ancestry_panel_seeds) {
  stop("panel_seed is not in ancestry_panel_seeds")
}
dir.create(ancestry_worker_directory, recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(
  ancestry_worker_directory,
  sprintf("results-T%d-seed%d.csv", T_generations, panel_seed)
)

evaluate_ancestry_fit <- function(object, causal) {
  sets <- susie_get_reported_cs(object)
  covered <- vapply(causal, function(j) {
    any(vapply(sets, function(set) j %in% set, logical(1L)))
  }, logical(1L))
  true_sets <- sum(vapply(sets, function(set) {
    any(set %in% causal)
  }, logical(1L)))
  raw_fit <- iw_unwrap_fit(object)
  data.frame(
    selected_lambda = if (!is.null(object$selected_lambda)) {
      object$selected_lambda
    } else NA_real_,
    raw_eta0 = if (!is.null(object$raw_eta0)) {
      object$raw_eta0
    } else if (!is.null(object$selected_raw_eta0)) {
      object$selected_raw_eta0
    } else NA_real_,
    eta0_used = if (!is.null(object$eta0_used)) {
      object$eta0_used
    } else if (!is.null(object$selected_eta0)) {
      object$selected_eta0
    } else NA_real_,
    selected_model_L = if (!is.null(object$selected_L)) {
      object$selected_L
    } else if (!is.null(object$L)) {
      object$L
    } else if (!is.null(raw_fit$V)) {
      length(raw_fit$V)
    } else NA_integer_,
    reported_L = length(sets),
    true_sets = true_sets,
    false_sets = length(sets) - true_sets,
    causal_power = mean(covered),
    all_causal_covered = all(covered),
    mean_causal_pip = mean(iw_pip(object)[causal]),
    fallback = isTRUE(object$postprocess$fallback_used),
    converged = isTRUE(object$converged),
    stringsAsFactors = FALSE
  )
}

message(sprintf(
  "Simulating two populations: T=%d generations, panel seed=%d",
  T_generations, panel_seed
))
panel <- simulate_two_population_ld(
  J = ancestry_J, T_generations = T_generations,
  panel_seed = panel_seed
)
simulation <- simulate_ancestry_summary_statistics(
  panel$R, T_generations = T_generations, panel_seed = panel_seed
)
v <- simulation$v
R <- panel$R

in_sample_time <- system.time({
  in_sample_fit <- fit_susie_rss(
    v = v, R = R, N = ancestry_N, L = ancestry_maximum_L
  )
})[["elapsed"]]
in_sample_summary <- evaluate_ancestry_fit(in_sample_fit, simulation$causal)

rows <- vector("list", length(ancestry_N0_grid) * 4L)
row_index <- 0L
for (N0 in ancestry_N0_grid) {
  message(sprintf("  N0=%d", N0))
  R0 <- panel$R0[[as.character(N0)]]
  mismatch <- ancestry_ld_metrics(R, R0)

  reference_time <- system.time({
    reference_fit <- fit_susie_rss(
      v = v, R = R0, N = ancestry_N, L = ancestry_maximum_L
    )
  })[["elapsed"]]
  collapsed_time <- system.time({
    collapsed_fit <- susie_iw(
      v = v, R0 = R0, N = ancestry_N, L = ancestry_maximum_L,
      lambda_grid = ancestry_lambda_grid,
      eta0_multiplier = ancestry_collapsed_multiplier,
      scaled_prior_variance = 0.2, estimate_prior_variance = TRUE,
      coverage = ancestry_coverage, min_abs_corr = ancestry_purity,
      ser_fallback = TRUE, partition_ser = FALSE,
      max_iter = 200L, tol = 1e-4, verbose = FALSE
    )
  })[["elapsed"]]
  aiwn_time <- system.time({
    aiwn_fit <- suppressMessages(susie_aiw(
      v = v, R0 = R0, N = ancestry_N, L = ancestry_maximum_L,
      approximation = "N", lambda_grid = ancestry_lambda_grid,
      eta0_multiplier = ancestry_aiwn_multiplier,
      scaled_prior_variance = 0.2, estimate_prior_variance = TRUE,
      coverage = ancestry_coverage, min_abs_corr = ancestry_purity,
      ser_fallback = TRUE, partition_ser = FALSE,
      max_iter = 500L, tol = 1e-3, verbose = FALSE
    ))
  })[["elapsed"]]

  fits <- list(
    "SuSiE-RSS: in-sample R" = in_sample_fit,
    "SuSiE-RSS: reference R0" = reference_fit,
    "collapsed-R" = collapsed_fit,
    "AIW-N" = aiwn_fit
  )
  runtimes <- c(in_sample_time, reference_time, collapsed_time, aiwn_time)
  for (method_index in seq_along(fits)) {
    row_index <- row_index + 1L
    fit_summary <- if (method_index == 1L) {
      in_sample_summary
    } else {
      evaluate_ancestry_fit(fits[[method_index]], simulation$causal)
    }
    rows[[row_index]] <- cbind(
      data.frame(
        T_generations = T_generations,
        N0 = N0,
        panel_seed = panel_seed,
        J = ancestry_J,
        N = ancestry_N,
        true_L = ancestry_true_L,
        h2 = ancestry_h2,
        method = names(fits)[method_index],
        fit_seconds = runtimes[method_index],
        causal_1 = simulation$causal[1L],
        causal_2 = simulation$causal[2L],
        causal_LD = simulation$causal_LD,
        expected_z_stronger = max(abs(simulation$expected_z)),
        expected_z_weaker = min(abs(simulation$expected_z)),
        stringsAsFactors = FALSE
      ),
      mismatch,
      fit_summary
    )
  }
}

results <- do.call(rbind, rows)
stopifnot(nrow(results) == length(ancestry_N0_grid) * 4L)
temporary_file <- paste0(output_file, ".tmp")
utils::write.csv(results, temporary_file, row.names = FALSE)
if (!file.rename(temporary_file, output_file)) {
  stop("Could not atomically move worker output into place")
}
message("Wrote ", output_file)
