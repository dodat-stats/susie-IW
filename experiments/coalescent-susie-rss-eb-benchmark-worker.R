#!/usr/bin/env Rscript

# Exploratory extension of the coalescent SuSiE-IW benchmark. This worker
# regenerates each frozen simulation deterministically and fits only the new
# susieR finite-reference + empirical-Bayes mismatch method. It does not alter
# the benchmark files used by the manuscript.

library(data.table)
source(file.path("experiments", "simulation-iw-collapsed-R.R"))
source(file.path("code", "03-SER-fallback-and-plot.R"))
source(file.path("code", "credible-set-distance.R"))

numeric_grid <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.numeric(default))
  as.numeric(strsplit(value, ",", fixed = TRUE)[[1L]])
}

scalar_or_na <- function(value, mode = c("numeric", "logical", "character")) {
  mode <- match.arg(mode)
  if (is.null(value) || !length(value)) {
    return(switch(mode, numeric = NA_real_, logical = NA, character = NA_character_))
  }
  switch(
    mode,
    numeric = as.numeric(value[1L]),
    logical = as.logical(value[1L]),
    character = as.character(value[1L])
  )
}

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) {
  stop(
    "Usage: Rscript experiments/coalescent-susie-rss-eb-benchmark-worker.R ",
    "J PANEL_SEED"
  )
}
J <- as.integer(arguments[1L])
panel_seed <- as.integer(arguments[2L])

true_eta0_grid <- numeric_grid(
  "COAL_BENCH_TRUE_ETA0", c(20, 50, 100, 200, 500, 1000, 2000, 5000)
)
true_L_grid <- as.integer(numeric_grid("COAL_BENCH_TRUE_L", 1:3))
true_lambda <- as.numeric(Sys.getenv("COAL_BENCH_TRUE_LAMBDA", "0.001"))
N <- as.integer(Sys.getenv("COAL_BENCH_N", "200000"))
N0 <- as.integer(Sys.getenv("COAL_BENCH_N0", "2000"))
coalescent_N <- as.integer(Sys.getenv("COAL_BENCH_COAL_N", "2000"))
h2 <- as.numeric(Sys.getenv("COAL_BENCH_H2", "0.001"))
maximum_L <- as.integer(Sys.getenv("COAL_BENCH_MAX_L", "5"))
output_directory <- Sys.getenv(
  "COAL_EB_OUTPUT_DIRECTORY",
  file.path("experiments", "coalescent-susie-rss-eb-benchmark-20rep-workers")
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(
  output_directory, sprintf("results-J%d-panel%d.csv", J, panel_seed)
)

method_name <- "SuSiE-RSS: finite R0 + EB mismatch"
required_arguments <- c("R_finite", "R_mismatch")
missing_arguments <- setdiff(
  required_arguments, names(formals(susieR::susie_rss))
)
if (length(missing_arguments)) {
  stop(
    "The installed susieR ", as.character(packageVersion("susieR")),
    " does not support ", paste(missing_arguments, collapse = ", "),
    ". Install the current development version before running this benchmark."
  )
}
stopifnot(
  J >= 20L, panel_seed > 0L, all(true_eta0_grid > 0),
  all(true_L_grid %in% 1:3), maximum_L >= max(true_L_grid),
  N > 0L, N0 > 1L, h2 > 0
)

R0 <- simulate_coalescent_reference_ld(
  coalescent_sample_size = coalescent_N,
  reference_sample_size = N0,
  number_of_variants = J,
  seed = panel_seed,
  heritability = h2
)

choose_nested_weak_ld_variants <- function(R) {
  candidates <- unique(as.integer(round(seq(
    max(2, 0.08 * nrow(R)), min(nrow(R) - 1, 0.92 * nrow(R)),
    length.out = min(250L, nrow(R) - 2L)
  ))))
  minimum_separation <- max(8L, floor(nrow(R) / 12))
  first <- candidates[which.min(abs(candidates - round(nrow(R) / 4)))]
  available <- candidates[abs(candidates - first) >= minimum_separation]
  second <- available[which.min(abs(R[first, available]))]
  available <- candidates[
    abs(candidates - first) >= minimum_separation &
      abs(candidates - second) >= minimum_separation
  ]
  third <- available[which.min(vapply(available, function(j) {
    max(abs(R[j, c(first, second)]))
  }, numeric(1L)))]
  c(first, second, third)
}

simulate_summary_statistics <- function(R, causal, true_L, noise) {
  weights <- c(7, 5, 3)[seq_len(true_L)]
  active <- causal[seq_len(true_L)]
  beta <- numeric(nrow(R))
  beta[active] <- weights * sqrt(
    h2 / drop(crossprod(weights, R[active, active, drop = FALSE] %*% weights))
  )
  list(
    v = as.numeric(R %*% beta + noise),
    causal = active,
    expected_z = sqrt(N) * as.numeric(R[active, , drop = FALSE] %*% beta)
  )
}

baseline_details <- function(fit, R) {
  credible_sets <- susieR::susie_get_cs(
    fit, Xcorr = R, coverage = 0.95, min_abs_corr = 0.5
  )
  list(
    sets = iw_cs_list(credible_sets),
    components = iw_reported_components(fit, credible_sets),
    pip = as.numeric(fit$pip),
    raw_L = length(iw_cs_list(credible_sets)),
    model_L = sum(as.numeric(fit$V) > 0),
    converged = isTRUE(fit$converged)
  )
}

summarize_details <- function(details, causal, reference_components, seconds) {
  causal_union <- unique(unlist(details$sets, use.names = FALSE))
  causal_sets <- if (length(details$sets)) {
    vapply(details$sets, function(set) any(set %in% causal), logical(1L))
  } else logical()
  distance <- iw_component_distance(details$components, reference_components)
  data.table(
    method = method_name,
    model_L = details$model_L,
    raw_L = details$raw_L,
    reported_L = length(details$sets),
    fallback = FALSE,
    partitioned_SER = FALSE,
    selected_lambda = NA_real_,
    raw_eta0 = NA_real_,
    eta0_used = NA_real_,
    converged = details$converged,
    fit_seconds = seconds,
    causal_coverage = mean(causal %in% causal_union),
    all_causal_covered = all(causal %in% causal_union),
    false_sets = sum(!causal_sets),
    true_sets = sum(causal_sets),
    mean_set_size = if (length(details$sets)) mean(lengths(details$sets)) else NA_real_,
    mean_causal_pip = mean(details$pip[causal]),
    d_cs = distance$distance,
    d_cs_normalized = distance$normalized_distance
  )
}

existing <- if (file.exists(output_file)) fread(output_file) else data.table()

for (true_eta0 in true_eta0_grid) {
  eta0_value <- true_eta0
  simulation_seed <- as.integer(2000000 + 10000 * panel_seed + J + true_eta0)
  set.seed(simulation_seed)
  Rbar_true <- make_stabilized_reference_ld(R0, true_lambda)
  R <- standardize_correlation(sample_inverse_wishart(
    degrees_of_freedom = true_eta0 + J + 1,
    scale_matrix = true_eta0 * Rbar_true
  ))
  causal_candidates <- choose_nested_weak_ld_variants(R)
  set.seed(simulation_seed + 500000L)
  noise <- as.numeric(t(stable_cholesky(R)) %*% rnorm(J)) / sqrt(N)

  for (true_L in true_L_grid) {
    L_value <- true_L
    already_complete <- if (nrow(existing)) {
      nrow(existing[
        abs(get("true_eta0") - eta0_value) < 1e-12 &
          get("true_L") == L_value & method == method_name
      ]) == 1L
    } else FALSE
    if (already_complete) next

    message(sprintf(
      "J=%d panel=%d true_eta0=%g true_L=%d", J, panel_seed,
      true_eta0, true_L
    ))
    simulation <- simulate_summary_statistics(
      R, causal_candidates, true_L, noise
    )
    causal <- simulation$causal

    # This inexpensive reference fit is used only to retain the benchmark's
    # component-distance metric. Its runtime is not charged to the EB method.
    in_fit <- fit_susie_rss(simulation$v, R, N, L = maximum_L)
    reference_components <- baseline_details(in_fit, R)$components

    z <- iw_z_from_v(simulation$v, N)
    eb_time <- system.time(eb_fit <- suppressMessages(susieR::susie_rss(
      z = z,
      R = R0,
      n = N,
      L = maximum_L,
      scaled_prior_variance = 0.2,
      estimate_prior_variance = TRUE,
      coverage = 0.95,
      min_abs_corr = 0.5,
      max_iter = 500L,
      tol = 1e-3,
      R_finite = N0,
      R_mismatch = "eb",
      track_fit = FALSE
    )))
    diagnostics <- eb_fit$R_finite_diagnostics
    row <- summarize_details(
      baseline_details(eb_fit, R0), causal, reference_components,
      unname(eb_time["elapsed"])
    )
    row[, `:=`(
      susieR_version = as.character(packageVersion("susieR")),
      R_finite = scalar_or_na(diagnostics$B),
      R_mismatch = scalar_or_na(diagnostics$R_mismatch, "character"),
      mismatch_estimator = scalar_or_na(eb_fit$R_mismatch_method, "character"),
      lambda_bias = scalar_or_na(eb_fit$lambda_bias),
      corrected_reference_size = scalar_or_na(eb_fit$B_corrected),
      Q_art = scalar_or_na(diagnostics$Q_art),
      residual_artifact_flag = scalar_or_na(
        diagnostics$artifact_flag, "logical"
      ),
      overall_reliability_flag = scalar_or_na(
        diagnostics$R_reliability_flag, "logical"
      ),
      convergence_reason = scalar_or_na(
        eb_fit$convergence_reason, "character"
      ),
      iterations = scalar_or_na(eb_fit$niter)
    )]
    causal_ld <- R[causal, causal, drop = FALSE]
    row[, `:=`(
      J = J, panel_seed = panel_seed, true_eta0 = true_eta0,
      true_L = L_value, true_lambda = true_lambda, h2 = h2,
      N = N, N0 = N0, causal = paste(causal, collapse = ";"),
      maximum_abs_causal_ld = if (true_L == 1L) 0 else {
        max(abs(causal_ld[upper.tri(causal_ld)]))
      },
      minimum_expected_abs_z = min(abs(simulation$expected_z)),
      maximum_expected_abs_z = max(abs(simulation$expected_z))
    )]
    setcolorder(row, c(
      "J", "panel_seed", "true_eta0", "true_L", "true_lambda",
      "h2", "N", "N0", "causal", "maximum_abs_causal_ld",
      "minimum_expected_abs_z", "maximum_expected_abs_z", "method"
    ))
    if (nrow(existing)) {
      existing <- existing[!(
        abs(get("true_eta0") - eta0_value) < 1e-12 &
          get("true_L") == L_value
      )]
    }
    existing <- rbindlist(list(existing, row), fill = TRUE)
    setorder(existing, true_eta0, true_L)
    fwrite(existing, output_file)
    rm(in_fit, eb_fit)
    gc()
  }
  rm(R)
  gc()
}

stopifnot(
  nrow(existing) == length(true_eta0_grid) * length(true_L_grid),
  all(existing$method == method_name)
)
message("Completed ", basename(output_file))
