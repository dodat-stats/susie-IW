#!/usr/bin/env Rscript

library(data.table)
source(file.path("experiments", "simulation-iw-collapsed-R.R"))
source(file.path("code", "01-collapsed-R.R"))
source(file.path("code", "02-AIW-and-AIW-N.R"))
source(file.path("code", "03-SER-fallback-and-plot.R"))
source(file.path("code", "credible-set-distance.R"))

numeric_grid <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.numeric(default))
  as.numeric(strsplit(value, ",", fixed = TRUE)[[1L]])
}

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) {
  stop("Usage: Rscript experiments/coalescent-susie-iw-benchmark-worker.R J PANEL_SEED")
}
J <- as.integer(arguments[1L])
panel_seed <- as.integer(arguments[2L])

true_eta0_grid <- numeric_grid(
  "COAL_BENCH_TRUE_ETA0", c(20, 50, 100, 200, 500, 1000, 2000, 5000)
)
true_L_grid <- as.integer(numeric_grid("COAL_BENCH_TRUE_L", 1:3))
lambda_grid <- numeric_grid(
  "COAL_BENCH_LAMBDA_GRID",
  c(1e-5, 1e-4, 1e-3, 2e-3, 4e-3, 6e-3, 8e-3, 1e-2)
)
true_lambda <- as.numeric(Sys.getenv("COAL_BENCH_TRUE_LAMBDA", "0.001"))
N <- as.integer(Sys.getenv("COAL_BENCH_N", "200000"))
N0 <- as.integer(Sys.getenv("COAL_BENCH_N0", "2000"))
coalescent_N <- as.integer(Sys.getenv("COAL_BENCH_COAL_N", "2000"))
h2 <- as.numeric(Sys.getenv("COAL_BENCH_H2", "0.001"))
maximum_L <- as.integer(Sys.getenv("COAL_BENCH_MAX_L", "5"))
collapsed_multiplier <- as.numeric(Sys.getenv(
  "COAL_BENCH_COLLAPSED_MULTIPLIER", "1.677"
))
aiwn_multiplier <- as.numeric(Sys.getenv(
  "COAL_BENCH_AIWN_MULTIPLIER", "1.424"
))
output_directory <- Sys.getenv(
  "COAL_BENCH_OUTPUT_DIRECTORY",
  file.path("experiments", "coalescent-susie-iw-benchmark-20rep-workers")
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(
  output_directory, sprintf("results-J%d-panel%d.csv", J, panel_seed)
)

stopifnot(
  J >= 20L, panel_seed > 0L, all(true_eta0_grid > 0),
  all(true_L_grid %in% 1:3), maximum_L >= max(true_L_grid),
  true_lambda %in% lambda_grid, N > 0L, N0 > 1L, h2 > 0
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
    beta = beta,
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
    fallback = FALSE,
    partition = FALSE,
    model_L = sum(as.numeric(fit$V) > 0),
    raw_eta0 = NA_real_, eta0_used = NA_real_, lambda = NA_real_,
    converged = isTRUE(fit$converged)
  )
}

proposed_details <- function(fit) {
  postprocess <- fit$postprocess
  list(
    sets = susie_get_reported_cs(fit),
    components = iw_reported_components(
      fit$raw_fit, postprocess$raw_credible_sets, postprocess
    ),
    pip = as.numeric(fit$pip),
    raw_L = length(fit$raw_sets),
    fallback = isTRUE(postprocess$fallback_used),
    partition = isTRUE(postprocess$partition_used),
    model_L = if (!is.null(fit$selected_L)) {
      fit$selected_L
    } else {
      sum(as.numeric(fit$raw_fit$V) > 0)
    },
    raw_eta0 = if (!is.null(fit$selected_raw_eta0)) {
      fit$selected_raw_eta0
    } else fit$raw_eta0,
    eta0_used = fit$eta0_used,
    lambda = fit$selected_lambda,
    converged = isTRUE(fit$raw_fit$converged)
  )
}

summarize_details <- function(
    details, method, causal, reference_components, seconds) {
  causal_union <- unique(unlist(details$sets, use.names = FALSE))
  causal_sets <- if (length(details$sets)) {
    vapply(details$sets, function(set) any(set %in% causal), logical(1L))
  } else logical()
  distance <- iw_component_distance(details$components, reference_components)
  data.table(
    method = method,
    model_L = details$model_L,
    raw_L = details$raw_L,
    reported_L = length(details$sets),
    fallback = details$fallback,
    partitioned_SER = details$partition,
    selected_lambda = details$lambda,
    raw_eta0 = details$raw_eta0,
    eta0_used = details$eta0_used,
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
methods <- c("SuSiE-RSS: in-sample R", "SuSiE-RSS: reference R0",
             "collapsed-R", "AIW-N")

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
    completed <- if (nrow(existing)) {
      existing[
        abs(get("true_eta0") - eta0_value) < 1e-12 &
          get("true_L") == L_value,
        unique(method)
      ]
    } else character()
    if (setequal(completed, methods)) next

    message(sprintf(
      "J=%d panel=%d true_eta0=%g true_L=%d", J, panel_seed,
      true_eta0, true_L
    ))
    simulation <- simulate_summary_statistics(
      R, causal_candidates, true_L, noise
    )
    causal <- simulation$causal

    in_time <- system.time(in_fit <- fit_susie_rss(
      simulation$v, R, N, L = maximum_L
    ))
    out_time <- system.time(out_fit <- fit_susie_rss(
      simulation$v, R0, N, L = maximum_L
    ))
    collapsed_time <- system.time(collapsed_fit <- susie_iw(
      v = simulation$v, R0 = R0, N = N, L = maximum_L,
      lambda_grid = lambda_grid, eta0_multiplier = collapsed_multiplier,
      scaled_prior_variance = 0.2, estimate_prior_variance = TRUE,
      coverage = 0.95, min_abs_corr = 0.5, ser_fallback = TRUE,
      partition_ser = FALSE,
      max_iter = 200L, tol = 1e-4, verbose = FALSE
    ))
    aiwn_time <- system.time(aiwn_fit <- suppressMessages(susie_aiw(
      v = simulation$v, R0 = R0, N = N, L = maximum_L,
      approximation = "N", lambda_grid = lambda_grid,
      eta0_multiplier = aiwn_multiplier,
      scaled_prior_variance = 0.2, estimate_prior_variance = TRUE,
      coverage = 0.95, min_abs_corr = 0.5, ser_fallback = TRUE,
      partition_ser = FALSE,
      max_iter = 500L, tol = 1e-3, verbose = FALSE
    )))

    in_details <- baseline_details(in_fit, R)
    reference_components <- in_details$components
    rows <- rbindlist(list(
      summarize_details(
        in_details, methods[1L], causal, reference_components,
        unname(in_time["elapsed"])
      ),
      summarize_details(
        baseline_details(out_fit, R0), methods[2L], causal,
        reference_components, unname(out_time["elapsed"])
      ),
      summarize_details(
        proposed_details(collapsed_fit), methods[3L], causal,
        reference_components, unname(collapsed_time["elapsed"])
      ),
      summarize_details(
        proposed_details(aiwn_fit), methods[4L], causal,
        reference_components, unname(aiwn_time["elapsed"])
      )
    ), fill = TRUE)
    causal_ld <- R[causal, causal, drop = FALSE]
    rows[, `:=`(
      J = J, panel_seed = panel_seed, true_eta0 = true_eta0,
      true_L = true_L, true_lambda = true_lambda, h2 = h2,
      N = N, N0 = N0, causal = paste(causal, collapse = ";"),
      maximum_abs_causal_ld = if (true_L == 1L) 0 else {
        max(abs(causal_ld[upper.tri(causal_ld)]))
      },
      minimum_expected_abs_z = min(abs(simulation$expected_z)),
      maximum_expected_abs_z = max(abs(simulation$expected_z))
    )]
    setcolorder(rows, c(
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
    existing <- rbindlist(list(existing, rows), fill = TRUE)
    setorder(existing, true_eta0, true_L, method)
    fwrite(existing, output_file)
    rm(in_fit, out_fit, collapsed_fit, aiwn_fit)
    gc()
  }
  rm(R)
  gc()
}

stopifnot(
  nrow(existing) == length(true_eta0_grid) * length(true_L_grid) * length(methods)
)
message("Completed ", basename(output_file))
