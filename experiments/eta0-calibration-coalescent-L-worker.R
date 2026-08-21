#!/usr/bin/env Rscript

library(data.table)
source(file.path("experiments", "simulation-iw-collapsed-R.R"))
source(file.path("code", "01-collapsed-R.R"))
source(file.path("code", "02-AIW-and-AIW-N.R"))
source(file.path("code", "03-SER-fallback-and-plot.R"))

numeric_grid <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.numeric(default))
  as.numeric(strsplit(value, ",", fixed = TRUE)[[1L]])
}

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) {
  stop("Usage: Rscript experiments/eta0-calibration-coalescent-L-worker.R J PANEL_SEED")
}
J <- as.integer(arguments[1L])
panel_seed <- as.integer(arguments[2L])

true_eta0_grid <- numeric_grid(
  "ETA0_L_TRUE_GRID", c(20, 50, 100, 200, 500, 1000, 2000, 5000)
)
true_L_grid <- as.integer(numeric_grid("ETA0_L_TRUE_L_GRID", 1:3))
lambda_grid <- numeric_grid(
  "ETA0_L_LAMBDA_GRID",
  c(1e-5, 1e-4, 1e-3, 2e-3, 4e-3, 6e-3, 8e-3, 1e-2)
)
true_lambda <- as.numeric(Sys.getenv("ETA0_L_TRUE_LAMBDA", "0.001"))
sample_size <- as.integer(Sys.getenv("ETA0_L_N", "200000"))
heritability <- as.numeric(Sys.getenv("ETA0_L_H2", "0.001"))
coalescent_sample_size <- as.integer(Sys.getenv("ETA0_L_N_COAL", "2000"))
reference_sample_size <- as.integer(Sys.getenv("ETA0_L_N0", "2000"))
maximum_fitted_L <- as.integer(Sys.getenv("ETA0_L_MAX_FITTED_L", "5"))
output_directory <- Sys.getenv(
  "ETA0_L_OUTPUT_DIRECTORY",
  file.path("experiments", "eta0-calibration-coalescent-L-workers")
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(
  output_directory, sprintf("results-J%d-panel%d.csv", J, panel_seed)
)

stopifnot(
  J >= 20L, panel_seed > 0L, all(true_eta0_grid > 0),
  setequal(true_L_grid, 1:3), true_lambda %in% lambda_grid,
  maximum_fitted_L >= max(true_L_grid), heritability > 0
)

R0 <- simulate_coalescent_reference_ld(
  coalescent_sample_size = coalescent_sample_size,
  reference_sample_size = reference_sample_size,
  number_of_variants = J,
  seed = panel_seed,
  heritability = heritability
)

# Select three separated variants greedily so that the L=1,2,3 designs are
# nested and all selected pairs have weak LD.
choose_nested_weak_ld_variants <- function(R) {
  J <- nrow(R)
  candidates <- unique(as.integer(round(seq(
    max(2, 0.08 * J), min(J - 1, 0.92 * J),
    length.out = min(250L, J - 2L)
  ))))
  minimum_separation <- max(8L, floor(J / 12))
  first <- candidates[which.min(abs(candidates - round(J / 4)))]
  available <- candidates[abs(candidates - first) >= minimum_separation]
  second <- available[which.min(abs(R[first, available]))]
  available <- candidates[
    abs(candidates - first) >= minimum_separation &
      abs(candidates - second) >= minimum_separation
  ]
  score <- vapply(available, function(j) {
    max(abs(R[j, c(first, second)]))
  }, numeric(1L))
  third <- available[which.min(score)]
  c(first, second, third)
}

simulate_v <- function(R, causal, true_L, noise) {
  weights <- c(7, 5, 3)[seq_len(true_L)]
  active <- causal[seq_len(true_L)]
  active_R <- R[active, active, drop = FALSE]
  scale <- sqrt(heritability / drop(crossprod(
    weights, active_R %*% weights
  )))
  beta <- numeric(nrow(R))
  beta[active] <- weights * scale
  list(
    v = as.numeric(R %*% beta + noise),
    beta = beta,
    causal = active,
    achieved_h2 = drop(crossprod(beta, R %*% beta)),
    expected_z = sqrt(sample_size) * as.numeric(R[active, , drop = FALSE] %*% beta)
  )
}

existing <- if (file.exists(output_file)) fread(output_file) else data.table()

for (true_eta0 in true_eta0_grid) {
  eta0_value <- true_eta0
  eta_seed <- as.integer(1000000 + 10000 * panel_seed + J + true_eta0)
  set.seed(eta_seed)
  Rbar_true <- make_stabilized_reference_ld(R0, true_lambda)
  R <- standardize_correlation(sample_inverse_wishart(
    degrees_of_freedom = true_eta0 + J + 1,
    scale_matrix = true_eta0 * Rbar_true
  ))
  causal <- choose_nested_weak_ld_variants(R)
  causal_ld <- R[causal, causal, drop = FALSE]
  set.seed(eta_seed + 500000L)
  noise <- as.numeric(t(stable_cholesky(R)) %*% rnorm(J)) /
    sqrt(sample_size)

  for (true_L in true_L_grid) {
    L_value <- true_L
    if (nrow(existing[
      abs(get("true_eta0") - eta0_value) < 1e-12 &
        get("true_L") == L_value
    ]) == 2L) next

    message(sprintf(
      "J=%d panel=%d true_eta0=%g true_L=%d",
      J, panel_seed, true_eta0, true_L
    ))
    simulated <- simulate_v(R, causal, true_L, noise)

    collapsed_time <- system.time(collapsed <- susie_iw(
      v = simulated$v, R0 = R0, N = sample_size,
      L = maximum_fitted_L, lambda_grid = lambda_grid,
      eta0_multiplier = 1, scaled_prior_variance = 0.2,
      estimate_prior_variance = TRUE, ser_fallback = FALSE,
      max_iter = 200L, tol = 1e-4, verbose = FALSE
    ))
    aiwn_time <- system.time(aiwn <- suppressMessages(susie_aiw(
      v = simulated$v, R0 = R0, N = sample_size,
      L = maximum_fitted_L, approximation = "N",
      lambda_grid = lambda_grid, eta0_multiplier = 1,
      scaled_prior_variance = 0.2, estimate_prior_variance = TRUE,
      ser_fallback = FALSE, max_iter = 500L, tol = 1e-3,
      verbose = FALSE
    )))

    common <- list(
      J = J, panel_seed = panel_seed, true_eta0 = true_eta0,
      true_L = true_L, true_lambda = true_lambda,
      h2 = heritability, N = sample_size, N0 = reference_sample_size,
      causal = paste(simulated$causal, collapse = ";"),
      maximum_abs_causal_ld = if (true_L == 1L) 0 else max(abs(
        causal_ld[seq_len(true_L), seq_len(true_L), drop = FALSE][
          upper.tri(causal_ld[seq_len(true_L), seq_len(true_L), drop = FALSE])
        ]
      )),
      minimum_expected_abs_z = min(abs(simulated$expected_z)),
      maximum_expected_abs_z = max(abs(simulated$expected_z)),
      achieved_h2 = simulated$achieved_h2
    )
    new_rows <- rbindlist(list(
      as.data.table(c(common, list(
        method = "collapsed-R",
        raw_eta0 = collapsed$raw_eta0,
        selected_lambda = collapsed$selected_lambda,
        selected_model_L = maximum_fitted_L,
        raw_pure_L = length(collapsed$raw_sets),
        eta0_boundary = if (identical(
          collapsed$raw_fit$nu0_boundary, "upper"
        )) {
          "normal_limit"
        } else if (identical(collapsed$raw_fit$nu0_boundary, "lower")) {
          "optimization_boundary"
        } else "interior",
        converged = isTRUE(collapsed$raw_fit$converged),
        fit_seconds = unname(collapsed_time["elapsed"])
      ))),
      as.data.table(c(common, list(
        method = "AIW-N",
        raw_eta0 = aiwn$selected_raw_eta0,
        selected_lambda = aiwn$selected_lambda,
        selected_model_L = aiwn$selected_L,
        raw_pure_L = length(aiwn$raw_sets),
        eta0_boundary = if (is.infinite(aiwn$selected_raw_eta0)) {
          "normal_limit"
        } else "interior",
        converged = isTRUE(aiwn$raw_fit$converged),
        fit_seconds = unname(aiwn_time["elapsed"])
      )))
    ), fill = TRUE)
    existing <- existing[!(
      abs(get("true_eta0") - eta0_value) < 1e-12 &
        get("true_L") == L_value
    )]
    existing <- rbindlist(list(existing, new_rows), fill = TRUE)
    setorder(existing, true_eta0, true_L, method)
    fwrite(existing, output_file)
    rm(collapsed, aiwn)
    gc()
  }
  rm(R)
  gc()
}

stopifnot(nrow(existing) == length(true_eta0_grid) *
  length(true_L_grid) * 2L)
message("Completed ", basename(output_file))
