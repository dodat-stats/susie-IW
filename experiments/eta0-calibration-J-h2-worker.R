#!/usr/bin/env Rscript

library(data.table)
source(file.path("experiments", "simulation-iw-collapsed-R.R"))
source(file.path("code", "02-AIW-and-AIW-N.R"))

numeric_grid <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.numeric(default))
  as.numeric(strsplit(value, ",", fixed = TRUE)[[1L]])
}

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) {
  stop("Usage: Rscript experiments/eta0-calibration-J-h2-worker.R J PANEL_SEED")
}
J <- as.integer(arguments[1L])
panel_seed <- as.integer(arguments[2L])
if (!is.finite(J) || J < 20L || !is.finite(panel_seed)) {
  stop("J and PANEL_SEED must be valid integers.")
}

true_eta0_grid <- numeric_grid(
  "ETA0_CAL_TRUE_GRID", c(20, 50, 100, 200, 500, 1000)
)
h2_grid <- numeric_grid(
  "ETA0_CAL_H2_GRID", c(2e-4, 1e-3, 5e-3)
)
lambda_grid <- numeric_grid(
  "ETA0_CAL_LAMBDA_GRID",
  c(1e-5, 1e-4, 1e-3, 2e-3, 4e-3, 6e-3, 8e-3, 1e-2)
)
true_lambda <- as.numeric(Sys.getenv(
  "ETA0_CAL_TRUE_LAMBDA", unset = "0.001"
))
sample_size <- as.integer(Sys.getenv("ETA0_CAL_N", unset = "200000"))
coalescent_sample_size <- as.integer(Sys.getenv(
  "ETA0_CAL_N_COAL", unset = "2000"
))
reference_sample_size <- as.integer(Sys.getenv(
  "ETA0_CAL_N0", unset = "2000"
))
output_directory <- Sys.getenv(
  "ETA0_CAL_OUTPUT_DIRECTORY",
  unset = file.path(
    susie_iw_root, "analysis", "eta0-calibration-J-h2-workers"
  )
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(
  output_directory, sprintf("results-J%d-panel%d.csv", J, panel_seed)
)

stopifnot(
  true_lambda %in% lambda_grid,
  all(true_eta0_grid > 0),
  all(h2_grid > 0),
  all(lambda_grid > 0 & lambda_grid <= 1e-2),
  !anyDuplicated(lambda_grid)
)

existing <- if (file.exists(output_file)) {
  fread(output_file)
} else {
  data.table()
}

R0 <- simulate_coalescent_reference_ld(
  coalescent_sample_size = coalescent_sample_size,
  reference_sample_size = reference_sample_size,
  number_of_variants = J,
  seed = panel_seed,
  heritability = 1e-3
)

reference_models <- lapply(lambda_grid, function(lambda) {
  Rbar <- make_stabilized_reference_ld(R0, lambda)
  chol_Rbar <- chol(Rbar)
  list(
    lambda = lambda,
    Rbar = Rbar,
    chol_Rbar = chol_Rbar,
    logdetRbar = 2 * sum(log(diag(chol_Rbar)))
  )
})

choose_weak_ld_pair <- function(R) {
  candidates <- unique(as.integer(round(seq(
    max(2, 0.1 * nrow(R)), min(nrow(R) - 1, 0.9 * nrow(R)),
    length.out = min(200L, nrow(R) - 2L)
  ))))
  absolute_ld <- abs(R[candidates, candidates, drop = FALSE])
  separation <- abs(outer(candidates, candidates, "-"))
  absolute_ld[separation < max(10L, floor(nrow(R) / 10))] <- NA_real_
  diag(absolute_ld) <- NA_real_
  pair <- which(
    absolute_ld == min(absolute_ld, na.rm = TRUE),
    arr.ind = TRUE
  )[1L, ]
  sort(candidates[pair])
}

simulate_marginal_vector <- function(R, causal, h2, noise) {
  R_causal <- R[causal, causal, drop = FALSE]
  effect_template <- c(7, sign(R_causal[1L, 2L]) * 3)
  template_variance <- drop(crossprod(
    effect_template, R_causal %*% effect_template
  ))
  beta_causal <- effect_template * sqrt(h2 / template_variance)
  beta <- numeric(nrow(R))
  beta[causal] <- beta_causal
  list(
    v = as.numeric(R %*% beta + noise),
    beta = beta,
    achieved_h2 = drop(crossprod(
      beta_causal, R_causal %*% beta_causal
    )),
    expected_z_1 = sqrt(sample_size) *
      drop(R_causal[1L, ] %*% beta_causal),
    expected_z_2 = sqrt(sample_size) *
      drop(R_causal[2L, ] %*% beta_causal)
  )
}

select_aiwn_row <- function(scores, subset_rows) {
  candidates <- scores[subset_rows, , drop = FALSE]
  best <- max(candidates$profiled_objective)
  eligible <- which(candidates$profiled_objective >= best - 1e-3)
  candidates[eligible[
    order(candidates$L[eligible], candidates$lambda[eligible])[1L]
  ], , drop = FALSE]
}

make_aiwn_row <- function(
    score, stage, true_eta0, h2, causal_correlation, simulation,
    elapsed_seconds) {
  data.table(
    method = "AIW-N",
    estimation_stage = stage,
    J = J,
    h2 = h2,
    true_eta0 = true_eta0,
    true_lambda = true_lambda,
    panel_seed = panel_seed,
    raw_eta0 = score$raw_eta0,
    selected_lambda = score$lambda,
    selected_L = score$L,
    eta0_finite = is.finite(score$raw_eta0) && score$raw_eta0 > 0,
    eta0_boundary = if (is.infinite(score$raw_eta0)) {
      "normal_limit"
    } else if (is.na(score$raw_eta0)) {
      "missing"
    } else {
      "interior"
    },
    objective = score$profiled_objective,
    converged = TRUE,
    causal_correlation = causal_correlation,
    achieved_h2 = simulation$achieved_h2,
    expected_z_1 = simulation$expected_z_1,
    expected_z_2 = simulation$expected_z_2,
    fit_seconds = elapsed_seconds
  )
}

fit_collapsed_profiles <- function(v) {
  rows <- vector("list", length(reference_models))
  for (index in seq_along(reference_models)) {
    model <- reference_models[[index]]
    inverse_times_v <- backsolve(
      model$chol_Rbar,
      forwardsolve(t(model$chol_Rbar), v)
    )
    pre <- list(
      r0 = as.numeric(crossprod(v, inverse_times_v)),
      logdetRbar = model$logdetRbar
    )
    timing <- system.time(fit <-
      fit_susie_collapse_r_shared_veb(
        x = v,
        Rbar = model$Rbar,
        N = sample_size,
        L = 5L,
        sigma2 = 0.2,
        nu0_init = 100,
        nu0_bounds = c(1, 1e6),
        nu0_coarse_grid_size = 15L,
        estimate_sigma2 = TRUE,
        max_iter = 200L,
        tol = 1e-4,
        verbose = FALSE,
        pre = pre
      ))
    rows[[index]] <- data.table(
      lambda = model$lambda,
      raw_eta0 = fit$nu0,
      objective = fit$lower_bound,
      converged = isTRUE(fit$converged),
      eta0_boundary = if (isTRUE(fit$nu0_boundary)) {
        "optimization_boundary"
      } else {
        "interior"
      },
      selected_L = length(susieR::susie_get_cs(
        fit, Xcorr = model$Rbar, coverage = 0.95,
        min_abs_corr = 0.5, dedup = TRUE
      )$cs),
      fit_seconds = unname(timing["elapsed"])
    )
  }
  rbindlist(rows)
}

make_collapsed_row <- function(
    score, stage, true_eta0, h2, causal_correlation, simulation,
    total_seconds) {
  data.table(
    method = "collapsed-R",
    estimation_stage = stage,
    J = J,
    h2 = h2,
    true_eta0 = true_eta0,
    true_lambda = true_lambda,
    panel_seed = panel_seed,
    raw_eta0 = score$raw_eta0,
    selected_lambda = score$lambda,
    selected_L = score$selected_L,
    eta0_finite = is.finite(score$raw_eta0) && score$raw_eta0 > 0,
    eta0_boundary = score$eta0_boundary,
    objective = score$objective,
    converged = score$converged,
    causal_correlation = causal_correlation,
    achieved_h2 = simulation$achieved_h2,
    expected_z_1 = simulation$expected_z_1,
    expected_z_2 = simulation$expected_z_2,
    fit_seconds = total_seconds
  )
}

for (true_eta0 in true_eta0_grid) {
  eta0_value <- true_eta0
  complete_h2 <- if (nrow(existing)) {
    existing[
      true_eta0 == eta0_value &
        method == "AIW-N" &
        estimation_stage == "actual_selection",
      unique(h2)
    ]
  } else {
    numeric()
  }
  if (setequal(complete_h2, h2_grid)) {
    message(sprintf(
      "J=%d panel=%d eta0=%g already complete", J, panel_seed, true_eta0
    ))
    next
  }

  true_Rbar <- reference_models[[
    match(true_lambda, lambda_grid)
  ]]$Rbar
  covariance_seed <- as.integer(
    1000000 + 10000 * panel_seed + J + round(true_eta0)
  )
  set.seed(covariance_seed)
  in_sample_ld <- standardize_correlation(sample_inverse_wishart(
    degrees_of_freedom = true_eta0 + J + 1,
    scale_matrix = true_eta0 * true_Rbar
  ))
  causal <- choose_weak_ld_pair(in_sample_ld)
  causal_correlation <- in_sample_ld[causal[1L], causal[2L]]
  noise_seed <- covariance_seed + 500000L
  set.seed(noise_seed)
  noise <- as.numeric(
    t(stable_cholesky(in_sample_ld)) %*% rnorm(J)
  ) / sqrt(sample_size)

  for (h2 in h2_grid) {
    h2_value <- h2
    already_complete <- nrow(existing) && nrow(existing[
      true_eta0 == eta0_value & abs(get("h2") - h2_value) < 1e-15
    ]) == 5L
    if (already_complete) next

    message(sprintf(
      "J=%d panel=%d eta0=%g h2=%g", J, panel_seed, true_eta0, h2
    ))
    simulation <- simulate_marginal_vector(
      in_sample_ld, causal, h2, noise
    )

    aiwn_timing <- system.time(aiwn_fit <- suppressMessages(
      fit_aiwn_lambda_grid(
        v = simulation$v,
        R0 = R0,
        N = sample_size,
        lambda_grid = lambda_grid,
        L_grid = 1:5,
        scaled_prior_variance = 0.2,
        selection_tolerance = 1e-3,
        cs_coverage = 0.95,
        cs_min_abs_corr = 0.5,
        susie_args = list(max_iter = 500L, tol = 1e-3),
        reference_models = reference_models,
        eta0_multiplier = 1
      )
    ))
    aiwn_scores <- as.data.table(aiwn_fit$scores)
    oracle_aiwn <- aiwn_scores[
      abs(lambda - true_lambda) < 1e-15 & L == 2L
    ]
    lambda_selected_aiwn <- select_aiwn_row(
      aiwn_scores, aiwn_scores$L == 2L
    )
    actual_aiwn <- select_aiwn_row(
      aiwn_scores, rep(TRUE, nrow(aiwn_scores))
    )
    aiwn_rows <- rbindlist(list(
      make_aiwn_row(
        oracle_aiwn, "oracle_lambda_oracle_L", true_eta0, h2,
        causal_correlation, simulation, unname(aiwn_timing["elapsed"])
      ),
      make_aiwn_row(
        lambda_selected_aiwn, "selected_lambda_oracle_L", true_eta0, h2,
        causal_correlation, simulation, unname(aiwn_timing["elapsed"])
      ),
      make_aiwn_row(
        actual_aiwn, "actual_selection", true_eta0, h2,
        causal_correlation, simulation, unname(aiwn_timing["elapsed"])
      )
    ))

    collapsed_profiles <- fit_collapsed_profiles(simulation$v)
    oracle_collapsed <- collapsed_profiles[
      abs(lambda - true_lambda) < 1e-15
    ]
    selected_collapsed <- collapsed_profiles[which.max(objective)]
    collapsed_rows <- rbindlist(list(
      make_collapsed_row(
        oracle_collapsed, "oracle_lambda", true_eta0, h2,
        causal_correlation, simulation,
        sum(collapsed_profiles$fit_seconds)
      ),
      make_collapsed_row(
        selected_collapsed, "actual_selection", true_eta0, h2,
        causal_correlation, simulation,
        sum(collapsed_profiles$fit_seconds)
      )
    ))

    new_rows <- rbindlist(list(aiwn_rows, collapsed_rows), fill = TRUE)
    existing <- existing[!(
      true_eta0 == eta0_value & abs(get("h2") - h2_value) < 1e-15
    )]
    existing <- rbindlist(list(existing, new_rows), fill = TRUE)
    setorder(existing, true_eta0, h2, method, estimation_stage)
    fwrite(existing, output_file)
    rm(aiwn_fit, aiwn_scores, collapsed_profiles)
    gc()
  }
  rm(in_sample_ld)
  gc()
}

expected_rows <- length(true_eta0_grid) * length(h2_grid) * 5L
stopifnot(nrow(existing) == expected_rows)
message("Completed ", basename(output_file))
