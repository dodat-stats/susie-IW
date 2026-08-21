#!/usr/bin/env Rscript

# Joint lambda--eta0 recovery in the inverse-Wishart collapsed-R simulation.
#
# This preserves the easy simulation design in
# experiments/simulation-iw-collapsed-R.R. For each simulated data set, collapsed-R
# is fitted at every candidate lambda while eta0 is profiled by VEB. We also
# refit every candidate at a calibrated multiple of its profiled eta0 and
# select lambda using the corrected-fit ELBO. SuSiE-RSS using in-sample and
# out-of-sample LD are included as benchmarks. Intermediate results are
# checkpointed after every simulated data set.

source(file.path("experiments", "simulation-iw-collapsed-R.R"))

get_numeric_grid <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.numeric(default))
  as.numeric(strsplit(value, ",", fixed = TRUE)[[1L]])
}

get_integer_value <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.integer(default))
  as.integer(value)
}

get_character_value <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) default else value
}

lambda_grid <- get_numeric_grid(
  "JOINT_LAMBDA_GRID",
  c(1e-5, 1e-4, 1e-3, 2e-3, 4e-3, 6e-3, 8e-3, 1e-2)
)
true_lambda_grid <- get_numeric_grid(
  "JOINT_TRUE_LAMBDA_GRID", c(1e-3, 1e-2)
)
true_eta0_grid <- get_numeric_grid(
  "JOINT_TRUE_ETA0_GRID",
  unique(round(exp(seq(log(20), log(10000), length.out = 9))))
)
reference_panel_seeds <- as.integer(get_numeric_grid(
  "JOINT_R0_SEEDS", c(3L, 4L, 5L)
))
number_of_replicates <- get_integer_value("JOINT_N_REP", 5L)
number_of_variants <- get_integer_value("JOINT_J", 500L)
sample_size <- get_integer_value("JOINT_N", 20000L)
coalescent_sample_size <- get_integer_value("JOINT_N_COAL", 2000L)
reference_sample_size <- get_integer_value("JOINT_N0", 2000L)
maximum_iterations <- get_integer_value("JOINT_MAX_ITER", 500L)
simulation_seed <- get_integer_value("JOINT_SIMULATION_SEED", 1L)
output_prefix <- get_character_value(
  "JOINT_OUTPUT_PREFIX",
  file.path(
    susie_iw_root, "analysis",
    "simulation-iw-collapsed-R-joint-lambda-eta0-capped-1.537"
  )
)

heritability <- 0.01
true_number_of_effects <- 2L
fitted_number_of_effects <- 5L
coverage <- 0.95
minimum_absolute_correlation <- 0.5
eta0_bounds <- c(2, 20000)
eta0_initial_value <- 100
eta0_coarse_grid_size <- 21L
convergence_tolerance <- 1e-7
eta0_multiplier <- 1.537

stopifnot(
  all(true_lambda_grid %in% lambda_grid),
  all(lambda_grid > 0 & lambda_grid < 1),
  !anyDuplicated(lambda_grid),
  !anyDuplicated(true_lambda_grid),
  all(true_eta0_grid > 0),
  eta0_multiplier >= 1,
  number_of_replicates >= 1L
)

results_file <- paste0(output_prefix, "-results.csv")
profiles_file <- paste0(output_prefix, "-profiles.csv")
postprocessed_results_file <- paste0(
  output_prefix, "-postprocessed-results.csv"
)
raw_summary_file <- paste0(output_prefix, "-raw-summary.csv")
summary_file <- paste0(output_prefix, "-summary.csv")
eta0_calibration_file <- paste0(
  output_prefix, "-eta0-calibration.csv"
)
output_file <- paste0(output_prefix, "-summary.rds")
figure_prefix <- file.path(
  dirname(output_prefix), "figures", basename(output_prefix)
)
dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(figure_prefix), recursive = TRUE, showWarnings = FALSE)

empty_results <- function() {
  data.frame(
    true_lambda = numeric(), true_eta0 = numeric(), r0_seed = integer(),
    replicate = integer(), seed = integer(), fit_type = character(),
    selected_lambda = numeric(), raw_eta0_hat = numeric(),
    eta0_used = numeric(), eta0_multiplier = numeric(),
    raw_log_eta0_error = numeric(), used_log_eta0_error = numeric(),
    log_lambda_error = numeric(),
    eta0_boundary = character(), lambda_boundary = character(),
    lower_bound = numeric(), elbo_gain_over_matching_oracle = numeric(),
    selected_L = integer(), coverage = numeric(), power = numeric(),
    all_causal_covered = logical(), mean_causal_pip = numeric(),
    minimum_causal_pip = numeric(), converged = logical(),
    iterations = integer(), fit_seconds = numeric(), causal = character(),
    stringsAsFactors = FALSE
  )
}

empty_profiles <- function() {
  data.frame(
    true_lambda = numeric(), true_eta0 = numeric(), r0_seed = integer(),
    replicate = integer(), seed = integer(), profile_type = character(),
    candidate_lambda = numeric(), raw_eta0_hat = numeric(),
    eta0_used = numeric(), eta0_multiplier = numeric(),
    raw_log_eta0_error = numeric(), used_log_eta0_error = numeric(),
    log_lambda_error = numeric(), eta0_boundary = character(),
    lower_bound = numeric(), elbo_below_best = numeric(),
    selected_L = integer(), coverage = numeric(), power = numeric(),
    all_causal_covered = logical(), mean_causal_pip = numeric(),
    minimum_causal_pip = numeric(), converged = logical(),
    iterations = integer(), fit_seconds = numeric(), causal = character(),
    stringsAsFactors = FALSE
  )
}

results <- if (file.exists(results_file)) {
  utils::read.csv(results_file, stringsAsFactors = FALSE)
} else {
  empty_results()
}
profiles <- if (file.exists(profiles_file)) {
  utils::read.csv(profiles_file, stringsAsFactors = FALSE)
} else {
  empty_profiles()
}

dataset_is_complete <- function(true_lambda, true_eta0, r0_seed, replicate) {
  if (nrow(profiles) == 0L) return(FALSE)
  rows <- profiles$true_lambda == true_lambda &
    profiles$true_eta0 == true_eta0 &
    profiles$r0_seed == r0_seed &
    profiles$replicate == replicate
  dataset_profiles <- profiles[rows, , drop = FALSE]
  expected_profile_types <- c("profiled_eta0", "calibrated_eta0")
  profile_complete <- nrow(dataset_profiles) ==
      length(lambda_grid) * length(expected_profile_types) &&
    setequal(dataset_profiles$candidate_lambda, lambda_grid) &&
    setequal(dataset_profiles$profile_type, expected_profile_types) &&
    all(table(dataset_profiles$profile_type) == length(lambda_grid))
  if (!profile_complete || nrow(results) == 0L) return(FALSE)
  result_rows <- results$true_lambda == true_lambda &
    results$true_eta0 == true_eta0 &
    results$r0_seed == r0_seed &
    results$replicate == replicate
  expected_fit_types <- c(
    "susie_rss_in_sample", "susie_rss_out_of_sample",
    "joint_lambda_profiled_eta0", "joint_lambda_calibrated_eta0"
  )
  dataset_results <- results[result_rows, , drop = FALSE]
  nrow(dataset_results) == length(expected_fit_types) &&
    setequal(dataset_results$fit_type, expected_fit_types)
}

prepare_reference_models <- function(R0) {
  lapply(lambda_grid, function(lambda) {
    Rbar <- make_stabilized_reference_ld(R0, lambda = lambda)
    precision <- solve(Rbar)
    list(
      lambda = lambda,
      Rbar = Rbar,
      precision = precision,
      logdet = scr_logdet_spd(Rbar)
    )
  })
}

summarize_fitted_model <- function(fit, Rbar, true_effects) {
  credible_sets <- susieR::susie_get_cs(
    fit,
    Xcorr = Rbar,
    coverage = coverage,
    min_abs_corr = minimum_absolute_correlation,
    dedup = TRUE
  )$cs
  metrics <- credible_set_metrics_from_list(credible_sets, true_effects)
  causal_variants <- which(true_effects != 0)
  c(
    selected_L = metrics$selected_L,
    coverage = metrics$coverage,
    power = metrics$power,
    all_causal_covered = as.numeric(metrics$all_causal_covered),
    mean_causal_pip = mean(fit$pip[causal_variants]),
    minimum_causal_pip = min(fit$pip[causal_variants])
  )
}

fit_susie_rss_baseline <- function(
    marginal_effects, ld_matrix, true_effects) {
  elapsed <- system.time(fit <- susieR::susie_suff_stat(
    XtX = sample_size * ld_matrix,
    Xty = sample_size * marginal_effects,
    n = sample_size,
    yty = sample_size,
    L = fitted_number_of_effects,
    max_iter = maximum_iterations,
    tol = convergence_tolerance
  ))
  list(
    fit = fit,
    metrics = summarize_fitted_model(fit, ld_matrix, true_effects),
    seconds = unname(elapsed["elapsed"])
  )
}

fit_candidate_lambda <- function(
    marginal_effects, reference_model, true_effects) {
  pre <- list(
    r0 = as.numeric(crossprod(
      marginal_effects,
      reference_model$precision %*% marginal_effects
    )),
    logdetRbar = reference_model$logdet
  )
  profile_elapsed <- system.time(raw_fit <- fit_susie_collapse_r_shared_veb(
    x = marginal_effects,
    Rbar = reference_model$Rbar,
    N = sample_size,
    L = fitted_number_of_effects,
    sigma2 = 0.2,
    nu0_init = eta0_initial_value,
    nu0_bounds = eta0_bounds,
    nu0_coarse_grid_size = eta0_coarse_grid_size,
    estimate_sigma2 = TRUE,
    max_iter = maximum_iterations,
    tol = convergence_tolerance,
    verbose = FALSE,
    pre = pre
  ))
  eta0_used <- min(eta0_bounds[2L], eta0_multiplier * raw_fit$nu0)
  corrected_elapsed <- system.time(corrected_fit <-
    fit_susie_collapse_r_shared_fixed_nu0(
      x = marginal_effects,
      Rbar = reference_model$Rbar,
      N = sample_size,
      L = fitted_number_of_effects,
      sigma2 = 0.2,
      nu0 = eta0_used,
      estimate_sigma2 = TRUE,
      max_iter = maximum_iterations,
      tol = convergence_tolerance,
      verbose = FALSE,
      pre = pre
    ))
  list(
    raw_fit = raw_fit,
    raw_metrics = summarize_fitted_model(
      raw_fit, reference_model$Rbar, true_effects
    ),
    profile_seconds = unname(profile_elapsed["elapsed"]),
    corrected_fit = corrected_fit,
    corrected_metrics = summarize_fitted_model(
      corrected_fit, reference_model$Rbar, true_effects
    ),
    corrected_seconds = unname(corrected_elapsed["elapsed"]),
    eta0_used = eta0_used
  )
}

make_profile_row <- function(
    candidate, true_lambda, true_eta0, r0_seed, replicate, seed,
    causal, profile_type = c("profiled_eta0", "calibrated_eta0")) {
  profile_type <- match.arg(profile_type)
  corrected <- profile_type == "calibrated_eta0"
  fit <- if (corrected) candidate$corrected_fit else candidate$raw_fit
  metrics <- if (corrected) {
    candidate$corrected_metrics
  } else {
    candidate$raw_metrics
  }
  raw_eta0_hat <- candidate$raw_fit$nu0
  eta0_used <- if (corrected) candidate$eta0_used else raw_eta0_hat
  multiplier <- if (corrected) eta0_multiplier else 1
  data.frame(
    true_lambda = true_lambda,
    true_eta0 = true_eta0,
    r0_seed = r0_seed,
    replicate = replicate,
    seed = seed,
    profile_type = profile_type,
    candidate_lambda = candidate$lambda,
    raw_eta0_hat = raw_eta0_hat,
    eta0_used = eta0_used,
    eta0_multiplier = multiplier,
    raw_log_eta0_error = log(raw_eta0_hat) - log(true_eta0),
    used_log_eta0_error = log(eta0_used) - log(true_eta0),
    log_lambda_error = log(candidate$lambda) - log(true_lambda),
    eta0_boundary = candidate$raw_fit$nu0_boundary,
    lower_bound = fit$lower_bound,
    elbo_below_best = NA_real_,
    selected_L = as.integer(metrics["selected_L"]),
    coverage = metrics["coverage"],
    power = metrics["power"],
    all_causal_covered = as.logical(metrics["all_causal_covered"]),
    mean_causal_pip = metrics["mean_causal_pip"],
    minimum_causal_pip = metrics["minimum_causal_pip"],
    converged = fit$converged,
    iterations = length(fit$elbo),
    fit_seconds = if (corrected) {
      candidate$profile_seconds + candidate$corrected_seconds
    } else {
      candidate$profile_seconds
    },
    causal = causal,
    stringsAsFactors = FALSE
  )
}

make_result_row <- function(profile_row, fit_type, oracle_elbo) {
  data.frame(
    true_lambda = profile_row$true_lambda,
    true_eta0 = profile_row$true_eta0,
    r0_seed = profile_row$r0_seed,
    replicate = profile_row$replicate,
    seed = profile_row$seed,
    fit_type = fit_type,
    selected_lambda = profile_row$candidate_lambda,
    raw_eta0_hat = profile_row$raw_eta0_hat,
    eta0_used = profile_row$eta0_used,
    eta0_multiplier = profile_row$eta0_multiplier,
    raw_log_eta0_error = profile_row$raw_log_eta0_error,
    used_log_eta0_error = profile_row$used_log_eta0_error,
    log_lambda_error = profile_row$log_lambda_error,
    eta0_boundary = profile_row$eta0_boundary,
    lambda_boundary = if (profile_row$candidate_lambda == min(lambda_grid)) {
      "lower"
    } else if (profile_row$candidate_lambda == max(lambda_grid)) {
      "upper"
    } else {
      "interior"
    },
    lower_bound = profile_row$lower_bound,
    elbo_gain_over_matching_oracle = profile_row$lower_bound - oracle_elbo,
    selected_L = profile_row$selected_L,
    coverage = profile_row$coverage,
    power = profile_row$power,
    all_causal_covered = profile_row$all_causal_covered,
    mean_causal_pip = profile_row$mean_causal_pip,
    minimum_causal_pip = profile_row$minimum_causal_pip,
    converged = profile_row$converged,
    iterations = profile_row$iterations,
    fit_seconds = profile_row$fit_seconds,
    causal = profile_row$causal,
    stringsAsFactors = FALSE
  )
}

make_susie_result_row <- function(
    baseline, fit_type, true_lambda, true_eta0, r0_seed,
    replicate, seed, causal) {
  fit <- baseline$fit
  metrics <- baseline$metrics
  data.frame(
    true_lambda = true_lambda,
    true_eta0 = true_eta0,
    r0_seed = r0_seed,
    replicate = replicate,
    seed = seed,
    fit_type = fit_type,
    selected_lambda = NA_real_,
    raw_eta0_hat = NA_real_,
    eta0_used = NA_real_,
    eta0_multiplier = NA_real_,
    raw_log_eta0_error = NA_real_,
    used_log_eta0_error = NA_real_,
    log_lambda_error = NA_real_,
    eta0_boundary = NA_character_,
    lambda_boundary = NA_character_,
    lower_bound = tail(fit$elbo, 1L),
    elbo_gain_over_matching_oracle = NA_real_,
    selected_L = as.integer(metrics["selected_L"]),
    coverage = metrics["coverage"],
    power = metrics["power"],
    all_causal_covered = as.logical(metrics["all_causal_covered"]),
    mean_causal_pip = metrics["mean_causal_pip"],
    minimum_causal_pip = metrics["minimum_causal_pip"],
    converged = isTRUE(fit$converged),
    iterations = length(fit$elbo),
    fit_seconds = baseline$seconds,
    causal = causal,
    stringsAsFactors = FALSE
  )
}

total_datasets <- length(true_lambda_grid) * length(reference_panel_seeds) *
  length(true_eta0_grid) * number_of_replicates
dataset_index <- 0L
reference_ld_by_seed <- list()

for (r0_seed in reference_panel_seeds) {
  message("Preparing reference LD seed ", r0_seed)
  R0 <- simulate_coalescent_reference_ld(
    coalescent_sample_size = coalescent_sample_size,
    reference_sample_size = reference_sample_size,
    number_of_variants = number_of_variants,
    seed = r0_seed,
    heritability = heritability
  )
  reference_ld_by_seed[[as.character(r0_seed)]] <- R0
  reference_models <- prepare_reference_models(R0)

  for (true_lambda in true_lambda_grid) {
    true_model_index <- match(true_lambda, lambda_grid)
    true_reference_ld <- reference_models[[true_model_index]]$Rbar

    for (true_eta0 in true_eta0_grid) {
      for (replicate in seq_len(number_of_replicates)) {
        dataset_index <- dataset_index + 1L
        replicate_seed <- simulation_seed + 100000L +
          (match(r0_seed, reference_panel_seeds) - 1L) *
            length(true_eta0_grid) * number_of_replicates +
          (match(true_eta0, true_eta0_grid) - 1L) * number_of_replicates +
          replicate

        if (dataset_is_complete(
          true_lambda, true_eta0, r0_seed, replicate
        )) {
          message(sprintf(
            "[%d/%d] already complete: lambda=%g eta0=%g panel=%d rep=%d",
            dataset_index, total_datasets, true_lambda, true_eta0,
            r0_seed, replicate
          ))
          next
        }

        message(sprintf(
          "[%d/%d] true lambda=%g eta0=%g panel=%d rep=%d seed=%d",
          dataset_index, total_datasets, true_lambda, true_eta0,
          r0_seed, replicate, replicate_seed
        ))
        set.seed(replicate_seed)
        in_sample_ld <- sample_inverse_wishart(
          degrees_of_freedom = true_eta0 + number_of_variants + 1,
          scale_matrix = true_eta0 * true_reference_ld
        )
        in_sample_ld <- standardize_correlation(in_sample_ld)
        true_effects <- simulate_effects(
          in_sample_ld,
          heritability = heritability,
          true_number_of_effects = true_number_of_effects
        )
        marginal_effects <- simulate_marginal_effects(
          in_sample_ld,
          effects = true_effects,
          sample_size = sample_size
        )
        causal <- paste(which(true_effects != 0), collapse = ", ")

        susie_in_sample <- fit_susie_rss_baseline(
          marginal_effects, in_sample_ld, true_effects
        )
        susie_out_of_sample <- fit_susie_rss_baseline(
          marginal_effects, R0, true_effects
        )

        candidate_fits <- lapply(reference_models, function(model) {
          candidate <- fit_candidate_lambda(
            marginal_effects, model, true_effects
          )
          candidate$lambda <- model$lambda
          candidate
        })
        raw_profile <- do.call(rbind, lapply(candidate_fits, function(x) {
          make_profile_row(
            x, true_lambda, true_eta0, r0_seed, replicate,
            replicate_seed, causal, profile_type = "profiled_eta0"
          )
        }))
        corrected_profile <- do.call(
          rbind, lapply(candidate_fits, function(x) {
            make_profile_row(
              x, true_lambda, true_eta0, r0_seed, replicate,
              replicate_seed, causal, profile_type = "calibrated_eta0"
            )
          })
        )
        raw_best_index <- which.max(raw_profile$lower_bound)
        corrected_best_index <- which.max(corrected_profile$lower_bound)
        raw_oracle_index <- match(
          true_lambda, raw_profile$candidate_lambda
        )
        corrected_oracle_index <- match(
          true_lambda, corrected_profile$candidate_lambda
        )
        raw_profile$elbo_below_best <- raw_profile$lower_bound -
          raw_profile$lower_bound[raw_best_index]
        corrected_profile$elbo_below_best <-
          corrected_profile$lower_bound -
          corrected_profile$lower_bound[corrected_best_index]
        dataset_profile <- rbind(raw_profile, corrected_profile)

        key_rows <- if (nrow(profiles) == 0L) logical() else {
          profiles$true_lambda == true_lambda &
            profiles$true_eta0 == true_eta0 &
            profiles$r0_seed == r0_seed &
            profiles$replicate == replicate
        }
        if (length(key_rows) && any(key_rows)) profiles <- profiles[!key_rows, ]
        profiles <- rbind(profiles, dataset_profile)

        raw_oracle_elbo <- raw_profile$lower_bound[raw_oracle_index]
        corrected_oracle_elbo <-
          corrected_profile$lower_bound[corrected_oracle_index]
        dataset_results <- rbind(
          make_susie_result_row(
            susie_in_sample, "susie_rss_in_sample",
            true_lambda, true_eta0, r0_seed, replicate,
            replicate_seed, causal
          ),
          make_susie_result_row(
            susie_out_of_sample, "susie_rss_out_of_sample",
            true_lambda, true_eta0, r0_seed, replicate,
            replicate_seed, causal
          ),
          make_result_row(
            raw_profile[raw_best_index, , drop = FALSE],
            "joint_lambda_profiled_eta0", raw_oracle_elbo
          ),
          make_result_row(
            corrected_profile[corrected_best_index, , drop = FALSE],
            "joint_lambda_calibrated_eta0", corrected_oracle_elbo
          )
        )
        result_key_rows <- if (nrow(results) == 0L) logical() else {
          results$true_lambda == true_lambda &
            results$true_eta0 == true_eta0 &
            results$r0_seed == r0_seed &
            results$replicate == replicate
        }
        if (length(result_key_rows) && any(result_key_rows)) {
          results <- results[!result_key_rows, ]
        }
        results <- rbind(results, dataset_results)

        profiles <- profiles[order(
          profiles$true_lambda, profiles$r0_seed, profiles$true_eta0,
          profiles$replicate, profiles$profile_type,
          profiles$candidate_lambda
        ), , drop = FALSE]
        results <- results[order(
          results$true_lambda, results$r0_seed, results$true_eta0,
          results$replicate, results$fit_type
        ), , drop = FALSE]
        utils::write.csv(profiles, profiles_file, row.names = FALSE)
        utils::write.csv(results, results_file, row.names = FALSE)
      }
    }
  }
}

expected_datasets <- length(true_lambda_grid) *
  length(reference_panel_seeds) * length(true_eta0_grid) *
  number_of_replicates
stopifnot(
  nrow(profiles) == expected_datasets * length(lambda_grid) * 2L,
  nrow(results) == expected_datasets * 4L
)

single_effect_credible_set <- function(alpha, coverage = 0.95) {
  posterior_probabilities <- as.numeric(alpha[1L, ])
  ordered_variants <- order(posterior_probabilities, decreasing = TRUE)
  cutoff <- which(
    cumsum(posterior_probabilities[ordered_variants]) >= coverage
  )[1L]
  ordered_variants[seq_len(cutoff)]
}

fit_ser_fallback <- function(marginal_effects) {
  elapsed <- system.time(fit <- susieR::susie_suff_stat(
    XtX = sample_size * diag(1, length(marginal_effects)),
    Xty = sample_size * marginal_effects,
    n = sample_size,
    yty = sample_size,
    L = 1L,
    scaled_prior_variance = 0.2,
    estimate_prior_variance = TRUE,
    estimate_residual_variance = FALSE
  ))
  list(
    fit = fit,
    credible_set = single_effect_credible_set(fit$alpha, coverage),
    seconds = unname(elapsed["elapsed"])
  )
}

parse_causal_variants <- function(causal) {
  as.integer(strsplit(as.character(causal), ", ", fixed = TRUE)[[1L]])
}

reconstruct_ser_fallback <- function(dataset_row) {
  R0 <- reference_ld_by_seed[[as.character(dataset_row$r0_seed)]]
  true_reference_ld <- make_stabilized_reference_ld(
    R0, lambda = dataset_row$true_lambda
  )
  set.seed(as.integer(dataset_row$seed))
  in_sample_ld <- sample_inverse_wishart(
    degrees_of_freedom = dataset_row$true_eta0 +
      number_of_variants + 1,
    scale_matrix = dataset_row$true_eta0 * true_reference_ld
  )
  in_sample_ld <- standardize_correlation(in_sample_ld)
  true_effects <- simulate_effects(
    in_sample_ld,
    heritability = heritability,
    true_number_of_effects = true_number_of_effects
  )
  marginal_effects <- simulate_marginal_effects(
    in_sample_ld,
    effects = true_effects,
    sample_size = sample_size
  )
  causal_variants <- which(true_effects != 0)
  stopifnot(identical(
    causal_variants, parse_causal_variants(dataset_row$causal)
  ))
  ser <- fit_ser_fallback(marginal_effects)
  data.frame(
    true_lambda = dataset_row$true_lambda,
    true_eta0 = dataset_row$true_eta0,
    r0_seed = dataset_row$r0_seed,
    replicate = dataset_row$replicate,
    seed = dataset_row$seed,
    ser_coverage = as.numeric(any(
      ser$credible_set %in% causal_variants
    )),
    ser_power = mean(causal_variants %in% ser$credible_set),
    ser_all_causal_covered = all(
      causal_variants %in% ser$credible_set
    ),
    ser_mean_causal_pip = mean(ser$fit$pip[causal_variants]),
    ser_minimum_causal_pip = min(ser$fit$pip[causal_variants]),
    ser_credible_set_size = length(ser$credible_set),
    ser_fit_seconds = ser$seconds,
    stringsAsFactors = FALSE
  )
}

apply_ser_fallback <- function(results) {
  dataset_columns <- c(
    "true_lambda", "true_eta0", "r0_seed", "replicate", "seed", "causal"
  )
  datasets <- unique(results[, dataset_columns, drop = FALSE])
  fallback_metrics <- do.call(
    rbind, lapply(seq_len(nrow(datasets)), function(index) {
      reconstruct_ser_fallback(datasets[index, , drop = FALSE])
    })
  )
  make_key <- function(data) paste(
    data$true_lambda, data$true_eta0, data$r0_seed, data$replicate,
    data$seed, sep = "|"
  )
  fallback_index <- match(
    make_key(results), make_key(fallback_metrics)
  )
  stopifnot(!anyNA(fallback_index))

  postprocessed <- results
  postprocessed$raw_selected_L <- results$selected_L
  postprocessed$raw_coverage <- results$coverage
  postprocessed$raw_power <- results$power
  postprocessed$raw_all_causal_covered <- results$all_causal_covered
  postprocessed$raw_mean_causal_pip <- results$mean_causal_pip
  postprocessed$raw_minimum_causal_pip <- results$minimum_causal_pip
  proposed_method <- results$fit_type %in% c(
    "joint_lambda_profiled_eta0", "joint_lambda_calibrated_eta0"
  )
  postprocessed$fallback_used <-
    proposed_method & results$selected_L <= 1L
  postprocessed$ser_credible_set_size <-
    fallback_metrics$ser_credible_set_size[fallback_index]
  postprocessed$ser_fit_seconds <-
    fallback_metrics$ser_fit_seconds[fallback_index]

  rows <- postprocessed$fallback_used
  matched <- fallback_index[rows]
  postprocessed$selected_L[rows] <- 1L
  postprocessed$coverage[rows] <- fallback_metrics$ser_coverage[matched]
  postprocessed$power[rows] <- fallback_metrics$ser_power[matched]
  postprocessed$all_causal_covered[rows] <-
    fallback_metrics$ser_all_causal_covered[matched]
  postprocessed$mean_causal_pip[rows] <-
    fallback_metrics$ser_mean_causal_pip[matched]
  postprocessed$minimum_causal_pip[rows] <-
    fallback_metrics$ser_minimum_causal_pip[matched]
  postprocessed
}

median_or_missing <- function(x) {
  if (all(is.na(x))) NA_real_ else stats::median(x, na.rm = TRUE)
}

quantile_or_missing <- function(x, probability) {
  if (all(is.na(x))) {
    NA_real_
  } else {
    unname(stats::quantile(x, probability, na.rm = TRUE))
  }
}

mean_or_missing <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

summarize_results <- function(results) {
  groups <- split(
    results,
    interaction(
      results$true_lambda, results$true_eta0, results$fit_type,
      drop = TRUE
    )
  )
  rows <- lapply(groups, function(group) {
    data.frame(
      true_lambda = group$true_lambda[1L],
      true_eta0 = group$true_eta0[1L],
      fit_type = group$fit_type[1L],
      fits = nrow(group),
      median_lambda_hat = median_or_missing(group$selected_lambda),
      lambda_exact_frequency = mean_or_missing(
        group$selected_lambda == group$true_lambda
      ),
      median_raw_eta0_hat = median_or_missing(group$raw_eta0_hat),
      raw_eta0_q10 = quantile_or_missing(group$raw_eta0_hat, 0.1),
      raw_eta0_q90 = quantile_or_missing(group$raw_eta0_hat, 0.9),
      median_eta0_used = median_or_missing(group$eta0_used),
      eta0_used_q10 = quantile_or_missing(group$eta0_used, 0.1),
      eta0_used_q90 = quantile_or_missing(group$eta0_used, 0.9),
      median_raw_log_eta0_error = median_or_missing(
        group$raw_log_eta0_error
      ),
      median_used_log_eta0_error = median_or_missing(
        group$used_log_eta0_error
      ),
      median_log_lambda_error = median_or_missing(group$log_lambda_error),
      raw_eta0_underestimate_frequency = mean_or_missing(
        group$raw_eta0_hat < group$true_eta0
      ),
      used_eta0_underestimate_frequency = mean_or_missing(
        group$eta0_used < group$true_eta0
      ),
      eta0_boundary_frequency = mean_or_missing(
        group$eta0_boundary != "interior"
      ),
      lambda_boundary_frequency = mean_or_missing(
        group$lambda_boundary != "interior"
      ),
      mean_elbo_gain_over_matching_oracle = mean_or_missing(
        group$elbo_gain_over_matching_oracle
      ),
      median_elbo_gain_over_matching_oracle = median_or_missing(
        group$elbo_gain_over_matching_oracle
      ),
      mean_selected_L = mean(group$selected_L),
      mean_raw_selected_L = if ("raw_selected_L" %in% names(group)) {
        mean(group$raw_selected_L)
      } else {
        mean(group$selected_L)
      },
      fallback_frequency = if ("fallback_used" %in% names(group)) {
        mean(group$fallback_used)
      } else {
        0
      },
      mean_coverage = mean(group$coverage),
      mean_power = mean(group$power),
      all_causal_covered_frequency = mean(group$all_causal_covered),
      mean_causal_pip = mean(group$mean_causal_pip),
      convergence_frequency = mean(group$converged),
      median_iterations = stats::median(group$iterations),
      median_fit_seconds = stats::median(group$fit_seconds),
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, rows)
  rownames(summary) <- NULL
  summary[order(
    summary$true_lambda, summary$true_eta0, summary$fit_type
  ), , drop = FALSE]
}

raw_summary <- summarize_results(results)
utils::write.csv(raw_summary, raw_summary_file, row.names = FALSE)
postprocessed_results <- apply_ser_fallback(results)
stopifnot(
  nrow(postprocessed_results) == nrow(results),
  all(postprocessed_results$selected_L[
    postprocessed_results$fit_type %in% c(
      "joint_lambda_profiled_eta0", "joint_lambda_calibrated_eta0"
    )
  ] >= 1L),
  all(
    postprocessed_results$fallback_used ==
      (
        postprocessed_results$fit_type %in% c(
          "joint_lambda_profiled_eta0", "joint_lambda_calibrated_eta0"
        ) & postprocessed_results$raw_selected_L <= 1L
      )
  )
)
utils::write.csv(
  postprocessed_results, postprocessed_results_file, row.names = FALSE
)
summary <- summarize_results(postprocessed_results)
utils::write.csv(summary, summary_file, row.names = FALSE)

summarize_eta0_calibration <- function(data, stratum) {
  fit <- stats::lm(log(true_eta0) ~ log(raw_eta0_hat), data = data)
  coefficients <- stats::coef(fit)
  fixed_slope_log_shift <- mean(
    log(data$true_eta0) - log(data$raw_eta0_hat)
  )
  data.frame(
    stratum = stratum,
    fits = nrow(data),
    intercept = unname(coefficients[1L]),
    exp_intercept = exp(unname(coefficients[1L])),
    slope = unname(coefficients[2L]),
    r_squared = base::summary(fit)$r.squared,
    fixed_slope_one_log_shift = fixed_slope_log_shift,
    fixed_slope_one_multiplier = exp(fixed_slope_log_shift),
    stringsAsFactors = FALSE
  )
}

eta0_calibration_data <- results[
  results$fit_type == "joint_lambda_profiled_eta0", , drop = FALSE
]
eta0_calibration <- rbind(
  summarize_eta0_calibration(eta0_calibration_data, "pooled"),
  do.call(rbind, lapply(sort(unique(eta0_calibration_data$true_lambda)),
    function(lambda) {
      summarize_eta0_calibration(
        eta0_calibration_data[
          eta0_calibration_data$true_lambda == lambda, , drop = FALSE
        ],
        paste0("true_lambda=", format(lambda, scientific = FALSE))
      )
    }
  ))
)
utils::write.csv(
  eta0_calibration, eta0_calibration_file, row.names = FALSE
)
settings <- list(
  lambda_grid = lambda_grid,
  true_lambda_grid = true_lambda_grid,
  true_eta0_grid = true_eta0_grid,
  reference_panel_seeds = reference_panel_seeds,
  number_of_replicates = number_of_replicates,
  number_of_variants = number_of_variants,
  sample_size = sample_size,
  coalescent_sample_size = coalescent_sample_size,
  reference_sample_size = reference_sample_size,
  heritability = heritability,
  true_number_of_effects = true_number_of_effects,
  fitted_number_of_effects = fitted_number_of_effects,
  coverage = coverage,
  minimum_absolute_correlation = minimum_absolute_correlation,
  eta0_bounds = eta0_bounds,
  maximum_iterations = maximum_iterations,
  convergence_tolerance = convergence_tolerance,
  eta0_multiplier = eta0_multiplier,
  simulation_seed = simulation_seed
)
saveRDS(
  list(
    results = postprocessed_results,
    raw_results = results,
    profiles = profiles,
    summary = summary,
    raw_summary = raw_summary,
    eta0_calibration = eta0_calibration,
    settings = settings,
    reporting_rule = paste(
      "For the two joint collapsed-R methods, use SER-RSS without an",
      "LD-purity filter when raw purity-filtered L is zero or one"
    )
  ),
  output_file
)

draw_eta0_recovery <- function() {
  colors <- c(
    joint_lambda_profiled_eta0 = "#762A83",
    joint_lambda_calibrated_eta0 = "#1B7837"
  )
  labels <- c(
    joint_lambda_profiled_eta0 = expression(paste("profiled ", eta[0])),
    joint_lambda_calibrated_eta0 = expression(paste("1.537", eta[0]))
  )
  old_par <- graphics::par(mfrow = c(1, length(true_lambda_grid)), mar = c(4.5, 4.6, 3, 1))
  on.exit(graphics::par(old_par))
  for (true_lambda in true_lambda_grid) {
    panel <- summary[summary$true_lambda == true_lambda, ]
    limits <- range(c(panel$true_eta0, panel$median_eta0_used), finite = TRUE)
    graphics::plot(
      NA, xlim = range(panel$true_eta0), ylim = limits, log = "xy",
      xlab = expression(paste("true ", eta[0])),
      ylab = expression(paste("median ", eta[0], " used")),
      main = bquote(lambda[true] == .(true_lambda))
    )
    graphics::grid()
    graphics::abline(0, 1, lty = 2, col = "gray45")
    for (fit_type in names(colors)) {
      rows <- panel$fit_type == fit_type
      graphics::lines(
        panel$true_eta0[rows], panel$median_eta0_used[rows],
        type = "b", pch = match(fit_type, names(colors)) + 15L,
        lwd = 2, col = colors[fit_type]
      )
    }
    graphics::legend(
      "topleft", legend = labels, col = colors, pch = 16:17,
      lty = 1, lwd = 2, bty = "n"
    )
  }
}

draw_lambda_recovery <- function() {
  old_par <- graphics::par(mfrow = c(1, length(true_lambda_grid)), mar = c(4.5, 4.6, 3, 1))
  on.exit(graphics::par(old_par))
  joint_types <- c(
    "joint_lambda_profiled_eta0", "joint_lambda_calibrated_eta0"
  )
  joint_colors <- c(
    joint_lambda_profiled_eta0 = "#762A83",
    joint_lambda_calibrated_eta0 = "#1B7837"
  )
  joint <- results[results$fit_type %in% joint_types, ]
  for (true_lambda in true_lambda_grid) {
    panel <- joint[joint$true_lambda == true_lambda, ]
    graphics::plot(
      NA, xlim = range(panel$true_eta0), ylim = range(lambda_grid),
      log = "xy",
      xlab = expression(paste("true ", eta[0])),
      ylab = expression(paste("selected ", lambda)),
      main = bquote(lambda[true] == .(true_lambda)),
      yaxs = "i"
    )
    graphics::grid()
    graphics::abline(h = true_lambda, lty = 2, col = "gray35")
    for (fit_type in joint_types) {
      rows <- panel[panel$fit_type == fit_type, ]
      graphics::points(
        jitter(rows$true_eta0, factor = 0.04), rows$selected_lambda,
        pch = if (fit_type == joint_types[1L]) 16 else 17,
        col = grDevices::adjustcolor(
          joint_colors[fit_type], alpha.f = 0.25
        )
      )
      medians <- stats::aggregate(
        selected_lambda ~ true_eta0, rows, stats::median
      )
      graphics::lines(
        medians$true_eta0, medians$selected_lambda,
        type = "b", pch = if (fit_type == joint_types[1L]) 16 else 17,
        lwd = 2.5, col = joint_colors[fit_type]
      )
    }
    graphics::legend(
      "topleft", c("profiled eta0", "1.537x eta0"),
      col = joint_colors, pch = c(16, 17), lty = 1, lwd = 2,
      bty = "n"
    )
  }
}

draw_elbo_profile <- function() {
  eta_examples <- unique(c(
    min(true_eta0_grid),
    true_eta0_grid[ceiling(length(true_eta0_grid) / 2)],
    max(true_eta0_grid)
  ))
  old_par <- graphics::par(
    mfrow = c(length(true_lambda_grid), length(eta_examples)),
    mar = c(4.2, 4.3, 2.8, 1)
  )
  on.exit(graphics::par(old_par))
  for (true_lambda in true_lambda_grid) {
    for (true_eta0 in eta_examples) {
      panel <- profiles[
        profiles$true_lambda == true_lambda &
          profiles$true_eta0 == true_eta0,
      ]
      raw <- panel[panel$profile_type == "profiled_eta0", ]
      corrected <- panel[panel$profile_type == "calibrated_eta0", ]
      aggregate_raw <- stats::aggregate(
        cbind(elbo_below_best, eta0_used) ~ candidate_lambda,
        raw, stats::median
      )
      aggregate_corrected <- stats::aggregate(
        cbind(elbo_below_best, eta0_used) ~ candidate_lambda,
        corrected, stats::median
      )
      graphics::plot(
        aggregate_raw$candidate_lambda,
        aggregate_raw$elbo_below_best,
        type = "b", pch = 16, log = "x", col = "#762A83", lwd = 2,
        xlab = expression(lambda),
        ylab = "Median ELBO - best ELBO",
        main = bquote(paste(
          lambda[true] == .(true_lambda), ", ", eta[0] == .(true_eta0)
        ))
      )
      graphics::grid()
      graphics::abline(v = true_lambda, lty = 2, col = "gray35")
      graphics::lines(
        aggregate_corrected$candidate_lambda,
        aggregate_corrected$elbo_below_best,
        type = "b", pch = 17, col = "#1B7837", lwd = 2
      )
      graphics::legend(
        "bottomleft", c("profiled eta0", "1.537x eta0"),
        col = c("#762A83", "#1B7837"), pch = c(16, 17),
        lty = 1, lwd = 2, bty = "n"
      )
    }
  }
}

draw_selected_L <- function() {
  colors <- c(
    susie_rss_in_sample = "#2166AC",
    susie_rss_out_of_sample = "#E66101",
    joint_lambda_profiled_eta0 = "#762A83",
    joint_lambda_calibrated_eta0 = "#1B7837"
  )
  labels <- c(
    "SuSiE-RSS: in-sample R", "SuSiE-RSS: out-of-sample R0",
    "Collapsed-R: profiled eta0", "Collapsed-R: 1.537x eta0"
  )
  old_par <- graphics::par(
    mfrow = c(1, length(true_lambda_grid)), mar = c(4.5, 4.6, 3, 1)
  )
  on.exit(graphics::par(old_par))
  for (true_lambda in true_lambda_grid) {
    panel <- summary[summary$true_lambda == true_lambda, ]
    graphics::plot(
      NA, xlim = range(panel$true_eta0), ylim = c(0, fitted_number_of_effects),
      log = "x", xlab = expression(paste("true ", eta[0])),
      ylab = "Mean reported L after SER fallback",
      main = bquote(lambda[true] == .(true_lambda))
    )
    graphics::grid()
    graphics::abline(h = true_number_of_effects, lty = 2, col = "gray35")
    for (fit_type in names(colors)) {
      rows <- panel$fit_type == fit_type
      graphics::lines(
        panel$true_eta0[rows], panel$mean_selected_L[rows],
        type = "b", pch = match(fit_type, names(colors)) + 14L,
        lwd = 2, col = colors[fit_type]
      )
    }
    graphics::legend(
      "topleft", labels, col = colors, pch = 15:18,
      lty = 1, lwd = 2, bty = "n", cex = 0.8
    )
  }
}

draw_summary_metrics <- function() {
  colors <- c(
    susie_rss_in_sample = "#2166AC",
    susie_rss_out_of_sample = "#E66101",
    joint_lambda_profiled_eta0 = "#762A83",
    joint_lambda_calibrated_eta0 = "#1B7837"
  )
  old_par <- graphics::par(
    mfrow = c(length(true_lambda_grid), 2L), mar = c(4.3, 4.5, 3, 1)
  )
  on.exit(graphics::par(old_par))
  for (true_lambda in true_lambda_grid) {
    panel <- summary[summary$true_lambda == true_lambda, ]
    for (metric in c("mean_coverage", "mean_power")) {
      graphics::plot(
        NA, xlim = range(panel$true_eta0), ylim = c(0, 1), log = "x",
        xlab = expression(paste("true ", eta[0])),
        ylab = if (metric == "mean_coverage") "CS coverage" else "Causal power",
        main = bquote(lambda[true] == .(true_lambda))
      )
      graphics::grid()
      if (metric == "mean_coverage") {
        graphics::abline(h = coverage, lty = 2, col = "gray35")
      }
      for (fit_type in names(colors)) {
        rows <- panel$fit_type == fit_type
        graphics::lines(
          panel$true_eta0[rows], panel[[metric]][rows],
          type = "b", pch = match(fit_type, names(colors)) + 14L,
          lwd = 2, col = colors[fit_type]
        )
      }
      if (true_lambda == true_lambda_grid[1L] &&
          metric == "mean_coverage") {
        graphics::legend(
          "bottomright",
          c(
            "SuSiE-RSS: in-sample R", "SuSiE-RSS: out-of-sample R0",
            "Collapsed-R: profiled eta0", "Collapsed-R: 1.537x eta0"
          ),
          col = colors, pch = 15:18, lty = 1, lwd = 2,
          bty = "n", cex = 0.75
        )
      }
    }
  }
}

save_figure <- function(name, draw, width, height) {
  grDevices::png(
    paste0(figure_prefix, "-", name, ".png"),
    width = width * 150, height = height * 150, res = 150
  )
  draw()
  grDevices::dev.off()
  grDevices::pdf(
    paste0(figure_prefix, "-", name, ".pdf"),
    width = width, height = height
  )
  draw()
  grDevices::dev.off()
}

save_figure("eta0-recovery", draw_eta0_recovery, 11, 5.5)
save_figure("lambda-recovery", draw_lambda_recovery, 11, 5.5)
save_figure("elbo-profile", draw_elbo_profile, 12, 7.5)
save_figure("selected-L", draw_selected_L, 11, 5.5)
save_figure("summary-metrics", draw_summary_metrics, 11, 9)

print(summary)
