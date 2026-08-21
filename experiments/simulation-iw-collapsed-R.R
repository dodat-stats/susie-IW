#!/usr/bin/env Rscript

# Inverse-Wishart simulation comparing the collapsed-R model
# with SuSiE-RSS using in-sample and out-of-sample LD.

library(susieR)

find_susie_iw_root <- function(path = getwd()) {
  path <- normalizePath(path, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "code", "01-collapsed-R.R")) &&
        file.exists(file.path(path, "code", "sim_2pop.R"))) {
      return(path)
    }
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not find the susie-IW repository root")
    }
    path <- parent
  }
}

susie_iw_root <- find_susie_iw_root()
source(file.path(susie_iw_root, "code", "01-collapsed-R.R"))
source(file.path(susie_iw_root, "code", "sim_2pop.R"))

standardize_correlation <- function(R) {
  R <- (R + t(R)) / 2
  R <- stats::cov2cor(R)
  R <- (R + t(R)) / 2
  diag(R) <- 1
  R
}

make_stabilized_reference_ld <- function(R0, lambda = 1e-3) {
  J <- ncol(R0)
  standardize_correlation((1 - lambda) * R0 + lambda * diag(1, J))
}

stable_cholesky <- function(A, ridge = 1e-10, maximum_attempts = 8L) {
  A <- (A + t(A)) / 2
  J <- nrow(A)
  for (attempt in seq_len(maximum_attempts)) {
    factor <- try(chol(A + ridge * diag(1, J)), silent = TRUE)
    if (!inherits(factor, "try-error")) {
      return(factor)
    }
    ridge <- ridge * 10
  }
  stop("Cholesky factorization failed after adding a ridge")
}

sample_inverse_wishart <- function(degrees_of_freedom, scale_matrix) {
  scale_matrix <- (scale_matrix + t(scale_matrix)) / 2
  precision <- stats::rWishart(
    1,
    df = degrees_of_freedom,
    Sigma = solve(scale_matrix)
  )[, , 1]
  solve((precision + t(precision)) / 2)
}

simulate_coalescent_reference_ld <- function(
    coalescent_sample_size, reference_sample_size, number_of_variants,
    seed, heritability = 0.01, chromosome = "chr22",
    start = 20e6, end = 21e6, split_time = 10L) {
  message(sprintf(
    paste0(
      "Simulating coalescent reference panel: coalescent N=%d, ",
      "reference N=%d, J=%d, seed=%d"
    ),
    coalescent_sample_size, reference_sample_size, number_of_variants, seed
  ))
  simulation <- sim_2pop(
    chrom = chromosome,
    start = start,
    end = end,
    T_split = split_time,
    n_gwas = coalescent_sample_size,
    n_ref = reference_sample_size,
    J = number_of_variants,
    h2 = heritability,
    seed = seed
  )
  standardize_correlation(simulation$R0)
}

choose_causal_variants <- function(in_sample_ld, true_number_of_effects = 2L) {
  number_of_variants <- ncol(in_sample_ld)
  true_number_of_effects <- as.integer(true_number_of_effects)
  stopifnot(
    true_number_of_effects >= 1L,
    true_number_of_effects <= number_of_variants
  )
  sort(sample.int(number_of_variants, true_number_of_effects))
}

simulate_effects <- function(in_sample_ld, heritability = 0.01,
                             true_number_of_effects = 2L,
                             causal_variants = NULL) {
  number_of_variants <- ncol(in_sample_ld)
  if (is.null(causal_variants)) {
    causal_variants <- choose_causal_variants(
      in_sample_ld,
      true_number_of_effects = true_number_of_effects
    )
  }

  effects <- numeric(number_of_variants)
  effect_magnitudes <- if (length(causal_variants) == 2L) {
    c(
      stats::runif(1L, min = 1, max = 2),
      stats::runif(1L, min = 2.5, max = 3)
    )
  } else {
    stats::runif(length(causal_variants), min = 1, max = 2)
  }
  effect_signs <- sample(c(-1, 1), length(causal_variants), replace = TRUE)
  effects[causal_variants] <- effect_signs * effect_magnitudes

  explained_variance <- as.numeric(crossprod(
    effects,
    in_sample_ld %*% effects
  ))
  if (!is.finite(explained_variance) || explained_variance <= 0) {
    stop("Could not scale effects because their explained variance is invalid")
  }
  effects * sqrt(heritability / explained_variance)
}

simulate_marginal_effects <- function(in_sample_ld, effects, sample_size) {
  mean_marginal_effects <- as.numeric(in_sample_ld %*% effects)
  noise <- as.numeric(
    t(stable_cholesky(in_sample_ld)) %*%
      stats::rnorm(ncol(in_sample_ld))
  ) / sqrt(sample_size)
  mean_marginal_effects + noise
}

make_credible_sets_from_alpha <- function(
    alpha, ld_for_purity, coverage = 0.95,
    minimum_absolute_correlation = 0.5,
    minimum_top_alpha = 1e-4) {
  alpha <- as.matrix(alpha)
  credible_sets <- list()
  purity <- numeric()

  for (effect_index in seq_len(nrow(alpha))) {
    posterior_weights <- as.numeric(alpha[effect_index, ])
    if (max(posterior_weights) < minimum_top_alpha) {
      next
    }

    ordered_variants <- order(posterior_weights, decreasing = TRUE)
    cutoff <- which(
      cumsum(posterior_weights[ordered_variants]) >= coverage
    )[1]
    variants <- ordered_variants[seq_len(cutoff)]
    minimum_correlation <- if (length(variants) == 1L) {
      1
    } else {
      absolute_ld <- abs(ld_for_purity[variants, variants, drop = FALSE])
      min(absolute_ld[upper.tri(absolute_ld)])
    }

    if (is.finite(minimum_correlation) &&
        minimum_correlation >= minimum_absolute_correlation) {
      name <- paste0("L", effect_index)
      credible_sets[[name]] <- variants
      purity[name] <- minimum_correlation
    }
  }

  list(cs = credible_sets, purity = purity)
}

credible_set_metrics_from_list <- function(credible_sets, true_effects) {
  causal_variants <- which(true_effects != 0)
  number_of_credible_sets <- length(credible_sets)
  contains_causal <- if (number_of_credible_sets == 0L) {
    logical()
  } else {
    vapply(
      credible_sets,
      function(variants) any(variants %in% causal_variants),
      logical(1)
    )
  }
  credible_set_union <- if (number_of_credible_sets == 0L) {
    integer()
  } else {
    unique(unlist(credible_sets))
  }

  data.frame(
    selected_L = number_of_credible_sets,
    coverage = if (number_of_credible_sets == 0L) {
      0
    } else {
      mean(contains_causal)
    },
    power = mean(causal_variants %in% credible_set_union),
    all_causal_covered = all(causal_variants %in% credible_set_union),
    n_causal = length(causal_variants)
  )
}

credible_set_metrics_from_alpha <- function(
    alpha, ld_for_purity, true_effects, coverage = 0.95,
    minimum_absolute_correlation = 0.5) {
  credible_sets <- make_credible_sets_from_alpha(
    alpha = alpha,
    ld_for_purity = ld_for_purity,
    coverage = coverage,
    minimum_absolute_correlation = minimum_absolute_correlation
  )$cs
  credible_set_metrics_from_list(credible_sets, true_effects)
}

summarize_standard_fit <- function(
    method, fit, ld_for_purity, true_effects, seed,
    replicate_index, true_eta0, fit_seconds,
    coverage = 0.95, minimum_absolute_correlation = 0.5) {
  metrics <- credible_set_metrics_from_alpha(
    alpha = fit$alpha,
    ld_for_purity = ld_for_purity,
    true_effects = true_effects,
    coverage = coverage,
    minimum_absolute_correlation = minimum_absolute_correlation
  )

  transform(
    metrics,
    method = method,
    seed = seed,
    rep = replicate_index,
    nu0_true = true_eta0,
    lambda = NA_real_,
    nu0_hat = NA_real_,
    final_elbo = NA_real_,
    fit_seconds = fit_seconds,
    causal = paste(which(true_effects != 0), collapse = ", "),
    causal_pip = paste(
      sprintf("%.3f", fit$pip[true_effects != 0]),
      collapse = ", "
    ),
    top = paste(
      apply(as.matrix(fit$alpha), 1, which.max),
      collapse = ", "
    ),
    top_alpha = paste(
      sprintf("%.3f", apply(as.matrix(fit$alpha), 1, max)),
      collapse = ", "
    )
  )
}

summarize_collapsed_fit <- function(
    fit, reference_ld, true_effects, seed, replicate_index,
    true_eta0, lambda, fit_seconds, coverage = 0.95,
    minimum_absolute_correlation = 0.5) {
  credible_set_object <- susieR::susie_get_cs(
    fit,
    Xcorr = reference_ld,
    coverage = coverage,
    min_abs_corr = minimum_absolute_correlation,
    dedup = TRUE
  )
  metrics <- credible_set_metrics_from_list(
    credible_set_object$cs,
    true_effects
  )

  list(
    summary = transform(
      metrics,
      method = "collapsed_r",
      seed = seed,
      rep = replicate_index,
      nu0_true = true_eta0,
      lambda = lambda,
      nu0_hat = fit$nu0,
      final_elbo = fit$lower_bound,
      fit_seconds = fit_seconds,
      causal = paste(which(true_effects != 0), collapse = ", "),
      causal_pip = paste(
        sprintf("%.3f", fit$pip[true_effects != 0]),
        collapse = ", "
      ),
      top = paste(apply(as.matrix(fit$alpha), 1, which.max), collapse = ", "),
      top_alpha = paste(
        sprintf("%.3f", apply(as.matrix(fit$alpha), 1, max)),
        collapse = ", "
      )
    ),
    credible_sets = credible_set_object
  )
}

simulate_one_collapsed_replicate <- function(
    true_eta0, replicate_index, seed, R0, reference_ld,
    sample_size, heritability, true_number_of_effects,
    fitted_number_of_effects, fitting_eta0_grid,
    coverage, minimum_absolute_correlation, lambda,
    collapsed_maximum_iterations = 500L,
    collapsed_tolerance = 1e-7, verbose = TRUE) {
  set.seed(seed)
  number_of_variants <- ncol(reference_ld)

  if (verbose) {
    message(sprintf(
      "true eta0=%g rep=%d seed=%d: simulate and fit",
      true_eta0, replicate_index, seed
    ))
  }

  inverse_wishart_time <- system.time({
    in_sample_ld <- sample_inverse_wishart(
      degrees_of_freedom = true_eta0 + number_of_variants + 1,
      scale_matrix = true_eta0 * reference_ld
    )
    in_sample_ld <- standardize_correlation(in_sample_ld)
  })
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

  in_sample_time <- system.time({
    in_sample_fit <- susieR::susie_suff_stat(
      XtX = sample_size * in_sample_ld,
      Xty = sample_size * marginal_effects,
      n = sample_size,
      yty = sample_size,
      L = fitted_number_of_effects
    )
  })
  out_of_sample_time <- system.time({
    out_of_sample_fit <- susieR::susie_suff_stat(
      XtX = sample_size * R0,
      Xty = sample_size * marginal_effects,
      n = sample_size,
      yty = sample_size,
      L = fitted_number_of_effects
    )
  })
  collapsed_time <- system.time({
    collapsed_fit <- fit_susie_collapse_r_shared(
      x = marginal_effects,
      Rbar = reference_ld,
      N = sample_size,
      L = fitted_number_of_effects,
      sigma2 = 0.2,
      nu0_grid = fitting_eta0_grid,
      estimate_sigma2 = TRUE,
      max_iter = collapsed_maximum_iterations,
      tol = collapsed_tolerance,
      verbose = FALSE
    )
  })

  collapsed_summary <- summarize_collapsed_fit(
    fit = collapsed_fit,
    reference_ld = reference_ld,
    true_effects = true_effects,
    seed = seed,
    replicate_index = replicate_index,
    true_eta0 = true_eta0,
    lambda = lambda,
    fit_seconds = unname(collapsed_time["elapsed"]),
    coverage = coverage,
    minimum_absolute_correlation = minimum_absolute_correlation
  )

  results <- rbind(
    summarize_standard_fit(
      method = "susie_insample_R",
      fit = in_sample_fit,
      ld_for_purity = in_sample_ld,
      true_effects = true_effects,
      seed = seed,
      replicate_index = replicate_index,
      true_eta0 = true_eta0,
      fit_seconds = unname(in_sample_time["elapsed"]),
      coverage = coverage,
      minimum_absolute_correlation = minimum_absolute_correlation
    ),
    summarize_standard_fit(
      method = "susie_outsample_R0",
      fit = out_of_sample_fit,
      ld_for_purity = R0,
      true_effects = true_effects,
      seed = seed,
      replicate_index = replicate_index,
      true_eta0 = true_eta0,
      fit_seconds = unname(out_of_sample_time["elapsed"]),
      coverage = coverage,
      minimum_absolute_correlation = minimum_absolute_correlation
    ),
    collapsed_summary$summary
  )

  results$purity_filtered_L <- NA_integer_
  results$purity_filtered_L[results$method == "collapsed_r"] <-
    collapsed_summary$summary$selected_L
  results$eta0_grid_boundary <- NA_character_
  collapsed_boundary <- if (collapsed_fit$nu0 == min(fitting_eta0_grid)) {
    "lower"
  } else if (collapsed_fit$nu0 == max(fitting_eta0_grid)) {
    "upper"
  } else {
    "interior"
  }
  results$eta0_grid_boundary[results$method == "collapsed_r"] <-
    collapsed_boundary
  results$estimated_prior_variances <- NA_character_
  results$estimated_prior_variances[results$method == "collapsed_r"] <-
    paste(format(collapsed_fit$sigma2, digits = 5), collapse = ", ")
  results$iw_seconds <- unname(inverse_wishart_time["elapsed"])
  results$N <- sample_size
  results$p <- number_of_variants
  results$L_true <- true_number_of_effects
  results$L_fit <- fitted_number_of_effects
  results$h2 <- heritability

  profile <- data.frame(
    eta0 = collapsed_fit$nu0_grid,
    lower_bound = collapsed_fit$grid_lower_bound,
    seed = seed,
    rep = replicate_index,
    nu0_true = true_eta0
  )

  list(results = results, profile = profile)
}

mean_or_missing <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

median_or_missing <- function(x) {
  if (all(is.na(x))) NA_real_ else stats::median(x, na.rm = TRUE)
}

summarize_collapsed_simulation_results <- function(results) {
  results$nu0_true <- as.numeric(as.character(results$nu0_true))
  means <- stats::aggregate(
    cbind(selected_L, coverage, power, all_causal_covered, fit_seconds) ~
      method + nu0_true,
    data = results,
    FUN = mean
  )
  eta0_means <- stats::aggregate(
    nu0_hat ~ method + nu0_true,
    data = results,
    FUN = mean_or_missing,
    na.action = na.pass
  )
  names(eta0_means)[3] <- "mean_nu0_hat"
  eta0_medians <- stats::aggregate(
    nu0_hat ~ method + nu0_true,
    data = results,
    FUN = median_or_missing,
    na.action = na.pass
  )
  names(eta0_medians)[3] <- "median_nu0_hat"

  summary <- Reduce(
    function(x, y) merge(x, y, by = c("method", "nu0_true"), all = TRUE),
    list(means, eta0_means, eta0_medians)
  )
  summary$nu0_true <- as.numeric(as.character(summary$nu0_true))
  summary[order(summary$nu0_true, summary$method), , drop = FALSE]
}

resolve_simulation_output_prefix <- function(output_prefix) {
  if (grepl("^/", output_prefix)) {
    output_prefix
  } else {
    file.path(susie_iw_root, output_prefix)
  }
}

run_collapsed_iw_simulation <- function(
    number_of_replicates = 5L,
    true_eta0_grid = unique(round(exp(seq(log(20), log(10000), length.out = 9)))),
    fitting_eta0_grid = NULL,
    sample_size = 20000L,
    coalescent_sample_size = 2000L,
    reference_sample_size = 2000L,
    number_of_variants = 500L,
    heritability = 0.01,
    true_number_of_effects = 2L,
    fitted_number_of_effects = 5L,
    lambda = 1e-3,
    coverage = 0.95,
    minimum_absolute_correlation = 0.5,
    simulation_seed = 1L,
    reference_panel_seeds = c(3L, 4L, 5L),
    output_prefix = "experiments/simulation-iw-collapsed-R",
    collapsed_maximum_iterations = 500L,
    collapsed_tolerance = 1e-7,
    save_plots = TRUE,
    verbose = TRUE) {
  if (is.null(fitting_eta0_grid)) {
    fitting_eta0_grid <- sort(unique(c(
      exp(seq(log(2), log(20000), length.out = 31)),
      true_eta0_grid
    )))
  }
  output_prefix <- resolve_simulation_output_prefix(output_prefix)
  figure_prefix <- file.path(
    dirname(output_prefix),
    "figures",
    basename(output_prefix)
  )
  dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(figure_prefix), recursive = TRUE, showWarnings = FALSE)

  number_of_replicates <- as.integer(number_of_replicates)
  reference_panel_seeds <- as.integer(reference_panel_seeds)
  total_runs <- length(reference_panel_seeds) *
    length(true_eta0_grid) * number_of_replicates
  runs <- vector("list", total_runs)
  profiles <- vector("list", total_runs)
  run_index <- 0L

  for (reference_panel_seed in reference_panel_seeds) {
    R0 <- simulate_coalescent_reference_ld(
      coalescent_sample_size = coalescent_sample_size,
      reference_sample_size = reference_sample_size,
      number_of_variants = number_of_variants,
      seed = reference_panel_seed,
      heritability = heritability
    )
    reference_ld <- make_stabilized_reference_ld(R0, lambda = lambda)

    for (true_eta0 in true_eta0_grid) {
      for (replicate_index in seq_len(number_of_replicates)) {
        run_index <- run_index + 1L
        replicate_seed <- simulation_seed + 100000L + run_index
        run <- simulate_one_collapsed_replicate(
          true_eta0 = true_eta0,
          replicate_index = replicate_index,
          seed = replicate_seed,
          R0 = R0,
          reference_ld = reference_ld,
          sample_size = sample_size,
          heritability = heritability,
          true_number_of_effects = true_number_of_effects,
          fitted_number_of_effects = fitted_number_of_effects,
          fitting_eta0_grid = fitting_eta0_grid,
          coverage = coverage,
          minimum_absolute_correlation = minimum_absolute_correlation,
          lambda = lambda,
          collapsed_maximum_iterations = collapsed_maximum_iterations,
          collapsed_tolerance = collapsed_tolerance,
          verbose = verbose
        )
        run$results$r0_seed <- reference_panel_seed
        run$results$N0 <- reference_sample_size
        run$results$N_coal <- coalescent_sample_size
        run$profile$r0_seed <- reference_panel_seed
        runs[[run_index]] <- run$results
        profiles[[run_index]] <- run$profile

        partial_results <- do.call(rbind, runs[seq_len(run_index)])
        partial_results <- partial_results[
          order(
            partial_results$r0_seed,
            partial_results$nu0_true,
            partial_results$rep,
            partial_results$method
          ),
          ,
          drop = FALSE
        ]
        utils::write.csv(
          partial_results,
          paste0(output_prefix, "-results.csv"),
          row.names = FALSE
        )
      }
    }
  }

  results <- do.call(rbind, runs)
  profiles <- do.call(rbind, profiles)
  results$nu0_true <- as.numeric(as.character(results$nu0_true))
  method_order <- c("susie_insample_R", "susie_outsample_R0", "collapsed_r")
  results$method <- factor(results$method, levels = method_order)
  results <- results[
    order(results$r0_seed, results$nu0_true, results$rep, results$method),
    ,
    drop = FALSE
  ]
  summary <- summarize_collapsed_simulation_results(results)

  settings <- list(
    number_of_replicates = number_of_replicates,
    true_eta0_grid = true_eta0_grid,
    fitting_eta0_grid = fitting_eta0_grid,
    sample_size = sample_size,
    coalescent_sample_size = coalescent_sample_size,
    reference_sample_size = reference_sample_size,
    number_of_variants = number_of_variants,
    heritability = heritability,
    true_number_of_effects = true_number_of_effects,
    fitted_number_of_effects = fitted_number_of_effects,
    lambda = lambda,
    coverage = coverage,
    minimum_absolute_correlation = minimum_absolute_correlation,
    simulation_seed = simulation_seed,
    reference_panel_seeds = reference_panel_seeds,
    collapsed_maximum_iterations = collapsed_maximum_iterations,
    collapsed_tolerance = collapsed_tolerance
  )
  output <- list(
    results = results,
    summary = summary,
    profiles = profiles,
    settings = settings
  )

  utils::write.csv(results, paste0(output_prefix, "-results.csv"), row.names = FALSE)
  utils::write.csv(summary, paste0(output_prefix, "-summary.csv"), row.names = FALSE)
  saveRDS(output, paste0(output_prefix, "-summary.rds"))

  if (save_plots) {
    plot_collapsed_iw_simulation(
      results = results,
      summary = summary,
      figure_prefix = figure_prefix,
      true_number_of_effects = true_number_of_effects
    )
  }
  output
}

plot_collapsed_iw_simulation <- function(
    results, summary,
    figure_prefix = file.path(
      susie_iw_root,
      "analysis",
      "figures",
      "simulation-iw-collapsed-R"
    ),
    true_number_of_effects = 2L) {
  methods <- c("susie_insample_R", "susie_outsample_R0", "collapsed_r")
  labels <- c(
    susie_insample_R = "SuSiE in-sample R",
    susie_outsample_R0 = "SuSiE out-of-sample R0",
    collapsed_r = "Collapsed-R"
  )
  colors <- c(
    susie_insample_R = "#1f78b4",
    susie_outsample_R0 = "#33a02c",
    collapsed_r = "#e31a1c"
  )
  results$method <- as.character(results$method)
  summary$method <- as.character(summary$method)
  results$nu0_true <- as.numeric(as.character(results$nu0_true))
  summary$nu0_true <- as.numeric(as.character(summary$nu0_true))
  dir.create(dirname(figure_prefix), recursive = TRUE, showWarnings = FALSE)

  draw_selected_L <- function() {
    plot(
      range(results$nu0_true),
      range(c(0, results$selected_L, true_number_of_effects), finite = TRUE),
      type = "n",
      log = "x",
      xlab = expression("true " * eta[0]),
      ylab = "Purity-filtered number of CSs",
      main = expression("Collapsed-R reported signals against true " * eta[0])
    )
    abline(h = true_number_of_effects, lty = 2, col = "gray40")
    for (method in methods) {
      rows <- results$method == method
      set.seed(100 + match(method, methods))
      jittered_eta0 <- results$nu0_true[rows] *
        exp(stats::runif(sum(rows), -0.025, 0.025))
      points(
        jittered_eta0,
        results$selected_L[rows],
        pch = 16,
        col = grDevices::adjustcolor(colors[method], alpha.f = 0.35)
      )
      method_summary <- summary[summary$method == method, , drop = FALSE]
      method_summary <- method_summary[
        order(method_summary$nu0_true),
        ,
        drop = FALSE
      ]
      lines(
        method_summary$nu0_true,
        method_summary$selected_L,
        type = "b",
        pch = 17,
        lwd = 2,
        col = colors[method]
      )
    }
    grid()
    legend(
      "topright",
      legend = c(
        unname(labels[methods]),
        sprintf("true L = %d", true_number_of_effects)
      ),
      col = c(colors[methods], "gray40"),
      pch = c(rep(17, length(methods)), NA),
      lty = c(rep(1, length(methods)), 2),
      bty = "n"
    )
  }

  draw_metrics <- function() {
    old_parameters <- par(no.readonly = TRUE)
    on.exit(par(old_parameters), add = TRUE)
    par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
    for (metric in c("coverage", "power")) {
      title <- if (metric == "coverage") "CS coverage" else "Power"
      plot(
        range(summary$nu0_true),
        c(0, 1),
        type = "n",
        log = "x",
        xlab = expression("true " * eta[0]),
        ylab = title,
        main = title
      )
      for (method in methods) {
        method_summary <- summary[summary$method == method, , drop = FALSE]
        method_summary <- method_summary[
          order(method_summary$nu0_true),
          ,
          drop = FALSE
        ]
        lines(
          method_summary$nu0_true,
          method_summary[[metric]],
          type = "b",
          pch = 16,
          lwd = 2,
          col = colors[method]
        )
      }
      grid()
      if (metric == "coverage") {
        legend(
          "bottomright",
          legend = unname(labels[methods]),
          col = colors[methods],
          pch = 16,
          lty = 1,
          bty = "n",
          cex = 0.8
        )
      }
    }
  }

  draw_estimated_eta0 <- function() {
    collapsed_results <- results[
      results$method == "collapsed_r",
      ,
      drop = FALSE
    ]
    collapsed_summary <- summary[
      summary$method == "collapsed_r",
      ,
      drop = FALSE
    ]
    finite_eta0 <- collapsed_results$nu0_hat[
      is.finite(collapsed_results$nu0_hat) & collapsed_results$nu0_hat > 0
    ]
    plot(
      range(collapsed_results$nu0_true),
      range(c(collapsed_results$nu0_true, finite_eta0)),
      type = "n",
      log = "xy",
      xlab = expression("true " * eta[0]),
      ylab = expression("selected " * eta[0]),
      main = expression("Collapsed-R selected " * eta[0])
    )
    boundary <- collapsed_results$eta0_grid_boundary != "interior"
    points(
      collapsed_results$nu0_true[!boundary],
      collapsed_results$nu0_hat[!boundary],
      pch = 16,
      col = grDevices::adjustcolor(colors["collapsed_r"], alpha.f = 0.4)
    )
    points(
      collapsed_results$nu0_true[boundary],
      collapsed_results$nu0_hat[boundary],
      pch = 24,
      bg = grDevices::adjustcolor(colors["collapsed_r"], alpha.f = 0.4),
      col = colors["collapsed_r"]
    )
    abline(0, 1, lty = 2, col = "gray40")
    lines(
      collapsed_summary$nu0_true,
      collapsed_summary$median_nu0_hat,
      type = "b",
      pch = 17,
      lwd = 2,
      col = colors["collapsed_r"]
    )
    grid()
    legend(
      "topleft",
      legend = c("interior grid point", "grid boundary", "median", "y = x"),
      col = c(
        grDevices::adjustcolor(colors["collapsed_r"], alpha.f = 0.4),
        colors["collapsed_r"], colors["collapsed_r"], "gray40"
      ),
      pch = c(16, 24, 17, NA),
      pt.bg = c(
        NA,
        grDevices::adjustcolor(colors["collapsed_r"], alpha.f = 0.4),
        NA,
        NA
      ),
      lty = c(NA, NA, 1, 2),
      bty = "n"
    )
  }

  grDevices::png(
    paste0(figure_prefix, "-selected-L.png"),
    width = 1400,
    height = 850,
    res = 150
  )
  draw_selected_L()
  grDevices::dev.off()
  grDevices::pdf(
    paste0(figure_prefix, "-selected-L.pdf"),
    width = 8,
    height = 5.5
  )
  draw_selected_L()
  grDevices::dev.off()

  grDevices::png(
    paste0(figure_prefix, "-summary-metrics.png"),
    width = 1400,
    height = 650,
    res = 150
  )
  draw_metrics()
  grDevices::dev.off()
  grDevices::pdf(
    paste0(figure_prefix, "-summary-metrics.pdf"),
    width = 11,
    height = 4.5
  )
  draw_metrics()
  grDevices::dev.off()

  grDevices::png(
    paste0(figure_prefix, "-estimated-eta0.png"),
    width = 1400,
    height = 850,
    res = 150
  )
  draw_estimated_eta0()
  grDevices::dev.off()
  grDevices::pdf(
    paste0(figure_prefix, "-estimated-eta0.pdf"),
    width = 7,
    height = 5.5
  )
  draw_estimated_eta0()
  grDevices::dev.off()

  invisible(TRUE)
}
