#!/usr/bin/env Rscript

# Combine the frozen four-method benchmark with the exploratory susieR
# finite-reference + empirical-Bayes mismatch fits. All output names are new;
# this script never overwrites the manuscript benchmark or its figures.

library(data.table)

original_results_file <- Sys.getenv(
  "COAL_EB_ORIGINAL_RESULTS",
  file.path("experiments", "coalescent-susie-iw-benchmark-results.csv")
)
worker_directory <- Sys.getenv(
  "COAL_EB_OUTPUT_DIRECTORY",
  file.path("experiments", "coalescent-susie-rss-eb-benchmark-20rep-workers")
)
output_prefix <- Sys.getenv(
  "COAL_EB_SUMMARY_PREFIX",
  file.path("experiments", "coalescent-susie-iw-plus-rss-eb-benchmark")
)
figure_prefix <- Sys.getenv(
  "COAL_EB_FIGURE_PREFIX",
  file.path("experiments", "figures", "coalescent-susie-iw-plus-rss-eb")
)
bootstrap_replicates <- as.integer(Sys.getenv("COAL_EB_BOOTSTRAPS", "2000"))

if (!file.exists(original_results_file)) {
  stop("Frozen benchmark results not found: ", original_results_file)
}
files <- list.files(
  worker_directory,
  pattern = "^results-J[0-9]+-panel[0-9]+[.]csv$",
  full.names = TRUE
)
if (!length(files)) stop("No completed EB benchmark workers found")

original <- fread(original_results_file)
eb_results <- rbindlist(lapply(files, fread), fill = TRUE)
new_method <- "SuSiE-RSS: finite R0 + EB mismatch"
original_methods <- c(
  "SuSiE-RSS: in-sample R", "SuSiE-RSS: reference R0",
  "collapsed-R", "AIW-N"
)
expected_methods <- c(original_methods[1:2], new_method, original_methods[3:4])

stopifnot(
  nrow(eb_results) > 0L,
  all(eb_results$method == new_method),
  setequal(unique(original$method), original_methods)
)

# The frozen result table supplies the replication labels and the four existing
# methods. Subsetting it to exactly the configurations completed by the EB
# workers also makes small smoke-test runs summarizable.
task_map <- unique(original[, .(J, panel_seed, replication, evaluation_panel)])
if ("replication" %in% names(eb_results)) eb_results[, replication := NULL]
if ("evaluation_panel" %in% names(eb_results)) {
  eb_results[, evaluation_panel := NULL]
}
eb_results <- merge(
  eb_results, task_map, by = c("J", "panel_seed"), all.x = TRUE
)
if (anyNA(eb_results$replication)) {
  stop("At least one EB worker does not match a frozen benchmark panel")
}
eb_results[, configuration_cs_coverage := fifelse(
  reported_L > 0, true_sets / reported_L, NA_real_
)]

configuration_keys <- unique(
  eb_results[, .(J, panel_seed, true_eta0, true_L)]
)
original_subset <- original[
  configuration_keys,
  on = .(J, panel_seed, true_eta0, true_L),
  nomatch = 0L
]
results <- rbindlist(list(original_subset, eb_results), fill = TRUE)
setorder(results, true_L, true_eta0, method, replication)

true_eta0_values <- sort(unique(eb_results$true_eta0))
true_L_values <- sort(unique(eb_results$true_L))
expected_configurations <- unique(eb_results[, .(
  J, panel_seed, true_eta0, true_L
)])
method_counts <- results[, .N, by = .(J, panel_seed, true_eta0, true_L, method)]
stopifnot(
  all(method_counts$N == 1L),
  setequal(results$method, expected_methods),
  nrow(results) == nrow(expected_configurations) * length(expected_methods)
)

results[, method := factor(method, levels = expected_methods)]
results[, evaluation_panel := as.character(replication)]

metric_functions <- list(
  mean_reported_L = function(x) mean(x$reported_L),
  causal_power = function(x) mean(x$causal_coverage),
  all_causal_power = function(x) mean(x$all_causal_covered),
  cs_coverage = function(x) {
    if (sum(x$reported_L) == 0) return(NA_real_)
    sum(x$true_sets) / sum(x$reported_L)
  },
  false_sets = function(x) mean(x$false_sets),
  fallback_frequency = function(x) mean(x$fallback),
  partition_frequency = function(x) mean(x$partitioned_SER),
  mean_set_size = function(x) mean(x$mean_set_size, na.rm = TRUE),
  normalized_d_cs = function(x) mean(x$d_cs_normalized, na.rm = TRUE),
  mean_causal_pip = function(x) mean(x$mean_causal_pip),
  convergence_frequency = function(x) mean(x$converged),
  median_seconds = function(x) median(x$fit_seconds),
  median_raw_eta0 = function(x) {
    values <- x$raw_eta0[is.finite(x$raw_eta0)]
    if (length(values)) median(values) else NA_real_
  },
  median_eta0_used = function(x) {
    values <- x$eta0_used[is.finite(x$eta0_used)]
    if (length(values)) median(values) else NA_real_
  },
  median_lambda = function(x) {
    values <- x$selected_lambda[is.finite(x$selected_lambda)]
    if (length(values)) median(values) else NA_real_
  }
)

summarize_cell <- function(data) {
  estimates <- vapply(metric_functions, function(f) f(data), numeric(1L))
  interval_metrics <- c(
    "mean_reported_L", "causal_power", "all_causal_power", "cs_coverage",
    "false_sets", "fallback_frequency", "mean_set_size", "normalized_d_cs"
  )
  intervals <- matrix(
    NA_real_, nrow = length(metric_functions), ncol = 2L,
    dimnames = list(names(metric_functions), c("lower", "upper"))
  )
  set.seed(
    730000L + 1000L * data$true_L[1L] + data$true_eta0[1L] +
      as.integer(data$method[1L])
  )
  indices <- matrix(
    replicate(
      bootstrap_replicates,
      sample.int(nrow(data), nrow(data), replace = TRUE)
    ),
    nrow = nrow(data)
  )
  bootstrap_mean <- function(values, na.rm = FALSE) {
    sampled <- matrix(values[indices], nrow = nrow(indices))
    colMeans(sampled, na.rm = na.rm)
  }
  bootstrap_values <- list(
    mean_reported_L = bootstrap_mean(data$reported_L),
    causal_power = bootstrap_mean(data$causal_coverage),
    all_causal_power = bootstrap_mean(data$all_causal_covered),
    cs_coverage = {
      sampled_true <- matrix(data$true_sets[indices], nrow = nrow(indices))
      sampled_reported <- matrix(
        data$reported_L[indices], nrow = nrow(indices)
      )
      denominator <- colSums(sampled_reported)
      fifelse(denominator > 0, colSums(sampled_true) / denominator, NA_real_)
    },
    false_sets = bootstrap_mean(data$false_sets),
    fallback_frequency = bootstrap_mean(data$fallback),
    mean_set_size = bootstrap_mean(data$mean_set_size, na.rm = TRUE),
    normalized_d_cs = bootstrap_mean(data$d_cs_normalized, na.rm = TRUE)
  )
  for (metric in interval_metrics) {
    intervals[metric, ] <- quantile(
      bootstrap_values[[metric]], c(0.025, 0.975),
      na.rm = TRUE, names = FALSE
    )
  }
  output <- as.list(estimates)
  for (metric in interval_metrics) {
    output[[paste0(metric, "_lower")]] <- intervals[metric, "lower"]
    output[[paste0(metric, "_upper")]] <- intervals[metric, "upper"]
  }
  as.data.table(output)
}

summary_rows <- list()
counter <- 0L
for (L_value in true_L_values) {
  for (eta_value in true_eta0_values) {
    for (method_value in expected_methods) {
      counter <- counter + 1L
      cell_data <- results[
        true_L == L_value & true_eta0 == eta_value &
          as.character(method) == method_value
      ]
      row <- summarize_cell(cell_data)
      row[, `:=`(
        true_L = L_value, true_eta0 = eta_value, method = method_value,
        configurations = nrow(cell_data)
      )]
      summary_rows[[counter]] <- row
    }
  }
}
summary <- rbindlist(summary_rows, fill = TRUE)
summary[, method := factor(method, levels = expected_methods)]
setcolorder(summary, c("true_L", "true_eta0", "method", "configurations"))

summarize_average_cell <- function(data) {
  average_over_L <- function(x, metric) {
    mean(vapply(true_L_values, function(L_value) {
      metric_functions[[metric]](x[true_L == L_value])
    }, numeric(1L)))
  }
  estimate <- c(
    causal_power = average_over_L(data, "causal_power"),
    cs_coverage = average_over_L(data, "cs_coverage")
  )
  panels <- unique(data$evaluation_panel)
  ordered <- copy(data)
  ordered[, panel_order := match(evaluation_panel, panels)]
  setorder(ordered, panel_order, true_L)
  causal_matrix <- matrix(
    ordered$causal_coverage, nrow = length(panels),
    ncol = length(true_L_values), byrow = TRUE
  )
  true_set_matrix <- matrix(
    ordered$true_sets, nrow = length(panels),
    ncol = length(true_L_values), byrow = TRUE
  )
  reported_matrix <- matrix(
    ordered$reported_L, nrow = length(panels),
    ncol = length(true_L_values), byrow = TRUE
  )
  set.seed(810000L + data$true_eta0[1L] + as.integer(data$method[1L]))
  sampled_indices <- matrix(
    replicate(
      bootstrap_replicates,
      sample.int(length(panels), length(panels), replace = TRUE)
    ),
    nrow = length(panels)
  )
  bootstrap <- vapply(seq_len(bootstrap_replicates), function(index) {
    sampled <- sampled_indices[, index]
    c(
      causal_power = mean(causal_matrix[sampled, , drop = FALSE]),
      cs_coverage = mean(
        colSums(true_set_matrix[sampled, , drop = FALSE]) /
          colSums(reported_matrix[sampled, , drop = FALSE])
      )
    )
  }, numeric(2L))
  data.table(
    causal_power = estimate["causal_power"],
    causal_power_lower = quantile(
      bootstrap["causal_power", ], 0.025, names = FALSE
    ),
    causal_power_upper = quantile(
      bootstrap["causal_power", ], 0.975, names = FALSE
    ),
    cs_coverage = estimate["cs_coverage"],
    cs_coverage_lower = quantile(
      bootstrap["cs_coverage", ], 0.025, names = FALSE
    ),
    cs_coverage_upper = quantile(
      bootstrap["cs_coverage", ], 0.975, names = FALSE
    )
  )
}

average_rows <- list()
counter <- 0L
for (eta_value in true_eta0_values) {
  for (method_value in expected_methods) {
    counter <- counter + 1L
    cell_data <- results[
      true_eta0 == eta_value & as.character(method) == method_value
    ]
    row <- summarize_average_cell(cell_data)
    row[, `:=`(
      true_eta0 = eta_value, method = method_value,
      configurations = nrow(cell_data)
    )]
    average_rows[[counter]] <- row
  }
}
average_summary <- rbindlist(average_rows, fill = TRUE)
average_summary[, method := factor(method, levels = expected_methods)]
setcolorder(average_summary, c("true_eta0", "method", "configurations"))

eb_diagnostics <- eb_results[, .(
  configurations = .N,
  convergence_frequency = mean(converged),
  reliability_flag_frequency = mean(overall_reliability_flag, na.rm = TRUE),
  artifact_flag_frequency = mean(residual_artifact_flag, na.rm = TRUE),
  median_lambda_bias = median(lambda_bias, na.rm = TRUE),
  median_corrected_reference_size = median(
    corrected_reference_size, na.rm = TRUE
  ),
  median_Q_art = median(Q_art, na.rm = TRUE),
  median_iterations = median(iterations, na.rm = TRUE),
  median_seconds = median(fit_seconds, na.rm = TRUE)
), by = .(true_L, true_eta0)]
setorder(eb_diagnostics, true_L, true_eta0)

dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(figure_prefix), recursive = TRUE, showWarnings = FALSE)
fwrite(results, paste0(output_prefix, "-results.csv"))
fwrite(summary, paste0(output_prefix, "-summary.csv"))
fwrite(average_summary, paste0(output_prefix, "-average-over-L-summary.csv"))
fwrite(eb_diagnostics, paste0(output_prefix, "-eb-diagnostics.csv"))
saveRDS(
  list(
    results = results, summary = summary,
    average_summary = average_summary, eb_diagnostics = eb_diagnostics
  ),
  paste0(output_prefix, ".rds")
)

method_colors <- c("#222222", "#d73027", "#1f78b4", "#2ca02c", "#6a3d9a")
method_shapes <- c(15, 17, 8, 16, 18)
names(method_colors) <- names(method_shapes) <- expected_methods

draw_panel <- function(data, metric, ylab, main, ylim = NULL,
                       legend = FALSE, legend_position = "topright",
                       truth = NULL) {
  x <- true_eta0_values
  values <- range(
    c(data[[paste0(metric, "_lower")]], data[[paste0(metric, "_upper")]],
      truth),
    na.rm = TRUE
  )
  if (is.null(ylim)) {
    padding <- max(0.04, diff(values) * 0.08)
    ylim <- values + c(-padding, padding)
  }
  plot(
    x, rep(NA_real_, length(x)), log = if (length(x) > 1L) "x" else "",
    xlim = if (length(x) > 1L) range(x) else x + c(-1, 1), ylim = ylim,
    xlab = expression("true mismatch concentration " * eta[0]),
    ylab = ylab, main = main, xaxt = "n"
  )
  axis(1, at = x, labels = format(x, big.mark = ",", scientific = FALSE))
  grid(col = "grey88", lty = 3)
  if (!is.null(truth)) abline(h = truth, col = "grey35", lty = 2)
  for (method_value in expected_methods) {
    method_data <- data[as.character(method) == method_value]
    lines(
      method_data$true_eta0, method_data[[metric]],
      col = method_colors[method_value], lwd = 2,
      type = "b", pch = method_shapes[method_value]
    )
    lower <- method_data[[paste0(metric, "_lower")]]
    upper <- method_data[[paste0(metric, "_upper")]]
    show_interval <- is.finite(lower) & is.finite(upper) & lower != upper
    if (any(show_interval)) {
      arrows(
        method_data$true_eta0[show_interval], lower[show_interval],
        method_data$true_eta0[show_interval], upper[show_interval],
        angle = 90, code = 3, length = 0.025,
        col = method_colors[method_value]
      )
    }
  }
  if (legend) legend(
    legend_position, legend = expected_methods, col = method_colors,
    pch = method_shapes, lwd = 2, bty = "n", cex = 0.65
  )
}

draw_selected_L <- function() {
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  par(
    mfrow = c(1, length(true_L_values)), mar = c(4.4, 4.4, 3, 1),
    las = 1
  )
  for (L_value in true_L_values) {
    draw_panel(
      summary[true_L == L_value], "mean_reported_L",
      "reported credible sets", paste("true L =", L_value),
      ylim = c(0.5, 5.2), legend = L_value == max(true_L_values),
      truth = L_value
    )
  }
}

draw_performance_average <- function() {
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  par(mfrow = c(1, 2), mar = c(4.4, 4.4, 3, 1), las = 1)
  for (metric in c("cs_coverage", "causal_power")) {
    draw_panel(
      average_summary, metric,
      if (metric == "cs_coverage") {
        "credible-set coverage"
      } else {
        "causal-variant power"
      },
      if (metric == "cs_coverage") "Coverage" else "Power",
      ylim = c(0, 1.03), legend = metric == "causal_power",
      legend_position = "bottomright"
    )
  }
}

draw_performance_for_L <- function(L_value) {
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  par(mfrow = c(1, 2), mar = c(4.4, 4.4, 3, 1), las = 1)
  for (metric in c("cs_coverage", "causal_power")) {
    draw_panel(
      summary[true_L == L_value], metric,
      if (metric == "cs_coverage") {
        "credible-set coverage"
      } else {
        "causal-variant power"
      },
      if (metric == "cs_coverage") "Coverage" else "Power",
      ylim = c(0, 1.03), legend = metric == "causal_power",
      legend_position = "bottomright"
    )
  }
}

write_figure <- function(stem, width, height, draw) {
  png(
    paste0(stem, ".png"), width = width * 180, height = height * 180,
    res = 180
  )
  draw()
  dev.off()
  pdf(paste0(stem, ".pdf"), width = width, height = height)
  draw()
  dev.off()
}

selected_width <- max(4.5, 4 * length(true_L_values))
write_figure(
  paste0(figure_prefix, "-selected-L"), selected_width, 4.5, draw_selected_L
)
write_figure(
  paste0(figure_prefix, "-performance"), 10, 4.6,
  draw_performance_average
)
for (L_value in true_L_values) {
  local({
    current_L <- L_value
    write_figure(
      paste0(figure_prefix, "-performance-L", current_L),
      10, 4.6, function() draw_performance_for_L(current_L)
    )
  })
}

message("Combined ", nrow(results), " rows across ", length(expected_methods),
        " methods")
message("Wrote exploratory outputs with prefix ", output_prefix)
message("Wrote exploratory figures with prefix ", figure_prefix)
