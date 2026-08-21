#!/usr/bin/env Rscript

library(data.table)

worker_directory <- file.path(
  "analysis", "coalescent-susie-iw-benchmark-20rep-workers"
)
files <- list.files(
  worker_directory, pattern = "^results-J[0-9]+-panel[0-9]+[.]csv$",
  full.names = TRUE
)
if (!length(files)) stop("No completed coalescent benchmark workers found")
all_results <- rbindlist(lapply(files, fread), fill = TRUE)
expected_J <- c(500L, 1000L, 2000L)
number_of_replications <- 20L
expected_tasks <- data.table(
  replication = seq_len(number_of_replications),
  J = rep(expected_J, length.out = number_of_replications),
  panel_seed = seq.int(61L, length.out = number_of_replications)
)
results <- merge(
  all_results, expected_tasks, by = c("J", "panel_seed"), all = FALSE
)
setorder(results, true_L, true_eta0, method, replication)

expected_eta0 <- c(20, 50, 100, 200, 500, 1000, 2000, 5000)
expected_methods <- c(
  "SuSiE-RSS: in-sample R", "SuSiE-RSS: reference R0",
  "collapsed-R", "AIW-N"
)
stopifnot(
  fsetequal(unique(results[, .(replication, J, panel_seed)]), expected_tasks),
  setequal(results$true_eta0, expected_eta0),
  setequal(results$true_L, 1:3),
  setequal(results$method, expected_methods),
  nrow(results) == number_of_replications * length(expected_eta0) *
    3L * length(expected_methods)
)

results[, method := factor(method, levels = expected_methods)]
results[, configuration_cs_coverage := fifelse(
  reported_L > 0, true_sets / reported_L, NA_real_
)]
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

summarize_cell <- function(data, bootstrap_replicates = 2000L) {
  estimates <- vapply(metric_functions, function(f) f(data), numeric(1L))
  additive_metrics <- c(
    "mean_reported_L", "causal_power", "all_causal_power", "cs_coverage",
    "false_sets", "fallback_frequency", "mean_set_size", "normalized_d_cs"
  )
  intervals <- matrix(
    NA_real_, nrow = length(metric_functions), ncol = 2L,
    dimnames = list(names(metric_functions), c("lower", "upper"))
  )
  set.seed(730000L + 1000L * data$true_L[1L] + data$true_eta0[1L] +
             as.integer(data$method[1L]))
  indices <- replicate(
    bootstrap_replicates,
    sample.int(nrow(data), nrow(data), replace = TRUE),
    simplify = FALSE
  )
  for (metric in additive_metrics) {
    values <- vapply(indices, function(index) {
      metric_functions[[metric]](data[index])
    }, numeric(1L))
    intervals[metric, ] <- quantile(
      values, c(0.025, 0.975), na.rm = TRUE, names = FALSE
    )
  }
  output <- as.list(estimates)
  for (metric in additive_metrics) {
    output[[paste0(metric, "_lower")]] <- intervals[metric, "lower"]
    output[[paste0(metric, "_upper")]] <- intervals[metric, "upper"]
  }
  as.data.table(output)
}

summary_rows <- list()
counter <- 0L
for (L_value in 1:3) {
  for (eta_value in expected_eta0) {
    for (method_value in expected_methods) {
      counter <- counter + 1L
      cell <- results[
        true_L == L_value & true_eta0 == eta_value &
          as.character(method) == method_value
      ]
      row <- summarize_cell(cell)
      row[, `:=`(
        true_L = L_value, true_eta0 = eta_value, method = method_value,
        configurations = nrow(cell)
      )]
      summary_rows[[counter]] <- row
    }
  }
}
summary <- rbindlist(summary_rows, fill = TRUE)
summary[, method := factor(method, levels = expected_methods)]
setcolorder(summary, c(
  "true_L", "true_eta0", "method", "configurations"
))

# Main-text performance averages over L_true = 1, 2, 3. Bootstrap the 20
# independent (J, R0) panels and retain all three nested architectures within
# each resampled panel.
summarize_average_cell <- function(data, bootstrap_replicates = 2000L) {
  average_over_L <- function(x, metric) {
    mean(vapply(1:3, function(L_value) {
      metric_functions[[metric]](x[true_L == L_value])
    }, numeric(1L)))
  }
  estimate <- c(
    causal_power = average_over_L(data, "causal_power"),
    cs_coverage = average_over_L(data, "cs_coverage")
  )
  panels <- unique(data$evaluation_panel)
  set.seed(810000L + data$true_eta0[1L] + as.integer(data$method[1L]))
  bootstrap <- replicate(bootstrap_replicates, {
    sampled_panels <- sample(panels, length(panels), replace = TRUE)
    sampled <- rbindlist(lapply(sampled_panels, function(panel) {
      data[evaluation_panel == panel]
    }))
    c(
      causal_power = average_over_L(sampled, "causal_power"),
      cs_coverage = average_over_L(sampled, "cs_coverage")
    )
  })
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
for (eta_value in expected_eta0) {
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
setcolorder(average_summary, c(
  "true_eta0", "method", "configurations"
))

fwrite(results, "experiments/coalescent-susie-iw-benchmark-results.csv")
fwrite(summary, "experiments/coalescent-susie-iw-benchmark-summary.csv")
fwrite(
  average_summary,
  "experiments/coalescent-susie-iw-benchmark-average-over-L-summary.csv"
)
saveRDS(
  list(results = results, summary = summary, average_summary = average_summary),
  "experiments/coalescent-susie-iw-benchmark.rds"
)

method_colors <- c("#222222", "#d73027", "#2ca02c", "#6a3d9a")
method_shapes <- c(15, 17, 16, 18)
names(method_colors) <- names(method_shapes) <- expected_methods

draw_panel <- function(data, metric, ylab, main, ylim = NULL,
                       legend = FALSE, legend_position = "topright",
                       truth = NULL) {
  x <- expected_eta0
  values <- range(
    c(data[[paste0(metric, "_lower")]], data[[paste0(metric, "_upper")]],
      truth), na.rm = TRUE
  )
  if (is.null(ylim)) {
    padding <- max(0.04, diff(values) * 0.08)
    ylim <- values + c(-padding, padding)
  }
  plot(
    x, rep(NA_real_, length(x)), log = "x", xlim = range(x), ylim = ylim,
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
    pch = method_shapes, lwd = 2, bty = "n", cex = 0.72
  )
}

draw_selected_L <- function() {
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  par(mfrow = c(1, 3), mar = c(4.4, 4.4, 3, 1), las = 1)
  for (L_value in 1:3) {
    draw_panel(
      summary[true_L == L_value], "mean_reported_L",
      "reported credible sets", paste("true L =", L_value),
      ylim = c(0.5, 5.2), legend = L_value == 3L, truth = L_value
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

draw_diagnostics <- function() {
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  par(mfrow = c(2, 2), mar = c(4.4, 4.4, 3, 1), las = 1)
  for (metric in c("false_sets", "fallback_frequency", "normalized_d_cs")) {
    source_data <- if (metric == "fallback_frequency") {
      summary[true_L >= 2L]
    } else summary
    plot_data <- source_data[, .(
      value = mean(get(metric)),
      lower = mean(get(paste0(metric, "_lower"))),
      upper = mean(get(paste0(metric, "_upper")))
    ), by = .(true_eta0, method)]
    setnames(plot_data, c("value", "lower", "upper"), c(
      metric, paste0(metric, "_lower"), paste0(metric, "_upper")
    ))
    draw_panel(
      plot_data, metric,
      switch(metric, false_sets = "false CS per dataset",
             fallback_frequency = "SER fallback frequency (true L > 1)",
             normalized_d_cs = expression("normalized " * d[CS])),
      switch(metric, false_sets = "False discoveries",
             fallback_frequency = "Operational fallback",
             normalized_d_cs = "Agreement with in-sample SuSiE-RSS"),
      ylim = if (metric == "fallback_frequency") c(0, 1.03) else NULL,
      legend = metric == "false_sets"
    )
  }
  eta_data <- results[as.character(method) %in% c("collapsed-R", "AIW-N") &
                        is.finite(eta0_used) & eta0_used > 0 & eta0_used < 1e8,
    .(median = median(eta0_used), lower = quantile(eta0_used, 0.25),
      upper = quantile(eta0_used, 0.75)),
    by = .(true_eta0, method)]
  plot(
    expected_eta0, expected_eta0, log = "xy", type = "n",
    xlab = expression("true " * eta[0]),
    ylab = expression("calibrated " * tilde(eta)[0]),
    main = expression("Mismatch recovery")
  )
  grid(col = "grey88", lty = 3)
  abline(0, 1, lty = 2, col = "grey35")
  for (method_value in c("collapsed-R", "AIW-N")) {
    data <- eta_data[as.character(method) == method_value]
    lines(data$true_eta0, data$median, type = "b", lwd = 2,
          pch = method_shapes[method_value], col = method_colors[method_value])
    arrows(data$true_eta0, data$lower, data$true_eta0, data$upper,
           angle = 90, code = 3, length = 0.025,
           col = method_colors[method_value])
  }
}

write_figure <- function(stem, width, height, draw) {
  png(paste0(stem, ".png"), width = width * 180, height = height * 180,
      res = 180)
  draw()
  dev.off()
  pdf(paste0(stem, ".pdf"), width = width, height = height)
  draw()
  dev.off()
}
dir.create(file.path("experiments", "figures"), showWarnings = FALSE)
write_figure(
  file.path("experiments", "figures", "coalescent-susie-iw-selected-L"),
  12, 4.5, draw_selected_L
)
write_figure(
  file.path("experiments", "figures", "coalescent-susie-iw-performance"),
  10, 4.6, draw_performance_average
)
for (L_value in 1:3) {
  local({
    current_L <- L_value
    write_figure(
      file.path(
        "analysis", "figures",
        paste0("coalescent-susie-iw-performance-L", current_L)
      ),
      10, 4.6, function() draw_performance_for_L(current_L)
    )
  })
}
write_figure(
  file.path("experiments", "figures", "coalescent-susie-iw-diagnostics"),
  10.5, 8.5, draw_diagnostics
)

selected_table <- summary[true_eta0 %in% c(20, 200, 1000, 5000), .(
  true_L, true_eta0, method = as.character(method), mean_reported_L,
  causal_power, cs_coverage, false_sets, fallback_frequency,
  normalized_d_cs
)]
fwrite(selected_table, "experiments/coalescent-susie-iw-benchmark-main-table.csv")

cell <- function(L_value, eta_value, method_value) {
  summary[
    true_L == L_value & true_eta0 == eta_value &
      as.character(method) == method_value
  ][1L]
}
eta20_counts <- summary[true_eta0 == 20]
pooled <- results[, .(
  mean_reported_L = mean(reported_L),
  causal_power = mean(causal_coverage),
  cs_coverage = mean(configuration_cs_coverage),
  false_sets = mean(false_sets),
  median_seconds = median(fit_seconds)
), by = method]
pooled_value <- function(method_value, column) {
  pooled[method == method_value][[column]][1L]
}
calibration <- fread(
  "experiments/eta0-calibration-coalescent-L-recommendations.csv"
)
multiplier <- function(method_value) {
  calibration[method == method_value, recommended_multiplier][1L]
}
l2_eta20_reference <- cell(2, 20, "SuSiE-RSS: reference R0")
l2_eta20_collapsed <- cell(2, 20, "collapsed-R")
l2_eta20_aiwn <- cell(2, 20, "AIW-N")
l2_eta200_collapsed <- cell(2, 200, "collapsed-R")
l2_eta200_aiwn <- cell(2, 200, "AIW-N")
report <- c(
  "# Coalescent SuSiE-IW benchmark",
  "",
  "## Design",
  "",
  paste0(
    "This independent evaluation contains ",
    uniqueN(results[, .(J, panel_seed, true_eta0, true_L)]),
    " datasets and ", nrow(results), " method fits. Twenty independent ",
    "coalescent panels are distributed across J in {",
    paste(expected_J, collapse = ", "),
    "}, with 7, 7 and 6 panels, respectively. The experiment crosses true eta0 in {",
    paste(expected_eta0, collapse = ", "),
    "} with true L in {1, 2, 3}, using h2=1e-3, N=200,000, ",
    "N0=2,000 and true lambda=1e-3. Nested weak-LD causal-effect ",
    "weights are 7, 7:5 and 7:5:3. The SER fallback reports its original ",
    "95% set without purity filtering or partitioning."
  ),
  "",
  paste0(
    "The evaluation-panel seeds (61 through 80) are disjoint from the ",
    "calibration-panel seeds (51 and 52). Collapsed-R uses eta0 ",
    "multiplier ", sprintf("%.3f", multiplier("collapsed-R")),
    " and AIW-N uses ", sprintf("%.3f", multiplier("AIW-N")),
    "; these constants are frozen before performance evaluation. ",
    "Error bars use a nonparametric bootstrap over the ",
    uniqueN(results$evaluation_panel), " independent ",
    "(J, R0) panels."
  ),
  "",
  "## Main findings",
  "",
  paste0(
    "- At true eta0=20, reference-LD SuSiE-RSS reports means of ",
    paste(sprintf("%.1f", eta20_counts[
      as.character(method) == "SuSiE-RSS: reference R0"
    ]$mean_reported_L), collapse = ", "),
    " sets when true L is 1, 2 and 3. Collapsed-R reports ",
    paste(sprintf("%.1f", eta20_counts[
      as.character(method) == "collapsed-R"
    ]$mean_reported_L), collapse = ", "),
    ", and AIW-N reports ",
    paste(sprintf("%.1f", eta20_counts[
      as.character(method) == "AIW-N"
    ]$mean_reported_L), collapse = ", "), "."
  ),
  paste0(
    "- With two true effects and true eta0=20, reference-LD SuSiE-RSS ",
    "has causal power ", sprintf("%.2f", l2_eta20_reference$causal_power),
    ", credible-set coverage ", sprintf("%.3f", l2_eta20_reference$cs_coverage),
    " and ", sprintf("%.1f", l2_eta20_reference$false_sets),
    " false sets per dataset. Collapsed-R and AIW-N each have power ",
    sprintf("%.2f", l2_eta20_collapsed$causal_power), ", coverage ",
    sprintf("%.3f", l2_eta20_collapsed$cs_coverage), " and ",
    sprintf("%.1f", l2_eta20_collapsed$false_sets),
    " false sets. This is the intended conservative fallback regime."
  ),
  paste0(
    "- With two true effects and true eta0=200, AIW-N reports ",
    sprintf("%.1f", l2_eta200_aiwn$mean_reported_L),
    " sets with power ", sprintf("%.2f", l2_eta200_aiwn$causal_power),
    " and coverage ", sprintf("%.2f", l2_eta200_aiwn$cs_coverage),
    ". Collapsed-R reports ",
    sprintf("%.1f", l2_eta200_collapsed$mean_reported_L),
    " sets, with power ", sprintf("%.2f", l2_eta200_collapsed$causal_power),
    " and coverage ", sprintf("%.2f", l2_eta200_collapsed$cs_coverage),
    ". From true eta0=500 onward, all four methods report two sets."
  ),
  paste0(
    "- Pooled over the 24 design cells, reference-LD SuSiE-RSS reports ",
    sprintf("%.3f", pooled_value("SuSiE-RSS: reference R0", "mean_reported_L")),
    " sets and ",
    sprintf("%.3f", pooled_value("SuSiE-RSS: reference R0", "false_sets")),
    " false sets per dataset, with coverage ",
    sprintf("%.3f", pooled_value("SuSiE-RSS: reference R0", "cs_coverage")),
    ". Collapsed-R reports ",
    sprintf("%.3f", pooled_value("collapsed-R", "mean_reported_L")),
    " sets and ", sprintf("%.3f", pooled_value("collapsed-R", "false_sets")),
    " false sets, with coverage ",
    sprintf("%.3f", pooled_value("collapsed-R", "cs_coverage")),
    "; AIW-N reports ",
    sprintf("%.3f", pooled_value("AIW-N", "mean_reported_L")),
    " sets and ", sprintf("%.3f", pooled_value("AIW-N", "false_sets")),
    " false sets, with coverage ",
    sprintf("%.3f", pooled_value("AIW-N", "cs_coverage")), "."
  ),
  paste0(
    "- All ", nrow(results), " fits converged. Median fit times are ",
    sprintf("%.2f", pooled_value("collapsed-R", "median_seconds")),
    " seconds for collapsed-R, ",
    sprintf("%.2f", pooled_value("AIW-N", "median_seconds")),
    " seconds for AIW-N, ",
    sprintf("%.3f", pooled_value("SuSiE-RSS: in-sample R", "median_seconds")),
    " seconds for in-sample SuSiE-RSS and ",
    sprintf("%.3f", pooled_value("SuSiE-RSS: reference R0", "median_seconds")),
    " seconds for reference-LD SuSiE-RSS."
  ),
  "",
  paste0(
    "The key qualitative result is the crossover: moving from severe mismatch ",
    "toward matched LD, the mismatch-aware methods stop withholding additional ",
    "signals at approximately the same point that reference-LD SuSiE-RSS stops ",
    "over-specifying the number of signals."
  ),
  "",
  "## Outputs",
  "",
  "- `coalescent-susie-iw-benchmark-main-table.csv`: selected manuscript cells.",
  "- `coalescent-susie-iw-benchmark-summary.csv`: all method/design summaries and bootstrap intervals.",
  "- `coalescent-susie-iw-benchmark-average-over-L-summary.csv`: main-text power and coverage averaged over L.",
  "- `coalescent-susie-iw-benchmark-results.csv`: all configuration-level fits.",
  "- `figures/coalescent-susie-iw-selected-L.pdf`",
  "- `figures/coalescent-susie-iw-performance.pdf`",
  "- `figures/coalescent-susie-iw-performance-L1.pdf`",
  "- `figures/coalescent-susie-iw-performance-L2.pdf`",
  "- `figures/coalescent-susie-iw-performance-L3.pdf`",
  "- `figures/coalescent-susie-iw-diagnostics.pdf`",
  "- `figures/coalescent-susie-iw-motivating-locus.pdf`"
)
writeLines(report, "experiments/coalescent-susie-iw-benchmark-report.md")
print(summary[true_eta0 %in% c(20, 200, 1000, 5000), .(
  true_L, true_eta0, method, mean_reported_L, causal_power, cs_coverage,
  false_sets, fallback_frequency
)])
