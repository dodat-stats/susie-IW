#!/usr/bin/env Rscript

library(data.table)
source(file.path("experiments", "simulation-iw-collapsed-R.R"))
source(file.path("code", "01-collapsed-R.R"))
source(file.path("code", "02-AIW-and-AIW-N.R"))
source(file.path("code", "03-SER-fallback-and-plot.R"))

result_file <- "experiments/coalescent-susie-iw-benchmark-results.csv"
if (!file.exists(result_file)) stop("Run and summarize the benchmark first")
results <- fread(result_file)

wide <- dcast(
  results[true_L == 2L & J == 1000L],
  J + panel_seed + true_eta0 + true_L ~ method,
  value.var = c("reported_L", "false_sets", "causal_coverage")
)
preferred <- wide[
  get("reported_L_SuSiE-RSS: in-sample R") == 2L &
    get("false_sets_SuSiE-RSS: in-sample R") == 0L &
    get("causal_coverage_SuSiE-RSS: in-sample R") == 1 &
    get("reported_L_SuSiE-RSS: reference R0") >= 3L &
    get("false_sets_SuSiE-RSS: reference R0") >= 1L &
    get("reported_L_collapsed-R") <= 2L &
    get("false_sets_collapsed-R") == 0L &
    get("causal_coverage_collapsed-R") == 1 &
    get("reported_L_AIW-N") <= 2L &
    get("false_sets_AIW-N") == 0L &
    get("causal_coverage_AIW-N") == 1
]
if (!nrow(preferred)) stop("No J=1000 configuration met the motivating-example rule")

preferred[, priority := abs(log(true_eta0 / 200))]
setorder(preferred, priority, panel_seed)
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments)) {
  requested_panel_seed <- as.integer(arguments[1L])
  chosen <- preferred[panel_seed == requested_panel_seed]
  if (!nrow(chosen)) {
    stop("The requested panel seed does not meet the motivating-example rule")
  }
  chosen <- chosen[1L]
} else {
  chosen <- preferred[1L]
}

J <- 1000L
panel_seed <- chosen$panel_seed
true_eta0 <- chosen$true_eta0
true_L <- 2L
N <- 200000L
N0 <- 2000L
h2 <- 1e-3
true_lambda <- 1e-3
maximum_L <- 5L
lambda_grid <- c(1e-5, 1e-4, 1e-3, 2e-3, 4e-3, 6e-3, 8e-3, 1e-2)

R0 <- simulate_coalescent_reference_ld(
  coalescent_sample_size = 2000L,
  reference_sample_size = N0,
  number_of_variants = J,
  seed = panel_seed,
  heritability = h2
)
Rbar_true <- make_stabilized_reference_ld(R0, true_lambda)
simulation_seed <- as.integer(2000000 + 10000 * panel_seed + J + true_eta0)
set.seed(simulation_seed)
R <- standardize_correlation(sample_inverse_wishart(
  degrees_of_freedom = true_eta0 + J + 1,
  scale_matrix = true_eta0 * Rbar_true
))

candidates <- unique(as.integer(round(seq(
  max(2, 0.08 * J), min(J - 1, 0.92 * J),
  length.out = min(250L, J - 2L)
))))
minimum_separation <- max(8L, floor(J / 12))
first <- candidates[which.min(abs(candidates - round(J / 4)))]
available <- candidates[abs(candidates - first) >= minimum_separation]
second <- available[which.min(abs(R[first, available]))]
causal <- c(first, second)

weights <- c(7, 5)
beta <- numeric(J)
beta[causal] <- weights * sqrt(
  h2 / drop(crossprod(weights, R[causal, causal, drop = FALSE] %*% weights))
)
set.seed(simulation_seed + 500000L)
noise <- as.numeric(t(stable_cholesky(R)) %*% rnorm(J)) / sqrt(N)
v <- as.numeric(R %*% beta + noise)

fits <- list(
  "SuSiE-RSS: in-sample LD" = fit_susie_rss(v, R, N, L = maximum_L),
  "SuSiE-RSS: reference LD" = fit_susie_rss(v, R0, N, L = maximum_L),
  "Collapsed-R" = susie_iw(
    v = v, R0 = R0, N = N, L = maximum_L,
    lambda_grid = lambda_grid, eta0_multiplier = 1.677,
    verbose = FALSE
  ),
  "AIW-N" = suppressMessages(susie_aiw(
    v = v, R0 = R0, N = N, L = maximum_L,
    approximation = "N", lambda_grid = lambda_grid,
    eta0_multiplier = 1.424, verbose = FALSE
  ))
)
plot_items <- list(
  list(object = fits[[1L]], R = R),
  list(object = fits[[2L]], R = R0),
  list(object = fits[[3L]], R = fits[[3L]]$selected_Rbar),
  list(object = fits[[4L]], R = fits[[4L]]$selected_Rbar)
)
names(plot_items) <- names(fits)
output_stem <- sprintf(
  "experiments/figures/coalescent-susie-iw-motivating-locus-J1000-seed%d",
  panel_seed
)

draw <- function() {
  plot_iw_comparison(
    plot_items, causal = causal, ncol = 2L,
    main = sprintf(
      "Coalescent SuSiE-IW example: J=%d, true eta0=%g, true L=2",
      J, true_eta0
    ),
    point_cex = 0.75,
    credible_set_cex = 2.15,
    credible_set_lwd = 3.25,
    legend_cex = 0.82,
    legend_point_cex = 1.55,
    causal_cex = 1.4
  )
}
png(paste0(output_stem, ".png"), width = 1900, height = 1600, res = 180)
draw()
dev.off()
pdf(paste0(output_stem, ".pdf"), width = 10.5, height = 8.9)
draw()
dev.off()

summary <- rbindlist(lapply(seq_along(fits), function(i) {
  row <- as.data.table(summarize_iw_fit(fits[[i]], causal = causal))
  credible_sets <- susie_get_reported_cs(fits[[i]])
  set_size_values <- lengths(credible_sets)
  pip <- iw_pip(fits[[i]])
  row[, method := names(fits)[i]]
  row[, credible_set_sizes := paste(set_size_values, collapse = ";")]
  row[, mean_credible_set_size := if (length(set_size_values)) {
    mean(set_size_values)
  } else {
    NA_real_
  }]
  row[, causal_pips := paste(sprintf("%.4f", pip[causal]), collapse = ";")]
  row
}), fill = TRUE)
fwrite(
  summary,
  sprintf(
    "experiments/coalescent-susie-iw-motivating-locus-J1000-seed%d-summary.csv",
    panel_seed
  )
)
writeLines(c(
  sprintf("J: %d", J),
  sprintf("panel seed: %d", panel_seed),
  sprintf("simulation seed: %d", simulation_seed),
  sprintf("true eta0: %g", true_eta0),
  sprintf("causal variants: %s", paste(causal, collapse = ", ")),
  sprintf("causal LD: %.6f", R[causal[1L], causal[2L]])
), sprintf(
  "experiments/coalescent-susie-iw-motivating-locus-J1000-seed%d-settings.txt",
  panel_seed
))
print(summary)
