#!/usr/bin/env Rscript

library(data.table)

worker_directory <- file.path(
  "analysis", "eta0-calibration-coalescent-L-workers"
)
files <- list.files(
  worker_directory, pattern = "^results-J[0-9]+-panel[0-9]+[.]csv$",
  full.names = TRUE
)
if (!length(files)) stop("No coalescent L-calibration results were found")
results <- rbindlist(lapply(files, fread), fill = TRUE)
setorder(results, method, J, panel_seed, true_eta0, true_L)
# A collapsed-R profile at the upper optimization boundary represents the
# matched normal limit, just as an infinite AIW-N profile does. Some worker
# output predates the explicit boundary label, so values within 10% of the
# configured upper bound (1e8) are classified here as normal-limit fits.
results[method == "collapsed-R" & raw_eta0 >= 0.9e8,
        eta0_boundary := "normal_limit"]
results[, finite_eta0 := is.finite(raw_eta0) & raw_eta0 > 0 &
  eta0_boundary != "normal_limit"]
results[, log_ratio := log(true_eta0) - log(raw_eta0)]
results[, log_J_centered := log(J / 500)]
results[, calibration_group := paste(J, panel_seed, sep = ":")]

expected_J <- c(100L, 250L, 500L, 1000L, 2000L)
expected_panels <- c(51L, 52L)
expected_eta0 <- c(20, 50, 100, 200, 500, 1000, 2000, 5000)
stopifnot(
  setequal(results$J, expected_J),
  setequal(results$panel_seed, expected_panels),
  setequal(results$true_eta0, expected_eta0),
  setequal(results$true_L, 1:3),
  nrow(results) == length(expected_J) * length(expected_panels) *
    length(expected_eta0) * 3L * 2L
)

model_definitions <- list(
  constant = log_ratio ~ 1,
  J = log_ratio ~ log_J_centered,
  L = log_ratio ~ factor(true_L),
  J_and_L = log_ratio ~ log_J_centered + factor(true_L)
)

coefficient_value <- function(fit, name, default = NA_real_) {
  values <- coef(fit)
  if (!name %in% names(values)) return(default)
  unname(values[[name]])
}

fit_rows <- list()
cv_rows <- list()
fit_counter <- 0L
cv_counter <- 0L
for (method_value in c("AIW-N", "collapsed-R")) {
  data <- results[method == method_value & finite_eta0]
  for (model_name in names(model_definitions)) {
    fit <- lm(model_definitions[[model_name]], data = data)
    interval <- suppressMessages(confint(fit, "(Intercept)", level = 0.95))
    fit_counter <- fit_counter + 1L
    fit_rows[[fit_counter]] <- data.table(
      method = method_value,
      model = model_name,
      finite_fits = nrow(data),
      intercept = coefficient_value(fit, "(Intercept)"),
      multiplier_at_J_500_L_1 = exp(coefficient_value(fit, "(Intercept)")),
      multiplier_CI_lower = exp(interval[1L]),
      multiplier_CI_upper = exp(interval[2L]),
      J_exponent = coefficient_value(fit, "log_J_centered", 0),
      J_exponent_SE = if ("log_J_centered" %in% rownames(summary(fit)$coef)) {
        summary(fit)$coef["log_J_centered", "Std. Error"]
      } else 0,
      L2_log_shift = coefficient_value(fit, "factor(true_L)2", 0),
      L3_log_shift = coefficient_value(fit, "factor(true_L)3", 0),
      adjusted_R_squared = summary(fit)$adj.r.squared,
      AIC = AIC(fit)
    )

    groups <- unique(data$calibration_group)
    for (group in groups) {
      training <- data[calibration_group != group]
      testing <- data[calibration_group == group]
      fold_fit <- lm(model_definitions[[model_name]], data = training)
      prediction <- log(testing$raw_eta0) + predict(
        fold_fit, newdata = testing
      )
      error <- prediction - log(testing$true_eta0)
      cv_counter <- cv_counter + 1L
      cv_rows[[cv_counter]] <- data.table(
        method = method_value, model = model_name,
        held_out_group = group, observations = nrow(testing),
        squared_error_sum = sum(error^2),
        absolute_error_sum = sum(abs(error)),
        error_sum = sum(error)
      )
    }
  }
}
model_fits <- rbindlist(fit_rows)
cv_folds <- rbindlist(cv_rows)
cv_summary <- cv_folds[, .(
  held_out_fits = sum(observations),
  CV_RMSE = sqrt(sum(squared_error_sum) / sum(observations)),
  CV_MAE = sum(absolute_error_sum) / sum(observations),
  CV_bias = sum(error_sum) / sum(observations)
), by = .(method, model)]
model_comparison <- merge(model_fits, cv_summary, by = c("method", "model"))
setorder(model_comparison, method, CV_RMSE)

by_L <- results[finite_eta0 == TRUE, .(
  finite_fits = .N,
  multiplier = exp(mean(log_ratio)),
  mean_log_ratio = mean(log_ratio),
  SE_log_ratio = sd(log_ratio) / sqrt(.N),
  RMSE_after_stratum_multiplier = sqrt(mean(
    (log(raw_eta0) + mean(log_ratio) - log(true_eta0))^2
  ))
), by = .(method, true_L)]
by_L[, `:=`(
  multiplier_CI_lower = exp(mean_log_ratio - 1.96 * SE_log_ratio),
  multiplier_CI_upper = exp(mean_log_ratio + 1.96 * SE_log_ratio)
)]

by_J <- results[finite_eta0 == TRUE, .(
  finite_fits = .N,
  multiplier = exp(mean(log_ratio)),
  mean_log_ratio = mean(log_ratio),
  SE_log_ratio = sd(log_ratio) / sqrt(.N)
), by = .(method, J)]
by_J[, `:=`(
  multiplier_CI_lower = exp(mean_log_ratio - 1.96 * SE_log_ratio),
  multiplier_CI_upper = exp(mean_log_ratio + 1.96 * SE_log_ratio)
)]

boundary_summary <- results[, .(
  fits = .N,
  finite_fits = sum(finite_eta0),
  finite_frequency = mean(finite_eta0),
  convergence_frequency = mean(converged),
  normal_limit_frequency = mean(eta0_boundary == "normal_limit"),
  optimization_boundary_frequency = mean(
    eta0_boundary == "optimization_boundary"
  ),
  selected_L_mean = mean(selected_model_L),
  selected_L_correct = mean(selected_model_L == true_L)
), by = .(method, true_L)]

recommendations <- model_comparison[model == "constant", .(
  method,
  recommended_multiplier = multiplier_at_J_500_L_1,
  multiplier_CI_lower,
  multiplier_CI_upper,
  finite_fits,
  CV_RMSE,
  CV_MAE,
  CV_bias
)]
recommendations <- merge(
  recommendations,
  model_comparison[model == "J", .(
    method, J_exponent, J_exponent_SE,
    J_model_CV_RMSE = CV_RMSE,
    J_model_CV_MAE = CV_MAE
  )], by = "method"
)
recommendations[, J_exponent_CI_lower := J_exponent - 1.96 * J_exponent_SE]
recommendations[, J_exponent_CI_upper := J_exponent + 1.96 * J_exponent_SE]
recommendations[, J_CV_RMSE_change := J_model_CV_RMSE - CV_RMSE]

fwrite(results, "experiments/eta0-calibration-coalescent-L-results.csv")
fwrite(model_comparison,
       "experiments/eta0-calibration-coalescent-L-model-comparison.csv")
fwrite(by_L, "experiments/eta0-calibration-coalescent-L-by-L.csv")
fwrite(by_J, "experiments/eta0-calibration-coalescent-L-by-J.csv")
fwrite(boundary_summary,
       "experiments/eta0-calibration-coalescent-L-boundary-summary.csv")
fwrite(recommendations,
       "experiments/eta0-calibration-coalescent-L-recommendations.csv")

figure_png <- "experiments/figures/eta0-calibration-coalescent-L.png"
figure_pdf <- "experiments/figures/eta0-calibration-coalescent-L.pdf"
dir.create(dirname(figure_png), recursive = TRUE, showWarnings = FALSE)

draw_figure <- function() {
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  par(mfrow = c(2, 2), mar = c(4.2, 4.4, 2.6, 1), las = 1)
  colors <- c("#1f78b4", "#33a02c", "#e31a1c")
  shapes <- c(16, 17, 15)
  for (method_value in c("AIW-N", "collapsed-R")) {
    data <- results[method == method_value & finite_eta0]
    recommendation <- recommendations[method == method_value]
    plot(
      data$true_eta0, data$raw_eta0,
      log = "xy", col = colors[data$true_L], pch = shapes[data$true_L],
      xlab = expression("true " * eta[0]),
      ylab = expression("raw estimated " * hat(eta)[0]),
      main = method_value, cex = 0.7
    )
    abline(0, 1, lty = 2, col = "grey45")
    abline(
      a = -log(recommendation$recommended_multiplier), b = 1,
      col = "black", lwd = 2
    )
    legend(
      "topleft", legend = paste("true L =", 1:3),
      col = colors, pch = shapes, bty = "n", cex = 0.75
    )
  }
  for (method_value in c("AIW-N", "collapsed-R")) {
    data <- by_J[method == method_value]
    plot(
      data$J, data$multiplier, type = "b", log = "x", pch = 16,
      ylim = range(c(data$multiplier_CI_lower, data$multiplier_CI_upper)),
      xlab = "number of variants J", ylab = "fixed-slope multiplier",
      main = paste(method_value, "by locus size")
    )
    arrows(
      data$J, data$multiplier_CI_lower,
      data$J, data$multiplier_CI_upper,
      angle = 90, code = 3, length = 0.04
    )
    abline(
      h = recommendations[method == method_value]$recommended_multiplier,
      lty = 2, col = "grey35"
    )
  }
}

grDevices::png(figure_png, width = 1800, height = 1500, res = 180)
draw_figure()
grDevices::dev.off()
grDevices::pdf(figure_pdf, width = 10, height = 8.3)
draw_figure()
grDevices::dev.off()

format_row <- function(row) {
  sprintf(
    "| %s | %.4f | %.4f--%.4f | %.4f | %.4f | %.4f--%.4f |",
    row$method, row$recommended_multiplier,
    row$multiplier_CI_lower, row$multiplier_CI_upper,
    row$CV_RMSE, row$J_exponent, row$J_exponent_CI_lower,
    row$J_exponent_CI_upper
  )
}
report <- c(
  "# Coalescent SuSiE-IW calibration across true effect counts",
  "",
  paste0(
    "The design contains ", uniqueN(results[, .(J, panel_seed, true_eta0, true_L)]),
    " datasets: J in {", paste(expected_J, collapse = ", "),
    "}, two independent coalescent reference panels per J, true eta0 in {",
    paste(expected_eta0, collapse = ", "),
    "}, and true L in {1, 2, 3}. The generating model fixes h2=1e-3, ",
    "N=200,000, N0=2,000, and lambda_true=1e-3."
  ),
  "",
  "The primary model fixes the coefficient of log(raw eta0) to one:",
  "",
  "log(true eta0) = log(multiplier) + log(raw estimated eta0).",
  "",
  "| method | multiplier | 95% CI | CV RMSE | J exponent | J exponent 95% CI |",
  "| --- | ---: | ---: | ---: | ---: | ---: |",
  vapply(seq_len(nrow(recommendations)), function(i) {
    format_row(recommendations[i])
  }, character(1L)),
  "",
  "The J model adds gamma log(J/500). Its cross-validated RMSE change relative to the constant model is:",
  "",
  paste0(
    "- ", recommendations$method, ": ",
    sprintf("%+.4f", recommendations$J_CV_RMSE_change)
  ),
  "",
  "Multipliers stratified by true L:",
  "",
  "| method | true L | multiplier | 95% CI | finite fits |",
  "| --- | ---: | ---: | ---: | ---: |",
  vapply(seq_len(nrow(by_L)), function(i) sprintf(
    "| %s | %d | %.4f | %.4f--%.4f | %d |",
    by_L$method[i], by_L$true_L[i], by_L$multiplier[i],
    by_L$multiplier_CI_lower[i], by_L$multiplier_CI_upper[i],
    by_L$finite_fits[i]
  ), character(1L)),
  "",
  "![Calibration across true L](figures/eta0-calibration-coalescent-L.png)"
)
writeLines(report, "experiments/eta0-calibration-coalescent-L-report.md")
print(recommendations)
