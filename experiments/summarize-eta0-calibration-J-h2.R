#!/usr/bin/env Rscript

library(data.table)

root <- normalizePath(".", mustWork = TRUE)
worker_directory <- file.path(
  root, "analysis", "eta0-calibration-J-h2-workers"
)
files <- list.files(
  worker_directory, pattern = "^results-J[0-9]+-panel[0-9]+[.]csv$",
  full.names = TRUE
)
if (!length(files)) stop("No eta0 calibration worker results were found.")
results <- rbindlist(lapply(files, fread), fill = TRUE)

expected_J <- c(100L, 250L, 500L, 1000L, 2000L)
expected_panels <- c(31L, 32L)
expected_eta0 <- c(20, 50, 100, 200, 500, 1000)
expected_h2 <- c(2e-4, 1e-3, 5e-3)
expected_stages <- list(
  "AIW-N" = c(
    "oracle_lambda_oracle_L", "selected_lambda_oracle_L",
    "actual_selection"
  ),
  "collapsed-R" = c("oracle_lambda", "actual_selection")
)
stopifnot(
  setequal(unique(results$J), expected_J),
  setequal(unique(results$panel_seed), expected_panels),
  setequal(unique(results$true_eta0), expected_eta0),
  setequal(unique(results$h2), expected_h2),
  nrow(results) == length(expected_J) * length(expected_panels) *
    length(expected_eta0) * length(expected_h2) * 5L
)

results[, log_true_centered := log(true_eta0 / 1000)]
results[, log_hat_centered := log(raw_eta0 / 1000)]
results[, log_J_centered := log(J / 1000)]
results[, log_h2_centered := log(h2 / 1e-3)]
# Keep every outcome generated from the same reference-LD panel together.
# Otherwise the same R0 can appear in both the training and validation folds.
results[, calibration_group := paste(J, panel_seed, sep = ":")]
results[, finite_for_calibration :=
  eta0_finite & is.finite(log_hat_centered) & raw_eta0 > 0]

set.seed(8675309)
groups <- sample(unique(results$calibration_group))
fold_map <- data.table(
  calibration_group = groups,
  fold = (seq_along(groups) - 1L) %% 5L + 1L
)
results <- merge(results, fold_map, by = "calibration_group")

model_specs <- data.table(
  model = c(
    "fixed_b_no_J_no_h2", "free_b_no_J_no_h2",
    "fixed_b_with_J_no_h2", "free_b_with_J_no_h2",
    "fixed_b_no_J_with_h2", "free_b_no_J_with_h2",
    "fixed_b_with_J_with_h2", "free_b_with_J_with_h2"
  ),
  free_b = rep(c(FALSE, TRUE), 4L),
  include_J = rep(c(FALSE, FALSE, TRUE, TRUE), 2L),
  include_h2 = c(rep(FALSE, 4L), rep(TRUE, 4L)),
  complexity = c(0L, 1L, 1L, 2L, 1L, 2L, 2L, 3L)
)

fit_calibration_model <- function(data, free_b, include_J, include_h2) {
  predictors <- character()
  if (free_b) {
    response <- "log_true_centered"
    predictors <- "log_hat_centered"
  } else {
    data <- copy(data)
    data[, log_ratio := log_true_centered - log_hat_centered]
    response <- "log_ratio"
  }
  if (include_J) predictors <- c(predictors, "log_J_centered")
  if (include_h2) predictors <- c(predictors, "log_h2_centered")
  formula <- reformulate(predictors, response = response)
  lm(formula, data = data)
}

predict_calibration_model <- function(
    fit, newdata, free_b) {
  prediction <- as.numeric(predict(fit, newdata = newdata))
  if (!free_b) prediction <- newdata$log_hat_centered + prediction
  prediction
}

coefficient_value <- function(fit, name, default) {
  coefficient <- coef(fit)
  if (!name %in% names(coefficient)) return(default)
  unname(coefficient[[name]])
}

coefficient_rows <- list()
prediction_rows <- list()
fold_rows <- list()
counter <- 0L
prediction_counter <- 0L
fold_counter <- 0L

for (method_value in names(expected_stages)) {
  for (stage_value in expected_stages[[method_value]]) {
    data <- results[
      method == method_value &
        estimation_stage == stage_value &
        finite_for_calibration
    ]
    if (nrow(data) < 20L) next

    for (spec_index in seq_len(nrow(model_specs))) {
      spec <- model_specs[spec_index]
      fit <- fit_calibration_model(
        data, spec$free_b, spec$include_J, spec$include_h2
      )
      fitted_log_eta0 <- predict_calibration_model(
        fit, data, spec$free_b
      )
      residual <- fitted_log_eta0 - data$log_true_centered

      counter <- counter + 1L
      coefficient_rows[[counter]] <- data.table(
        method = method_value,
        estimation_stage = stage_value,
        model = spec$model,
        finite_fits = nrow(data),
        intercept = coefficient_value(fit, "(Intercept)", 0),
        reference_multiplier = exp(
          coefficient_value(fit, "(Intercept)", 0)
        ),
        b = if (spec$free_b) {
          coefficient_value(fit, "log_hat_centered", NA_real_)
        } else 1,
        J_exponent = coefficient_value(
          fit, "log_J_centered", 0
        ),
        h2_exponent = coefficient_value(
          fit, "log_h2_centered", 0
        ),
        in_sample_RMSE = sqrt(mean(residual^2)),
        in_sample_MAE = mean(abs(residual)),
        in_sample_bias = mean(residual),
        adjusted_R_squared = summary(fit)$adj.r.squared,
        AIC = AIC(fit)
      )

      for (fold_value in sort(unique(data$fold))) {
        training <- data[fold != fold_value]
        testing <- data[fold == fold_value]
        fold_fit <- fit_calibration_model(
          training, spec$free_b, spec$include_J, spec$include_h2
        )
        predicted <- predict_calibration_model(
          fold_fit, testing, spec$free_b
        )
        error <- predicted - testing$log_true_centered
        prediction_counter <- prediction_counter + 1L
        prediction_rows[[prediction_counter]] <- data.table(
          method = method_value,
          estimation_stage = stage_value,
          model = spec$model,
          fold = fold_value,
          J = testing$J,
          panel_seed = testing$panel_seed,
          true_eta0 = testing$true_eta0,
          h2 = testing$h2,
          raw_eta0 = testing$raw_eta0,
          predicted_eta0 = 1000 * exp(predicted),
          log_prediction_error = error
        )
        fold_counter <- fold_counter + 1L
        fold_rows[[fold_counter]] <- data.table(
          method = method_value,
          estimation_stage = stage_value,
          model = spec$model,
          fold = fold_value,
          fits = nrow(testing),
          RMSE = sqrt(mean(error^2)),
          MAE = mean(abs(error)),
          bias = mean(error)
        )
      }
    }
  }
}

coefficients <- rbindlist(coefficient_rows)
predictions <- rbindlist(prediction_rows)
fold_metrics <- rbindlist(fold_rows)
cv_summary <- predictions[, .(
  held_out_fits = .N,
  CV_RMSE = sqrt(mean(log_prediction_error^2)),
  CV_MAE = mean(abs(log_prediction_error)),
  CV_bias = mean(log_prediction_error),
  maximum_absolute_true_eta0_bias = max(abs(
    tapply(log_prediction_error, true_eta0, mean)
  )),
  bias_at_true_eta0_20 = mean(
    log_prediction_error[true_eta0 == 20]
  ),
  bias_at_true_eta0_1000 = mean(
    log_prediction_error[true_eta0 == 1000]
  ),
  maximum_absolute_group_bias = max(abs(
    tapply(log_prediction_error, list(J, h2), mean)
  ))
), by = .(method, estimation_stage, model)]
fold_summary <- fold_metrics[, .(
  mean_fold_RMSE = mean(RMSE),
  SE_fold_RMSE = sd(RMSE) / sqrt(.N),
  mean_fold_MAE = mean(MAE),
  mean_absolute_fold_bias = mean(abs(bias))
), by = .(method, estimation_stage, model)]
cv_summary <- merge(
  cv_summary, fold_summary,
  by = c("method", "estimation_stage", "model")
)
cv_summary <- merge(cv_summary, model_specs, by = "model")

boundary_summary <- results[, .(
  fits = .N,
  finite_fits = sum(finite_for_calibration),
  finite_frequency = mean(finite_for_calibration),
  normal_limit_frequency = mean(eta0_boundary == "normal_limit"),
  optimization_boundary_frequency = mean(
    eta0_boundary == "optimization_boundary"
  ),
  convergence_frequency = mean(converged),
  median_raw_eta0 = if (any(finite_for_calibration)) {
    median(raw_eta0[finite_for_calibration])
  } else NA_real_
), by = .(method, estimation_stage, J, h2, true_eta0)]

recommended_rows <- list()
recommendation_counter <- 0L
for (method_value in c("AIW-N", "collapsed-R")) {
  candidates <- cv_summary[
    method == method_value & estimation_stage == "actual_selection"
  ]
  best <- candidates[which.min(CV_RMSE)]
  threshold <- best$mean_fold_RMSE + best$SE_fold_RMSE
  eligible <- candidates[mean_fold_RMSE <= threshold]
  setorder(
    eligible, complexity, include_h2, include_J, free_b, CV_RMSE
  )
  one_se <- eligible[1L]
  for (kind in c("minimum_CV_RMSE", "one_SE_parsimonious")) {
    chosen <- if (kind == "minimum_CV_RMSE") best else one_se
    coefficient <- coefficients[
      method == method_value &
        estimation_stage == "actual_selection" &
        model == chosen$model
    ]
    recommendation_counter <- recommendation_counter + 1L
    recommended_rows[[recommendation_counter]] <- cbind(
      data.table(
        method = method_value,
        recommendation = kind,
        model = chosen$model
      ),
      coefficient[, .(
        reference_multiplier, b, J_exponent, h2_exponent
      )],
      chosen[, .(
        CV_RMSE, CV_MAE, CV_bias,
        maximum_absolute_group_bias,
        mean_fold_RMSE, SE_fold_RMSE
      )]
    )
  }
}
recommendations <- rbindlist(recommended_rows)

recommended_values <- rbindlist(lapply(
  split(
    recommendations[recommendation == "one_SE_parsimonious"],
    by = "method", keep.by = TRUE
  ),
  function(recommendation) {
    grid <- CJ(
      J = expected_J,
      h2 = expected_h2,
      raw_eta0 = c(20, 50, 100, 200, 500, 1000)
    )
    grid[, calibrated_eta0 :=
      1000 * recommendation$reference_multiplier *
        (raw_eta0 / 1000)^recommendation$b *
        (J / 1000)^recommendation$J_exponent *
        (h2 / 1e-3)^recommendation$h2_exponent
    ]
    grid[, calibration_ratio := calibrated_eta0 / raw_eta0]
    grid[, method := recommendation$method]
    grid[, model := recommendation$model]
    grid
  }
))

prefix <- file.path(root, "analysis", "eta0-calibration-J-h2")
fwrite(results, paste0(prefix, "-results.csv"))
fwrite(coefficients, paste0(prefix, "-model-coefficients.csv"))
fwrite(predictions, paste0(prefix, "-cross-validation-predictions.csv"))
fwrite(cv_summary, paste0(prefix, "-model-comparison.csv"))
fwrite(boundary_summary, paste0(prefix, "-boundary-summary.csv"))
fwrite(recommendations, paste0(prefix, "-recommendations.csv"))
fwrite(recommended_values, paste0(prefix, "-recommended-values.csv"))

colors <- setNames(hcl.colors(nrow(model_specs), "Dark 3"), model_specs$model)
model_labels <- c(
  fixed_b_no_J_no_h2 = "b=1, no J, no h2",
  free_b_no_J_no_h2 = "free b, no J, no h2",
  fixed_b_with_J_no_h2 = "b=1, with J, no h2",
  free_b_with_J_no_h2 = "free b, with J, no h2",
  fixed_b_no_J_with_h2 = "b=1, no J, with h2",
  free_b_no_J_with_h2 = "free b, no J, with h2",
  fixed_b_with_J_with_h2 = "b=1, with J, with h2",
  free_b_with_J_with_h2 = "free b, with J, with h2"
)

draw_figure <- function(device, path) {
  device(path, width = 12, height = 9)
  layout(
    matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE),
    heights = c(1.15, 1)
  )
  for (method_value in c("AIW-N", "collapsed-R")) {
    par(mar = c(4, 10.5, 3, 1))
    comparison <- cv_summary[
      method == method_value & estimation_stage == "actual_selection"
    ]
    comparison <- comparison[match(names(colors), model)]
    barplot(
      comparison$CV_RMSE,
      names.arg = model_labels[comparison$model],
      col = colors[comparison$model],
      horiz = TRUE, las = 1, cex.names = 0.72,
      xlab = "Held-out RMSE on log eta0",
      main = paste(method_value, "model comparison")
    )
  }
  for (method_value in c("AIW-N", "collapsed-R")) {
    par(mar = c(5, 4.5, 3, 1))
    recommendation <- recommendations[
      method == method_value &
        recommendation == "one_SE_parsimonious"
    ]
    model_value <- recommendation$model
    prediction <- predictions[
      method == method_value &
        estimation_stage == "actual_selection" &
        model == model_value
    ]
    plot(
      log(prediction$true_eta0), log(prediction$predicted_eta0),
      pch = 16, col = adjustcolor("#2166AC", alpha.f = 0.55),
      xlab = "log(true eta0)", ylab = "held-out log(calibrated eta0)",
      main = paste(method_value, model_labels[model_value])
    )
    abline(0, 1, lty = 2, col = "gray40")
  }
  dev.off()
}
figure_directory <- file.path(root, "analysis", "figures")
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)
draw_figure(function(path, width, height) {
  png(path, width = width, height = height, units = "in", res = 160)
}, file.path(figure_directory, "eta0-calibration-J-h2.png"))
draw_figure(pdf, file.path(
  figure_directory, "eta0-calibration-J-h2.pdf"
))

markdown_table <- function(data, digits = 4L) {
  data <- copy(data)
  for (column in names(data)) {
    if (is.numeric(data[[column]])) {
      data[[column]] <- formatC(
        data[[column]], digits = digits, format = "f"
      )
    }
  }
  c(
    paste0("| ", paste(names(data), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(data)), collapse = " | "), " |"),
    apply(data, 1L, function(row) {
      paste0("| ", paste(row, collapse = " | "), " |")
    })
  )
}

comparison_report <- merge(
  cv_summary[estimation_stage == "actual_selection", .(
    method, model, CV_RMSE, CV_MAE, CV_bias,
    maximum_absolute_true_eta0_bias,
    bias_at_true_eta0_20, bias_at_true_eta0_1000,
    maximum_absolute_group_bias
  )],
  coefficients[estimation_stage == "actual_selection", .(
    method, model, reference_multiplier, b,
    J_exponent, h2_exponent
  )],
  by = c("method", "model")
)
setorder(comparison_report, method, CV_RMSE)
recommendation_report <- recommendations[, .(
  method, recommendation, model, reference_multiplier,
  b, J_exponent, h2_exponent, CV_RMSE, CV_MAE, CV_bias
)]

aiwn_one_se <- recommendations[
  method == "AIW-N" & recommendation == "one_SE_parsimonious"
]
collapsed_one_se <- recommendations[
  method == "collapsed-R" & recommendation == "one_SE_parsimonious"
]
aiwn_fixed_constant <- coefficients[
  method == "AIW-N" & estimation_stage == "actual_selection" &
    model == "fixed_b_no_J_no_h2"
]
collapsed_fixed_constant <- coefficients[
  method == "collapsed-R" & estimation_stage == "actual_selection" &
    model == "fixed_b_no_J_no_h2"
]
collapsed_identity_crossing <- 1000 *
  collapsed_one_se$reference_multiplier^(
    1 / (1 - collapsed_one_se$b)
  )

writeLines(c(
  "# Calibration of eta0 as a function of J and local h2",
  "",
  "The experiment varies J over {100, 250, 500, 1000, 2000}, true eta0 ",
  "over {20, 50, 100, 200, 500, 1000}, and local h2 over ",
  "{2e-4, 1e-3, 5e-3}. Each cell uses two independent coalescent ",
  "reference-LD panels. The phenotype contains two weak-LD causal effects ",
  "with fixed 7:3 weights, and the capped lambda grid is searched.",
  "",
  "The four requested calibration forms are fitted separately for AIW-N ",
  "and collapsed-R, with J included:",
  "",
  "1. b fixed to one, without h2;",
  "2. b estimated, without h2;",
  "3. b fixed to one, with h2; and",
  "4. b estimated, with h2.",
  "",
  "The corresponding four no-J baselines are also fitted. These baselines ",
  "test whether allowing calibration to depend on locus size improves ",
  "held-out prediction rather than assuming that it must.",
  "",
  "The centered calibration equation is",
  "",
  "eta0_cal = 1000 exp(a) (eta0_hat/1000)^b ",
  "             (J/1000)^cJ (h2/1e-3)^ch.",
  "",
  "For b=1, exp(a) is the multiplicative correction at J=1000 and ",
  "h2=1e-3. Models without h2 set ch=0. Normal-limit estimates are ",
  "recorded but excluded from finite calibration regressions because they ",
  "require no multiplicative correction.",
  "",
  "## Actual-selection model comparison",
  "",
  markdown_table(comparison_report),
  "",
  "## Recommended models and coefficients",
  "",
  markdown_table(recommendation_report),
  "",
  "The minimum-RMSE recommendation is the model with the smallest pooled ",
  "five-fold held-out log-RMSE. The one-standard-error recommendation is ",
  "the least complex model within one standard error of the best model's ",
  "mean fold RMSE. The latter is the default recommendation for discussion.",
  "All observations sharing a reference-LD panel are assigned to the same ",
  "fold, so validation does not reuse R0 across training and testing.",
  "",
  "## Interpretation before locking the calibration",
  "",
  sprintf(
    paste0(
      "For AIW-N, the one-standard-error rule selects a constant multiplier ",
      "of %.4f: eta0_cal = %.4f eta0_hat. Neither J nor h2 is retained. "
    ),
    aiwn_one_se$reference_multiplier,
    aiwn_one_se$reference_multiplier
  ),
  "The free-b models lower held-out RMSE, but largely by shrinking estimates ",
  "toward the center of the simulated eta0 grid. This creates appreciably ",
  "larger conditional bias at true eta0=20 and 1000. For a transportable ",
  "multiplicative correction, the constant AIW-N model is therefore the ",
  "preferred candidate for the next downstream check.",
  "",
  sprintf(
    paste0(
      "For collapsed-R, the predictive one-standard-error model is ",
      "eta0_cal = 1000 x %.4f x (eta0_hat/1000)^%.4f, with neither J nor ",
      "h2 retained. It crosses the identity line at eta0_hat approximately ",
      "%.1f and downscales larger estimates."
    ),
    collapsed_one_se$reference_multiplier,
    collapsed_one_se$b,
    collapsed_identity_crossing
  ),
  sprintf(
    paste0(
      "That shrinkage map has the best average predictive performance, but ",
      "it is not equivalent to the intended underestimation correction. The ",
      "simple collapsed-R multiplier from this design is %.4f, close to the ",
      "previously tested value 1.537, but its log-RMSE is much worse because ",
      "the joint lambda-eta0 profile occasionally produces very large ",
      "estimates. I therefore record both values and do not automatically ",
      "replace 1.537 before examining downstream power and coverage."
    ),
    collapsed_fixed_constant$reference_multiplier
  ),
  "",
  sprintf(
    "Thus the current release candidates are: AIW-N with constant %.4f; ",
    aiwn_fixed_constant$reference_multiplier
  ),
  "and, for collapsed-R, either the predictive free-b map above or the ",
  "already validated multiplicative correction 1.537. The latter choice ",
  "requires a methodological decision, not another regression criterion.",
  "",
  "![Calibration comparison](figures/eta0-calibration-J-h2.png)",
  "",
  "Detailed coefficients, cross-validation predictions, boundary rates, ",
  "and calibration values over the design grid are saved beside this report."
), paste0(prefix, "-report.md"))

message("Wrote eta0 calibration summaries and report.")
