# Maintained AIW and AIW-N implementations for the SuSiE-IW project.
#
# Public functions:
#   susie_aiw(v, R0, N, approximation = c("N", "t"), ...)
#   fit_aiwn(v, R0, N, ...)  # legacy standardized-input interface
#   fit_aiw(v, R0, N, ...)   # legacy standardized-input interface
#
# Both methods always estimate eta0. eta0_multiplier changes each candidate's
# profiled eta0 before the candidate objectives are compared; it does not fix
# eta0. AIW-N defaults to the frozen calibration constant 1.424. Exact AIW
# has not been separately calibrated, so its default multiplier is 1.
#
# This file requires only the CRAN package susieR.

aiwn_solve_chol <- function(chol_A, b) {
  backsolve(chol_A, forwardsolve(t(chol_A), b))
}

aiwn_logdet_chol <- function(chol_A) {
  2 * sum(log(diag(chol_A)))
}

# susieR <= 0.14 exports susie_suff_stat(); susieR >= 0.16 exports the same
# sufficient-statistics interface as susie_ss(). Support both so AIW-N can be
# compared with the newer susie_rss LD-mismatch implementation.
aiwn_susie_ss_function <- function() {
  exports <- getNamespaceExports("susieR")
  if ("susie_suff_stat" %in% exports) {
    return(getExportedValue("susieR", "susie_suff_stat"))
  }
  if ("susie_ss" %in% exports) {
    return(getExportedValue("susieR", "susie_ss"))
  }
  stop("susieR does not export susie_suff_stat() or susie_ss()")
}

aiwn_dominant_distinct_configuration <- function(alpha) {
  alpha <- as.matrix(alpha)
  L <- nrow(alpha)
  J <- ncol(alpha)

  if (L < 1L || L > J) {
    stop("alpha must have between 1 and J rows")
  }
  if (any(!is.finite(alpha)) || any(alpha < 0)) {
    stop("alpha must be finite and nonnegative")
  }

  log_alpha <- log(alpha)
  n_states <- bitwShiftL(1L, L)
  scores <- rep(-Inf, n_states)
  paths <- vector("list", n_states)
  scores[1L] <- 0
  paths[[1L]] <- integer(L)

  for (j in seq_len(J)) {
    next_scores <- scores
    next_paths <- paths

    for (mask in 0:(n_states - 1L)) {
      state <- mask + 1L
      if (!is.finite(scores[state])) {
        next
      }

      for (ell in seq_len(L)) {
        bit <- bitwShiftL(1L, ell - 1L)
        if (bitwAnd(mask, bit) != 0L || !is.finite(log_alpha[ell, j])) {
          next
        }

        next_mask <- bitwOr(mask, bit)
        next_state <- next_mask + 1L
        candidate <- scores[state] + log_alpha[ell, j]

        if (candidate > next_scores[next_state]) {
          path <- paths[[state]]
          path[ell] <- j
          next_scores[next_state] <- candidate
          next_paths[[next_state]] <- path
        }
      }
    }

    scores <- next_scores
    paths <- next_paths
  }

  gamma <- paths[[n_states]]
  if (is.null(gamma) || any(gamma == 0L) || anyDuplicated(gamma)) {
    stop("failed to find a distinct configuration")
  }

  list(gamma = gamma, log_probability = scores[n_states])
}

aiwn_score_fit <- function(fit, v, Rbar, N, chol_Rbar = NULL,
                            logdet_Rbar = NULL, full_quadratic = NULL) {
  v <- as.numeric(v)
  Rbar <- as.matrix(Rbar)
  J <- length(v)
  L <- nrow(as.matrix(fit$alpha))

  if (!all(dim(Rbar) == c(J, J))) {
    stop("Rbar must be a J by J matrix")
  }
  if (L >= J) {
    stop("AIW-N requires L < J")
  }

  if (is.null(chol_Rbar)) {
    chol_Rbar <- chol(Rbar)
  }
  if (is.null(logdet_Rbar)) {
    logdet_Rbar <- aiwn_logdet_chol(chol_Rbar)
  }
  if (is.null(full_quadratic)) {
    full_quadratic <- drop(crossprod(
      v, aiwn_solve_chol(chol_Rbar, v)
    ))
  }

  config <- aiwn_dominant_distinct_configuration(fit$alpha)
  gamma <- config$gamma
  D_gamma <- Rbar[gamma, gamma, drop = FALSE]
  chol_D <- chol(D_gamma)
  logdet_D <- aiwn_logdet_chol(chol_D)
  v_gamma <- v[gamma]

  A_gamma <- drop(crossprod(
    v_gamma, aiwn_solve_chol(chol_D, v_gamma)
  ))
  r_gamma <- max(full_quadratic - A_gamma, 0)
  logdet_S_gamma <- logdet_Rbar - logdet_D

  index <- cbind(seq_len(L), gamma)
  beta_mean <- as.numeric(fit$mu[index])
  beta_variance <- pmax(
    as.numeric(fit$mu2[index]) - beta_mean^2,
    0
  )
  active_residual <- v_gamma - drop(D_gamma %*% beta_mean)
  active_quadratic <- drop(crossprod(
    active_residual,
    aiwn_solve_chol(chol_D, active_residual)
  ))
  active_variance <- sum(diag(D_gamma) * beta_variance)

  active_loglik <- -L / 2 * log(2 * pi) +
    L / 2 * log(N) -
    logdet_D / 2 -
    N / 2 * (active_quadratic + active_variance)

  if (is.null(fit$KL) || any(!is.finite(fit$KL))) {
    stop("the SuSiE fit does not contain a finite per-effect KL")
  }
  kl <- sum(pmax(as.numeric(fit$KL), 0))
  fixed_term <- active_loglik - kl

  average_inactive_residual <- r_gamma / (J - L)
  h_hat <- max(average_inactive_residual, 1 / N)
  eta_hat <- if (A_gamma <= sqrt(.Machine$double.eps)) {
    NA_real_
  } else if (average_inactive_residual <= 1 / N) {
    Inf
  } else {
    A_gamma / (average_inactive_residual - 1 / N)
  }

  inactive_loglik <- -(J - L) / 2 * log(2 * pi) -
    logdet_S_gamma / 2 -
    (J - L) / 2 * log(h_hat) -
    r_gamma / (2 * h_hat)

  data.frame(
    L = L,
    inactive_dimension = J - L,
    objective = fixed_term + inactive_loglik,
    eta0 = eta_hat,
    h_hat = h_hat,
    A_gamma = A_gamma,
    r_gamma = r_gamma,
    residual_per_snp = average_inactive_residual,
    logdet_S_gamma = logdet_S_gamma,
    active_loglik = active_loglik,
    inactive_loglik = inactive_loglik,
    KL = kl,
    gamma_hat = paste(gamma, collapse = ","),
    distinct_log_probability = config$log_probability,
    stringsAsFactors = FALSE
  )
}

# Apply a multiplicative bias correction to the profiled eta0 before using the
# inactive block to compare (lambda, L) models. A multiplier of one preserves
# the original profiled objective exactly.
aiwn_rescore_with_eta0_multiplier <- function(score, N, eta0_multiplier = 1) {
  if (length(eta0_multiplier) != 1L || !is.finite(eta0_multiplier) ||
      eta0_multiplier <= 0) {
    stop("eta0_multiplier must be one positive finite number")
  }
  score <- as.data.frame(score)
  raw_eta0 <- score$eta0
  raw_objective <- score$objective
  score$raw_eta0 <- raw_eta0
  score$eta0_used <- raw_eta0
  score$eta0_multiplier <- eta0_multiplier
  score$profiled_objective <- raw_objective
  score$objective_penalty_from_eta0_multiplier <- 0

  if (eta0_multiplier == 1 || is.na(raw_eta0) || is.infinite(raw_eta0)) {
    return(score)
  }

  eta0_used <- eta0_multiplier * raw_eta0
  inactive_dimension <- score$inactive_dimension
  h_used <- 1 / N + score$A_gamma / eta0_used
  inactive_loglik <- -inactive_dimension / 2 * log(2 * pi) -
    score$logdet_S_gamma / 2 -
    inactive_dimension / 2 * log(h_used) -
    score$r_gamma / (2 * h_used)
  corrected_objective <- score$active_loglik - score$KL + inactive_loglik

  score$objective <- corrected_objective
  score$eta0 <- eta0_used
  score$eta0_used <- eta0_used
  score$h_hat <- h_used
  score$inactive_loglik <- inactive_loglik
  score$objective_penalty_from_eta0_multiplier <-
    corrected_objective - raw_objective
  score
}

fit_aiwn_path <- function(v, Rbar, N, L_grid = 1:5,
                          scaled_prior_variance = 0.2,
                          selection_tolerance = 1e-3,
                          cs_coverage = 0.95,
                          cs_min_abs_corr = 0.5,
                          susie_args = list(),
                          chol_Rbar = NULL,
                          logdet_Rbar = NULL,
                          eta0_multiplier = 1) {
  if (!requireNamespace("susieR", quietly = TRUE)) {
    stop("The susieR package is required")
  }

  v <- as.numeric(v)
  Rbar <- as.matrix(Rbar)
  J <- length(v)
  L_grid <- as.integer(L_grid)

  if (any(L_grid < 1L) || any(L_grid >= J) || anyDuplicated(L_grid)) {
    stop("L_grid must contain distinct integers between 1 and J - 1")
  }

  if (is.null(chol_Rbar)) chol_Rbar <- chol(Rbar)
  if (is.null(logdet_Rbar)) {
    logdet_Rbar <- aiwn_logdet_chol(chol_Rbar)
  }
  full_quadratic <- drop(crossprod(
    v, aiwn_solve_chol(chol_Rbar, v)
  ))

  fits <- vector("list", length(L_grid))
  names(fits) <- as.character(L_grid)
  scores <- vector("list", length(L_grid))

  base_args <- list(
    XtX = N * Rbar,
    Xty = N * v,
    yty = N,
    n = N,
    residual_variance = 1,
    estimate_residual_variance = FALSE,
    estimate_prior_variance = TRUE,
    scaled_prior_variance = scaled_prior_variance,
    standardize = FALSE,
    check_prior = FALSE,
    max_iter = 200,
    tol = 1e-3
  )
  fit_args <- utils::modifyList(base_args, susie_args)

  for (i in seq_along(L_grid)) {
    fit_args$L <- L_grid[i]
    fit <- do.call(aiwn_susie_ss_function(), fit_args)
    if (!isTRUE(fit$converged)) {
      warning(sprintf("SuSiE-RSS did not converge for L = %d", L_grid[i]))
    }

    fits[[i]] <- fit
    scores[[i]] <- tryCatch(
      aiwn_score_fit(
        fit = fit,
        v = v,
        Rbar = Rbar,
        N = N,
        chol_Rbar = chol_Rbar,
        logdet_Rbar = logdet_Rbar,
        full_quadratic = full_quadratic
      ),
      error = function(error) {
        if (!grepl(
          "failed to find a distinct configuration",
          conditionMessage(error), fixed = TRUE
        )) {
          stop(error)
        }
        warning(sprintf(
          paste0(
            "AIW-N could not find a positive-probability distinct ",
            "configuration for L = %d; excluding this L"
          ),
          L_grid[i]
        ))
        data.frame(
          L = L_grid[i],
          inactive_dimension = J - L_grid[i],
          objective = -Inf,
          eta0 = NA_real_,
          h_hat = NA_real_,
          A_gamma = NA_real_,
          r_gamma = NA_real_,
          residual_per_snp = NA_real_,
          logdet_S_gamma = NA_real_,
          active_loglik = NA_real_,
          inactive_loglik = NA_real_,
          KL = NA_real_,
          gamma_hat = NA_character_,
          distinct_log_probability = -Inf,
          stringsAsFactors = FALSE
        )
      }
    )
    scores[[i]] <- aiwn_rescore_with_eta0_multiplier(
      scores[[i]], N = N, eta0_multiplier = eta0_multiplier
    )
  }

  score_table <- do.call(rbind, scores)
  rownames(score_table) <- NULL
  best_objective <- max(score_table$objective)
  if (!is.finite(best_objective)) {
    stop("AIW-N could not score any fitted value of L")
  }
  eligible <- which(
    score_table$objective >= best_objective - selection_tolerance
  )
  best_index <- eligible[which.min(score_table$L[eligible])]
  selected_fit <- fits[[best_index]]
  selected_cs <- susieR::susie_get_cs(
    selected_fit,
    Xcorr = Rbar,
    coverage = cs_coverage,
    min_abs_corr = cs_min_abs_corr
  )
  purity_filtered_L <- length(selected_cs$cs)

  structure(
    list(
      scores = score_table,
      fits = fits,
      selected_L = score_table$L[best_index],
      selected_eta0 = score_table$eta0[best_index],
      selected_raw_eta0 = score_table$raw_eta0[best_index],
      selected_eta0_multiplier = eta0_multiplier,
      selected_fit = selected_fit,
      # Keep the literal pure-CS count, but do not turn a selected SER into L=0.
      purity_filtered_L = purity_filtered_L,
      reported_L = max(1L, purity_filtered_L),
      credible_sets = selected_cs,
      selection_tolerance = selection_tolerance,
      cs_coverage = cs_coverage,
      cs_min_abs_corr = cs_min_abs_corr,
      selected_gamma = as.integer(strsplit(
        score_table$gamma_hat[best_index], ",", fixed = TRUE
      )[[1L]]),
      eta0_multiplier = eta0_multiplier
    ),
    class = "aiwn_path"
  )
}

aiwn_score_fit_at_eta0 <- function(fit, v, Rbar, N, eta0,
                                    chol_Rbar = NULL,
                                    logdet_Rbar = NULL,
                                    full_quadratic = NULL) {
  if (length(eta0) != 1L || is.na(eta0) || eta0 <= 0) {
    stop("eta0 must be one positive value or Inf")
  }

  profiled_score <- aiwn_score_fit(
    fit = fit,
    v = v,
    Rbar = Rbar,
    N = N,
    chol_Rbar = chol_Rbar,
    logdet_Rbar = logdet_Rbar,
    full_quadratic = full_quadratic
  )
  J <- length(v)
  L <- profiled_score$L
  inactive_dimension <- J - L
  h_eta0 <- if (is.infinite(eta0)) {
    1 / N
  } else {
    1 / N + profiled_score$A_gamma / eta0
  }
  inactive_loglik <- -inactive_dimension / 2 * log(2 * pi) -
    profiled_score$logdet_S_gamma / 2 -
    inactive_dimension / 2 * log(h_eta0) -
    profiled_score$r_gamma / (2 * h_eta0)

  profiled_score$profiled_objective <- profiled_score$objective
  profiled_score$profiled_eta0 <- profiled_score$eta0
  profiled_score$profiled_h_hat <- profiled_score$h_hat
  profiled_score$objective <- profiled_score$active_loglik -
    profiled_score$KL +
    inactive_loglik
  profiled_score$eta0 <- eta0
  profiled_score$h_hat <- h_eta0
  profiled_score$inactive_loglik <- inactive_loglik
  profiled_score
}

fit_aiwn_fixed_eta0_path <- function(
    v, Rbar, N, eta0, L_grid = 1:5,
    scaled_prior_variance = 0.2,
    selection_tolerance = 1e-3,
    cs_coverage = 0.95,
    cs_min_abs_corr = 0.5,
    susie_args = list(),
    aiwn_path = NULL) {
  v <- as.numeric(v)
  Rbar <- as.matrix(Rbar)
  L_grid <- as.integer(L_grid)

  if (is.null(aiwn_path)) {
    aiwn_path <- fit_aiwn_path(
      v = v,
      Rbar = Rbar,
      N = N,
      L_grid = L_grid,
      scaled_prior_variance = scaled_prior_variance,
      selection_tolerance = selection_tolerance,
      cs_coverage = cs_coverage,
      cs_min_abs_corr = cs_min_abs_corr,
      susie_args = susie_args
    )
  }
  if (!identical(as.integer(names(aiwn_path$fits)), L_grid)) {
    stop("aiwn_path fits do not match L_grid")
  }

  chol_Rbar <- chol(Rbar)
  logdet_Rbar <- aiwn_logdet_chol(chol_Rbar)
  full_quadratic <- drop(crossprod(
    v, aiwn_solve_chol(chol_Rbar, v)
  ))
  scores <- vector("list", length(L_grid))

  for (i in seq_along(L_grid)) {
    scores[[i]] <- aiwn_score_fit_at_eta0(
      fit = aiwn_path$fits[[i]],
      v = v,
      Rbar = Rbar,
      N = N,
      eta0 = eta0,
      chol_Rbar = chol_Rbar,
      logdet_Rbar = logdet_Rbar,
      full_quadratic = full_quadratic
    )
  }

  score_table <- do.call(rbind, scores)
  rownames(score_table) <- NULL
  best_objective <- max(score_table$objective)
  eligible <- which(
    score_table$objective >= best_objective - selection_tolerance
  )
  best_index <- eligible[which.min(score_table$L[eligible])]
  selected_fit <- aiwn_path$fits[[best_index]]
  selected_cs <- susieR::susie_get_cs(
    selected_fit,
    Xcorr = Rbar,
    coverage = cs_coverage,
    min_abs_corr = cs_min_abs_corr
  )
  purity_filtered_L <- length(selected_cs$cs)

  structure(
    list(
      scores = score_table,
      fits = aiwn_path$fits,
      selected_L = score_table$L[best_index],
      selected_eta0 = eta0,
      selected_fit = selected_fit,
      purity_filtered_L = purity_filtered_L,
      reported_L = max(1L, purity_filtered_L),
      credible_sets = selected_cs,
      selection_tolerance = selection_tolerance,
      cs_coverage = cs_coverage,
      cs_min_abs_corr = cs_min_abs_corr,
      selected_gamma = as.integer(strsplit(
        score_table$gamma_hat[best_index], ",", fixed = TRUE
      )[[1L]]),
      profiled_path = aiwn_path
    ),
    class = c("aiwn_fixed_eta0_path", "aiwn_path")
  )
}

fit_aiwn_lambda_grid <- function(
    v, R0, N,
    lambda_grid = c(1e-5, 1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5),
    L_grid = 1:5,
    scaled_prior_variance = 0.2,
    selection_tolerance = 1e-3,
    cs_coverage = 0.95,
    cs_min_abs_corr = 0.5,
    susie_args = list(),
    reference_models = NULL,
    eta0_multiplier = 1) {
  v <- as.numeric(v)
  R0 <- as.matrix(R0)
  lambda_grid <- as.numeric(lambda_grid)
  L_grid <- as.integer(L_grid)
  J <- length(v)

  if (!all(dim(R0) == c(J, J))) {
    stop("R0 must be a J by J matrix")
  }
  if (any(!is.finite(lambda_grid)) ||
      any(lambda_grid < 0 | lambda_grid > 1) ||
      anyDuplicated(lambda_grid)) {
    stop("lambda_grid must contain distinct finite values between 0 and 1")
  }
  if (!is.null(reference_models)) {
    if (length(reference_models) != length(lambda_grid) ||
        !isTRUE(all.equal(
          vapply(reference_models, `[[`, numeric(1), "lambda"),
          lambda_grid
        ))) {
      stop("reference_models must match lambda_grid in the same order")
    }
  }

  paths <- vector("list", length(lambda_grid))
  names(paths) <- format(lambda_grid, scientific = FALSE, trim = TRUE)
  score_tables <- vector("list", length(lambda_grid))

  for (i in seq_along(lambda_grid)) {
    lambda <- lambda_grid[i]
    if (is.null(reference_models)) {
      Rbar <- (1 - lambda) * R0 + lambda * diag(1, J)
      Rbar <- stats::cov2cor((Rbar + t(Rbar)) / 2)
      chol_Rbar <- NULL
      logdet_Rbar <- NULL
    } else {
      Rbar <- reference_models[[i]]$Rbar
      chol_Rbar <- reference_models[[i]]$chol_Rbar
      logdet_Rbar <- reference_models[[i]]$logdetRbar
    }
    path <- fit_aiwn_path(
      v = v,
      Rbar = Rbar,
      N = N,
      L_grid = L_grid,
      scaled_prior_variance = scaled_prior_variance,
      selection_tolerance = selection_tolerance,
      cs_coverage = cs_coverage,
      cs_min_abs_corr = cs_min_abs_corr,
      susie_args = susie_args,
      chol_Rbar = chol_Rbar,
      logdet_Rbar = logdet_Rbar,
      eta0_multiplier = eta0_multiplier
    )
    paths[[i]] <- path
    score_tables[[i]] <- transform(path$scores, lambda = lambda)
  }

  scores <- do.call(rbind, score_tables)
  rownames(scores) <- NULL
  best_objective <- max(scores$objective)
  eligible <- which(scores$objective >= best_objective - selection_tolerance)
  best_index <- eligible[
    order(scores$L[eligible], scores$lambda[eligible])[1L]
  ]
  selected_lambda <- scores$lambda[best_index]
  selected_L <- scores$L[best_index]
  lambda_index <- match(selected_lambda, lambda_grid)
  L_index <- match(selected_L, L_grid)
  selected_fit <- paths[[lambda_index]]$fits[[L_index]]
  if (is.null(reference_models)) {
    selected_Rbar <- (1 - selected_lambda) * R0 +
      selected_lambda * diag(1, J)
    selected_Rbar <- stats::cov2cor(
      (selected_Rbar + t(selected_Rbar)) / 2
    )
  } else {
    selected_Rbar <- reference_models[[lambda_index]]$Rbar
  }
  selected_cs <- susieR::susie_get_cs(
    selected_fit,
    Xcorr = selected_Rbar,
    coverage = cs_coverage,
    min_abs_corr = cs_min_abs_corr
  )
  purity_filtered_L <- length(selected_cs$cs)

  lambda_profile <- do.call(rbind, lapply(seq_along(lambda_grid), function(i) {
    rows <- scores$lambda == lambda_grid[i]
    lambda_scores <- scores[rows, , drop = FALSE]
    lambda_best <- max(lambda_scores$objective)
    lambda_eligible <- which(
      lambda_scores$objective >= lambda_best - selection_tolerance
    )
    chosen <- lambda_eligible[
      which.min(lambda_scores$L[lambda_eligible])
    ]
    data.frame(
      lambda = lambda_grid[i],
      objective = lambda_scores$objective[chosen],
      selected_L = lambda_scores$L[chosen],
      eta0 = lambda_scores$eta0[chosen],
      raw_eta0 = lambda_scores$raw_eta0[chosen],
      eta0_multiplier = eta0_multiplier,
      objective_penalty_from_eta0_multiplier =
        lambda_scores$objective_penalty_from_eta0_multiplier[chosen],
      stringsAsFactors = FALSE
    )
  }))

  structure(
    list(
      scores = scores,
      lambda_profile = lambda_profile,
      paths = paths,
      selected_lambda = selected_lambda,
      selected_L = selected_L,
      selected_eta0 = scores$eta0[best_index],
      selected_raw_eta0 = scores$raw_eta0[best_index],
      selected_eta0_multiplier = eta0_multiplier,
      selected_objective = scores$objective[best_index],
      selected_fit = selected_fit,
      selected_Rbar = selected_Rbar,
      purity_filtered_L = purity_filtered_L,
      reported_L = max(1L, purity_filtered_L),
      credible_sets = selected_cs,
      selected_gamma = as.integer(strsplit(
        scores$gamma_hat[best_index], ",", fixed = TRUE
      )[[1L]]),
      lambda_grid = lambda_grid,
      L_grid = L_grid,
      selection_tolerance = selection_tolerance,
      cs_coverage = cs_coverage,
      cs_min_abs_corr = cs_min_abs_corr,
      eta0_multiplier = eta0_multiplier
    ),
    class = "aiwn_lambda_grid"
  )
}

# -------------------------------------------------------------------------
# Student-t AIW score
# -------------------------------------------------------------------------

aiw_lgamma_increment <- function(base, increment) {
  if (!is.finite(base) || base <= 0 ||
      !is.finite(increment) || increment <= 0) {
    stop("base and increment must be finite and positive")
  }
  lgamma(increment) - lbeta(base, increment)
}

aiw_inactive_loglik <- function(eta0, A_gamma, r_gamma, logdet_S_gamma,
                                J, L, N) {
  if (!is.finite(eta0) || eta0 <= 0) {
    stop("eta0 must be finite and positive")
  }
  inactive_dimension <- J - L
  if (inactive_dimension < 1L) stop("AIW requires L < J")
  df <- eta0 + L + 2
  total_df <- eta0 + J + 2
  active_scale <- A_gamma + eta0 / N
  aiw_lgamma_increment(df / 2, inactive_dimension / 2) -
    inactive_dimension / 2 * log(pi) -
    logdet_S_gamma / 2 -
    inactive_dimension / 2 * log(active_scale) -
    total_df / 2 * log1p(r_gamma / active_scale)
}

aiw_inactive_normal_limit <- function(r_gamma, logdet_S_gamma, J, L, N) {
  inactive_dimension <- J - L
  -inactive_dimension / 2 * log(2 * pi) -
    logdet_S_gamma / 2 +
    inactive_dimension / 2 * log(N) -
    N * r_gamma / 2
}

aiw_profile_eta0 <- function(A_gamma, r_gamma, logdet_S_gamma, J, L, N,
                             eta0_bounds = c(1e-4, 1e12),
                             grid_size = 161L,
                             boundary_tolerance = 1e-6) {
  eta0_bounds <- as.numeric(eta0_bounds)
  if (length(eta0_bounds) != 2L || any(!is.finite(eta0_bounds)) ||
      eta0_bounds[1L] <= 0 || eta0_bounds[2L] <= eta0_bounds[1L]) {
    stop("eta0_bounds must contain two increasing positive values")
  }
  grid_size <- as.integer(grid_size)
  if (grid_size < 5L) stop("grid_size must be at least 5")

  log_grid <- seq(log(eta0_bounds[1L]), log(eta0_bounds[2L]),
                  length.out = grid_size)
  objective_on_log_scale <- function(log_eta0) {
    aiw_inactive_loglik(
      eta0 = exp(log_eta0), A_gamma = A_gamma, r_gamma = r_gamma,
      logdet_S_gamma = logdet_S_gamma, J = J, L = L, N = N
    )
  }
  grid_objective <- vapply(log_grid, objective_on_log_scale, numeric(1))
  grid_best <- which.max(grid_objective)
  bracket <- c(max(1L, grid_best - 1L), min(grid_size, grid_best + 1L))
  refined <- stats::optimize(
    objective_on_log_scale,
    interval = log_grid[bracket], maximum = TRUE, tol = 1e-10
  )
  finite_log_eta0 <- c(log_grid[1L], refined$maximum, log_grid[grid_size])
  finite_objective <- c(
    grid_objective[1L], refined$objective, grid_objective[grid_size]
  )
  finite_best <- which.max(finite_objective)
  normal_limit <- aiw_inactive_normal_limit(
    r_gamma, logdet_S_gamma, J, L, N
  )
  if (normal_limit >= finite_objective[finite_best] - boundary_tolerance) {
    return(list(
      eta0 = Inf, inactive_loglik = normal_limit, boundary = "infinity",
      finite_eta0 = exp(finite_log_eta0[finite_best]),
      finite_inactive_loglik = finite_objective[finite_best],
      normal_limit = normal_limit
    ))
  }
  eta0 <- exp(finite_log_eta0[finite_best])
  boundary <- c("lower", "interior", "upper")[finite_best]
  list(
    eta0 = eta0,
    inactive_loglik = finite_objective[finite_best],
    boundary = boundary,
    finite_eta0 = eta0,
    finite_inactive_loglik = finite_objective[finite_best],
    normal_limit = normal_limit
  )
}

aiw_score_fit <- function(fit, v, Rbar, N, chol_Rbar = NULL,
                          logdet_Rbar = NULL, full_quadratic = NULL,
                          eta0_bounds = c(1e-4, 1e12),
                          eta0_grid_size = 161L,
                          eta0_boundary_tolerance = 1e-6) {
  normal_score <- aiwn_score_fit(
    fit, v, Rbar, N, chol_Rbar, logdet_Rbar, full_quadratic
  )
  eta_profile <- aiw_profile_eta0(
    A_gamma = normal_score$A_gamma,
    r_gamma = normal_score$r_gamma,
    logdet_S_gamma = normal_score$logdet_S_gamma,
    J = length(v), L = normal_score$L, N = N,
    eta0_bounds = eta0_bounds,
    grid_size = eta0_grid_size,
    boundary_tolerance = eta0_boundary_tolerance
  )
  data.frame(
    L = normal_score$L,
    objective = normal_score$active_loglik - normal_score$KL +
      eta_profile$inactive_loglik,
    eta0 = eta_profile$eta0,
    eta0_boundary = eta_profile$boundary,
    finite_eta0 = eta_profile$finite_eta0,
    A_gamma = normal_score$A_gamma,
    r_gamma = normal_score$r_gamma,
    residual_per_snp = normal_score$residual_per_snp,
    logdet_S_gamma = normal_score$logdet_S_gamma,
    active_loglik = normal_score$active_loglik,
    inactive_loglik = eta_profile$inactive_loglik,
    normal_limit = eta_profile$normal_limit,
    KL = normal_score$KL,
    gamma_hat = normal_score$gamma_hat,
    distinct_log_probability = normal_score$distinct_log_probability,
    stringsAsFactors = FALSE
  )
}

fit_aiw_path <- function(
    v, Rbar, N, L_grid = 1:5,
    scaled_prior_variance = 0.2,
    selection_tolerance = 1e-3,
    cs_coverage = 0.95,
    cs_min_abs_corr = 0.5,
    eta0_bounds = c(1e-4, 1e12),
    eta0_grid_size = 161L,
    eta0_boundary_tolerance = 1e-6,
    susie_args = list(), aiwn_path = NULL) {
  if (is.null(aiwn_path)) {
    aiwn_path <- fit_aiwn_path(
      v, Rbar, N, L_grid, scaled_prior_variance,
      selection_tolerance, cs_coverage, cs_min_abs_corr,
      susie_args = susie_args, eta0_multiplier = 1
    )
  }
  if (!identical(as.integer(names(aiwn_path$fits)), as.integer(L_grid))) {
    stop("aiwn_path fits do not match L_grid")
  }
  chol_Rbar <- chol(Rbar)
  logdet_Rbar <- aiwn_logdet_chol(chol_Rbar)
  full_quadratic <- drop(crossprod(
    v, aiwn_solve_chol(chol_Rbar, v)
  ))
  scores <- lapply(aiwn_path$fits, function(fit) {
    aiw_score_fit(
      fit, v, Rbar, N, chol_Rbar, logdet_Rbar, full_quadratic,
      eta0_bounds, eta0_grid_size, eta0_boundary_tolerance
    )
  })
  score_table <- do.call(rbind, scores)
  rownames(score_table) <- NULL
  best <- max(score_table$objective)
  eligible <- which(score_table$objective >= best - selection_tolerance)
  best_index <- eligible[which.min(score_table$L[eligible])]
  selected_fit <- aiwn_path$fits[[best_index]]
  selected_cs <- susieR::susie_get_cs(
    selected_fit, Xcorr = Rbar, coverage = cs_coverage,
    min_abs_corr = cs_min_abs_corr
  )
  structure(list(
    scores = score_table,
    fits = aiwn_path$fits,
    selected_L = score_table$L[best_index],
    selected_eta0 = score_table$eta0[best_index],
    selected_fit = selected_fit,
    credible_sets = selected_cs,
    selected_gamma = as.integer(strsplit(
      score_table$gamma_hat[best_index], ",", fixed = TRUE
    )[[1L]]),
    aiwn_path = aiwn_path
  ), class = "aiw_path")
}

# -------------------------------------------------------------------------
# Public collaborator-facing interfaces
# -------------------------------------------------------------------------

aiw_default_lambda_grid <- function() {
  c(1e-5, 1e-4, 1e-3, 2e-3, 4e-3, 6e-3, 8e-3, 1e-2)
}

aiw_validate_public_inputs <- function(v, R0, N, L, lambda_grid,
                                       eta0_multiplier) {
  v <- as.numeric(v)
  R0 <- as.matrix(R0)
  J <- length(v)
  if (J < 2L || any(!is.finite(v)) ||
      !all(dim(R0) == c(J, J)) || any(!is.finite(R0))) {
    stop("v and R0 must be a finite vector and matching square matrix")
  }
  if (length(N) != 1L || !is.finite(N) || N <= 0) {
    stop("N must be one positive finite number")
  }
  L <- as.integer(L)
  if (length(L) != 1L || is.na(L) || L < 1L || L >= J) {
    stop("L must be one integer between 1 and length(v) - 1")
  }
  lambda_grid <- as.numeric(lambda_grid)
  if (!length(lambda_grid) || any(!is.finite(lambda_grid)) ||
      any(lambda_grid < 0 | lambda_grid >= 1) ||
      anyDuplicated(lambda_grid)) {
    stop("lambda_grid must contain distinct values in [0, 1)")
  }
  if (length(eta0_multiplier) != 1L || !is.finite(eta0_multiplier) ||
      eta0_multiplier <= 0) {
    stop("eta0_multiplier must be one positive finite number")
  }
  list(
    v = v, R0 = (R0 + t(R0)) / 2, N = as.numeric(N), L = L,
    lambda_grid = lambda_grid, eta0_multiplier = eta0_multiplier
  )
}

# Fit AIW-N while comparing candidate effect counts 1, ..., L and all lambda
# values. eta0 is profiled in closed form for every (lambda, effect-count)
# candidate, multiplied by eta0_multiplier, and then used to recompute that
# candidate's objective before the final comparison.
fit_aiwn <- function(
    v, R0, N, L = 5L,
    lambda_grid = aiw_default_lambda_grid(),
    eta0_multiplier = 1.424,
    prior_variance = 0.2,
    estimate_prior_variance = TRUE,
    prior_weights = NULL,
    coverage = 0.95,
    min_abs_corr = 0.5,
    selection_tolerance = 1e-3,
    max_iter = 500L,
    tol = 1e-3,
    verbose = FALSE,
    track_fit = FALSE) {
  checked <- aiw_validate_public_inputs(
    v, R0, N, L, lambda_grid, eta0_multiplier
  )
  if (verbose) message("AIW-N: fitting lambda-by-L path")
  result <- suppressMessages(fit_aiwn_lambda_grid(
    v = checked$v, R0 = checked$R0, N = checked$N,
    lambda_grid = checked$lambda_grid,
    L_grid = seq_len(checked$L),
    scaled_prior_variance = prior_variance,
    selection_tolerance = selection_tolerance,
    cs_coverage = coverage,
    cs_min_abs_corr = min_abs_corr,
    susie_args = list(
      max_iter = max_iter, tol = tol,
      estimate_prior_variance = estimate_prior_variance,
      prior_weights = prior_weights,
      verbose = verbose, track_fit = track_fit
    ),
    eta0_multiplier = checked$eta0_multiplier
  ))
  result$method <- "AIW-N"
  result$pip <- as.numeric(result$selected_fit$pip)
  result$eta0_used <- result$selected_eta0
  class(result) <- c("aiwn_fit", class(result))
  if (verbose) print(result)
  result
}

# Fit the Student-t AIW working model. As in AIW-N, L is the maximum candidate
# effect count and eta0 is always learned separately for every candidate.
fit_aiw <- function(
    v, R0, N, L = 5L,
    lambda_grid = aiw_default_lambda_grid(),
    eta0_multiplier = 1,
    prior_variance = 0.2,
    estimate_prior_variance = TRUE,
    prior_weights = NULL,
    coverage = 0.95,
    min_abs_corr = 0.5,
    selection_tolerance = 1e-3,
    eta0_bounds = c(1e-4, 1e12),
    eta0_grid_size = 161L,
    max_iter = 500L,
    tol = 1e-3,
    verbose = FALSE,
    track_fit = FALSE) {
  checked <- aiw_validate_public_inputs(
    v, R0, N, L, lambda_grid, eta0_multiplier
  )
  paths <- vector("list", length(checked$lambda_grid))
  score_tables <- vector("list", length(checked$lambda_grid))
  J <- length(checked$v)
  for (i in seq_along(checked$lambda_grid)) {
    lambda <- checked$lambda_grid[i]
    if (verbose) message(sprintf("AIW: lambda=%g", lambda))
    Rbar <- (1 - lambda) * checked$R0 + lambda * diag(J)
    Rbar <- stats::cov2cor((Rbar + t(Rbar)) / 2)
    path <- fit_aiw_path(
      v = checked$v, Rbar = Rbar, N = checked$N,
      L_grid = seq_len(checked$L),
      scaled_prior_variance = prior_variance,
      selection_tolerance = selection_tolerance,
      cs_coverage = coverage, cs_min_abs_corr = min_abs_corr,
      eta0_bounds = eta0_bounds,
      eta0_grid_size = eta0_grid_size,
      susie_args = list(
        max_iter = max_iter, tol = tol,
        estimate_prior_variance = estimate_prior_variance,
        prior_weights = prior_weights,
        verbose = verbose, track_fit = track_fit
      )
    )
    scores <- path$scores
    scores$J <- J
    scores$lambda <- lambda
    scores$raw_eta0 <- scores$eta0
    scores$eta0_used <- ifelse(
      is.infinite(scores$eta0), Inf,
      checked$eta0_multiplier * scores$eta0
    )
    for (row in seq_len(nrow(scores))) {
      inactive <- if (is.infinite(scores$eta0_used[row])) {
        scores$normal_limit[row]
      } else {
        aiw_inactive_loglik(
          scores$eta0_used[row], scores$A_gamma[row],
          scores$r_gamma[row], scores$logdet_S_gamma[row],
          J, scores$L[row], checked$N
        )
      }
      scores$inactive_loglik[row] <- inactive
      scores$objective[row] <- scores$active_loglik[row] -
        scores$KL[row] + inactive
    }
    scores$eta0_multiplier <- checked$eta0_multiplier
    paths[[i]] <- path
    paths[[i]]$Rbar <- Rbar
    score_tables[[i]] <- scores
  }
  scores <- do.call(rbind, score_tables)
  rownames(scores) <- NULL
  best <- max(scores$objective)
  eligible <- which(scores$objective >= best - selection_tolerance)
  best_index <- eligible[
    order(scores$L[eligible], scores$lambda[eligible])[1L]
  ]
  lambda_index <- match(scores$lambda[best_index], checked$lambda_grid)
  L_index <- match(scores$L[best_index], seq_len(checked$L))
  selected_fit <- paths[[lambda_index]]$fits[[L_index]]
  selected_Rbar <- paths[[lambda_index]]$Rbar
  credible_sets <- susieR::susie_get_cs(
    selected_fit, Xcorr = selected_Rbar,
    coverage = coverage, min_abs_corr = min_abs_corr
  )
  result <- list(
    method = "AIW",
    selected_fit = selected_fit,
    pip = as.numeric(selected_fit$pip),
    credible_sets = credible_sets,
    selected_Rbar = selected_Rbar,
    selected_lambda = scores$lambda[best_index],
    selected_L = scores$L[best_index],
    selected_raw_eta0 = scores$raw_eta0[best_index],
    selected_eta0 = scores$eta0_used[best_index],
    eta0_used = scores$eta0_used[best_index],
    eta0_multiplier = checked$eta0_multiplier,
    selected_objective = scores$objective[best_index],
    scores = scores,
    paths = paths,
    call = match.call()
  )
  class(result) <- "aiw_fit"
  if (verbose) print(result)
  result
}

print.aiwn_fit <- function(x, ...) {
  cat("AIW-N fit\n")
  cat("  lambda:", format(x$selected_lambda),
      " selected L:", x$selected_L, "\n")
  cat("  raw eta0:", format(x$selected_raw_eta0),
      " eta0 used:", format(x$selected_eta0),
      " multiplier:", format(x$selected_eta0_multiplier), "\n")
  invisible(x)
}

print.aiw_fit <- function(x, ...) {
  cat("AIW fit\n")
  cat("  lambda:", format(x$selected_lambda),
      " selected L:", x$selected_L, "\n")
  cat("  raw eta0:", format(x$selected_raw_eta0),
      " eta0 used:", format(x$selected_eta0),
      " multiplier:", format(x$eta0_multiplier), "\n")
  invisible(x)
}

# Paper-facing interface to AIW and AIW-N. The primary inputs are the marginal
# association v = X'y/N, reference LD R0, and GWAS sample size N; z-scores are
# accepted as an explicitly named alternative. approximation="N" uses the
# normal approximation and is the default; approximation="t" uses the
# Student-t working likelihood. Both always estimate eta0 for every
# (lambda, L) candidate before applying the requested multiplier and comparing
# the resulting objectives.
susie_aiw <- function(
    v = NULL, R0, N, z = NULL,
    L = min(10L, ncol(R0)),
    approximation = c("N", "t"),
    lambda_grid = aiw_default_lambda_grid(),
    eta0_multiplier = NULL,
    scaled_prior_variance = 0.2,
    estimate_prior_variance = TRUE,
    prior_weights = NULL,
    coverage = 0.95,
    min_abs_corr = 0.5,
    ser_fallback = TRUE,
    partition_ser = FALSE,
    max_partition_sets = 3L,
    min_partition_probability = 0.1,
    selection_tolerance = 1e-3,
    eta0_bounds = c(1e-4, 1e12),
    eta0_grid_size = 161L,
    max_iter = 500L,
    tol = 1e-3,
    verbose = FALSE,
    track_fit = FALSE) {
  if (!exists("iw_prepare_summary_inputs", mode = "function") ||
      !exists("iw_finalize_public_fit", mode = "function")) {
    stop("Source 03-SER-fallback-and-plot.R before calling susie_aiw()")
  }
  supplied_call <- match.call()
  if (length(approximation) > 1L) approximation <- approximation[1L]
  if (tolower(approximation) == "n") approximation <- "N"
  if (tolower(approximation) == "t") approximation <- "t"
  approximation <- match.arg(approximation, c("N", "t"))
  if (is.null(eta0_multiplier)) {
    eta0_multiplier <- if (approximation == "N") 1.424 else 1
  }
  summary_data <- iw_prepare_summary_inputs(
    v = v, R0 = R0, N = N, z = z
  )

  if (approximation == "N") {
    method_fit <- fit_aiwn(
      v = summary_data$v, R0 = summary_data$R0,
      N = summary_data$effective_N, L = L,
      lambda_grid = lambda_grid,
      eta0_multiplier = eta0_multiplier,
      prior_variance = scaled_prior_variance,
      estimate_prior_variance = estimate_prior_variance,
      prior_weights = prior_weights,
      coverage = coverage, min_abs_corr = min_abs_corr,
      selection_tolerance = selection_tolerance,
      max_iter = max_iter, tol = tol, verbose = verbose,
      track_fit = track_fit
    )
    raw_eta0 <- method_fit$selected_raw_eta0
    eta0_used <- method_fit$selected_eta0
    profile <- method_fit$scores
    method_label <- "AIW-N"
  } else {
    method_fit <- fit_aiw(
      v = summary_data$v, R0 = summary_data$R0,
      N = summary_data$effective_N, L = L,
      lambda_grid = lambda_grid,
      eta0_multiplier = eta0_multiplier,
      prior_variance = scaled_prior_variance,
      estimate_prior_variance = estimate_prior_variance,
      prior_weights = prior_weights,
      coverage = coverage, min_abs_corr = min_abs_corr,
      selection_tolerance = selection_tolerance,
      eta0_bounds = eta0_bounds, eta0_grid_size = eta0_grid_size,
      max_iter = max_iter, tol = tol, verbose = verbose,
      track_fit = track_fit
    )
    raw_eta0 <- method_fit$selected_raw_eta0
    eta0_used <- method_fit$selected_eta0
    profile <- method_fit$scores
    method_label <- "AIW-t"
  }

  diagnostics <- list(
    approximation = approximation,
    selected_lambda = method_fit$selected_lambda,
    selected_Rbar = method_fit$selected_Rbar,
    raw_eta0 = raw_eta0,
    eta0 = eta0_used,
    eta0_multiplier = eta0_multiplier,
    selected_L = method_fit$selected_L,
    profile = profile,
    paths = if (isTRUE(track_fit)) method_fit$paths else NULL,
    input = summary_data$input,
    supplied_N = summary_data$supplied_N,
    effective_N = summary_data$effective_N
  )
  result <- iw_finalize_public_fit(
    method_fit = method_fit, v = summary_data$v,
    N = summary_data$effective_N,
    R = method_fit$selected_Rbar,
    method = method_label, subclass = "susie_aiw",
    coverage = coverage, min_abs_corr = min_abs_corr,
    ser_fallback = ser_fallback,
    partition_ser = partition_ser,
    max_partition_sets = max_partition_sets,
    min_partition_probability = min_partition_probability,
    prior_variance = scaled_prior_variance,
    max_iter = max_iter, tol = tol, call = supplied_call,
    diagnostics = diagnostics
  )
  result$approximation <- approximation
  result$selected_lambda <- method_fit$selected_lambda
  result$selected_Rbar <- method_fit$selected_Rbar
  result$selected_raw_eta0 <- raw_eta0
  result$selected_eta0 <- eta0_used
  result$eta0 <- eta0_used
  result$eta0_used <- eta0_used
  result$eta0_multiplier <- eta0_multiplier
  result$selected_L <- method_fit$selected_L
  result
}

print.susie_aiw <- function(x, ...) {
  cat(if (identical(x$approximation, "N")) "AIW-N fit\n" else "AIW-t fit\n")
  cat("  lambda:", format(x$selected_lambda),
      " selected L:", x$selected_L, "\n")
  cat("  raw eta0:", format(x$selected_raw_eta0),
      " eta0 used:", format(x$selected_eta0),
      " multiplier:", format(x$eta0_multiplier), "\n")
  cat("  reported CS:", length(susie_get_reported_cs(x)),
      " SER fallback:", isTRUE(x$postprocess$fallback_used), "\n")
  invisible(x)
}
