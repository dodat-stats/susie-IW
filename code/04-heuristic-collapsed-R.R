# Heuristic collapsed-R inference using a correlation-standardized surrogate.
#
# Source code/01-collapsed-R.R and code/03-SER-fallback-and-plot.R first.
# The exact collapsed-R model, q(omega), eta0 profile, prior-variance update,
# and reported ELBO are unchanged. Only the q(b) update is modified:
#
#   Rhat   = E_q(1 / omega) * (eta0 * Rbar + N * tcrossprod(v))
#   Rtilde = cov2cor(Rhat)
#
# and Rtilde is supplied to the IBSS updates. Because E_q(1 / omega) is a
# positive scalar, it cancels from Rtilde. The resulting algorithm is not
# coordinate ascent for the reported collapsed-R ELBO.

hiw_require_core <- function(require_postprocessing = FALSE) {
  required <- c(
    "scr_initialize_effect", "scr_update_effect",
    "scr_effect_second_moment", "scr_optimize_profiled_nu0",
    "scr_gig_log_integral", "scr_gig_inv_mean_stable",
    "scr_gig_chisq_kl_from_log_integral", "scr_lower_bound",
    "scr_logdet_spd", "collapsed_r_default_lambda_grid",
    "collapsed_r_validate_inputs"
  )
  missing <- required[!vapply(required, exists, logical(1), mode = "function")]
  if (length(missing)) {
    stop(
      "Source code/01-collapsed-R.R before code/04-heuristic-collapsed-R.R"
    )
  }
  if (require_postprocessing &&
      (!exists("iw_prepare_summary_inputs", mode = "function") ||
       !exists("iw_finalize_public_fit", mode = "function"))) {
    stop(
      "Source code/03-SER-fallback-and-plot.R before calling ",
      "susie_iw_heuristic()"
    )
  }
  invisible(TRUE)
}

hiw_make_surrogate <- function(A, inverse_mean) {
  A <- (as.matrix(A) + t(as.matrix(A))) / 2
  diagonal <- diag(A)
  if (any(!is.finite(diagonal)) || any(diagonal <= 0)) {
    stop("The collapsed surrogate must have a strictly positive diagonal")
  }
  Rhat <- inverse_mean * A
  Rtilde <- stats::cov2cor(Rhat)
  Rtilde <- (Rtilde + t(Rtilde)) / 2
  diag(Rtilde) <- 1
  list(Rhat = Rhat, Rtilde = Rtilde)
}

hiw_alpha_matrix <- function(effects) {
  do.call(rbind, lapply(effects, function(effect) {
    as.numeric(effect$alpha)
  }))
}

hiw_fit_core <- function(
    x, Rbar, N, L = 2L, sigma2 = 0.2,
    nu0 = 100, estimate_nu0 = TRUE,
    nu0_bounds = c(1, 1e8), nu0_coarse_grid_size = 15L,
    pi_prior = NULL, estimate_sigma2 = TRUE,
    max_iter = 200L, tol = 1e-4, stop_when_stable = TRUE,
    verbose = FALSE, pre = NULL, track_elbo_updates = FALSE) {
  hiw_require_core()
  x <- as.numeric(x)
  Rbar <- as.matrix(Rbar)
  J <- length(x)
  L <- as.integer(L)
  max_iter <- as.integer(max_iter)
  if (length(sigma2) == 1L) sigma2 <- rep(as.numeric(sigma2), L)
  stopifnot(
    nrow(Rbar) == J, ncol(Rbar) == J,
    N > 0, L >= 1L, length(sigma2) == L, all(sigma2 > 0),
    length(nu0) == 1L, is.finite(nu0), nu0 > 0,
    length(nu0_bounds) == 2L, nu0_bounds[1L] > 0,
    nu0_bounds[2L] > nu0_bounds[1L],
    max_iter >= 1L, tol >= 0
  )

  if (is.null(pi_prior)) pi_prior <- rep(1 / J, J)
  pi_prior <- as.numeric(pi_prior)
  pi_prior <- pi_prior / sum(pi_prior)
  if (length(pi_prior) != J || any(!is.finite(pi_prior)) ||
      any(pi_prior <= 0)) {
    stop("All prior weights must be finite and positive")
  }

  if (is.null(pre)) {
    chol_Rbar <- chol((Rbar + t(Rbar)) / 2)
    inverse_x <- backsolve(
      chol_Rbar, forwardsolve(t(chol_Rbar), x)
    )
    pre <- list(
      r0 = as.numeric(crossprod(x, inverse_x)),
      logdetRbar = 2 * sum(log(diag(chol_Rbar)))
    )
  }

  nu0 <- min(max(nu0, nu0_bounds[1L]), nu0_bounds[2L])
  initial_nu0 <- nu0
  A <- nu0 * Rbar + N * tcrossprod(x)
  gig_shape <- (nu0 + 3) / 2
  inverse_mean <- 1 / (nu0 + 1)
  chi <- 0
  omega_kl <- 0
  surrogate <- hiw_make_surrogate(A, inverse_mean)
  effects <- lapply(seq_len(L), function(ell) {
    scr_initialize_effect(J, sigma2[ell], pi_prior)
  })

  elbo <- numeric()
  nu0_trace <- numeric()
  nu0_boundary_trace <- character()
  objective_evaluations <- integer()
  state_trace <- data.frame(
    iteration = integer(),
    alpha_change = numeric(),
    eta0_log_change = numeric(),
    sigma2_log_change = numeric(),
    state_change = numeric(),
    elbo_change = numeric(),
    stringsAsFactors = FALSE
  )
  elbo_updates <- data.frame(
    update_index = integer(),
    iteration = integer(),
    update = character(),
    effect = integer(),
    nu0 = numeric(),
    lower_bound = numeric(),
    stringsAsFactors = FALSE
  )
  record_elbo_update <- function(iteration, update, effect = NA_integer_) {
    if (!isTRUE(track_elbo_updates)) return(invisible(NULL))
    elbo_updates <<- rbind(
      elbo_updates,
      data.frame(
        update_index = nrow(elbo_updates) + 1L,
        iteration = as.integer(iteration),
        update = update,
        effect = as.integer(effect),
        nu0 = nu0,
        lower_bound = scr_lower_bound(
          x, Rbar, N, nu0, sigma2, effects, A, inverse_mean,
          omega_kl, pi_prior, pre
        ),
        stringsAsFactors = FALSE
      )
    )
    invisible(NULL)
  }
  record_elbo_update(0L, "initial")

  converged <- FALSE
  for (iter in seq_len(max_iter)) {
    old_alpha <- hiw_alpha_matrix(effects)
    old_nu0 <- nu0
    old_sigma2 <- sigma2

    # This is the only heuristic step. Rtilde has unit diagonal, so passing
    # eta=1 makes scr_update_effect use Rtilde as its working LD matrix.
    m_sum <- Reduce(`+`, lapply(effects, `[[`, "m"))
    for (ell in seq_len(L)) {
      m_minus <- m_sum - effects[[ell]]$m
      old_m <- effects[[ell]]$m
      effects[[ell]] <- scr_update_effect(
        effect = effects[[ell]], m_minus = m_minus, x = x,
        A = surrogate$Rtilde, diagA = rep(1, J), N = N,
        sigma2 = sigma2[ell], eta = 1, pi_prior = pi_prior
      )
      m_sum <- m_sum - old_m + effects[[ell]]$m
      record_elbo_update(iter, "heuristic_effect", ell)
    }

    if (estimate_sigma2) {
      sigma2 <- vapply(effects, function(effect) {
        max(sum(as.numeric(effect$alpha) * effect$m2), 1e-12)
      }, numeric(1))
      record_elbo_update(iter, "prior_variance")
    }

    if (estimate_nu0) {
      profiled <- scr_optimize_profiled_nu0(
        current_nu0 = nu0, x = x, Rbar = Rbar, N = N,
        effects = effects, pre = pre, nu0_bounds = nu0_bounds,
        coarse_grid_size = nu0_coarse_grid_size
      )
      nu0 <- profiled$nu0
      gig_shape <- profiled$lambda
      chi <- profiled$chi
      inverse_mean <- scr_gig_inv_mean_stable(gig_shape, chi)
      omega_kl <- scr_gig_chisq_kl_from_log_integral(
        lambda = gig_shape, chi = chi, inverse_mean = inverse_mean,
        log_integral = profiled$log_integral
      )
      boundary <- profiled$boundary
      evaluations <- profiled$number_of_objective_evaluations
      update_name <- "eta0_omega"
    } else {
      gig_shape <- (nu0 + 3) / 2
      chi <- N * scr_effect_second_moment(effects, A)
      log_integral <- scr_gig_log_integral(gig_shape, chi)
      inverse_mean <- scr_gig_inv_mean_stable(gig_shape, chi)
      omega_kl <- scr_gig_chisq_kl_from_log_integral(
        lambda = gig_shape, chi = chi, inverse_mean = inverse_mean,
        log_integral = log_integral
      )
      boundary <- "fixed"
      evaluations <- 0L
      update_name <- "omega"
    }

    A <- nu0 * Rbar + N * tcrossprod(x)
    surrogate <- hiw_make_surrogate(A, inverse_mean)
    nu0_trace[iter] <- nu0
    nu0_boundary_trace[iter] <- boundary
    objective_evaluations[iter] <- evaluations
    elbo[iter] <- scr_lower_bound(
      x, Rbar, N, nu0, sigma2, effects, A, inverse_mean,
      omega_kl, pi_prior, pre
    )
    record_elbo_update(iter, update_name)

    new_alpha <- hiw_alpha_matrix(effects)
    alpha_change <- max(abs(new_alpha - old_alpha))
    eta0_log_change <- abs(log(nu0) - log(old_nu0))
    sigma2_log_change <- max(abs(log(sigma2) - log(old_sigma2)))
    state_change <- max(
      alpha_change, eta0_log_change, sigma2_log_change,
      na.rm = TRUE
    )
    elbo_change <- if (iter == 1L) NA_real_ else elbo[iter] - elbo[iter - 1L]
    state_trace <- rbind(
      state_trace,
      data.frame(
        iteration = iter, alpha_change = alpha_change,
        eta0_log_change = eta0_log_change,
        sigma2_log_change = sigma2_log_change,
        state_change = state_change, elbo_change = elbo_change
      )
    )

    if (verbose && (iter == 1L || iter %% 10L == 0L)) {
      message(sprintf(
        paste0(
          "heuristic iter=%d ELBO=%.6f eta0=%.5g ",
          "state_change=%.3g"
        ),
        iter, elbo[iter], nu0, state_change
      ))
    }
    if (isTRUE(stop_when_stable) && iter > 1L &&
        is.finite(state_change) && state_change < tol) {
      converged <- TRUE
      break
    }
  }

  alpha <- hiw_alpha_matrix(effects)
  rownames(alpha) <- paste0("effect", seq_len(L))
  colnames(alpha) <- paste0("j", seq_len(J))
  pip <- 1 - apply(1 - alpha, 2, prod)
  completed <- length(elbo)

  list(
    alpha = alpha,
    pip = pip,
    gamma_hat = apply(alpha, 1, which.max),
    nu0 = nu0,
    nu0_init = initial_nu0,
    nu0_bounds = nu0_bounds,
    nu0_boundary = tail(nu0_boundary_trace, 1L),
    nu0_trace = nu0_trace[seq_len(completed)],
    eta = inverse_mean,
    chi = chi,
    omega_kl = omega_kl,
    lower_bound = tail(elbo, 1L),
    elbo = elbo,
    elbo_updates = if (isTRUE(track_elbo_updates)) elbo_updates else NULL,
    state_trace = state_trace,
    effects = effects,
    sigma2 = sigma2,
    N = N,
    Rbar = Rbar,
    x = x,
    A = A,
    Rhat = surrogate$Rhat,
    Rtilde = surrogate$Rtilde,
    converged = converged,
    stop_when_stable = isTRUE(stop_when_stable),
    objective_evaluations = objective_evaluations[seq_len(completed)],
    heuristic = TRUE
  )
}

fit_collapsed_r_heuristic <- function(
    v, R0, N, L = 5L,
    lambda_grid = collapsed_r_default_lambda_grid(),
    eta0 = NULL, eta0_multiplier = 1.677,
    eta0_bounds = c(1, 1e8), eta0_init = 100,
    eta0_coarse_grid_size = 15L,
    prior_variance = 0.2, estimate_prior_variance = TRUE,
    prior_weights = NULL, max_iter = 200L, tol = 1e-4,
    stop_when_stable = TRUE, verbose = FALSE,
    track_elbo_updates = FALSE, keep_candidate_fits = FALSE) {
  hiw_require_core()
  checked <- collapsed_r_validate_inputs(v, R0, N, L, lambda_grid)
  v <- checked$v
  R0 <- checked$R0
  N <- checked$N
  L <- checked$L
  lambda_grid <- checked$lambda_grid
  J <- length(v)

  if (length(eta0_multiplier) != 1L || !is.finite(eta0_multiplier) ||
      eta0_multiplier <= 0) {
    stop("eta0_multiplier must be one positive finite number")
  }
  eta0_bounds <- as.numeric(eta0_bounds)
  if (length(eta0_bounds) != 2L || any(!is.finite(eta0_bounds)) ||
      eta0_bounds[1L] <= 0 || eta0_bounds[2L] <= eta0_bounds[1L]) {
    stop("eta0_bounds must contain two increasing positive values")
  }
  estimate_eta0 <- is.null(eta0)
  if (!estimate_eta0 &&
      (length(eta0) != 1L || !is.finite(eta0) || eta0 <= 0)) {
    stop("eta0 must be NULL or one positive finite number")
  }
  if (!is.null(prior_weights)) {
    prior_weights <- as.numeric(prior_weights)
    if (length(prior_weights) != J || any(!is.finite(prior_weights)) ||
        any(prior_weights <= 0)) {
      stop("prior_weights must contain one positive value per variant")
    }
    prior_weights <- prior_weights / sum(prior_weights)
  }

  candidate_fits <- vector("list", length(lambda_grid))
  profile_rows <- vector("list", length(lambda_grid))
  for (i in seq_along(lambda_grid)) {
    lambda <- lambda_grid[i]
    Rbar <- (1 - lambda) * R0 + lambda * diag(J)
    Rbar <- (Rbar + t(Rbar)) / 2
    chol_Rbar <- tryCatch(chol(Rbar), error = function(error) NULL)
    if (is.null(chol_Rbar)) {
      stop(sprintf(
        "Rbar is not positive definite at lambda=%g; increase lambda",
        lambda
      ))
    }
    inverse_v <- backsolve(
      chol_Rbar, forwardsolve(t(chol_Rbar), v)
    )
    pre <- list(
      r0 = as.numeric(crossprod(v, inverse_v)),
      logdetRbar = 2 * sum(log(diag(chol_Rbar)))
    )

    if (verbose) {
      message(sprintf(
        "heuristic collapsed-R: lambda=%g (%d/%d)",
        lambda, i, length(lambda_grid)
      ))
    }
    if (estimate_eta0) {
      raw_fit <- hiw_fit_core(
        x = v, Rbar = Rbar, N = N, L = L,
        sigma2 = prior_variance, nu0 = eta0_init,
        estimate_nu0 = TRUE, nu0_bounds = eta0_bounds,
        nu0_coarse_grid_size = eta0_coarse_grid_size,
        pi_prior = prior_weights,
        estimate_sigma2 = estimate_prior_variance,
        max_iter = max_iter, tol = tol,
        stop_when_stable = stop_when_stable, verbose = verbose,
        pre = pre, track_elbo_updates = track_elbo_updates
      )
      raw_eta0 <- raw_fit$nu0
      eta0_used <- min(eta0_bounds[2L], eta0_multiplier * raw_eta0)
      if (isTRUE(all.equal(eta0_used, raw_eta0))) {
        selected_fit <- raw_fit
      } else {
        selected_fit <- hiw_fit_core(
          x = v, Rbar = Rbar, N = N, L = L,
          sigma2 = prior_variance, nu0 = eta0_used,
          estimate_nu0 = FALSE, nu0_bounds = eta0_bounds,
          nu0_coarse_grid_size = eta0_coarse_grid_size,
          pi_prior = prior_weights,
          estimate_sigma2 = estimate_prior_variance,
          max_iter = max_iter, tol = tol,
          stop_when_stable = stop_when_stable, verbose = verbose,
          pre = pre, track_elbo_updates = track_elbo_updates
        )
      }
    } else {
      raw_fit <- NULL
      raw_eta0 <- as.numeric(eta0)
      eta0_used <- as.numeric(eta0)
      selected_fit <- hiw_fit_core(
        x = v, Rbar = Rbar, N = N, L = L,
        sigma2 = prior_variance, nu0 = eta0_used,
        estimate_nu0 = FALSE, nu0_bounds = eta0_bounds,
        nu0_coarse_grid_size = eta0_coarse_grid_size,
        pi_prior = prior_weights,
        estimate_sigma2 = estimate_prior_variance,
        max_iter = max_iter, tol = tol,
        stop_when_stable = stop_when_stable, verbose = verbose,
        pre = pre, track_elbo_updates = track_elbo_updates
      )
    }

    candidate_fits[[i]] <- list(
      fit = selected_fit, raw_fit = raw_fit, Rbar = Rbar
    )
    profile_rows[[i]] <- data.frame(
      lambda = lambda,
      raw_eta0 = raw_eta0,
      eta0_used = eta0_used,
      eta0_multiplier = if (estimate_eta0) eta0_multiplier else 1,
      eta0_estimated = estimate_eta0,
      objective = selected_fit$lower_bound,
      converged = isTRUE(selected_fit$converged),
      iterations = length(selected_fit$elbo),
      stringsAsFactors = FALSE
    )
  }

  profile <- do.call(rbind, profile_rows)
  best_objective <- max(profile$objective)
  eligible <- which(profile$objective >= best_objective - 1e-8)
  best <- eligible[which.min(profile$lambda[eligible])]
  chosen <- candidate_fits[[best]]
  profile$selected <- seq_len(nrow(profile)) == best

  result <- list(
    method = "heuristic collapsed-R",
    fit = chosen$fit,
    pip = as.numeric(chosen$fit$pip),
    alpha = chosen$fit$alpha,
    selected_lambda = profile$lambda[best],
    selected_Rbar = chosen$Rbar,
    selected_Rhat = chosen$fit$Rhat,
    selected_Rtilde = chosen$fit$Rtilde,
    raw_eta0 = profile$raw_eta0[best],
    eta0_used = profile$eta0_used[best],
    eta0_multiplier = profile$eta0_multiplier[best],
    eta0_estimated = estimate_eta0,
    L = L,
    profile = profile,
    converged = isTRUE(chosen$fit$converged),
    call = match.call()
  )
  if (keep_candidate_fits) result$candidate_fits <- candidate_fits
  # Inherit from collapsed_r_fit so the shared post-processing helpers unwrap
  # the internal fit exactly as they do for the principal method.
  class(result) <- c("heuristic_collapsed_r_fit", "collapsed_r_fit")
  result
}

print.heuristic_collapsed_r_fit <- function(x, ...) {
  cat("Heuristic collapsed-R fit (cov2cor surrogate)\n")
  cat("  lambda:", format(x$selected_lambda), "\n")
  cat("  eta0:", format(x$eta0_used),
      if (x$eta0_estimated) "(estimated and calibrated)" else "(fixed)",
      "\n")
  if (x$eta0_estimated) {
    cat("  raw eta0:", format(x$raw_eta0),
        " multiplier:", format(x$eta0_multiplier), "\n")
  }
  cat("  L:", x$L, " state-converged:", x$converged, "\n")
  invisible(x)
}

susie_iw_heuristic <- function(
    v = NULL, R0, N, z = NULL,
    L = min(10L, ncol(R0)),
    lambda_grid = collapsed_r_default_lambda_grid(),
    eta0 = NULL, eta0_multiplier = 1.677,
    eta0_bounds = c(1, 1e8), eta0_init = 100,
    eta0_coarse_grid_size = 15L,
    scaled_prior_variance = 0.2,
    estimate_prior_variance = TRUE, prior_weights = NULL,
    coverage = 0.95, min_abs_corr = 0.5,
    ser_fallback = TRUE, partition_ser = FALSE,
    max_partition_sets = 3L, min_partition_probability = 0.1,
    max_iter = 200L, tol = 1e-4, stop_when_stable = TRUE,
    verbose = FALSE, track_fit = FALSE,
    keep_candidate_fits = FALSE) {
  hiw_require_core(require_postprocessing = TRUE)
  supplied_call <- match.call()
  summary_data <- iw_prepare_summary_inputs(v = v, R0 = R0, N = N, z = z)
  method_fit <- fit_collapsed_r_heuristic(
    v = summary_data$v, R0 = summary_data$R0,
    N = summary_data$effective_N, L = L,
    lambda_grid = lambda_grid, eta0 = eta0,
    eta0_multiplier = eta0_multiplier,
    eta0_bounds = eta0_bounds, eta0_init = eta0_init,
    eta0_coarse_grid_size = eta0_coarse_grid_size,
    prior_variance = scaled_prior_variance,
    estimate_prior_variance = estimate_prior_variance,
    prior_weights = prior_weights,
    max_iter = max_iter, tol = tol,
    stop_when_stable = stop_when_stable, verbose = verbose,
    track_elbo_updates = track_fit,
    keep_candidate_fits = keep_candidate_fits
  )
  diagnostics <- list(
    approximation = "heuristic collapsed-R with cov2cor surrogate",
    selected_lambda = method_fit$selected_lambda,
    selected_Rbar = method_fit$selected_Rbar,
    selected_Rhat = method_fit$selected_Rhat,
    selected_Rtilde = method_fit$selected_Rtilde,
    raw_eta0 = method_fit$raw_eta0,
    eta0 = method_fit$eta0_used,
    eta0_multiplier = method_fit$eta0_multiplier,
    eta0_estimated = method_fit$eta0_estimated,
    L = method_fit$L,
    profile = method_fit$profile,
    candidate_fits = method_fit$candidate_fits,
    input = summary_data$input,
    supplied_N = summary_data$supplied_N,
    effective_N = summary_data$effective_N,
    elbo_is_diagnostic = TRUE
  )
  result <- iw_finalize_public_fit(
    method_fit = method_fit, v = summary_data$v,
    N = summary_data$effective_N,
    R = method_fit$selected_Rtilde,
    method = "Heuristic SuSiE-IW (cov2cor collapsed-R)",
    subclass = "susie_iw_heuristic",
    coverage = coverage, min_abs_corr = min_abs_corr,
    ser_fallback = ser_fallback, partition_ser = partition_ser,
    max_partition_sets = max_partition_sets,
    min_partition_probability = min_partition_probability,
    prior_variance = scaled_prior_variance,
    max_iter = max_iter, tol = tol, call = supplied_call,
    diagnostics = diagnostics
  )
  result$selected_lambda <- method_fit$selected_lambda
  result$selected_Rbar <- method_fit$selected_Rbar
  result$selected_Rhat <- method_fit$selected_Rhat
  result$selected_Rtilde <- method_fit$selected_Rtilde
  result$raw_eta0 <- method_fit$raw_eta0
  result$eta0 <- method_fit$eta0_used
  result$eta0_used <- method_fit$eta0_used
  result$eta0_multiplier <- method_fit$eta0_multiplier
  result$L <- method_fit$L
  result
}

print.susie_iw_heuristic <- function(x, ...) {
  cat("Heuristic SuSiE-IW fit (cov2cor collapsed-R)\n")
  cat("  lambda:", format(x$selected_lambda), "\n")
  cat("  raw eta0:", format(x$raw_eta0),
      " eta0 used:", format(x$eta0),
      " multiplier:", format(x$eta0_multiplier), "\n")
  cat("  reported CS:", length(susie_get_reported_cs(x)),
      " SER fallback:", isTRUE(x$postprocess$fallback_used), "\n")
  cat("  state-converged:", isTRUE(x$converged), "\n")
  invisible(x)
}
