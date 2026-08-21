# Maintained collapsed-R implementation for the SuSiE-IW project.
#
# Public functions:
#   susie_iw(v, R0, N, ...)
#   fit_collapsed_r(v, R0, N, ...)  # legacy standardized-input interface
#
# Inputs use the standardized RSS convention
#   v = X' y / N, R0 = X0' X0 / N0, var(y) = 1.
# The implementation otherwise uses only base R.
#
# SuSiE collapsed-R with q(omega).
#
# This implements the multi-effect scalar-augmented approximation: IBSS-style
# coordinate updates for L single effects and a GIG q(omega). The
# public fit profiles eta0 by variational empirical Bayes, or holds it fixed
# when the caller supplies eta0.

scr_logdet_spd <- function(A) {
  as.numeric(determinant(A, logarithm = TRUE)$modulus)
}

scr_softmax <- function(x) {
  z <- x - max(x)
  exp_z <- exp(z)
  exp_z / sum(exp_z)
}

scr_log_bessel_k <- function(x, nu) {
  log(besselK(x, nu, expon.scaled = TRUE)) - x
}

scr_small_chi <- function(lambda, chi) {
  chi <= pmax(1e-10, 1e-4 * lambda^2)
}

scr_gig_inv_mean <- function(lambda, chi) {
  small <- scr_small_chi(lambda, chi)
  if (small) {
    return(1 / (2 * (lambda - 1)))
  }

  x <- sqrt(chi)
  val <- exp(
    scr_log_bessel_k(x, lambda - 1) -
      scr_log_bessel_k(x, lambda) -
      0.5 * log(chi)
  )
  if (is.finite(val)) val else 1 / (2 * (lambda - 1))
}

scr_gig_chisq_kl <- function(lambda, chi, eta) {
  if (scr_small_chi(lambda, chi)) {
    return(0)
  }

  x <- sqrt(chi)
  val <- lgamma(lambda) +
    (lambda - 1) * log(2) -
    lambda / 2 * log(chi) -
    scr_log_bessel_k(x, lambda) -
    chi * eta / 2

  if (is.finite(val)) max(val, 0) else 0
}

# Evaluate the GIG normalizing integral
#
#   E_{omega ~ chi-square(2 lambda)} exp{-chi / (2 omega)}.
#
# The Bessel formula is fast, but base R can overflow when both the order and
# eta0 are large. The numerical fallback integrates over log(omega / 2), with
# the integrand shifted by its analytic mode.
scr_gig_numeric_summary <- function(lambda, chi, relative_tolerance = 1e-9) {
  stopifnot(lambda > 1, chi >= 0)
  if (chi == 0) {
    return(list(log_integral = 0, inverse_mean = 1 / (2 * (lambda - 1))))
  }

  integrate_log_kernel <- function(power) {
    log_integrand <- function(log_y) {
      y <- exp(log_y)
      power * log_y - y - chi / (4 * y)
    }
    mode_y <- (power + sqrt(power^2 + chi)) / 2
    mode_log_y <- log(mode_y)
    maximum <- log_integrand(mode_log_y)
    target <- maximum - 45

    find_bound <- function(direction) {
      bound <- mode_log_y
      step <- 0.5
      for (attempt in seq_len(100L)) {
        candidate <- bound + direction * step
        value <- log_integrand(candidate)
        bound <- candidate
        if (!is.finite(value) || value <= target) return(bound)
        step <- min(4, step * 1.35)
      }
      stop("Could not bracket the GIG log-integral")
    }
    lower <- find_bound(-1)
    upper <- find_bound(1)
    scaled <- function(log_y) {
      exp(vapply(log_y, log_integrand, numeric(1)) - maximum)
    }
    value <- stats::integrate(
      scaled,
      lower = lower,
      upper = upper,
      rel.tol = relative_tolerance,
      subdivisions = 500L,
      stop.on.error = TRUE
    )$value
    log(value) + maximum
  }
  log_denominator <- integrate_log_kernel(lambda)
  log_numerator_without_half <- integrate_log_kernel(lambda - 1)

  list(
    log_integral = log_denominator - lgamma(lambda),
    inverse_mean = 0.5 * exp(
      log_numerator_without_half - log_denominator
    )
  )
}

scr_gig_log_integral <- function(lambda, chi) {
  stopifnot(lambda > 1, chi >= 0)
  if (chi == 0) return(0)

  # When chi is small relative to lambda^2, direct evaluation of the
  # high-order Bessel function loses precision.  With
  # omega ~ chi-square(2 * lambda), a second-order cumulant expansion of
  # log E exp{-chi / (2 * omega)} is accurate in this regime and has the
  # correct normal-limit behavior as lambda grows.
  if (scr_small_chi(lambda, chi)) {
    first_order <- -chi / (4 * (lambda - 1))
    second_order <- if (lambda > 2) {
      chi^2 / (32 * (lambda - 1)^2 * (lambda - 2))
    } else {
      0
    }
    return(min(first_order + second_order, 0))
  }

  value <- log(2) + lambda / 2 * log(chi) +
    scr_log_bessel_k(sqrt(chi), lambda) -
    lambda * log(2) - lgamma(lambda)
  if (is.finite(value) && value <= 1e-8) return(min(value, 0))
  scr_gig_numeric_summary(lambda, chi)$log_integral
}

scr_gig_inv_mean_stable <- function(lambda, chi) {
  stopifnot(lambda > 1, chi >= 0)
  if (chi == 0) return(1 / (2 * (lambda - 1)))
  if (scr_small_chi(lambda, chi)) {
    return(1 / (2 * (lambda - 1)))
  }

  value <- exp(
    scr_log_bessel_k(sqrt(chi), lambda - 1) -
      scr_log_bessel_k(sqrt(chi), lambda) -
      0.5 * log(chi)
  )
  if (is.finite(value) && value > 0) return(value)
  scr_gig_numeric_summary(lambda, chi)$inverse_mean
}

scr_gig_chisq_kl_from_log_integral <- function(
    lambda, chi, inverse_mean, log_integral = NULL) {
  if (chi == 0) return(0)
  if (is.null(log_integral)) {
    log_integral <- scr_gig_log_integral(lambda, chi)
  }
  # At the coordinate optimum,
  # log integral = -chi E(1 / omega) / 2 - KL(q || prior).
  max(0, -chi * inverse_mean / 2 - log_integral)
}

scr_log_p0 <- function(x, Rbar, N, nu0, pre = NULL) {
  J <- length(x)
  if (is.null(pre)) {
    Omega <- solve(Rbar)
    r0 <- as.numeric(crossprod(x, Omega %*% x))
    logdetRbar <- scr_logdet_spd(Rbar)
  } else {
    r0 <- pre$r0
    logdetRbar <- pre$logdetRbar
  }

  logdetSigma <- J * log(nu0 / (N * (nu0 + 2))) + logdetRbar
  lgamma((nu0 + J + 2) / 2) -
    lgamma((nu0 + 2) / 2) -
    J / 2 * log((nu0 + 2) * pi) -
    0.5 * logdetSigma -
    (nu0 + J + 2) / 2 * log1p(N * r0 / nu0)
}

scr_initialize_effect <- function(J, sigma2, pi_prior) {
  alpha <- pi_prior
  mu <- rep(0, J)
  s2 <- rep(sigma2, J)
  m2 <- mu^2 + s2
  m <- alpha * mu

  list(
    alpha = alpha,
    mu = mu,
    s2 = s2,
    m2 = m2,
    m = m,
    score = rep(NA_real_, J)
  )
}

scr_update_effect <- function(effect, m_minus, x, A, diagA, N, sigma2,
                              eta, pi_prior) {
  cvec <- as.numeric(A %*% m_minus)
  natural_mean <- x - eta * cvec
  s2 <- 1 / (1 / sigma2 + N * eta * diagA)
  mu <- s2 * N * natural_mean
  m2 <- mu^2 + s2

  score <- 0.5 * (log(s2 / sigma2) + 1 - m2 / sigma2) +
    N * mu * natural_mean -
    0.5 * N * eta * diagA * m2
  alpha <- scr_softmax(log(pi_prior) + score)
  names(alpha) <- paste0("j", seq_along(alpha))

  list(
    alpha = alpha,
    mu = mu,
    s2 = s2,
    m2 = m2,
    m = as.numeric(alpha) * mu,
    score = score
  )
}

scr_effect_second_moment <- function(effects, A) {
  m_total <- Reduce(`+`, lapply(effects, `[[`, "m"))
  sq <- as.numeric(crossprod(m_total, A %*% m_total))
  diagA <- diag(A)

  for (eff in effects) {
    sq <- sq +
      sum(as.numeric(eff$alpha) * eff$m2 * diagA) -
      as.numeric(crossprod(eff$m, A %*% eff$m))
  }

  sq
}

scr_effect_rank_one_second_moment <- function(effects, x) {
  x <- as.numeric(x)
  m_total <- Reduce(`+`, lapply(effects, `[[`, "m"))
  sq <- sum(x * m_total)^2

  for (eff in effects) {
    sq <- sq +
      sum(as.numeric(eff$alpha) * eff$m2 * x^2) -
      sum(x * eff$m)^2
  }

  max(0, sq)
}

scr_profiled_omega_block <- function(
    nu0, x, Rbar, N, effects, pre, reference_second_moment = NULL,
    marginal_second_moment = NULL) {
  if (is.null(reference_second_moment)) {
    reference_second_moment <- scr_effect_second_moment(effects, Rbar)
  }
  if (is.null(marginal_second_moment)) {
    marginal_second_moment <- scr_effect_rank_one_second_moment(effects, x)
  }
  chi <- N * (
    nu0 * reference_second_moment + N * marginal_second_moment
  )
  chi <- max(0, chi)
  lambda <- (nu0 + 3) / 2
  log_integral <- scr_gig_log_integral(lambda, chi)

  list(
    objective = scr_log_p0(x, Rbar, N, nu0, pre) + log_integral,
    lambda = lambda,
    chi = chi,
    log_integral = log_integral,
    reference_second_moment = reference_second_moment,
    marginal_second_moment = marginal_second_moment
  )
}

scr_optimize_profiled_nu0 <- function(
    current_nu0, x, Rbar, N, effects, pre,
    nu0_bounds = c(1, 1e5), coarse_grid_size = 21L) {
  stopifnot(
    length(nu0_bounds) == 2L,
    all(is.finite(nu0_bounds)),
    nu0_bounds[1] > 0,
    nu0_bounds[2] > nu0_bounds[1],
    coarse_grid_size >= 5L
  )
  current_nu0 <- min(max(current_nu0, nu0_bounds[1]), nu0_bounds[2])
  reference_second_moment <- scr_effect_second_moment(effects, Rbar)
  marginal_second_moment <- scr_effect_rank_one_second_moment(effects, x)
  cache <- new.env(parent = emptyenv())

  evaluate <- function(log_nu0) {
    nu0 <- exp(log_nu0)
    key <- sprintf("%.16g", nu0)
    if (!exists(key, envir = cache, inherits = FALSE)) {
      assign(
        key,
        scr_profiled_omega_block(
          nu0 = nu0,
          x = x,
          Rbar = Rbar,
          N = N,
          effects = effects,
          pre = pre,
          reference_second_moment = reference_second_moment,
          marginal_second_moment = marginal_second_moment
        ),
        envir = cache
      )
    }
    get(key, envir = cache, inherits = FALSE)
  }

  log_grid <- sort(unique(c(
    seq(
      log(nu0_bounds[1]),
      log(nu0_bounds[2]),
      length.out = as.integer(coarse_grid_size)
    ),
    log(current_nu0)
  )))
  grid_objective <- vapply(
    log_grid,
    function(value) evaluate(value)$objective,
    numeric(1)
  )
  candidate_log_nu0 <- log_grid
  candidate_objective <- grid_objective

  local_maximum <- which(
    grid_objective >= c(-Inf, head(grid_objective, -1L)) &
      grid_objective >= c(tail(grid_objective, -1L), -Inf)
  )
  for (index in local_maximum) {
    lower_index <- max(1L, index - 1L)
    upper_index <- min(length(log_grid), index + 1L)
    if (lower_index == upper_index) next
    refined <- stats::optimize(
      function(log_nu0) -evaluate(log_nu0)$objective,
      interval = log_grid[c(lower_index, upper_index)],
      tol = 1e-7
    )
    candidate_log_nu0 <- c(candidate_log_nu0, refined$minimum)
    candidate_objective <- c(candidate_objective, -refined$objective)
  }

  best <- which.max(candidate_objective)
  selected_log_nu0 <- candidate_log_nu0[best]
  selected <- evaluate(selected_log_nu0)
  selected$nu0 <- exp(selected_log_nu0)
  boundary_tolerance <- 1e-4
  selected$boundary <- if (
    abs(selected_log_nu0 - log(nu0_bounds[1])) < boundary_tolerance
  ) {
    "lower"
  } else if (
    abs(selected_log_nu0 - log(nu0_bounds[2])) < boundary_tolerance
  ) {
    "upper"
  } else {
    "interior"
  }
  selected$coarse_grid <- data.frame(
    nu0 = exp(log_grid),
    objective = grid_objective
  )
  selected$number_of_objective_evaluations <- length(ls(cache))
  selected
}

scr_lower_bound <- function(x, Rbar, N, nu0, sigma2, effects, A, eta,
                            omega_kl, pi_prior, pre) {
  m_total <- Reduce(`+`, lapply(effects, `[[`, "m"))
  S_Q <- scr_effect_second_moment(effects, A)
  out <- scr_log_p0(x, Rbar, N, nu0, pre) +
    N * sum(x * m_total) -
    0.5 * N * eta * S_Q -
    omega_kl

  for (ell in seq_along(effects)) {
    eff <- effects[[ell]]
    active <- as.numeric(eff$alpha) > 0
    normal_term <- 0.5 * (
      log(eff$s2 / sigma2[ell]) + 1 - eff$m2 / sigma2[ell]
    )
    out <- out + sum(as.numeric(eff$alpha)[active] * (
      log(pi_prior[active]) - log(as.numeric(eff$alpha)[active]) +
        normal_term[active]
    ))
  }

  out
}

fit_susie_collapse_r_shared_fixed_nu0 <- function(
    x, Rbar, N, L = 2, sigma2 = 0.2^2, nu0,
    pi_prior = NULL, estimate_sigma2 = FALSE,
    max_iter = 100, tol = 1e-6, verbose = FALSE, pre = NULL,
    track_elbo_updates = FALSE) {
  x <- as.numeric(x)
  Rbar <- as.matrix(Rbar)
  J <- length(x)
  if (length(sigma2) == 1) {
    sigma2 <- rep(sigma2, L)
  }
  stopifnot(
    nrow(Rbar) == J, ncol(Rbar) == J,
    N > 0, L >= 1, length(sigma2) == L,
    all(sigma2 > 0), nu0 > 0,
    max_iter >= 1, tol > 0
  )

  if (is.null(pi_prior)) {
    pi_prior <- rep(1 / J, J)
  }
  pi_prior <- as.numeric(pi_prior)
  pi_prior <- pi_prior / sum(pi_prior)
  if (any(pi_prior <= 0)) {
    stop("All prior weights must be positive.")
  }

  if (is.null(pre)) {
    Omega <- solve(Rbar)
    pre <- list(
      r0 = as.numeric(crossprod(x, Omega %*% x)),
      logdetRbar = scr_logdet_spd(Rbar)
    )
  }

  A <- nu0 * Rbar + N * tcrossprod(x)
  diagA <- diag(A)
  lambda <- (nu0 + 3) / 2
  eta <- 1 / (nu0 + 1)
  chi <- NA_real_
  omega_kl <- 0
  effects <- lapply(seq_len(L), function(ell) {
    scr_initialize_effect(J, sigma2[ell], pi_prior)
  })
  elbo <- numeric()
  elbo_updates <- data.frame(
    update_index = integer(),
    iteration = integer(),
    update = character(),
    effect = integer(),
    lower_bound = numeric()
  )
  record_elbo_update <- function(iteration, update, effect = NA_integer_) {
    if (!isTRUE(track_elbo_updates)) {
      return(invisible(NULL))
    }
    elbo_updates <<- rbind(
      elbo_updates,
      data.frame(
        update_index = nrow(elbo_updates) + 1L,
        iteration = as.integer(iteration),
        update = update,
        effect = as.integer(effect),
        lower_bound = scr_lower_bound(
          x, Rbar, N, nu0, sigma2, effects, A, eta,
          omega_kl, pi_prior, pre
        )
      )
    )
    invisible(NULL)
  }
  record_elbo_update(iteration = 0L, update = "initial")

  for (iter in seq_len(max_iter)) {
    m_sum <- Reduce(`+`, lapply(effects, `[[`, "m"))

    for (ell in seq_len(L)) {
      m_minus <- m_sum - effects[[ell]]$m
      old_m <- effects[[ell]]$m
      effects[[ell]] <- scr_update_effect(
        effects[[ell]], m_minus, x, A, diagA, N, sigma2[ell], eta, pi_prior
      )
      m_sum <- m_sum - old_m + effects[[ell]]$m
      record_elbo_update(
        iteration = iter,
        update = "effect",
        effect = ell
      )
    }

    if (estimate_sigma2) {
      sigma2 <- vapply(effects, function(eff) {
        max(sum(as.numeric(eff$alpha) * eff$m2), 1e-12)
      }, numeric(1))
      record_elbo_update(iteration = iter, update = "prior_variance")
    }

    S_Q <- scr_effect_second_moment(effects, A)
    chi <- N * S_Q
    log_integral <- scr_gig_log_integral(lambda, chi)
    eta <- scr_gig_inv_mean_stable(lambda, chi)
    omega_kl <- scr_gig_chisq_kl_from_log_integral(
      lambda = lambda,
      chi = chi,
      inverse_mean = eta,
      log_integral = log_integral
    )
    elbo[iter] <- scr_lower_bound(
      x, Rbar, N, nu0, sigma2, effects, A, eta, omega_kl, pi_prior, pre
    )
    record_elbo_update(iteration = iter, update = "omega")

    if (verbose && (iter == 1 || iter %% 10 == 0)) {
      tops <- vapply(effects, function(eff) which.max(eff$alpha), integer(1))
      pips <- vapply(effects, function(eff) max(eff$alpha), numeric(1))
      message(sprintf(
        "nu0=%.4f iter=%d elbo=%.6f eta=%.4g top=(%s) pi=(%s)",
        nu0, iter, elbo[iter], eta,
        paste(tops, collapse = ","),
        paste(sprintf("%.3f", pips), collapse = ",")
      ))
    }

    elbo_diff <- abs(elbo[iter] - elbo[iter - 1])
    if (iter > 1 && is.finite(elbo_diff) &&
        elbo_diff < tol * (1 + abs(elbo[iter - 1]))) {
      elbo <- elbo[seq_len(iter)]
      break
    }
  }

  alpha <- do.call(rbind, lapply(effects, function(eff) as.numeric(eff$alpha)))
  rownames(alpha) <- paste0("effect", seq_len(L))
  colnames(alpha) <- paste0("j", seq_len(J))
  pip <- 1 - apply(1 - alpha, 2, prod)

  list(
    alpha = alpha,
    pip = pip,
    gamma_hat = apply(alpha, 1, which.max),
    nu0 = nu0,
    eta = eta,
    chi = chi,
    omega_kl = omega_kl,
    lower_bound = tail(elbo, 1),
    elbo = elbo,
    elbo_updates = if (isTRUE(track_elbo_updates)) elbo_updates else NULL,
    effects = effects,
    sigma2 = sigma2,
    N = N,
    Rbar = Rbar,
    x = x,
    converged = length(elbo) < max_iter
  )
}

fit_susie_collapse_r_shared_veb <- function(
    x, Rbar, N, L = 2, sigma2 = 0.2^2,
    nu0_init = 100, nu0_bounds = c(1, 1e5),
    nu0_coarse_grid_size = 21L,
    pi_prior = NULL, estimate_sigma2 = FALSE,
    max_iter = 100, tol = 1e-6, verbose = FALSE,
    track_elbo_updates = FALSE, pre = NULL) {
  x <- as.numeric(x)
  Rbar <- as.matrix(Rbar)
  J <- length(x)
  if (length(sigma2) == 1L) sigma2 <- rep(sigma2, L)
  stopifnot(
    nrow(Rbar) == J, ncol(Rbar) == J,
    N > 0, L >= 1, length(sigma2) == L,
    all(sigma2 > 0), nu0_init > 0,
    length(nu0_bounds) == 2L,
    nu0_bounds[1] > 0,
    nu0_bounds[2] > nu0_bounds[1],
    max_iter >= 1, tol > 0
  )

  if (is.null(pi_prior)) pi_prior <- rep(1 / J, J)
  pi_prior <- as.numeric(pi_prior)
  pi_prior <- pi_prior / sum(pi_prior)
  if (any(pi_prior <= 0)) stop("All prior weights must be positive.")

  if (is.null(pre)) {
    Omega <- solve(Rbar)
    pre <- list(
      r0 = as.numeric(crossprod(x, Omega %*% x)),
      logdetRbar = scr_logdet_spd(Rbar)
    )
  }
  nu0 <- min(max(nu0_init, nu0_bounds[1]), nu0_bounds[2])
  A <- nu0 * Rbar + N * tcrossprod(x)
  diagA <- diag(A)
  lambda <- (nu0 + 3) / 2
  inverse_mean <- 1 / (nu0 + 1)
  chi <- 0
  omega_kl <- 0
  effects <- lapply(seq_len(L), function(ell) {
    scr_initialize_effect(J, sigma2[ell], pi_prior)
  })
  elbo <- numeric()
  nu0_trace <- numeric()
  nu0_boundary_trace <- character()
  objective_evaluations <- integer()
  elbo_updates <- data.frame(
    update_index = integer(),
    iteration = integer(),
    update = character(),
    effect = integer(),
    nu0 = numeric(),
    lower_bound = numeric()
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
        )
      )
    )
    invisible(NULL)
  }
  record_elbo_update(iteration = 0L, update = "initial")

  for (iter in seq_len(max_iter)) {
    m_sum <- Reduce(`+`, lapply(effects, `[[`, "m"))
    for (ell in seq_len(L)) {
      m_minus <- m_sum - effects[[ell]]$m
      old_m <- effects[[ell]]$m
      effects[[ell]] <- scr_update_effect(
        effects[[ell]], m_minus, x, A, diagA, N, sigma2[ell],
        inverse_mean, pi_prior
      )
      m_sum <- m_sum - old_m + effects[[ell]]$m
      record_elbo_update(iteration = iter, update = "effect", effect = ell)
    }

    if (estimate_sigma2) {
      sigma2 <- vapply(effects, function(effect) {
        max(sum(as.numeric(effect$alpha) * effect$m2), 1e-12)
      }, numeric(1))
      record_elbo_update(iteration = iter, update = "prior_variance")
    }

    profiled <- scr_optimize_profiled_nu0(
      current_nu0 = nu0,
      x = x,
      Rbar = Rbar,
      N = N,
      effects = effects,
      pre = pre,
      nu0_bounds = nu0_bounds,
      coarse_grid_size = nu0_coarse_grid_size
    )
    nu0 <- profiled$nu0
    lambda <- profiled$lambda
    chi <- profiled$chi
    inverse_mean <- scr_gig_inv_mean_stable(lambda, chi)
    omega_kl <- scr_gig_chisq_kl_from_log_integral(
      lambda = lambda,
      chi = chi,
      inverse_mean = inverse_mean,
      log_integral = profiled$log_integral
    )
    A <- nu0 * Rbar + N * tcrossprod(x)
    diagA <- diag(A)
    nu0_trace[iter] <- nu0
    nu0_boundary_trace[iter] <- profiled$boundary
    objective_evaluations[iter] <-
      profiled$number_of_objective_evaluations
    elbo[iter] <- scr_lower_bound(
      x, Rbar, N, nu0, sigma2, effects, A, inverse_mean,
      omega_kl, pi_prior, pre
    )
    record_elbo_update(iteration = iter, update = "eta0_omega")

    if (verbose && (iter == 1L || iter %% 10L == 0L)) {
      tops <- vapply(effects, function(effect) {
        which.max(effect$alpha)
      }, integer(1))
      message(sprintf(
        paste0(
          "iter=%d elbo=%.6f eta0=%.5g E(1/omega)=%.5g ",
          "boundary=%s top=(%s)"
        ),
        iter, elbo[iter], nu0, inverse_mean, profiled$boundary,
        paste(tops, collapse = ",")
      ))
    }

    if (iter > 1L) {
      elbo_difference <- abs(elbo[iter] - elbo[iter - 1L])
      if (is.finite(elbo_difference) &&
          elbo_difference < tol * (1 + abs(elbo[iter - 1L]))) {
        elbo <- elbo[seq_len(iter)]
        nu0_trace <- nu0_trace[seq_len(iter)]
        nu0_boundary_trace <- nu0_boundary_trace[seq_len(iter)]
        objective_evaluations <- objective_evaluations[seq_len(iter)]
        break
      }
    }
  }

  alpha <- do.call(rbind, lapply(effects, function(effect) {
    as.numeric(effect$alpha)
  }))
  rownames(alpha) <- paste0("effect", seq_len(L))
  colnames(alpha) <- paste0("j", seq_len(J))
  pip <- 1 - apply(1 - alpha, 2, prod)

  list(
    alpha = alpha,
    pip = pip,
    gamma_hat = apply(alpha, 1, which.max),
    nu0 = nu0,
    nu0_init = nu0_init,
    nu0_bounds = nu0_bounds,
    nu0_boundary = tail(nu0_boundary_trace, 1L),
    nu0_trace = nu0_trace,
    eta = inverse_mean,
    chi = chi,
    omega_kl = omega_kl,
    lower_bound = tail(elbo, 1L),
    elbo = elbo,
    elbo_updates = if (isTRUE(track_elbo_updates)) elbo_updates else NULL,
    effects = effects,
    sigma2 = sigma2,
    N = N,
    Rbar = Rbar,
    x = x,
    converged = length(elbo) < max_iter,
    objective_evaluations = objective_evaluations
  )
}

fit_susie_collapse_r_shared <- function(
    x, Rbar, N, L = 2, sigma2 = 0.2^2,
    nu0_grid = exp(seq(log(5), log(5000), length.out = 61)),
    pi_prior = NULL, estimate_sigma2 = FALSE,
    max_iter = 100, tol = 1e-6, verbose = FALSE,
    track_elbo_updates = FALSE) {
  x <- as.numeric(x)
  Rbar <- as.matrix(Rbar)
  J <- length(x)
  stopifnot(
    nrow(Rbar) == J, ncol(Rbar) == J,
    N > 0, L >= 1,
    length(nu0_grid) >= 1,
    all(nu0_grid > 0)
  )

  Omega <- solve(Rbar)
  pre <- list(
    r0 = as.numeric(crossprod(x, Omega %*% x)),
    logdetRbar = scr_logdet_spd(Rbar)
  )

  fits <- lapply(sort(unique(as.numeric(nu0_grid))), function(nu0) {
    fit_susie_collapse_r_shared_fixed_nu0(
      x = x,
      Rbar = Rbar,
      N = N,
      L = L,
      sigma2 = sigma2,
      nu0 = nu0,
      pi_prior = pi_prior,
      estimate_sigma2 = estimate_sigma2,
      max_iter = max_iter,
      tol = tol,
      verbose = verbose,
      pre = pre,
      track_elbo_updates = track_elbo_updates
    )
  })

  bounds <- vapply(fits, `[[`, numeric(1), "lower_bound")
  best <- which.max(bounds)
  fit <- fits[[best]]
  fit$nu0_grid <- vapply(fits, `[[`, numeric(1), "nu0")
  fit$grid_lower_bound <- bounds
  fit$all_grid_fits <- fits
  fit
}

scr_sample_inverse_wishart <- function(df, scale) {
  W <- stats::rWishart(1, df = df, Sigma = solve(scale))[, , 1]
  solve(W)
}

# -------------------------------------------------------------------------
# Public collaborator-facing interface
# -------------------------------------------------------------------------

collapsed_r_default_lambda_grid <- function() {
  c(1e-5, 1e-4, 1e-3, 2e-3, 4e-3, 6e-3, 8e-3, 1e-2)
}

collapsed_r_validate_inputs <- function(v, R0, N, L, lambda_grid) {
  v <- as.numeric(v)
  R0 <- as.matrix(R0)
  J <- length(v)
  if (J < 2L || any(!is.finite(v))) {
    stop("v must be a finite vector containing at least two variants")
  }
  if (!all(dim(R0) == c(J, J)) || any(!is.finite(R0))) {
    stop("R0 must be a finite length(v) by length(v) matrix")
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
  R0 <- (R0 + t(R0)) / 2
  list(v = v, R0 = R0, N = as.numeric(N), L = L,
       lambda_grid = lambda_grid)
}

# Fit collapsed-R while selecting lambda by the calibrated ELBO.
#
# eta0 = NULL:
#   estimate eta0 separately for every lambda, multiply each estimate by
#   eta0_multiplier, refit at the corrected value, and compare the refitted
#   ELBOs. The default multiplier 1.677 is the frozen calibration constant.
#
# eta0 = a positive number:
#   hold eta0 fixed at exactly that number for every lambda. No eta0 update or
#   multiplicative correction is performed.
#
# L is the fixed number of single-effect components. lambda_grid may be
# replaced by the user. Set verbose=TRUE for lambda- and iteration-level ELBO
# progress, and keep_candidate_fits=TRUE only when the complete profile fits
# are needed for diagnostics.
fit_collapsed_r <- function(
    v, R0, N, L = 5L,
    lambda_grid = collapsed_r_default_lambda_grid(),
    eta0 = NULL,
    eta0_multiplier = 1.677,
    eta0_bounds = c(1, 1e8),
    eta0_init = 100,
    eta0_coarse_grid_size = 15L,
    prior_variance = 0.2,
    estimate_prior_variance = TRUE,
    prior_weights = NULL,
    max_iter = 200L,
    tol = 1e-4,
    verbose = FALSE,
    track_elbo_updates = FALSE,
    keep_candidate_fits = FALSE) {
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
        "collapsed-R: lambda=%g (%d/%d)", lambda, i, length(lambda_grid)
      ))
    }
    if (estimate_eta0) {
      raw_fit <- fit_susie_collapse_r_shared_veb(
        x = v, Rbar = Rbar, N = N, L = L,
        sigma2 = prior_variance,
        nu0_init = eta0_init,
        nu0_bounds = eta0_bounds,
        nu0_coarse_grid_size = eta0_coarse_grid_size,
        pi_prior = prior_weights,
        estimate_sigma2 = estimate_prior_variance,
        max_iter = max_iter, tol = tol, verbose = verbose,
        track_elbo_updates = track_elbo_updates, pre = pre
      )
      raw_eta0 <- raw_fit$nu0
      eta0_used <- min(eta0_bounds[2L], eta0_multiplier * raw_eta0)
      if (isTRUE(all.equal(eta0_used, raw_eta0))) {
        selected_fit <- raw_fit
      } else {
        selected_fit <- fit_susie_collapse_r_shared_fixed_nu0(
          x = v, Rbar = Rbar, N = N, L = L,
          sigma2 = prior_variance,
          nu0 = eta0_used,
          pi_prior = prior_weights,
          estimate_sigma2 = estimate_prior_variance,
          max_iter = max_iter, tol = tol, verbose = verbose,
          track_elbo_updates = track_elbo_updates, pre = pre
        )
      }
    } else {
      raw_fit <- NULL
      raw_eta0 <- as.numeric(eta0)
      eta0_used <- as.numeric(eta0)
      selected_fit <- fit_susie_collapse_r_shared_fixed_nu0(
        x = v, Rbar = Rbar, N = N, L = L,
        sigma2 = prior_variance,
        nu0 = eta0_used,
        pi_prior = prior_weights,
        estimate_sigma2 = estimate_prior_variance,
        max_iter = max_iter, tol = tol, verbose = verbose,
        track_elbo_updates = track_elbo_updates, pre = pre
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
    if (verbose) {
      message(sprintf(
        "  objective=%.6f raw_eta0=%.4g eta0_used=%.4g",
        selected_fit$lower_bound, raw_eta0, eta0_used
      ))
    }
  }

  profile <- do.call(rbind, profile_rows)
  best_objective <- max(profile$objective)
  eligible <- which(profile$objective >= best_objective - 1e-8)
  best <- eligible[which.min(profile$lambda[eligible])]
  chosen <- candidate_fits[[best]]
  profile$selected <- seq_len(nrow(profile)) == best

  result <- list(
    method = "collapsed-R",
    fit = chosen$fit,
    pip = as.numeric(chosen$fit$pip),
    alpha = chosen$fit$alpha,
    selected_lambda = profile$lambda[best],
    selected_Rbar = chosen$Rbar,
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
  class(result) <- "collapsed_r_fit"
  result
}

print.collapsed_r_fit <- function(x, ...) {
  cat("Collapsed-R fit\n")
  cat("  lambda:", format(x$selected_lambda), "\n")
  cat("  eta0:", format(x$eta0_used),
      if (x$eta0_estimated) "(estimated and calibrated)" else "(fixed)",
      "\n")
  if (x$eta0_estimated) {
    cat("  raw eta0:", format(x$raw_eta0),
        " multiplier:", format(x$eta0_multiplier), "\n")
  }
  cat("  L:", x$L, " converged:", x$converged, "\n")
  invisible(x)
}

# Paper-facing interface to collapsed-R. The primary inputs are the marginal
# association v = X'y/N, reference LD R0, and GWAS sample size N. Marginal
# z-scores are accepted as an explicitly named alternative. The top-level
# posterior fields are the final solution after the manuscript SER fallback;
# the original collapsed-R fit is retained unchanged in $raw_fit.
susie_iw <- function(
    v = NULL, R0, N, z = NULL,
    L = min(10L, ncol(R0)),
    lambda_grid = collapsed_r_default_lambda_grid(),
    eta0 = NULL,
    eta0_multiplier = 1.677,
    eta0_bounds = c(1, 1e8),
    eta0_init = 100,
    eta0_coarse_grid_size = 15L,
    scaled_prior_variance = 0.2,
    estimate_prior_variance = TRUE,
    prior_weights = NULL,
    coverage = 0.95,
    min_abs_corr = 0.5,
    ser_fallback = TRUE,
    partition_ser = FALSE,
    max_partition_sets = 3L,
    min_partition_probability = 0.1,
    max_iter = 200L,
    tol = 1e-4,
    verbose = FALSE,
    track_fit = FALSE,
    keep_candidate_fits = FALSE) {
  if (!exists("iw_prepare_summary_inputs", mode = "function") ||
      !exists("iw_finalize_public_fit", mode = "function")) {
    stop("Source 03-SER-fallback-and-plot.R before calling susie_iw()")
  }
  supplied_call <- match.call()
  summary_data <- iw_prepare_summary_inputs(
    v = v, R0 = R0, N = N, z = z
  )
  method_fit <- fit_collapsed_r(
    v = summary_data$v, R0 = summary_data$R0,
    N = summary_data$effective_N, L = L,
    lambda_grid = lambda_grid, eta0 = eta0,
    eta0_multiplier = eta0_multiplier,
    eta0_bounds = eta0_bounds, eta0_init = eta0_init,
    eta0_coarse_grid_size = eta0_coarse_grid_size,
    prior_variance = scaled_prior_variance,
    estimate_prior_variance = estimate_prior_variance,
    prior_weights = prior_weights,
    max_iter = max_iter, tol = tol, verbose = verbose,
    track_elbo_updates = track_fit,
    keep_candidate_fits = keep_candidate_fits
  )
  diagnostics <- list(
    approximation = "collapsed-R",
    selected_lambda = method_fit$selected_lambda,
    selected_Rbar = method_fit$selected_Rbar,
    raw_eta0 = method_fit$raw_eta0,
    eta0 = method_fit$eta0_used,
    eta0_multiplier = method_fit$eta0_multiplier,
    eta0_estimated = method_fit$eta0_estimated,
    L = method_fit$L,
    profile = method_fit$profile,
    candidate_fits = method_fit$candidate_fits,
    input = summary_data$input,
    supplied_N = summary_data$supplied_N,
    effective_N = summary_data$effective_N
  )
  result <- iw_finalize_public_fit(
    method_fit = method_fit, v = summary_data$v,
    N = summary_data$effective_N,
    R = method_fit$selected_Rbar,
    method = "SuSiE-IW (collapsed-R)", subclass = "susie_iw",
    coverage = coverage, min_abs_corr = min_abs_corr,
    ser_fallback = ser_fallback,
    partition_ser = partition_ser,
    max_partition_sets = max_partition_sets,
    min_partition_probability = min_partition_probability,
    prior_variance = scaled_prior_variance,
    max_iter = max_iter, tol = tol, call = supplied_call,
    diagnostics = diagnostics
  )
  result$selected_lambda <- method_fit$selected_lambda
  result$selected_Rbar <- method_fit$selected_Rbar
  result$raw_eta0 <- method_fit$raw_eta0
  result$eta0 <- method_fit$eta0_used
  result$eta0_used <- method_fit$eta0_used
  result$eta0_multiplier <- method_fit$eta0_multiplier
  result$L <- method_fit$L
  result
}

print.susie_iw <- function(x, ...) {
  cat("SuSiE-IW (collapsed-R) fit\n")
  cat("  lambda:", format(x$selected_lambda), "\n")
  cat("  raw eta0:", format(x$raw_eta0),
      " eta0 used:", format(x$eta0),
      " multiplier:", format(x$eta0_multiplier), "\n")
  cat("  reported CS:", length(susie_get_reported_cs(x)),
      " SER fallback:", isTRUE(x$postprocess$fallback_used), "\n")
  cat("  converged:", isTRUE(x$converged), "\n")
  invisible(x)
}
