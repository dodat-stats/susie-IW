# Shared SuSiE-RSS, SER fallback, and plotting utilities
#
# This file depends only on susieR. It can be sourced after either method file,
# but does not require them: objects are recognized by their common fields.

iw_require_susieR <- function() {
  if (!requireNamespace("susieR", quietly = TRUE)) {
    stop("Install the susieR package before fitting these models")
  }
}

# susieR <= 0.14 uses susie_suff_stat() for this interface; susieR >= 0.16
# renamed the exported function to susie_ss().
iw_susie_ss_function <- function() {
  iw_require_susieR()
  exports <- getNamespaceExports("susieR")
  if ("susie_suff_stat" %in% exports) {
    return(getExportedValue("susieR", "susie_suff_stat"))
  }
  if ("susie_ss" %in% exports) {
    return(getExportedValue("susieR", "susie_ss"))
  }
  stop("susieR does not export susie_suff_stat() or susie_ss()")
}

# Validate either the paper's marginal association v = X'y/N or a marginal
# z-score vector. The v interface uses N exactly. The z interface mirrors the
# finite-sample adjustment in susieR::susie_rss and therefore has working
# sufficient statistics XtX = (N - 1) R0 and yty = N - 1.
iw_prepare_summary_inputs <- function(v = NULL, R0, N, z = NULL) {
  iw_require_susieR()
  R0 <- as.matrix(R0)
  has_v <- !is.null(v)
  has_z <- !is.null(z)
  if (has_v == has_z) {
    stop("Provide exactly one of v and z")
  }
  summary_statistic <- if (has_v) as.numeric(v) else as.numeric(z)
  if (!length(summary_statistic) || any(!is.finite(summary_statistic))) {
    stop(if (has_v) "v must be a nonempty finite vector" else
      "z must be a nonempty finite vector")
  }
  if (!all(dim(R0) == c(length(summary_statistic),
                        length(summary_statistic))) ||
      any(!is.finite(R0))) {
    stop("R0 must be a finite square matrix matching the summary statistic")
  }
  if (length(N) != 1L || !is.finite(N) || N <= 0) {
    stop("N must be one positive finite number")
  }

  if (has_v) {
    internal_v <- summary_statistic
    adjusted_z <- NULL
    effective_N <- as.numeric(N)
    input <- "v"
  } else {
    if (N <= 2) stop("N must be greater than 2 when z is supplied")
    adjusted_z <- summary_statistic * sqrt(
      (N - 1) / (summary_statistic^2 + N - 2)
    )
    effective_N <- as.numeric(N - 1)
    internal_v <- adjusted_z / sqrt(effective_N)
    input <- "z"
  }

  list(
    input = input,
    v = internal_v,
    z = if (has_z) summary_statistic else NULL,
    adjusted_z = adjusted_z,
    R0 = (R0 + t(R0)) / 2,
    supplied_N = as.numeric(N),
    effective_N = effective_N
  )
}

# Backward-compatible helper used by earlier collaborator scripts.
iw_prepare_rss_inputs <- function(
    z = NULL, R, n, bhat = NULL, shat = NULL, var_y = NULL) {
  if (is.null(z)) {
    if (is.null(bhat) || is.null(shat)) {
      stop("Provide z or both bhat and shat")
    }
    bhat <- as.numeric(bhat)
    shat <- as.numeric(shat)
    if (length(shat) == 1L) shat <- rep(shat, length(bhat))
    if (length(bhat) != length(shat) || any(!is.finite(bhat)) ||
        any(!is.finite(shat)) || any(shat <= 0)) {
      stop("bhat and shat must be finite, equally sized, and shat > 0")
    }
    z <- bhat / shat
  }
  prepared <- iw_prepare_summary_inputs(z = z, R0 = R, N = n)
  list(
    z = prepared$z,
    adjusted_z = prepared$adjusted_z,
    R = prepared$R0,
    n = prepared$supplied_N,
    N = prepared$effective_N,
    v = prepared$v,
    var_y = if (is.null(var_y)) 1 else as.numeric(var_y)
  )
}

# Inverse of the finite-sample conversion above. This is mainly useful for
# migrating analyses that stored v = X'y/(n - 1) instead of the original
# marginal z-scores.
iw_z_from_v <- function(v, n) {
  v <- as.numeric(v)
  if (length(n) != 1L || !is.finite(n) || n <= 2 ||
      any(!is.finite(v)) || any(abs(v) >= 1)) {
    stop("v must be finite with abs(v) < 1, and n must be greater than 2")
  }
  adjusted_z <- sqrt(n - 1) * v
  sign(adjusted_z) * sqrt(
    adjusted_z^2 * (n - 2) / (n - 1 - adjusted_z^2)
  )
}

# Standard SuSiE-RSS fit from the marginal association vector v and LD R.
fit_susie_rss <- function(
    v, R, N, L = 5L, prior_variance = 0.2,
    max_iter = 500L, tol = 1e-3) {
  iw_require_susieR()
  v <- as.numeric(v)
  R <- as.matrix(R)
  if (!all(dim(R) == c(length(v), length(v)))) {
    stop("R must be a square matrix matching v")
  }
  iw_susie_ss_function()(
    XtX = N * R, Xty = N * v, yty = N, n = N, L = as.integer(L),
    scaled_prior_variance = prior_variance,
    estimate_prior_variance = TRUE,
    residual_variance = 1,
    estimate_residual_variance = FALSE,
    standardize = FALSE,
    check_prior = FALSE,
    max_iter = as.integer(max_iter), tol = tol
  )
}

iw_unwrap_fit <- function(object) {
  if (inherits(object, "collapsed_r_fit")) return(object$fit)
  if (inherits(object, c("aiwn_fit", "aiw_fit"))) {
    return(object$selected_fit)
  }
  object
}

iw_selected_ld <- function(object, R = NULL) {
  if (!is.null(R)) return(as.matrix(R))
  if (!is.null(object$selected_Rtilde)) {
    return(as.matrix(object$selected_Rtilde))
  }
  if (!is.null(object$selected_Rbar)) return(as.matrix(object$selected_Rbar))
  NULL
}

iw_pip <- function(object) {
  fit <- iw_unwrap_fit(object)
  if (!is.null(fit$pip)) return(as.numeric(fit$pip))
  alpha <- as.matrix(fit$alpha)
  if (!length(alpha)) stop("The object does not contain PIPs or alpha")
  1 - apply(1 - alpha, 2L, prod)
}

iw_cs_list <- function(credible_sets) {
  if (is.null(credible_sets)) return(list())
  if (is.list(credible_sets) && "cs" %in% names(credible_sets)) {
    credible_sets <- credible_sets$cs
  }
  if (is.null(credible_sets)) return(list())
  if (!length(credible_sets)) return(list())
  unname(lapply(credible_sets, as.integer))
}

iw_reported_sets_object <- function(sets, fit, R, coverage = 0.95) {
  sets <- iw_cs_list(sets)
  if (!length(sets)) {
    return(list(
      cs = NULL, purity = NULL, cs_index = NULL, coverage = NULL,
      requested_coverage = coverage
    ))
  }
  if (is.null(names(sets)) || any(!nzchar(names(sets)))) {
    names(sets) <- paste0("L", seq_along(sets))
  }
  purity <- t(vapply(sets, function(variants) {
    if (length(variants) <= 1L) return(rep(1, 3L))
    values <- abs(R[variants, variants, drop = FALSE])
    values <- values[upper.tri(values)]
    c(min(values), mean(values), stats::median(values))
  }, numeric(3L)))
  colnames(purity) <- c("min.abs.corr", "mean.abs.corr", "median.abs.corr")
  rownames(purity) <- names(sets)

  alpha <- as.matrix(fit$alpha)
  set_coverage <- vapply(sets, function(variants) {
    if (nrow(alpha) == 1L) return(sum(alpha[1L, variants]))
    max(rowSums(alpha[, variants, drop = FALSE]))
  }, numeric(1L))
  list(
    cs = sets,
    purity = as.data.frame(purity),
    cs_index = seq_along(sets),
    coverage = set_coverage,
    requested_coverage = coverage
  )
}

# Convert an internal collapsed-R fit to the core layout returned by susieR.
# A genuine susie fit is returned unchanged apart from refreshing its PIPs.
iw_as_susie_object <- function(fit) {
  if (inherits(fit, "susie")) {
    fit$pip <- susieR::susie_get_pip(fit)
    return(fit)
  }
  if (is.null(fit$alpha) || is.null(fit$effects)) {
    stop("Cannot convert the internal fit to a susie-compatible object")
  }
  result <- list(
    alpha = as.matrix(fit$alpha),
    mu = do.call(rbind, lapply(fit$effects, `[[`, "mu")),
    mu2 = do.call(rbind, lapply(fit$effects, `[[`, "m2")),
    V = as.numeric(fit$sigma2),
    pip = as.numeric(fit$pip),
    elbo = as.numeric(fit$elbo),
    converged = isTRUE(fit$converged),
    niter = length(fit$elbo),
    sets = NULL
  )
  class(result) <- "susie"
  result
}

iw_get_credible_sets <- function(
    object, R = NULL, coverage = 0.95, min_abs_corr = 0.5) {
  iw_require_susieR()
  fit <- iw_unwrap_fit(object)
  R <- iw_selected_ld(object, R)
  if (is.null(R)) {
    stop("Supply R, or use a method object that stores selected_Rbar")
  }
  susieR::susie_get_cs(
    fit, Xcorr = R, coverage = coverage, min_abs_corr = min_abs_corr
  )
}

fit_ser_rss <- function(
    v, N, coverage = 0.95, prior_variance = 0.2,
    max_iter = 500L, tol = 1e-3) {
  fit <- fit_susie_rss(
    v = v, R = diag(length(v)), N = N, L = 1L,
    prior_variance = prior_variance, max_iter = max_iter, tol = tol
  )
  probability <- as.numeric(fit$alpha[1L, ])
  ordered <- order(probability, decreasing = TRUE)
  cutoff <- which(cumsum(probability[ordered]) >= coverage)[1L]
  list(fit = fit, cs = list(ordered[seq_len(cutoff)]))
}

iw_minimum_set_correlation <- function(variants, R) {
  variants <- as.integer(variants)
  if (length(variants) <= 1L) return(1)
  within_set <- abs(R[variants, variants, drop = FALSE])
  min(within_set[upper.tri(within_set)])
}

# Rare rescue for an impure SER set. Complete-linkage clustering is accepted
# only if every child set is pure, carries enough SER probability, and has a
# representative separated from the representatives of the other child sets.
partition_ser_set <- function(
    fit, ser_set, R, min_abs_corr = 0.5, max_sets = 3L,
    min_probability = 0.1) {
  probability <- as.numeric(fit$alpha[1L, ])
  ser_set <- unique(as.integer(ser_set))
  R <- as.matrix(R)
  original_purity <- iw_minimum_set_correlation(ser_set, R)
  no_partition <- function(reason) list(
    sets = list(ser_set), partition_used = FALSE, number_of_sets = 1L,
    original_purity = original_purity, reason = reason
  )
  if (length(ser_set) < 2L) return(no_partition("set_has_fewer_than_two_variants"))
  if (original_purity >= min_abs_corr) {
    return(no_partition("original_set_is_already_pure"))
  }

  absolute_ld <- abs(R[ser_set, ser_set, drop = FALSE])
  absolute_ld[!is.finite(absolute_ld)] <- 0
  diag(absolute_ld) <- 1
  distance <- 1 - absolute_ld
  distance[distance < 0] <- 0
  clustering <- stats::hclust(
    stats::as.dist(distance), method = "complete"
  )
  largest_candidate <- min(as.integer(max_sets), length(ser_set))
  for (number_of_sets in seq.int(2L, largest_candidate)) {
    membership <- stats::cutree(clustering, k = number_of_sets)
    candidate_sets <- unname(split(ser_set, membership))
    candidate_purity <- vapply(
      candidate_sets, iw_minimum_set_correlation, numeric(1), R = R
    )
    if (any(candidate_purity < min_abs_corr)) next
    candidate_probability <- vapply(candidate_sets, function(variants) {
      sum(probability[variants])
    }, numeric(1))
    if (any(candidate_probability < min_probability)) next
    representatives <- vapply(candidate_sets, function(variants) {
      variants[which.max(probability[variants])]
    }, integer(1))
    representative_ld <- abs(R[representatives, representatives, drop = FALSE])
    if (any(representative_ld[upper.tri(representative_ld)] >= min_abs_corr)) {
      next
    }
    names(candidate_sets) <- paste0("SER", seq_along(candidate_sets))
    return(list(
      sets = candidate_sets, partition_used = TRUE,
      number_of_sets = number_of_sets, original_purity = original_purity,
      partition_purity = candidate_purity,
      partition_probability = candidate_probability,
      representatives = representatives,
      reason = "partition_passed_all_checks"
    ))
  }
  no_partition("no_two_or_three_set_partition_passed_all_checks")
}

# Apply the manuscript rule: when a method reports at most one pure CS, refit
# SER-RSS and report its 95% set without purity filtering. Partitioning an
# impure SER set is retained as an optional exploratory postprocessing step,
# but is disabled by default because it does not change the reported union.
apply_ser_fallback <- function(
    object, v, N, R = NULL, coverage = 0.95, min_abs_corr = 0.5,
    partition_ser = FALSE,
    max_partition_sets = 3L, min_partition_probability = 0.1,
    prior_variance = 0.2, max_iter = 500L, tol = 1e-3) {
  fit <- iw_unwrap_fit(object)
  R <- iw_selected_ld(object, R)
  if (is.null(R)) stop("Supply the LD matrix R used to assess purity")
  credible_sets <- iw_get_credible_sets(
    object, R = R, coverage = coverage, min_abs_corr = min_abs_corr
  )
  raw_sets <- iw_cs_list(credible_sets)
  if (length(raw_sets) > 1L) {
    return(structure(list(
      raw_sets = raw_sets, reported_sets = raw_sets,
      fallback_used = FALSE, partition_used = FALSE,
      partition = NULL, reported_fit = fit, pip = iw_pip(fit),
      raw_credible_sets = credible_sets
    ), class = "iw_postprocessed_fit"))
  }

  ser <- fit_ser_rss(
    v = v, N = N, coverage = coverage, prior_variance = prior_variance,
    max_iter = max_iter, tol = tol
  )
  partition <- if (isTRUE(partition_ser)) {
    partition_ser_set(
      fit = ser$fit, ser_set = ser$cs[[1L]], R = R,
      min_abs_corr = min_abs_corr, max_sets = max_partition_sets,
      min_probability = min_partition_probability
    )
  } else {
    list(
      sets = ser$cs, partition_used = FALSE, number_of_sets = 1L,
      original_purity = iw_minimum_set_correlation(ser$cs[[1L]], R),
      partition_purity = NULL, partition_probability = NULL,
      representatives = NULL, reason = "partition_disabled"
    )
  }
  structure(list(
    raw_sets = raw_sets, reported_sets = partition$sets,
    fallback_used = TRUE, partition_used = partition$partition_used,
    partition = partition, reported_fit = ser$fit,
    pip = iw_pip(ser$fit), raw_credible_sets = credible_sets
  ), class = "iw_postprocessed_fit")
}

# Build the public object returned by susie_iw() and susie_aiw(). The fields
# familiar from susieR describe the final, post-processed solution. The exact
# pre-fallback method fit is retained under raw_fit, with raw_alpha and raw_pip
# supplied as convenient aliases.
iw_finalize_public_fit <- function(
    method_fit, v, N, R, method, subclass,
    coverage = 0.95, min_abs_corr = 0.5,
    ser_fallback = TRUE, partition_ser = FALSE, max_partition_sets = 3L,
    min_partition_probability = 0.1,
    prior_variance = 0.2, max_iter = 500L, tol = 1e-3,
    call = NULL, diagnostics = list()) {
  raw_fit <- iw_unwrap_fit(method_fit)
  raw_susie <- iw_as_susie_object(raw_fit)
  raw_credible_sets <- susieR::susie_get_cs(
    raw_susie, Xcorr = R, coverage = coverage,
    min_abs_corr = min_abs_corr
  )
  raw_sets <- iw_cs_list(raw_credible_sets)

  if (isTRUE(ser_fallback)) {
    postprocessed <- apply_ser_fallback(
      method_fit, v = v, N = N, R = R,
      coverage = coverage, min_abs_corr = min_abs_corr,
      partition_ser = partition_ser,
      max_partition_sets = max_partition_sets,
      min_partition_probability = min_partition_probability,
      prior_variance = prior_variance, max_iter = max_iter, tol = tol
    )
  } else {
    postprocessed <- structure(list(
      raw_sets = raw_sets,
      reported_sets = raw_sets,
      fallback_used = FALSE,
      partition_used = FALSE,
      partition = NULL,
      reported_fit = raw_fit,
      pip = as.numeric(raw_susie$pip),
      raw_credible_sets = raw_credible_sets
    ), class = "iw_postprocessed_fit")
  }

  final_fit <- iw_as_susie_object(postprocessed$reported_fit)
  final_fit$pip <- as.numeric(postprocessed$pip)
  final_fit$sets <- iw_reported_sets_object(
    postprocessed$reported_sets, final_fit, R, coverage
  )
  final_fit$method <- method
  final_fit$raw_fit <- raw_fit
  final_fit$raw_alpha <- as.matrix(raw_susie$alpha)
  final_fit$raw_pip <- as.numeric(raw_susie$pip)
  final_fit$raw_sets <- raw_sets
  final_fit$reported_sets <- postprocessed$reported_sets
  final_fit$postprocess <- postprocessed
  final_fit$iw <- diagnostics
  final_fit$call <- call
  class(final_fit) <- unique(c(subclass, "susie", class(final_fit)))
  final_fit
}

susie_get_reported_cs <- function(fit) {
  if (is.null(fit$sets) || is.null(fit$sets$cs)) return(list())
  iw_cs_list(fit$sets$cs)
}

iw_method_label <- function(object) {
  if (!is.null(object$method)) return(as.character(object$method))
  if (inherits(object, "susie")) return("SuSiE-RSS")
  class(object)[1L]
}

summarize_iw_fit <- function(object, postprocessed = NULL, causal = NULL) {
  raw_fit <- iw_unwrap_fit(object)
  sets <- if (is.null(postprocessed)) {
    if (!is.null(object$sets)) iw_cs_list(object$sets) else list()
  } else {
    postprocessed$reported_sets
  }
  result <- data.frame(
    method = iw_method_label(object),
    selected_lambda = if (!is.null(object$selected_lambda)) object$selected_lambda else NA_real_,
    raw_eta0 = if (!is.null(object$raw_eta0)) object$raw_eta0 else if (!is.null(object$selected_raw_eta0)) object$selected_raw_eta0 else NA_real_,
    eta0_used = if (!is.null(object$eta0_used)) object$eta0_used else if (!is.null(object$selected_eta0)) object$selected_eta0 else NA_real_,
    model_L = if (!is.null(object$selected_L)) {
      object$selected_L
    } else if (!is.null(object$L)) {
      object$L
    } else if (!is.null(raw_fit$V)) {
      length(raw_fit$V)
    } else {
      NA_integer_
    },
    reported_L = length(sets),
    fallback = if (is.null(postprocessed)) {
      isTRUE(object$postprocess$fallback_used)
    } else postprocessed$fallback_used,
    partitioned_SER = if (is.null(postprocessed)) {
      isTRUE(object$postprocess$partition_used)
    } else postprocessed$partition_used,
    stringsAsFactors = FALSE
  )
  if (!is.null(causal)) {
    covered <- vapply(as.integer(causal), function(j) {
      any(vapply(sets, function(set) j %in% set, logical(1)))
    }, logical(1))
    result$causal_coverage <- mean(covered)
    result$all_causal_covered <- all(covered)
    result$false_sets <- sum(!vapply(sets, function(set) {
      any(set %in% causal)
    }, logical(1)))
  }
  result
}

plot_iw_pip <- function(
    object, postprocessed = NULL, sets = NULL, causal = integer(),
    R = NULL, main = iw_method_label(object),
    xlab = "Variant index", ylab = "PIP",
    point_cex = 0.7, title_cex = 0.9,
    credible_set_cex = 1.5, credible_set_lwd = 2.5,
    legend_cex = 0.62, legend_point_cex = 1.2,
    causal_cex = 1.15,
    credible_set_colors = c(
      "dodgerblue2", "green4", "purple3", "darkorange2", "goldenrod2"
    )) {
  pip <- if (is.null(postprocessed)) iw_pip(object) else postprocessed$pip
  if (is.null(sets) && !is.null(postprocessed)) sets <- postprocessed$reported_sets
  if (is.null(sets) && !is.null(object$sets)) sets <- iw_cs_list(object$sets)
  if (is.null(sets)) sets <- list()
  R <- iw_selected_ld(object, R)
  J <- length(pip)
  graphics::plot(
    seq_len(J), pip, type = "n", ylim = c(0, 1),
    xlab = xlab, ylab = ylab, main = main,
    cex.main = title_cex, las = 1
  )
  graphics::grid(col = "grey82", lty = 3)
  graphics::points(seq_len(J), pip, pch = 16, cex = point_cex, col = "black")

  legend_labels <- character(length(sets))
  legend_colors <- character(length(sets))
  for (i in seq_along(sets)) {
    variants <- sets[[i]]
    color <- credible_set_colors[
      (i - 1L) %% length(credible_set_colors) + 1L
    ]
    graphics::points(
      variants, pip[variants], pch = 1, lwd = credible_set_lwd,
      cex = credible_set_cex,
      col = color
    )
    purity <- if (is.null(R)) NA_real_ else {
      iw_minimum_set_correlation(variants, R)
    }
    set_name <- if (!is.null(names(sets)) && nzchar(names(sets)[i])) {
      names(sets)[i]
    } else {
      paste0("L", i)
    }
    legend_labels[i] <- if (is.finite(purity)) {
      sprintf("%s: size=%d, purity=%.3f", set_name, length(variants), purity)
    } else {
      sprintf("%s: size=%d", set_name, length(variants))
    }
    legend_colors[i] <- color
  }
  if (length(causal)) {
    graphics::points(
      causal, pip[causal], pch = 16, cex = causal_cex, col = "firebrick"
    )
    graphics::rug(causal, col = "firebrick", lwd = 1.5)
  }
  if (length(sets)) {
    graphics::legend(
      "topright", legend = legend_labels, col = legend_colors,
      pch = 1, pt.cex = legend_point_cex, pt.lwd = credible_set_lwd,
      cex = legend_cex, bty = "n"
    )
  }
  invisible(list(pip = pip, sets = sets))
}

plot_iw_comparison <- function(
    results, causal = integer(), ncol = 3L,
    main = NULL, point_cex = 0.65,
    credible_set_cex = 1.5, credible_set_lwd = 2.5,
    legend_cex = 0.62, legend_point_cex = 1.2,
    causal_cex = 1.15) {
  if (is.null(names(results))) names(results) <- paste("Method", seq_along(results))
  nrow <- ceiling(length(results) / ncol)
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old), add = TRUE)
  graphics::par(mfrow = c(nrow, ncol), mar = c(4.2, 4.2, 3.1, 0.8),
                oma = c(0, 0, if (is.null(main)) 0 else 2, 0))
  for (i in seq_along(results)) {
    item <- results[[i]]
    if (is.list(item) && !is.null(item$object)) {
      plot_iw_pip(
        item$object, postprocessed = item$postprocessed,
        causal = causal, R = item$R,
        main = names(results)[i], point_cex = point_cex,
        credible_set_cex = credible_set_cex,
        credible_set_lwd = credible_set_lwd,
        legend_cex = legend_cex,
        legend_point_cex = legend_point_cex,
        causal_cex = causal_cex
      )
    } else {
      plot_iw_pip(
        item, causal = causal, main = names(results)[i], point_cex = point_cex,
        credible_set_cex = credible_set_cex,
        credible_set_lwd = credible_set_lwd,
        legend_cex = legend_cex,
        legend_point_cex = legend_point_cex,
        causal_cex = causal_cex
      )
    }
  }
  if (!is.null(main)) graphics::mtext(main, outer = TRUE, line = 0.4, cex = 1.15)
  invisible(results)
}
