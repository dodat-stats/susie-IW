# Generic helpers for comparing fine-mapping results. SuSiE stores alpha with
# components in rows and variants in columns; the functions below expose
# reported components as a list of probability vectors.

iw_reported_components <- function(fit, credible_sets, fallback = NULL) {
  alpha <- as.matrix(fit$alpha)
  if (is.null(fallback) || !isTRUE(fallback$fallback_used)) {
    indices <- credible_sets$cs_index
    if (is.null(indices) || !length(indices)) return(list())
    return(lapply(as.integer(indices), function(index) {
      as.numeric(alpha[index, ])
    }))
  }

  ser_probability <- as.numeric(fallback$reported_fit$alpha[1L, ])
  if (!isTRUE(fallback$partition_used)) return(list(ser_probability))

  # A partitioned SER rescue is rare.  Represent each reported partition by
  # restricting the SER allocation probabilities to that partition and
  # renormalizing them.  This makes each reported alpha column sum to one.
  lapply(fallback$reported_sets, function(variants) {
    probability <- numeric(length(ser_probability))
    probability[variants] <- ser_probability[variants]
    total <- sum(probability)
    if (total <= 0) {
      probability[variants] <- 1 / length(variants)
    } else {
      probability <- probability / total
    }
    probability
  })
}

iw_component_assignment <- function(cost) {
  cost <- as.matrix(cost)
  if (!nrow(cost)) return(list(cost = 0, assignment = integer()))
  if (nrow(cost) > ncol(cost)) stop("The cost matrix must have rows <= columns.")

  search <- function(row, available) {
    if (row > nrow(cost)) return(list(cost = 0, assignment = integer()))
    candidates <- lapply(available, function(column) {
      remainder <- search(row + 1L, setdiff(available, column))
      list(
        cost = cost[row, column] + remainder$cost,
        assignment = c(column, remainder$assignment)
      )
    })
    candidates[[which.min(vapply(candidates, `[[`, numeric(1), "cost"))]]
  }
  search(1L, seq_len(ncol(cost)))
}

iw_component_distance <- function(components, reference_components) {
  L <- length(components)
  L_reference <- length(reference_components)
  if (!L && !L_reference) {
    return(list(
      distance = 0, matched_distance = 0, count_distance = 0,
      normalized_distance = 0, assignment = integer()
    ))
  }

  smaller_is_method <- L <= L_reference
  smaller <- if (smaller_is_method) components else reference_components
  larger <- if (smaller_is_method) reference_components else components
  if (!length(smaller)) {
    matched <- 0
    assignment <- integer()
  } else {
    cost <- outer(seq_along(smaller), seq_along(larger), Vectorize(
      function(i, j) sum(abs(smaller[[i]] - larger[[j]]))
    ))
    solution <- iw_component_assignment(cost)
    matched <- solution$cost
    assignment <- solution$assignment
  }
  count <- abs(L - L_reference)
  total <- matched + count
  list(
    distance = total,
    matched_distance = matched,
    count_distance = count,
    normalized_distance = total / (L + L_reference),
    assignment = assignment,
    smaller_is_method = smaller_is_method
  )
}

iw_matched_set_agreement <- function(sets, reference_sets) {
  L <- length(sets)
  L_reference <- length(reference_sets)
  if (!L && !L_reference) {
    return(list(mean_jaccard = 1, matched_jaccard_sum = 0))
  }
  smaller <- if (L <= L_reference) sets else reference_sets
  larger <- if (L <= L_reference) reference_sets else sets
  if (!length(smaller)) {
    return(list(mean_jaccard = 0, matched_jaccard_sum = 0))
  }
  jaccard <- outer(seq_along(smaller), seq_along(larger), Vectorize(
    function(i, j) {
      denominator <- length(union(smaller[[i]], larger[[j]]))
      if (!denominator) 1 else
        length(intersect(smaller[[i]], larger[[j]])) / denominator
    }
  ))
  solution <- iw_component_assignment(1 - jaccard)
  matched_sum <- sum(jaccard[cbind(seq_along(smaller), solution$assignment)])
  list(
    mean_jaccard = matched_sum / max(L, L_reference),
    matched_jaccard_sum = matched_sum
  )
}

iw_validate_component_distance <- function() {
  a <- c(0.8, 0.2, 0)
  b <- c(0, 0.3, 0.7)
  identical_permuted <- iw_component_distance(list(a, b), list(b, a))
  one_unmatched <- iw_component_distance(list(a), list(a, b))
  disjoint <- iw_component_distance(list(c(1, 0)), list(c(0, 1)))
  reverse <- iw_component_distance(list(a, b), list(a))
  stopifnot(
    identical_permuted$distance < 1e-12,
    abs(one_unmatched$distance - 1) < 1e-12,
    abs(disjoint$distance - 2) < 1e-12,
    abs(reverse$distance - one_unmatched$distance) < 1e-12,
    abs(one_unmatched$normalized_distance - 1 / 3) < 1e-12
  )
  invisible(TRUE)
}
