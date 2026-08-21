# Shared settings and simulation helpers for the out-of-model ancestry-
# mismatch experiment. GWAS and reference genotypes are sampled from two
# populations that split T generations ago. Their empirical LD matrices are
# used directly; no inverse-Wishart relationship is imposed between them.

ancestry_get_integer <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.integer(default))
  as.integer(value)
}

ancestry_J <- ancestry_get_integer("ANCESTRY_J", 500L)
ancestry_panel_seeds <- seq.int(
  ancestry_get_integer("ANCESTRY_FIRST_SEED", 101L),
  length.out = ancestry_get_integer("ANCESTRY_REPLICATIONS", 8L)
)
ancestry_T_grid <- c(1L, 50L, 200L, 1000L, 5000L)
ancestry_N0_grid <- c(100L, 500L, 2000L)
ancestry_gwas_ld_sample_size <- 2000L
ancestry_N <- 200000L
ancestry_h2 <- 1e-3
ancestry_true_L <- 2L
ancestry_effect_weights <- c(7, 5)
ancestry_maximum_L <- 5L
ancestry_lambda_grid <- c(
  1e-5, 1e-4, 1e-3, 2e-3, 4e-3, 6e-3, 8e-3, 1e-2
)
ancestry_collapsed_multiplier <- 1.677
ancestry_aiwn_multiplier <- 1.424
ancestry_coverage <- 0.95
ancestry_purity <- 0.5
ancestry_worker_directory <- file.path(
  "analysis", "ancestry-mismatch-coalescent-workers"
)
ancestry_output_stem <- file.path(
  "analysis", "ancestry-mismatch-coalescent"
)

ancestry_load_python <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Install the R package reticulate")
  }
  if (!nzchar(Sys.getenv("RETICULATE_PYTHON"))) {
    python <- Sys.which("python3")
    if (nzchar(python)) {
      Sys.setenv(RETICULATE_PYTHON = python)
      reticulate::use_python(python, required = FALSE)
    }
  }
  list(
    stdpopsim = reticulate::import("stdpopsim"),
    msprime = reticulate::import("msprime")
  )
}

ancestry_to_diploid <- function(genotype_haploid, sample_ids) {
  columns <- as.integer(sample_ids) + 1L
  haplotypes <- genotype_haploid[, columns, drop = FALSE]
  genotypes <- haplotypes[, seq(1L, ncol(haplotypes), by = 2L)] +
    haplotypes[, seq(2L, ncol(haplotypes), by = 2L)]
  t(genotypes)
}

ancestry_ld <- function(genotypes) {
  R <- stats::cor(genotypes)
  if (any(!is.finite(R))) stop("Non-finite empirical LD matrix")
  R <- (R + t(R)) / 2
  diag(R) <- 1
  R
}

simulate_two_population_ld <- function(
    J, T_generations, panel_seed,
    gwas_sample_size = ancestry_gwas_ld_sample_size,
    reference_sample_sizes = ancestry_N0_grid,
    chromosome = "chr22", start = 20e6, end = 21e6) {
  python <- ancestry_load_python()
  stdpopsim <- python$stdpopsim
  msprime <- python$msprime
  maximum_reference_size <- max(reference_sample_sizes)

  species <- stdpopsim$get_species("HomSap")
  contig <- species$get_contig(
    chromosome, left = as.integer(start), right = as.integer(end)
  )
  demography <- msprime$Demography()
  demography$add_population(name = "ANC", initial_size = 20000L)
  demography$add_population(name = "GWAS", initial_size = 100000L)
  demography$add_population(name = "REF", initial_size = 100000L)
  demography$add_population_split(
    time = as.integer(T_generations),
    derived = list("GWAS", "REF"), ancestral = "ANC"
  )
  models <- list(
    msprime$DiscreteTimeWrightFisher(duration = 100L),
    msprime$StandardCoalescent()
  )
  ancestry_seed <- as.integer(
    100000L + 10000L * match(T_generations, ancestry_T_grid) + panel_seed
  )
  mutation_seed <- as.integer(200000L + ancestry_seed)
  tree_sequence <- msprime$sim_ancestry(
    samples = reticulate::dict(
      GWAS = as.integer(gwas_sample_size),
      REF = as.integer(maximum_reference_size)
    ),
    demography = demography,
    recombination_rate = contig$recombination_map,
    model = models,
    random_seed = ancestry_seed
  )
  tree_sequence <- msprime$sim_mutations(
    tree_sequence, rate = 1.29e-8, random_seed = mutation_seed
  )

  population_ids <- list()
  for (population in reticulate::iterate(tree_sequence$populations())) {
    population_ids[[population$metadata$name]] <- population$id
  }
  genotype_haploid <- tree_sequence$genotype_matrix()
  genotype_gwas <- ancestry_to_diploid(
    genotype_haploid,
    tree_sequence$samples(population = as.integer(population_ids$GWAS))
  )
  genotype_reference <- ancestry_to_diploid(
    genotype_haploid,
    tree_sequence$samples(population = as.integer(population_ids$REF))
  )

  set.seed(300000L + ancestry_seed)
  reference_order <- sample.int(nrow(genotype_reference))
  genotype_reference <- genotype_reference[reference_order, , drop = FALSE]
  smallest_reference <- genotype_reference[
    seq_len(min(reference_sample_sizes)), , drop = FALSE
  ]
  frequency_gwas <- colMeans(genotype_gwas) / 2
  frequency_reference <- colMeans(genotype_reference) / 2
  eligible <- frequency_gwas > 0.05 & frequency_gwas < 0.95 &
    frequency_reference > 0.05 & frequency_reference < 0.95 &
    apply(genotype_gwas, 2L, stats::sd) > 0 &
    apply(smallest_reference, 2L, stats::sd) > 0
  eligible_columns <- which(eligible)
  if (length(eligible_columns) < J) {
    stop(sprintf(
      "Only %d eligible variants for J=%d (T=%d, seed=%d)",
      length(eligible_columns), J, T_generations, panel_seed
    ))
  }

  # Select a seeded contiguous block among variants that pass common-SNP QC in
  # both populations and remain polymorphic in the smallest reference panel.
  set.seed(400000L + ancestry_seed)
  first <- if (length(eligible_columns) == J) {
    1L
  } else {
    sample.int(length(eligible_columns) - J + 1L, 1L)
  }
  columns <- eligible_columns[first + seq_len(J) - 1L]
  genotype_gwas <- genotype_gwas[, columns, drop = FALSE]
  genotype_reference <- genotype_reference[, columns, drop = FALSE]
  R <- ancestry_ld(genotype_gwas)
  R0 <- lapply(reference_sample_sizes, function(N0) {
    ancestry_ld(genotype_reference[seq_len(N0), , drop = FALSE])
  })
  names(R0) <- as.character(reference_sample_sizes)

  list(
    R = R, R0 = R0,
    gwas_frequencies = colMeans(genotype_gwas) / 2,
    reference_frequencies = lapply(reference_sample_sizes, function(N0) {
      colMeans(genotype_reference[seq_len(N0), , drop = FALSE]) / 2
    }),
    ancestry_seed = ancestry_seed,
    mutation_seed = mutation_seed,
    selected_columns = columns
  )
}

choose_ancestry_weak_ld_variants <- function(R) {
  candidates <- unique(as.integer(round(seq(
    max(2, 0.08 * nrow(R)), min(nrow(R) - 1, 0.92 * nrow(R)),
    length.out = min(250L, nrow(R) - 2L)
  ))))
  minimum_separation <- max(8L, floor(nrow(R) / 12))
  first <- candidates[which.min(abs(candidates - round(nrow(R) / 4)))]
  available <- candidates[abs(candidates - first) >= minimum_separation]
  second <- available[which.min(abs(R[first, available]))]
  c(first, second)
}

simulate_ancestry_summary_statistics <- function(
    R, T_generations, panel_seed, N = ancestry_N, h2 = ancestry_h2) {
  causal <- choose_ancestry_weak_ld_variants(R)
  scale <- sqrt(
    h2 / drop(crossprod(
      ancestry_effect_weights,
      R[causal, causal, drop = FALSE] %*% ancestry_effect_weights
    ))
  )
  beta <- numeric(nrow(R))
  beta[causal] <- scale * ancestry_effect_weights
  set.seed(
    500000L + 10000L * match(T_generations, ancestry_T_grid) + panel_seed
  )
  eigendecomposition <- eigen(R, symmetric = TRUE)
  noise <- as.numeric(
    eigendecomposition$vectors %*%
      (sqrt(pmax(eigendecomposition$values, 0)) * stats::rnorm(nrow(R)))
  ) / sqrt(N)
  v <- as.numeric(R %*% beta + noise)
  list(
    v = v, beta = beta, causal = causal,
    causal_LD = R[causal[1L], causal[2L]],
    expected_z = sqrt(N) * as.numeric(R[causal, , drop = FALSE] %*% beta)
  )
}

ancestry_ld_metrics <- function(R, R0) {
  upper <- upper.tri(R)
  difference <- R0[upper] - R[upper]
  data.frame(
    ld_rmse = sqrt(mean(difference^2)),
    ld_mae = mean(abs(difference)),
    ld_correlation = stats::cor(R[upper], R0[upper]),
    maximum_ld_error = max(abs(difference))
  )
}
