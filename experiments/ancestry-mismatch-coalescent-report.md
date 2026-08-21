# Out-of-model ancestry-mismatch coalescent experiment

Two empirical populations split T generations ago. The experiment uses J=500, 8 independent panels, N0 in {100, 500, 2000}, T in {1, 50, 200, 1000, 5000}, true L=2, 7:5 effect weights, h2=0.001 and GWAS N=200,000.

The inverse-Wishart relation is not used to generate R from R0; both are empirical LD matrices from the two populations.

At N0=100, reference-LD SuSiE-RSS reported 4.25--4.75 sets across the divergence grid, with credible-set coverage 0.342--0.395. Both proposed methods remained in the one-set fallback regime, with coverage one and power 0.5. Thus finite-panel noise dominates this row of the design.

At N0=2,000 and T<=200, reference-LD SuSiE-RSS and AIW-N both reported two sets with coverage and power one. At T=1,000, reference SuSiE-RSS reported 2.375 sets with coverage 0.842; AIW-N reported 1.750 sets with coverage 1.000 and power 0.875. At T=5,000, reference SuSiE-RSS reported 3.375 sets with coverage 0.556, whereas both proposed methods reported one set with coverage one and power 0.5.

Off-diagonal LD RMSE increased with T at every N0. The fitted eta0 values decreased along the same gradient, showing that both working models detect structured ancestry mismatch even though the inverse-Wishart model did not generate either empirical LD matrix.

All 400 distinct model fits converged across 120 reference configurations.

Generated files:

- `experiments/ancestry-mismatch-coalescent-results.csv`
- `experiments/ancestry-mismatch-coalescent-summary.csv`
- `experiments/figures/ancestry-mismatch-performance.{pdf,png}`
- `experiments/figures/ancestry-mismatch-selected-L.{pdf,png}`
- `experiments/figures/ancestry-mismatch-ld-error.{pdf,png}`
- `experiments/figures/ancestry-mismatch-estimated-eta0.{pdf,png}`
