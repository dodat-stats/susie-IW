# Coalescent SuSiE-IW benchmark

## Design

This independent evaluation contains 480 datasets and 1920 method fits. Twenty independent coalescent panels are distributed across J in {500, 1000, 2000}, with 7, 7 and 6 panels, respectively. The experiment crosses true eta0 in {20, 50, 100, 200, 500, 1000, 2000, 5000} with true L in {1, 2, 3}, using h2=1e-3, N=200,000, N0=2,000 and true lambda=1e-3. Nested weak-LD causal-effect weights are 7, 7:5 and 7:5:3. The SER fallback reports its original 95% set without purity filtering or partitioning.

The evaluation-panel seeds (61 through 80) are disjoint from the calibration-panel seeds (51 and 52). Collapsed-R uses eta0 multiplier 1.677 and AIW-N uses 1.424; these constants are frozen before performance evaluation. Error bars use a nonparametric bootstrap over the 20 independent (J, R0) panels.

## Main findings

- At true eta0=20, reference-LD SuSiE-RSS reports means of 4.5, 4.5, 4.6 sets when true L is 1, 2 and 3. Collapsed-R reports 1.0, 1.0, 1.0, and AIW-N reports 1.0, 1.0, 1.0.
- With two true effects and true eta0=20, reference-LD SuSiE-RSS has causal power 0.62, credible-set coverage 0.278 and 3.2 false sets per dataset. Collapsed-R and AIW-N each have power 0.47, coverage 0.950 and 0.1 false sets. This is the intended conservative fallback regime.
- With two true effects and true eta0=200, AIW-N reports 2.0 sets with power 1.00 and coverage 1.00. Collapsed-R reports 1.9 sets, with power 0.95 and coverage 1.00. From true eta0=500 onward, all four methods report two sets.
- Pooled over the 24 design cells, reference-LD SuSiE-RSS reports 2.515 sets and 0.892 false sets per dataset, with coverage 0.770. Collapsed-R reports 1.485 sets and 0.029 false sets, with coverage 0.973; AIW-N reports 1.587 sets and 0.044 false sets, with coverage 0.970.
- All 1920 fits converged. Median fit times are 1.48 seconds for collapsed-R, 3.44 seconds for AIW-N, 0.040 seconds for in-sample SuSiE-RSS and 0.049 seconds for reference-LD SuSiE-RSS.

The key qualitative result is the crossover: moving from severe mismatch toward matched LD, the mismatch-aware methods stop withholding additional signals at approximately the same point that reference-LD SuSiE-RSS stops over-specifying the number of signals.

## Outputs

- `coalescent-susie-iw-benchmark-main-table.csv`: selected manuscript cells.
- `coalescent-susie-iw-benchmark-summary.csv`: all method/design summaries and bootstrap intervals.
- `coalescent-susie-iw-benchmark-average-over-L-summary.csv`: main-text power and coverage averaged over L.
- `coalescent-susie-iw-benchmark-results.csv`: all configuration-level fits.
- `figures/coalescent-susie-iw-selected-L.pdf`
- `figures/coalescent-susie-iw-performance.pdf`
- `figures/coalescent-susie-iw-performance-L1.pdf`
- `figures/coalescent-susie-iw-performance-L2.pdf`
- `figures/coalescent-susie-iw-performance-L3.pdf`
- `figures/coalescent-susie-iw-diagnostics.pdf`
- `figures/coalescent-susie-iw-motivating-locus.pdf`
