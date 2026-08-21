# Joint lambda and eta0 experiment with 1.537 calibration

## Main result

The requested experiment supports using the fitted constant correction

$$
\widetilde\eta_0=\min\{20000,1.537\widehat\eta_0\}.
$$

Relative to the uncorrected profiled estimate, the correction raises mean
reported L from 1.544 to 1.570 and mean causal power from 0.757 to 0.770.
Mean credible-set coverage is unchanged at 0.978. In paired dataset-level
comparisons, the correction increases reported L and power in 7 of 270
datasets, never decreases either one, and never decreases coverage.

The benefit is concentrated in the intended moderate-mismatch regime. At true
eta0 = 206, mean reported L rises from 1.50 to 1.60 and power rises from 0.733
to 0.783, while coverage remains 0.983. At true eta0 = 447, power rises from
0.933 to 0.950 with coverage fixed at 1.000. Under severe mismatch, the method
still uses the conservative SER-RSS fallback; under weak mismatch, the raw and
corrected fits are essentially identical.

## Experimental design

The experiment uses the controlled inverse-Wishart data-generating process:

- J = 500 variants, GWAS sample size N = 20,000, and reference-panel size
  N0 = 2,000;
- two causal effects, fitted maximum L = 5, and regional heritability h2 =
  0.01;
- three reference LD matrices and five phenotype replicates, giving 15
  datasets per parameter setting;
- true eta0 in `{20, 43, 95, 206, 447, 972, 2115, 4599, 10000}`;
- true lambda equal to 0.001 or 0.01.

For each dataset,

$$
\bar R=(1-\lambda_{\mathrm{true}})R_0+\lambda_{\mathrm{true}}I,
$$

and

$$
R=\operatorname{cov2cor}\left[
  \operatorname{IW}\left(
    \eta_{0,\mathrm{true}}+J+1,
    \eta_{0,\mathrm{true}}\bar R
  \right)
\right].
$$

The fitted lambda grid is

$$
\lambda\in\{10^{-5},10^{-4},10^{-3},0.002,0.004,0.006,0.008,0.01\}.
$$

At every candidate lambda, collapsed-R first profiles eta0 by variational
empirical Bayes. The calibrated method then refits at
$1.537\widehat\eta_0$, capped at 20,000, and selects lambda and L using that
refit's ELBO. The comparison contains only four reported methods:

1. SuSiE-RSS with the in-sample LD R;
2. SuSiE-RSS with the out-of-sample LD R0;
3. collapsed-R with jointly selected lambda and profiled eta0;
4. collapsed-R with jointly selected lambda and calibrated
   $1.537\widehat\eta_0$.

No oracle-lambda method is included.

## Postprocessing

For each proposed collapsed-R fit, raw L is the number of deduplicated 95%
credible sets passing the minimum absolute-correlation purity threshold of
0.5. The project's postprocessing rule is then applied:

1. if raw L is zero or one, refit SER-RSS using only the unit diagonal and
   report its single 95% credible set without an LD-purity filter;
2. if raw L is larger than one, report the collapsed-R credible sets.

The SuSiE-RSS benchmarks are reported directly and do not receive this
fallback. The output records both raw and postprocessed results and a fallback
indicator.

## Comparison with SuSiE-RSS

The entries below are mean reported L / credible-set coverage / causal power,
averaged over the two true lambda values, three reference panels, and five
phenotype replicates.

| True eta0 | In-sample SuSiE-RSS | Out-of-sample SuSiE-RSS | Collapsed-R: raw eta0 | Collapsed-R: 1.537x eta0 |
|---:|:---:|:---:|:---:|:---:|
| 20 | 1.77 / 1.000 / 0.883 | 4.63 / 0.240 / 0.550 | 1.00 / 0.933 / 0.467 | 1.00 / 0.933 / 0.467 |
| 43 | 1.97 / 1.000 / 0.983 | 4.17 / 0.358 / 0.733 | 1.00 / 1.000 / 0.500 | 1.03 / 1.000 / 0.517 |
| 95 | 1.97 / 1.000 / 0.983 | 3.23 / 0.540 / 0.817 | 1.13 / 0.933 / 0.533 | 1.20 / 0.933 / 0.567 |
| 206 | 1.93 / 0.933 / 0.900 | 2.43 / 0.768 / 0.833 | 1.50 / 0.983 / 0.733 | 1.60 / 0.983 / 0.783 |
| 447 | 2.00 / 1.000 / 1.000 | 2.00 / 1.000 / 1.000 | 1.87 / 1.000 / 0.933 | 1.90 / 1.000 / 0.950 |
| 972 | 1.90 / 0.989 / 0.933 | 1.93 / 0.978 / 0.933 | 1.87 / 1.000 / 0.933 | 1.87 / 1.000 / 0.933 |
| 2,115 | 1.90 / 0.967 / 0.917 | 1.83 / 0.967 / 0.883 | 1.80 / 0.983 / 0.883 | 1.80 / 0.983 / 0.883 |
| 4,599 | 1.93 / 0.983 / 0.950 | 1.93 / 0.967 / 0.933 | 1.90 / 0.967 / 0.917 | 1.90 / 0.967 / 0.917 |
| 10,000 | 1.80 / 1.000 / 0.900 | 1.83 / 1.000 / 0.917 | 1.83 / 1.000 / 0.917 | 1.83 / 1.000 / 0.917 |

Across all 270 datasets:

| Method | Mean reported L | Coverage | Power | Mean causal PIP | Fallback frequency |
|:---|---:|---:|---:|---:|---:|
| SuSiE-RSS: in-sample R | 1.907 | 0.986 | 0.939 | 0.744 | 0 |
| SuSiE-RSS: out-of-sample R0 | 2.667 | 0.757 | 0.844 | 0.672 | 0 |
| Collapsed-R: profiled eta0 | 1.544 | 0.978 | 0.757 | 0.574 | 0.456 |
| Collapsed-R: 1.537x eta0 | 1.570 | 0.978 | 0.770 | 0.584 | 0.430 |

This makes the main trade-off especially clear. Out-of-sample SuSiE-RSS has
higher average power than collapsed-R, but at low eta0 it strongly
overspecifies L and its credible-set coverage is very poor. The proposed
method sacrifices power under severe mismatch, where it deliberately reports
the SER fallback, but retains high coverage. From true eta0 around 447 onward,
all methods behave similarly and the proposed fit approaches the in-sample
benchmark.

## Eta0 recovery

The following table summarizes the jointly selected calibrated fit, averaged
over both generating lambda values. Eta0 columns are medians; performance
columns are means.

| True eta0 | Raw eta0 estimate | Eta0 used | Reported L | Coverage | Power | Fallback frequency |
|---:|---:|---:|---:|---:|---:|---:|
| 20 | 3.8 | 5.8 | 1.00 | 0.933 | 0.467 | 1.00 |
| 43 | 25.5 | 39.1 | 1.03 | 1.000 | 0.517 | 0.97 |
| 95 | 57.6 | 88.5 | 1.20 | 0.933 | 0.567 | 0.80 |
| 206 | 147.6 | 226.8 | 1.60 | 0.983 | 0.783 | 0.40 |
| 447 | 353.9 | 544.0 | 1.90 | 1.000 | 0.950 | 0.10 |
| 972 | 820.2 | 1,260.6 | 1.87 | 1.000 | 0.933 | 0.13 |
| 2,115 | 2,253.9 | 3,464.3 | 1.80 | 0.983 | 0.883 | 0.20 |
| 4,599 | 1,999.3 | 3,072.9 | 1.90 | 0.967 | 0.917 | 0.10 |
| 10,000 | 5,247.4 | 8,065.3 | 1.83 | 1.000 | 0.917 | 0.17 |

The 1.537 multiplier corrects much of the underestimation through the middle
of the range, although it cannot fully remove the nonlinear recovery pattern.
In particular, recovery remains noisy at true eta0 = 4,599 and 10,000. This is
not consequential for fine-mapping here because the fitted results have
already reached their weak-mismatch behavior.

For reference, the previous two-fold correction had mean L = 1.596, coverage
= 0.976, and power = 0.780. Thus 1.537x is slightly more conservative than 2x:
it gives up about 0.009 mean power and gains about 0.002 mean coverage. Unlike
2x, however, 1.537 is the multiplier directly estimated by the fixed-slope
calibration model.

## Numerical checks

- All 1,080 selected fits converged.
- Of 4,320 candidate profiles, 4,319 converged. The one nonconverged unused
  candidate was never selected.
- All proposed postprocessed results have reported L at least one.
- Fallback is used exactly when a proposed method's raw purity-filtered L is
  zero or one: 123 times for raw profiled eta0 and 116 times for calibrated
  eta0.
- Fallback is never applied to either SuSiE-RSS benchmark.
- The result files contain 4,320 candidate profiles, 1,080 selected fits,
  1,080 postprocessed fits, and 72 aggregate rows.

## Output files

- Analysis code:
  `experiments/simulation-iw-collapsed-R-joint-lambda-eta0.R`
- Raw selected fits:
  `experiments/simulation-iw-collapsed-R-joint-lambda-eta0-capped-1.537-results.csv`
- Postprocessed selected fits:
  `experiments/simulation-iw-collapsed-R-joint-lambda-eta0-capped-1.537-postprocessed-results.csv`
- Postprocessed summary used by the figures:
  `experiments/simulation-iw-collapsed-R-joint-lambda-eta0-capped-1.537-summary.csv`
- Candidate lambda profiles:
  `experiments/simulation-iw-collapsed-R-joint-lambda-eta0-capped-1.537-profiles.csv`
- Eta0 calibration coefficients:
  `experiments/simulation-iw-collapsed-R-joint-lambda-eta0-capped-1.537-eta0-calibration.csv`
- Figures:
  `experiments/figures/simulation-iw-collapsed-R-joint-lambda-eta0-capped-1.537-*.{png,pdf}`
