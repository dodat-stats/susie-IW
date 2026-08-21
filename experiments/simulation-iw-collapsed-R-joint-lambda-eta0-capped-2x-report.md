# Capped lambda, SER fallback, and eta0 calibration

## Main conclusions

The earlier low coverage at true eta0 = 20 was caused by reporting raw
purity-filtered credible sets instead of applying the project's SER-RSS
fallback. After applying the fallback whenever raw purity-filtered L is zero or
one, both generating-lambda settings have mean reported L = 1, coverage =
0.933, and power = 0.467 at true eta0 = 20.

The two-fold eta0 correction modestly improves power at intermediate eta0
without causing overspecification. Across all 270 paired joint fits, mean power
increases from 0.757 to 0.780, while coverage is essentially unchanged (0.978
versus 0.976). Neither method ever reports more than the true L = 2.

The requested calibration regression is

$$
\log(\eta_{0,\mathrm{true}})
=1.5888+0.7957\log(\widehat\eta_0),
$$

with $R^2=0.873$. Therefore,

$$
\exp(\text{intercept})=4.898.
$$

Because the fitted slope is 0.796 rather than one, 4.898 is not a single
multiplicative correction. The fitted calibration curve is

$$
\widehat\eta_{0,\mathrm{calibrated}}
=4.898\widehat\eta_0^{0.7957}.
$$

If a single constant multiplier is required, the regression should instead
fix the slope at one. That estimate is 1.537, not 4.898.

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

The fitting grid is capped at

$$
\lambda\in\{10^{-5},10^{-4},10^{-3},0.002,0.004,0.006,0.008,0.01\}.
$$

At every candidate lambda, collapsed-R profiles eta0 by variational empirical
Bayes. For the two-fold comparison, the model is refitted at

$$
\widetilde\eta_0(\lambda)
=\min\{20000,2\widehat\eta_0(\lambda)\},
$$

and lambda is selected using the corrected-fit ELBO.

## Reporting rule

Let raw L denote the number of deduplicated 95% credible sets passing the
minimum absolute-correlation purity threshold of 0.5. The reported result is:

1. if raw L is zero or one, fit SER-RSS using only the unit diagonal and report
   its single 95% credible set without an LD-purity filter;
2. if raw L is larger than one, report the collapsed-R credible sets.

The output now stores raw L and reported L separately, together with the
fallback indicator. Thus, raw L may be zero, but reported L is always at least
one in these associated-locus simulations.

## Results after SER fallback

The raw eta0 column below is the profiled estimate at the lambda selected by
the corrected method. Each row summarizes 15 datasets.

| True lambda | True eta0 | Lambda: raw to corrected | Raw eta0 | Eta0 used | Reported L: raw to corrected | Fallback frequency: raw to corrected | Coverage: raw to corrected | Power: raw to corrected |
|---:|---:|:---:|---:|---:|:---:|:---:|:---:|:---:|
| 0.001 | 20 | 0.006 to 0.006 | 3.9 | 7.9 | 1.00 to 1.00 | 1.00 to 1.00 | 0.93 to 0.93 | 0.47 to 0.47 |
| 0.001 | 43 | 0.004 to 0.004 | 29.2 | 58.4 | 1.00 to 1.07 | 1.00 to 0.93 | 1.00 to 1.00 | 0.50 to 0.53 |
| 0.001 | 95 | 0.004 to 0.002 | 60.0 | 119.9 | 1.13 to 1.40 | 0.87 to 0.60 | 0.93 to 1.00 | 0.53 to 0.70 |
| 0.001 | 206 | 0.002 to 0.002 | 120.9 | 241.9 | 1.53 to 1.67 | 0.47 to 0.33 | 0.97 to 0.93 | 0.73 to 0.77 |
| 0.001 | 447 | 0.002 to 0.002 | 374.7 | 749.5 | 1.87 to 1.93 | 0.13 to 0.07 | 1.00 to 1.00 | 0.93 to 0.97 |
| 0.001 | 972 | 0.001 to 0.002 | 878.6 | 1757.1 | 1.87 to 1.87 | 0.13 to 0.13 | 1.00 to 1.00 | 0.93 to 0.93 |
| 0.001 | 2,115 | 0.001 to 0.001 | 1963.1 | 3926.3 | 1.80 to 1.80 | 0.20 to 0.20 | 1.00 to 0.97 | 0.90 to 0.87 |
| 0.001 | 4,599 | 0.001 to 0.001 | 1429.2 | 2858.5 | 1.87 to 1.93 | 0.13 to 0.07 | 0.97 to 0.97 | 0.90 to 0.93 |
| 0.001 | 10,000 | 0.001 to 0.001 | 5702.7 | 11405.4 | 1.87 to 1.87 | 0.13 to 0.13 | 1.00 to 1.00 | 0.93 to 0.93 |
| 0.01 | 20 | 0.01 to 0.01 | 2.4 | 4.8 | 1.00 to 1.00 | 1.00 to 1.00 | 0.93 to 0.93 | 0.47 to 0.47 |
| 0.01 | 43 | 0.01 to 0.01 | 23.3 | 46.6 | 1.00 to 1.07 | 1.00 to 0.93 | 1.00 to 1.00 | 0.50 to 0.53 |
| 0.01 | 95 | 0.01 to 0.01 | 55.0 | 110.1 | 1.13 to 1.20 | 0.87 to 0.80 | 0.93 to 0.93 | 0.53 to 0.57 |
| 0.01 | 206 | 0.01 to 0.01 | 143.2 | 286.3 | 1.47 to 1.60 | 0.53 to 0.40 | 1.00 to 0.97 | 0.73 to 0.77 |
| 0.01 | 447 | 0.01 to 0.01 | 355.5 | 711.1 | 1.87 to 1.93 | 0.13 to 0.07 | 1.00 to 1.00 | 0.93 to 0.97 |
| 0.01 | 972 | 0.01 to 0.01 | 860.9 | 1721.8 | 1.87 to 1.87 | 0.13 to 0.13 | 1.00 to 1.00 | 0.93 to 0.93 |
| 0.01 | 2,115 | 0.01 to 0.01 | 2389.2 | 4778.5 | 1.80 to 1.80 | 0.20 to 0.20 | 0.97 to 0.97 | 0.87 to 0.87 |
| 0.01 | 4,599 | 0.01 to 0.01 | 2529.5 | 5059.0 | 1.93 to 1.93 | 0.07 to 0.07 | 0.97 to 0.97 | 0.93 to 0.93 |
| 0.01 | 10,000 | 0.01 to 0.01 | 4128.0 | 8256.1 | 1.80 to 1.80 | 0.20 to 0.20 | 1.00 to 1.00 | 0.90 to 0.90 |

Across all paired joint fits:

| Metric | Raw profiled eta0 | Two-fold eta0 |
|:---|---:|---:|
| Mean reported L | 1.544 | 1.596 |
| Mean coverage | 0.978 | 0.976 |
| Mean power | 0.757 | 0.780 |
| Mean causal PIP | 0.574 | 0.593 |
| Fallback frequency | 0.456 | 0.404 |

The correction increases reported L in 14 datasets and never decreases it.
Power increases in 14 datasets, decreases in 3, and is unchanged in 253.
Coverage increases once, decreases three times, and is unchanged in 266. These
small changes are consistent with essentially preserved coverage and a modest
power gain.

## Eta0 calibration

The calibration regression uses all 270 uncorrected jointly selected fits:

| Stratum | Fits | Intercept | exp(intercept) | Slope | R-squared |
|:---|---:|---:|---:|---:|---:|
| Pooled | 270 | 1.5888 | 4.898 | 0.7957 | 0.873 |
| True lambda = 0.001 | 135 | 1.5968 | 4.937 | 0.7941 | 0.861 |
| True lambda = 0.01 | 135 | 1.5810 | 4.860 | 0.7972 | 0.884 |

The nearly identical stratified coefficients suggest that the calibration
curve is stable across these two generating lambda values. However, the slope
is clearly below one. Rewriting the pooled curve as a correction factor gives

$$
\frac{\widehat\eta_{0,\mathrm{calibrated}}}{\widehat\eta_0}
=4.898\widehat\eta_0^{-0.2043}.
$$

Thus, the implied multiplier decreases with estimated eta0. It is approximately
3.69 at estimated eta0 = 4, 2.54 at 25, 1.91 at 100, and 1.19 at 1000. A
fixed two-fold correction is therefore a rough approximation around the middle
of the useful range, rather than a globally estimated constant.

For comparison, forcing the slope to one estimates

$$
\log(\eta_{0,\mathrm{true}})
-\log(\widehat\eta_0)=0.4299,
$$

so the fitted single multiplier is

$$
\exp(0.4299)=1.537.
$$

Before replacing the two-fold correction, the natural next comparison is
therefore between the nonlinear calibration
$4.898\widehat\eta_0^{0.7957}$ and the genuine constant-multiplier correction
$1.537\widehat\eta_0$, evaluated on independently simulated datasets.

## Numerical checks

- All 1,080 selected oracle and joint fits converged.
- All postprocessed results have reported L at least one.
- Fallback is used exactly when raw purity-filtered L is zero or one.
- At true eta0 = 20, fallback frequency is one for all four oracle/joint and
  raw/corrected comparisons.
- The result files contain 4,320 candidate profiles, 1,080 raw selected/oracle
  fits, 1,080 postprocessed fits, and 72 aggregate rows.

## Output files

- Analysis code:
  `experiments/simulation-iw-collapsed-R-joint-lambda-eta0.R`
- Raw dataset-level selected fits:
  `experiments/simulation-iw-collapsed-R-joint-lambda-eta0-capped-2x-results.csv`
- Postprocessed dataset-level fits:
  `experiments/simulation-iw-collapsed-R-joint-lambda-eta0-capped-2x-postprocessed-results.csv`
- Raw summary:
  `experiments/simulation-iw-collapsed-R-joint-lambda-eta0-capped-2x-raw-summary.csv`
- Postprocessed summary used by the figures:
  `experiments/simulation-iw-collapsed-R-joint-lambda-eta0-capped-2x-summary.csv`
- Eta0 calibration coefficients:
  `experiments/simulation-iw-collapsed-R-joint-lambda-eta0-capped-2x-eta0-calibration.csv`
- Figures:
  `experiments/figures/simulation-iw-collapsed-R-joint-lambda-eta0-capped-2x-*.{png,pdf}`
