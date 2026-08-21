---
editor_options: 
  markdown: 
    wrap: 72
---

# Joint recovery of lambda and eta0 in the inverse-Wishart simulation

## Main conclusion

Jointly learning $\lambda$ and $\eta_0$ is reasonable as an operational
collapsed-R fitting procedure, but the two parameters are not separately
identified when the mismatch is very large. When $\eta_0$ is small, the
ELBO generally prefers a value of $\lambda$ above the generating value
and compensates with a larger estimate of $\eta_0$. When $\eta_0$ is
moderate or large, the selected $\lambda$ moves back to the generating
value or a nearby grid point.

Thus, the experiment confirms the expected $\lambda$-$\eta_0$ trade-off.
It also shows that the trade-off is mostly a strong-mismatch phenomenon.
In the range where the LD relationship is informative, joint fitting
recovers the generating mean structure well enough for fine-mapping.

## Experimental design

This experiment preserves the easy inverse-Wishart simulation in
`local-only/analysis-archive/early-simulations/experiments-simulation-iw-collapsed-R.Rmd`:

-   $J=500$ variants, GWAS sample size $N=20{,}000$, and reference-panel
    size $N_0=2{,}000$;
-   two causal effects, fitted maximum $L=5$, and regional heritability
    $h^2=0.01$;
-   three independently simulated reference LD matrices and five
    replicates for each parameter setting, giving 15 datasets per
    setting;
-   true $\eta_0\in\{20,43,95,206,447,972,2115,4599,10000\}$;
-   $\lambda_{\mathrm{true}}=10^{-3}$ or $10^{-2}$.

For each dataset, the generating prior mean is

$$
\bar R=(1-\lambda_{\mathrm{true}})R_0+
\lambda_{\mathrm{true}}I.
$$

The in-sample LD matrix is sampled from the inverse-Wishart model used
by the original experiment and then standardized to a correlation
matrix. The method receives $R_0$, but not the generating $\bar R$.

For each candidate

$$
\lambda\in\{10^{-5},10^{-4},10^{-3},0.003,0.005,0.008,
0.01,0.015,0.02,0.025,0.03\},
$$

collapsed-R profiles $\eta_0$ by VEB over $[2,20000]$. The joint fit is
the candidate with the largest final ELBO. The fit evaluated at
$\lambda_{\mathrm{true}}$ is retained as the oracle-$\lambda$
comparison. The full experiment contains 270 datasets and 2,970
candidate-$\lambda$ fits.

## Parameter recovery

### Generating lambda = 0.001

| True eta0 | Joint median lambda | Exact lambda | Oracle median eta0 | Joint median eta0 | Oracle median log-error | Joint median log-error |
|----------:|----------:|----------:|----------:|----------:|----------:|----------:|
| 20 | 0.008 | 0% | 2.0 | 10.1 | -2.303 | -0.684 |
| 43 | 0.003 | 7% | 11.0 | 25.8 | -1.366 | -0.511 |
| 95 | 0.003 | 20% | 30.0 | 53.4 | -1.153 | -0.575 |
| 206 | 0.001 | 53% | 86.1 | 120.9 | -0.872 | -0.533 |
| 447 | 0.003 | 40% | 198.2 | 381.6 | -0.813 | -0.158 |
| 972 | 0.001 | 93% | 495.4 | 527.9 | -0.674 | -0.610 |
| 2,115 | 0.001 | 100% | 1,505.0 | 1,505.0 | -0.340 | -0.340 |
| 4,599 | 0.001 | 100% | 1,429.2 | 1,429.2 | -1.169 | -1.169 |
| 10,000 | 0.001 | 100% | 5,702.7 | 5,702.7 | -0.562 | -0.562 |

Here and below, log-error means
$\log(\widehat\eta_0)-\log(\eta_{0,\mathrm{true}})$. Joint fitting
strongly reduces the low-$\eta_0$ underestimation, but it does so by
selecting a larger $\lambda$, not by recovering the complete generating
parameter pair. At $\eta_0\geq 972$, $\lambda$ is recovered reliably and
the joint and oracle estimates of $\eta_0$ become nearly or exactly
identical.

### Generating lambda = 0.01

| True eta0 | Joint median lambda | Exact lambda | Oracle median eta0 | Joint median eta0 | Oracle median log-error | Joint median log-error |
|----------:|----------:|----------:|----------:|----------:|----------:|----------:|
| 20 | 0.015 | 13% | 2.4 | 14.2 | -2.110 | -0.342 |
| 43 | 0.015 | 33% | 23.3 | 31.7 | -0.613 | -0.305 |
| 95 | 0.015 | 47% | 55.2 | 77.3 | -0.544 | -0.206 |
| 206 | 0.010 | 53% | 152.0 | 189.6 | -0.304 | -0.083 |
| 447 | 0.010 | 60% | 355.5 | 376.3 | -0.229 | -0.172 |
| 972 | 0.010 | 53% | 883.0 | 883.0 | -0.096 | -0.096 |
| 2,115 | 0.010 | 87% | 2,389.2 | 2,882.8 | 0.122 | 0.310 |
| 4,599 | 0.010 | 80% | 2,529.5 | 3,033.9 | -0.598 | -0.416 |
| 10,000 | 0.010 | 67% | 8,679.7 | 5,968.6 | -0.142 | -0.516 |

The median selected $\lambda$ equals the generating value from
$\eta_0=206$ onward. Exact-grid recovery is lower than in the
$\lambda_{\mathrm{true}}=0.001$ experiment because the fitting grid is
much denser around 0.01: nearby values 0.008 and 0.015 often have very
similar ELBOs. Therefore, the median trend and the ELBO profile are more
informative than exact selection frequency here.

## The ELBO profiles show two regimes

For $\lambda_{\mathrm{true}}=0.001$, the median joint-versus-oracle ELBO
gain is 43.5, 20.8, and 13.5 at true $\eta_0=20,43,95$, respectively.
These are not flat-profile accidents: under strong mismatch, the
variational objective clearly prefers the compensating larger-$\lambda$,
larger-$\eta_0$ solution. By $\eta_0=972$, the median gain is zero and
the true $\lambda$ is selected in 93% of datasets.

For $\lambda_{\mathrm{true}}=0.01$, the corresponding median gains are
much smaller: 3.61, 1.31, and 0.05 at true $\eta_0=20,43,95$, and zero
from $\eta_0=206$ onward. Around 0.01, adjacent grid values commonly lie
within two ELBO units of the winner. This is a shallow local profile
rather than a large displacement of the inferred prior mean.

The ELBO profiles are shown in
`experiments/figures/simulation-iw-collapsed-R-joint-lambda-eta0-elbo-profile.png`.

## Fine-mapping behavior

Relative to fixing the generating $\lambda$:

-   with $\lambda_{\mathrm{true}}=0.001$, joint fitting increased power
    in 23 of 135 datasets and left it unchanged in the other 112;
-   with $\lambda_{\mathrm{true}}=0.01$, it increased power in 5 of 135
    datasets and left it unchanged in the other 130;
-   it did not reduce power, reported-CS coverage, or the
    purity-filtered number of credible sets in any dataset;
-   mean causal PIP increased by 0.064 for
    $\lambda_{\mathrm{true}}=0.001$ and by 0.022 for
    $\lambda_{\mathrm{true}}=0.01$;
-   the fine-mapping differences were concentrated at true
    $\eta_0\leq 447$. At $\eta_0\geq 972$, average selected $L$, power,
    and coverage were unchanged.

These results support joint fitting operationally in this easy
simulation: the method gains power in some strong-mismatch datasets and
agrees with the oracle fit once the LD relationship becomes informative.
They do not by themselves establish calibration outside this
model-generated setting.

## Numerical checks

-   All 270 selected joint fits converged.
-   No selected $\lambda$ was at either end of the candidate grid.
-   2,969 of 2,970 candidate fits converged. The only exception used
    $\lambda=10^{-5}$ and was 737 ELBO units below its dataset's winner,
    so it cannot affect the selected fit.
-   The joint fit reached an $\eta_0$ boundary in 27 of 270 datasets,
    compared with 37 of 270 oracle-$\lambda$ fits. Most boundary cases
    occurred at the smallest or largest generating $\eta_0$.

## Interpretation and recommended use

The experiment supports learning $\lambda$ and $\eta_0$ jointly, with
one important qualification:

1.  At moderate or small mismatch, $\lambda$ is recovered at or near its
    generating value, and joint fitting behaves like the oracle method.
2.  At very large mismatch, the data do not identify the exact prior
    mean well enough to separate $\lambda$ from $\eta_0$. The selected
    pair should be interpreted as an ELBO-optimal effective mismatch
    model, not as two independently calibrated physical parameters.
3.  Joint fitting reduces the severe low-$\eta_0$ underestimation but
    does not eliminate the previously observed high-$\eta_0$
    underestimation. For example, at $\lambda_{mathrm{true}}=0.001$ and
    true $\eta_0=4599$, both procedures have median estimate 1429.

There is also a simulation-specific caveat. As in the original
experiment, an inverse-Wishart covariance draw is converted to a
correlation matrix before the summary statistics are simulated.
Consequently, after this normalization,
$(\lambda_{\mathrm{true}},\eta_{0,\mathrm{true}})$ is not necessarily
the exact finite-sample optimum of the collapsed likelihood. This effect
should be strongest at small $\eta_0$, precisely where the upward shift
in $\lambda$ is observed. A small literal-model check could distinguish
this normalization effect from intrinsic non-identifiability, but it is
not needed to conclude that joint fitting is operationally useful.

## Output files

-   Analysis code:
    `experiments/simulation-iw-collapsed-R-joint-lambda-eta0.R`
-   Dataset-level oracle and joint results:
    `experiments/simulation-iw-collapsed-R-joint-lambda-eta0-results.csv`
-   Complete candidate profiles:
    `experiments/simulation-iw-collapsed-R-joint-lambda-eta0-profiles.csv`
-   Aggregated results:
    `experiments/simulation-iw-collapsed-R-joint-lambda-eta0-summary.csv`
-   Recovery and ELBO figures:
    `experiments/figures/simulation-iw-collapsed-R-joint-lambda-eta0-*.png`
