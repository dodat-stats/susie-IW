# Method implementations

This directory contains the maintained R implementation used by the
manuscript and workflowr analyses.

- `01-collapsed-R.R`: collapsed-R inference for the SuSiE-IW model. The main
  public entry point is `susie_iw()`.
- `02-AIW-and-AIW-N.R`: AIW and AIW-N working approximations. The main public
  entry point is `susie_aiw()` with `approximation = "t"` or `"N"`.
- `03-SER-fallback-and-plot.R`: shared SuSiE-RSS fitting, operational SER-RSS
  fallback, credible-set extraction, summaries, and PIP plotting.
- `04-heuristic-collapsed-R.R`: experimental collapsed-R inference with a
  correlation-standardized posterior working matrix in each IBSS update. Its
  reported collapsed-R ELBO is diagnostic rather than guaranteed monotone.
- `credible-set-distance.R`: generic helpers for the credible-set allocation
  distance used in method comparisons.
- `sim_2pop.R`: coalescent two-population simulation helpers.

Source files 01--03, in that order, before fitting both primary methods. Source
file 04 afterward when using the experimental heuristic:

```r
source("code/01-collapsed-R.R")
source("code/02-AIW-and-AIW-N.R")
source("code/03-SER-fallback-and-plot.R")
source("code/04-heuristic-collapsed-R.R")
```

The paper-facing simulations are under `experiments/`; workflowr pages are
under `analysis/`.
