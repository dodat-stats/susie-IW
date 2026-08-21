# Reproducible experiments

This directory contains scripts, intermediate results, reports, and figures
that support the current manuscript. The workflowr source pages are kept
separately under `analysis/`.

The active experiment families are:

- the final coalescent concentration calibrations;
- the independent coalescent benchmark;
- the two-population ancestry-mismatch stress test;
- and the joint LD-stabilization and concentration experiment used in S1
  Appendix.

Earlier development experiments, private-data analyses, and their development
outputs are stored separately in a local-only archive excluded from version
control. Publication-ready aggregate figures required by the manuscript are
kept under `tex/`.

The main coalescent workflow is:

```sh
Rscript --vanilla experiments/run-paper-coalescent-experiments.R
```
