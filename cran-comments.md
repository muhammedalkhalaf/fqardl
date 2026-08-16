# cran-comments.md — fqardl 1.0.3

## Purpose of this submission

This is a **correctness release**. Version 1.0.2 estimated the wrong equation and reported
long-run coefficients with the wrong sign, together with a cointegration statistic that was
not the statistic its critical values were simulated for. Results change materially. The
release is being submitted promptly because the current CRAN version can produce sign-reversed
inference in published work.

The defect was reported by a user, Prof. Chi-Lu (Edward) Peng, National Kaohsiung University
of Science and Technology.

## What changed

1. `build_ardl_design()`, `estimate_nardl()` and `estimate_mtnardl()` used a conditional
   error-correction right-hand side with the dependent variable in levels rather than in first
   differences. The coefficient reported as the error-correction term was therefore `1 + rho`
   instead of `rho`, positive at every quantile, and the long-run multipliers `-theta/phi`
   carried the wrong sign.
2. The bounds "F statistic" was `mean(t^2)` over the lagged level terms, an average of squared
   marginal t ratios, compared against tables that Pesaran, Shin and Smith (2001) simulated for
   a Wald statistic. It is now the Wald form of their equation (21).
3. The bounds verdict now reports the inconclusive region, and is computed at every quantile
   rather than only at the median.
4. `case` was accepted and never used; any value other than 3 now raises an error.
5. `std_errors` and `t_statistics` were unnamed while `coefficients` was renamed, so lookups
   such as `std_errors["y_lag1"]` returned `NA` silently.

Full details are in `NEWS.md`.

## Reverse dependencies

None. `revdepcheck` was not needed.

## Test environments

- Ubuntu 22.04, R 4.4.3, quantreg 6.1.

`R CMD check --as-cran` gives no ERRORs and no WARNINGs. The NOTEs are:

- *Packages suggested but not available for checking: 'testthat', 'knitr', 'rmarkdown',
  'plotly', 'covr'.* These have compiled dependencies that could not be installed in the
  check environment. The `testthat` suite is included in the tarball and passes; see below.
- *unable to verify current time.* Environment artefact, no network clock in the container.

## Tests

`tests/testthat/test-regression-1.0.3.R` was added specifically so that each defect fixed here
cannot silently return. Eleven test blocks, thirty-two assertions, all passing. In particular:

- the design's left-hand side is asserted to equal `Delta y_t`;
- the package's bounds F is asserted to equal a hand-computed Wald statistic;
- the estimator is asserted to recover a known long-run coefficient (`0.80`) and a known
  adjustment speed (`-0.35`) from simulated data, over fifteen replications;
- the shipped Case III critical values are asserted to match Tables CI(iii) and CII(iii) at
  the 5 percent level for every `k` from 1 to 10.

## Known limitations retained and now documented

Three routines have limitations that are documented in `NEWS.md` and that now emit a warning at
the point of use rather than reporting a number that cannot be justified: the bootstrap does
not generate pseudo-samples under the null; the cross-quantile Wald test assumes independence
across quantiles and is conservative; and the quantile-specific bounds verdicts use critical
values simulated for the conditional mean. These are scheduled for version 2.0.0.

## Note to the CRAN team

If you would prefer that version 1.0.2 be archived immediately rather than superseded, please
let me know and I will request it. The concern is that the current version reports long-run
coefficients with a reversed sign, which a user may not notice.
