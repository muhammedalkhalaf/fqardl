## fqardl 1.1.0

This release corrects the computations below; the version on CRAN is 1.0.6.

It follows 1.0.6 within days because a re-audit found the bootstrap bounds tests and the analytical verdicts for models with Fourier terms invalid and the fnardl and mtnardl bounds F statistic wrong when p or q exceeds 1, so results change; mtnardl() is deprecated in favour of ardlverse::mtnardl().

### Breaking changes

* `fqardl(bootstrap = TRUE)` / `bootstrap_bounds_test()` and `fnardl(bootstrap = TRUE)` / `bootstrap_nardl()` use a new bootstrap (below). Bootstrap p-values, critical values and decisions change. The `decision` element of the bootstrap result is now a code (`"COINTEGRATION"`, `"NO_COINTEGRATION"`, `"DEGENERATE_1"`, `"DEGENERATE_2"`, `"UNDETERMINED"`) with the text in `label`; failed replications are `NA` and counted in `n_valid` instead of being dropped. The `tau` argument of `bootstrap_bounds_test()` is no longer used. `bootstrap_nardl()` (internal) has new arguments: the original regressors, their names and the decomposed variables.
* `perform_bounds_test()` and the bounds test of `fnardl()` give no analytical verdict when the regression contains Fourier terms, which is always the case for `fqardl()` and `fnardl()`; the decision now reads "No analytical verdict: the PSS bounds are for the model without Fourier terms and are not valid with Fourier terms", followed by "; use the bootstrap (bootstrap = TRUE)" when `fqardl()` or `fnardl()` was called without the bootstrap. The per-quantile `decision_F` and `decision_t` read "not available (Fourier terms)" and the `decision` column of `bounds_test$summary` carries the new text. The F and t statistics are unchanged (see below for the fnardl F statistic); the PSS bounds are still returned, labelled "Bounds for the model without Fourier terms; not valid with Fourier terms", with the new elements `bounds_valid` and `bounds_note`. The evidence that the PSS (2001) bounds do not apply with Fourier terms is the Monte Carlo below, not a published result.
* `bootstrap_bounds_test()` no longer emits its "indicative only" warning; its arguments `case` and `tau` now have the defaults 3 and `NULL`.
* The fnardl and mtnardl bounds tests read the PSS table at k equal to the number of level regressors, so `cv_5` and `t_cv_5` change whenever the old count (all `_lag1` terms minus one) differed, and the mtnardl verdict can change with the corrected F and bounds.
* `orig_t` of the fnardl bootstrap result keeps its name `y_lag1` as in 1.0.6; its value is now the OLS t statistic of the same model.
* The header printed by `print()` for mtnardl objects is now "Multi-Threshold NARDL Model (fqardl, deprecated)". `plot()` on an fqardl mtnardl object now stops with an informative error (method `plot.fqardl_mtnardl`) instead of dispatching through the inherited class to the plot method of ardlverse, which failed with "argument is of length zero".
* `mtnardl()` is deprecated and will be removed in fqardl 2.0.0; use `ardlverse::mtnardl()`. The first call in a session gives a warning that also states that the conventions differ (p in ardlverse equals p - 1 here; thresholds are a single numeric vector there, a named list here). Its results change: the bounds F statistic is corrected when the selected p or q exceeds 1 (example below), thresholds other than 0 and cases other than 3 give an error, and the object has class `c("fqardl_mtnardl", "mtnardl")`.
* `decompose_multi_threshold()` gives an error for thresholds other than 0.

### Bootstrap bounds tests (fqardl and fnardl)

* Up to 1.0.6 each replication regressed a new Gaussian random walk on the fixed observed regressors and Fourier terms (`fqardl`: quantile regression at the first quantile; `fnardl`: OLS, t test only, without any warning), and `fqardl` declared cointegration when the F or the t test rejected. This is not a bootstrap under the null of the estimated model and it was oversized (Monte Carlo below).
* The pseudo-samples now come from the shared recursive bootstrap engine `R/ardl_boot_engine.R` (engine version 1.1.0, copied verbatim from the package ardlverse; a test checks its version and checksum). Under each null the restricted conditional ECM and a marginal model for the differenced regressors are estimated by OLS, and y* and x* are generated recursively from resampled, recentred residuals. `boot_scheme = "bvz"` (default) follows Bertelli, Vacca and Zoia (2022): a separate restricted model for each statistic, marginal VECM, recentring after each draw, initial values from a random block, with the levels of the regressors and partial sums continuing from the block; the x equation has as many lags as the y equation has lagged differences of y (BVZ eq. 19). `boot_scheme = "mcnown"` applies the Fov null for all statistics (McNown, Sam and Goh 2018, Steps 1-8) to the conditional ECM (package choice; MSG write the y equation in unconditional form): the null of the overall F test generates the samples for all statistics, unrestricted x equation, residuals mean-centred once.
* The statistics are OLS statistics of the package's own design (`build_ardl_design()` for fqardl, `build_nardl_design()` for fnardl): the overall F test on the lagged levels of y and the regressors (Fov), the t ratio on the lagged level of y (t) and the F test on the lagged levels of the regressors (Find). For fnardl the positive and negative partial sums are rebuilt from x* with the package's own decomposition (starting at zero). The Fourier terms are held fixed in the data generating process.
* The Fourier frequency and the lag orders are re-selected in every replication (a package choice: Bertelli, Vacca and Zoia 2022, step 5(c), re-estimate the unrestricted model on each bootstrap sample but do not discuss re-selection) with the criteria of `select_fourier_frequency()` and `select_optimal_lags()` (computed with `qr()` for speed; a test checks that the choices are identical). The new argument `fourier_k` of `fqardl()` and `fnardl()` fixes the frequency in advance; it is then held fixed in the bootstrap.
* Decision: AND rule. Cointegration only if Fov, t and Find all reject at 5 percent (bootstrap order-statistic critical values); Fov rejecting without t, and Fov and t rejecting without Find, are reported as the degenerate cases of the first and second type. Cointegration is never declared from a single test (no OR rule).
* The quantile ARDL of Cho, Kim and Shin (2015) has no bounds test and no bootstrap for the null of no cointegration, so no quantile-specific bootstrap statistic is computed; the quantile Wald bounds statistics of `perform_bounds_test()` are descriptive.
* New arguments: `fqardl(fourier_k, boot_scheme)`, `fnardl(fourier_k, boot_scheme, seed)`, `bootstrap_bounds_test(reselect, scheme, level, seed)`. A `seed` given to the bootstrap functions is used locally and the caller's random number state is restored. `fqardl(seed)` still calls `set.seed()` at the start, as before.

### Bounds statistics

* `fnardl()`: the bounds F statistic tested every coefficient whose name ended in `_lag1`, so with p > 1 or q > 1 it also tested `dy_lag1` and the lagged differences of the regressors. It now tests only the lagged level of y and the lagged levels of the regressors, with k equal to the number of level regressors. With p = q = 1 the value is unchanged. Examples with data `set.seed(s); n <- 120; x1 <- cumsum(rnorm(n)); x2 <- cumsum(rnorm(n)); y <- 0.5 * x1 + cumsum(rnorm(n, sd = 0.5)); d <- data.frame(y, x1, x2)` and `fnardl(y ~ x1 + x2, d, decompose = "x1", max_k = 3, max_p = 3, max_q = 3)` (selected p = 2, q = 1): s = 1, F 4.9160 becomes 6.1315; s = 2, 2.2814 becomes 2.8463; s = 42, 6.2361 becomes 7.5356. The bounds test also returns `t_cv_5`, `Find_stat` and `k`.
* `mtnardl()`: same correction (defect F1 of the 2026-10-03 audit); with the data above, s = 2, and `mtnardl(y ~ x1, d, max_p = 3, max_q = 3)` (selected p = 2, q = 1), F 1.5051 becomes 1.9533; k for the PSS table is the number of level regressors (partial sums and undecomposed variables), a package choice pending verification against Shin, Yu and Greenwood-Nimmo (2014) and Pal and Mitra (2016). The decision is the F-based three-way verdict, as before.

### mtnardl (deprecated)

* F2: the multi-threshold decomposition was mathematically wrong. With a nonzero threshold (0 is always added) the regime increments did not sum to the change of the variable (for example, with thresholds -2, 0 and 2 the changes -3, -1, 1, 3 and 0.5 gave regime increments summing to -5, 1, 3, 5 and 2.5; with the single threshold 2 the changes -3 and -1 gave -5 and -3). Thresholds other than 0 now give an informative error in `mtnardl()` and `decompose_multi_threshold()`; with threshold 0 (the default) the decomposition is the standard pair of partial sums and results other than the bounds F are unchanged.
* F3: `case` was ignored in the bounds test and forced to 3 in the lag selection; cases other than 3 now give an error.
* The print and summary methods are registered for class `"fqardl_mtnardl"` only, so they no longer mask the methods of ardlverse for its `"mtnardl"` objects when both packages are loaded.

### Printing

* `summary()` of an fnardl object printed the long-run multipliers as `NA` (the elements of `multipliers$long_run_positive` and `long_run_negative` are named `<var>.<var>_pos_lag1` and were indexed by `<var>`). The printout now shows the values; the stored values are unchanged.
* Author names in the documentation and in the headers printed by the Fourier unit root tests are joined with "and" instead of "&"; the help page of `get_pss_critical_values()` now describes the function as implemented since 1.0.4.

### Monte Carlo (size at the 5 percent level)

y and x independent random walks, n = 100, 200 replications, B = 199, default calls (BIC, max_p = max_q = 4, max_k = 3); rejection rates with their Monte Carlo standard errors sqrt(r(1 - r)/200) in parentheses. Before: fqardl 1.0.6 (audit of 2026-10-03, same data seeds).

| Statistic | 1.0.6 | 1.1.0 |
|---|---|---|
| fqardl bootstrap, Fov | 0.085 (0.020) | 0.065 (0.017) |
| fqardl bootstrap, t | 0.120 (0.023) | 0.100 (0.021) |
| fqardl bootstrap, Find | not computed | 0.075 (0.019) |
| fqardl bootstrap, decision | 0.130 (0.024) (F or t) | 0.035 (0.013) (Fov and t and Find) |
| fnardl bootstrap, Fov | not computed | 0.060 (0.017) |
| fnardl bootstrap, t | 0.170 (0.027) | 0.105 (0.022) |
| fnardl bootstrap, Find | not computed | 0.100 (0.021) |
| fnardl bootstrap, decision | none | 0.030 (0.012) (Fov and t and Find) |
| fqardl analytical decision (Fourier) | 0.100 (0.021) | no verdict |
| fnardl analytical decision (Fourier) | 0.475 (0.035) | no verdict |

The table shows the run with engine 1.1.0. A run with engine 1.0.0 on the same data seeds gave identical fqardl rates and, for fnardl, Fov 0.060 (0.017), t 0.120 (0.023), Find 0.105 (0.022), decision 0.035 (0.013); an independent review run gave for fnardl a decision rate of 0.055 (0.016) and t and Find rates of 0.120 (0.023). Across these runs the combined decisions are close to the nominal level (0.030 to 0.055). The individual t tests and the fnardl Find test remain above 5 percent (0.100 to 0.120) with the frequency and lags re-selected in every replication; read them through the combined decision.

### Unchanged

* QARDL point estimates, quantile standard errors and covariances, long-run and short-run multipliers, the quantile Wald bounds F and t statistics, `quantile_wald_test()`, the fnardl and mtnardl coefficients, multipliers and asymmetry tests, the fnardl and mtnardl t statistics, the Fourier frequency and lag selection, and the Fourier unit root tests are identical to 1.0.6 (checked on five simulated data sets).

## Test environments

* Ubuntu 24.04, R 4.3.3 and R-devel, R CMD check --as-cran

## R CMD check results

0 errors | 0 warnings | 0 notes
