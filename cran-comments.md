## fqardl 1.1.1

CRAN version: 1.1.0.

### Changes

* test_asymmetry(): accepts the object returned by fnardl() as well as its element nardl_result; the long-run test is the Wald test of L(+) = L(-) with the delta-method variance; short-run symmetry tests (pairwise and additive) added to its result.
* New function dynamic_multipliers(): cumulative dynamic multipliers of Shin, Yu and Greenwood-Nimmo (2014, eq. 9.13).
* plot_dynamic_multipliers(), plot_cumulative_multipliers(): plot the multipliers of dynamic_multipliers().
* compute_asymmetric_multipliers(): elements named by the decomposed variable.
* summary() of fnardl objects prints the short-run symmetry tests.
* Man pages and tests updated. tests/testthat/test-bootstrap-1.1.0.R: the comparison of two bootstrap runs with the same seed is made with expect_equal(tolerance = 1e-10), not expect_identical().

### Test environments

* Local: R 4.3.3 on Linux (x86_64).

### R CMD check results

R CMD check --as-cran on R 4.3.3 (Linux, LANG=C.UTF-8, _R_CHECK_SYSTEM_CLOCK_=0): 0 errors, 0 warnings, 1 NOTE. The NOTE reports the GitHub URLs in DESCRIPTION with status 403; the build environment reaches GitHub only through a proxy, and the URLs are those of the published 1.1.0.
