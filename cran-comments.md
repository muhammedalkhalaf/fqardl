## fqardl 1.0.6

This release follows 1.0.5 closely. I apologise for the quick resubmission: an audit of the unit root functions against the published tables found critical values without a source and p-values interpolated between critical values, which give wrong test decisions. The corrections are listed below.

* Fourier ADF test: critical values now follow Enders and Lee (2012, Table 1a) by frequency (k = 1 to 5) and sample size, as tabulated in TSPDLIB. The previous values for k = 2 and 3 had no source, `fourier_adf()` used the k = 1 values for every frequency, and approximate p-values were interpolated between critical values; p-values are no longer reported.
* F test for the Fourier terms: critical values from Enders and Lee (2012, Table 1b) by model and sample size; the p-value from the standard F distribution, which does not apply under the unit root null, is no longer reported.
* Fourier KPSS test: critical values from Becker, Enders and Lee (2006) by frequency and sample size; the previous table mixed values of the two models for k = 2 and had no source for k = 3 or for the trend model. Approximate p-values removed.
* `fourier_adf()` now calls `fourier_adf_test()` (constant model); a duplicate `print.fadf()` method was removed.
* The frequency is limited to k <= 5, the range of the published tables.

## Test environments

* Ubuntu 24.04, R 4.3.3 and R-devel, R CMD check --as-cran

## R CMD check results

0 errors | 0 warnings | 0 notes
