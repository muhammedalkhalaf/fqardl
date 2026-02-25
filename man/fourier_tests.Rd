\name{fourier_adf_test}
\alias{fourier_adf_test}
\alias{fourier_kpss_test}
\alias{fourier_unit_root_analysis}
\alias{fadf_f_test}
\title{Fourier Unit Root Tests}
\description{
Fourier ADF and KPSS unit root tests with structural breaks following 
Enders & Lee (2012) and Becker, Enders & Lee (2006).
}
\usage{
fourier_adf_test(y, model = c("c", "ct"), max_freq = 3, max_lag = NULL)
fourier_kpss_test(y, model = c("c", "ct"), max_freq = 3)
fourier_unit_root_analysis(y, name = "series")
fadf_f_test(y, model = c("c", "ct"), max_freq = 3)
}
\arguments{
  \item{y}{Numeric vector of the time series.}
  \item{model}{Model type: "c" (constant) or "ct" (constant + trend).}
  \item{max_freq}{Maximum Fourier frequency to consider.}
  \item{max_lag}{Maximum lag order (NULL for automatic selection).}
  \item{name}{Name for the series in output.}
}
\value{
Test results with statistics and critical values.
}
\references{
Enders, W., & Lee, J. (2012). The flexible Fourier form and Dickey-Fuller 
type unit root tests. Economics Letters.

Becker, R., Enders, W., & Lee, J. (2006). A stationarity test in the 
presence of an unknown number of smooth breaks.
}
\author{Muhammad Alkhalaf}
