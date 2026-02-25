\name{perform_bounds_test}
\alias{perform_bounds_test}
\alias{bootstrap_bounds_test}
\alias{get_pss_critical_values}
\title{Bounds Testing for Cointegration}
\description{
PSS bounds testing following Pesaran, Shin & Smith (2001) with bootstrap option.
}
\usage{
perform_bounds_test(model, case = 3)
bootstrap_bounds_test(model, n_boot = 1000, case = 3)
get_pss_critical_values(k, case = 3, level = 0.05)
}
\arguments{
  \item{model}{Estimated ARDL model object.}
  \item{case}{PSS case (1-5), default is 3 (unrestricted intercept, no trend).}
  \item{n_boot}{Number of bootstrap replications.}
  \item{k}{Number of regressors.}
  \item{level}{Significance level.}
}
\value{
Bounds test results with F-statistic and critical values.
}
\references{
Pesaran, M. H., Shin, Y., & Smith, R. J. (2001). Bounds testing approaches 
to the analysis of level relationships.
}
\author{Muhammad Alkhalaf}
