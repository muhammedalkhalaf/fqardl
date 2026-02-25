\name{decompose_variables}
\alias{decompose_variables}
\alias{decompose_multi_threshold}
\alias{generate_fourier_terms}
\alias{select_fourier_frequency}
\alias{compute_asymmetric_multipliers}
\alias{test_asymmetry}
\title{Helper Functions for ARDL Analysis}
\description{
Helper functions for variable decomposition, Fourier term generation, 
and asymmetry testing.
}
\usage{
decompose_variables(data, variables)
decompose_multi_threshold(x, thresholds)
generate_fourier_terms(n, k, cumulative = FALSE)
select_fourier_frequency(y, X, max_k = 3, criterion = c("BIC", "AIC", "HQ"))
compute_asymmetric_multipliers(nardl_result, decompose, x_names)
test_asymmetry(nardl_result, decompose)
}
\arguments{
  \item{data}{Data frame containing variables.}
  \item{variables}{Variable names to decompose.}
  \item{x}{Numeric vector to decompose.}
  \item{thresholds}{Threshold values for decomposition.}
  \item{n}{Number of observations.}
  \item{k}{Fourier frequency.}
  \item{cumulative}{Logical; use cumulative Fourier terms?}
  \item{y}{Dependent variable.}
  \item{X}{Matrix of independent variables.}
  \item{max_k}{Maximum frequency to consider.}
  \item{criterion}{Information criterion for selection.}
  \item{nardl_result}{Estimated NARDL model object.}
  \item{decompose}{Variable name that was decomposed.}
  \item{x_names}{Names of independent variables.}
}
\value{
Processed data or test results depending on function.
}
\author{Muhammad Alkhalaf}
