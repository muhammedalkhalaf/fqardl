\name{fqardl}
\alias{fqardl}
\title{Fourier Quantile ARDL Estimation}
\description{
Estimates Fourier Quantile ARDL models for analyzing relationships across 
the conditional distribution with Fourier terms to capture structural breaks.
}
\usage{
fqardl(formula, data, tau = c(0.25, 0.5, 0.75), max_p = 4, max_q = 4, 
       max_k = 3, ic = c("AIC", "BIC", "HQ"), bootstrap = FALSE, 
       n_boot = 1000, ...)
}
\arguments{
  \item{formula}{A formula specifying the model (y ~ x1 + x2).}
  \item{data}{A data frame containing the variables.}
  \item{tau}{Vector of quantiles to estimate (default: 0.25, 0.5, 0.75).}
  \item{max_p}{Maximum lag order for dependent variable.}
  \item{max_q}{Maximum lag order for independent variables.}
  \item{max_k}{Maximum Fourier frequency.}
  \item{ic}{Information criterion for lag selection.}
  \item{bootstrap}{Logical; perform bootstrap inference?}
  \item{n_boot}{Number of bootstrap replications.}
  \item{...}{Additional arguments.}
}
\value{
An object of class "fqardl" containing estimation results.
}
\references{
Cho, J. S., Greenwood-Nimmo, M., & Shin, Y. (2023). Quantile Cointegration.
}
\author{Muhammad Alkhalaf}
\examples{
\dontrun{
  result <- fqardl(y ~ x1 + x2, data = mydata, tau = c(0.1, 0.5, 0.9))
  summary(result)
}
}
