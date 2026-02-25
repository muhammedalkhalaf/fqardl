\name{fnardl}
\alias{fnardl}
\title{Fourier Nonlinear ARDL Estimation}
\description{
Estimates Fourier Nonlinear ARDL models with asymmetric cointegration 
following Shin, Yu & Greenwood-Nimmo (2014).
}
\usage{
fnardl(formula, data, decompose, max_p = 4, max_q = 4, max_k = 3,
       ic = c("AIC", "BIC", "HQ"), ...)
}
\arguments{
  \item{formula}{A formula specifying the model.}
  \item{data}{A data frame containing the variables.}
  \item{decompose}{Character vector of variables to decompose into positive/negative.}
  \item{max_p}{Maximum lag order for dependent variable.}
  \item{max_q}{Maximum lag order for independent variables.}
  \item{max_k}{Maximum Fourier frequency.}
  \item{ic}{Information criterion for selection.}
  \item{...}{Additional arguments.}
}
\value{
An object of class "fnardl" containing estimation results.
}
\references{
Shin, Y., Yu, B., & Greenwood-Nimmo, M. (2014). Modelling Asymmetric 
Cointegration and Dynamic Multipliers in a Nonlinear ARDL Framework.
}
\author{Muhammad Alkhalaf}
\examples{
\dontrun{
  result <- fnardl(y ~ x, data = mydata, decompose = "x")
  summary(result)
}
}
