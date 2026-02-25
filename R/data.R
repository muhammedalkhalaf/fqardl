#' Simulated Macroeconomic Data with Structural Break
#'
#' A simulated quarterly dataset containing GDP, inflation, and interest rate
#' with a structural break at observation 100.
#'
#' @format A data frame with 200 rows and 4 variables:
#' \describe{
#'   \item{date}{Date of observation (quarterly)}
#'   \item{gdp}{Simulated GDP index}
#'   \item{inflation}{Simulated inflation rate with Fourier component}
#'   \item{interest_rate}{Simulated interest rate}
#' }
#'
#' @details
#' The data is generated with a cointegrating relationship between variables
#' that changes at the break point (observation 100). This makes it suitable
#' for demonstrating FQARDL estimation with structural breaks.
#'
#' @source Simulated data for package demonstration
#'
#' @examples
#' data(macro_data)
#' head(macro_data)
"macro_data"


#' Simulated Oil Price and GDP Data with Asymmetric Effects
#'
#' A simulated quarterly dataset where GDP responds asymmetrically to
#' oil price changes - negative oil shocks have larger effects than
#' positive shocks.
#'
#' @format A data frame with 200 rows and 3 variables:
#' \describe{
#'   \item{date}{Date of observation (quarterly)}
#'   \item{gdp}{Simulated GDP index}
#'   \item{oil_price}{Simulated oil price}
#' }
#'
#' @details
#' The data is generated such that:
#' - Positive oil price increases reduce GDP by factor -0.3
#' - Negative oil price decreases increase GDP by factor +0.7 (asymmetric)
#'
#' This makes the data suitable for demonstrating FNARDL estimation
#' with asymmetric effects.
#'
#' @source Simulated data for package demonstration
#'
#' @examples
#' data(oil_gdp_data)
#' head(oil_gdp_data)
"oil_gdp_data"
