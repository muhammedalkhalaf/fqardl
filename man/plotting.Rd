\name{plot_dynamic_multipliers}
\alias{plot_dynamic_multipliers}
\alias{plot_cumulative_multipliers}
\alias{generate_fnardl_report}
\title{Plotting Functions for ARDL Models}
\description{
Create visualizations for dynamic and cumulative multipliers from ARDL models.
}
\usage{
plot_dynamic_multipliers(obj, variable, horizon = 20)
plot_cumulative_multipliers(obj, variable, horizon = 20)
generate_fnardl_report(obj, file = "fnardl_report.html", horizon = 20)
}
\arguments{
  \item{obj}{Estimated FNARDL or FQARDL model object.}
  \item{variable}{Variable name for multiplier plots.}
  \item{horizon}{Number of periods for multiplier computation.}
  \item{file}{Output file path for report.}
}
\value{
A ggplot object for plot functions, or file path for report.
}
\author{Muhammad Alkhalaf}
