# Global variables to avoid R CMD check notes
# These are column names used in dplyr/tidyr/ggplot2 pipelines

utils::globalVariables(c(
  # Column names used in plotting
  "value", "type", "estimate", "se", "tau", "variable", "coefficient",
  "multiplier", "horizon", "persistence", "positive", "negative", 
  "asymmetric", "density",
  

  # Other symbols
  "."
))
