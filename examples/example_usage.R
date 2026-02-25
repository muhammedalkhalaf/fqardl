#' =============================================================================
#' FQARDL & FNARDL Package - Usage Examples
#' =============================================================================

# Install dependencies
# install.packages(c("quantreg", "ggplot2", "dplyr", "tidyr", "zoo"))

# Load the package (development mode)
# devtools::load_all("C:/Users/acad_/.openclaw/workspace/fqardl-r")

# =============================================================================
# Example 1: Generate Simulated Data
# =============================================================================

set.seed(123)
n <- 200

# Generate cointegrated time series with structural break
t <- 1:n
break_point <- 100

# Independent variables
x1 <- cumsum(rnorm(n)) + 0.5 * sin(2 * pi * t / n)  # With smooth break
x2 <- cumsum(rnorm(n))

# Dependent variable (cointegrated with x1 and x2)
# Different relationship before and after break
y <- numeric(n)
y[1] <- 0
for (i in 2:n) {
  if (i < break_point) {
    y[i] <- 0.3 * y[i-1] + 0.5 * x1[i] + 0.3 * x2[i] + rnorm(1, sd = 0.5)
  } else {
    y[i] <- 0.3 * y[i-1] + 0.8 * x1[i] + 0.2 * x2[i] + rnorm(1, sd = 0.5)
  }
}

# Create data frame
sim_data <- data.frame(
  y = y,
  x1 = x1,
  x2 = x2
)

# =============================================================================
# Example 2: FQARDL Estimation
# =============================================================================

cat("\n========================================\n")
cat("FQARDL Example\n")
cat("========================================\n\n")

# Estimate FQARDL model
fqardl_result <- fqardl(
  formula = y ~ x1 + x2,
  data = sim_data,
  tau = c(0.25, 0.50, 0.75),  # Lower, median, upper quantiles
  max_k = 3,                   # Max Fourier frequency
  max_p = 4,                   # Max lag for y
  max_q = 4,                   # Max lag for X
  criterion = "BIC",
  bootstrap = FALSE            # Set TRUE for bootstrap test
)

# View summary
summary(fqardl_result)

# Plot coefficients across quantiles
plot(fqardl_result, type = "coefficients")

# Plot multipliers
plot(fqardl_result, type = "multipliers")

# Plot persistence profile
plot_persistence(fqardl_result, horizons = 30)

# =============================================================================
# Example 3: FNARDL Estimation (Asymmetric Effects)
# =============================================================================

cat("\n========================================\n")
cat("FNARDL Example\n")
cat("========================================\n\n")

# Generate data with asymmetric effects
set.seed(456)
n <- 200

# Oil price with both positive and negative shocks
oil_price <- cumsum(rnorm(n, sd = 2))

# GDP responds asymmetrically to oil price
# Negative oil shocks have larger effect
gdp <- numeric(n)
gdp[1] <- 100
for (i in 2:n) {
  d_oil <- oil_price[i] - oil_price[i-1]
  if (d_oil > 0) {
    # Positive oil shock: smaller negative effect on GDP
    gdp[i] <- 0.95 * gdp[i-1] - 0.3 * d_oil + rnorm(1, sd = 0.5)
  } else {
    # Negative oil shock: larger positive effect on GDP
    gdp[i] <- 0.95 * gdp[i-1] - 0.7 * d_oil + rnorm(1, sd = 0.5)
  }
}

asym_data <- data.frame(
  gdp = gdp,
  oil_price = oil_price
)

# Estimate FNARDL model
fnardl_result <- fnardl(
  formula = gdp ~ oil_price,
  data = asym_data,
  decompose = c("oil_price"),  # Decompose oil price into + and -
  max_k = 2,
  max_p = 4,
  max_q = 4,
  criterion = "BIC",
  bootstrap = FALSE
)

# View summary
summary(fnardl_result)

# Plot asymmetric multipliers
plot(fnardl_result, type = "asymmetry")

# Plot dynamic multipliers
plot(fnardl_result, type = "dynamic", variable = "oil_price", horizon = 30)

# Plot cumulative multipliers
plot(fnardl_result, type = "cumulative", variable = "oil_price", horizon = 30)

# =============================================================================
# Example 4: Fourier Unit Root Test
# =============================================================================

cat("\n========================================\n")
cat("Fourier ADF Test\n")
cat("========================================\n\n")

# Test for unit root with structural breaks
fadf_result <- fourier_adf(
  y = sim_data$y,
  max_k = 3,
  max_lag = 8,
  criterion = "BIC"
)

print(fadf_result)

# =============================================================================
# Example 5: Practical Application Template
# =============================================================================

# Load your own data:
# my_data <- read.csv("your_data.csv")
# 
# # FQARDL for quantile analysis
# result <- fqardl(
#   formula = dependent ~ independent1 + independent2,
#   data = my_data,
#   tau = c(0.1, 0.25, 0.5, 0.75, 0.9),
#   max_k = 3,
#   bootstrap = TRUE,
#   n_boot = 1000
# )
# 
# # FNARDL for asymmetric analysis
# result_asym <- fnardl(
#   formula = dependent ~ independent1 + independent2,
#   data = my_data,
#   decompose = c("independent1"),
#   max_k = 3,
#   bootstrap = TRUE
# )

cat("\n========================================\n")
cat("Examples completed successfully!\n")
cat("========================================\n")
