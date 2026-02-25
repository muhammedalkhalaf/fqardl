# Script to create example datasets for fqardl package

# =============================================================================
# Example 1: Cointegrated series with structural break
# =============================================================================

set.seed(2026)
n <- 200

t <- 1:n
break_point <- 100

# Independent variables
x1 <- cumsum(rnorm(n)) + 0.5 * sin(2 * pi * t / n)
x2 <- cumsum(rnorm(n))

# Dependent variable with regime change
y <- numeric(n)
y[1] <- 0
for (i in 2:n) {
  if (i < break_point) {
    y[i] <- 0.3 * y[i-1] + 0.5 * x1[i] + 0.3 * x2[i] + rnorm(1, sd = 0.5)
  } else {
    y[i] <- 0.3 * y[i-1] + 0.8 * x1[i] + 0.2 * x2[i] + rnorm(1, sd = 0.5)
  }
}

macro_data <- data.frame(
  date = seq(as.Date("2010-01-01"), by = "quarter", length.out = n),
  gdp = y,
  inflation = x1,
  interest_rate = x2
)

# Save
usethis::use_data(macro_data, overwrite = TRUE)

# =============================================================================
# Example 2: Oil price and GDP with asymmetric effects
# =============================================================================

set.seed(2025)
n <- 200

oil_price <- cumsum(rnorm(n, sd = 2)) + 50

gdp <- numeric(n)
gdp[1] <- 100

for (i in 2:n) {
  d_oil <- oil_price[i] - oil_price[i-1]
  if (d_oil > 0) {
    # Positive oil shock: smaller negative effect
    gdp[i] <- 0.95 * gdp[i-1] - 0.3 * d_oil + rnorm(1, sd = 0.5)
  } else {
    # Negative oil shock: larger positive effect
    gdp[i] <- 0.95 * gdp[i-1] - 0.7 * d_oil + rnorm(1, sd = 0.5)
  }
}

oil_gdp_data <- data.frame(
  date = seq(as.Date("2010-01-01"), by = "quarter", length.out = n),
  gdp = gdp,
  oil_price = oil_price
)

# Save
usethis::use_data(oil_gdp_data, overwrite = TRUE)

cat("Example datasets created successfully!\n")
