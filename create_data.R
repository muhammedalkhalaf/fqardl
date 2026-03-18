# Create datasets for fqardl
setwd("C:/Users/acad_/.openclaw/workspace/fqardl-r")
set.seed(123)

# Generate macro_data for fqardl examples
n <- 100
x1 <- cumsum(rnorm(n, sd = 0.3))
x2 <- cumsum(rnorm(n, sd = 0.2))
t_idx <- 1:n
fourier <- 0.5 * sin(2 * pi * t_idx / n) + 0.3 * cos(2 * pi * t_idx / n)
y_lr <- 3 + 0.4 * x1 + 0.25 * x2 + fourier
y <- rep(0, n)
y[1] <- y_lr[1]
for (i in 2:n) {
  ec <- y[i-1] - y_lr[i-1]
  dy <- -0.15 * ec + 0.2 * (x1[i] - x1[i-1]) + 0.1 * (x2[i] - x2[i-1]) + rnorm(1, sd = 0.05)
  y[i] <- y[i-1] + dy
}
inflation <- 2 + 0.3 * diff(c(0, y)) + rnorm(n, sd = 0.2)
interest_rate <- 3 + 0.5 * inflation + rnorm(n, sd = 0.3)
macro_data <- data.frame(
  gdp = y, 
  inflation = inflation, 
  interest_rate = interest_rate,
  oil_price = x1,
  exchange_rate = x2
)

# Save
if (!dir.exists("data")) dir.create("data")
save(macro_data, file = "data/macro_data.rda", compress = "xz")
cat("Dataset created!\n")
