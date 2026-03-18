# Create oil_gdp_data for fqardl
setwd("C:/Users/acad_/.openclaw/workspace/fqardl-r")
set.seed(456)

n <- 200
d_oil <- rnorm(n, mean = 0, sd = 1)
oil_price <- 50 + cumsum(d_oil)
d_oil_pos <- pmax(d_oil, 0)
d_oil_neg <- pmin(d_oil, 0)
gdp <- rep(100, n)
for (i in 2:n) {
  gdp[i] <- gdp[i-1] - 0.3 * d_oil_pos[i] + 0.7 * abs(d_oil_neg[i]) + rnorm(1, sd = 0.5)
}
date_seq <- seq.Date(as.Date("1970-01-01"), by = "quarter", length.out = n)
oil_gdp_data <- data.frame(date = date_seq, gdp = gdp, oil_price = oil_price)

save(oil_gdp_data, file = "data/oil_gdp_data.rda", compress = "xz")
cat("oil_gdp_data created!\n")
