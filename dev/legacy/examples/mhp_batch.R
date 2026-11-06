### Name: mhp_batch
### Title: Batch Modified HP Filter
### Aliases: mhp_batch

### ** Examples

## No test: 
# Example 1: Multiple macroeconomic series
set.seed(456)
n_periods <- 60
n_countries <- 5
gdp_matrix <- matrix(nrow = n_periods, ncol = n_countries)
colnames(gdp_matrix) <- c("USA", "UK", "Germany", "France", "Japan")

# Generate series with different characteristics
for (i in 1:n_countries) {
  trend <- cumsum(rnorm(n_periods, mean = 0.5, sd = 0.3))
  cycle <- rnorm(n_periods, sd = 1 + i * 0.2) # Increasing volatility
  gdp_matrix[, i] <- trend + cycle
}

# Apply batch Modified HP filter
results <- mhp_batch(gdp_matrix, max_lambda = 10000)

# Extract optimal lambdas
lambdas <- attr(results, "lambdas")
print(lambdas)

# Example 2: Sectoral data
set.seed(789)
n_time <- 120
n_sectors <- 8
sector_names <- c(
  "Agriculture", "Mining", "Manufacturing", "Construction",
  "Trade", "Transport", "Finance", "Services"
)
sector_data <- matrix(rnorm(n_time * n_sectors), nrow = n_time)

# Add sector-specific trends and cycles
for (i in 1:n_sectors) {
  trend_growth <- runif(1, 0.2, 1.0)
  cycle_amplitude <- runif(1, 0.5, 3.0)
  sector_data[, i] <- cumsum(rnorm(n_time, mean = trend_growth / 4, sd = 0.3)) +
    cycle_amplitude * sin(2 * pi * (1:n_time) / (20 + i * 5))
}
colnames(sector_data) <- sector_names

sector_results <- mhp_batch(sector_data, max_lambda = 50000)

# View results for first few periods
head(sector_results)
## End(No test)




