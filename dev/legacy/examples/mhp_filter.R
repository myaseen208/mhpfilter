### Name: mhp_filter
### Title: Modified Hodrick-Prescott Filter
### Aliases: mhp_filter

### ** Examples

# Simulate a trend + cycle series
set.seed(42)
n <- 100
trend <- cumsum(c(0, rnorm(n - 1, mean = 0.5, sd = 0.2)))
cycle <- 2 * sin(2 * pi * (1:n) / 20) + rnorm(n, sd = 0.5)
y <- trend + cycle

# Apply Modified HP filter
result <- mhp_filter(y, max_lambda = 10000)

# Extract optimal lambda
get_lambda(result)

# Extract GCV value
get_gcv(result)

# Print summary
print(result)

# Plot with ggplot2
if (require(ggplot2)) {
  autoplot(mhp_filter(y, max_lambda = 10000, as_dt = FALSE))
}




