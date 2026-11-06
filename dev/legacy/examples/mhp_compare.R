### Name: mhp_compare
### Title: Compare HP vs Modified HP Filter
### Aliases: mhp_compare

### ** Examples

# Example 1: Quarterly GDP-like series
set.seed(789)
n <- 100
gdp <- cumsum(rnorm(n, mean = 0.7, sd = 0.5)) + 2 * cos(2 * pi * (1:n) / 32)
comparison <- mhp_compare(gdp, frequency = "quarterly", max_lambda = 10000)
print(comparison)

# Example 2: Annual series
set.seed(101)
n_annual <- 50
annual_series <- cumsum(rnorm(n_annual, mean = 2.0, sd = 1.0)) +
  3 * sin(2 * pi * (1:n_annual) / 10)
annual_comparison <- mhp_compare(annual_series, frequency = "annual", max_lambda = 5000)
print(annual_comparison)

# Example 3: Visual comparison
set.seed(2023)
test_series <- cumsum(rnorm(120, mean = 0.5, sd = 0.4)) +
  runif(1, 1, 3) * sin(2 * pi * (1:120) / 30)

comp_result <- mhp_compare(test_series, frequency = "quarterly", max_lambda = 20000)

if (require(ggplot2)) {
  # Create visualization
  hp_result <- hp_filter(test_series, lambda = 1600, as_dt = FALSE)
  mhp_result <- mhp_filter(test_series, max_lambda = 20000, as_dt = FALSE)

  plot_comparison(test_series, frequency = "quarterly", max_lambda = 20000)
}




