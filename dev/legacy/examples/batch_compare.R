### Name: batch_compare
### Title: Batch Comparison of HP vs Modified HP
### Aliases: batch_compare

### ** Examples

# Example 1: Country GDP comparison
set.seed(101)
n <- 80
countries <- c("USA", "UK", "Japan", "Germany", "France", "Italy", "Canada", "Australia")
gdp_data <- sapply(countries, function(ctry) {
  # Varying volatility and persistence
  vol <- runif(1, 0.5, 2.5)
  persist <- runif(1, 0.6, 0.95)
  trend <- cumsum(rnorm(n, 0.5, 0.3))
  cycle <- arima.sim(list(ar = persist), n, sd = vol)
  trend + cycle
})

results <- batch_compare(gdp_data, frequency = "quarterly", max_lambda = 10000)
print(results)

# Example 2: Sectoral analysis with visualization
set.seed(2024)
n_time <- 100
sectors <- c("Tech", "Finance", "Energy", "Healthcare", "Consumer")
sector_returns <- matrix(rnorm(n_time * length(sectors)), nrow = n_time)

# Add sector-specific characteristics
for (i in 1:length(sectors)) {
  drift <- runif(1, -0.1, 0.3)
  volatility <- runif(1, 0.5, 2.0)
  sector_returns[, i] <- cumsum(rnorm(n_time, mean = drift / 100, sd = volatility / 100)) +
    runif(1, 0.5, 2) * sin(2 * pi * (1:n_time) / (20 + i * 3))
}
colnames(sector_returns) <- sectors

sector_comparison <- batch_compare(sector_returns, frequency = "quarterly", max_lambda = 5000)

if (require(ggplot2)) {
  # Plot lambda comparison
  lambda_plot <- ggplot2::ggplot(
    sector_comparison,
    ggplot2::aes(x = series, y = mhp_lambda)
  ) +
    ggplot2::geom_col(fill = "steelblue", alpha = 0.7) +
    ggplot2::geom_hline(yintercept = 1600, linetype = "dashed", color = "red") +
    ggplot2::labs(
      title = "Modified HP Optimal Lambdas by Sector",
      subtitle = "Red line shows fixed HP lambda (1600)",
      x = "Sector", y = "Optimal Lambda"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  print(lambda_plot)
}




