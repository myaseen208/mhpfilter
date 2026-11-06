### Name: hp_filter
### Title: Standard Hodrick-Prescott Filter
### Aliases: hp_filter

### ** Examples

# Example 1: Simple random walk with cycle
set.seed(123)
n <- 80
y <- cumsum(rnorm(n)) + sin((1:n) * pi / 10)
result <- hp_filter(y, lambda = 1600)
head(result)

# Example 2: GDP-like series
set.seed(456)
gdp <- cumsum(rnorm(100, mean = 0.5, sd = 0.3)) + 2 * cos(2 * pi * (1:100) / 40)
gdp_decomp <- hp_filter(gdp, lambda = 1600)

# Plot the decomposition
if (require(ggplot2)) {
  plot_data <- data.table::data.table(
    t = 1:length(gdp),
    Original = gdp,
    Trend = gdp_decomp$trend,
    Cycle = gdp_decomp$cycle
  )
  plot_data_long <- data.table::melt(plot_data, id.vars = "t")

  ggplot2::ggplot(plot_data_long, ggplot2::aes(x = t, y = value, color = variable)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::facet_wrap(~variable, ncol = 1, scales = "free_y") +
    ggplot2::labs(
      title = "HP Filter Decomposition (lambda = 1600)",
      x = "Time", y = "Value"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
}




