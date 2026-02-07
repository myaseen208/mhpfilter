# Batch Modified HP Filter

Process multiple time series efficiently using the Modified HP Filter.
Optimized for processing large collections of time series (e.g., panel
data).

## Usage

``` r
mhp_batch(X, max_lambda = 100000L)
```

## Arguments

- X:

  Matrix or data.frame. Each column is a separate time series. Rows
  represent time periods, columns represent different series.

- max_lambda:

  Integer. Maximum lambda value to search. Default 100000.

## Value

A data.table in long format with columns:

- series:

  Series identifier (column name or V1, V2, ...)

- t:

  Time index (1, 2, ..., T)

- original:

  Original series values

- trend:

  Estimated trend component

- cycle:

  Cyclical component

With attribute `"lambdas"` containing a data.table of optimal lambdas
and GCV values for each series:

- series:

  Series identifier

- lambda:

  Optimal lambda for the series

- gcv:

  GCV value at optimal lambda

## Details

This function efficiently processes multiple series by: 1. Pre-computing
the \\K'K\\ matrix once for all series 2. Performing parallelizable grid
search for each series 3. Using optimized C++ routines via RcppArmadillo

## Examples

``` r
# \donttest{
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
#>     series lambda      gcv
#>     <char>  <int>    <num>
#> 1:     USA   2385 1.294325
#> 2:      UK   2678 1.707017
#> 3: Germany  10000 1.775680
#> 4:  France   5393 3.038703
#> 5:   Japan   3065 3.196394

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
#>         series     t  original    trend        cycle
#>         <char> <int>     <num>    <num>        <num>
#> 1: Agriculture     1 0.6635266 1.127298 -0.463771849
#> 2: Agriculture     2 1.2963652 1.484541 -0.188175872
#> 3: Agriculture     3 1.6167897 1.825792 -0.209001765
#> 4: Agriculture     4 2.1350412 2.128569  0.006472319
#> 5: Agriculture     5 2.5815389 2.363185  0.218353706
#> 6: Agriculture     6 3.1001873 2.500176  0.600011407
# }
```
