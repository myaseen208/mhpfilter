# Extract Optimal Lambda

Extract the optimal smoothing parameter lambda from Modified HP filter
results.

## Usage

``` r
get_lambda(x)
```

## Arguments

- x:

  Result from
  [`mhp_filter`](https://myaseen208.com/mhpfilter/reference/mhp_filter.md)
  (either data.table or mhp object).

## Value

Integer. The optimal lambda value.

## Examples

``` r
set.seed(123)
result <- mhp_filter(cumsum(rnorm(100)), max_lambda = 10000)
get_lambda(result)
#> [1] 562

# With mhp object
result_obj <- mhp_filter(cumsum(rnorm(100)), max_lambda = 10000, as_dt = FALSE)
get_lambda(result_obj)
#> [1] 896
```
