# Extract GCV Value

Extract the generalized cross-validation (GCV) value from Modified HP
filter results.

## Usage

``` r
get_gcv(x)
```

## Arguments

- x:

  Result from
  [`mhp_filter`](https://myaseen208.com/mhpfilter/reference/mhp_filter.md)
  (either data.table or mhp object).

## Value

Numeric. The GCV value at the optimal lambda.

## Examples

``` r
set.seed(123)
result <- mhp_filter(cumsum(rnorm(100)), max_lambda = 10000)
get_gcv(result)
#> [1] 1.346882

# With mhp object
result_obj <- mhp_filter(cumsum(rnorm(100)), max_lambda = 10000, as_dt = FALSE)
get_gcv(result_obj)
#> [1] 1.14044
```
