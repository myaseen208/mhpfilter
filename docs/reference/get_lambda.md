# Extract the smoothing parameter or the criterion of a result

`get_lambda()` returns the smoothing parameter and `get_gcv()` the value
of McDermott's criterion of a result of
[`hp_filter()`](https://myaseen208.com/mhpfilter/reference/hp_filter.md),
[`mhp_filter()`](https://myaseen208.com/mhpfilter/reference/mhp_filter.md)
or
[`fmhp_filter()`](https://myaseen208.com/mhpfilter/reference/fmhp_filter.md),
whether it is a `data.table` (`as_dt = TRUE`) or a list
(`as_dt = FALSE`).

## Usage

``` r
get_lambda(x)

# Default S3 method
get_lambda(x)

get_gcv(x)

# Default S3 method
get_gcv(x)
```

## Arguments

- x:

  a filter result.

## Value

`get_lambda()`: the smoothing parameter; `get_gcv()`: the criterion at
that value.

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
data(us_macro)
m <-
  mhp_filter(
      x          = log(us_macro$gdp)
    , max_lambda = 100000
    , as_dt      = TRUE
    , log_scale  = FALSE
    , criterion  = c("approximate", "exact")
    )
get_lambda(x = m)
#> [1] 683
get_gcv(x = m)
#> [1] 0.0002819926
```
