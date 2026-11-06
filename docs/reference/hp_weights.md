# Weights of an HP-type filter

`hp_weights()` gives the matrix \\\mathbf{B}\\ of the trend
\\\hat{\mathbf{g}} = \mathbf{B}\mathbf{x}\\: row \\t\\ holds the weights
of the estimate at time \\t\\. Symmetric in the middle of the series,
the rows become asymmetric near the ends, the source of the end-point
bias (Hanif et al. 2017, Figures 1 and 1a).

## Usage

``` r
hp_weights(x, lambda = 1600, form = c("rows", "differences"), rows = NULL)

# Default S3 method
hp_weights(x, lambda = 1600, form = c("rows", "differences"), rows = NULL)
```

## Arguments

- x:

  the length \\T\\ of a series, or a result of
  [`hp_filter()`](https://myaseen208.com/mhpfilter/reference/hp_filter.md),
  [`mhp_filter()`](https://myaseen208.com/mhpfilter/reference/mhp_filter.md)
  or
  [`fmhp_filter()`](https://myaseen208.com/mhpfilter/reference/fmhp_filter.md)
  (with `as_dt = FALSE`).

- lambda:

  a smoothing parameter, or one per observation (`form = "rows"`) or per
  second difference (`form = "differences"`); ignored for a result.

- form:

  `"rows"` (default) or `"differences"`; see
  [`fmhp_filter()`](https://myaseen208.com/mhpfilter/reference/fmhp_filter.md).

- rows:

  the rows drawn by
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) (default: the
  first, third, two middle and last two odd ones, as in the paper's
  figures).

## Value

An object of class `c("mhp_weights", "mhp_result", "list")`: `weights`
(the \\T \times T\\ matrix), `lambda_t`, `form`, `rows`.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
chosen rows.

## See also

[`hp_gain`](https://myaseen208.com/mhpfilter/reference/hp_gain.md),
[`hp_loss`](https://myaseen208.com/mhpfilter/reference/hp_loss.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
w <-
  hp_weights(
      x      = 50
    , lambda = 100
    , form   = "rows"
    , rows   = NULL
    )
round(x = w$weights[1:3, 1:6], digits = 3)
#>       [,1]  [,2]  [,3]  [,4]  [,5]  [,6]
#> [1,] 0.362 0.282 0.208 0.145 0.093 0.052
#> [2,] 0.282 0.237 0.190 0.144 0.104 0.070
#> [3,] 0.208 0.190 0.169 0.143 0.114 0.087
plot(x = w)
```
