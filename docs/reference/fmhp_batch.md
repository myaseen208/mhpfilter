# Fully modified HP filter of many series

`fmhp_batch()` applies
[`fmhp_filter()`](https://myaseen208.com/mhpfilter/reference/fmhp_filter.md)
to every column of a matrix or data frame, as
[`mhp_batch()`](https://myaseen208.com/mhpfilter/reference/mhp_batch.md)
does for the modified HP filter.

## Usage

``` r
fmhp_batch(
  X,
  max_lambda = 3000L,
  k_linear = 15L,
  k_quadratic = NULL,
  alpha_linear = seq(0, 2000, 10),
  alpha_quadratic = seq(0, 2000, 5),
  schemes = c("linear", "quadratic"),
  form = c("rows", "differences"),
  lambda = NULL,
  log_scale = TRUE
)

# Default S3 method
fmhp_batch(
  X,
  max_lambda = 3000L,
  k_linear = 15L,
  k_quadratic = NULL,
  alpha_linear = seq(0, 2000, 10),
  alpha_quadratic = seq(0, 2000, 5),
  schemes = c("linear", "quadratic"),
  form = c("rows", "differences"),
  lambda = NULL,
  log_scale = TRUE
)
```

## Arguments

- X:

  a numeric matrix or data frame, one series per column; missing values
  only at the start or the end.

- max_lambda:

  the largest \\\lambda\\ of the searches over \\1, 2, \dots,\\
  `max_lambda` (default 3000, as in the published code).

- k_linear:

  number of end values with linear weights (default 15).

- k_quadratic:

  number of end values with quadratic weights (default `NULL`:
  \\\mathrm{round}(\min(22, T/2 - 5))\\).

- alpha_linear, alpha_quadratic:

  the grids of the weight \\\alpha\\ (defaults `seq(0, 2000, 10)` and
  `seq(0, 2000, 5)`).

- schemes:

  the schemes compared: `c("linear", "quadratic")` (default) or one of
  them.

- form:

  `"rows"` (default, the published code) or `"differences"`; see
  Details.

- lambda:

  `NULL` (default: chosen from the data) or a fixed smoothing parameter.

- log_scale:

  filter the logarithm of the series (default `TRUE`, as in the
  published code).

## Value

A `data.table` with the columns `series`, `t`, `original`, `trend`,
`cycle` and the attribute `lambdas`, a `data.table` with `series`,
`lambda`, `gcv`, `lambda_mhp`, `scheme`, `alpha`, `k` and `loss`.

## See also

[`fmhp_filter`](https://myaseen208.com/mhpfilter/reference/fmhp_filter.md),
[`mhp_batch`](https://myaseen208.com/mhpfilter/reference/mhp_batch.md),
[`plot_batch`](https://myaseen208.com/mhpfilter/reference/plot_batch.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
data(us_macro)
X <- as.matrix(us_macro[, c("gdp", "consumption", "investment")])
fb <-
  fmhp_batch(
      X               = X
    , max_lambda      = 3000
    , k_linear        = 15
    , k_quadratic     = NULL
    , alpha_linear    = seq(from = 0, to = 2000, by = 10)
    , alpha_quadratic = seq(from = 0, to = 2000, by = 5)
    , schemes         = c("linear", "quadratic")
    , form            = "rows"
    , lambda          = NULL
    , log_scale       = TRUE
    )
attr(x = fb, which = "lambdas")
#>         series lambda          gcv lambda_mhp    scheme alpha     k      loss
#>         <char>  <int>        <num>      <int>    <char> <num> <int>     <num>
#> 1:         gdp   1448 0.0003335623        683 quadratic    15    22 0.6163442
#> 2: consumption   1261 0.0002248243        534 quadratic    15    22 0.5850792
#> 3:  investment   2961 0.0073465762       1009 quadratic    20    22 0.6728679
```
