# Modified HP filter of many series

`mhp_batch()` applies
[`mhp_filter()`](https://myaseen208.com/mhpfilter/reference/mhp_filter.md)
to every column of a matrix or data frame and returns one long table,
with the smoothing parameter of each series in an attribute. Series of
different lengths can share the matrix: leading and trailing missing
values are dropped series by series.

## Usage

``` r
mhp_batch(X, max_lambda = 100000L)

# Default S3 method
mhp_batch(X, max_lambda = 100000L)
```

## Arguments

- X:

  a numeric matrix or data frame, one series per column; missing values
  only at the start or the end.

- max_lambda:

  the largest \\\lambda\\ of the search over \\1, 2, \dots,\\
  `max_lambda` (default 100000).

## Value

A `data.table` with the columns `series`, `t`, `original`, `trend`,
`cycle` and the attribute `lambdas`, a `data.table` with `series`,
`lambda` and `gcv` (as in version 0.1.0).

## See also

[`mhp_filter`](https://myaseen208.com/mhpfilter/reference/mhp_filter.md),
[`fmhp_batch`](https://myaseen208.com/mhpfilter/reference/fmhp_batch.md),
[`batch_compare`](https://myaseen208.com/mhpfilter/reference/mhp_compare.md),
[`plot_batch`](https://myaseen208.com/mhpfilter/reference/plot_batch.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
data(us_macro)
X <- log(as.matrix(us_macro[, c("gdp", "consumption", "investment")]))
mb <-
  mhp_batch(
      X          = X
    , max_lambda = 100000
    )
attr(x = mb, which = "lambdas")
#>         series lambda          gcv
#>         <char>  <int>        <num>
#> 1:         gdp    683 0.0002819926
#> 2: consumption    534 0.0001753795
#> 3:  investment   1009 0.0063778004
```
