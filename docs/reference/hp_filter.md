# Hodrick-Prescott filter with a fixed smoothing parameter

`hp_filter()` splits a time series into a smooth trend and a cycle with
the filter of Hodrick and Prescott (1997) and a smoothing parameter
\\\lambda\\ chosen by the user: 1600 for quarterly data, 100 for annual
data and 129600 for monthly data are the conventional values (Ravn and
Uhlig 2002).

## Usage

``` r
hp_filter(x, lambda = 1600, as_dt = TRUE, log_scale = FALSE)

# Default S3 method
hp_filter(x, lambda = 1600, as_dt = TRUE, log_scale = FALSE)
```

## Arguments

- x:

  a numeric vector or a `ts` object, without missing values; at least 5
  observations.

- lambda:

  the smoothing parameter \\\lambda \> 0\\ (default 1600, for quarterly
  data).

- as_dt:

  `TRUE` (default) returns a `data.table`; `FALSE` a list.

- log_scale:

  filter the logarithm of the series (default `FALSE`, as in version
  0.1.0).

## Value

With `as_dt = TRUE`, a `data.table` with the columns `original`,
`trend`, `cycle` (and `time` for a `ts`) and the attribute `lambda`.
With `as_dt = FALSE`, a list of class
`c("mhp_hp", "mhp", "mhp_result", "list")` with `original`, `trend`,
`cycle`, `lambda`, then `log_scale`, `trend_log` and `cycle_log` (with
`log_scale = TRUE`) and `time` and `frequency` (for a `ts`). Methods:
[`print()`](https://rdrr.io/r/base/print.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html),
[`ggplot2::autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html),
[`mhp_interpret()`](https://myaseen208.com/mhpfilter/reference/mhp_interpret.md).

## Details

The trend \\\mathbf{g} = (g_1, \dots, g_T)^\top\\ of the series
\\\mathbf{x}\\ minimises \$\$\sum\_{t=1}^{T}(x_t - g_t)^2 +
\lambda\sum\_{t=2}^{T-1}\left\[(g\_{t+1} - g_t) - (g_t -
g\_{t-1})\right\]^2,\$\$ so that \\\hat{\mathbf{g}} = (\mathbf{I} +
\lambda\mathbf{A})^{-1}\mathbf{x} = \mathbf{B}\mathbf{x}\\ with
\\\mathbf{A} = \mathbf{K}^\top\mathbf{K}\\ and \\\mathbf{K}\\ the \\(T -
2) \times T\\ matrix of second differences (rows \\(1, -2, 1)\\). The
cycle is \\\mathbf{c} = \mathbf{x} - \hat{\mathbf{g}}\\. The system is
pentadiagonal, symmetric and positive definite and is solved in \\O(T)\\
operations.

With `log_scale = TRUE` the filter is applied to \\\log x_t\\; the trend
is returned as \\\exp(\hat{g}\_t)\\ and the cycle as \\x_t -
\exp(\hat{g}\_t)\\, and the components on the log scale are added to the
list result.

## References

Hodrick, R. J. and Prescott, E. C. (1997). Postwar U.S. business cycles:
an empirical investigation. *Journal of Money, Credit and Banking*,
**29**(1), 1–16. [doi:10.2307/2953682](https://doi.org/10.2307/2953682)

Ravn, M. O. and Uhlig, H. (2002). On adjusting the Hodrick-Prescott
filter for the frequency of observations. *The Review of Economics and
Statistics*, **84**(2), 371–376.
[doi:10.1162/003465302317411604](https://doi.org/10.1162/003465302317411604)

## See also

[`mhp_filter`](https://myaseen208.com/mhpfilter/reference/mhp_filter.md),
[`fmhp_filter`](https://myaseen208.com/mhpfilter/reference/fmhp_filter.md),
[`mhp_compare`](https://myaseen208.com/mhpfilter/reference/mhp_compare.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
data(us_macro)
gdp <- ts(data = us_macro$gdp, start = c(1959, 1), frequency = 4)
hp <-
  hp_filter(
      x         = log(gdp)
    , lambda    = 1600
    , as_dt     = TRUE
    , log_scale = FALSE
    )
hp
#>      original    trend        cycle    time
#>         <num>    <num>        <num>   <num>
#>   1: 7.904833 7.896154  0.008678366 1959.00
#>   2: 7.929775 7.905529  0.024246310 1959.25
#>   3: 7.928582 7.914908  0.013673747 1959.50
#>   4: 7.932076 7.924314  0.007762667 1959.75
#>   5: 7.954267 7.933774  0.020492108 1960.00
#>  ---                                       
#> 199: 9.497367 9.490038  0.007328946 2008.50
#> 200: 9.483562 9.492102 -0.008539432 2008.75
#> 201: 9.466950 9.494061 -0.027110867 2009.00
#> 202: 9.465099 9.495969 -0.030869902 2009.25
#> 203: 9.471961 9.497861 -0.025899315 2009.50
attr(x = hp, which = "lambda")
#> [1] 1600
```
