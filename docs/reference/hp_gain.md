# Gain functions of an HP-type filter

`hp_gain()` gives the gain of the estimate at every time \\t\\ as a
function of the frequency (Hanif et al. 2017, eq. 6, Figures 2 and 2b):
the factor by which the filter scales an oscillation of that frequency.
In the middle of a series the gain of the trend falls to zero at high
frequencies; at the ends it does not, so the trend there keeps part of
the cycle.

## Usage

``` r
hp_gain(
  x,
  lambda = 1600,
  form = c("rows", "differences"),
  theta = seq(0, pi, length.out = 101),
  rows = NULL
)

# Default S3 method
hp_gain(
  x,
  lambda = 1600,
  form = c("rows", "differences"),
  theta = seq(0, pi, length.out = 101),
  rows = NULL
)
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

- theta:

  the frequencies (default 101 values from 0 to \\\pi\\).

- rows:

  the rows drawn by
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) (default: the
  first, third, two middle and last two odd ones, as in the paper's
  figures).

## Value

An object of class `c("mhp_gain", "mhp_result", "list")`: `gain` (a
`data.table` with `t`, `theta`, `gain` for every time and frequency),
`lambda_t`, `form`, `rows`.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the gains
of the chosen rows.

## Details

With \\b\_{tj}\\ the weights of
[`hp_weights()`](https://myaseen208.com/mhpfilter/reference/hp_weights.md),
\\P_t(\omega) = \|\sum_j b\_{tj} e^{i\omega(j - t)}\|\\. It is computed
from two solves of the smoothing system per frequency.

## See also

[`hp_loss`](https://myaseen208.com/mhpfilter/reference/hp_loss.md),
[`hp_weights`](https://myaseen208.com/mhpfilter/reference/hp_weights.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
g <-
  hp_gain(
      x      = 50
    , lambda = 100
    , form   = "rows"
    , theta  = seq(from = 0, to = pi, length.out = 101)
    , rows   = NULL
    )
plot(x = g)
```
