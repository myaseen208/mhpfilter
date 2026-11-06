# End-point loss of an HP-type filter

`hp_loss()` measures the end-point bias of a filter by the loss of
Bloechl (2014) used by Hanif, Iqbal and Choudhary (2017): for every time
\\t\\, the squared distance of the gain of its estimate from the gain at
the centre of the series. The loss is zero in the middle and grows
towards the ends; the fully modified HP filter chooses its end-point
weights to make the total small.

## Usage

``` r
hp_loss(
  x,
  lambda = 1600,
  form = c("rows", "differences"),
  theta = seq(0, 3.1, 0.1),
  weight = 0.1
)

# Default S3 method
hp_loss(
  x,
  lambda = 1600,
  form = c("rows", "differences"),
  theta = seq(0, 3.1, 0.1),
  weight = 0.1
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

  the frequencies (default \\0, 0.1, \dots, 3.1\\).

- weight:

  the factor \\\theta\\ (default 0.1).

## Value

An object of class `c("mhp_loss", "mhp_result", "list")`: `loss` (a
`data.table` with `t` and `loss`), `total` (\\L\\), `lambda_t`, `form`.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws \\l(t)\\
against \\t\\ (Figures 3 and 3b of the paper).

## Details

\$\$l(t) = \theta\sum\_{\omega}\left\[P_c(\omega) -
P_t(\omega)\right\]^2, \qquad L = \sum\_{t=1}^{T} l(t),\$\$ with the
gains \\P_t\\ of
[`hp_gain()`](https://myaseen208.com/mhpfilter/reference/hp_gain.md),
\\\omega = 0, 0.1, \dots, 3.1\\, \\\theta = 0.1\\ (their spacing) and
\\c = \mathrm{round}(T/2)\\ (Hanif et al. 2017, eqs. 7 and 8, and their
MATLAB code). The loss depends on \\T\\ and the smoothing parameters
only, not on the data: for the HP filter with \\\lambda = 100\\ it is
1.731 to 1.735 for 41 to 55 observations, the values of Table 3a of the
paper.

## References

Bloechl, A. (2014). Reducing the excess variability of the
Hodrick-Prescott filter by flexible penalization. *Discussion Paper*
2014-1, University of Munich.

Hanif, M. N., Iqbal, J. and Choudhary, M. A. (2017). Fully modified HP
filter. *SBP Working Paper Series*, No. 88, State Bank of Pakistan.

## See also

[`fmhp_filter`](https://myaseen208.com/mhpfilter/reference/fmhp_filter.md),
[`hp_gain`](https://myaseen208.com/mhpfilter/reference/hp_gain.md),
[`fmhp_compare`](https://myaseen208.com/mhpfilter/reference/fmhp_compare.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
hl <-
  hp_loss(
      x      = 55
    , lambda = 100
    , form   = "rows"
    , theta  = seq(from = 0, to = 3.1, by = 0.1)
    , weight = 0.1
    )
hl$total
#> [1] 1.734501
plot(x = hl)
```
