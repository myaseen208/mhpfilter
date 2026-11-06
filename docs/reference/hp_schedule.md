# Smoothing parameters with raised end values

`hp_schedule()` gives the smoothing parameter of every observation (or
second difference) under the end-point schemes of Hanif, Iqbal and
Choudhary (2017, eqs. 9a and 9b): \\\lambda\\ in the middle and
\\\lambda + \alpha j^{i}\\ for the \\j\\-th of the \\k\\ values at each
end, counted from the middle (\\j = k\\ at the first and last
observation).

## Usage

``` r
hp_schedule(T, lambda, alpha, k, power = 1, form = c("rows", "differences"))
```

## Arguments

- T:

  length of the series.

- lambda:

  the base smoothing parameter.

- alpha:

  the weight \\\alpha \ge 0\\.

- k:

  number of raised values at each end.

- power:

  \\i\\: 1 (linear, default) or 2 (quadratic); any positive number is
  accepted.

- form:

  `"rows"` (default: \\T\\ values, one per observation, as in the
  published code) or `"differences"` (\\T - 2\\ values, one per second
  difference).

## Value

A numeric vector of smoothing parameters.

## See also

[`fmhp_filter`](https://myaseen208.com/mhpfilter/reference/fmhp_filter.md),
[`hp_loss`](https://myaseen208.com/mhpfilter/reference/hp_loss.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
hp_schedule(
    T      = 50
  , lambda = 100
  , alpha  = 10
  , k      = 20
  , power  = 1
  , form   = "rows"
  )
#>  [1] 300 290 280 270 260 250 240 230 220 210 200 190 180 170 160 150 140 130 120
#> [20] 110 100 100 100 100 100 100 100 100 100 100 110 120 130 140 150 160 170 180
#> [39] 190 200 210 220 230 240 250 260 270 280 290 300
```
