# Compare the HP and the modified HP filter

`mhp_compare()` filters one series with the HP filter at the
conventional \\\lambda\\ of its frequency and with the modified HP
filter, and compares the cycles; `batch_compare()` does the same for
every column of a matrix.

## Usage

``` r
mhp_compare(
  x,
  frequency = c("quarterly", "annual", "monthly"),
  max_lambda = 100000L
)

# Default S3 method
mhp_compare(
  x,
  frequency = c("quarterly", "annual", "monthly"),
  max_lambda = 100000L
)

batch_compare(
  X,
  frequency = c("quarterly", "annual", "monthly"),
  max_lambda = 100000L
)

# Default S3 method
batch_compare(
  X,
  frequency = c("quarterly", "annual", "monthly"),
  max_lambda = 100000L
)
```

## Arguments

- x:

  a numeric vector or a `ts` object, without missing values; at least 5
  observations.

- frequency:

  `"quarterly"` (default), `"annual"` or `"monthly"`.

- max_lambda:

  the largest \\\lambda\\ of the search over \\1, 2, \dots,\\
  `max_lambda` (default 100000).

- X:

  a numeric matrix or data frame, one series per column.

## Value

`mhp_compare()`: a `data.table` with one row per method (`HP`,
`Modified HP`) and the columns `method`, `lambda`, `cycle_sd`,
`cycle_mean`, `ar1`, `cycle_range`, `gcv`. `batch_compare()`: a
`data.table` with one row per series and the columns `series`,
`hp_lambda`, `mhp_lambda`, `hp_cycle_sd`, `mhp_cycle_sd`, `sd_diff`,
`hp_ar1`, `mhp_ar1`, `ar1_diff`, `relative_sd`.

## Details

The conventional values are \\\lambda = 1600\\ (quarterly), 100 (annual)
and 129600 (monthly, Ravn and Uhlig 2002). `ar1` is the first-order
autocorrelation of the cycle, \\\sum_t c_t c\_{t+1} / \sum_t c_t^2\\.

## See also

[`fmhp_compare`](https://myaseen208.com/mhpfilter/reference/fmhp_compare.md),
[`plot_comparison`](https://myaseen208.com/mhpfilter/reference/plot_comparison.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
data(us_macro)
cmp <-
  mhp_compare(
      x          = log(us_macro$gdp)
    , frequency  = "quarterly"
    , max_lambda = 100000
    )
cmp
#>         method lambda   cycle_sd    cycle_mean       ar1 cycle_range
#>         <char>  <num>      <num>         <num>     <num>       <num>
#> 1:          HP   1600 0.01543904  4.942713e-12 0.8668207  0.08590516
#> 2: Modified HP    683 0.01333176 -1.634336e-13 0.8280266  0.07801575
#>             gcv
#>           <num>
#> 1:           NA
#> 2: 0.0002819926
bc <-
  batch_compare(
      X          = log(as.matrix(us_macro[, c("gdp", "consumption", "investment")]))
    , frequency  = "quarterly"
    , max_lambda = 100000
    )
bc
#>         series hp_lambda mhp_lambda hp_cycle_sd mhp_cycle_sd      sd_diff
#>         <char>     <num>      <int>       <num>        <num>        <num>
#> 1:         gdp      1600        683  0.01543904   0.01333176 -0.002107273
#> 2: consumption      1600        534  0.01241982   0.01000618 -0.002413639
#> 3:  investment      1600       1009  0.07189806   0.06760443 -0.004293628
#>       hp_ar1   mhp_ar1    ar1_diff relative_sd
#>        <num>     <num>       <num>       <num>
#> 1: 0.8668207 0.8280266 -0.03879413   0.8635101
#> 2: 0.8779254 0.8231665 -0.05475890   0.8056623
#> 3: 0.8143287 0.7914367 -0.02289200   0.9402817
```
