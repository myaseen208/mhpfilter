# Compare the HP, modified HP and fully modified HP filters

`fmhp_compare()` filters one series with the three filters and compares
their cycles and their end-point loss (Hanif et al. 2017, Tables 3a and
3b): the HP filter at the conventional \\\lambda\\ of the frequency, the
modified HP filter and the fully modified HP filter.

## Usage

``` r
fmhp_compare(
  x,
  frequency = c("quarterly", "annual", "monthly"),
  max_lambda = 3000L,
  log_scale = TRUE,
  end_share = 0.1
)

# Default S3 method
fmhp_compare(
  x,
  frequency = c("quarterly", "annual", "monthly"),
  max_lambda = 3000L,
  log_scale = TRUE,
  end_share = 0.1
)
```

## Arguments

- x:

  a numeric vector or a `ts` object, without missing values; at least 5
  observations.

- frequency:

  `"quarterly"` (default), `"annual"` or `"monthly"`.

- max_lambda:

  the largest \\\lambda\\ of the searches over \\1, 2, \dots,\\
  `max_lambda` (default 3000, as in the published code).

- log_scale:

  filter the logarithm of the series (default `TRUE`, as in the
  published code).

- end_share:

  share of the observations counted as ends on each side (default 0.1,
  as in Table 2 of the paper).

## Value

A `data.table` with one row per filter (`HP`, `Modified HP`,
`Fully modified HP`) and the columns `method`, `lambda`, `cycle_sd`,
`cycle_mean`, `ar1`, `cycle_range`, `gcv`, `loss`, `sd_ends`,
`sd_middle`.

## Details

Besides the moments of
[`mhp_compare()`](https://myaseen208.com/mhpfilter/reference/mhp_compare.md),
the table gives the total loss \\L\\ of each filter
([`hp_loss()`](https://myaseen208.com/mhpfilter/reference/hp_loss.md)),
its standard deviation of the cycle at the ends (the first and last
`end_share` of the observations) and in the middle. A filter without
end-point bias has a loss near zero and about the same cycle variability
at the ends as in the middle.

## See also

[`fmhp_filter`](https://myaseen208.com/mhpfilter/reference/fmhp_filter.md),
[`mhp_compare`](https://myaseen208.com/mhpfilter/reference/mhp_compare.md),
[`hp_loss`](https://myaseen208.com/mhpfilter/reference/hp_loss.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
data(us_macro)
fc <-
  fmhp_compare(
      x          = us_macro$gdp
    , frequency  = "quarterly"
    , max_lambda = 3000
    , log_scale  = TRUE
    , end_share  = 0.1
    )
fc
#>               method lambda   cycle_sd    cycle_mean       ar1 cycle_range
#>               <char>  <num>      <num>         <num>     <num>       <num>
#> 1:                HP   1600 0.01543904  4.942713e-12 0.8668207  0.08590516
#> 2:       Modified HP    683 0.01333176 -1.634336e-13 0.8280266  0.07801575
#> 3: Fully modified HP   1448 0.01618038 -4.145448e-05 0.8815215  0.08713978
#>             gcv      loss    sd_ends  sd_middle
#>           <num>     <num>      <num>      <num>
#> 1:           NA 1.7617778 0.01458541 0.01568213
#> 2: 0.0002819926 1.7652928 0.01202396 0.01366805
#> 3: 0.0003335623 0.7898641 0.01930994 0.01538061
```
