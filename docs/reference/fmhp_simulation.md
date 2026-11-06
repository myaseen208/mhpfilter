# Power comparison of the filters on simulated series

`fmhp_simulation()` repeats the simulation study of Hanif, Iqbal and
Choudhary (2017, Tables 1 and 2): for each design it simulates series
with
[`mhp_simulate()`](https://myaseen208.com/mhpfilter/reference/mhp_simulate.md),
filters them with the HP, modified HP and fully modified HP filters, and
compares the root mean square error of the estimated cycles against the
true one, over the whole series, its middle and its ends.

## Usage

``` r
fmhp_simulation(
  designs = mhp_designs(),
  n_rep = 100L,
  T = 200L,
  trend = c("linear", "nonlinear"),
  aggregate = c("none", "systematic", "sum", "average"),
  lambda = NULL,
  max_lambda = 3000L,
  seed = NULL,
  form = c("rows", "differences")
)
```

## Arguments

- designs:

  designs as in
  [`mhp_designs()`](https://myaseen208.com/mhpfilter/reference/mhp_designs.md)
  (default: all 30).

- n_rep:

  replicates per design (default 100; the paper used 1000).

- T:

  number of quarterly observations (default 200).

- trend, aggregate:

  as in
  [`mhp_simulate()`](https://myaseen208.com/mhpfilter/reference/mhp_simulate.md).

- lambda:

  the HP smoothing parameter (default `NULL`: 1600, or 100 for
  aggregated series).

- max_lambda:

  the largest \\\lambda\\ of the searches (default 3000).

- seed:

  seed of the random numbers (`NULL`: the current state).

- form:

  the form of the end-point weights of the fully modified HP filter
  (default `"rows"`); see
  [`fmhp_filter()`](https://myaseen208.com/mhpfilter/reference/fmhp_filter.md).

## Value

An object of class `c("mhp_simulation", "mhp_result", "list")`: `rmse`
(a `data.table` with one row per design and filter: `model`, `sd_ratio`,
`ar1`, `ar2`, `method`, `rmse`, `rmse_middle`, `rmse_ends`), `power`
(one row per design: the share of replicates in which FMHP beats HP and
MHP) and `settings`.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the RMSE
by design.

## Details

The RMSE of a filter is \\\\T^{-1}\sum_t(\hat{c}\_t - c_t)^2\\^{1/2}\\;
the middle is the central 80% of the observations and the ends the first
and last 10%. `fmhp_better` is the share of replicates in which the FMHP
filter has a smaller RMSE than the HP filter (Table 1 of the paper
reports it in percent). The series are filtered on their own scale
(`log_scale = FALSE`).

## References

Hanif, M. N., Iqbal, J. and Choudhary, M. A. (2017). Fully modified HP
filter. *SBP Working Paper Series*, No. 88, State Bank of Pakistan.

## See also

[`mhp_simulate`](https://myaseen208.com/mhpfilter/reference/mhp_simulate.md),
[`fmhp_filter`](https://myaseen208.com/mhpfilter/reference/fmhp_filter.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
sim <-
  fmhp_simulation(
      designs    = mhp_designs()[c(3, 13, 23)]
    , n_rep      = 2
    , T          = 120
    , trend      = "linear"
    , aggregate  = "none"
    , lambda     = NULL
    , max_lambda = 3000
    , seed       = 1
    , form       = "rows"
    )
sim$power
#> Key: <model>
#>    model sd_ratio   ar1   ar2 fmhp_better_hp fmhp_better_mhp
#>    <int>    <num> <num> <num>          <num>           <num>
#> 1:     3     10.0   1.2  -0.4            0.0             0.0
#> 2:    13      2.0   1.2  -0.4            1.0             0.5
#> 3:    23      0.5   1.2  -0.4            0.5             0.5
```
