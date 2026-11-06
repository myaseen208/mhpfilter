# Plot the HP and modified HP filters of one series

`plot_comparison()` draws the series with the trends of the HP filter
(at the conventional \\\lambda\\ of its frequency) and of the modified
HP filter, and optionally the cycles; with `fmhp = TRUE` the fully
modified HP filter is added.

## Usage

``` r
plot_comparison(
  x,
  frequency = c("quarterly", "annual", "monthly"),
  max_lambda = 100000L,
  show_cycle = FALSE,
  fmhp = FALSE
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

- show_cycle:

  draw the cycles too (default `FALSE`).

- fmhp:

  add the fully modified HP filter (default `FALSE`; it uses
  `max_lambda = 3000` and the other defaults of
  [`fmhp_filter()`](https://myaseen208.com/mhpfilter/reference/fmhp_filter.md)
  on the scale of the series).

## Value

A `ggplot` object.

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
data(us_macro)
plot_comparison(
    x          = log(us_macro$gdp)
  , frequency  = "quarterly"
  , max_lambda = 100000
  , show_cycle = TRUE
  , fmhp       = TRUE
  )
```
