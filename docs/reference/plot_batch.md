# Plot the results of many series

`plot_batch()` draws the cycles (or trends) of a result of
[`mhp_batch()`](https://myaseen208.com/mhpfilter/reference/mhp_batch.md)
or
[`fmhp_batch()`](https://myaseen208.com/mhpfilter/reference/fmhp_batch.md),
one panel per series or all in one panel, with the smoothing parameter
of each series in its label.

## Usage

``` r
plot_batch(x, show = c("cycle", "trend"), facet = TRUE, highlight = NULL)
```

## Arguments

- x:

  a result of
  [`mhp_batch()`](https://myaseen208.com/mhpfilter/reference/mhp_batch.md)
  or
  [`fmhp_batch()`](https://myaseen208.com/mhpfilter/reference/fmhp_batch.md).

- show:

  `"cycle"` (default) or `"trend"`.

- facet:

  one panel per series (default `TRUE`).

- highlight:

  names of series drawn in full colour, the others faded (default
  `NULL`: all).

## Value

A `ggplot` object.

## Details

**Change from version 0.1.0.** The faceted plot passed `linewidth` twice
(a warning) and ignored `highlight`; both are corrected.

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
data(us_macro)
mb <-
  mhp_batch(
      X          = log(as.matrix(us_macro[, c("gdp", "consumption", "investment")]))
    , max_lambda = 100000
    )
plot_batch(
    x         = mb
  , show      = "cycle"
  , facet     = TRUE
  , highlight = NULL
  )
```
