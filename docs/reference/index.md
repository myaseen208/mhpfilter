# Package index

## Core Filtering Functions

Main functions for applying HP and Modified HP filters to time series
data. The Modified HP filter automatically selects the optimal smoothing
parameter lambda using generalized cross-validation (GCV).

- [`mhp_filter()`](https://myaseen208.com/mhpfilter/reference/mhp_filter.md)
  : Modified Hodrick-Prescott Filter
- [`hp_filter()`](https://myaseen208.com/mhpfilter/reference/hp_filter.md)
  : Standard Hodrick-Prescott Filter
- [`mhp_compare()`](https://myaseen208.com/mhpfilter/reference/mhp_compare.md)
  : Compare HP vs Modified HP Filter

## Batch Processing

Efficient batch processing functions for multiple time series. Ideal for
cross-country comparisons and analyzing multiple economic indicators
simultaneously.

- [`mhp_batch()`](https://myaseen208.com/mhpfilter/reference/mhp_batch.md)
  : Batch Modified HP Filter
- [`batch_compare()`](https://myaseen208.com/mhpfilter/reference/batch_compare.md)
  : Batch Comparison of HP vs Modified HP

## Visualization Functions

Professional visualization tools using ggplot2. Create publication-ready
plots of trend-cycle decompositions with minimal code.

- [`autoplot(`*`<mhp>`*`)`](https://myaseen208.com/mhpfilter/reference/autoplot.mhp.md)
  : Plot Method for mhp Objects
- [`plot_comparison()`](https://myaseen208.com/mhpfilter/reference/plot_comparison.md)
  : Plot Comparison of HP vs Modified HP
- [`plot_batch()`](https://myaseen208.com/mhpfilter/reference/plot_batch.md)
  : Plot Batch Results

## Utility Functions

Helper functions for extracting results and information from filter
objects.

- [`get_lambda()`](https://myaseen208.com/mhpfilter/reference/get_lambda.md)
  : Extract Optimal Lambda
- [`get_gcv()`](https://myaseen208.com/mhpfilter/reference/get_gcv.md) :
  Extract GCV Value
- [`print(`*`<mhp>`*`)`](https://myaseen208.com/mhpfilter/reference/print.mhp.md)
  : Print Method for mhp Objects

## Package Documentation

Overall package documentation and information.

- [`mhpfilter-package`](https://myaseen208.com/mhpfilter/reference/mhpfilter-package.md)
  [`mhpfilter`](https://myaseen208.com/mhpfilter/reference/mhpfilter-package.md)
  : mhpfilter: Fast Modified Hodrick-Prescott Filter
