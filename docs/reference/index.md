# Package index

## Filters

HP, modified HP (Choudhary et al. 2014) and fully modified HP (Hanif et
al. 2017) filters.

- [`hp_filter()`](https://myaseen208.com/mhpfilter/reference/hp_filter.md)
  : Hodrick-Prescott filter with a fixed smoothing parameter
- [`mhp_filter()`](https://myaseen208.com/mhpfilter/reference/mhp_filter.md)
  : Modified Hodrick-Prescott filter: the smoothing parameter chosen
  from the data
- [`fmhp_filter()`](https://myaseen208.com/mhpfilter/reference/fmhp_filter.md)
  : Fully modified Hodrick-Prescott filter

## Many series and comparisons

- [`mhp_batch()`](https://myaseen208.com/mhpfilter/reference/mhp_batch.md)
  : Modified HP filter of many series
- [`fmhp_batch()`](https://myaseen208.com/mhpfilter/reference/fmhp_batch.md)
  : Fully modified HP filter of many series
- [`mhp_compare()`](https://myaseen208.com/mhpfilter/reference/mhp_compare.md)
  [`batch_compare()`](https://myaseen208.com/mhpfilter/reference/mhp_compare.md)
  : Compare the HP and the modified HP filter
- [`fmhp_compare()`](https://myaseen208.com/mhpfilter/reference/fmhp_compare.md)
  : Compare the HP, modified HP and fully modified HP filters

## End-point diagnostics

Weights, gains and end-point loss of HP-type filters with any smoothing
parameters.

- [`hp_schedule()`](https://myaseen208.com/mhpfilter/reference/hp_schedule.md)
  : Smoothing parameters with raised end values
- [`hp_weights()`](https://myaseen208.com/mhpfilter/reference/hp_weights.md)
  : Weights of an HP-type filter
- [`hp_gain()`](https://myaseen208.com/mhpfilter/reference/hp_gain.md) :
  Gain functions of an HP-type filter
- [`hp_loss()`](https://myaseen208.com/mhpfilter/reference/hp_loss.md) :
  End-point loss of an HP-type filter

## Simulation

- [`mhp_designs()`](https://myaseen208.com/mhpfilter/reference/mhp_designs.md)
  : The 30 simulation designs of Hanif, Iqbal and Choudhary (2017)
- [`mhp_simulate()`](https://myaseen208.com/mhpfilter/reference/mhp_simulate.md)
  : Simulate a series with a known trend and cycle
- [`fmhp_simulation()`](https://myaseen208.com/mhpfilter/reference/fmhp_simulation.md)
  : Power comparison of the filters on simulated series

## Results, plots and interpretation

- [`get_lambda()`](https://myaseen208.com/mhpfilter/reference/get_lambda.md)
  [`get_gcv()`](https://myaseen208.com/mhpfilter/reference/get_lambda.md)
  : Extract the smoothing parameter or the criterion of a result
- [`print(`*`<mhp>`*`)`](https://myaseen208.com/mhpfilter/reference/print.mhp.md)
  : Print a filter result
- [`autoplot.mhp()`](https://myaseen208.com/mhpfilter/reference/autoplot.mhp.md)
  : Plot the decomposition of a series
- [`plot_comparison()`](https://myaseen208.com/mhpfilter/reference/plot_comparison.md)
  : Plot the HP and modified HP filters of one series
- [`plot_batch()`](https://myaseen208.com/mhpfilter/reference/plot_batch.md)
  : Plot the results of many series
- [`mhp_interpret()`](https://myaseen208.com/mhpfilter/reference/mhp_interpret.md)
  [`print(`*`<mhp_interpretation>`*`)`](https://myaseen208.com/mhpfilter/reference/mhp_interpret.md)
  [`format(`*`<mhp_interpretation>`*`)`](https://myaseen208.com/mhpfilter/reference/mhp_interpret.md)
  : Plain-language interpretation of a result
- [`print(`*`<mhp_result>`*`)`](https://myaseen208.com/mhpfilter/reference/mhp-methods.md)
  [`summary(`*`<mhp_result>`*`)`](https://myaseen208.com/mhpfilter/reference/mhp-methods.md)
  [`print(`*`<summary.mhp_result>`*`)`](https://myaseen208.com/mhpfilter/reference/mhp-methods.md)
  [`format(`*`<mhp_result>`*`)`](https://myaseen208.com/mhpfilter/reference/mhp-methods.md)
  [`as.data.table(`*`<mhp_result>`*`)`](https://myaseen208.com/mhpfilter/reference/mhp-methods.md)
  [`as.data.frame(`*`<mhp_result>`*`)`](https://myaseen208.com/mhpfilter/reference/mhp-methods.md)
  [`plot(`*`<mhp_result>`*`)`](https://myaseen208.com/mhpfilter/reference/mhp-methods.md)
  : Methods for mhpfilter results

## Application

- [`mhp_app()`](https://myaseen208.com/mhpfilter/reference/mhp_app.md)
  [`mhp_app_object()`](https://myaseen208.com/mhpfilter/reference/mhp_app.md)
  : Trend and cycle analysis in a web browser

## Data

- [`us_macro`](https://myaseen208.com/mhpfilter/reference/us_macro.md) :
  United States real GDP, consumption and investment, 1959 to 2009
  (quarterly)

## Utilities

- [`mhp_cite()`](https://myaseen208.com/mhpfilter/reference/mhp_cite.md)
  : References of a result or of the package
- [`mhp_theme()`](https://myaseen208.com/mhpfilter/reference/mhp_theme.md)
  : ggplot2 theme of mhpfilter
- [`mhp_save()`](https://myaseen208.com/mhpfilter/reference/mhp_save.md)
  : Save a figure (vector PDF by default)
- [`mhp_notation()`](https://myaseen208.com/mhpfilter/reference/mhp_notation.md)
  : Notation of mhpfilter
- [`mhp_reference()`](https://myaseen208.com/mhpfilter/reference/mhp_reference.md)
  : Published values used by the tests and vignettes
- [`mhpfilter`](https://myaseen208.com/mhpfilter/reference/mhpfilter-package.md)
  [`mhpfilter-package`](https://myaseen208.com/mhpfilter/reference/mhpfilter-package.md)
  : mhpfilter: Modified and Fully Modified Hodrick-Prescott Filters
