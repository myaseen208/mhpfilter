# mhpfilter: Modified and Fully Modified Hodrick-Prescott Filters

mhpfilter splits a macroeconomic time series into a trend and a cycle
with the filter of Hodrick and Prescott (1997) and with two versions
whose settings come from the data: the modified HP filter of McDermott
(1997), which chooses the smoothing parameter of each series by
generalised cross-validation (Choudhary, Hanif and Iqbal 2014), and the
fully modified HP filter of Hanif, Iqbal and Choudhary (2017), which
also reduces the bias of the HP filter at the ends of the series by
raising the smoothing parameter there, with weights that minimise the
end-point loss of Bloechl (2014).

## Details

**Functions**, in the order of an analysis:

- [`hp_filter`](https://myaseen208.com/mhpfilter/reference/hp_filter.md),
  [`mhp_filter`](https://myaseen208.com/mhpfilter/reference/mhp_filter.md),
  [`fmhp_filter`](https://myaseen208.com/mhpfilter/reference/fmhp_filter.md):
  the three filters.

- [`mhp_batch`](https://myaseen208.com/mhpfilter/reference/mhp_batch.md),
  [`fmhp_batch`](https://myaseen208.com/mhpfilter/reference/fmhp_batch.md):
  many series at once.

- [`mhp_compare`](https://myaseen208.com/mhpfilter/reference/mhp_compare.md),
  [`batch_compare`](https://myaseen208.com/mhpfilter/reference/mhp_compare.md),
  [`fmhp_compare`](https://myaseen208.com/mhpfilter/reference/fmhp_compare.md):
  the filters compared, by the moments of the cycles and the end-point
  loss.

- [`hp_weights`](https://myaseen208.com/mhpfilter/reference/hp_weights.md),
  [`hp_gain`](https://myaseen208.com/mhpfilter/reference/hp_gain.md),
  [`hp_loss`](https://myaseen208.com/mhpfilter/reference/hp_loss.md),
  [`hp_schedule`](https://myaseen208.com/mhpfilter/reference/hp_schedule.md):
  the weights, gains and end-point loss of any HP-type filter.

- [`mhp_designs`](https://myaseen208.com/mhpfilter/reference/mhp_designs.md),
  [`mhp_simulate`](https://myaseen208.com/mhpfilter/reference/mhp_simulate.md),
  [`fmhp_simulation`](https://myaseen208.com/mhpfilter/reference/fmhp_simulation.md):
  simulated series and the power comparison of the paper.

- [`autoplot.mhp`](https://myaseen208.com/mhpfilter/reference/autoplot.mhp.md),
  [`plot_comparison`](https://myaseen208.com/mhpfilter/reference/plot_comparison.md),
  [`plot_batch`](https://myaseen208.com/mhpfilter/reference/plot_batch.md):
  graphics;
  [`get_lambda`](https://myaseen208.com/mhpfilter/reference/get_lambda.md),
  [`get_gcv`](https://myaseen208.com/mhpfilter/reference/get_lambda.md):
  accessors.

- [`mhp_interpret`](https://myaseen208.com/mhpfilter/reference/mhp_interpret.md),
  [`mhp_cite`](https://myaseen208.com/mhpfilter/reference/mhp_cite.md):
  a plain-language interpretation and the references of a result;
  [`mhp_app`](https://myaseen208.com/mhpfilter/reference/mhp_app.md):
  the built-in application.

Every filter is an S3 generic with a `.default` method. With
`as_dt = TRUE` (the default) the filters return a `data.table`, as in
version 0.1.0; with `as_dt = FALSE` a list with a class and the methods
of
[`mhp-methods`](https://myaseen208.com/mhpfilter/reference/mhp-methods.md).
Utilities:
[`mhp_theme`](https://myaseen208.com/mhpfilter/reference/mhp_theme.md),
[`mhp_save`](https://myaseen208.com/mhpfilter/reference/mhp_save.md),
[`mhp_notation`](https://myaseen208.com/mhpfilter/reference/mhp_notation.md),
[`mhp_reference`](https://myaseen208.com/mhpfilter/reference/mhp_reference.md).

**Data.**
[`us_macro`](https://myaseen208.com/mhpfilter/reference/us_macro.md):
quarterly real GDP, consumption and investment of the United States,
1959 to 2009.

**Vignettes.** Start with
[`vignette("mhpfilter_01_getting_started", package = "mhpfilter")`](https://myaseen208.com/mhpfilter/articles/mhpfilter_01_getting_started.md);
`browseVignettes("mhpfilter")` lists all of them, in HTML and PDF.

## Notation

Scalars are upper-case letters in normal face, indices lower-case
italic, vectors lower-case bold and matrices upper-case bold; Greek
letters follow the same rule. The full table is
[`mhp_notation()`](https://myaseen208.com/mhpfilter/reference/mhp_notation.md).

|  |  |
|----|----|
| \\T\\ | number of observations |
| \\\mathbf{x}\\, \\\mathbf{g}\\, \\\mathbf{c}\\ | series, trend and cycle |
| \\\mathbf{K}\\, \\\mathbf{A} = \mathbf{K}^\top\mathbf{K}\\ | second differences and the penalty matrix |
| \\\mathbf{B}\\ | smoother (filter weights) |
| \\\lambda\\, \\\lambda_t\\ | smoothing parameter(s) |
| \\\alpha\\, \\k\\ | end-point weight and number of end values |
| \\P_t(\omega)\\, \\l(t)\\, \\L\\ | gain, loss and total loss |

## References

Bloechl, A. (2014). Reducing the excess variability of the
Hodrick-Prescott filter by flexible penalization. *Discussion Paper*
2014-1, University of Munich.

Choudhary, M. A., Hanif, M. N. and Iqbal, J. (2014). On smoothing
macroeconomic time series using the modified HP filter. *Applied
Economics*, **46**(19), 2205–2214.
[doi:10.1080/00036846.2014.894631](https://doi.org/10.1080/00036846.2014.894631)

Hanif, M. N., Iqbal, J. and Choudhary, M. A. (2017). Fully modified HP
filter. *SBP Working Paper Series*, No. 88, State Bank of Pakistan.

Hodrick, R. J. and Prescott, E. C. (1997). Postwar U.S. business cycles:
an empirical investigation. *Journal of Money, Credit and Banking*,
**29**(1), 1–16. [doi:10.2307/2953682](https://doi.org/10.2307/2953682)

McDermott, C. J. (1997). An automatic method for choosing the smoothing
parameter in the HP filter. Unpublished, International Monetary Fund,
Washington DC.

## See also

[`mhp-methods`](https://myaseen208.com/mhpfilter/reference/mhp-methods.md),
[`mhp_notation`](https://myaseen208.com/mhpfilter/reference/mhp_notation.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>), ORCID 0000-0002-5923-1714.
