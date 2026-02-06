# mhpfilter: Fast Modified Hodrick-Prescott Filter

High-performance implementation of the Modified HP Filter for
decomposing time series into trend and cyclical components. Based on the
methodology of Choudhary, Hanif & Iqbal (2014) which uses generalized
cross-validation to automatically select the optimal smoothing parameter
lambda.

## Details

The standard Hodrick-Prescott (1997) filter decomposes a time series
\\y_t\\ into trend \\g_t\\ and cycle \\c_t\\ components by minimizing:

\$\$\sum\_{t=1}^{T}(y_t - g_t)^2 + \lambda
\sum\_{t=1}^{T-2}\[(g\_{t+2} - g\_{t+1}) - (g\_{t+1} - g_t)\]^2\$\$

where \\\lambda\\ is the smoothing parameter that controls the trade-off
between trend smoothness and cycle fit.

The Modified HP Filter (McDermott, 1997) selects \\\lambda\\ optimally
using generalized cross-validation (GCV). The GCV criterion is:

\$\$GCV(\lambda) = \frac{SSR(\lambda)}{T} \left\[ 1 + \frac{2}{T\lambda}
\right\]\$\$

where \\SSR(\lambda)\\ is the sum of squared residuals. The optimal
\\\lambda\\ minimizes this GCV criterion.

## Main Functions

- [`mhp_filter`](https://myaseen208.com/mhpfilter/reference/mhp_filter.md):
  Apply Modified HP filter to single series

- [`hp_filter`](https://myaseen208.com/mhpfilter/reference/hp_filter.md):
  Apply standard HP filter with fixed lambda

- [`mhp_batch`](https://myaseen208.com/mhpfilter/reference/mhp_batch.md):
  Batch process multiple series efficiently

- [`mhp_compare`](https://myaseen208.com/mhpfilter/reference/mhp_compare.md):
  Compare HP vs Modified HP for single series

- [`batch_compare`](https://myaseen208.com/mhpfilter/reference/batch_compare.md):
  Compare HP vs Modified HP for multiple series

- [`autoplot.mhp`](https://myaseen208.com/mhpfilter/reference/autoplot.mhp.md):
  ggplot2 visualization for mhp objects

- [`plot_comparison`](https://myaseen208.com/mhpfilter/reference/plot_comparison.md):
  Compare HP and Modified HP visually

- [`plot_batch`](https://myaseen208.com/mhpfilter/reference/plot_batch.md):
  Visualize batch processing results

- [`get_lambda`](https://myaseen208.com/mhpfilter/reference/get_lambda.md):
  Extract optimal lambda from results

- [`get_gcv`](https://myaseen208.com/mhpfilter/reference/get_gcv.md):
  Extract GCV value from results

## References

Choudhary, M.A., Hanif, M.N., & Iqbal, J. (2014). On smoothing
macroeconomic time series using the modified HP filter. *Applied
Economics*, 46(19), 2205-2214.

Hodrick, R.J., & Prescott, E.C. (1997). Postwar US business cycles: An
empirical investigation. *Journal of Money, Credit and Banking*, 29(1),
1-16.

McDermott, C.J. (1997). Note on the modified Hodrick-Prescott filter.
*IMF Working Paper No. 97/108*.

## See also

Useful links:

- <https://github.com/myaseen208/mhpfilter>

- Report bugs at <https://github.com/myaseen208/mhpfilter/issues>

## Author

**Maintainer**: Muhammad Yaseen <myaseen208@gmail.com>
([ORCID](https://orcid.org/0000-0002-5923-1714))

Other contributors:

- Javed Iqbal <Javed.iqbal6@sbp.org.pk> (Original methodology author)
  \[contributor\]

- M. Nadim Hanif <Nadeem.hanif@sbp.org.pk> (Original methodology author)
  \[contributor\]
