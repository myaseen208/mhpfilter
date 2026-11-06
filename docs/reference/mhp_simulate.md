# Simulate a series with a known trend and cycle

`mhp_simulate()` generates a series \\x_t = g_t + c_t\\ from the data
generating process of Hanif, Iqbal and Choudhary (2017, eqs. 10a and
10b), after Harvey and Jaeger (1993) and Guay and St-Amant (2005): a
random walk trend with drift and an AR(2) cycle, \$\$g_t = \mu +
\tau_t + g\_{t-1} + \varepsilon_t, \qquad c_t = \phi_1 c\_{t-1} + \phi_2
c\_{t-2} + \delta_t,\$\$ with \\\varepsilon_t \sim N(0,
\sigma\_\varepsilon^2)\\ and \\\delta_t \sim N(0, \sigma\_\delta^2)\\.
Quarterly series can be aggregated to annual ones as in the paper.

## Usage

``` r
mhp_simulate(
  T = 200L,
  sd_ratio = 1,
  ar = c(1.2, -0.4),
  sd_cycle = 1,
  drift = 0.5,
  trend = c("linear", "nonlinear"),
  aggregate = c("none", "systematic", "sum", "average"),
  seed = NULL
)
```

## Arguments

- T:

  number of (quarterly) observations (default 200, as in the paper).

- sd_ratio:

  \\\sigma\_\varepsilon / \sigma\_\delta\\ (default 1).

- ar:

  \\(\phi_1, \phi_2)\\ of the cycle, stationary (default
  `c(1.2, -0.4)`).

- sd_cycle:

  \\\sigma\_\delta\\ (default 1).

- drift:

  \\\mu\\ (default 0.5).

- trend:

  `"linear"` (default) or `"nonlinear"`.

- aggregate:

  `"none"` (default), `"systematic"`, `"sum"` or `"average"`.

- seed:

  seed of the random numbers (`NULL`: the current state).

## Value

A `data.table` with the columns `t`, `x`, `trend`, `cycle`.

## Details

The paper does not give the drift \\\mu\\ or the form of the nonlinear
trend term \\\tau_t\\. Here \\\tau_t = 0\\ for `trend = "linear"`, and
\\\tau_t = \mu\sin(2\pi t/T)\\ for `trend = "nonlinear"`, a drift that
rises and falls once over the sample. The cycle starts from a burn-in of
100 values. Aggregation to annual data takes every fourth value
(`"systematic"`), the sum of four consecutive values (`"sum"`) or their
average (`"average"`), for the series and both components.

## References

Guay, A. and St-Amant, P. (2005). Do the Hodrick-Prescott and
Baxter-King filters provide a good approximation of business cycles?
*Annales d'Economie et de Statistique*, **77**, 133–155.
[doi:10.2307/20079119](https://doi.org/10.2307/20079119)

Harvey, A. C. and Jaeger, A. (1993). Detrending, stylized facts and the
business cycle. *Journal of Applied Econometrics*, **8**(3), 231–247.
[doi:10.1002/jae.3950080302](https://doi.org/10.1002/jae.3950080302)

## See also

[`mhp_designs`](https://myaseen208.com/mhpfilter/reference/mhp_designs.md),
[`fmhp_simulation`](https://myaseen208.com/mhpfilter/reference/fmhp_simulation.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
s <-
  mhp_simulate(
      T         = 200
    , sd_ratio  = 1
    , ar        = c(1.2, -0.4)
    , sd_cycle  = 1
    , drift     = 0.5
    , trend     = "linear"
    , aggregate = "none"
    , seed      = 1
    )
s
#>          t          x       trend      cycle
#>      <int>      <num>       <num>      <num>
#>   1:     1   1.679742  -0.1264538  1.8061955
#>   2:     2   1.037125   0.5571895  0.4799353
#>   3:     3   2.046342   0.2215609  1.8247815
#>   4:     4   3.930973   2.3168417  1.6141316
#>   5:     5   6.007540   3.1463495  2.8611906
#>  ---                                        
#> 196:   196 108.256252 104.6517201  3.6045314
#> 197:   197 109.635655 106.5928778  3.0427775
#> 198:   198 107.179656 106.0770304  1.1026259
#> 199:   199 107.442699 106.9890051  0.4536938
#> 200:   200 106.338047 107.1079290 -0.7698824
```
