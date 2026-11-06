# Plain-language interpretation of a result

`mhp_interpret()` turns a result of the package into a short report: a
headline, findings, recommendations and cautions, computed from its
numbers, with the sources of the method. It covers the filters
([`hp_filter()`](https://myaseen208.com/mhpfilter/reference/hp_filter.md),
[`mhp_filter()`](https://myaseen208.com/mhpfilter/reference/mhp_filter.md),
[`fmhp_filter()`](https://myaseen208.com/mhpfilter/reference/fmhp_filter.md),
as lists or as `data.table`s), the end-point loss
([`hp_loss()`](https://myaseen208.com/mhpfilter/reference/hp_loss.md))
and the power comparison
([`fmhp_simulation()`](https://myaseen208.com/mhpfilter/reference/fmhp_simulation.md)).

## Usage

``` r
mhp_interpret(x, top = 3L, ...)

# Default S3 method
mhp_interpret(x, top = 3L, ...)

# S3 method for class 'mhp_interpretation'
print(x, ...)

# S3 method for class 'mhp_interpretation'
format(x, style = c("text", "markdown"), width = getOption("width", 80L), ...)
```

## Arguments

- x:

  a result of a function of the package.

- top:

  number of largest positive and negative cycle values named (default
  3).

- ...:

  unused.

- style:

  `"text"` (default) or `"markdown"`.

- width:

  line width of the text (default `getOption("width")`).

## Value

An object of class `"mhp_interpretation"`: a list with `analysis`,
`title`, `headline`, `findings`, `recommendations`, `cautions`, `facts`
(the numbers behind the sentences), `guide` (how to read the method) and
`sources`. [`format()`](https://rdrr.io/r/base/format.html) gives it as
text or Markdown.

## Details

For a filter the report compares the chosen smoothing parameter with the
conventional one of the frequency (1600 quarterly, 100 annual, 129600
monthly; known when the series is a `ts`), names the largest booms and
recessions of the cycle, and warns when a search ended at its limit or
the series is short. For the fully modified HP filter it also gives the
reduction of the end-point loss against the HP filter with the same
\\\lambda\\.

## See also

[`mhp_cite`](https://myaseen208.com/mhpfilter/reference/mhp_cite.md),
[`mhp_app`](https://myaseen208.com/mhpfilter/reference/mhp_app.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
data(us_macro)
gdp <- ts(data = us_macro$gdp, start = c(1959, 1), frequency = 4)
f <-
  fmhp_filter(
      x               = gdp
    , max_lambda      = 3000
    , k_linear        = 15
    , k_quadratic     = NULL
    , alpha_linear    = seq(from = 0, to = 2000, by = 10)
    , alpha_quadratic = seq(from = 0, to = 2000, by = 5)
    , schemes         = c("linear", "quadratic")
    , form            = "rows"
    , lambda          = NULL
    , log_scale       = TRUE
    , as_dt           = FALSE
    )
mhp_interpret(
    x   = f
  , top = 3
  )
#> Fully modified HP filter
#> 
#> The fully modified HP filter chose lambda = 1,448 (the modified HP filter
#> alone: 683) and quadratic end-point weights on 22 values at each end, which
#> cut the end-point loss of the HP filter with the same lambda by 65%.
#> 
#> Findings:
#> - Cycle (log scale, about the share of the trend): standard deviation 0.0162,
#>   first-order autocorrelation 0.88, range 0.0871.
#> - Largest positive cycle values (booms) at 2000Q2, 1978Q4, 1999Q4; largest
#>   negative ones (recessions) at 2009Q2, 2009Q3, 2009Q1.
#> - Linear weights: loss 0.855 (alpha = 240); quadratic weights: loss 0.616
#>   (alpha = 15); the quadratic scheme was kept.
#> 
#> Recommendations:
#> - Compare the filters on this series with fmhp_compare() (moments of the
#>   cycles and end-point loss).
#> 
#> Cautions:
#> - The last values of the cycle are revised as new observations arrive: read
#>   them with care.
#> 
#> Sources:
#> - Hanif, M. N., Iqbal, J. and Choudhary, M. A. (2017). Fully modified HP
#>   filter. SBP Working Paper Series, No. 88, State Bank of Pakistan.
#> - Hanif, M. N. (2017). Fully modified HP filter function (MATLAB code).
#>   MATLAB Central File Exchange, 63198.
#> - Bloechl, A. (2014). Reducing the excess variability of the Hodrick-Prescott
#>   filter by flexible penalization. Discussion Paper 2014-1, University of
#>   Munich.
#> - McDermott, C. J. (1997). An automatic method for choosing the smoothing
#>   parameter in the HP filter. Unpublished, International Monetary Fund,
#>   Washington DC.
#> - Hodrick, R. J. and Prescott, E. C. (1997). Postwar U.S. business cycles: an
#>   empirical investigation. Journal of Money, Credit and Banking, 29(1), 1-16.
#>   doi:10.2307/2953682
```
