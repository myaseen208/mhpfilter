# Print a filter result

[`print()`](https://rdrr.io/r/base/print.html) of a result of
[`hp_filter()`](https://myaseen208.com/mhpfilter/reference/hp_filter.md),
[`mhp_filter()`](https://myaseen208.com/mhpfilter/reference/mhp_filter.md)
or
[`fmhp_filter()`](https://myaseen208.com/mhpfilter/reference/fmhp_filter.md)
(with `as_dt = FALSE`) shows the smoothing parameter, the criterion and
summaries of the cycle, as in version 0.1.0, and the sources of the
method.

## Usage

``` r
# S3 method for class 'mhp'
print(x, ...)
```

## Arguments

- x:

  a result of
  [`hp_filter()`](https://myaseen208.com/mhpfilter/reference/hp_filter.md),
  [`mhp_filter()`](https://myaseen208.com/mhpfilter/reference/mhp_filter.md)
  or
  [`fmhp_filter()`](https://myaseen208.com/mhpfilter/reference/fmhp_filter.md)
  with `as_dt = FALSE`.

- ...:

  further arguments (`digits`).

## Value

`x`, invisibly.

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
data(us_macro)
m <-
  mhp_filter(
      x          = log(us_macro$gdp)
    , max_lambda = 100000
    , as_dt      = FALSE
    , log_scale  = FALSE
    , criterion  = c("approximate", "exact")
    )
print(x = m)
#> Modified HP Filter Result
#> ==============================
#> Observations: 203 
#> Optimal lambda: 683 
#> GCV: 2.82e-04 
#> 
#> Cycle Statistics:
#>   Mean: -1.634e-13 
#>   SD:   0.01333 
#>   AR1:  0.828 
#>   Range: 0.07802 
#> 
#> Sources: Hodrick and Prescott (1997); McDermott (1997); Coe and McDermott
#> (1997); Choudhary et al. (2014); Choudhary et al. (2014, MATLAB code). Full
#> references: mhp_cite(x).
```
