# The 30 simulation designs of Hanif, Iqbal and Choudhary (2017)

`mhp_designs()` lists the 30 combinations of Table 1 of the paper: five
ratios \\\sigma\_\varepsilon / \sigma\_\delta\\ of the standard
deviations of the trend and cycle shocks (10, 5, 2, 1, 0.5) crossed with
six AR(2) cycles \\(\phi_1, \phi_2)\\: (0.9, 0.01), (1.2, -0.25), (1.2,
-0.4), (1.2, -0.55), (1.2, -0.75) and (0.8, 0). The designs of
Choudhary, Hanif and Iqbal (2014, Table 1) have the same cycles and the
ratios 10, 5, 1, 0.5 and 0.01.

## Usage

``` r
mhp_designs(paper = c("2017", "2014"))
```

## Arguments

- paper:

  `"2017"` (default, Hanif et al.) or `"2014"` (Choudhary et al.).

## Value

A `data.table` with the columns `model` (1 to 30, as in the paper),
`sd_ratio`, `ar1`, `ar2`.

## References

Choudhary, M. A., Hanif, M. N. and Iqbal, J. (2014). On smoothing
macroeconomic time series using the modified HP filter. *Applied
Economics*, **46**(19), 2205–2214.
[doi:10.1080/00036846.2014.894631](https://doi.org/10.1080/00036846.2014.894631)

Hanif, M. N., Iqbal, J. and Choudhary, M. A. (2017). Fully modified HP
filter. *SBP Working Paper Series*, No. 88, State Bank of Pakistan.

## See also

[`mhp_simulate`](https://myaseen208.com/mhpfilter/reference/mhp_simulate.md),
[`fmhp_simulation`](https://myaseen208.com/mhpfilter/reference/fmhp_simulation.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
mhp_designs(
    paper = "2017"
  )
#>     model sd_ratio   ar1   ar2
#>     <int>    <num> <num> <num>
#>  1:     1     10.0   0.9  0.01
#>  2:     2     10.0   1.2 -0.25
#>  3:     3     10.0   1.2 -0.40
#>  4:     4     10.0   1.2 -0.55
#>  5:     5     10.0   1.2 -0.75
#>  6:     6      5.0   0.9  0.01
#>  7:     7      5.0   1.2 -0.25
#>  8:     8      5.0   1.2 -0.40
#>  9:     9      5.0   1.2 -0.55
#> 10:    10      5.0   1.2 -0.75
#> 11:    11      2.0   0.9  0.01
#> 12:    12      2.0   1.2 -0.25
#> 13:    13      2.0   1.2 -0.40
#> 14:    14      2.0   1.2 -0.55
#> 15:    15      2.0   1.2 -0.75
#> 16:    16      1.0   0.9  0.01
#> 17:    17      1.0   1.2 -0.25
#> 18:    18      1.0   1.2 -0.40
#> 19:    19      1.0   1.2 -0.55
#> 20:    20      1.0   1.2 -0.75
#> 21:    21      0.5   0.9  0.01
#> 22:    22      0.5   1.2 -0.25
#> 23:    23      0.5   1.2 -0.40
#> 24:    24      0.5   1.2 -0.55
#> 25:    25      0.5   1.2 -0.75
#> 26:    26     10.0   0.8  0.00
#> 27:    27      5.0   0.8  0.00
#> 28:    28      2.0   0.8  0.00
#> 29:    29      1.0   0.8  0.00
#> 30:    30      0.5   0.8  0.00
#>     model sd_ratio   ar1   ar2
#>     <int>    <num> <num> <num>
```
