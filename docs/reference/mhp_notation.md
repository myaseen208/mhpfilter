# Notation of mhpfilter

The symbols used in the help pages and vignettes of mhpfilter, with
their meaning and the argument or result component that holds them. The
notation follows one rule throughout: scalars are upper-case letters in
normal face (\\P\\, \\R\\), indices lower-case italic (\\i, j\\),
vectors lower-case bold (\\\mathbf{y}\\, \\\boldsymbol{\beta}\\) and
matrices upper-case bold (\\\mathbf{X}\\, \\\boldsymbol{\Sigma}\\).

## Usage

``` r
mhp_notation(format = c("table", "path", "tex"))
```

## Arguments

- format:

  `"table"` (default: a `data.table` with columns `symbol` (LaTeX),
  `text` (plain text), `kind` (scalar, index, vector, matrix, set,
  function, operator), `meaning` and `argument`), `"path"` (the path of
  `notation.csv`) or `"tex"` (the path of `Notations.tex`, the notation
  appendix used in the articles).

## Value

See `format`.

## See also

Other mhpfilter utilities:
[`mhp_reference()`](https://myaseen208.com/mhpfilter/reference/mhp_reference.md),
[`mhp_save()`](https://myaseen208.com/mhpfilter/reference/mhp_save.md),
[`mhp_theme()`](https://myaseen208.com/mhpfilter/reference/mhp_theme.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
mhp_notation(
    format = c("table", "path", "tex")[1]
  )
#>                                    symbol            text   kind
#>                                    <char>          <char> <char>
#>  1:                                     T               T scalar
#>  2:                                  t, j             t j  index
#>  3:                           \\mathbf{x}               x vector
#>  4:                           \\mathbf{g}               g vector
#>  5:                           \\mathbf{c}               c vector
#>  6:                           \\mathbf{K}               K matrix
#>  7:                           \\mathbf{A}               A matrix
#>  8:                           \\mathbf{I}               I matrix
#>  9:                           \\mathbf{B}               B matrix
#> 10:                              \\lambda          lambda scalar
#> 11:                            \\lambda_t        lambda_t vector
#> 12:                               \\alpha           alpha scalar
#> 13:                                     k               k scalar
#> 14:                               \\omega           omega scalar
#> 15:                          P_t(\\omega)      P_t(omega) scalar
#> 16:                               l(t), L          l(t) L scalar
#> 17:                               \\theta           theta scalar
#> 18:                                \\mu_i            mu_i scalar
#> 19:                      \\phi_1, \\phi_2     phi_1 phi_2 scalar
#> 20: \\sigma_\\varepsilon, \\sigma_\\delta sigma_e sigma_d scalar
#>                                    symbol            text   kind
#>                                    <char>          <char> <char>
#>                                                              meaning
#>                                                               <char>
#>  1:                                           number of observations
#>  2:     time and the observations weighted in the estimate at time t
#>  3:                 the series (its logarithm when log_scale = TRUE)
#>  4:                                                            trend
#>  5:                                                     cycle, x - g
#>  6:                         (T - 2) x T matrix of second differences
#>  7:                                               penalty matrix K'K
#>  8:                                                  identity matrix
#>  9:  smoother: g = B x; row t holds the weights of the estimate at t
#> 10:                                              smoothing parameter
#> 11: smoothing parameter of observation t (or of second difference t)
#> 12:                                                 end-point weight
#> 13:            number of end values with raised smoothing parameters
#> 14:                                              frequency (radians)
#> 15:                                        gain of the estimate at t
#> 16:               end-point loss of the estimate at t, and its total
#> 17:                     spacing of the frequencies in the loss (0.1)
#> 18:                                     eigenvalues of A (exact GCV)
#> 19:                          AR(2) coefficients of a simulated cycle
#> 20:                standard deviations of the trend and cycle shocks
#>                                                              meaning
#>                                                               <char>
#>         argument
#>           <char>
#>  1:            x
#>  2:             
#>  3:            x
#>  4:             
#>  5:             
#>  6:             
#>  7:             
#>  8:             
#>  9:             
#> 10:       lambda
#> 11:     lambda_t
#> 12: alpha_linear
#> 13:     k_linear
#> 14:        theta
#> 15:             
#> 16:             
#> 17:       weight
#> 18:             
#> 19:           ar
#> 20:     sd_ratio
#>         argument
#>           <char>
```
