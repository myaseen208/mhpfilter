# Published values used by the tests and vignettes

Reads `inst/extdata/reference.csv`: one row per published number, with
its source, the R call that recomputes it with this package (`call`),
the expression that extracts the number from the result (`extract`, in
terms of `x`), the printed value and the digits to compare. The same
file drives the test of published values and the vignette that
reproduces the published examples. Each distinct `call` is evaluated
once.

## Usage

``` r
mhp_reference(source = NULL, evaluate = FALSE, level = 1)
```

## Arguments

- source:

  optional regular expression; keep rows whose `source` matches.

- evaluate:

  logical; evaluate the rows and add columns `computed` and `agrees`
  (default `FALSE`).

- level:

  largest `level` evaluated (1 = fast rows, 2 = slow rows; default 1);
  rows above it get `NA`.

## Value

A `data.table` with columns `id`, `source`, `quantity`, `call`,
`extract`, `value`, `digits`, `level`, `note` (and `computed`,
`agrees`).

## See also

Other mhpfilter utilities:
[`mhp_notation()`](https://myaseen208.com/mhpfilter/reference/mhp_notation.md),
[`mhp_save()`](https://myaseen208.com/mhpfilter/reference/mhp_save.md),
[`mhp_theme()`](https://myaseen208.com/mhpfilter/reference/mhp_theme.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
ref <-
  mhp_reference(
      source   = NULL
    , evaluate = FALSE
    , level    = 1
    )
ref
#>                          id
#>                      <char>
#>  1:    hanif2017_t3a_hp_T55
#>  2:    hanif2017_t3a_hp_T54
#>  3:    hanif2017_t3a_hp_T40
#>  4:    hanif2017_t3b_hp_T72
#>  5:   hanif2017_t3b_hp_T103
#>  6:  hanif2017_t3a_ours_T55
#>  7:  hanif2017_t3a_ours_T54
#>  8:  hanif2017_t3a_ours_T50
#>  9: hanif2017_fig1_hp_range
#> 10:   hanif2017_fig1a_range
#> 11:   hanif2017_t3b_ours_US
#>                                                         source
#>                                                         <char>
#>  1:                Hanif, Iqbal and Choudhary (2017), Table 3a
#>  2:                Hanif, Iqbal and Choudhary (2017), Table 3a
#>  3:                Hanif, Iqbal and Choudhary (2017), Table 3a
#>  4:                Hanif, Iqbal and Choudhary (2017), Table 3b
#>  5:                Hanif, Iqbal and Choudhary (2017), Table 3b
#>  6:                Hanif, Iqbal and Choudhary (2017), Table 3a
#>  7:                Hanif, Iqbal and Choudhary (2017), Table 3a
#>  8:                Hanif, Iqbal and Choudhary (2017), Table 3a
#>  9:  Hanif, Iqbal and Choudhary (2017), Section 5 and Figure 1
#> 10: Hanif, Iqbal and Choudhary (2017), Section 5 and Figure 1a
#> 11:                Hanif, Iqbal and Choudhary (2017), Table 3b
#>                                                                                                                  quantity
#>                                                                                                                    <char>
#>  1:                              end-point loss of the HP filter (lambda = 100) for 55 annual observations (1960 to 2014)
#>  2:                                             end-point loss of the HP filter (lambda = 100) for 54 annual observations
#>  3:                       end-point loss of the HP filter (lambda = 100) for 40 annual observations (the shortest series)
#>  4:         end-point loss of the HP filter (lambda = 1600) for 72 quarterly observations (the shortest span of 18 years)
#>  5:       end-point loss of the HP filter (lambda = 1600) for 103 quarterly observations (the largest value of the table)
#>  6: end-point loss of the HP filter (lambda = 100) with quadratic weights (k = 20; alpha = 10) for 55 annual observations
#>  7: end-point loss of the HP filter (lambda = 100) with quadratic weights (k = 20; alpha = 10) for 54 annual observations
#>  8: end-point loss of the HP filter (lambda = 100) with quadratic weights (k = 20; alpha = 10) for 50 annual observations
#>  9:                                    range of the weights of the first estimate of the HP filter (lambda = 100; T = 50)
#> 10:          range of the weights of the first estimate with quadratic weights (k = 20; alpha = 10; lambda = 100; T = 50)
#> 11:                  end-point loss of the HP filter (lambda = 1600) with the paper's weights for the US quarterly series
#>                                                                                              call
#>                                                                                            <char>
#>  1:                                                                 hp_loss(x = 55, lambda = 100)
#>  2:                                                                 hp_loss(x = 54, lambda = 100)
#>  3:                                                                 hp_loss(x = 40, lambda = 100)
#>  4:                                                                hp_loss(x = 72, lambda = 1600)
#>  5:                                                               hp_loss(x = 103, lambda = 1600)
#>  6:    hp_loss(x = 55, lambda = hp_schedule(T = 55, lambda = 100, alpha = 10, k = 20, power = 2))
#>  7:    hp_loss(x = 54, lambda = hp_schedule(T = 54, lambda = 100, alpha = 10, k = 20, power = 2))
#>  8:    hp_loss(x = 50, lambda = hp_schedule(T = 50, lambda = 100, alpha = 10, k = 20, power = 2))
#>  9:                                                              hp_weights(x = 50, lambda = 100)
#> 10: hp_weights(x = 50, lambda = hp_schedule(T = 50, lambda = 100, alpha = 10, k = 20, power = 2))
#> 11: hp_loss(x = 200, lambda = hp_schedule(T = 200, lambda = 1600, alpha = 10, k = 20, power = 2))
#>                         extract value digits level
#>                          <char> <num>  <int> <int>
#>  1:                     x$total 1.735      3     1
#>  2:                     x$total 1.733      3     1
#>  3:                     x$total 1.731      3     1
#>  4:                     x$total 1.757      3     1
#>  5:                     x$total 1.770      3     1
#>  6:                     x$total 0.626      3     1
#>  7:                     x$total 0.627      3     1
#>  8:                     x$total 0.623      3     1
#>  9: diff(range(x$weights[1, ])) 0.385      3     1
#> 10: diff(range(x$weights[1, ])) 0.182      3     1
#> 11:                     x$total 0.797      3     1
#>                                                                                                                                                                                                              note
#>                                                                                                                                                                                                            <char>
#>  1:                                                                                                                                                                                                              
#>  2:                                                                                                                                                                                                              
#>  3:                                                                                                                                                                                                              
#>  4:                                                                                                                                                                                                              
#>  5:                                                                                                                                                                                                              
#>  6:                                                                                                                                                                                                              
#>  7:                                                                                                                                                                                                              
#>  8:                                                                                                                                                                                                              
#>  9:                                                                                                                                                                                                              
#> 10:                                                                                                                                                                                                              
#> 11: not reproduced: the quarterly column (0.770 to 0.827) is not given by the weights of the annual column (quadratic, k = 20, alpha = 10: about 1.00); the weights and the length of the series are not reported
# \donttest{
mhp_reference(
    source   = "."
  , evaluate = TRUE
  , level    = 1
  )
#>                          id
#>                      <char>
#>  1:    hanif2017_t3a_hp_T55
#>  2:    hanif2017_t3a_hp_T54
#>  3:    hanif2017_t3a_hp_T40
#>  4:    hanif2017_t3b_hp_T72
#>  5:   hanif2017_t3b_hp_T103
#>  6:  hanif2017_t3a_ours_T55
#>  7:  hanif2017_t3a_ours_T54
#>  8:  hanif2017_t3a_ours_T50
#>  9: hanif2017_fig1_hp_range
#> 10:   hanif2017_fig1a_range
#> 11:   hanif2017_t3b_ours_US
#>                                                         source
#>                                                         <char>
#>  1:                Hanif, Iqbal and Choudhary (2017), Table 3a
#>  2:                Hanif, Iqbal and Choudhary (2017), Table 3a
#>  3:                Hanif, Iqbal and Choudhary (2017), Table 3a
#>  4:                Hanif, Iqbal and Choudhary (2017), Table 3b
#>  5:                Hanif, Iqbal and Choudhary (2017), Table 3b
#>  6:                Hanif, Iqbal and Choudhary (2017), Table 3a
#>  7:                Hanif, Iqbal and Choudhary (2017), Table 3a
#>  8:                Hanif, Iqbal and Choudhary (2017), Table 3a
#>  9:  Hanif, Iqbal and Choudhary (2017), Section 5 and Figure 1
#> 10: Hanif, Iqbal and Choudhary (2017), Section 5 and Figure 1a
#> 11:                Hanif, Iqbal and Choudhary (2017), Table 3b
#>                                                                                                                  quantity
#>                                                                                                                    <char>
#>  1:                              end-point loss of the HP filter (lambda = 100) for 55 annual observations (1960 to 2014)
#>  2:                                             end-point loss of the HP filter (lambda = 100) for 54 annual observations
#>  3:                       end-point loss of the HP filter (lambda = 100) for 40 annual observations (the shortest series)
#>  4:         end-point loss of the HP filter (lambda = 1600) for 72 quarterly observations (the shortest span of 18 years)
#>  5:       end-point loss of the HP filter (lambda = 1600) for 103 quarterly observations (the largest value of the table)
#>  6: end-point loss of the HP filter (lambda = 100) with quadratic weights (k = 20; alpha = 10) for 55 annual observations
#>  7: end-point loss of the HP filter (lambda = 100) with quadratic weights (k = 20; alpha = 10) for 54 annual observations
#>  8: end-point loss of the HP filter (lambda = 100) with quadratic weights (k = 20; alpha = 10) for 50 annual observations
#>  9:                                    range of the weights of the first estimate of the HP filter (lambda = 100; T = 50)
#> 10:          range of the weights of the first estimate with quadratic weights (k = 20; alpha = 10; lambda = 100; T = 50)
#> 11:                  end-point loss of the HP filter (lambda = 1600) with the paper's weights for the US quarterly series
#>                                                                                              call
#>                                                                                            <char>
#>  1:                                                                 hp_loss(x = 55, lambda = 100)
#>  2:                                                                 hp_loss(x = 54, lambda = 100)
#>  3:                                                                 hp_loss(x = 40, lambda = 100)
#>  4:                                                                hp_loss(x = 72, lambda = 1600)
#>  5:                                                               hp_loss(x = 103, lambda = 1600)
#>  6:    hp_loss(x = 55, lambda = hp_schedule(T = 55, lambda = 100, alpha = 10, k = 20, power = 2))
#>  7:    hp_loss(x = 54, lambda = hp_schedule(T = 54, lambda = 100, alpha = 10, k = 20, power = 2))
#>  8:    hp_loss(x = 50, lambda = hp_schedule(T = 50, lambda = 100, alpha = 10, k = 20, power = 2))
#>  9:                                                              hp_weights(x = 50, lambda = 100)
#> 10: hp_weights(x = 50, lambda = hp_schedule(T = 50, lambda = 100, alpha = 10, k = 20, power = 2))
#> 11: hp_loss(x = 200, lambda = hp_schedule(T = 200, lambda = 1600, alpha = 10, k = 20, power = 2))
#>                         extract value digits level
#>                          <char> <num>  <int> <int>
#>  1:                     x$total 1.735      3     1
#>  2:                     x$total 1.733      3     1
#>  3:                     x$total 1.731      3     1
#>  4:                     x$total 1.757      3     1
#>  5:                     x$total 1.770      3     1
#>  6:                     x$total 0.626      3     1
#>  7:                     x$total 0.627      3     1
#>  8:                     x$total 0.623      3     1
#>  9: diff(range(x$weights[1, ])) 0.385      3     1
#> 10: diff(range(x$weights[1, ])) 0.182      3     1
#> 11:                     x$total 0.797      3     1
#>                                                                                                                                                                                                              note
#>                                                                                                                                                                                                            <char>
#>  1:                                                                                                                                                                                                              
#>  2:                                                                                                                                                                                                              
#>  3:                                                                                                                                                                                                              
#>  4:                                                                                                                                                                                                              
#>  5:                                                                                                                                                                                                              
#>  6:                                                                                                                                                                                                              
#>  7:                                                                                                                                                                                                              
#>  8:                                                                                                                                                                                                              
#>  9:                                                                                                                                                                                                              
#> 10:                                                                                                                                                                                                              
#> 11: not reproduced: the quarterly column (0.770 to 0.827) is not given by the weights of the annual column (quadratic, k = 20, alpha = 10: about 1.00); the weights and the length of the series are not reported
#>      computed agrees
#>         <num> <lgcl>
#>  1: 1.7345013   TRUE
#>  2: 1.7333085   TRUE
#>  3: 1.7312665   TRUE
#>  4: 1.7571952   TRUE
#>  5: 1.7701498   TRUE
#>  6: 0.6256808   TRUE
#>  7: 0.6270290   TRUE
#>  8: 0.6230345   TRUE
#>  9: 0.3853298   TRUE
#> 10: 0.1817440   TRUE
#> 11: 0.9998254  FALSE
# }
```
