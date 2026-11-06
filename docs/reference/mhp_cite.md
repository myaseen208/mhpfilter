# References of a result or of the package

`mhp_cite()` gives the references of the methods behind a result of the
package, or the references of every method in the package, as text or
BibTeX. `citation("mhpfilter")` cites the package and its main methods.

## Usage

``` r
mhp_cite(x = NULL, style = c("text", "bibtex"))

# Default S3 method
mhp_cite(x = NULL, style = c("text", "bibtex"))
```

## Arguments

- x:

  a result of a function of the package (with `as_dt = FALSE` for the
  filters), or `NULL` (default) for every reference of the package.

- style:

  `"text"` (default) or `"bibtex"`.

## Value

An object of class `bibentry` (printed in the chosen style).

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
mhp_cite(
    x     = NULL
  , style = "text"
  )
#> Hodrick, R. J. and Prescott, E. C. (1997). Postwar U.S. business cycles: an empirical investigation. Journal of Money, Credit and Banking, 29(1), 1-16. doi:10.2307/2953682
#> 
#> Ravn, M. O. and Uhlig, H. (2002). On adjusting the Hodrick-Prescott filter for the frequency of observations. The Review of Economics and Statistics, 84(2), 371-376. doi:10.1162/003465302317411604
#> 
#> McDermott, C. J. (1997). An automatic method for choosing the smoothing parameter in the HP filter. Unpublished, International Monetary Fund, Washington DC.
#> 
#> Coe, D. T. and McDermott, C. J. (1997). Does the gap model work in Asia? IMF Staff Papers, 44(1), 59-80. doi:10.2307/3867497
#> 
#> Choudhary, M. A., Hanif, M. N. and Iqbal, J. (2014). On smoothing macroeconomic time series using the modified HP filter. Applied Economics, 46(19), 2205-2214. doi:10.1080/00036846.2014.894631
#> 
#> Choudhary, M. A., Hanif, M. N. and Iqbal, J. (2014). modifiedhpfilternew.m: MATLAB code of the modified HP filter of Choudhary, Hanif and Iqbal (2014).
#> 
#> Craven, P. and Wahba, G. (1979). Smoothing noisy data with spline functions. Numerische Mathematik, 31, 377-403. doi:10.1007/BF01404567
#> 
#> Silverman, B. W. (1984). A fast and efficient cross-validation method for smoothing parameter choice in spline regression. Journal of the American Statistical Association, 79(387), 584-589. doi:10.1080/01621459.1984.10478084
#> 
#> Bloechl, A. (2014). Reducing the excess variability of the Hodrick-Prescott filter by flexible penalization. Discussion Paper 2014-1, University of Munich.
#> 
#> Hanif, M. N., Iqbal, J. and Choudhary, M. A. (2017). Fully modified HP filter. SBP Working Paper Series, No. 88, State Bank of Pakistan.
#> 
#> Hanif, M. N. (2017). Fully modified HP filter function (MATLAB code). MATLAB Central File Exchange, 63198.
#> 
#> Harvey, A. C. and Jaeger, A. (1993). Detrending, stylized facts and the business cycle. Journal of Applied Econometrics, 8(3), 231-247. doi:10.1002/jae.3950080302
#> 
#> Guay, A. and St-Amant, P. (2005). Do the Hodrick-Prescott and Baxter-King filters provide a good approximation of business cycles? Annales d'Economie et de Statistique, 77, 133-155. doi:10.2307/20079119 
data(us_macro)
f <-
  mhp_filter(
      x          = log(us_macro$gdp)
    , max_lambda = 100000
    , as_dt      = FALSE
    , log_scale  = FALSE
    , criterion  = "approximate"
    )
mhp_cite(
    x     = f
  , style = "bibtex"
  )
#> @Article{HodrickPrescott1997,
#>   title = {Postwar {U.S.} Business Cycles: An Empirical Investigation},
#>   author = {Robert J. Hodrick and Edward C. Prescott},
#>   journal = {Journal of Money, Credit and Banking},
#>   year = {1997},
#>   volume = {29},
#>   number = {1},
#>   pages = {1--16},
#>   doi = {10.2307/2953682},
#> }
#> 
#> @Unpublished{McDermott1997,
#>   title = {An Automatic Method for Choosing the Smoothing Parameter in the {HP} Filter},
#>   author = {C. John McDermott},
#>   note = {International Monetary Fund, Washington DC},
#>   year = {1997},
#> }
#> 
#> @Article{CoeMcDermott1997,
#>   title = {Does the Gap Model Work in {Asia}?},
#>   author = {David T. Coe and C. John McDermott},
#>   journal = {IMF Staff Papers},
#>   year = {1997},
#>   volume = {44},
#>   number = {1},
#>   pages = {59--80},
#>   doi = {10.2307/3867497},
#> }
#> 
#> @Article{Choudhary2014,
#>   title = {On Smoothing Macroeconomic Time Series Using the Modified {HP} Filter},
#>   author = {M. Ali Choudhary and M. Nadim Hanif and Javed Iqbal},
#>   journal = {Applied Economics},
#>   year = {2014},
#>   volume = {46},
#>   number = {19},
#>   pages = {2205--2214},
#>   doi = {10.1080/00036846.2014.894631},
#> }
#> 
#> @Misc{Choudhary2014code,
#>   title = {{modifiedhpfilternew.m}: {MATLAB} Code of the Modified {HP} Filter},
#>   author = {M. Ali Choudhary and M. Nadim Hanif and Javed Iqbal},
#>   year = {2014},
#> }
```
