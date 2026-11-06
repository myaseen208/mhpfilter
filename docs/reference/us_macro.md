# United States real GDP, consumption and investment, 1959 to 2009 (quarterly)

Quarterly real gross domestic product, real personal consumption
expenditures and real gross private domestic investment of the United
States from the first quarter of 1959 to the third quarter of 2009: the
three series (income, consumption and investment) of the empirical
analysis of Hanif, Iqbal and Choudhary (2017), for one country.

## Usage

``` r
us_macro
```

## Format

A data frame with 203 rows and 5 columns:

- `year`:

  year, 1959 to 2009.

- `quarter`:

  quarter, 1 to 4.

- `gdp`:

  real gross domestic product (billions of chained 2005 US dollars,
  seasonally adjusted annual rate).

- `consumption`:

  real personal consumption expenditures (same units).

- `investment`:

  real gross private domestic investment (same units).

## Source

Federal Reserve Economic Data (FRED), Federal Reserve Bank of St. Louis,
accessed 15 December 2009, as compiled in the data set `macrodata` of
the Python library statsmodels (public domain).

## References

Hanif, M. N., Iqbal, J. and Choudhary, M. A. (2017). Fully modified HP
filter. *SBP Working Paper Series*, No. 88, State Bank of Pakistan.

## Examples

``` r
data(us_macro)
head(us_macro)
#>   year quarter      gdp consumption investment
#> 1 1959       1 2710.349      1707.4    286.898
#> 2 1959       2 2778.801      1733.7    310.859
#> 3 1959       3 2775.488      1751.8    289.226
#> 4 1959       4 2785.204      1753.7    299.356
#> 5 1960       1 2847.699      1770.5    331.722
#> 6 1960       2 2834.390      1792.9    298.152
gdp <- ts(data = us_macro$gdp, start = c(1959, 1), frequency = 4)
plot(log(gdp))
```
