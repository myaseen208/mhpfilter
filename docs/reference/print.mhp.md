# Print Method for mhp Objects

Print a summary of Modified HP filter results.

## Usage

``` r
# S3 method for class 'mhp'
print(x, ...)
```

## Arguments

- x:

  An object of class "mhp" from `mhp_filter(as_dt = FALSE)`.

- ...:

  Additional arguments (ignored).

## Value

Invisibly returns the input object.

## Details

Prints a formatted summary including: - Number of observations - Optimal
lambda value - GCV criterion value - Cycle statistics (mean, SD, AR1)

Suitable for quick inspection of filter results.

## Examples

``` r
set.seed(42)
y <- cumsum(rnorm(100)) + sin((1:100) * pi / 20)
result <- mhp_filter(y, max_lambda = 10000, as_dt = FALSE)
print(result)
#> Modified HP Filter Result
#> ============================== 
#> Observations: 100 
#> Optimal lambda: 592 
#> GCV: 2.317e+00 
#> 
#> Cycle Statistics:
#>   Mean: -3.858e-14 
#>   SD:   1.323 
#>   AR1:  0.7203 
#>   Range: 6.537 
```
