### Name: print.mhp
### Title: Print Method for mhp Objects
### Aliases: print.mhp

### ** Examples

set.seed(42)
y <- cumsum(rnorm(100)) + sin((1:100) * pi / 20)
result <- mhp_filter(y, max_lambda = 10000, as_dt = FALSE)
print(result)




