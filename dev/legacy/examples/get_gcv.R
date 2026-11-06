### Name: get_gcv
### Title: Extract GCV Value
### Aliases: get_gcv

### ** Examples

set.seed(123)
result <- mhp_filter(cumsum(rnorm(100)), max_lambda = 10000)
get_gcv(result)

# With mhp object
result_obj <- mhp_filter(cumsum(rnorm(100)), max_lambda = 10000, as_dt = FALSE)
get_gcv(result_obj)




