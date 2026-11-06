### Name: get_lambda
### Title: Extract Optimal Lambda
### Aliases: get_lambda

### ** Examples

set.seed(123)
result <- mhp_filter(cumsum(rnorm(100)), max_lambda = 10000)
get_lambda(result)

# With mhp object
result_obj <- mhp_filter(cumsum(rnorm(100)), max_lambda = 10000, as_dt = FALSE)
get_lambda(result_obj)




