# Cases timed by tools/bench.R: calls that run in version 0.1.0 and in this version, on simulated series (the
# released version has no data set). Each case is a quoted call evaluated in the global environment.
set.seed(2017)
x200 <- cumsum(rnorm(200, 0.5)) + 100
x500 <- cumsum(rnorm(500, 0.5)) + 100
X10  <- sapply(seq_len(10), function(j) cumsum(rnorm(200, 0.5)) + 100)
list(
    hp_T200      = quote(hp_filter(x = x200, lambda = 1600, as_dt = FALSE))
  , mhp_T200     = quote(mhp_filter(x = x200, max_lambda = 3000, as_dt = FALSE))
  , mhp_T200_big = quote(mhp_filter(x = x200, max_lambda = 100000, as_dt = FALSE))
  , mhp_T500     = quote(mhp_filter(x = x500, max_lambda = 3000, as_dt = FALSE))
  , batch_10     = quote(mhp_batch(X = X10, max_lambda = 3000))
  )
