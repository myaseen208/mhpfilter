# The power comparison of the vignette "Simulation": all 30 designs of Hanif, Iqbal and Choudhary (2017), quarterly
# series (T = 200) and series averaged to annual data (T = 50), both forms of the end-point weights.
# Writes inst/extdata/simulation_rmse.csv and inst/extdata/simulation_power.csv. Run from the package folder:
#   Rscript tools/simulation-study.R [n_rep]
args  <- commandArgs(trailingOnly = TRUE)
n_rep <- if (length(args)) as.integer(args[1L]) else 50L
pkgload::load_all(".", quiet = TRUE)
runs <- data.table::CJ(aggregate = c("none", "average"), trend = c("linear", "nonlinear"),
                       form = c("rows", "differences"), sorted = FALSE)
res <- lapply(seq_len(nrow(runs)), function(i) {
  r <- runs[i]
  t0 <- proc.time()[["elapsed"]]
  s <- fmhp_simulation(designs = mhp_designs(), n_rep = n_rep, T = 200L, trend = r$trend, aggregate = r$aggregate,
                       form = r$form, seed = 2017L + i)
  message(sprintf("%s %s %s: %.0f s", r$aggregate, r$trend, r$form, proc.time()[["elapsed"]] - t0))
  list(rmse = cbind(r, s$rmse), power = cbind(r, s$power))
})
data.table::fwrite(data.table::rbindlist(lapply(res, `[[`, "rmse")), "inst/extdata/simulation_rmse.csv")
data.table::fwrite(data.table::rbindlist(lapply(res, `[[`, "power")), "inst/extdata/simulation_power.csv")
