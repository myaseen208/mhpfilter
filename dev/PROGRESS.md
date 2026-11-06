# mhpfilter: progress (read this first)

Prefix `mhp`. Standard: SKILL.qmd. Released version audited in dev/AUDIT.md; legacy snapshots in dev/legacy/.

## Phase

| phase | state | note |
|---|---|---|
| P0 audit and snapshots | done | 0.1.0 (CRAN); dev/AUDIT.md, dev/legacy |
| P1 repository and scaffold | done | common layer from baystability, prefix mhp |
| P2 cores (computation) | done | src/mhp_core.cpp: pentadiagonal LDL', O(T) per lambda |
| P3 interface and results | done | S3 generics with .default; data.table results kept (as_dt = TRUE) |
| P4 help pages | done | LaTeX maths, every argument in the examples |
| P5 tests | done | legacy, dense, MATLAB ports, published values |
| P6 vignettes (HTML + PDF) | done | seven vignettes |
| P7 README, NEWS, site | done | README.qmd, NEWS.qmd, _pkgdown.yml |
| P8 checks and release | done | tools/check-all.R |

## Functions -> cores -> tests -> vignette

| exported function | core(s) | test file | vignette | state |
|---|---|---|---|---|
| hp_filter | .mhp_trend_cpp | test-filters, test-legacy | 01, 07 | done |
| mhp_filter | .mhp_crit_cpp, .mhp_gcv_exact_cpp | test-core, test-filters, test-legacy | 01, 02, 03 | done |
| fmhp_filter | .mhp_scheme_loss_cpp, .mhp_crit_cpp | test-core, test-filters | 01, 02, 04 | done |
| mhp_batch, fmhp_batch | as the filters | test-batch, test-legacy | 03, 04 | done |
| mhp_compare, batch_compare, fmhp_compare | as the filters | test-batch, test-legacy | 01, 03 | done |
| hp_schedule | .mhp_schedule_cpp | test-diagnostics | 02, 04 | done |
| hp_weights, hp_gain, hp_loss | .mhp_weights_cpp, .mhp_gain_cpp, .mhp_loss_cpp | test-core, test-diagnostics, test-reference | 02, 04 | done |
| mhp_designs, mhp_simulate, fmhp_simulation | as the filters | test-simulate | 03, 05 | done |
| mhp_interpret | | test-interpret | 01 | done |
| mhp_cite | | test-references | 01 | done |
| mhp_app, mhp_app_object | | test-app | 06 | done |
| plot_comparison, plot_batch, autoplot.mhp | | test-filters, test-batch | 01, 03 | done |

## Corrections of the released version

| function | result | old value | new value | reason | in NEWS |
|---|---|---|---|---|---|
| mhp_filter (help page) | GCV factor | 1 + 2/(T lambda) | 1 + 2T/lambda | the code and eq. 6 use the latter | yes |
| plot_batch | linewidth, highlight | warning; highlight ignored with facet | corrected | bug | yes |
| mhp_batch, batch_compare | missing values | error | leading and trailing NA dropped | series of different lengths | yes |

No numerical result of 0.1.0 changes (test-legacy, 1e-8).

## Last test summary

See the output of tools/check-all.R at release.

## Decisions specific to this package

* The fully modified HP filter offers both forms of the end-point weights; `form = "rows"` (the default) reproduces
  the MATLAB code, `"differences"` is the penalised least-squares form of eq. 9.
* Defaults of `fmhp_filter()` follow the MATLAB code (k = 15 linear, round(min(22, T/2 - 5)) quadratic), not the
  k = 20 of the text.
* MATLAB rounding (halves up) for the centre of the loss and the quadratic k: needed to reproduce Table 3a/3b.
* The annual "our scheme" column of Table 3a is the quadratic scheme with k = 20, alpha = 10 (reference.csv);
  the quarterly column is not reproduced and is annotated.
* The simulation study ships precomputed (tools/simulation-study.R, 20 replicates); the paper's 100% claim is not
  confirmed and the vignette reports what the runs give.
* `as_dt = TRUE` remains the default of the filters, as in 0.1.0.
