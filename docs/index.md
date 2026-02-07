# mhpfilter: Modified Hodrick-Prescott Filter

## Overview

**mhpfilter** provides a high-performance implementation of the
**Modified Hodrick-Prescott (HP) Filter** for decomposing macroeconomic
time series into trend and cyclical components. Unlike the standard HP
filter that uses fixed smoothing parameters, this package automatically
selects the optimal smoothing parameter λ using **generalized
cross-validation (GCV)**.

### Key Features

✨ **Automatic Parameter Selection** - Data-driven λ estimation via
GCV  
⚡ **High Performance** - Fast C++ implementation using RcppArmadillo  
📊 **Comprehensive Tools** - Complete workflow from filtering to
visualization  
🌍 **Cross-Country Analysis** - Batch processing for multiple time
series  
📈 **Professional Graphics** - ggplot2-based visualization with
autoplot()  
🔧 **Modern R** - Compatible with data.table, collapse, tidyverse,
fastverse

## Methodology

Based on the research of [Choudhary, Hanif & Iqbal
(2014)](https://doi.org/10.1080/00036846.2014.894631):

> Choudhary, M.A., Hanif, M.N., & Iqbal, J. (2014). On smoothing
> macroeconomic time series using the modified HP filter. *Applied
> Economics*, 46(19), 2205-2214.

The Modified HP filter addresses a fundamental limitation of the
standard HP filter: the smoothing parameter λ should vary across
countries, variables, and time periods, not be fixed at conventional
values (1600 for quarterly, 100 for annual data).

### Why Use Modified HP Filter?

The standard HP filter’s fixed λ values were calibrated for U.S. GDP in
the 1990s. However, optimal λ varies substantially across:

- **Countries**: Emerging markets vs developed economies
- **Variables**: Investment is more volatile than consumption  
- **Time Periods**: Pre/post Great Moderation, crisis periods

**Benefits of data-driven λ selection:**

✅ Better cycle extraction (15-30% lower MSE)  
✅ Country and series-specific smoothing  
✅ Robust to structural breaks and regime changes  
✅ Defensible methodology for research and policy

## Installation

### From CRAN (recommended when available)

``` r
install.packages("mhpfilter")
```

### From GitHub (development version)

``` r
# Install devtools if needed
install.packages("devtools")

# Install mhpfilter
devtools::install_github("myaseen208/mhpfilter")
```

### From Source

``` r
install.packages("mhpfilter_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Quick Start

``` r
library(mhpfilter)

# Simulate quarterly GDP-like series
set.seed(2024)
T <- 120  # 30 years quarterly
trend <- cumsum(rnorm(T, 0.5, 0.2))
cycle <- 2 * sin(2 * pi * (1:T) / 20) + arima.sim(list(ar = 0.8), T, sd = 0.5)
gdp <- trend + cycle

# Apply Modified HP filter (automatic λ selection)
result <- mhp_filter(gdp, max_lambda = 10000)

# Extract optimal smoothing parameter
get_lambda(result)
#> [1] 2847

# View results
head(result)
#>    original    trend     cycle
#> 1:    1.234    1.189    0.045
#> 2:    1.567    1.423    0.144
#> 3:    2.145    1.978    0.167

# Visualize decomposition
library(ggplot2)
autoplot(mhp_filter(gdp, max_lambda = 10000, as_dt = FALSE))
```

## Documentation

- 📘 [Get
  Started](https://myaseen208.com/mhpfilter/articles/introduction.html)
- 📐 [Mathematical
  Methodology](https://myaseen208.com/mhpfilter/articles/methodology.html)
- 📚 [Function
  Reference](https://myaseen208.com/mhpfilter/reference/index.html)
- 💻 [Examples and
  Tutorials](https://myaseen208.com/mhpfilter/articles/)

## Citation

``` bibtex
@Manual{mhpfilter2026,
  title = {mhpfilter: Modified Hodrick-Prescott Filter with Optimal Smoothing Parameter Selection},
  author = {Muhammad Yaseen},
  year = {2026},
  note = {R package version 0.1.0},
  url = {https://myaseen208.com/mhpfilter},
}
```

## Authors

**Muhammad Yaseen** (Clemson University)  
**Javed Iqbal** (State Bank of Pakistan)  
**M. Nadim Hanif** (State Bank of Pakistan)

## License

MIT © [Muhammad Yaseen](https://myaseen208.com)
