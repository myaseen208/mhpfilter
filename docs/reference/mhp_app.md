# Trend and cycle analysis in a web browser

`mhp_app()` starts the application of mhpfilter: load a series (the
example
[`us_macro`](https://myaseen208.com/mhpfilter/reference/us_macro.md) or
a CSV file), filter it with the HP, modified HP or fully modified HP
filter, compare the three filters, inspect their weights, gains and
end-point loss, filter all the series of a file at once and run the
power comparison of Hanif, Iqbal and Choudhary (2017). Every result page
shows an interpretation computed from the results
([`mhp_interpret()`](https://myaseen208.com/mhpfilter/reference/mhp_interpret.md)),
the output as R prints it, formatted tables, plots and the R code that
repeats the analysis, each with a Copy button.

## Usage

``` r
mhp_app(launch.browser = interactive(), ...)

mhp_app_object()
```

## Arguments

- launch.browser:

  passed to
  [`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html): `TRUE`
  (default in an interactive session) opens the application in the web
  browser.

- ...:

  further arguments to
  [`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html), such
  as `port`.

## Value

`mhp_app()` is called for its side effect and returns the value of
[`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html) when the
application is closed; `mhp_app_object()` returns the application object
(for [`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html),
[`shiny::testServer()`](https://rdrr.io/pkg/shiny/man/testServer.html)
or deployment on a server).

## Details

The application needs the suggested packages shiny, bslib and DT. Plots
are sent to the browser as WebP images when magick is installed (smaller
than PNG; PNG otherwise) and can be downloaded as PNG, PDF or WebP. Text
is set in IBM Plex Sans and output in IBM Plex Mono when the browser can
reach Google Fonts, in the system fonts otherwise.

## See also

[`mhp_filter`](https://myaseen208.com/mhpfilter/reference/mhp_filter.md),
[`fmhp_filter`](https://myaseen208.com/mhpfilter/reference/fmhp_filter.md),
[`fmhp_compare`](https://myaseen208.com/mhpfilter/reference/fmhp_compare.md),
[`mhp_interpret`](https://myaseen208.com/mhpfilter/reference/mhp_interpret.md)

## Author

Muhammad Yaseen (<myaseen208@gmail.com>)

## Examples

``` r
if (interactive()) {
  mhp_app(
      launch.browser = TRUE
    )
}
```
