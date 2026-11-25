# Run several GBLUP scenarios from `cv_scenarios()`

Run several GBLUP scenarios from `cv_scenarios()`

## Usage

``` r
asr_scenarios(
  sc,
  ginv,
  scenarios = c("st_sp", "st_mp", "mt_sp", "mt_mp"),
  dir = ".",
  tag = deparse(substitute(ginv)),
  ...
)
```

## Arguments

- sc:

  Output of `mintyr::cv_scenarios()`.

- ginv:

  Inverse genomic relationship matrix.

- scenarios:

  Scenarios to run. `"mt_bt"` (one trait as a different trait per group,
  two groups) is fitted with `residual_cov = "diag"` and without rp,
  since no individual has records in both groups.

- dir:

  Folder for the workbooks `<scenario>_<tag>.xlsx`; `NULL` writes
  nothing.

- tag:

  Label used in the file names; defaults to the name of `ginv`.

- ...:

  Further arguments for
  [`asr_gblup()`](https://tony2015116.github.io/mintyrasr/reference/asr_gblup.md)
  (e.g. `fixed`, `inbreeding`).

## Value

A named list of
[`asr_gblup()`](https://tony2015116.github.io/mintyrasr/reference/asr_gblup.md)
results.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- asr_scenarios(sc, ginv, fixed = ~ hhys, inbreeding = inbreeding)
} # }
```
