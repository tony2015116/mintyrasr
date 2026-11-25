# GBLUP parameter estimation and cross-validation with ASReml-R

`asr_gblup()` runs the whole GBLUP workflow for one scenario directly in
R, so no input files have to be exported for ASReml:

1.  fits the model to every data set in `param`: variance components,
    heritability and, for two-trait models, genetic / phenotypic
    correlations plus an LRT of rg = 0 (`idh(trait)` vs the full genetic
    structure, `boundary = FALSE`);

2.  refits the model to every training set in `cv`, with the variance
    components handled as set by `cv_vc`;

3.  predicts GEBVs of the validation individuals, computes their
    reliability `1 - PEV / ((1 + F) * Va)`, drops GEBVs with reliability
    below `min_rel`, and returns per fold the predictive ability
    `pa = r(y, GEBV)`, the accuracy `pa / sqrt(h2)` (h2 from the full
    data, Dekkers 2007) and the slope `b(y ~ GEBV)`;

4.  summarises them over folds (n, mean, sd, se) and, optionally, writes
    everything to one Excel workbook.

h2, rg and rp are built from the **names** of the variance components,
not from fixed `V` positions, so they stay correct when random effects
are added (`random = ~ ide(id) + litter`) or covariance structures
change (`us`, `idh`/`diag`, `corh`, `corgh`): h2 = Va / (all variances
of the trait), rg = Cov_a / sqrt(Va1 Va2) (the correlation itself for
`corh`/`corgh`), rp = sum of all covariances / sqrt(Vp1 Vp2). The
`vc_index` table shows which `V` each component is and `param$formulas`
the formulas used. `V` numbers follow `fit$vparameters`, the order
`vpredict()` uses. Components fixed at the boundary (bound `B`) are
reported with a warning and in `param$boundary`. `fa()` / `rr()` terms
and models with several genetic terms (e.g. maternal
[`str()`](https://rdrr.io/r/utils/str.html) models) need your own
formulas in `vpredict`.

Single-trait or two-trait models are chosen from the nested data: a
`value` column gives a single-trait model, `value1` + `value2` a
two-trait model (layout of `mintyr::w2l_nest()` / `mintyr::c2p_nest()`).

## Usage

``` r
asr_gblup(
  param,
  cv,
  ginv,
  fixed = ~1,
  random = NULL,
  id = "id",
  traits = NULL,
  genetic_cov = "us",
  random_cov = "us",
  residual_cov = "us",
  cv_vc = c("fixed", "warm", "estimate"),
  vpredict = NULL,
  inbreeding = NULL,
  min_rel = 0.05,
  neg_to_na = TRUE,
  test_rg = TRUE,
  test_rp = FALSE,
  rp = TRUE,
  na_x = c("include", "fail", "omit"),
  dense = TRUE,
  maxit = 100L,
  workspace = "30Gb",
  max_update = 0L,
  asreml_args = list(),
  on_error = c("warn", "stop"),
  keep_pred = TRUE,
  path = NULL,
  verbose = TRUE
)
```

## Arguments

- param:

  Nested table used for parameter estimation, with a list-column `data`
  (e.g. `cv_scenarios(...)$st_sp$param`).

- cv:

  Cross-validation table with list-columns `train` and `validate` (e.g.
  `cv_scenarios(...)$st_sp$cv` or the output of `mintyr::nest_cv()`).
  Each row is matched to the row of `param` that shares its outer
  columns.

- ginv:

  Inverse genomic (or pedigree) relationship matrix as used by `vm()`,
  e.g. `G.inverse(G, sparseform = TRUE)$Ginv.sparse`. It must carry
  `attr(, "INVERSE") = TRUE`, otherwise `vm()` takes it as G itself (a
  warning is given).

- fixed:

  One-sided formula of fixed effects, e.g. `~ hhys`. For two-trait
  models every term is nested in trait (`trait + trait:hhys`).

- random:

  Optional one-sided formula of further random effects, e.g. `~ ide(id)`
  (permanent environment, repeated records), `~ litter + pen`. In
  two-trait models each term gets the structure `random_cov`.

- id:

  Name of the individual ID column in the nested data.

- traits:

  Optional response column(s), one or two, when the nested data do not
  use `value` / `value1` + `value2`, e.g. `c("bfs", "nls")` for one
  trait recorded in two herds (`cv_scenarios(...)$mt_bt`).

- genetic_cov, random_cov, residual_cov:

  Two-trait models: variance structure across traits for the genomic
  term, the extra random terms and the residual: `"us"`, `"idh"`,
  `"diag"`, `"corv"`, `"corh"` or `"corgh"` (residual: `"us"`, `"idh"`,
  `"diag"`).

- cv_vc:

  How the variance components are handled in each fold: `"fixed"`
  (default, fast): fixed at the full-data estimates (`Constraint = "F"`,
  `maxit = 1`; slight information leakage because validation records
  were used to estimate them); `"warm"`: re-estimated from the training
  set, starting from the full-data values; `"estimate"`: re-estimated
  from the training set from default start values (no leakage, slowest).

- vpredict:

  Optional named list of extra `vpredict()` formulas, e.g.
  `list(rep = rep ~ (V1 + V2) / (V1 + V2 + V3))`. Check the `V` numbers
  in the `vc_index` table of a first run: they depend on the model.

- inbreeding:

  Inbreeding coefficients: a data.frame whose first column is the ID and
  second column F, or a named numeric vector. `NULL` uses F = 0. For
  GBLUP, `F = diag(G) - 1` is common.

- min_rel:

  GEBVs with reliability below this value are set to `NA` before
  computing PA / accuracy. Use `-Inf` to keep all.

- neg_to_na:

  Set negative PA / accuracy / slope to `NA` before the summary (the
  behaviour of the original scripts).

- test_rg:

  Two-trait models only: LRT of rg = 0 (`idh(trait)` genetic structure,
  other structures unchanged, `boundary = FALSE`).

- test_rp:

  Two-trait models only: LRT of *all* covariances = 0 (`idh` for the
  genetic, extra random and residual terms; `df` = number of covariance
  / correlation parameters). This tests whether the traits are
  correlated at all, not rp = 0 as such.

- rp:

  Compute the phenotypic correlation. Set `FALSE` when no individual has
  records for both traits (e.g. one trait in two herds): the residual
  covariance is then not estimable and rp is not meaningful.

- na_x:

  `na.method(x = )`: `"include"` (default, as in the original scripts)
  treats missing covariates / factor levels as 0 in the design matrix -
  centre numeric covariates first; ASReml's own default is `"fail"`.

- dense:

  Put the genomic term in `dense` (faster for GBLUP).

- maxit, workspace:

  Passed to `asreml()` (maxit is 1 for `cv_vc = "fixed"`). ASReml's own
  defaults are 13 and `"128mb"`; two-trait GBLUP usually needs more.

- max_update:

  Number of [`update()`](https://rdrr.io/r/stats/update.html) calls on a
  fit that has not converged.

- asreml_args:

  Named list of further `asreml()` arguments, e.g.
  `list(ai.sing = TRUE)`.

- on_error:

  `"warn"` (default): a failing fit gives a warning and `NA` results,
  the loop continues. `"stop"`: stop at the first failure.

- keep_pred:

  Return the validation predictions (`pred`).

- path:

  Optional `.xlsx` file to write the results to.

- verbose:

  Print progress.

## Value

A list of class `asr_gblup` with `varcomp`, `vc_index`, `param`,
`accuracy` (one row per fold and trait), `summary` and, if `keep_pred`,
`pred`. Fitted model objects are not kept, which keeps memory use low.

## See also

[`asr_scenarios()`](https://tony2015116.github.io/mintyrasr/reference/asr_scenarios.md)
to run several scenarios from `mintyr::cv_scenarios()` in one call.

## Examples

``` r
if (FALSE) { # \dontrun{
library(asreml)
sc <- mintyr::cv_scenarios(pheno, cols = c("adg", "bf"), id = "id",
                           group = "test_herd", v = 5, repeats = 5,
                           seed = 2026, out_type = "df")
res <- asr_gblup(sc$st_sp$param, sc$st_sp$cv, ginv = ginv,
                 fixed = ~ hhys, inbreeding = inbreeding,
                 path = "st_sp_ginv.xlsx")
res$vc_index   # which V is which component
res$summary
} # }
```
