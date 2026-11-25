# Get started with mintyrasr

[中文版](https://tony2015116.github.io/mintyrasr/articles/mintyrasr-zh.html)

ASReml-R is licensed software, so the code on this page is not run when
the site is built. Every chunk runs as shown on a computer with ASReml-R
4.

## What the package does

mintyr prepares the data for the usual comparison of genomic prediction
models: for every scenario a nested table for **parameter estimation**
(all records) and a table of **training / validation sets** (one row per
repeat x fold), all validated on the same individuals. mintyrasr takes
both tables and, for every row:

1.  fits the model to the full data: variance components, h2 and, for
    two traits, rg, rp and an LRT of rg = 0;
2.  refits it to every training set (variance components fixed,
    warm-started or re-estimated);
3.  predicts the GEBVs of the validation individuals and computes their
    reliability, the predictive ability, the accuracy and the slope;
4.  summarises over folds and writes one Excel workbook per scenario.

| function | use |
|----|----|
| `asr_gblup(param, cv, ginv, ...)` | one scenario (or any nested table + CV table) |
| `asr_scenarios(sc, ginv, ...)` | all scenarios of `mintyr::cv_scenarios()` in one call |

## Installation

``` r

# ASReml-R 4 must be installed and licensed first (https://vsni.co.uk/)
install.packages("mintyr")
pak::pak("tony2015116/mintyrasr")
```

## 1. Prepare the inputs

You need three objects:

- **phenotypes**: one row per individual, with the ID, the grouping
  column (herd, line, …), fixed effects and the traits;
- **`ginv`**: the inverse genomic relationship matrix in a form `vm()`
  accepts. It must carry `attr(ginv, "INVERSE") <- TRUE`, otherwise
  `vm()` treats it as G itself (mintyrasr warns);
- **inbreeding** (optional): a data.frame `(id, F)` for the reliability
  `1 - PEV / ((1 + F) Va)`. For GBLUP, `F = diag(G) - 1` is common.
  Without it, F = 0 is used.

``` r

library(asreml)
library(ASRgenomics)
library(mintyr)
library(mintyrasr)

pheno <- import_xlsx("pheno.xlsx")             # id, test_herd, hhys, adg, bf, fcr

G    <- G.matrix(M = snp, method = "VanRaden")$G
ginv <- G.inverse(G = G, sparseform = TRUE)$Ginv.sparse   # carries INVERSE = TRUE
inbreeding <- data.frame(id = rownames(G), F = diag(G) - 1)
```

## 2. Build the scenarios with mintyr

``` r

sc <- cv_scenarios(pheno, cols = c("adg", "bf", "fcr"), id = "id",
                   group = "test_herd", v = 5, repeats = 5, seed = 2026,
                   out_type = "df")      # ASReml-R needs data.frames
sc
```

| scenario | model | parameter data | training set |
|----|----|----|----|
| `st_sp` | single-trait, single group | trait x group | target group |
| `st_mp` | single-trait, groups combined | trait | target group + all other groups |
| `mt_sp` | two-trait, single group | trait pair x group | target group |
| `mt_mp` | two-trait, groups combined | trait pair | target group + all other groups |
| `mt_bt` | one trait as a different trait in each group | trait | target group + all other groups |

## 3. Run all scenarios

``` r

res <- asr_scenarios(sc, ginv, fixed = ~ hhys, inbreeding = inbreeding,
                     dir = "results")
# results/st_sp_ginv.xlsx, st_mp_ginv.xlsx, mt_sp_ginv.xlsx, mt_mp_ginv.xlsx

res$st_sp$summary
res$mt_sp$param
```

Add `"mt_bt"` to `scenarios` for one trait treated as a different trait
in each of two groups (G x E). It is fitted with a `diag` residual and
without rp, because no individual has records in both groups:

``` r

res_bt <- asr_scenarios(sc, ginv, scenarios = "mt_bt", fixed = ~ hhys,
                        genetic_cov = "corgh", inbreeding = inbreeding)
res_bt$mt_bt$param        # h2 per group, rg between groups (+ SE), rg_p
```

## 4. One scenario at a time

[`asr_gblup()`](https://tony2015116.github.io/mintyrasr/reference/asr_gblup.md)
takes a nested parameter table and a CV table. They are matched through
their shared outer columns (`name`, `pairs`, `test_herd`, …).

``` r

st <- asr_gblup(sc$st_sp$param, sc$st_sp$cv, ginv = ginv,
                fixed = ~ hhys, inbreeding = inbreeding,
                path = "results/st_sp_ginv.xlsx")
st
```

Tables built by hand with `w2l_nest()` / `c2p_nest()` + `nest_cv()` work
too. When the response columns are not called `value` / `value1` +
`value2`, name them in `traits`:

``` r

nested <- w2l_nest(pheno, cols = c("adg", "bf"), by = "test_herd", out_type = "df")
cvtab  <- nest_cv(nested, v = 5, repeats = 2, seed = 1, out_type = "df")
asr_gblup(nested, cvtab, ginv, fixed = ~ hhys)

asr_gblup(all_nest5, all_cv5, ginv, fixed = ~ hhys, traits = c("bfs", "nls"),
          genetic_cov = "corgh", residual_cov = "diag", rp = FALSE)
```

## 5. What you get

Each result is a list; every element becomes one sheet of the workbook.

| element | content |
|----|----|
| `varcomp` | `summary(fit)$varcomp` of every parameter fit, with its `V` number |
| `vc_index` | which `V` is which component: term, role (genetic / random / residual), kind (var / cov / cor), traits, used |
| `param` | h2 (two traits: per trait), rg, rp with SE; `formulas` used; `converge`, `loglik`, `boundary`; `rg_lrt`, `rg_p` (and `cov_*` with `test_rp = TRUE`) |
| `accuracy` | per repeat x fold x trait: `n`, `pa`, `accuracy`, `slope`, `va`, `h2_full`, `n_obs` |
| `summary` | `pa`, `accuracy`, `slope`: n, mean, sd, se over folds |
| `pred` | validation individuals: `obs`, `gebv`, `pev`, `inb` (F), `rel` |

- `pa` (predictive ability) = r(y, GEBV) on the validation individuals;
- `accuracy` = pa / sqrt(h2), with h2 from the full data (Dekkers 2007);
- `slope` = b(y ~ GEBV); 1 means no inflation or deflation of the GEBVs.

GEBVs with reliability below `min_rel` (default 0.05) are left out, and
with `neg_to_na = TRUE` (default) negative values are set to `NA` before
the summary.

## 6. Variance components in cross-validation

| `cv_vc` | in each fold | note |
|----|----|----|
| `"fixed"` (default) | fixed at the full-data estimates (`Constraint = "F"`, `maxit = 1`) | fastest; slight leakage, because the validation records were used to estimate them |
| `"warm"` | re-estimated from the training set, starting from the full-data values |  |
| `"estimate"` | re-estimated from default start values | no leakage; slowest |

The `va` used for the reliability is always taken from the model of that
fold.

## 7. Practical notes

- Fitted models are not kept: memory stays low even with many folds.
- A failing fit gives a warning and `NA` results; the loop carries on
  (`on_error = "stop"` to stop instead). `max_update = 2` calls
  [`update()`](https://rdrr.io/r/stats/update.html) on fits that have
  not converged.
- Two-trait GBLUP needs more than ASReml’s defaults (`maxit = 13`,
  `workspace = "128mb"`); mintyrasr uses `maxit = 100`,
  `workspace = "30Gb"`.
- `na_x = "include"` (default, as in many scripts) treats missing
  covariates as 0 - centre numeric covariates first. ASReml’s own
  default is `"fail"`.
- Further `asreml()` arguments: `asreml_args = list(ai.sing = TRUE)`.

Next: [Models and genetic
parameters](https://tony2015116.github.io/mintyrasr/articles/models.html) -
random effects, variance structures, how the `V` numbers are found,
tests and your own `vpredict()` formulas.
