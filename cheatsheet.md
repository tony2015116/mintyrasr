## Setup inputs

    pak::pak("tony2015116/mintyrasr")   # needs licensed ASReml-R 4
    library(asreml); library(mintyr); library(mintyrasr)

    ginv <- G.inverse(G, sparseform = TRUE)$Ginv.sparse
    inbreeding <- data.frame(id = rownames(G), F = diag(G) - 1)

| input | requirement |
|----|----|
| `pheno` | one row per individual: id, group, fixed effects, traits |
| `ginv` | inverse G (or A) with `attr(, "INVERSE") = TRUE`; otherwise `vm()` takes it as G (warning) |
| `inbreeding` | data.frame `(id, F)` or named vector; `NULL` → F = 0 |

## [Scenarios](https://tony2015116.github.io/mintyrasr/articles/mintyrasr.md) mintyr

    sc <- cv_scenarios(pheno, cols = c("adg", "bf"), id = "id",
                       group = "test_herd", v = 5, repeats = 5,
                       seed = 2026, out_type = "df")

| scenario | model                               | train           |
|----------|-------------------------------------|-----------------|
| `st_sp`  | 1 trait, 1 group                    | target group    |
| `st_mp`  | 1 trait, groups combined            | \+ other groups |
| `mt_sp`  | 2 traits, 1 group                   | target group    |
| `mt_mp`  | 2 traits, groups combined           | \+ other groups |
| `mt_bt`  | 1 trait per group as 2 traits (G×E) | \+ other group  |

## [Run](https://tony2015116.github.io/mintyrasr/reference/asr_scenarios.md) all scenarios or one

**[asr_scenarios](https://tony2015116.github.io/mintyrasr/reference/asr_scenarios.md)**(sc,
ginv, scenarios = c("st_sp", "st_mp", "mt_sp", "mt_mp"), dir = ".", tag
= deparse(substitute(ginv)), ...)

One workbook per scenario: `<dir>/<scenario>_<tag>.xlsx`. `"mt_bt"`:
diag residual, no rp.

**[asr_gblup](https://tony2015116.github.io/mintyrasr/reference/asr_gblup.md)**(param,
cv, ginv, fixed = ~ 1, random = NULL, id = "id", traits = NULL,
genetic_cov = "us", random_cov = "us", residual_cov = "us", cv_vc =
"fixed", vpredict = NULL, inbreeding = NULL, min_rel = 0.05, neg_to_na =
TRUE, test_rg = TRUE, test_rp = FALSE, rp = TRUE, na_x = "include",
dense = TRUE, maxit = 100, workspace = "30Gb", max_update = 0,
asreml_args = list(), on_error = "warn", keep_pred = TRUE, path = NULL,
verbose = TRUE)

    res <- asr_scenarios(sc, ginv, fixed = ~ hhys,
                         inbreeding = inbreeding, dir = "results")
    st <- asr_gblup(sc$st_sp$param, sc$st_sp$cv, ginv,
                    fixed = ~ hhys, path = "st_sp.xlsx")
    # own nested tables, other column names
    asr_gblup(all_nest, all_cv, ginv, fixed = ~ hhys,
              traits = c("bfs", "nls"))

`param` and `cv` rows are matched through shared outer columns (`name`,
`pairs`, `test_herd`…). Single / two-trait from `value` /
`value1`+`value2` or `traits`.

## [Models](https://tony2015116.github.io/mintyrasr/articles/models.md) random effects & structures

    # 1 trait:  value ~ fixed,  ~ vm(id, ginv) + random,  ~ idv(units)
    # 2 traits: cbind(a, b) ~ trait + trait:fixed,
    #   ~ genetic_cov(trait):vm(id, ginv) + random_cov(trait):random,
    #   ~ units:residual_cov(trait)
    asr_gblup(..., random = ~ ide(id))            # repeated records
    asr_gblup(..., random = ~ litter, random_cov = "idh")
    asr_gblup(..., genetic_cov = "corgh")         # us too heavy
    asr_gblup(..., residual_cov = "diag", rp = FALSE)

| argument                    | values                                     |
|-----------------------------|--------------------------------------------|
| `genetic_cov`, `random_cov` | `us` `idh`/`diag` `corv` `corh` `corgh`    |
| `residual_cov`              | `us` `idh`/`diag`                          |
| `na_x`                      | `include` (NA covariate = 0) `fail` `omit` |

## [Genetic parameters](https://tony2015116.github.io/mintyrasr/articles/models.md) by component name

Each varcomp row is matched to `names(fit$vparameters)` (= `V` of
`vpredict`) and classified as genetic / random / residual and var / cov
/ cor. Scale rows `units!R` (after `units!units`) and `units:trait!R`
are skipped.

| parameter | formula                                        |
|-----------|------------------------------------------------|
| h2_(i)    | Va_(i) / (all variances of trait i)            |
| rg        | Cov_(a) / √(Va₁Va₂); cor\* → the correlation   |
| rp        | Σ covariances / √(Vp₁Vp₂); cor term → r√(V₁V₂) |

    # us G + us R
    h2_value1 = V1/(V1 + V5); rg = V2/sqrt(V1 * V3)
    # corgh G + diag R
    h2_nls = V2/(V2 + V5); rg = V1
    # own formulas (check res$vc_index first)
    asr_gblup(..., random = ~ ide(id),
      vpredict = list(rep = rep ~ (V1 + V2)/(V1 + V2 + V3)))

fa() / rr() and several genetic terms (maternal `str()`): write formulas
in `vpredict`.

## Tests lrt.asreml(boundary = FALSE)

| argument         | reduced model        | output                     |
|------------------|----------------------|----------------------------|
| `test_rg = TRUE` | genetic `idh(trait)` | `rg_lrt` `rg_p`            |
| `test_rp = TRUE` | all terms `idh`      | `cov_df` `cov_lrt` `cov_p` |

`test_rp`: are the traits correlated at all (all covariances = 0), not
rp = 0 alone.

## [Cross-validation](https://tony2015116.github.io/mintyrasr/articles/mintyrasr.md) per fold

| `cv_vc` | variance components in each fold |
|----|----|
| `"fixed"` | full-data values, `Constraint = "F"`, `maxit = 1` — fast, slight leakage |
| `"warm"` | re-estimated, starting from full-data values |
| `"estimate"` | re-estimated from defaults — no leakage, slowest |

| measure    | definition                                                   |
|------------|--------------------------------------------------------------|
| `rel`      | 1 − PEV / ((1 + F) Va), Va of the fold; \< `min_rel` dropped |
| `pa`       | r(y, GEBV)                                                   |
| `accuracy` | pa / √h2, full-data h2                                       |
| `slope`    | b(y ~ GEBV), ideally 1                                       |

## Output one sheet each

| element | content |
|----|----|
| `varcomp` | summary(fit)\$varcomp + `V` |
| `vc_index` | V → term, role, kind, traits, used |
| `param` | h2, rg, rp (+ SE), `formulas`, converge, loglik, boundary, tests |
| `accuracy` | rep × fold × trait: n, pa, accuracy, slope, va, h2_full |
| `summary` | pa / accuracy / slope: n, mean, sd, se |
| `pred` | validation: obs, gebv, pev, inb, rel |

## Checks & tips

- warnings: no `INVERSE` attribute; IDs not in `ginv`; components at the
  boundary (`B`, also `param$boundary`); `gammaPar = TRUE`; no
  convergence
- ASReml defaults `maxit = 13`, `workspace = "128mb"` are too small for
  two-trait GBLUP
- `max_update = 2` → `update()` on non-converged fits
- `on_error = "warn"`: a failing fit gives `NA`, the loop goes on
- models are not kept → low memory
- `asreml_args = list(ai.sing = TRUE)` for other `asreml()` arguments
