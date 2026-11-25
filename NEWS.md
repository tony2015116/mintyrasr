# mintyrasr 0.3.0

* `traits =` names the response column(s) when the nested data do not use
  `value` / `value1` + `value2` (e.g. `traits = c("bfs", "nls")`).
* `asr_scenarios()` runs the `mt_bt` scenario of `cv_scenarios()` (one trait
  treated as a different trait in each of two groups), with a `diag`
  residual and without rp.
* `test_rp = TRUE`: joint LRT of all covariances = 0. `rp = FALSE` drops the
  phenotypic correlation.
* corgh correlation rows such as `...!trait!b:!trait!a.cor` are recognised.

# mintyrasr 0.2.1

* `V` numbers follow `fit$vparameters`, the order used by `vpredict()`.
* Warning when `ginv` has no `INVERSE = TRUE` attribute, when components are
  fixed at the boundary (`B`, also in `param$boundary`) and when
  `asreml.options(gammaPar = TRUE)`.
* New structure `corv`.

# mintyrasr 0.2.0

* LRT of rg = 0 uses `idh(trait)` and `boundary = FALSE`.
* New columns: `pa` = r(y, GEBV) (the former `accuracy`), `accuracy` =
  pa / sqrt(h2) with the full-data h2, `va`, `h2_full`.
* `cv_vc = "fixed"` (`maxit = 1`), `"warm"` or `"estimate"`.
* `genetic_cov`, `random_cov`: `us`, `idh` / `diag`, `corh`, `corgh`;
  `random = ~ ide(id)` etc.; `na_x`; `asreml_args`.

# mintyrasr 0.1.1

* h2, rg and rp are built from the names of the variance components instead
  of fixed `V1`, `V2`, ... positions; `vc_index` table; `vpredict =` for
  your own formulas; `random =`, `random_cov`, `residual_cov`.

# mintyrasr 0.1.0

* First version: `asr_gblup()` and `asr_scenarios()`.
