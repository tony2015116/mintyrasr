# mintyrasr 0.4.0

* The LRTs of rg = 0 and of all covariances use `asreml::lrt()`: ASReml-R 4.2
  no longer exports `lrt.asreml()`, so `rg_lrt` / `rg_p` were silently `NA`.
  If the asreml LRT fails, LR is computed from the log-likelihoods (warning).
* `min_rel` now defaults to `-Inf` (keep all GEBVs). With 0.05, validation
  individuals of small single-group scenarios (`st_sp`, `mt_sp`, `mt_bt`)
  were almost all dropped and PA was `NA`. A warning is given when the filter
  leaves fewer than 3 GEBVs.
* Progress with count, %, elapsed time and ETA; summary at the end;
  `asr_scenarios()` shows scenario k/n. `trace = FALSE` (default) hides
  ASReml's iteration log.
* Safe run: `max_update` now defaults to 3 (`update()` on fits that have not
  converged). New `on_nonconverge = "skip"` (default) / `"keep"`: a fit that
  still has not converged is warned about, recorded and nothing is built on
  it. New `log` table (also in the workbook) with stage, label, status,
  `n_update`, seconds and the ASReml warnings of every fit; new columns
  `param$n_update` and `accuracy$status`. The expected "not converged"
  warnings of `cv_vc = "fixed"` fits (`maxit = 1`) are no longer shown.

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
