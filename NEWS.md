# mintyrasr 0.3.2

From the second validation run with ASReml-R 4.2.0.481:

* An LRT where both models have the same number of free parameters (df = 0,
  e.g. a covariance restrained at a bound) gives `rg_p = NA` with a warning,
  instead of a misleading p = 0.
* New `min_n = 10`: a fold needs at least `min_n` validation individuals with
  a usable GEBV (after `min_rel`) for PA, accuracy and slope. With a low h2
  the reliability filter could leave 2-7 animals and PA / slope were noise.
* `cv_iter` is documented as the length of the iteration history (2 for
  fixed variances).
* Confirmed on real ASReml-R: `rg_p` equals `lrt(full, idh, boundary = FALSE)`;
  fixed-variance CV fits give no warnings; `cv_vc = "warm"` needs about half
  the iterations of `"estimate"` and gives the same results.

# mintyrasr 0.3.1

Fixes found by running the validation script with ASReml-R 4.2.0.481:

* The LRT of rg = 0 (and `test_rp`) uses `lrt()`, the name in ASReml-R 4.2,
  with a fallback to `lrt.asreml()` for older releases. Before, the call
  failed silently and `rg_p` was always `NA`. LRT failures now give a warning.
* Fixed-variance CV fits (`cv_vc = "fixed"`, `maxit = 1`) no longer pass on
  ASReml's expected "Log-likelihood not converged" warnings.
* `accuracy` gains `cv_converge` and `cv_iter` (iterations of each fold's fit).
* Confirmed on real ASReml-R: component names for idv, us, idh/diag, corgh,
  corh and corv, `vparameters` order, coefficient names, and h2 / rg / rp /
  repeatability identical to hand-written `vpredict()` formulas.

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
