# Changelog

## mintyrasr 0.4.0

### Cross-validation without a reliability filter

- `cv_min_rel = -Inf` (new, default): validation individuals are **not**
  filtered by reliability. They have no own record in the training set,
  so a low reliability is expected; a filter kept only well-connected
  individuals, biased PA and shrank the sample. `min_rel` (default 0.05)
  now applies to the final GEBVs of the full-data fit.
- New output `gebv`: GEBV, PEV, F, reliability and
  `reliable = rel >= min_rel` for every individual of each full-data fit
  (`keep_gebv = TRUE`). Tables over Excel’s row limit are written to
  `<path>_<name>.csv`.
- `neg_to_na` now defaults to `FALSE`: dropping negative PA / slopes
  biased the mean upwards. `neg_to_na = TRUE` gives the old behaviour.

### Long runs

- Progress: one line per fit with status, iterations, time and estimated
  time left; a summary line at the end. ASReml’s iteration monitor is
  off while mintyrasr runs (`asreml_trace = TRUE` to keep it).
- `max_update` now defaults to 2: a fit that has not converged is
  continued with [`update()`](https://rdrr.io/r/stats/update.html); the
  number used is in `param$updates` / `accuracy$cv_updates`.
- `nonconv = "skip"` (default): a parameter fit that still has not
  converged is reported with `status = "not_converged"` and its CV rows
  are skipped; a CV fit gives `NA` for its fold. `nonconv = "keep"` uses
  such fits.
- Safe running: errors, non-convergence, boundary estimates and invalid
  LRTs are recorded in the new `$log` (time, stage, label, type,
  message) instead of interrupting the loop; one warning at the end.
  `param` and `accuracy` gain `status`; `summary` uses only rows with
  `status == "ok"`.

## mintyrasr 0.3.2

From the second validation run with ASReml-R 4.2.0.481:

- An LRT where both models have the same number of free parameters (df =
  0, e.g. a covariance restrained at a bound) gives `rg_p = NA` with a
  warning, instead of a misleading p = 0.
- New `min_n = 10`: a fold needs at least `min_n` validation individuals
  with a usable GEBV (after `min_rel`) for PA, accuracy and slope. With
  a low h2 the reliability filter could leave 2-7 animals and PA / slope
  were noise.
- `cv_iter` is documented as the length of the iteration history (2 for
  fixed variances).
- Confirmed on real ASReml-R: `rg_p` equals
  `lrt(full, idh, boundary = FALSE)`; fixed-variance CV fits give no
  warnings; `cv_vc = "warm"` needs about half the iterations of
  `"estimate"` and gives the same results.

## mintyrasr 0.3.1

Fixes found by running the validation script with ASReml-R 4.2.0.481:

- The LRT of rg = 0 (and `test_rp`) uses `lrt()`, the name in ASReml-R
  4.2, with a fallback to `lrt.asreml()` for older releases. Before, the
  call failed silently and `rg_p` was always `NA`. LRT failures now give
  a warning.
- Fixed-variance CV fits (`cv_vc = "fixed"`, `maxit = 1`) no longer pass
  on ASReml’s expected “Log-likelihood not converged” warnings.
- `accuracy` gains `cv_converge` and `cv_iter` (iterations of each
  fold’s fit).
- Confirmed on real ASReml-R: component names for idv, us, idh/diag,
  corgh, corh and corv, `vparameters` order, coefficient names, and h2 /
  rg / rp / repeatability identical to hand-written `vpredict()`
  formulas.

## mintyrasr 0.3.0

- `traits =` names the response column(s) when the nested data do not
  use `value` / `value1` + `value2` (e.g. `traits = c("bfs", "nls")`).
- [`asr_scenarios()`](https://tony2015116.github.io/mintyrasr/reference/asr_scenarios.md)
  runs the `mt_bt` scenario of
  [`cv_scenarios()`](https://tony2015116.github.io/mintyr/reference/cv_scenarios.html)
  (one trait treated as a different trait in each of two groups), with a
  `diag` residual and without rp.
- `test_rp = TRUE`: joint LRT of all covariances = 0. `rp = FALSE` drops
  the phenotypic correlation.
- corgh correlation rows such as `...!trait!b:!trait!a.cor` are
  recognised.

## mintyrasr 0.2.1

- `V` numbers follow `fit$vparameters`, the order used by `vpredict()`.
- Warning when `ginv` has no `INVERSE = TRUE` attribute, when components
  are fixed at the boundary (`B`, also in `param$boundary`) and when
  `asreml.options(gammaPar = TRUE)`.
- New structure `corv`.

## mintyrasr 0.2.0

- LRT of rg = 0 uses `idh(trait)` and `boundary = FALSE`.
- New columns: `pa` = r(y, GEBV) (the former `accuracy`), `accuracy` =
  pa / sqrt(h2) with the full-data h2, `va`, `h2_full`.
- `cv_vc = "fixed"` (`maxit = 1`), `"warm"` or `"estimate"`.
- `genetic_cov`, `random_cov`: `us`, `idh` / `diag`, `corh`, `corgh`;
  `random = ~ ide(id)` etc.; `na_x`; `asreml_args`.

## mintyrasr 0.1.1

- h2, rg and rp are built from the names of the variance components
  instead of fixed `V1`, `V2`, … positions; `vc_index` table;
  `vpredict =` for your own formulas; `random =`, `random_cov`,
  `residual_cov`.

## mintyrasr 0.1.0

- First version:
  [`asr_gblup()`](https://tony2015116.github.io/mintyrasr/reference/asr_gblup.md)
  and
  [`asr_scenarios()`](https://tony2015116.github.io/mintyrasr/reference/asr_scenarios.md).
