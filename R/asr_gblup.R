#' GBLUP parameter estimation and cross-validation with ASReml-R
#'
#' @description
#' `asr_gblup()` runs the whole GBLUP workflow for one scenario directly in R,
#' so no input files have to be exported for ASReml:
#'
#' 1. fits the model to every data set in `param`: variance components,
#'    heritability and, for two-trait models, genetic / phenotypic
#'    correlations plus an LRT of \sb{r}{g} = 0 (`idh(trait)` vs the full genetic
#'    structure, `boundary = FALSE`);
#' 2. refits the model to every training set in `cv`, with the variance
#'    components handled as set by `cv_vc`;
#' 3. predicts GEBVs of the validation individuals (all of them: they have
#'    no own record, so a low reliability is expected and not filtered by
#'    default) and returns per fold the predictive ability
#'    `pa = r(y, GEBV)`, the accuracy pa / \sqrthsq{} (\hsq{} from the full data,
#'    Dekkers 2007) and the slope `b(y ~ GEBV)`;
#' 4. summarises them over folds (n, mean, sd, se), returns the GEBVs of all
#'    individuals from the full-data fit with their reliability (flag
#'    `reliable = rel >= min_rel`) and, optionally, writes everything to one
#'    Excel workbook.
#'
#' Long runs are safe: every error, non-convergence, boundary estimate and
#' invalid LRT is recorded in `$log` (with time, stage and label) instead of
#' stopping the loop; fits that do not converge are first continued with
#' `update()` (`max_update`), and are then skipped (`nonconv = "skip"`) or
#' kept. A progress line per fit shows status, iterations, time and the
#' estimated time left.
#'
#' \hsq{}, \sb{r}{g} and \sb{r}{p} are built from the **names** of the variance components,
#' not from fixed `V` positions, so they stay correct when random effects
#' are added (`random = ~ ide(id) + litter`) or covariance structures change
#' (`us`, `idh`/`diag`, `corh`, `corgh`):
#' \hsq{} = \sb{V}{a} / (all variances of the trait), \sb{r}{g} = \sb{Cov}{a} / \rt{\sb{V}{a,1} \sb{V}{a,2}}
#' (the correlation itself for `corh`/`corgh`), \sb{r}{p} = sum of all covariances /
#' \rt{\sb{V}{p,1} \sb{V}{p,2}}. The `vc_index` table shows which `V` each component is and
#' `param$formulas` the formulas used. `V` numbers follow `fit$vparameters`,
#' the order `vpredict()` uses. Components fixed at the boundary (bound `B`)
#' are reported with a warning and in `param$boundary`. `fa()` / `rr()` terms and models with
#' several genetic terms (e.g. maternal `str()` models) need your own
#' formulas in `vpredict`.
#'
#' Single-trait or two-trait models are chosen from the nested data:
#' a `value` column gives a single-trait model, `value1` + `value2` a
#' two-trait model (layout of [mintyr::w2l_nest()] / [mintyr::c2p_nest()]).
#'
#' @param param Nested table used for parameter estimation, with a
#'   list-column `data` (e.g. `cv_scenarios(...)$st_sp$param`).
#' @param cv Cross-validation table with list-columns `train` and `validate`
#'   (e.g. `cv_scenarios(...)$st_sp$cv` or the output of [mintyr::nest_cv()]).
#'   `NULL` (default): genetic parameters only, no cross-validation - e.g.
#'   for every group of [mintyr::w2l_nest()] or every trait pair of
#'   [mintyr::c2p_nest()].
#'   Each row is matched to the row of `param` that shares its outer columns.
#' @param ginv Inverse relationship matrix as used by `vm()`: \inv{A} for
#'   pedigree BLUP (`ainverse(ped)`), \inv{G} for GBLUP
#'   (`G.inverse(G, sparseform = TRUE)$Ginv.sparse`) or \inv{H} for single-step
#'   GBLUP (`H.inverse(...)`). It must carry `attr(, "INVERSE") = TRUE`,
#'   otherwise `vm()` takes it as the relationship matrix itself (a warning is
#'   given). With `ainverse()`, F is taken from its `inbreeding` attribute when
#'   `inbreeding` is `NULL`.
#' @param fixed One-sided formula of fixed effects, e.g. `~ hhys`. For
#'   two-trait models every term is nested in trait (`trait + trait:hhys`).
#' @param random Optional one-sided formula of further random effects, e.g.
#'   `~ ide(id)` (permanent environment, repeated records), `~ litter + pen`,
#'   `~ vm(dam, ginv) + ide(dam)` (maternal genetic and permanent maternal
#'   environment), `~ vm(id, ginv):herd` (genotype \by{} herd interaction; the
#'   factor is given the variance structure `idv(herd)` ASReml needs). In
#'   two-trait models each term gets the structure `random_cov`. Each term is
#'   also reported as a proportion of the phenotypic variance
#'   (`prop_<term>`), with `rep` (repeatability) for `ide(<id>)` and
#'   `rB_<X>` (type-B genetic correlation) for `vm(<id>, ginv):X`.
#' @param id Name of the individual ID column in the nested data.
#' @param traits Optional response column(s), one or two, when the nested
#'   data do not use `value` / `value1` + `value2`, e.g. `c("bfs", "nls")`
#'   for one trait recorded in two herds (`cv_scenarios(...)$mt_bt`).
#' @param genetic_cov,random_cov,residual_cov Two-trait models: variance
#'   structure across traits for the genomic term, the extra random terms
#'   and the residual: `"us"`, `"idh"`, `"diag"`, `"corv"`, `"corh"` or `"corgh"`
#'   (residual: `"us"`, `"idh"`, `"diag"`).
#' @param cv_vc How the variance components are handled in each fold:
#'   `"fixed"` (default, fast): fixed at the full-data estimates
#'   (`Constraint = "F"`, `maxit = 1`; slight information leakage because
#'   validation records were used to estimate them); `"warm"`: re-estimated
#'   from the training set, starting from the full-data values;
#'   `"estimate"`: re-estimated from the training set from default start
#'   values (no leakage, slowest).
#' @param vpredict Optional named list of extra `vpredict()` formulas, e.g.
#'   `list(rep = rep ~ (V1 + V2) / (V1 + V2 + V3))`. Check the `V` numbers in
#'   the `vc_index` table of a first run: they depend on the model.
#' @param inbreeding Inbreeding coefficients: a data.frame whose first column
#'   is the ID and second column F, or a named numeric vector. `NULL` uses
#'   F = 0. For GBLUP, `F = diag(G) - 1` is common.
#' @param cv_min_rel Reliability filter for the validation individuals in
#'   cross-validation. Default `-Inf`: none. Validation individuals have no
#'   own record, so low reliabilities are expected; a filter keeps only the
#'   best-connected individuals, biases PA and shrinks the sample.
#' @param min_rel Reliability threshold for the final GEBVs of the full-data
#'   fit (`$gebv`): sets `reliable = rel >= min_rel`; no rows are dropped.
#' @param min_n Minimum number of validation individuals with a GEBV and a
#'   record for PA, accuracy and slope of a fold; fewer gives `NA`
#'   (`n` still shows the number). With a low \hsq{} the reliability filter can
#'   leave only a handful of individuals, whose PA and slope are noise.
#' @param neg_to_na Set negative PA / accuracy / slope to `NA` before the
#'   summary. Default `FALSE`: dropping negative values biases the mean
#'   upwards; `TRUE` reproduces the original scripts.
#' @param test_rg Two-trait models only: LRT of \sb{r}{g} = 0 (`idh(trait)` genetic
#'   structure, other structures unchanged, `boundary = FALSE`). If both
#'   models end up with the same number of free parameters (df = 0, e.g. a
#'   restrained covariance), `rg_p` is `NA` with a warning.
#' @param test_rp Two-trait models only: LRT of *all* covariances = 0
#'   (`idh` for the genetic, extra random and residual terms; `df` = number
#'   of covariance / correlation parameters). This tests whether the traits
#'   are correlated at all, not \sb{r}{p} = 0 as such.
#' @param rp Compute the phenotypic correlation. Set `FALSE` when no
#'   individual has records for both traits (e.g. one trait in two herds):
#'   the residual covariance is then not estimable and \sb{r}{p} is not meaningful.
#' @param na_x `na.method(x = )`: `"include"` (default, as in the original
#'   scripts) treats missing covariates / factor levels as 0 in the design
#'   matrix - centre numeric covariates first; ASReml's own default is `"fail"`.
#' @param dense Put the relationship term in `dense`. `"auto"` (default):
#'   yes for a (nearly) full inverse such as \inv{G} or \inv{H}, no for a sparse
#'   one such as \inv{A} from `ainverse()`, where `dense` would cost much memory.
#' @param maxit Passed to `asreml()` (1 for `cv_vc = "fixed"`). ASReml's own
#'   default is 13; two-trait GBLUP usually needs more.
#' @param workspace Passed to `asreml()`. `"auto"` (default): estimated from the
#'   size of the dense part (individuals \by{} traits when the relationship term
#'   is in `dense`), at least `"512mb"`. ASReml reserves the whole workspace
#'   for every fit, so a large fixed value (e.g. `"30Gb"`) slows every fit
#'   down and can exhaust the memory of the computer.
#' @param max_update Number of `update()` calls (continue iterating from the
#'   current estimates) on a fit that has not converged; recorded in
#'   `param$updates` and `accuracy$cv_updates`.
#' @param nonconv What to do with a fit that has still not converged:
#'   `"skip"` (default) - a parameter fit is reported with
#'   `status = "not_converged"` and its CV rows are skipped, a CV fit gives `NA`
#'   for that fold; `"keep"` - use it anyway. Always recorded in `$log`.
#' @param asreml_trace Show ASReml's iteration monitor. Default `FALSE`: the
#'   progress lines of mintyrasr are shown instead.
#' @param asreml_args Named list of further `asreml()` arguments, e.g.
#'   `list(ai.sing = TRUE)`.
#' @param on_error `"warn"` (default): a failing fit is recorded in `$log`,
#'   gives `NA` results and the loop continues (one warning at the end).
#'   `"stop"`: stop at the first failure.
#' @param keep_pred Return the validation predictions (`pred`).
#' @param keep_gebv Return the GEBVs of all individuals from the full-data
#'   fits (`gebv`: gebv, pev, inb (F), rel, reliable).
#' @param path Optional `.xlsx` file to write the results to.
#' @param verbose Print progress.
#'
#' @return A list of class `asr_gblup`:
#'   * `varcomp`, `vc_index`: variance components and which `V` is which;
#'   * `param`: `h2`, `rg`, `rp` (+ SE), `formulas`, `status`, `converge`, `updates`,
#'     `iterations`, `loglik`, `boundary`, tests;
#'   * `accuracy`: one row per fold and trait: `status`, `n`, `pa`, `accuracy`,
#'     `slope`, `va`, `h2_full`, `cv_converge`, `cv_updates`, `cv_iter` (length
#'     of the iteration history: about the number of REML iterations, 2 for
#'     fixed variances);
#'   * `summary`: n, mean, sd, se over the folds with `status == "ok"`;
#'   * `pred` (if `keep_pred`): validation individuals;
#'   * `gebv` (if `keep_gebv`): all individuals of the full-data fits;
#'   * `log`: every error, non-convergence, boundary estimate and invalid LRT.
#'
#'   A `status` other than `"ok"` (`failed`, `not_converged`, `param_failed`,
#'   `param_not_converged`, `no_start_values`) marks rows left out of
#'   `summary`. Fitted model objects are not kept, which keeps memory use low.
#'
#' @seealso [asr_scenarios()] to run several scenarios from
#'   [mintyr::cv_scenarios()] in one call.
#' @export
#' @examples
#' \dontrun{
#' library(asreml)
#' sc <- mintyr::cv_scenarios(pheno, cols = c("adg", "bf"), id = "id",
#'                            group = "test_herd", v = 5, repeats = 5,
#'                            seed = 2026, out_type = "df")
#' res <- asr_gblup(sc$st_sp$param, sc$st_sp$cv, ginv = ginv,
#'                  fixed = ~ hhys, inbreeding = inbreeding,
#'                  path = "st_sp_ginv.xlsx")
#' res$vc_index   # which V is which component
#' res$summary
#' }
asr_gblup <- function(param, cv = NULL, ginv, fixed = ~ 1, random = NULL, id = "id", traits = NULL,
                      genetic_cov = "us", random_cov = "us", residual_cov = "us",
                      cv_vc = c("fixed", "warm", "estimate"), vpredict = NULL,
                      inbreeding = NULL, cv_min_rel = -Inf, min_rel = 0.05, min_n = 10L,
                      neg_to_na = FALSE, test_rg = TRUE, test_rp = FALSE, rp = TRUE,
                      na_x = c("include", "fail", "omit"), dense = "auto", maxit = 100L,
                      workspace = "auto", max_update = 2L, nonconv = c("skip", "keep"),
                      asreml_args = list(), asreml_trace = FALSE, on_error = c("warn", "stop"),
                      keep_pred = TRUE, keep_gebv = TRUE, path = NULL, verbose = TRUE) {
  on_error <- match.arg(on_error)
  cv_vc    <- match.arg(cv_vc)
  na_x     <- match.arg(na_x)
  nonconv  <- match.arg(nonconv)
  covs <- c("us", "idh", "diag", "corv", "corh", "corgh")
  genetic_cov  <- match.arg(genetic_cov, covs)
  random_cov   <- match.arg(random_cov, covs)
  residual_cov <- match.arg(residual_cov, c("us", "idh", "diag"))
  if (!requireNamespace("asreml", quietly = TRUE))
    stop("Package 'asreml' (ASReml-R 4) is required.", call. = FALSE)
  if (!inherits(fixed, "formula") || length(fixed) != 2L)
    stop("`fixed` must be a one-sided formula, e.g. ~ hhys.", call. = FALSE)
  if (!is.null(random) && (!inherits(random, "formula") || length(random) != 2L))
    stop("`random` must be a one-sided formula, e.g. ~ ide(id).", call. = FALSE)
  if (!is.null(vpredict) && (!is.list(vpredict) || is.null(names(vpredict))))
    stop("`vpredict` must be a named list of formulas.", call. = FALSE)
  if (!is.list(asreml_args) || (length(asreml_args) && is.null(names(asreml_args))))
    stop("`asreml_args` must be a named list.", call. = FALSE)

  P <- data.table::as.data.table(param)
  if (!"data" %in% names(P)) stop("`param` needs a list-column `data`.", call. = FALSE)
  C <- if (is.null(cv)) NULL else data.table::as.data.table(cv)   # NULL: parameters only
  if (!is.null(C) && !all(c("train", "validate") %in% names(C)))
    stop("`cv` needs list-columns `train` and `validate`.", call. = FALSE)

  p_outer <- .outer_cols(P)
  c_outer <- if (is.null(C)) character(0) else .outer_cols(C)
  keys      <- intersect(p_outer, c_outer)               # e.g. name / pairs (+ herd)
  fold_cols <- intersect(c_outer, c("rep", "fold", "id", "id2"))
  grp_cols  <- setdiff(c_outer, c(keys, fold_cols))      # validated group (mp scenarios)
  if (!is.null(C) && !length(keys) && nrow(P) != 1L)
    stop("`param` and `cv` share no outer column to match rows on.", call. = FALSE)

  d1   <- P$data[[1L]]
  if (!is.null(traits)) {
    traits <- as.character(traits)
    if (!length(traits) %in% 1:2) stop("`traits` must name one or two columns.", call. = FALSE)
    bad <- setdiff(traits, names(d1))
    if (length(bad)) stop("`traits` not in the nested data: ", paste(bad, collapse = ", "), call. = FALSE)
    resp <- traits
    type <- if (length(traits) == 1L) "st" else "mt"
  } else {
    type <- if (all(c("value1", "value2") %in% names(d1))) "mt"
            else if ("value" %in% names(d1)) "st"
            else stop("Nested data need `value` / `value1`+`value2`, or name the response ",
                      "column(s) in `traits`, e.g. traits = c(\"bfs\", \"nls\").", call. = FALSE)
    if (type == "mt" && "value3" %in% names(d1))
      stop("Only two-trait models are supported (c2p_nest(pairs_n = 2)).", call. = FALSE)
    resp <- if (type == "st") "value" else c("value1", "value2")
  }
  for (d in c(P$data, C$train[1L])) if (!id %in% names(d))
    stop(sprintf("ID column '%s' not found in the nested data.", id), call. = FALSE)

  .check_ginv(ginv)
  gp_opt <- tryCatch(asreml::asreml.options()$gammaPar, error = function(e) NULL)
  if (isTRUE(gp_opt))
    warning("asreml.options(gammaPar = TRUE): the automatic h2/rg/rp formulas assume the ",
            "default variance (sigma) parameterisation.", call. = FALSE)
  gids <- .ginv_ids(ginv)
  if (!is.null(gids)) {
    pid <- unique(as.character(unlist(lapply(P$data, `[[`, id))))
    miss <- setdiff(pid, gids)
    if (length(miss))
      warning(sprintf("%d of %d individuals in `param` are not in `ginv` (e.g. %s); they get no GEBV.",
                      length(miss), length(pid), paste(utils::head(miss, 3), collapse = ", ")),
              call. = FALSE)
  }
  Ftab <- .prep_F(inbreeding)
  if (is.null(Ftab) && !is.null(attr(ginv, "inbreeding")) && !is.null(gids) &&
      length(attr(ginv, "inbreeding")) == length(gids)) {
    Ftab <- data.table::data.table(id = gids, F = as.numeric(attr(ginv, "inbreeding")))
    if (verbose) message("F taken from attr(ginv, \"inbreeding\") (e.g. ainverse()).")
  }
  if (is.null(Ftab) && verbose) message("`inbreeding` is NULL: F = 0 is used for all individuals.")
  if (identical(dense, "auto")) dense <- .is_dense(ginv)
  if (!is.logical(dense) || length(dense) != 1L)
    stop("`dense` must be TRUE, FALSE or \"auto\".", call. = FALSE)
  if (identical(workspace, "auto")) workspace <- .auto_workspace(ginv, length(resp), dense)

  # ASReml's iteration monitor floods the console in long runs: off while we run
  if (!asreml_trace) {
    old_trace <- tryCatch(asreml::asreml.options()$trace, error = function(e) NULL)
    if (!is.null(old_trace)) {
      try(asreml::asreml.options(trace = FALSE), silent = TRUE)
      on.exit(try(asreml::asreml.options(trace = old_trace), silent = TRUE), add = TRUE)
    }
  }

  # ---- safe run: every error / warning goes to the log, the loop goes on -----
  issues <- list()
  note <- function(stage, label, type, msg) {
    issues[[length(issues) + 1L]] <<- data.table::data.table(
      time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), stage = stage, label = label,
      type = type, message = msg)
  }
  run <- function(expr, stage, label) {
    withCallingHandlers(
      tryCatch(expr, error = function(e) {
        note(stage, label, "error", conditionMessage(e))
        if (on_error == "stop")
          stop(sprintf("ASReml fit failed for %s: %s", label, conditionMessage(e)), call. = FALSE)
        NULL
      }),
      warning = function(w) {
        note(stage, label, "warning", conditionMessage(w))
        invokeRestart("muffleWarning")
      })
  }
  conv_ok <- function(fit) isTRUE(fit$converge)
  n_upd   <- function(fit) { k <- attr(fit, "n_update"); if (is.null(k)) 0L else as.integer(k) }

  spec <- list(type = type, resp = resp, fixed = fixed, random = random, id = id,
               genetic_cov = genetic_cov, random_cov = random_cov, residual_cov = residual_cov,
               dense = dense, na_x = na_x, workspace = workspace, max_update = max_update,
               asreml_args = asreml_args)
  fit1 <- function(d, ...) .asr_fit(d, spec = spec, relmat = ginv, ...)

  np <- nrow(P); nc <- if (is.null(C)) 0L else nrow(C)
  t_start <- Sys.time()
  say <- function(...) if (verbose) message(...)
  say(sprintf("asr_gblup: %s model | %d parameter fit(s) + %d CV fit(s)%s | dense = %s | workspace = %s",
              if (type == "mt") "two-trait" else "single-trait", np, nc,
              if (nc) sprintf(" | cv_vc = \"%s\"", cv_vc) else "", dense, workspace))
  fit_note <- function(fit, t1) {
    if (is.null(fit)) return("FAILED")
    st <- if (conv_ok(fit)) "ok" else "NOT CONVERGED"
    if (n_upd(fit)) st <- sprintf("%s after %d update()", st, n_upd(fit))
    it <- .n_iter(fit)
    sprintf("%s | %s%s", st, if (is.na(it)) "" else paste0(it, " it, "),
            .fmt_time(difftime(Sys.time(), t1, units = "secs")))
  }

  # ---- 1. parameter estimation (full data) -----------------------------------
  vc_rows <- par_rows <- idx_rows <- vc_cv <- h2_list <- gebv_rows <- vector("list", np)
  p_status <- rep("ok", np)
  for (i in seq_len(np)) {
    lab  <- .label(P[i], p_outer)
    t1   <- Sys.time()
    dat  <- P$data[[i]]
    meta <- P[i, p_outer, with = FALSE]
    fit  <- run(fit1(dat, maxit = maxit), "param", lab)
    say(sprintf("[param %d/%d] %s | %s", i, np, lab, fit_note(fit, t1)))
    if (is.null(fit)) {
      p_status[i] <- "failed"
      par_rows[[i]] <- cbind(meta, status = "failed", converge = NA)
      next
    }
    if (!conv_ok(fit)) {
      p_status[i] <- if (nonconv == "skip") "not_converged" else "not_converged_kept"
      note("param", lab, "not_converged", sprintf(
        "not converged after %d update(); %s", n_upd(fit),
        if (nonconv == "skip") "its CV rows are skipped" else "kept (nonconv = \"keep\")"))
    }
    vc  <- as.data.frame(summary(fit)$varcomp)
    idx <- .vc_index(vc, type, id, resp, names(fit$vparameters))
    bnd <- idx$item[idx$bound == "B"]
    if (length(bnd)) note("param", lab, "boundary",
                          paste("fixed at the boundary (B):", paste(unique(bnd), collapse = "; ")))
    vc_rows[[i]]  <- cbind(meta[rep(1L, nrow(vc))], V = idx$V[match(rownames(vc), idx$item)],
                           item = rownames(vc), vc)
    idx_rows[[i]] <- cbind(meta[rep(1L, nrow(idx))], idx[order(V), !"component"])
    gp <- run(.genetic_params(fit, c(.vp_formulas(idx, resp, rp = rp, id = id), vpredict)), "param", lab)
    if (is.null(gp)) gp <- data.table::data.table(formulas = NA_character_)
    h2_list[[i]] <- stats::setNames(
      vapply(resp, function(r) {
        nm <- if (type == "st") "h2" else paste0("h2_", r)
        if (nm %in% names(gp)) gp[[nm]] else NA_real_
      }, numeric(1L)), resp)
    par_rows[[i]] <- cbind(meta, status = p_status[i], gp, converge = conv_ok(fit),
                           updates = n_upd(fit), iterations = .n_iter(fit),
                           loglik = if (is.null(fit$loglik)) NA_real_ else fit$loglik,
                           boundary = paste(unique(bnd), collapse = "; "))
    va_full <- run(.va(idx, resp), "param", lab)
    if (keep_gebv && !is.null(va_full)) {
      g <- run(.gebv_table(fit, id, resp, va_full, Ftab), "param", lab)
      if (!is.null(g)) {
        g[, trait := .trait_names(P[i], keys, type, resp, given = !is.null(traits))[[1L]][match(trait, resp)]]
        g[, reliable := !is.na(rel) & rel >= min_rel]
        gebv_rows[[i]] <- cbind(meta[rep(1L, nrow(g))], g)
      }
    }
    lrt_test <- function(red_fit, what) {          # LRT against a reduced model
      if (is.null(red_fit)) return(NULL)
      if (!conv_ok(red_fit) && nonconv == "skip") {
        note("param", lab, "not_converged", paste(what, "not converged: test skipped"))
        return(NULL)
      }
      run(.lrt(fit, red_fit, lab), "param", lab)
    }
    if (type == "mt" && test_rg) {
      fit0 <- run(fit1(dat, maxit = maxit, genetic_cov = "idh"), "param", paste(lab, "(idh model)"))
      lr <- lrt_test(fit0, "idh model")
      par_rows[[i]][, `:=`(rg_lrt = if (is.null(lr)) NA_real_ else lr[[2L]][1L],
                            rg_p   = if (is.null(lr)) NA_real_ else lr[[ncol(lr)]][1L])]
      rm(fit0)
    }
    if (type == "mt" && test_rp) {   # all covariances = 0: idh genetic, random and residual
      fit00 <- run(fit1(dat, maxit = maxit, genetic_cov = "idh", random_cov = "idh",
                        residual_cov = "idh"), "param", paste(lab, "(all-idh model)"))
      lr <- lrt_test(fit00, "all-idh model")
      par_rows[[i]][, `:=`(cov_df  = if (is.null(lr)) NA_real_ else lr[[1L]][1L],
                            cov_lrt = if (is.null(lr)) NA_real_ else lr[[2L]][1L],
                            cov_p   = if (is.null(lr)) NA_real_ else lr[[ncol(lr)]][1L])]
      rm(fit00)
    }
    if (nc && cv_vc != "estimate" && p_status[i] != "not_converged") {   # start values by Component name
      init <- run(fit1(dat, maxit = maxit, start_values = TRUE), "param", paste(lab, "(start values)"))
      if (!is.null(init)) vc_cv[[i]] <- run(.vc_table(init, vc, lab, fix = cv_vc == "fixed"), "param", lab)
      rm(init)
    }
    rm(fit); gc(verbose = FALSE)
  }

  # ---- 2-3. cross-validation ------------------------------------------------
  acc_rows <- pred_rows <- list()
  if (nc) {
  pidx <- if (length(keys)) {
    m <- P[, c(keys), with = FALSE][, .pi := .I]
    m[C[, c(keys), with = FALSE], on = keys, .pi]
  } else rep(1L, nrow(C))
  if (anyNA(pidx)) stop("Some rows of `cv` have no matching row in `param`.", call. = FALSE)

  tnames <- .trait_names(P, keys, type, resp, given = !is.null(traits))
  acc_rows <- pred_rows <- vector("list", nc)
  t_cv <- Sys.time()
  for (j in seq_len(nc)) {
    i <- pidx[j]
    lab <- .label(C[j], c_outer)
    t1 <- Sys.time()
    meta <- C[j, c_outer, with = FALSE]
    tn <- tnames[[i]]
    pred <- va <- fit <- NULL
    conv <- NA; iter <- NA_integer_; upd <- NA_integer_
    status <- if (p_status[i] == "failed") "param_failed"
              else if (p_status[i] == "not_converged") "param_not_converged"
              else if (cv_vc != "estimate" && is.null(vc_cv[[i]])) "no_start_values"
              else "ok"
    if (status == "ok") {
      fit <- run(fit1(C$train[[j]], maxit = if (cv_vc == "fixed") 1L else maxit,
                      vc = vc_cv[[i]], check_converge = cv_vc != "fixed"), "cv", lab)
      if (is.null(fit)) status <- "failed"
      else {
        iter <- .n_iter(fit)
        if (cv_vc != "fixed") {
          conv <- conv_ok(fit); upd <- n_upd(fit)
          if (!conv) {
            note("cv", lab, "not_converged", sprintf("not converged after %d update()%s", upd,
                 if (nonconv == "skip") "; fold skipped" else ""))
            if (nonconv == "skip") status <- "not_converged"
          }
        }
        if (status == "ok") {
          va <- run(.va(.vc_index(as.data.frame(summary(fit)$varcomp), type, id, resp,
                                  names(fit$vparameters)), resp), "cv", lab)
          pred <- run(.validate(fit, C$validate[[j]], id, resp, va, Ftab, cv_min_rel), "cv", lab)
          if (is.null(pred)) status <- "failed"
        }
      }
    }
    if (verbose) {
      done <- j; el <- as.numeric(difftime(Sys.time(), t_cv, units = "secs"))
      eta <- el / done * (nc - done)
      what <- if (status == "ok" && !is.null(fit)) fit_note(fit, t1)
              else if (!is.null(fit)) sprintf("%s (skipped)", fit_note(fit, t1))
              else toupper(gsub("_", " ", status))
      message(sprintf("[cv %d/%d %3.0f%%] %s | %s | %s elapsed, ~%s left",
                      j, nc, 100 * j / nc, lab, what,
                      .fmt_time(difftime(Sys.time(), t_start, units = "secs")), .fmt_time(eta)))
    }
    rm(fit)
    acc_rows[[j]] <- cbind(meta[rep(1L, length(resp))], status = status,
                           .accuracy(pred, resp, tn, h2_list[[i]], va, min_n),
                           cv_converge = conv, cv_updates = upd, cv_iter = iter)
    if (keep_pred && !is.null(pred)) {
      pred[, trait := tn[match(trait, resp)]]
      pred_rows[[j]] <- cbind(meta[rep(1L, nrow(pred))], pred)
    }
    if (j %% 20L == 0L) gc(verbose = FALSE)
  }
  }

  # ---- 4. assemble ----------------------------------------------------------
  accuracy <- NULL
  if (nc) {
    accuracy <- data.table::rbindlist(acc_rows, use.names = TRUE, fill = TRUE)
    # traits without any record in a validation set (e.g. the other herd in mt_bt)
    accuracy <- accuracy[is.na(n_obs) | n_obs > 0]
    if (neg_to_na)
      accuracy[, `:=`(pa = .neg_na(pa), accuracy = .neg_na(accuracy), slope = .neg_na(slope))]
  }
  log <- if (length(issues)) data.table::rbindlist(issues)
         else data.table::data.table(time = character(), stage = character(), label = character(),
                                     type = character(), message = character())
  out <- list(
    varcomp  = data.table::rbindlist(vc_rows, use.names = TRUE, fill = TRUE),
    vc_index = data.table::rbindlist(idx_rows, use.names = TRUE, fill = TRUE),
    param    = data.table::rbindlist(par_rows, use.names = TRUE, fill = TRUE),
    accuracy = accuracy,
    summary  = if (nc) .summarise(accuracy[status == "ok"], c(keys, grp_cols, "trait")),
    gebv     = if (keep_gebv) data.table::rbindlist(gebv_rows, use.names = TRUE, fill = TRUE),
    pred     = if (nc && keep_pred) data.table::rbindlist(pred_rows, use.names = TRUE, fill = TRUE),
    log      = log
  )
  out <- Filter(Negate(is.null), out)

  cv_st <- if (nc) table(factor(unique(accuracy[, c(c_outer, "status"), with = FALSE])$status))
           else NULL
  msg <- sprintf("Done in %s | param: %d ok, %d not converged, %d failed | cv: %s | %d log entr%s",
                 .fmt_time(difftime(Sys.time(), t_start, units = "secs")),
                 sum(p_status == "ok"), sum(grepl("not_converged", p_status)), sum(p_status == "failed"),
                 if (nc) paste(names(cv_st), cv_st, sep = " ", collapse = ", ") else "none (cv = NULL)",
                 nrow(log), if (nrow(log) == 1L) "y" else "ies")
  say(msg)
  bad <- log[type %in% c("error", "not_converged")]
  if (nrow(bad))
    warning(sprintf("%d fit(s) failed or did not converge (%s); see `$log`.", nrow(bad),
                    paste(unique(bad$type), collapse = ", ")), call. = FALSE)

  if (!is.null(path)) .write_results(out, path, verbose)
  structure(out, class = c("asr_gblup", "list"), type = type, cv_vc = cv_vc, n_cv = nc)
}

#' Run several GBLUP scenarios from `cv_scenarios()`
#'
#' @param sc Output of [mintyr::cv_scenarios()].
#' @param ginv Inverse genomic relationship matrix.
#' @param scenarios Scenarios to run. `"mt_bt"` (one trait as a different
#'   trait per group, two groups) is fitted with `residual_cov = "diag"` and
#'   without \sb{r}{p}, since no individual has records in both groups.
#' @param dir Folder for the workbooks `<scenario>_<tag>.xlsx`; `NULL` writes nothing.
#' @param tag Label used in the file names; defaults to the name of `ginv`.
#' @param ... Further arguments for [asr_gblup()] (e.g. `fixed`, `inbreeding`).
#' @return A named list of [asr_gblup()] results.
#' @export
#' @examples
#' \dontrun{
#' res <- asr_scenarios(sc, ginv, fixed = ~ hhys, inbreeding = inbreeding)
#' }
asr_scenarios <- function(sc, ginv, scenarios = c("st_sp", "st_mp", "mt_sp", "mt_mp"),
                          dir = ".", tag = deparse(substitute(ginv)), ...) {
  force(tag)
  miss <- setdiff(scenarios, names(sc))
  if (length(miss)) stop("Not in `sc`: ", paste(miss, collapse = ", "), call. = FALSE)
  dots <- list(...)
  res <- lapply(scenarios, function(s) {
    message("==== ", s, " ====")
    args <- c(list(param = sc[[s]]$param, cv = sc[[s]]$cv, ginv = ginv,
                   path = if (is.null(dir)) NULL else file.path(dir, paste0(s, "_", tag, ".xlsx"))),
              dots)
    if (s == "mt_bt") {   # one trait as a different trait in each group
      grp <- attr(sc$folds, "cv_by")
      if (is.null(grp)) grp <- setdiff(.outer_cols(data.table::as.data.table(sc$mt_bt$cv)),
                                       c("name", "rep", "fold"))
      args$traits <- as.character(unique(sc$mt_bt$cv[[grp[1L]]]))
      if (length(args$traits) != 2L)
        stop("`mt_bt` is supported for two groups only.", call. = FALSE)
      args$residual_cov <- "diag"   # no individual has records in both groups
      args$rp <- FALSE
      args$test_rp <- FALSE
    }
    do.call(asr_gblup, args)
  })
  stats::setNames(res, scenarios)
}

#' @export
print.asr_gblup <- function(x, ...) {
  cat(sprintf("<asr_gblup> %s model | %d parameter fits | %d CV fits | CV variances: %s\n",
              if (identical(attr(x, "type"), "mt")) "two-trait" else "single-trait",
              nrow(x$param), as.integer(attr(x, "n_cv")), attr(x, "cv_vc")))
  if (!is.null(x$log) && nrow(x$log))
    cat(sprintf("log: %s\n", paste(names(table(x$log$type)), table(x$log$type), collapse = ", ")))
  print(if (!is.null(x$summary)) x$summary else x$param)
  invisible(x)
}
