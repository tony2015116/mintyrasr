#' GBLUP parameter estimation and cross-validation with ASReml-R
#'
#' @description
#' `asr_gblup()` runs the whole GBLUP workflow for one scenario directly in R,
#' so no input files have to be exported for ASReml:
#'
#' 1. fits the model to every data set in `param`: variance components,
#'    heritability and, for two-trait models, genetic / phenotypic
#'    correlations plus an LRT of rg = 0 (`idh(trait)` vs the full genetic
#'    structure, `boundary = FALSE`);
#' 2. refits the model to every training set in `cv`, with the variance
#'    components handled as set by `cv_vc`;
#' 3. predicts GEBVs of the validation individuals, computes their
#'    reliability `1 - PEV / ((1 + F) * Va)`, drops GEBVs with reliability
#'    below `min_rel`, and returns per fold the predictive ability
#'    `pa = r(y, GEBV)`, the accuracy `pa / sqrt(h2)` (h2 from the full data,
#'    Dekkers 2007) and the slope `b(y ~ GEBV)`;
#' 4. summarises them over folds (n, mean, sd, se) and, optionally, writes
#'    everything to one Excel workbook.
#'
#' h2, rg and rp are built from the **names** of the variance components,
#' not from fixed `V` positions, so they stay correct when random effects
#' are added (`random = ~ ide(id) + litter`) or covariance structures change
#' (`us`, `idh`/`diag`, `corh`, `corgh`):
#' h2 = Va / (all variances of the trait), rg = Cov_a / sqrt(Va1 Va2)
#' (the correlation itself for `corh`/`corgh`), rp = sum of all covariances /
#' sqrt(Vp1 Vp2). The `vc_index` table shows which `V` each component is and
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
#'   Each row is matched to the row of `param` that shares its outer columns.
#' @param ginv Inverse genomic (or pedigree) relationship matrix as used by
#'   `vm()`, e.g. `G.inverse(G, sparseform = TRUE)$Ginv.sparse`. It must carry
#'   `attr(, "INVERSE") = TRUE`, otherwise `vm()` takes it as G itself
#'   (a warning is given).
#' @param fixed One-sided formula of fixed effects, e.g. `~ hhys`. For
#'   two-trait models every term is nested in trait (`trait + trait:hhys`).
#' @param random Optional one-sided formula of further random effects, e.g.
#'   `~ ide(id)` (permanent environment, repeated records), `~ litter + pen`.
#'   In two-trait models each term gets the structure `random_cov`.
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
#' @param min_rel GEBVs with reliability below this value are set to `NA`
#'   before computing PA / accuracy. Use `-Inf` to keep all.
#' @param neg_to_na Set negative PA / accuracy / slope to `NA` before the
#'   summary (the behaviour of the original scripts).
#' @param test_rg Two-trait models only: LRT of rg = 0 (`idh(trait)` genetic
#'   structure, other structures unchanged, `boundary = FALSE`).
#' @param test_rp Two-trait models only: LRT of *all* covariances = 0
#'   (`idh` for the genetic, extra random and residual terms; `df` = number
#'   of covariance / correlation parameters). This tests whether the traits
#'   are correlated at all, not rp = 0 as such.
#' @param rp Compute the phenotypic correlation. Set `FALSE` when no
#'   individual has records for both traits (e.g. one trait in two herds):
#'   the residual covariance is then not estimable and rp is not meaningful.
#' @param na_x `na.method(x = )`: `"include"` (default, as in the original
#'   scripts) treats missing covariates / factor levels as 0 in the design
#'   matrix - centre numeric covariates first; ASReml's own default is `"fail"`.
#' @param dense Put the genomic term in `dense` (faster for GBLUP).
#' @param maxit,workspace Passed to `asreml()` (maxit is 1 for `cv_vc = "fixed"`).
#'   ASReml's own defaults are 13 and `"128mb"`; two-trait GBLUP usually
#'   needs more.
#' @param max_update Number of `update()` calls on a fit that has not converged.
#' @param asreml_args Named list of further `asreml()` arguments, e.g.
#'   `list(ai.sing = TRUE)`.
#' @param on_error `"warn"` (default): a failing fit gives a warning and `NA`
#'   results, the loop continues. `"stop"`: stop at the first failure.
#' @param keep_pred Return the validation predictions (`pred`).
#' @param path Optional `.xlsx` file to write the results to.
#' @param verbose Print progress.
#'
#' @return A list of class `asr_gblup` with `varcomp`, `vc_index`, `param`,
#'   `accuracy` (one row per fold and trait), `summary` and, if `keep_pred`,
#'   `pred`. Fitted model objects are not kept, which keeps memory use low.
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
asr_gblup <- function(param, cv, ginv, fixed = ~ 1, random = NULL, id = "id", traits = NULL,
                      genetic_cov = "us", random_cov = "us", residual_cov = "us",
                      cv_vc = c("fixed", "warm", "estimate"), vpredict = NULL,
                      inbreeding = NULL, min_rel = 0.05, neg_to_na = TRUE,
                      test_rg = TRUE, test_rp = FALSE, rp = TRUE, na_x = c("include", "fail", "omit"),
                      dense = TRUE, maxit = 100L, workspace = "30Gb", max_update = 0L,
                      asreml_args = list(), on_error = c("warn", "stop"),
                      keep_pred = TRUE, path = NULL, verbose = TRUE) {
  on_error <- match.arg(on_error)
  cv_vc    <- match.arg(cv_vc)
  na_x     <- match.arg(na_x)
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
  C <- data.table::as.data.table(cv)
  if (!"data" %in% names(P)) stop("`param` needs a list-column `data`.", call. = FALSE)
  if (!all(c("train", "validate") %in% names(C)))
    stop("`cv` needs list-columns `train` and `validate`.", call. = FALSE)

  p_outer <- .outer_cols(P)
  c_outer <- .outer_cols(C)
  keys      <- intersect(p_outer, c_outer)               # e.g. name / pairs (+ herd)
  fold_cols <- intersect(c_outer, c("rep", "fold", "id", "id2"))
  grp_cols  <- setdiff(c_outer, c(keys, fold_cols))      # validated group (mp scenarios)
  if (!length(keys) && nrow(P) != 1L)
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
  if (is.null(Ftab) && verbose) message("`inbreeding` is NULL: F = 0 is used for all individuals.")
  run <- function(expr, what) .try_fit(expr, what, on_error)
  spec <- list(type = type, resp = resp, fixed = fixed, random = random, id = id,
               genetic_cov = genetic_cov, random_cov = random_cov, residual_cov = residual_cov,
               dense = dense, na_x = na_x, workspace = workspace, max_update = max_update,
               asreml_args = asreml_args)
  fit1 <- function(d, ...) .asr_fit(d, spec = spec, ginv = ginv, ...)

  # ---- 1. parameter estimation (full data) -----------------------------------
  np <- nrow(P)
  vc_rows <- par_rows <- idx_rows <- vc_cv <- h2_list <- vector("list", np)
  for (i in seq_len(np)) {
    lab <- .label(P[i], p_outer)
    if (verbose) message(sprintf("[param %d/%d] %s", i, np, lab))
    dat  <- P$data[[i]]
    meta <- P[i, p_outer, with = FALSE]
    fit  <- run(fit1(dat, maxit = maxit), lab)
    if (is.null(fit)) {
      par_rows[[i]] <- cbind(meta, converge = NA)
      next
    }
    vc  <- as.data.frame(summary(fit)$varcomp)
    idx <- .vc_index(vc, type, id, resp, names(fit$vparameters))
    bnd <- idx$item[idx$bound == "B"]
    if (length(bnd)) warning(sprintf("[%s] fixed at the boundary (B): %s", lab,
                                     paste(unique(bnd), collapse = "; ")), call. = FALSE)
    vc_rows[[i]]  <- cbind(meta[rep(1L, nrow(vc))], V = idx$V[match(rownames(vc), idx$item)], item = rownames(vc), vc)
    idx_rows[[i]] <- cbind(meta[rep(1L, nrow(idx))], idx[order(V), !"component"])
    gp <- .genetic_params(fit, c(.vp_formulas(idx, resp, rp = rp), vpredict))
    h2_list[[i]] <- stats::setNames(
      vapply(resp, function(r) {
        nm <- if (type == "st") "h2" else paste0("h2_", r)
        if (nm %in% names(gp)) gp[[nm]] else NA_real_
      }, numeric(1L)), resp)
    par_rows[[i]] <- cbind(meta, gp, converge = isTRUE(fit$converge),
                           loglik = if (is.null(fit$loglik)) NA_real_ else fit$loglik,
                           boundary = paste(unique(bnd), collapse = "; "))
    if (type == "mt" && test_rg) {
      fit0 <- run(fit1(dat, maxit = maxit, genetic_cov = "idh"), paste(lab, "(idh model)"))
      lr <- if (is.null(fit0)) NULL else
        tryCatch(as.data.frame(asreml::lrt.asreml(fit, fit0, boundary = FALSE)),
                 error = function(e) NULL)
      par_rows[[i]][, `:=`(rg_lrt = if (is.null(lr)) NA_real_ else lr[[2L]][1L],
                            rg_p   = if (is.null(lr)) NA_real_ else lr[[ncol(lr)]][1L])]
      rm(fit0)
    }
    if (type == "mt" && test_rp) {   # all covariances = 0: idh genetic, random and residual
      fit00 <- run(fit1(dat, maxit = maxit, genetic_cov = "idh", random_cov = "idh",
                        residual_cov = "idh"), paste(lab, "(all-idh model)"))
      lr <- if (is.null(fit00)) NULL else
        tryCatch(as.data.frame(asreml::lrt.asreml(fit, fit00, boundary = FALSE)),
                 error = function(e) NULL)
      par_rows[[i]][, `:=`(cov_df  = if (is.null(lr)) NA_real_ else lr[[1L]][1L],
                            cov_lrt = if (is.null(lr)) NA_real_ else lr[[2L]][1L],
                            cov_p   = if (is.null(lr)) NA_real_ else lr[[ncol(lr)]][1L])]
      rm(fit00)
    }
    if (cv_vc != "estimate") {   # start-value table filled by Component name
      init <- run(fit1(dat, maxit = maxit, start_values = TRUE), paste(lab, "(start values)"))
      if (!is.null(init)) vc_cv[[i]] <- .vc_table(init, vc, lab, fix = cv_vc == "fixed")
      rm(init)
    }
    rm(fit); gc(verbose = FALSE)
  }

  # ---- 2-3. cross-validation ------------------------------------------------
  pidx <- if (length(keys)) {
    m <- P[, c(keys), with = FALSE][, .pi := .I]
    m[C[, c(keys), with = FALSE], on = keys, .pi]
  } else rep(1L, nrow(C))
  if (anyNA(pidx)) stop("Some rows of `cv` have no matching row in `param`.", call. = FALSE)

  tnames <- .trait_names(P, keys, type, resp, given = !is.null(traits))
  nc <- nrow(C)
  acc_rows <- pred_rows <- vector("list", nc)
  for (j in seq_len(nc)) {
    i <- pidx[j]
    lab <- .label(C[j], c_outer)
    if (verbose) message(sprintf("[cv %d/%d] %s", j, nc, lab))
    meta <- C[j, c_outer, with = FALSE]
    tn <- tnames[[i]]
    pred <- va <- NULL
    if (cv_vc == "estimate" || !is.null(vc_cv[[i]])) {
      fit <- run(fit1(C$train[[j]],
                      maxit = if (cv_vc == "fixed") 1L else maxit,
                      vc = vc_cv[[i]], check_converge = cv_vc != "fixed"), lab)
      if (!is.null(fit)) {
        va <- .va(.vc_index(as.data.frame(summary(fit)$varcomp), type, id, resp,
                            names(fit$vparameters)), resp)
        pred <- .validate(fit, C$validate[[j]], id, resp, va, Ftab, min_rel)
        rm(fit); gc(verbose = FALSE)
      }
    }
    acc_rows[[j]] <- cbind(meta[rep(1L, length(resp))], .accuracy(pred, resp, tn, h2_list[[i]], va))
    if (keep_pred && !is.null(pred)) {
      pred[, trait := tn[match(trait, resp)]]
      pred_rows[[j]] <- cbind(meta[rep(1L, nrow(pred))], pred)
    }
  }

  # ---- 4. assemble ----------------------------------------------------------
  accuracy <- data.table::rbindlist(acc_rows, use.names = TRUE, fill = TRUE)
  # traits without any record in a validation set (e.g. the other herd in mt_bt)
  accuracy <- accuracy[is.na(n_obs) | n_obs > 0]
  if (neg_to_na)
    accuracy[, `:=`(pa = .neg_na(pa), accuracy = .neg_na(accuracy), slope = .neg_na(slope))]
  out <- list(
    varcomp  = data.table::rbindlist(vc_rows, use.names = TRUE, fill = TRUE),
    vc_index = data.table::rbindlist(idx_rows, use.names = TRUE, fill = TRUE),
    param    = data.table::rbindlist(par_rows, use.names = TRUE, fill = TRUE),
    accuracy = accuracy,
    summary  = .summarise(accuracy, c(keys, grp_cols, "trait")),
    pred     = if (keep_pred) data.table::rbindlist(pred_rows, use.names = TRUE, fill = TRUE)
  )
  out <- Filter(Negate(is.null), out)
  if (!is.null(path)) {
    if (!requireNamespace("writexl", quietly = TRUE))
      stop("Package 'writexl' is needed to write `path`.", call. = FALSE)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writexl::write_xlsx(lapply(out, as.data.frame), path)
    if (verbose) message("Written: ", path)
  }
  structure(out, class = c("asr_gblup", "list"), type = type, cv_vc = cv_vc, n_cv = nc)
}

#' Run several GBLUP scenarios from `cv_scenarios()`
#'
#' @param sc Output of [mintyr::cv_scenarios()].
#' @param ginv Inverse genomic relationship matrix.
#' @param scenarios Scenarios to run. `"mt_bt"` (one trait as a different
#'   trait per group, two groups) is fitted with `residual_cov = "diag"` and
#'   without rp, since no individual has records in both groups.
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
  print(x$summary)
  invisible(x)
}
