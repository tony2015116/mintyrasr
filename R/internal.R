#' @import data.table
NULL

# Internal helpers for asr_gblup() ------------------------------------------

.outer_cols <- function(dt) names(dt)[!vapply(dt, is.list, logical(1L))]

.label <- function(row, cols) {
  if (!length(cols)) return("all")
  paste(paste0(cols, "=", vapply(cols, function(c) as.character(row[[c]][1L]), "")),
        collapse = " ")
}

# Elapsed time as "1h02m", "3m05s" or "12s"
.fmt_time <- function(sec) {
  if (is.null(sec) || !is.finite(sec)) return("--")
  sec <- round(sec)
  if (sec >= 3600) sprintf("%dh%02dm", sec %/% 3600, (sec %% 3600) %/% 60)
  else if (sec >= 60) sprintf("%dm%02ds", sec %/% 60, sec %% 60)
  else sprintf("%ds", sec)
}

# Progress line: [cv 12/60  20%] label | elapsed 3m05s | ETA 12m20s
.progress <- function(stage, k, n, t0, lab) {
  el  <- proc.time()[["elapsed"]] - t0
  eta <- if (k > 1L) el / (k - 1L) * (n - k + 1L) else NA_real_
  message(sprintf("[%s %d/%d %3.0f%%] %s | elapsed %s | ETA %s",
                  stage, k, n, 100 * (k - 1L) / n, lab, .fmt_time(el), .fmt_time(eta)))
}

# Run log: one row per fit (or skipped fit) with status and messages
.new_log <- function() {
  e <- new.env(parent = emptyenv())
  e$rows <- list()
  e$add <- function(stage, label, status, n_update = NA_integer_, seconds = NA_real_, message = "") {
    e$rows[[length(e$rows) + 1L]] <- data.table::data.table(
      stage = stage, label = label, status = status, n_update = as.integer(n_update),
      seconds = round(as.numeric(seconds), 2), message = message)
    invisible(NULL)
  }
  e$table <- function() {
    if (!length(e$rows)) return(data.table::data.table(stage = character(), label = character(),
      status = character(), n_update = integer(), seconds = numeric(), message = character()))
    data.table::rbindlist(e$rows)
  }
  e
}

# Run one fit safely: errors are caught (on_error), ASReml's own warnings are
# collected into the log instead of being printed one by one (muffle: regex
# of expected warnings to drop, e.g. "not converged" when all variances are
# fixed and maxit = 1), and the result carries attr "status":
# "ok" / "not_converged". Returns NULL on error.
.safe_fit <- function(expr, stage, lab, log, on_error, expect_converge = TRUE, muffle = NULL) {
  t1 <- proc.time()[["elapsed"]]
  warns <- character()
  fit <- tryCatch(
    withCallingHandlers(expr, warning = function(w) {
      m <- conditionMessage(w)
      if (is.null(muffle) || !grepl(muffle, m, ignore.case = TRUE)) warns <<- c(warns, m)
      invokeRestart("muffleWarning")
    }),
    error = function(e) e)
  secs <- proc.time()[["elapsed"]] - t1
  if (inherits(fit, "error")) {
    msg <- sprintf("ASReml fit failed for %s [%s]: %s", lab, stage, conditionMessage(fit))
    log$add(stage, lab, "failed", NA_integer_, secs, msg)
    if (on_error == "stop") stop(msg, call. = FALSE)
    warning(msg, call. = FALSE)
    return(NULL)
  }
  nu <- attr(fit, "n_update"); if (is.null(nu)) nu <- 0L
  status <- if (!expect_converge || isTRUE(fit$converge)) "ok" else "not_converged"
  log$add(stage, lab, status, nu, secs, paste(unique(warns), collapse = " | "))
  attr(fit, "status") <- status
  fit
}

# Likelihood-ratio test full vs reduced (REML, same fixed effects).
# ASReml-R 4.2 exports the generic lrt(); older versions lrt.asreml().
# If neither works, LR = 2 (logL_full - logL_reduced) with df = difference in
# the number of non-fixed variance parameters.
.asr_fun <- function(nms) {
  for (n in nms) {
    f <- tryCatch(getExportedValue("asreml", n), error = function(e) NULL)
    if (is.function(f)) return(f)
  }
  NULL
}
.n_free <- function(fit) {
  b <- as.character(summary(fit)$varcomp$bound)
  sum(!b %in% c("F", "Fixed"))
}
.lrt <- function(full, reduced, lab) {
  f <- .asr_fun(c("lrt", "lrt.asreml"))
  r <- if (is.null(f)) NULL else
    tryCatch(as.data.frame(f(full, reduced, boundary = FALSE)), error = function(e) e)
  if (is.data.frame(r) && nrow(r))
    return(list(df = r[[1L]][1L], lr = r[[2L]][1L], p = r[[ncol(r)]][1L]))
  why <- if (inherits(r, "error")) conditionMessage(r) else "asreml::lrt() not found"
  df <- .n_free(full) - .n_free(reduced)
  lr <- 2 * (full$loglik - reduced$loglik)
  if (!is.finite(lr) || df < 1L) {
    warning(sprintf("[%s] LRT not possible (%s).", lab, why), call. = FALSE)
    return(list(df = NA_real_, lr = NA_real_, p = NA_real_))
  }
  warning(sprintf("[%s] asreml LRT failed (%s); LR computed from the log-likelihoods (df = %d).",
                  lab, why, df), call. = FALSE)
  list(df = df, lr = lr, p = stats::pchisq(max(lr, 0), df, lower.tail = FALSE))
}

# data.frame with factors, as asreml() expects
.asr_df <- function(d, id) {
  d <- as.data.frame(d)
  for (nm in names(d)) if (is.character(d[[nm]])) d[[nm]] <- factor(d[[nm]])
  if (!is.factor(d[[id]])) d[[id]] <- factor(d[[id]])
  d
}

# One ASReml fit. Formulas are built here so that asreml finds `ginv`, `vc`
# ... in the formula / calling environment, and update() works.
# spec: type, fixed, random, id, genetic_cov, random_cov, residual_cov,
#       dense, na_x, workspace, max_update, asreml_args, trace, verbose
.asr_fit <- function(d, spec, ginv, maxit, vc = NULL, start_values = FALSE,
                     genetic_cov = spec$genetic_cov, random_cov = spec$random_cov,
                     residual_cov = spec$residual_cov, check_converge = TRUE) {
  id <- spec$id
  d  <- .asr_df(d, id)
  ft <- attr(stats::terms(spec$fixed), "term.labels")
  rt <- if (is.null(spec$random)) character(0) else attr(stats::terms(spec$random), "term.labels")
  env <- environment()
  if (spec$type == "st") {
    fx <- paste(spec$resp, "~", if (length(ft)) paste(ft, collapse = " + ") else "1")
    gt <- sprintf("vm(%s, ginv)", id)
    rd <- paste("~", paste(c(gt, rt), collapse = " + "))
    rs <- "~ idv(units)"
  } else {
    fx <- paste0("cbind(", paste(spec$resp, collapse = ", "), ") ~ trait",
                if (length(ft)) paste0("+ ", paste0("trait:", ft, collapse = " + ")) else "")
    gt <- sprintf("%s(trait):vm(%s, ginv)", genetic_cov, id)
    ot <- if (length(rt)) paste0(random_cov, "(trait):", rt) else character(0)
    rd <- paste("~", paste(c(gt, ot), collapse = " + "))
    rs <- sprintf("~ units:%s(trait)", residual_cov)
  }
  dn <- stats::as.formula(paste("~", gt), env = env)   # only the genomic term goes to `dense`
  fx <- stats::as.formula(fx, env = env)
  rd <- stats::as.formula(rd, env = env)
  rs <- stats::as.formula(rs, env = env)
  nam <- asreml::na.method(x = spec$na_x, y = "include")
  workspace <- spec$workspace

  cl <- quote(asreml::asreml(fixed = fx, random = rd, residual = rs, na.action = nam,
                             data = d, maxit = maxit, workspace = workspace))
  if (spec$dense) cl$dense <- quote(dn)
  if (!is.null(vc)) { cl$G.param <- quote(vc); cl$R.param <- quote(vc) }
  if (!is.null(spec$trace)) cl$trace <- spec$trace
  for (nm in names(spec$asreml_args)) cl[[nm]] <- spec$asreml_args[[nm]]
  if (start_values) { cl$start.values <- TRUE; return(eval(cl)) }
  fit <- eval(cl)
  attr(fit, "n_update") <- 0L
  if (!check_converge) return(fit)          # all variances fixed, maxit = 1: BLUP only
  # not converged: restart from the current estimates with update()
  k <- 0L
  while (!isTRUE(fit$converge) && k < spec$max_update) {
    k <- k + 1L
    if (isTRUE(spec$verbose))
      message(sprintf("    not converged: update() %d/%d", k, spec$max_update))
    fit <- stats::update(fit)
  }
  attr(fit, "n_update") <- k
  fit
}

# Start-value table from start.values = TRUE, with Value filled in by
# Component *name* from the full-data estimates; fix = TRUE sets all "F".
.vc_table <- function(init, vc, lab, fix) {
  tab <- as.data.frame(init$vparameters.table)
  val <- vc$component[match(tab$Component, rownames(vc))]
  miss <- is.na(val)
  if (any(miss))
    warning(sprintf("[%s] %d variance parameter(s) not matched to the estimates by name: %s",
                    lab, sum(miss), paste(tab$Component[miss], collapse = "; ")), call. = FALSE)
  tab$Value[!miss] <- val[!miss]
  if (fix) tab$Constraint <- "F"
  tab
}

# Regex for the genomic term vm(<id>, ginv) in component / coefficient names
.gpat <- function(id) paste0("vm\\(", id, ",\\s*ginv\\)")

# Map every row of summary(fit)$varcomp (= V1, V2, ... in vpredict) to its
# term, role (genetic / residual / random) and kind (var / cov / cor) for
# which trait(s). Trait names are found as tokens in the name, so us(),
# diag()/idh(), corh() and corgh() rows are all recognised, e.g.
#   "vm(id, ginv)", "ide(id)", "litter", "units!units", "units!R"  (single-trait)
#   "trait:vm(id, ginv)!trait_value2:value1"   us() covariance
#   "trait:litter!trait_value1"                 idh()/diag() variance
#   "...!trait!cor", "...value2...value1...cor"  corh() / corgh() correlation
#   "units:trait!R", "units!R" (with units!units)  scale fixed at 1 -> not used
# fa() / rr() terms are not interpreted (used = FALSE).
.vc_index <- function(vc, type, id, resp, vnames = NULL) {
  rn   <- rownames(vc)
  term <- sub("!.*$", "", rn)
  suf  <- ifelse(grepl("!", rn, fixed = TRUE), sub("^[^!]*!", "", rn), "")
  role <- ifelse(grepl(.gpat(id), term), "genetic",
                 ifelse(grepl("^units", term), "residual", "random"))
  bound <- if ("bound" %in% names(vc)) as.character(vc$bound) else rep("", length(rn))
  scale <- suf == "R" & (type == "mt" | bound %in% c("F", "Fixed") | "units!units" %in% rn)
  is_cor <- grepl("cor", suf, fixed = TRUE)
  fa_rr  <- grepl("(^|:)(fa|rr)\\(", term)
  n <- length(rn)
  kind <- rep(NA_character_, n); t1 <- t2 <- rep(NA_character_, n)
  if (type == "st") {
    kind <- ifelse(is_cor, "cor", "var"); t1 <- t2 <- rep(resp, n)
  } else {
    esc <- gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", resp)
    esc <- esc[order(-nchar(resp))]          # longest name first
    pat <- paste0("(?<![A-Za-z0-9])(", paste(esc, collapse = "|"), ")(?![A-Za-z0-9])")
    tok <- regmatches(suf, gregexpr(pat, suf, perl = TRUE))
    for (k in seq_len(n)) {
      x <- unique(tok[[k]])
      if (length(x) == 1L && !is_cor[k]) { kind[k] <- "var"; t1[k] <- t2[k] <- x }
      else if (length(x) == 2L) { kind[k] <- if (is_cor[k]) "cor" else "cov"; t1[k] <- x[1]; t2[k] <- x[2] }
      else if (!length(x) && is_cor[k]) { kind[k] <- "cor"; t1[k] <- resp[1]; t2[k] <- resp[2] }
      else if (!length(x) && grepl("(^|!)var$", suf[k])) kind[k] <- "hvar"   # corv(): one variance for both traits
    }
  }
  used <- !is.na(kind) & !scale & !fa_rr
  # corv(): the common variance is the variance of each trait
  if (any(kind == "hvar", na.rm = TRUE)) {
    h <- which(kind == "hvar")
    kind[h] <- "var"; t1[h] <- resp[1]; t2[h] <- resp[1]
    extra <- data.table::data.table(k = h)
  } else extra <- NULL
  # V numbers follow the order of fit$vparameters (what vpredict() uses)
  V <- seq_len(n)
  if (!is.null(vnames)) {
    m <- match(rn, vnames)
    if (anyNA(m)) warning("Variance components could not all be matched to fit$vparameters; ",
                          "V numbers follow summary()$varcomp order.", call. = FALSE)
    else V <- m
  }
  if (any(fa_rr)) warning("fa()/rr() terms are not interpreted for h2/rg/rp; ",
                          "give your own formulas in `vpredict`.", call. = FALSE)
  out <- data.table::data.table(V = V, item = rn, term = term, role = role, kind = kind,
                                trait1 = t1, trait2 = t2, used = used,
                                bound = bound, component = vc$component)
  if (!is.null(extra)) {                       # duplicate corv variance row for trait 2
    dup <- out[extra$k]; dup[, `:=`(trait1 = resp[2], trait2 = resp[2])]
    out <- rbind(out, dup)
  }
  out
}

# Genetic variance per response column, from the index table
.va <- function(idx, resp) {
  v <- vapply(resp, function(r) {
    hit <- idx$used & idx$role == "genetic" & idx$kind == "var" & idx$trait1 == r
    if (sum(hit) == 1L) idx$component[hit] else NA_real_
  }, numeric(1L))
  if (anyNA(v)) warning("Genetic variance not found for: ", paste(resp[is.na(v)], collapse = ", "),
                        call. = FALSE)
  stats::setNames(v, resp)
}

# vpredict() formulas built from the index table, so they follow the model:
#   h2 = Va / (all variances of the trait: genetic + other random + residual)
#   rg = Cov_a / sqrt(Va1 Va2)        (corh/corgh: the correlation itself)
#   rp = sum of all covariances / sqrt(Vp1 Vp2), covariance of a cor() term
#        written as r * sqrt(V1 V2)
.vp_formulas <- function(idx, resp, rp = TRUE) {
  u <- idx[idx$used, ]
  V <- function(k) paste0("V", k)
  var_rows <- function(r, roles = c("genetic", "residual", "random"), terms = NULL)
    sort(u$V[u$kind == "var" & u$trait1 == r & u$role %in% roles & (is.null(terms) | u$term %in% terms)])
  cov_expr <- function(k) {                 # covariance of row k as text
    row <- u[u$V == k, ][1L, ]
    if (row$kind == "cov") return(V(k))
    vi <- var_rows(row$trait1, terms = row$term); vj <- var_rows(row$trait2, terms = row$term)
    if (length(vi) != 1L || length(vj) != 1L) return(NA_character_)
    sprintf("%s * sqrt(%s * %s)", V(k), V(vi), V(vj))
  }
  mk <- function(lhs, txt) stats::as.formula(paste(lhs, "~", txt))
  f <- list()
  for (r in resp) {
    g <- var_rows(r, "genetic"); p <- var_rows(r)
    nm <- if (length(resp) == 1L) "h2" else paste0("h2_", r)
    if (length(g) == 1L) f[[nm]] <- mk(nm, sprintf("%s / (%s)", V(g), paste(V(p), collapse = " + ")))
  }
  if (length(resp) == 2L) {
    g1 <- var_rows(resp[1], "genetic"); g2 <- var_rows(resp[2], "genetic")
    gc <- u$V[u$role == "genetic" & u$kind %in% c("cov", "cor")]
    if (length(gc) == 1L && length(g1) == 1L && length(g2) == 1L)
      f$rg <- if (u$kind[u$V == gc][1L] == "cor") mk("rg", V(gc))
              else mk("rg", sprintf("%s / sqrt(%s * %s)", V(gc), V(g1), V(g2)))
    pc <- vapply(sort(unique(u$V[u$kind %in% c("cov", "cor")])), cov_expr, "")
    if (rp && length(pc) && !anyNA(pc))
      f$rp <- mk("rp", sprintf("(%s) / sqrt((%s) * (%s))", paste(pc, collapse = " + "),
                               paste(V(var_rows(resp[1])), collapse = " + "),
                               paste(V(var_rows(resp[2])), collapse = " + ")))
  }
  f
}

# Run the formulas (+ user formulas) through vpredict(); one row, est + se columns
.genetic_params <- function(fit, formulas) {
  out <- list()
  for (nm in names(formulas)) {
    r <- tryCatch(asreml::vpredict(fit, formulas[[nm]]), error = function(e) NULL)
    out[[nm]] <- if (is.null(r)) NA_real_ else r$Estimate[1L]
    out[[paste0(nm, "_se")]] <- if (is.null(r)) NA_real_ else r$SE[1L]
  }
  out$formulas <- paste(names(formulas), vapply(formulas, function(f) deparse(f[[3L]]), ""),
                        sep = " = ", collapse = "; ")
  data.table::as.data.table(out)
}

# vm() treats a matrix as an inverse only with attr(, "INVERSE") = TRUE
.check_ginv <- function(ginv) {
  inv <- attr(ginv, "INVERSE")
  if (!isTRUE(inv))
    warning("`ginv` has no attribute INVERSE = TRUE, so vm() treats it as a relationship ",
            "matrix (G), not its inverse. Use G.inverse(..., sparseform = TRUE)$Ginv.sparse or ",
            "set attr(ginv, \"INVERSE\") <- TRUE.", call. = FALSE)
  invisible(NULL)
}

# IDs of the relationship inverse (dense matrix or ASRgenomics/ainverse triplet)
.ginv_ids <- function(ginv) {
  ids <- attr(ginv, "rowNames")
  if (is.null(ids)) ids <- rownames(ginv)
  if (is.null(ids)) NULL else as.character(ids)
}

# GEBV, PEV and reliability of the validation individuals
.validate <- function(fit, val, id, resp, va, Ftab, min_rel) {
  cr <- as.data.frame(summary(fit, coef = TRUE)$coef.random)
  rn <- rownames(cr)
  g  <- grepl(paste0("^(trait_[^:]+:)?", .gpat(id), "_"), rn)   # genomic term only
  cr <- cr[g, , drop = FALSE]; rn <- rn[g]
  tr <- if (length(resp) == 1L) rep(resp, length(rn))
        else sub("^trait_([^:]+):.*$", "\\1", rn)
  ids <- sub(paste0("^(trait_[^:]+:)?", .gpat(id), "_"), "", rn)   # IDs may contain "_"
  gebv <- data.table::data.table(id = ids, trait = tr,
                                 gebv = cr$solution, pev = cr$std.error^2)
  v <- data.table::as.data.table(val)
  obs <- data.table::melt(
    data.table::data.table(id = as.character(v[[id]]), v[, resp, with = FALSE]),
    id.vars = "id", variable.name = "trait", value.name = "obs", variable.factor = FALSE)
  out <- gebv[obs, on = c("id", "trait")]
  out[, inb := if (is.null(Ftab)) 0 else Ftab$F[match(id, Ftab$id)]]
  out[, rel := 1 - pev / ((1 + inb) * va[trait])]
  if (is.finite(min_rel)) {
    had <- out[!is.na(gebv) & !is.na(obs), .N, by = trait]
    out[is.na(rel) | rel < min_rel, gebv := NA_real_]
    left <- out[!is.na(gebv) & !is.na(obs), .N, by = trait]
    gone <- setdiff(had$trait[had$N > 0], left$trait[left$N >= 3L])
    if (length(gone))
      warning(sprintf("min_rel = %g leaves fewer than 3 GEBVs for %s in this fold; PA is NA. ",
                      min_rel, paste(gone, collapse = ", ")),
              "Validation individuals have no own record, so low reliabilities are normal; ",
              "consider min_rel = -Inf.", call. = FALSE)
  }
  out[, c("id", "trait", "obs", "gebv", "pev", "inb", "rel"), with = FALSE]
}

# Predictive ability PA = r(y, GEBV); accuracy = PA / sqrt(h2) with the
# full-data h2 (Dekkers 2007); slope b(y ~ GEBV), ideally close to 1.
.accuracy <- function(pred, resp, tn, h2, va, status = "ok") {
  res <- lapply(resp, function(r) {
    na <- list(n = 0L, pa = NA_real_, slope = NA_real_,
               n_obs = if (is.null(pred)) NA_integer_ else sum(pred$trait == r & !is.na(pred$obs)))
    if (is.null(pred)) return(na)
    p  <- pred[trait == r & !is.na(obs) & !is.na(gebv)]
    n  <- nrow(p)
    if (n < 3L || stats::var(p$gebv) == 0) return(utils::modifyList(na, list(n = n)))
    list(n = n, pa = stats::cor(p$obs, p$gebv),
         slope = unname(stats::coef(stats::lm(obs ~ gebv, data = p))[2L]), n_obs = na$n_obs)
  })
  pa <- vapply(res, `[[`, numeric(1L), "pa")
  h  <- if (is.null(h2)) rep(NA_real_, length(resp)) else unname(h2[resp])
  data.table::data.table(trait = tn,
                         n = vapply(res, `[[`, integer(1L), "n"),
                         pa = pa,
                         accuracy = ifelse(!is.na(h) & h > 0, pa / sqrt(h), NA_real_),
                         slope = vapply(res, `[[`, numeric(1L), "slope"),
                         va = if (is.null(va)) NA_real_ else unname(va[resp]),
                         h2_full = h,
                         n_obs = vapply(res, `[[`, integer(1L), "n_obs"),
                         status = status)
}

.neg_na <- function(x) { x[!is.na(x) & x < 0] <- NA_real_; x }

.summarise <- function(acc, by) {
  long <- data.table::melt(acc, measure.vars = c("pa", "accuracy", "slope"),
                           variable.name = "type", value.name = "value",
                           variable.factor = FALSE)
  long[, list(n = sum(!is.na(value)),
              mean = if (all(is.na(value))) NA_real_ else mean(value, na.rm = TRUE),
              sd = stats::sd(value, na.rm = TRUE),
              se = stats::sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value)))),
       by = c(intersect(by, names(long)), "type")]
}

# Real trait names per param row (from `name` or `pairs = "a-b"`)
.trait_names <- function(P, keys, type, resp, given = FALSE) {
  lapply(seq_len(nrow(P)), function(i) {
    if (given) return(resp)
    if (type == "st" && "name" %in% names(P)) return(as.character(P$name[i]))
    if (type == "mt" && "pairs" %in% names(P)) {
      s <- strsplit(as.character(P$pairs[i]), "-", fixed = TRUE)[[1L]]
      if (length(s) == 2L) return(s)
    }
    resp
  })
}

.prep_F <- function(inbreeding) {
  if (is.null(inbreeding)) return(NULL)
  if (is.numeric(inbreeding) && !is.null(names(inbreeding)))
    return(data.table::data.table(id = names(inbreeding), F = unname(inbreeding)))
  if (is.data.frame(inbreeding) && ncol(inbreeding) >= 2L)
    return(data.table::data.table(id = as.character(inbreeding[[1L]]),
                                  F = as.numeric(inbreeding[[2L]])))
  stop("`inbreeding` must be a data.frame (ID, F) or a named numeric vector.", call. = FALSE)
}

utils::globalVariables(c(".pi", "trait", "obs", "gebv", "pev", "inb", "rel", "id",
                         "value", "pa", "accuracy", "slope", "n_obs", "V", ".", "status", "stage"))
