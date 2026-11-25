na.method <- function(x = "include", y = "include") list(x = x, y = y)
.vcnames <- function(mt, random, residual, rn = c("value1", "value2")) {
  tl <- attr(stats::terms(random), "term.labels")
  rs <- deparse(residual)
  if (!mt) return(c(tl, "units!units", "units!R"))     # idv(units): scale row at the end
  el <- function(term) {
    str  <- sub("\\(trait\\):.*$", "", term)
    base <- sub("^[a-z]+\\(trait\\):", "trait:", term)
    switch(str,
      us    = paste0(base, "!trait_", c(paste0(rn[1],":",rn[1]), paste0(rn[2],":",rn[1]), paste0(rn[2],":",rn[2]))),
      idh   = , diag = paste0(base, "!trait_", rn),
      corv  = paste0(base, c("!trait!cor", "!trait!var")),
      corh  = c(paste0(base, "!trait!cor"), paste0(base, "!trait_", rn)),
      corgh = c(paste0(base, "!trait!", rn[2], ":!trait!", rn[1], ".cor"), paste0(base, "!trait_", rn)))
  }
  G <- unlist(lapply(tl, el))
  R <- if (grepl("diag|idh", rs)) paste0("units:trait!trait_", rn)
       else paste0("units:trait!trait_", c(paste0(rn[1],":",rn[1]), paste0(rn[2],":",rn[1]), paste0(rn[2],":",rn[2])))
  c(G, "units:trait!R", R)
}
.comp <- function(nm) ifelse(grepl("!R$", nm), 1, ifelse(grepl("cor$", nm), 0.5,
                      ifelse(grepl("_([^:!]+):(?!\\1)[^:!]+$", nm, perl = TRUE), 0.1,
                      ifelse(grepl("units", nm), 0.6, ifelse(grepl("vm", nm), 0.3, 0.2)))))
asreml <- function(fixed, random, residual, dense = NULL, na.action, data, G.param = NULL,
                   R.param = NULL, start.values = FALSE, maxit = 30, workspace = NULL, trace = TRUE, ...) {
  stopifnot(is.data.frame(data), is.factor(data$id))
  for (v in all.vars(fixed)) if (!v %in% c(names(data), "trait")) stop("object '", v, "' not found")
  rtxt <- deparse(random)
  gname <- sub(".*vm\\(id, ([^)]+)\\).*", "\\1", rtxt)
  ginv <- get(gname, envir = environment(random))
  mt <- grepl("cbind", deparse(fixed[[2]]))
  diagG <- grepl("(diag|idh)\\(trait\\):vm", rtxt)
  rn <- if (mt) all.vars(fixed[[2]]) else NULL
  nm <- .vcnames(mt, random, residual, rn)
  if (start.values)
    return(list(vparameters.table = data.frame(Component = nm, Value = 0.1, Constraint = "P")))
  comp <- .comp(nm)
  if (!is.null(G.param)) {
    if (!all(nm %in% G.param$Component)) stop("G.param does not match the model")
    if (all(G.param$Constraint == "F") && maxit != 1) stop("mock: expected maxit = 1 for fixed")
    comp <- G.param$Value[match(nm, G.param$Component)]
    if (!all(G.param$Constraint == "F")) comp <- comp * 1.05   # re-estimated
  }
  vc <- data.frame(component = comp, std.error = 0.05, z.ratio = comp/0.05, bound = ifelse(comp == 1 & grepl("!R$", nm), "F", "P"), `%ch` = 0,
                   row.names = nm, check.names = FALSE)
  ids <- rownames(ginv)
  tbv <- if ("tbv" %in% names(data)) data$tbv[match(ids, as.character(data$id))] else NULL
  mk <- function(pref) {
    sol <- stats::rnorm(length(ids), 0, 0.3)
    known <- ids %in% as.character(data$id)
    sol[known] <- sol[known] + 1   # irrelevant offset
    # make predictions correlate with phenotypes stored in attribute
    tb <- attr(ginv, "tbv"); if (!is.null(tb)) sol <- tb[ids] + stats::rnorm(length(ids), 0, 0.5)
    data.frame(solution = sol, std.error = stats::runif(length(ids), 0.2, 0.5),
               z.ratio = 0, row.names = paste0(pref, sprintf("vm(id, %s)_", gname), ids))
  }
  cr <- if (!mt) mk("") else rbind(mk(paste0("trait_", rn[1], ":")), mk(paste0("trait_", rn[2], ":")))
  vp <- setNames(comp, nm)
  if (isTRUE(getOption("mock.shuffle"))) { vc <- vc[rev(seq_len(nrow(vc))), ] }   # varcomp order != vparameters
  if (isTRUE(getOption("mock.boundary")) && mt) vc$bound[grepl("litter", rownames(vc))] <- "B"
  if (isTRUE(getOption("mock.fail")) && is.null(G.param)) stop("mock: singular information matrix")
  if (!is.null(getOption("mock.warn")) && maxit == 1) warning(getOption("mock.warn"))
  conv <- !isTRUE(getOption("mock.nonconv")) || (!is.null(G.param) && all(G.param$Constraint == "F"))
  structure(list(converge = conv, vc = vc, cr = cr, vparameters = vp, loglik = -100 - diagG * 3), class = "asreml")
}
summary.asreml <- function(object, coef = FALSE, ...) {
  out <- list(varcomp = object$vc)
  if (coef) out$coef.random <- as.matrix(object$cr)
  out
}
update.asreml <- function(object, ...) {
  if (isTRUE(getOption("mock.update_fixes"))) object$converge <- TRUE
  object
}
vpredict <- function(object, transform) {
  V <- as.list(unname(object$vparameters)); names(V) <- paste0("V", seq_along(V))
  est <- eval(transform[[3]], V)
  if (!is.finite(est)) stop("mock vpredict: V index out of range")
  data.frame(Estimate = est, SE = 0.05, row.names = as.character(transform[[2]]))
}
lrt <- function(object, ..., boundary = TRUE) UseMethod("lrt")
lrt.asreml <- function(object, ..., boundary = TRUE) {
  if (boundary) stop("mock: rg test must use boundary = FALSE")
  o <- list(...)[[1]]; lr <- 2 * (object$loglik - o$loglik)
  data.frame(df = 1, `LR-statistic` = lr, `Pr(Chisq)` = stats::pchisq(lr, 1, lower.tail = FALSE),
             check.names = FALSE, row.names = "test")
}
asreml.options <- function(...) list(gammaPar = FALSE, maxit = 13, workspace = "128mb")
