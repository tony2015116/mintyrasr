#' Compare relationship matrices: pedigree, genomic and single-step
#'
#' Runs the same scenarios with several inverse relationship matrices - e.g.
#' A^-1 (pedigree BLUP), G^-1 (GBLUP) and H^-1 (single-step GBLUP) - and
#' stacks the genetic parameters and the cross-validation summaries with a
#' `method` column, so the methods can be compared directly. All methods use
#' the same folds, so the comparison is paired.
#'
#' @param sc Output of [mintyr::cv_scenarios()].
#' @param inverses Named list of inverse relationship matrices, e.g.
#'   `list(ABLUP = ainv, GBLUP = ginv, ssGBLUP = hinv)`. The names become
#'   the `method` column and the file tags.
#' @param scenarios Scenarios to run (see [asr_scenarios()]).
#' @param dir Folder for the workbooks: one per scenario and method
#'   (`<scenario>_<method>.xlsx`) plus `compare.xlsx`; `NULL` writes nothing.
#' @param inbreeding Optional named list of inbreeding tables, one per method
#'   (same names as `inverses`); missing entries use the defaults of
#'   [asr_gblup()] (e.g. F from `ainverse()`).
#' @param ... Further arguments for [asr_gblup()], e.g. `fixed = ~ hhys`.
#'
#' @return A list with `param` and `summary` (all methods and scenarios,
#'   with `method` and `scenario` columns) and `results`, the
#'   [asr_scenarios()] output of every method.
#' @export
#' @examples
#' \dontrun{
#' ainv <- ainverse(ped)
#' ginv <- G.inverse(G, sparseform = TRUE)$Ginv.sparse
#' cmp  <- asr_compare(sc, list(ABLUP = ainv, GBLUP = ginv),
#'                     scenarios = c("st_sp", "st_mp"), fixed = ~ hhys)
#' cmp$summary[type == "pa"]
#' }
asr_compare <- function(sc, inverses, scenarios = c("st_sp", "st_mp", "mt_sp", "mt_mp"),
                        dir = ".", inbreeding = NULL, ...) {
  if (!is.list(inverses) || is.null(names(inverses)) || any(!nzchar(names(inverses))))
    stop("`inverses` must be a named list, e.g. list(ABLUP = ainv, GBLUP = ginv).", call. = FALSE)
  res <- lapply(names(inverses), function(m) {
    message("######## method: ", m, " ########")
    args <- list(sc = sc, ginv = inverses[[m]], scenarios = scenarios, dir = dir, tag = m, ...)
    if (!is.null(inbreeding[[m]])) args$inbreeding <- inbreeding[[m]]
    do.call(asr_scenarios, args)
  })
  names(res) <- names(inverses)
  stack <- function(what) data.table::rbindlist(unlist(lapply(names(res), function(m)
    lapply(names(res[[m]]), function(s) {
      x <- res[[m]][[s]][[what]]
      if (is.null(x) || !nrow(x)) return(NULL)
      cbind(method = m, scenario = s, x)
    })), recursive = FALSE), use.names = TRUE, fill = TRUE)
  out <- list(param = stack("param"), summary = stack("summary"), results = res)
  if (!is.null(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    .write_results(out[c("param", "summary")], file.path(dir, "compare.xlsx"))
  }
  out
}
