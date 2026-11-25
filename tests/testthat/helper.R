# Model tests need the CI mock of asreml (tools/mock-asreml); they are skipped
# with the real ASReml-R, whose numbers are not fixed.
skip_if_no_mock <- function() {
  skip_if_not_installed("asreml")
  skip_if_not_installed("mintyr")
  title <- utils::packageDescription("asreml")$Title
  skip_if_not(isTRUE(grepl("Mock", title)), "real ASReml-R installed: model tests use the CI mock")
}

sim_data <- function(n = 120, groups = c("A", "B"), seed = 42) {
  set.seed(seed)
  ids <- paste0("ID_", seq_len(n))                      # underscores on purpose
  tbv <- stats::setNames(stats::rnorm(n), ids)
  pheno <- data.frame(id = ids, test_herd = rep(groups, each = n / length(groups)),
                      litter = sample(paste0("L", 1:20), n, TRUE),
                      hhys = sample(paste0("h", 1:4), n, TRUE),
                      adg = tbv + stats::rnorm(n), bf = 0.5 * tbv + stats::rnorm(n))
  ginv <- diag(n); rownames(ginv) <- colnames(ginv) <- ids
  attr(ginv, "tbv") <- tbv; attr(ginv, "INVERSE") <- TRUE
  list(pheno = pheno, ginv = ginv, F = data.frame(id = ids, F = 0.01))
}

vc_df <- function(names, comp = seq_along(names) / 10, bound = "P") {
  data.frame(component = comp, std.error = 0.01, bound = bound, row.names = names)
}
