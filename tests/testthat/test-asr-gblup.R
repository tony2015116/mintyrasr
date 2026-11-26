test_that("all four scenarios run and write workbooks", {
  skip_if_no_mock()
  skip_if_not_installed("writexl")
  s <- sim_data()
  sc <- mintyr::cv_scenarios(s$pheno[, c("id", "test_herd", "hhys", "adg", "bf")],
                             cols = c("adg", "bf"), id = "id", group = "test_herd",
                             v = 3, seed = 1, scenarios = c("st_sp", "st_mp", "mt_sp", "mt_mp"),
                             out_type = "df")
  dir <- tempfile("asr"); dir.create(dir)
  res <- asr_scenarios(sc, s$ginv, fixed = ~ hhys, inbreeding = s$F, dir = dir, verbose = FALSE)
  expect_named(res, c("st_sp", "st_mp", "mt_sp", "mt_mp"))
  expect_equal(length(list.files(dir, "\\.xlsx$")), 4L)
  for (r in res) {
    expect_true(all(c("varcomp", "vc_index", "param", "accuracy", "summary", "pred") %in% names(r)))
    expect_true(all(c("pa", "accuracy", "slope", "va", "h2_full") %in% names(r$accuracy)))
  }
  expect_true(all(c("rg", "rp", "rg_p", "h2_value1") %in% names(res$mt_sp$param)))
  expect_equal(nrow(res$st_sp$accuracy), nrow(sc$st_sp$cv))
  expect_equal(sort(unique(res$mt_sp$accuracy$trait)), c("adg", "bf"))
})

test_that("cv_vc modes and extra random effects", {
  skip_if_no_mock()
  s <- sim_data()
  sc <- mintyr::cv_scenarios(s$pheno, cols = c("adg", "bf"), id = "id", group = "test_herd",
                             v = 3, seed = 1, scenarios = "mt_mp", out_type = "df")
  for (m in c("fixed", "warm", "estimate")) {
    r <- asr_gblup(sc$mt_mp$param, sc$mt_mp$cv[1:2, ], s$ginv, fixed = ~ hhys,
                   cv_vc = m, verbose = FALSE)
    expect_equal(attr(r, "cv_vc"), m)
    expect_true(all(!is.na(r$accuracy$va)))
  }
  r <- asr_gblup(sc$mt_mp$param, sc$mt_mp$cv[1:2, ], s$ginv, fixed = ~ hhys,
                 random = ~ litter, random_cov = "idh", verbose = FALSE)
  expect_match(r$param$formulas, "h2_value1 = V1/(V1 + V4 + V7)", fixed = TRUE)
})

test_that("mt_bt, traits = and test_rp", {
  skip_if_no_mock()
  s <- sim_data(groups = c("nls", "bfs"))
  sc <- mintyr::cv_scenarios(s$pheno[, c("id", "test_herd", "hhys", "adg")], cols = "adg",
                             id = "id", group = "test_herd", v = 3, seed = 1,
                             scenarios = "mt_bt", out_type = "df")
  r <- asr_scenarios(sc, s$ginv, scenarios = "mt_bt", dir = NULL, fixed = ~ hhys,
                     genetic_cov = "corgh", verbose = FALSE)$mt_bt
  expect_match(r$param$formulas, "rg = V1", fixed = TRUE)
  expect_false("rp" %in% names(r$param))
  expect_true(all(r$accuracy$trait == r$accuracy$test_herd))

  d <- s$pheno; names(d)[names(d) == "bf"] <- "nl"
  nest <- data.table::data.table(name = "adg-nl", data = list(d))
  cv <- mintyr::nest_cv(nest, v = 3, seed = 1, out_type = "df")
  r2 <- asr_gblup(nest, cv, s$ginv, fixed = ~ hhys, traits = c("adg", "nl"),
                  test_rp = TRUE, verbose = FALSE)
  expect_true(all(c("h2_adg", "h2_nl", "rg", "rp", "cov_p") %in% names(r2$param)))
})

test_that("warnings for a matrix without INVERSE and for boundary components", {
  skip_if_no_mock()
  s <- sim_data()
  sc <- mintyr::cv_scenarios(s$pheno, cols = c("adg", "bf"), id = "id", group = "test_herd",
                             v = 3, seed = 1, scenarios = "mt_mp", out_type = "df")
  g <- s$ginv; attr(g, "INVERSE") <- NULL
  expect_warning(asr_gblup(sc$mt_mp$param, sc$mt_mp$cv[1, ], g, verbose = FALSE), "INVERSE")
  old <- options(mock.boundary = TRUE); on.exit(options(old))
  r <- asr_gblup(sc$mt_mp$param, sc$mt_mp$cv[1, ], s$ginv, random = ~ litter,
                 random_cov = "idh", verbose = FALSE)
  expect_match(r$param$boundary, "litter")
  expect_true("boundary" %in% r$log$type)
})

test_that("LRT via lrt() gives rg_p; fixed-variance CV fits are quiet", {
  skip_if_no_mock()
  s <- sim_data()
  sc <- mintyr::cv_scenarios(s$pheno, cols = c("adg", "bf"), id = "id", group = "test_herd",
                             v = 3, seed = 1, scenarios = "mt_mp", out_type = "df")
  expect_no_warning(r <- asr_gblup(sc$mt_mp$param, sc$mt_mp$cv[1:2, ], s$ginv,
                                   fixed = ~ hhys, verbose = FALSE))
  expect_true(is.finite(r$param$rg_p))
  expect_true(all(c("cv_converge", "cv_iter") %in% names(r$accuracy)))
  expect_true(all(r$accuracy$cv_iter >= 1L))
  r2 <- asr_gblup(sc$mt_mp$param, sc$mt_mp$cv[1:2, ], s$ginv, fixed = ~ hhys,
                  cv_vc = "estimate", verbose = FALSE)
  expect_true(all(r2$accuracy$cv_converge))
})

test_that("an LRT with df = 0 gives NA, not p = 0; min_n drops tiny folds", {
  skip_if_no_mock()
  s <- sim_data()
  sc <- mintyr::cv_scenarios(s$pheno, cols = c("adg", "bf"), id = "id", group = "test_herd",
                             v = 3, seed = 1, scenarios = "mt_mp", out_type = "df")
  old <- options(mock.lrt_df0 = TRUE); on.exit(options(old))
  r <- asr_gblup(sc$mt_mp$param, sc$mt_mp$cv[1, ], s$ginv, verbose = FALSE)
  expect_true(is.na(r$param$rg_p))
  expect_true(any(grepl("df = 0", r$log$message)))
  options(old)
  r2 <- asr_gblup(sc$mt_mp$param, sc$mt_mp$cv[1, ], s$ginv, min_n = 1e6, verbose = FALSE)
  expect_true(all(is.na(r2$accuracy$pa)) && all(r2$accuracy$n > 0))
})

mt_setup <- function() {
  s <- sim_data()
  sc <- mintyr::cv_scenarios(s$pheno, cols = c("adg", "bf"), id = "id", group = "test_herd",
                             v = 3, seed = 1, scenarios = "mt_mp", out_type = "df")
  list(s = s, sc = sc)
}

test_that("no reliability filter in CV by default; min_rel flags the final GEBVs", {
  skip_if_no_mock()
  x <- mt_setup()
  r <- asr_gblup(x$sc$mt_mp$param, x$sc$mt_mp$cv[1:2, ], x$s$ginv, inbreeding = x$s$F,
                 verbose = FALSE)
  expect_false(anyNA(r$pred$gebv))                       # nobody dropped by default
  r2 <- asr_gblup(x$sc$mt_mp$param, x$sc$mt_mp$cv[1:2, ], x$s$ginv, inbreeding = x$s$F,
                  cv_min_rel = 0.5, verbose = FALSE)
  expect_true(any(r2$pred$rel < 0.5))
  expect_true(all(is.na(r2$pred$gebv[r2$pred$rel < 0.5])))
  expect_equal(nrow(r$gebv), nrow(x$s$ginv) * 2L)          # all individuals x 2 traits
  expect_equal(r$gebv$reliable, r$gebv$rel >= 0.05)
  expect_setequal(unique(r$gebv$trait), c("adg", "bf"))
})

test_that("update() rescues a non-converged fit", {
  skip_if_no_mock()
  x <- mt_setup()
  old <- options(mock.nonconv = TRUE, mock.update_fixes = TRUE); on.exit(options(old))
  r <- asr_gblup(x$sc$mt_mp$param, x$sc$mt_mp$cv[1, ], x$s$ginv, verbose = FALSE)
  expect_true(r$param$converge)
  expect_equal(r$param$updates, 1L)
  expect_equal(unique(r$accuracy$status), "ok")
})

test_that("non-converged fits are skipped and logged, or kept", {
  skip_if_no_mock()
  x <- mt_setup()
  old <- options(mock.nonconv = TRUE); on.exit(options(old))
  expect_warning(r <- asr_gblup(x$sc$mt_mp$param, x$sc$mt_mp$cv[1:2, ], x$s$ginv,
                                verbose = FALSE), "see `\\$log`")
  expect_equal(r$param$status, "not_converged")
  expect_equal(r$param$updates, 2L)
  expect_true(all(r$accuracy$status == "param_not_converged"))
  expect_true(all(is.na(r$accuracy$pa)))
  expect_true("not_converged" %in% r$log$type)
  expect_true(nrow(r$summary) == 0L || all(r$summary$n == 0L))
  expect_warning(k <- asr_gblup(x$sc$mt_mp$param, x$sc$mt_mp$cv[1:2, ], x$s$ginv,
                                nonconv = "keep", verbose = FALSE))
  expect_equal(k$param$status, "not_converged_kept")
  expect_true(all(k$accuracy$status == "ok"))

  old2 <- options(mock.nonconv = TRUE)
  cvx <- suppressWarnings(asr_gblup(x$sc$mt_mp$param, x$sc$mt_mp$cv[1, ], x$s$ginv,
                                    nonconv = "keep", cv_vc = "estimate", verbose = FALSE))
  expect_true(all(cvx$accuracy$cv_updates == 2L))
})

test_that("a failing fit is logged and the loop goes on", {
  skip_if_no_mock()
  x <- mt_setup()
  cvt <- x$sc$mt_mp$cv[1:2, ]
  cvt$train[[2]]$id <- NULL                     # makes this fold fail
  expect_warning(r <- asr_gblup(x$sc$mt_mp$param, cvt, x$s$ginv, verbose = FALSE), "log")
  expect_equal(sort(unique(r$accuracy$status)), c("failed", "ok"))
  expect_true("error" %in% r$log$type)
  expect_error(asr_gblup(x$sc$mt_mp$param, cvt, x$s$ginv, on_error = "stop", verbose = FALSE))
})

test_that("progress lines", {
  skip_if_no_mock()
  x <- mt_setup()
  expect_message(asr_gblup(x$sc$mt_mp$param, x$sc$mt_mp$cv[1, ], x$s$ginv),
                 "\\[cv 1/1 +100%\\].*left")
})

test_that("parameters only (cv = NULL), F from ainverse(), sparse A^-1 -> no dense", {
  skip_if_no_mock()
  s <- sim_data()
  nest <- mintyr::w2l_nest(s$pheno, cols = c("adg", "bf"), by = "test_herd", out_type = "df")
  ids <- rownames(s$ginv)
  ainv <- cbind(seq_along(ids), seq_along(ids), 1)          # sparse triplet like ainverse()
  attr(ainv, "rowNames") <- ids; attr(ainv, "INVERSE") <- TRUE
  attr(ainv, "inbreeding") <- rep(0.03, length(ids)); attr(ainv, "tbv") <- attr(s$ginv, "tbv")
  expect_message(r <- asr_gblup(nest, ginv = ainv, fixed = ~ hhys), "dense = FALSE")
  expect_equal(nrow(r$param), 4L)                            # 2 traits x 2 herds
  expect_null(r$accuracy); expect_null(r$summary)
  expect_true(all(r$gebv$inb == 0.03))                       # F from attr "inbreeding"
  expect_s3_class(r, "asr_gblup")
})

test_that("asr_compare stacks methods", {
  skip_if_no_mock()
  skip_if_not_installed("writexl")
  s <- sim_data()
  sc <- mintyr::cv_scenarios(s$pheno[, c("id", "test_herd", "hhys", "adg")], cols = "adg",
                             id = "id", group = "test_herd", v = 3, seed = 1,
                             scenarios = "st_sp", out_type = "df")
  dir <- tempfile("cmp")
  cmp <- asr_compare(sc, list(GBLUP = s$ginv, other = s$ginv), scenarios = "st_sp",
                     dir = dir, fixed = ~ hhys, verbose = FALSE)
  expect_setequal(unique(cmp$summary$method), c("GBLUP", "other"))
  expect_true(file.exists(file.path(dir, "compare.xlsx")))
  expect_true(file.exists(file.path(dir, "st_sp_GBLUP.xlsx")))
})

test_that("random terms such as maternal effects run end to end", {
  skip_if_no_mock()
  s <- sim_data()
  d <- s$pheno; d$dam <- sample(d$id, nrow(d), TRUE)
  nest <- data.table::data.table(name = "adg", data = list(d[, c("id", "dam", "hhys", "adg")]))
  r <- asr_gblup(nest, ginv = s$ginv, traits = "adg", fixed = ~ hhys,
                 random = ~ vm(dam, ginv) + ide(dam), verbose = FALSE)
  expect_true(all(c("h2", "prop_vm_dam_ginv", "prop_ide_dam") %in% names(r$param)))
})
