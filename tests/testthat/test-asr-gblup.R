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
  expect_warning(r <- asr_gblup(sc$mt_mp$param, sc$mt_mp$cv[1, ], s$ginv, random = ~ litter,
                                random_cov = "idh", verbose = FALSE), "boundary")
  expect_match(r$param$boundary, "litter")
})

test_that("rg LRT uses asreml::lrt() (lrt.asreml is not exported in ASReml-R 4.2)", {
  skip_if_no_mock()
  s <- sim_data()
  sc <- mintyr::cv_scenarios(s$pheno, cols = c("adg", "bf"), id = "id", group = "test_herd",
                             v = 3, seed = 1, scenarios = "mt_mp", out_type = "df")
  r <- asr_gblup(sc$mt_mp$param, sc$mt_mp$cv[1, ], s$ginv, fixed = ~ hhys, verbose = FALSE)
  expect_false(anyNA(r$param$rg_lrt))
  expect_false(anyNA(r$param$rg_p))
})

test_that("min_rel defaults to -Inf and keeps all GEBVs", {
  skip_if_no_mock()
  s <- sim_data()
  sc <- mintyr::cv_scenarios(s$pheno, cols = c("adg", "bf"), id = "id", group = "test_herd",
                             v = 3, seed = 1, scenarios = "st_sp", out_type = "df")
  r <- asr_gblup(sc$st_sp$param, sc$st_sp$cv, s$ginv, verbose = FALSE)
  expect_false(anyNA(r$pred$gebv))
  expect_warning(asr_gblup(sc$st_sp$param, sc$st_sp$cv[1, ], s$ginv, min_rel = 2, verbose = FALSE),
                 "min_rel")
})

test_that("non-converged fits: update(), skip, keep and log", {
  skip_if_no_mock()
  s <- sim_data()
  sc <- mintyr::cv_scenarios(s$pheno, cols = c("adg", "bf"), id = "id", group = "test_herd",
                             v = 3, seed = 1, scenarios = "st_mp", out_type = "df")
  P <- sc$st_mp$param; C <- sc$st_mp$cv
  old <- options(mock.nonconv = TRUE, mock.update_fixes = NULL); on.exit(options(old))
  # skip: warning, nothing built on the fit, recorded in log
  w <- testthat::capture_warnings(r <- asr_gblup(P, C, s$ginv, max_update = 2, verbose = FALSE))
  expect_true(any(grepl("did not converge after 2 update", w)))
  expect_true(any(grepl("see the `log` table", w)))
  expect_true(all(r$accuracy$status == "skipped"))
  expect_true(all(is.na(r$accuracy$pa)))
  expect_true(all(r$param$converge == FALSE))
  expect_true(all(r$param$n_update == 2L))
  expect_true(all(c("param", "cv") %in% r$log$stage))
  expect_true(any(r$log$status == "not_converged"))
  # keep: results used anyway
  r2 <- suppressWarnings(asr_gblup(P, C, s$ginv, on_nonconverge = "keep", verbose = FALSE))
  expect_true(all(r2$accuracy$status == "ok"))
  expect_false(all(is.na(r2$accuracy$pa)))
  # update() rescues the fit
  options(mock.update_fixes = TRUE)
  r3 <- asr_gblup(P, C, s$ginv, verbose = FALSE)
  expect_true(all(r3$param$converge))
  expect_true(all(r3$param$n_update == 1L))
  expect_true(all(r3$log$status == "ok"))
})

test_that("failed fits are logged; expected warnings of fixed CV fits are muffled", {
  skip_if_no_mock()
  s <- sim_data()
  sc <- mintyr::cv_scenarios(s$pheno, cols = c("adg", "bf"), id = "id", group = "test_herd",
                             v = 3, seed = 1, scenarios = "st_mp", out_type = "df")
  old <- options(mock.warn = "Log-likelihood not converged", mock.fail = NULL); on.exit(options(old))
  expect_silent(r <- asr_gblup(sc$st_mp$param, sc$st_mp$cv, s$ginv, verbose = FALSE))
  expect_true(all(r$log$message == ""))
  options(mock.warn = NULL, mock.fail = TRUE)
  w <- testthat::capture_warnings(r <- asr_gblup(sc$st_mp$param, sc$st_mp$cv, s$ginv, verbose = FALSE))
  expect_true(any(grepl("ASReml fit failed", w)))
  expect_true(all(r$log$status[r$log$stage == "param"] == "failed"))
  expect_true(all(r$accuracy$status == "skipped"))
  expect_error(asr_gblup(sc$st_mp$param, sc$st_mp$cv, s$ginv, on_error = "stop", verbose = FALSE),
               "ASReml fit failed")
})

test_that("progress and summary messages", {
  skip_if_no_mock()
  s <- sim_data()
  sc <- mintyr::cv_scenarios(s$pheno, cols = c("adg", "bf"), id = "id", group = "test_herd",
                             v = 3, seed = 1, scenarios = "st_mp", out_type = "df")
  msgs <- testthat::capture_messages(asr_gblup(sc$st_mp$param, sc$st_mp$cv, s$ginv, inbreeding = s$F))
  expect_true(any(grepl("^\\[cv 1/\\d+ +0%\\].*elapsed .*ETA", msgs)))
  expect_true(any(grepl("^Done in ", msgs)))
})
