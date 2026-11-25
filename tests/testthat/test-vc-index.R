# Parser and formula builder: no ASReml needed

f_txt <- function(f) vapply(f, function(x) deparse(x[[3L]]), "")

test_that("single trait, idv(units) scale row is skipped", {
  vc <- vc_df(c("vm(id, ginv)", "ide(id)", "units!units", "units!R"), bound = c("P", "P", "P", "F"))
  idx <- mintyrasr:::.vc_index(vc, "st", "id", "value")
  expect_equal(idx$role, c("genetic", "random", "residual", "residual"))
  expect_equal(idx$used, c(TRUE, TRUE, TRUE, FALSE))
  expect_equal(f_txt(mintyrasr:::.vp_formulas(idx, "value")), c(h2 = "V1/(V1 + V2 + V3)"))
})

test_that("default residual units!R (not fixed) is the residual variance", {
  vc <- vc_df(c("vm(id, ginv)", "units!R"))
  idx <- mintyrasr:::.vc_index(vc, "st", "id", "value")
  expect_true(all(idx$used))
})

test_that("us G + us R reproduces the classic formulas", {
  nm <- c("trait:vm(id, ginv)!trait_value1:value1", "trait:vm(id, ginv)!trait_value2:value1",
          "trait:vm(id, ginv)!trait_value2:value2", "units:trait!R",
          "units:trait!trait_value1:value1", "units:trait!trait_value2:value1",
          "units:trait!trait_value2:value2")
  idx <- mintyrasr:::.vc_index(vc_df(nm, bound = c("P", "U", "P", "F", "P", "U", "P")),
                               "mt", "id", c("value1", "value2"))
  f <- f_txt(mintyrasr:::.vp_formulas(idx, c("value1", "value2")))
  expect_equal(unname(f["h2_value1"]), "V1/(V1 + V5)")
  expect_equal(unname(f["h2_value2"]), "V3/(V3 + V7)")
  expect_equal(unname(f["rg"]), "V2/sqrt(V1 * V3)")
  expect_equal(unname(f["rp"]), "(V2 + V6)/sqrt((V1 + V5) * (V3 + V7))")
})

test_that("corgh G + diag R with any trait names", {
  nm <- c("trait:vm(id, ginv)!trait!bfs:!trait!nls.cor", "trait:vm(id, ginv)!trait_nls",
          "trait:vm(id, ginv)!trait_bfs", "units:trait!R", "units:trait!trait_nls",
          "units:trait!trait_bfs")
  idx <- mintyrasr:::.vc_index(vc_df(nm), "mt", "id", c("nls", "bfs"))
  expect_equal(idx$kind[1], "cor")
  f <- f_txt(mintyrasr:::.vp_formulas(idx, c("nls", "bfs"), rp = FALSE))
  expect_equal(unname(f), c("V2/(V2 + V5)", "V3/(V3 + V6)", "V1"))
  expect_false("rp" %in% names(f))
})

test_that("extra random term enters the phenotypic variance", {
  nm <- c("trait:vm(id, ginv)!trait_a:a", "trait:vm(id, ginv)!trait_b:a", "trait:vm(id, ginv)!trait_b:b",
          "trait:litter!trait_a", "trait:litter!trait_b", "units:trait!R",
          "units:trait!trait_a:a", "units:trait!trait_b:a", "units:trait!trait_b:b")
  idx <- mintyrasr:::.vc_index(vc_df(nm), "mt", "id", c("a", "b"))
  f <- f_txt(mintyrasr:::.vp_formulas(idx, c("a", "b")))
  expect_equal(unname(f["h2_a"]), "V1/(V1 + V4 + V7)")
  expect_equal(unname(f["rp"]), "(V2 + V8)/sqrt((V1 + V4 + V7) * (V3 + V5 + V9))")
})

test_that("V numbers follow vparameters, not the varcomp row order", {
  nm <- c("vm(id, ginv)", "units!units", "units!R")
  vc <- vc_df(rev(nm), bound = c("F", "P", "P"))
  idx <- mintyrasr:::.vc_index(vc, "st", "id", "value", vnames = nm)
  expect_equal(f_txt(mintyrasr:::.vp_formulas(idx, "value")), c(h2 = "V1/(V1 + V2)"))
})

test_that("a pedigree term vm(dam, ginv) is not taken as the genomic term", {
  vc <- vc_df(c("vm(id, ginv)", "vm(dam, ginv)", "units!R"))
  idx <- mintyrasr:::.vc_index(vc, "st", "id", "value")
  expect_equal(idx$role, c("genetic", "random", "residual"))
})
