test_that("user-facing functions are exported", {
  ex <- getNamespaceExports("mintyrasr")
  expect_true(all(c("asr_gblup", "asr_scenarios", "asr_compare") %in% ex))
})
