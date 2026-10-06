test_that("the bundled engine is 0.6.1", {
  expect_identical(engine_version(), "0.6.1")
})

test_that("the CBLAS shim computes a row-major dgemm", {
  # A (2x3) %*% B (3x2), both row-major, returned row-major.
  A <- matrix(1:6, 2, 3, byrow = TRUE)
  B <- matrix(7:12, 3, 2, byrow = TRUE)
  expect_equal(sd_selftest_dgemm(), as.vector(t(A %*% B)))
})
