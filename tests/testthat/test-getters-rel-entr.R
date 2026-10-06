test_that("shape getters report dimensions and size", {
  x <- sd_variable(3L, 2L, 0L, 6L)
  expect_identical(sd_get_expr_dimensions(x), c(3L, 2L))
  expect_identical(sd_get_expr_size(x), 6L)
  s <- sd_sum(x, -1L)
  expect_identical(sd_get_expr_size(s), 1L)
})

## sd_rel_entr dispatches on operand size like sparsediffpy make_rel_entr. The
## engine requires both arguments to be distinct variable leaves (CVXPY's
## canonicalization guarantees it), so each case uses two sd_variable()s.
rel_entr_value <- function(nl, nr, u) {
  n <- nl + nr
  l <- sd_variable(nl, 1L, 0L, n)
  r <- sd_variable(nr, 1L, nl, n)
  z <- sd_rel_entr(l, r)
  prob <- sd_problem(sd_sum(z, -1L), list(z), FALSE)
  sd_init_jacobian_coo(prob)
  sd_init_hessian_coo(prob)
  sd_init_derivatives(prob)
  list(dims = sd_get_expr_dimensions(z), value = sd_constraint_forward(prob, u))
}

test_that("sd_rel_entr handles scalar-vector, vector-scalar and elementwise", {
  u <- c(2, 1, 3, 0.5)
  a <- rel_entr_value(1L, 3L, u)   # scalar l, vector r
  expect_identical(a$dims, c(3L, 1L))
  expect_equal(a$value, u[1] * log(u[1] / u[2:4]))
  b <- rel_entr_value(3L, 1L, u)   # vector l, scalar r
  expect_equal(b$value, u[1:3] * log(u[1:3] / u[4]))
  c2 <- rel_entr_value(2L, 2L, u)  # elementwise
  expect_equal(c2$value, u[1:2] * log(u[1:2] / u[3:4]))
})

test_that("relative-entropy bindings require two distinct variables", {
  x <- sd_variable(1L, 1L, 0L, 4L)
  y <- sd_variable(3L, 1L, 1L, 4L)
  expect_error(sd_rel_entr(sd_exp(x), y), "must be variables")
  expect_error(sd_rel_entr_first_scalar(x, sd_exp(y)), "must be variables")
  expect_error(sd_rel_entr_second_scalar(sd_exp(y), x), "must be variables")
  v <- sd_variable(2L, 1L, 0L, 2L)
  expect_error(sd_rel_entr(v, v), "different variables")
})
