## Engine fix #107: a matrix product whose parametric source is a composite of a
## parameter must use the updated parameter after sd_update_params(). The shape
## is the one CVXR emits for (g * A) %*% x with a scalar Parameter g and a
## constant A: vector_mult(promote(g), A) as the source of a dense left matmul.
## With engine 0.3.0 (sparsediff 0.4.0) the composite was never refreshed and
## the old value of g was used.

test_that("dense left matmul with a composite parametric source follows updates", {
  set.seed(5)
  m <- 2L
  n <- 3L
  A <- matrix(stats::rnorm(m * n), m, n)
  u <- stats::rnorm(n)

  x <- sd_variable(n, 1L, 0L, n)
  g <- sd_parameter(1L, 1L, 0L, n, 2)
  Ac <- sd_parameter(m, n, -1L, n, as.vector(A))  # fixed constant
  src <- sd_vector_mult(sd_promote(g, m, n), Ac)
  y <- sd_left_matmul_dense(src, x, m, n, numeric(0))
  prob <- sd_problem(sd_sum(y, -1L), list(y), FALSE)
  sd_register_params(prob, list(g))
  sd_init_jacobian_coo(prob)
  sd_init_derivatives(prob)

  expect_equal(sd_constraint_forward(prob, u), drop(2 * A %*% u))
  sd_update_params(prob, 3)
  expect_equal(sd_constraint_forward(prob, u), drop(3 * A %*% u))
  expect_equal(sd_objective_forward(prob, u), sum(3 * A %*% u))
})
