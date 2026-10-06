## Engine fix #105 (CVXPY gh-3442): gathers whose index arrays contain
## duplicates. With engine 0.3.0 (sparsediff 0.4.0) building the derivative
## structure overran a buffer and killed the R process, so this file is kept
## separate: on a broken engine it fails by crashing rather than by an
## expectation.

test_that("duplicate-index gathers give correct derivatives", {
  n <- 4L
  idx <- c(0L, 0L, 1L, 2L, 2L, 3L)  # 0-based, with repeats
  x <- sd_variable(n, 1L, 0L, n)
  g <- sd_index(x, length(idx), 1L, idx)
  y <- sd_elementwise_mult(g, g)
  prob <- sd_problem(sd_sum(y, -1L), list(y), FALSE)
  sd_init_jacobian_coo(prob)
  sd_init_hessian_coo(prob)
  sd_init_derivatives(prob)

  u <- c(1, 2, 3, 4)
  w <- rep(0.5, length(idx))
  expect_equal(sd_objective_forward(prob, u), sum(u[idx + 1L]^2))
  expect_equal(sd_constraint_forward(prob, u), u[idx + 1L]^2)
  expect_equal(sd_gradient(prob), as.vector(2 * tabulate(idx + 1L, n) * u))

  js <- sd_jacobian_sparsity(prob)
  J <- matrix(0, js$nrow, js$ncol)
  J[cbind(js$rows + 1L, js$cols + 1L)] <- sd_jacobian_values(prob)
  Jexp <- matrix(0, length(idx), n)
  Jexp[cbind(seq_along(idx), idx + 1L)] <- 2 * u[idx + 1L]
  expect_equal(J, Jexp)

  # Lagrangian Hessian: 1 * objective + sum_i w_i * y_i, both diagonal.
  hs <- sd_hessian_sparsity(prob)
  H <- lower_to_full(hs, sd_hessian_values(prob, 1, w), n)
  counts <- tabulate(idx + 1L, n)
  expect_equal(H, diag(2 * counts + 2 * 0.5 * counts))
})
