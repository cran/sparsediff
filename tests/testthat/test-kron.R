## Kronecker products (sparsediffpy make_left_kron / make_right_kron, used by
## CVXPY's kron converter). Values are checked against base::kronecker(); the
## atom is linear in the variable, so its Jacobian is checked by finite
## differences and its Hessian is zero.

kron_problem <- function(z) {
  prob <- sd_problem(sd_sum(z, -1L), list(z), FALSE)
  sd_init_jacobian_coo(prob)
  sd_init_hessian_coo(prob)
  sd_init_derivatives(prob)
  prob
}

dense_jacobian <- function(prob) {
  js <- sd_jacobian_sparsity(prob)
  J <- matrix(0, js$nrow, js$ncol)
  J[cbind(js$rows + 1L, js$cols + 1L)] <- sd_jacobian_values(prob)
  J
}

test_that("sd_left_kron matches kronecker(A, X) with structural zeros in A", {
  A <- matrix(c(1, 0, 2, 0, 3, -1), 2, 3)  # p x q = 2 x 3, two zeros
  r <- 2L; s <- 2L
  n <- r * s
  u <- c(0.5, -1, 2, 3)
  x <- sd_variable(r, s, 0L, n)
  Ac <- sd_parameter(2L, 3L, -1L, n, as.vector(A))
  active <- as.integer(which(A != 0) - 1L)  # 0-based column-major
  z <- sd_left_kron(Ac, x, 2L, 3L, r, s, active)
  expect_identical(sd_get_expr_dimensions(z), c(4L, 6L))

  prob <- kron_problem(z)
  K <- kronecker(A, matrix(u, r, s))
  expect_equal(sd_constraint_forward(prob, u), as.vector(K))
  expect_equal(sd_objective_forward(prob, u), sum(K))

  f <- function(v) as.vector(kronecker(A, matrix(v, r, s)))
  Jfd <- vapply(seq_len(n), function(j) {
    e <- numeric(n); e[j] <- 1e-6
    (f(u + e) - f(u - e)) / 2e-6
  }, numeric(length(K)))
  expect_equal(dense_jacobian(prob), Jfd, tolerance = 1e-8)
  hv <- sd_hessian_values(prob, 1, rep(1, length(K)))
  expect_true(all(hv == 0))
})

test_that("sd_right_kron matches kronecker(X, B)", {
  B <- matrix(c(2, -1, 0, 4), 2, 2)  # r x s = 2 x 2, one zero
  p <- 3L; q <- 1L
  n <- p * q
  u <- c(1, -2, 0.5)
  x <- sd_variable(p, q, 0L, n)
  Bc <- sd_parameter(2L, 2L, -1L, n, as.vector(B))
  active <- as.integer(which(B != 0) - 1L)
  prob <- kron_problem(sd_right_kron(Bc, x, p, q, 2L, 2L, active))
  expect_equal(sd_constraint_forward(prob, u), as.vector(kronecker(matrix(u, p, q), B)))
})

test_that("a parametric kron operand follows sd_update_params", {
  A1 <- matrix(c(1, 2, 3, 4), 2, 2)
  A2 <- matrix(c(-1, 0.5, 2, 0), 2, 2)
  n <- 2L
  u <- c(3, -1)
  x <- sd_variable(2L, 1L, 0L, n)
  P <- sd_parameter(2L, 2L, 0L, n, as.vector(A1))
  prob <- kron_problem(sd_left_kron(P, x, 2L, 2L, 2L, 1L, 0:3))  # all blocks
  sd_register_params(prob, list(P))
  expect_equal(sd_constraint_forward(prob, u), as.vector(kronecker(A1, matrix(u))))
  sd_update_params(prob, as.vector(A2))
  expect_equal(sd_constraint_forward(prob, u), as.vector(kronecker(A2, matrix(u))))
})

test_that("kron bindings validate their inputs", {
  n <- 4L
  x <- sd_variable(2L, 2L, 0L, n)
  A <- sd_parameter(2L, 3L, -1L, n, rep(1, 6))
  expect_error(sd_left_kron(A, x, 3L, 2L, 2L, 2L, 0:5), NA)  # 3 x 2 also has 6 entries
  expect_error(sd_left_kron(A, x, 2L, 2L, 2L, 2L, 0:3), "`param` has 6 entries")
  expect_error(sd_left_kron(A, x, 2L, 3L, 1L, 2L, 0:5), "`child` has 4 entries")
  expect_error(sd_left_kron(A, x, 2L, 3L, 2L, 2L, c(0, 1)), "integer vector")
  expect_error(sd_left_kron(A, x, 2L, 3L, 2L, 2L, c(0L, 6L)), "outside")
  expect_error(sd_left_kron(A, x, 2L, 3L, 2L, 2L, -1L), "outside")
  expect_error(sd_right_kron(A, x, 2L, 2L, 2L, 3L, 0:5), NA)
  expect_error(sd_right_kron(A, x, 2L, 2L, 2L, 3L, 6L), "outside")
  expect_error(sd_left_kron(A, x, 0L, 3L, 2L, 2L, 0L), "positive")
})
