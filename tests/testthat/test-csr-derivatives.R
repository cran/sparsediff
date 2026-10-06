## CSR derivative accessors (sparsediffpy problem_jacobian / get_jacobian /
## problem_init_hessian / problem_hessian / get_hessian), checked against the
## COO accessors on the same problem.

csr_to_dense <- function(csr) {
  M <- matrix(0, csr$shape[1L], csr$shape[2L])
  rows <- rep(seq_len(csr$shape[1L]), diff(csr$indptr))
  M[cbind(rows, csr$indices + 1L)] <- csr$data
  M
}

build_problem <- function() {
  n <- 3L
  x <- sd_variable(n, 1L, 0L, n)
  obj <- sd_sum(sd_exp(x), -1L)
  g1 <- sd_elementwise_mult(x, x)               # x_i^2
  g2 <- sd_sum(sd_sin(x), -1L)                  # sum sin(x_i)
  sd_problem(obj, list(g1, g2), FALSE)
}

test_that("CSR Jacobian and Hessian equal the COO forms", {
  u <- c(0.3, -1.2, 2)
  w <- c(0.5, -1, 2, 0.7)
  prob <- build_problem()
  sd_init_jacobian_coo(prob)
  sd_init_hessian_coo(prob)
  sd_init_jacobian(prob)
  sd_init_hessian(prob)
  sd_init_derivatives(prob)
  sd_objective_forward(prob, u)
  sd_constraint_forward(prob, u)

  J <- sd_jacobian(prob)
  expect_named(J, c("data", "indices", "indptr", "shape"))
  expect_identical(J$shape, c(4L, 3L))
  js <- sd_jacobian_sparsity(prob)
  Jcoo <- matrix(0, js$nrow, js$ncol)
  Jcoo[cbind(js$rows + 1L, js$cols + 1L)] <- sd_jacobian_values(prob)
  expect_equal(csr_to_dense(J), Jcoo)
  expect_equal(csr_to_dense(J), rbind(diag(2 * u), cos(u)))
  expect_identical(sd_get_jacobian(prob), J)

  H <- sd_hessian(prob, 1, w)
  expect_identical(H$shape, c(3L, 3L))
  hs <- sd_hessian_sparsity(prob)
  Hcoo <- lower_to_full(hs, sd_hessian_values(prob, 1, w), 3L)
  expect_equal(csr_to_dense(H), Hcoo)
  expect_equal(csr_to_dense(H), diag(exp(u) + 2 * w[1:3] - w[4] * sin(u)))
  expect_identical(sd_get_hessian(prob), H)
})

test_that("CSR accessors refuse to run before initialization", {
  prob <- build_problem()
  expect_error(sd_jacobian(prob), "sd_init_jacobian")
  expect_error(sd_get_jacobian(prob), "not initialized")
  expect_error(sd_hessian(prob, 1, rep(0, 4)), "sd_init_hessian")
  expect_error(sd_get_hessian(prob), "not initialized")
  sd_init_hessian(prob)
  expect_error(sd_hessian(prob, 1, rep(0, 3)), "length 3")
})

## The engine's Hessian init needs the Jacobian structure; called first on a
## fresh problem it used to segfault (sd_init_hessian_coo since 0.4.0).
test_that("the Hessian can be initialized before the Jacobian", {
  u <- c(0.3, -1.2, 2)
  w <- c(0.5, -1, 2, 0.7)
  expected <- diag(exp(u) + 2 * w[1:3] - w[4] * sin(u))

  prob <- build_problem()
  sd_init_hessian_coo(prob)
  sd_objective_forward(prob, u)
  sd_constraint_forward(prob, u)
  hs <- sd_hessian_sparsity(prob)
  expect_equal(lower_to_full(hs, sd_hessian_values(prob, 1, w), 3L), expected)

  prob2 <- build_problem()
  sd_init_hessian(prob2)
  sd_objective_forward(prob2, u)
  sd_constraint_forward(prob2, u)
  expect_equal(csr_to_dense(sd_hessian(prob2, 1, w)), expected)
})

test_that("sd_jacobian of an unconstrained problem is empty with the right shape", {
  x <- sd_variable(2L, 1L, 0L, 2L)
  prob <- sd_problem(sd_sum(sd_exp(x), -1L), list(), FALSE)
  J <- sd_jacobian(prob)
  expect_identical(J$shape, c(0L, 2L))
  expect_identical(J$indptr, 0L)
  expect_length(J$data, 0L)
})
