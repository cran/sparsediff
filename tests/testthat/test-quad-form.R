## Quadratic forms x' Q x: the sparse binding (sd_quad_form), the dense
## constant binding and the dense parametric binding (sd_quad_form_dense).

test_that("sparse and dense constant quad_form agree on value, gradient and Hessian", {
  set.seed(1)
  n <- 5L
  Q <- random_spd(n)
  u <- stats::rnorm(n)
  csr <- dense_csr(Q)

  sparse <- evaluate_problem(
    sd_quad_form(sd_variable(n, 1L, 0L, n), csr$p, csr$i, csr$x), n, u)
  dense <- evaluate_problem(
    sd_quad_form_dense(NULL, sd_variable(n, 1L, 0L, n), as.vector(t(Q))), n, u)

  expect_equal(sparse$value, drop(t(u) %*% Q %*% u))
  expect_equal(dense$value, sparse$value)
  expect_equal(dense$gradient, drop(2 * Q %*% u))
  expect_equal(sparse$gradient, dense$gradient)
  expect_equal(dense$hessian, 2 * Q)
  expect_equal(sparse$hessian, dense$hessian)
})

test_that("parametric quad_form follows sd_update_params", {
  set.seed(2)
  n <- 4L
  Q1 <- random_spd(n)
  Q2 <- random_spd(n)
  u <- stats::rnorm(n)

  x <- sd_variable(n, 1L, 0L, n)
  P <- sd_parameter(n, n, 0L, n, as.vector(Q1))
  res <- evaluate_problem(sd_quad_form_dense(P, x, numeric(0)), n, u)
  sd_register_params(res$prob, list(P))
  expect_equal(res$value, drop(t(u) %*% Q1 %*% u))

  sd_update_params(res$prob, as.vector(Q2))
  expect_equal(sd_objective_forward(res$prob, u), drop(t(u) %*% Q2 %*% u))
  expect_equal(sd_gradient(res$prob), drop(2 * Q2 %*% u))
  hs <- sd_hessian_sparsity(res$prob)
  expect_equal(lower_to_full(hs, sd_hessian_values(res$prob, 1, numeric(0)), n), 2 * Q2)
})

test_that("parametric quad_form with a composite source g * Q follows updates (engine #107)", {
  set.seed(3)
  n <- 3L
  Q <- random_spd(n)
  u <- stats::rnorm(n)

  x <- sd_variable(n, 1L, 0L, n)
  g <- sd_parameter(1L, 1L, 0L, n, 2)
  Qc <- sd_parameter(n, n, -1L, n, as.vector(Q))  # fixed constant
  qf <- sd_quad_form_dense(sd_vector_mult(sd_promote(g, n, n), Qc), x, numeric(0))
  res <- evaluate_problem(qf, n, u)
  sd_register_params(res$prob, list(g))
  expect_equal(res$value, 2 * drop(t(u) %*% Q %*% u))

  sd_update_params(res$prob, 3)
  expect_equal(sd_objective_forward(res$prob, u), 3 * drop(t(u) %*% Q %*% u))
  expect_equal(sd_gradient(res$prob), drop(6 * Q %*% u))
})

test_that("dense quad_form over a composite child has the right gradient", {
  set.seed(4)
  n <- 6L
  k <- 3L
  idx <- c(4L, 0L, 2L)  # 0-based
  Q <- random_spd(k)
  u <- stats::rnorm(n)

  build <- function() {
    x <- sd_variable(n, 1L, 0L, n)
    sd_quad_form_dense(NULL, sd_exp(sd_index(x, k, 1L, idx)), as.vector(t(Q)))
  }
  res <- evaluate_problem(build(), n, u)
  f <- function(v) { z <- exp(v[idx + 1L]); drop(t(z) %*% Q %*% z) }
  expect_equal(res$value, f(u))
  expect_equal(res$gradient, central_gradient(f, u), tolerance = 1e-6)
})

test_that("sd_quad_form_dense validates its inputs", {
  n <- 3L
  x <- sd_variable(n, 1L, 0L, n)
  Q <- diag(n)
  P <- sd_parameter(n, n, 0L, n, as.vector(Q))

  expect_error(sd_quad_form_dense(NULL, x, numeric(0)), "exactly one")
  expect_error(sd_quad_form_dense(P, x, as.vector(Q)), "exactly one")
  expect_error(sd_quad_form_dense(NULL, x, rep(1, 8)), "length 8")
  expect_error(sd_quad_form_dense(NULL, x, 1:9), "double vector")
  A <- Q; A[1, 2] <- 1
  expect_error(sd_quad_form_dense(NULL, x, as.vector(t(A))), "not symmetric")
  small <- sd_parameter(2L, 2L, 0L, n, rep(1, 4))
  expect_error(sd_quad_form_dense(small, x, numeric(0)), "4 entries")
  M <- sd_variable(2L, 2L, 0L, 4L)
  expect_error(sd_quad_form_dense(NULL, M, as.vector(diag(4))), "must be a vector")
})
