# Helpers shared by the tests: build a problem over a single n-vector variable,
# evaluate it at u, and return the objective, gradient and the full symmetric
# Lagrangian Hessian as dense matrices (the engine returns the lower triangle in
# COO form, 0-based).

evaluate_problem <- function(obj, n, u, constraints = list(), obj_w = 1, w = numeric(0)) {
  prob <- sd_problem(obj, constraints, FALSE)
  sd_init_jacobian_coo(prob)
  sd_init_hessian_coo(prob)
  sd_init_derivatives(prob)
  val <- sd_objective_forward(prob, u)
  grad <- sd_gradient(prob)
  if (length(constraints) > 0L) sd_constraint_forward(prob, u)
  hs <- sd_hessian_sparsity(prob)
  hv <- sd_hessian_values(prob, obj_w, w)
  list(prob = prob, value = val, gradient = grad, hessian = lower_to_full(hs, hv, n))
}

lower_to_full <- function(sparsity, values, n) {
  H <- matrix(0, n, n)
  H[cbind(sparsity$rows + 1L, sparsity$cols + 1L)] <- values
  H + t(H) - diag(diag(H), n)
}

# CSR arrays of a dense matrix, every entry stored.
dense_csr <- function(Q) {
  m <- nrow(Q); k <- ncol(Q)
  list(p = as.integer(seq(0L, m * k, by = k)),
       i = rep(0:(k - 1L), m),
       x = as.vector(t(Q)))
}

central_gradient <- function(f, u, eps = 1e-6) {
  vapply(seq_along(u), function(j) {
    up <- u; up[j] <- up[j] + eps
    um <- u; um[j] <- um[j] - eps
    (f(up) - f(um)) / (2 * eps)
  }, numeric(1))
}

random_spd <- function(n) {
  M <- matrix(stats::rnorm(n * n), n, n)
  crossprod(M) + diag(n)
}
