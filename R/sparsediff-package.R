#' sparsediff: R interface to the SparseDiffEngine differentiation backend
#'
#' R bindings to the 'SparseDiffEngine' C library --- the sparse Jacobian and
#' Hessian differentiation backend used by 'CVXPY' for its Disciplined
#' Nonlinear Programming (DNLP) extension. This package is the R analog of the
#' 'sparsediffpy' Python package and wraps the same C library (pinned at the
#' upstream v0.3.0 release).
#'
#' @useDynLib sparsediff, .registration = TRUE
#' @keywords internal
"_PACKAGE"
