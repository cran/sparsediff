# Documentation for the expression-graph constructors (the sd_* atoms). The
# functions are the cpp11-generated stubs in R/cpp11.R; these doc-only blocks
# group them into families via @name / @rdname / @aliases.
#
# Every constructor returns an *expression handle*: an external pointer into the
# SparseDiffEngine DAG that participates in reference-counted memory management
# and is freed when garbage collected. Handles are combined into a problem with
# sd_problem() and differentiated through the sparsediff-oracle interface.
#
# Indexing conventions (engine, 0-based): var_id / param_id are flat
# column-major offsets into the primal vector; sd_index() indices are 0-based
# column-major flat offsets; reduction axis is -1 (all entries), 0 (rows) or
# 1 (columns).

#' Leaf expressions: variables and parameters
#'
#' Create the leaves of an expression graph. A variable is a slice of the
#' differentiation vector; a parameter is fixed data that can be updated between
#' evaluations (see \code{\link{sd_register_params}}).
#'
#' @param d1,d2 row and column dimensions of the leaf.
#' @param var_id 0-based flat (column-major) offset of this variable's first
#'   entry within the primal vector \code{u}.
#' @param param_id 0-based parameter offset, analogous to \code{var_id}.
#' @param n_vars total number of variables in the problem.
#' @param values numeric data for the parameter (length \code{d1 * d2},
#'   column-major).
#' @return An expression handle.
#' @seealso \code{\link{sparsediff-elementwise}}, \code{\link{sd_problem}}
#' @usage
#' sd_variable(d1, d2, var_id, n_vars)
#' sd_parameter(d1, d2, param_id, n_vars, values)
#' @name sparsediff-leaves
#' @aliases sd_variable sd_parameter
NULL

#' Expression shape
#'
#' The dimensions and the number of entries of an expression node (the
#' 'sparsediffpy' \code{get_expr_dimensions} and \code{get_expr_size}).
#'
#' @param node an expression handle.
#' @return \code{sd_get_expr_dimensions}: an integer vector \code{c(d1, d2)},
#'   the node's rows and columns. \code{sd_get_expr_size}: an integer scalar,
#'   the number of entries \code{d1 * d2}.
#' @seealso \code{\link{sparsediff-leaves}}
#' @usage
#' sd_get_expr_dimensions(node)
#' sd_get_expr_size(node)
#' @name sparsediff-getters
#' @aliases sd_get_expr_dimensions sd_get_expr_size
NULL

#' Elementwise atoms
#'
#' Smooth elementwise functions of a single expression. Each returns an
#' expression of the same shape as its argument.
#'
#' @param child,c an expression handle (the argument).
#' @param p exponent for \code{sd_power}.
#' @return An expression handle.
#' @details
#' \describe{
#'   \item{\code{sd_exp}, \code{sd_log}}{exponential and natural logarithm.}
#'   \item{\code{sd_sin}, \code{sd_cos}, \code{sd_tan}}{trigonometric functions.}
#'   \item{\code{sd_sinh}, \code{sd_tanh}, \code{sd_asinh}, \code{sd_atanh}}{hyperbolic and inverse-hyperbolic functions.}
#'   \item{\code{sd_logistic}}{the logistic \eqn{\log(1 + e^x)}.}
#'   \item{\code{sd_xexp}}{\eqn{x e^x}.}
#'   \item{\code{sd_normal_cdf}}{the standard normal CDF.}
#'   \item{\code{sd_entr}}{the elementwise entropy \eqn{-x \log x}.}
#'   \item{\code{sd_power}}{the power \eqn{x^p}.}
#'   \item{\code{sd_neg}}{negation \eqn{-x}.}
#' }
#' @seealso \code{\link{sparsediff-affine}}, \code{\link{sparsediff-bivariate}}
#' @usage
#' sd_exp(child)
#' sd_log(c)
#' sd_sin(c)
#' sd_cos(c)
#' sd_tan(c)
#' sd_sinh(c)
#' sd_tanh(c)
#' sd_asinh(c)
#' sd_atanh(c)
#' sd_logistic(c)
#' sd_xexp(c)
#' sd_normal_cdf(c)
#' sd_entr(c)
#' sd_power(c, p)
#' sd_neg(child)
#' @name sparsediff-elementwise
#' @aliases sd_exp sd_log sd_sin sd_cos sd_tan sd_sinh sd_tanh sd_asinh sd_atanh sd_logistic sd_xexp sd_normal_cdf sd_entr sd_power sd_neg
NULL

#' Affine and shape atoms
#'
#' Affine combinations and shape manipulations of expressions. These have
#' constant (zero) second derivatives but participate in the Jacobian.
#'
#' @param left,right,child,c expression handles.
#' @param d1,d2 target row and column dimensions (for \code{sd_promote},
#'   \code{sd_reshape}, \code{sd_broadcast}, \code{sd_index}).
#' @param indices 0-based column-major flat indices selected by \code{sd_index}.
#' @param args a list of expression handles to stack
#'   (\code{sd_hstack}, \code{sd_vstack}).
#' @param n_vars total number of variables in the problem.
#' @param axis reduction axis for \code{sd_sum}: \code{-1} (all entries),
#'   \code{0} (down rows) or \code{1} (across columns).
#' @return An expression handle.
#' @details
#' \describe{
#'   \item{\code{sd_add}}{elementwise sum of two expressions.}
#'   \item{\code{sd_sum}}{sum reduction along \code{axis}.}
#'   \item{\code{sd_trace}}{matrix trace.}
#'   \item{\code{sd_transpose}}{matrix transpose.}
#'   \item{\code{sd_diag_vec}}{diagonal matrix from a vector.}
#'   \item{\code{sd_diag_mat}}{diagonal vector from a matrix.}
#'   \item{\code{sd_upper_tri}}{the strict upper-triangular entries.}
#'   \item{\code{sd_promote}}{promote a scalar to shape \code{d1 x d2}.}
#'   \item{\code{sd_reshape}}{reshape to \code{d1 x d2} (column-major).}
#'   \item{\code{sd_broadcast}}{broadcast to \code{d1 x d2}.}
#'   \item{\code{sd_index}}{select entries by 0-based flat \code{indices}.}
#'   \item{\code{sd_hstack}, \code{sd_vstack}}{horizontal / vertical stacking.}
#' }
#' @seealso \code{\link{sparsediff-elementwise}}, \code{\link{sparsediff-matrix}}
#' @usage
#' sd_add(left, right)
#' sd_sum(child, axis)
#' sd_trace(c)
#' sd_transpose(c)
#' sd_diag_vec(c)
#' sd_diag_mat(c)
#' sd_upper_tri(c)
#' sd_promote(c, d1, d2)
#' sd_reshape(c, d1, d2)
#' sd_broadcast(c, d1, d2)
#' sd_index(child, d1, d2, indices)
#' sd_hstack(args, n_vars)
#' sd_vstack(args, n_vars)
#' @name sparsediff-affine
#' @aliases sd_add sd_sum sd_trace sd_transpose sd_diag_vec sd_diag_mat sd_upper_tri sd_promote sd_reshape sd_broadcast sd_index sd_hstack sd_vstack
NULL

#' Bivariate atoms
#'
#' Functions of two expression arguments.
#'
#' @param l,r,x,y expression handles.
#' @return An expression handle.
#' @details
#' \describe{
#'   \item{\code{sd_elementwise_mult}}{elementwise (Hadamard) product.}
#'   \item{\code{sd_matmul}}{matrix product \eqn{x y}.}
#'   \item{\code{sd_quad_over_lin}}{the quadratic-over-linear \eqn{\lVert x \rVert^2 / y}.}
#'   \item{\code{sd_rel_entr}}{relative entropy \eqn{x \log(x / y)}. Like
#'     'sparsediffpy' \code{make_rel_entr}, it dispatches on operand size: a
#'     scalar \code{l} with a non-scalar \code{r} uses
#'     \code{sd_rel_entr_first_scalar}, the reverse uses
#'     \code{sd_rel_entr_second_scalar}, and otherwise it is elementwise.}
#'   \item{\code{sd_rel_entr_first_scalar}, \code{sd_rel_entr_second_scalar}}{relative
#'     entropy with a scalar first or second argument broadcast against the other.}
#' }
#' @seealso \code{\link{sparsediff-elementwise}}, \code{\link{sparsediff-reduction}}
#' @usage
#' sd_elementwise_mult(l, r)
#' sd_matmul(x, y)
#' sd_quad_over_lin(l, r)
#' sd_rel_entr(l, r)
#' sd_rel_entr_first_scalar(l, r)
#' sd_rel_entr_second_scalar(l, r)
#' @name sparsediff-bivariate
#' @aliases sd_elementwise_mult sd_matmul sd_quad_over_lin sd_rel_entr sd_rel_entr_first_scalar sd_rel_entr_second_scalar
NULL

#' Product-reduction atoms
#'
#' Multiplicative reductions of an expression.
#'
#' @param c an expression handle.
#' @return An expression handle.
#' @details
#' \describe{
#'   \item{\code{sd_prod}}{product of all entries.}
#'   \item{\code{sd_prod_axis_zero}}{column-wise products (reduce down rows).}
#'   \item{\code{sd_prod_axis_one}}{row-wise products (reduce across columns).}
#' }
#' @seealso \code{\link{sparsediff-affine}}
#' @usage
#' sd_prod(c)
#' sd_prod_axis_zero(c)
#' sd_prod_axis_one(c)
#' @name sparsediff-reduction
#' @aliases sd_prod sd_prod_axis_zero sd_prod_axis_one
NULL

#' Parameter- and constant-matrix atoms
#'
#' Operations that combine an expression with fixed data --- a registered
#' parameter node or a constant matrix --- so the data can flow through the
#' differentiated graph (and, for parameters, be updated between evaluations).
#'
#' @param param a parameter expression handle (see \code{\link{sd_parameter}}).
#' @param child an expression handle (the variable argument).
#' @param Qp,Qi,Qx the row-pointer, column-index and value arrays of a
#'   compressed-sparse-row (CSR) matrix \eqn{Q} for \code{sd_quad_form}'s
#'   \eqn{x^\top Q x}. \eqn{Q} must be symmetric, so its CSR arrays equal its
#'   compressed-sparse-column arrays and the \code{@p}, \code{@i}, \code{@x}
#'   slots of a \code{Matrix::dgCMatrix} holding \eqn{Q} can be passed directly.
#' @param Ap,Ai,Ax the row-pointer, column-index and value arrays of a constant
#'   matrix \eqn{A} in compressed-sparse-row (CSR) form, for the sparse matrix
#'   products. These are the \code{@p}, \code{@i}, \code{@x} slots of a
#'   \code{Matrix::dgCMatrix} holding \eqn{A^\top}, not \eqn{A}.
#' @param ncol number of columns of the sparse constant matrix \eqn{A}.
#' @param m,n row and column dimensions of the dense constant matrix.
#' @param p,q,r,s for the Kronecker products \eqn{Z = A \otimes B}: \eqn{A} is
#'   \eqn{p \times q} and \eqn{B} is \eqn{r \times s}.
#' @param active_blocks for the Kronecker products, an integer vector of 0-based
#'   column-major indices of the nonzero entries of the variable-free operand
#'   (\code{param}); only the output rows they cover are built. For a
#'   parametric operand pass every index, \code{0:(length - 1)}.
#' @param data the dense constant-matrix entries in row-major order (length
#'   \code{m * n} for the matrix products, \code{n * n} for
#'   \code{sd_quad_form_dense}); for an R matrix \code{M}, pass
#'   \code{as.vector(t(M))}. Pass \code{numeric(0)} when the matrix comes from
#'   \code{param}.
#' @return An expression handle.
#' @details
#' \describe{
#'   \item{\code{sd_scalar_mult}, \code{sd_vector_mult}}{multiply a child by a
#'     scalar / vector parameter.}
#'   \item{\code{sd_convolve}}{convolution of a parameter kernel with a child.}
#'   \item{\code{sd_quad_form}}{the quadratic form \eqn{x^\top Q x} with sparse
#'     constant \eqn{Q}.}
#'   \item{\code{sd_quad_form_dense}}{the quadratic form \eqn{x^\top Q x} with a
#'     dense symmetric \eqn{n \times n} \eqn{Q}, where \eqn{n} is the length of
#'     the vector \code{child}. Supply exactly one source: \code{param = NULL}
#'     and \code{data} for a constant \eqn{Q} (checked for symmetry), or a
#'     parameter handle of size \eqn{n^2} with \code{data = numeric(0)} for a
#'     parametric \eqn{Q} that follows \code{\link{sd_update_params}}. A
#'     parametric \eqn{Q} is not checked: keeping it symmetric is the caller's
#'     responsibility.}
#'   \item{\code{sd_left_matmul}, \code{sd_right_matmul}}{left / right product
#'     with a sparse constant matrix \eqn{A}.}
#'   \item{\code{sd_left_matmul_dense}, \code{sd_right_matmul_dense}}{left / right
#'     product with a dense constant or parametric matrix.}
#'   \item{\code{sd_left_kron}, \code{sd_right_kron}}{the Kronecker product
#'     \eqn{A \otimes B} with the variable-free operand \code{param} on the
#'     left (\eqn{A}) or on the right (\eqn{B}) and the variable operand
#'     \code{child} on the other side. \code{param} may be a parameter or a
#'     constant made with \code{sd_parameter(..., param_id = -1, ...)}.}
#' }
#' A parameter cannot be the source of the sparse products
#' \code{sd_left_matmul} / \code{sd_right_matmul} (the engine does not support
#' it; 'sparsediffpy' accepts the argument but the engine then aborts); use the
#' dense products for a parametric matrix.
#' @seealso \code{\link{sd_parameter}}, \code{\link{sd_register_params}}
#' @usage
#' sd_scalar_mult(param, child)
#' sd_vector_mult(param, child)
#' sd_convolve(param, child)
#' sd_quad_form(child, Qp, Qi, Qx)
#' sd_quad_form_dense(param, child, data)
#' sd_left_matmul(child, Ap, Ai, Ax, ncol)
#' sd_right_matmul(child, Ap, Ai, Ax, ncol)
#' sd_left_matmul_dense(param, child, m, n, data)
#' sd_right_matmul_dense(param, child, m, n, data)
#' sd_left_kron(param, child, p, q, r, s, active_blocks)
#' sd_right_kron(param, child, p, q, r, s, active_blocks)
#' @name sparsediff-matrix
#' @aliases sd_scalar_mult sd_vector_mult sd_convolve sd_quad_form sd_quad_form_dense sd_left_matmul sd_right_matmul sd_left_matmul_dense sd_right_matmul_dense sd_left_kron sd_right_kron
NULL
