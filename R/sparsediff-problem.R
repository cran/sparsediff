# Documentation for the high-level problem-assembly and derivative-oracle API.
# The functions themselves are the cpp11-generated stubs in R/cpp11.R; these are
# doc-only blocks attached by @name / @rdname.

#' Assemble a differentiable problem
#'
#' Combine an objective expression and a list of constraint expressions
#' (built with the \code{sd_*} atom constructors) into a single problem object
#' whose value and sparse derivatives can be evaluated repeatedly.
#'
#' @param objective an expression handle for the scalar objective.
#' @param constraints a list of expression handles, stacked vertically to form
#'   the constraint vector \eqn{g(x)} (may be empty).
#' @param verbose logical; if \code{TRUE}, print diagnostic information while the
#'   problem is assembled.
#' @param prob a problem handle returned by \code{sd_problem()}.
#' @param params a list of parameter expression handles (see \code{\link{sd_parameter}})
#'   to register for fast re-evaluation under changing data.
#' @param theta a numeric vector of new parameter values, concatenated in the
#'   order the parameters were registered.
#'
#' @return \code{sd_problem()} returns an external-pointer problem handle. The
#'   handle owns the underlying expression graph and frees it when garbage
#'   collected. \code{sd_register_params()} and \code{sd_update_params()} are
#'   called for their side effect and return \code{NULL} invisibly.
#'
#' @details
#' The typical lifecycle is: build expressions with the \code{sd_*} constructors,
#' assemble them with \code{sd_problem()}, initialise the derivative structures
#' once (\code{\link{sd_init_derivatives}} and friends), then evaluate the
#' objective/constraints and their sparse derivatives at as many primal points
#' as needed. When the problem has parameters, register them once with
#' \code{sd_register_params()} and push new values with \code{sd_update_params()}
#' to re-evaluate without rebuilding the graph (the DPP-style fast path).
#'
#' @seealso \code{\link{sd_init_derivatives}} for the evaluation oracle,
#'   \code{\link{sd_variable}} and the atom constructors for building expressions.
#' @name sparsediff-problem
#' @aliases sd_problem sd_register_params sd_update_params
NULL

#' Sparse derivative oracle
#'
#' Initialise and evaluate the value, gradient, sparse constraint Jacobian and
#' sparse lower-triangular Lagrangian Hessian of a problem built with
#' \code{\link{sd_problem}}.
#'
#' @param prob a problem handle from \code{\link{sd_problem}}.
#' @param u a numeric vector: the primal point at which to evaluate, laid out in
#'   the engine's column-major flat ordering (the same ordering used by the
#'   \code{var_id} offsets passed to \code{\link{sd_variable}}).
#' @param obj_w a scalar weight \eqn{\sigma} multiplying the objective Hessian.
#' @param w a numeric vector of constraint multipliers (length equal to the total
#'   constraint size) weighting the constraint Hessians.
#'
#' @return
#' \describe{
#'   \item{\code{sd_init_derivatives}, \code{sd_init_jacobian},
#'     \code{sd_init_jacobian_coo}, \code{sd_init_hessian_coo}}{called for their
#'     side effect; return \code{NULL} invisibly.}
#'   \item{\code{sd_objective_forward}}{the scalar objective value at \code{u}.}
#'   \item{\code{sd_constraint_forward}}{the constraint vector at \code{u}.}
#'   \item{\code{sd_gradient}}{the objective gradient, length \code{n_vars}.}
#'   \item{\code{sd_jacobian_sparsity}, \code{sd_hessian_sparsity}}{a list with
#'     integer \code{rows}/\code{cols} (0-based COO indices) and \code{nrow}/\code{ncol}.}
#'   \item{\code{sd_jacobian_values}}{the constraint-Jacobian nonzeros, matching
#'     the Jacobian sparsity order.}
#'   \item{\code{sd_hessian_values}}{the nonzeros of
#'     \eqn{\sigma\nabla^2 f + \sum_i w_i \nabla^2 g_i}, lower triangle, matching
#'     the Hessian sparsity order.}
#' }
#'
#' @details
#' Evaluation is ordered: a forward pass first
#' (\code{sd_objective_forward} / \code{sd_constraint_forward}) populates the
#' node values at \code{u}, after which \code{sd_gradient},
#' \code{sd_jacobian_values} and \code{sd_hessian_values} read them. Sparsity
#' patterns are structural --- fixed once the corresponding \code{sd_init_*}
#' routine has run --- so they are queried once and reused, while the values are
#' recomputed at each new point. Row and column indices are 0-based (the engine
#' convention; a higher-level modelling layer such as \pkg{CVXR} translates them
#' to 1-based as needed).
#'
#' @seealso \code{\link{sd_problem}}
#' @name sparsediff-oracle
#' @aliases sd_init_derivatives sd_init_jacobian sd_init_jacobian_coo
#'   sd_init_hessian_coo sd_objective_forward sd_constraint_forward sd_gradient
#'   sd_jacobian_sparsity sd_jacobian_values sd_hessian_sparsity sd_hessian_values
NULL
