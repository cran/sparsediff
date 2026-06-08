// sparsediff: cpp11 bindings to the SparseDiffEngine C library.
//
// The engine is pure C, so its headers are included inside extern "C". Engine
// `expr` (expression-DAG node) and `problem` objects are opaque foreign handles,
// so they cross the R boundary as SEXP external pointers (cpp11 cannot see the
// engine types in the generated cpp11.cpp, so we use SEXP rather than
// cpp11::external_pointer for these). cpp11 still handles the scalar/vector
// arguments, return vectors, error unwinding, and routine registration.
//
// Memory model. Engine nodes start at refcount 0 (init_expr); each parent and
// new_problem() expr_retain()s the nodes it references; free_expr() decrements
// and, at <= 0, frees the node and recurses into its children. So each R
// external pointer holds exactly ONE reference: wrap_expr() retains on wrap and
// the registered finalizer (free_expr) releases on GC -- a balanced +1/-1 per
// handle. Because engine constructors retain their own children, freeing the
// last owner cleans up the whole DAG exactly once, in any order.

#include <cpp11.hpp>
#include <string>
#include <vector>

extern "C" {
#include "problem.h"                            // pulls in expr.h
#include "atoms/affine.h"                       // variable, add, neg, sum, shape ops
#include "atoms/elementwise_full_dom.h"         // exp, sin, cos, ... power
#include "atoms/elementwise_restricted_dom.h"   // log, entr, atanh, tan
#include "atoms/bivariate_full_dom.h"           // elementwise_mult, matmul
#include "atoms/bivariate_restricted_dom.h"     // quad_over_lin, rel_entr*
#include "atoms/non_elementwise_full_dom.h"     // prod, prod_axis_*, quad_form
}
#include "sparsediff_cblas.h"             // self-guarded extern "C"

#ifndef DIFF_ENGINE_VERSION
#define DIFF_ENGINE_VERSION "unknown"
#endif

using namespace cpp11;
using namespace cpp11::literals;  // for "name"_nm = ... in list literals

// ---------------------------------------------------------------------------
//  external-pointer helpers for the opaque engine handles
// ---------------------------------------------------------------------------
static void expr_xp_finalizer(SEXP xp) {
  expr* p = static_cast<expr*>(R_ExternalPtrAddr(xp));
  if (p != nullptr) { free_expr(p); R_ClearExternalPtr(xp); }
}
static void problem_xp_finalizer(SEXP xp) {
  problem* p = static_cast<problem*>(R_ExternalPtrAddr(xp));
  if (p != nullptr) { free_problem(p); R_ClearExternalPtr(xp); }
}

static SEXP wrap_expr(expr* node) {
  if (node == nullptr) stop("sparsediff: engine returned a NULL expression");
  expr_retain(node);  // this R handle owns one reference
  SEXP xp = PROTECT(R_MakeExternalPtr(node, R_NilValue, R_NilValue));
  // Cast TRUE: on Windows/Rtools `TRUE` resolves to int (macro pollution), and
  // R_RegisterCFinalizerEx's 3rd arg is Rboolean -> int->Rboolean is a hard
  // error under g++ (-fpermissive). Explicit cast is portable across toolchains.
  R_RegisterCFinalizerEx(xp, expr_xp_finalizer, static_cast<Rboolean>(TRUE));
  UNPROTECT(1);
  return xp;
}
static SEXP wrap_problem(problem* p) {
  if (p == nullptr) stop("sparsediff: engine returned a NULL problem");
  SEXP xp = PROTECT(R_MakeExternalPtr(p, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(xp, problem_xp_finalizer, static_cast<Rboolean>(TRUE));
  UNPROTECT(1);
  return xp;
}
static expr* as_expr(SEXP xp) {
  expr* p = static_cast<expr*>(R_ExternalPtrAddr(xp));
  if (p == nullptr) stop("sparsediff: NULL or finalized expression pointer");
  return p;
}
static problem* as_problem(SEXP xp) {
  problem* p = static_cast<problem*>(R_ExternalPtrAddr(xp));
  if (p == nullptr) stop("sparsediff: NULL or finalized problem pointer");
  return p;
}
static expr* as_expr_or_null(SEXP xp) {
  return (xp == R_NilValue) ? nullptr : as_expr(xp);
}

// A CSR_matrix that ALIASES R's vector memory (zero copy). It is valid only for
// the duration of the call: the engine constructors copy it synchronously (the
// same one-time copy the Python binding incurs), so the aliased R vectors only
// need to survive the .Call, which they do. The engine's CSR uses int p/i and
// double x -- exactly R's INTSXP/REALSXP storage -- so no conversion is needed.
static CSR_matrix csr_view(SEXP p, SEXP i, SEXP x, int ncol) {
  CSR_matrix A;
  A.m = static_cast<int>(Rf_length(p)) - 1;
  A.n = ncol;
  A.nnz = static_cast<int>(Rf_length(x));
  A.p = INTEGER(p);
  A.i = INTEGER(i);
  A.x = REAL(x);
  return A;
}

// ---------------------------------------------------------------------------
//  diagnostics
// ---------------------------------------------------------------------------
[[cpp11::register]]
std::string sd_engine_version() { return std::string(DIFF_ENGINE_VERSION); }

// Self-test of the CBLAS->Fortran-BLAS shim: row-major C(2x2)=A(2x3)B(3x2).
[[cpp11::register]]
doubles sd_selftest_dgemm() {
  const double A[6] = {1, 2, 3, 4, 5, 6};
  const double B[6] = {7, 8, 9, 10, 11, 12};
  double C[4] = {0, 0, 0, 0};
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              2, 2, 3, 1.0, A, 3, B, 2, 0.0, C, 2);
  writable::doubles out(4);
  for (int i = 0; i < 4; i++) out[i] = C[i];
  return out;
}

// ---------------------------------------------------------------------------
//  expression-DAG constructors
//  (var_id is the 0-based offset of the variable's block in the flat input
//   vector u of length n_vars; axis is the engine convention: -1 = all,
//   0 = sum rows, 1 = sum columns.)
// ---------------------------------------------------------------------------
[[cpp11::register]]
SEXP sd_variable(int d1, int d2, int var_id, int n_vars) {
  return wrap_expr(new_variable(d1, d2, var_id, n_vars));
}

[[cpp11::register]]
SEXP sd_exp(SEXP child) { return wrap_expr(new_exp(as_expr(child))); }

[[cpp11::register]]
SEXP sd_neg(SEXP child) { return wrap_expr(new_neg(as_expr(child))); }

[[cpp11::register]]
SEXP sd_add(SEXP left, SEXP right) {
  return wrap_expr(new_add(as_expr(left), as_expr(right)));
}

[[cpp11::register]]
SEXP sd_sum(SEXP child, int axis) {
  return wrap_expr(new_sum(as_expr(child), axis));
}

// ---- elementwise, full domain (unary) ----
[[cpp11::register]] SEXP sd_sin(SEXP c)        { return wrap_expr(new_sin(as_expr(c))); }
[[cpp11::register]] SEXP sd_cos(SEXP c)        { return wrap_expr(new_cos(as_expr(c))); }
[[cpp11::register]] SEXP sd_sinh(SEXP c)       { return wrap_expr(new_sinh(as_expr(c))); }
[[cpp11::register]] SEXP sd_tanh(SEXP c)       { return wrap_expr(new_tanh(as_expr(c))); }
[[cpp11::register]] SEXP sd_asinh(SEXP c)      { return wrap_expr(new_asinh(as_expr(c))); }
[[cpp11::register]] SEXP sd_logistic(SEXP c)   { return wrap_expr(new_logistic(as_expr(c))); }
[[cpp11::register]] SEXP sd_xexp(SEXP c)       { return wrap_expr(new_xexp(as_expr(c))); }
[[cpp11::register]] SEXP sd_normal_cdf(SEXP c) { return wrap_expr(new_normal_cdf(as_expr(c))); }
[[cpp11::register]] SEXP sd_power(SEXP c, double p) { return wrap_expr(new_power(as_expr(c), p)); }

// ---- elementwise, restricted domain (unary) ----
[[cpp11::register]] SEXP sd_log(SEXP c)   { return wrap_expr(new_log(as_expr(c))); }
[[cpp11::register]] SEXP sd_entr(SEXP c)  { return wrap_expr(new_entr(as_expr(c))); }
[[cpp11::register]] SEXP sd_atanh(SEXP c) { return wrap_expr(new_atanh(as_expr(c))); }
[[cpp11::register]] SEXP sd_tan(SEXP c)   { return wrap_expr(new_tan(as_expr(c))); }

// ---- affine: shape / structural ----
[[cpp11::register]] SEXP sd_trace(SEXP c)     { return wrap_expr(new_trace(as_expr(c))); }
[[cpp11::register]] SEXP sd_transpose(SEXP c) { return wrap_expr(new_transpose(as_expr(c))); }
[[cpp11::register]] SEXP sd_diag_vec(SEXP c)  { return wrap_expr(new_diag_vec(as_expr(c))); }
[[cpp11::register]] SEXP sd_diag_mat(SEXP c)  { return wrap_expr(new_diag_mat(as_expr(c))); }
[[cpp11::register]] SEXP sd_upper_tri(SEXP c) { return wrap_expr(new_upper_tri(as_expr(c))); }
[[cpp11::register]] SEXP sd_promote(SEXP c, int d1, int d2)   { return wrap_expr(new_promote(as_expr(c), d1, d2)); }
[[cpp11::register]] SEXP sd_reshape(SEXP c, int d1, int d2)   { return wrap_expr(new_reshape(as_expr(c), d1, d2)); }
[[cpp11::register]] SEXP sd_broadcast(SEXP c, int d1, int d2) { return wrap_expr(new_broadcast(as_expr(c), d1, d2)); }

// indices are 0-based offsets into the (column-major) flattened child.
[[cpp11::register]]
SEXP sd_index(SEXP child, int d1, int d2, integers indices) {
  std::vector<int> idx(indices.begin(), indices.end());
  return wrap_expr(new_index(as_expr(child), d1, d2, idx.data(),
                             static_cast<int>(idx.size())));
}

[[cpp11::register]]
SEXP sd_hstack(list args, int n_vars) {
  std::vector<expr*> a;
  a.reserve(args.size());
  for (R_xlen_t i = 0; i < args.size(); i++) { SEXP s = args[i]; a.push_back(as_expr(s)); }
  return wrap_expr(new_hstack(a.data(), static_cast<int>(a.size()), n_vars));
}
[[cpp11::register]]
SEXP sd_vstack(list args, int n_vars) {
  std::vector<expr*> a;
  a.reserve(args.size());
  for (R_xlen_t i = 0; i < args.size(); i++) { SEXP s = args[i]; a.push_back(as_expr(s)); }
  return wrap_expr(new_vstack(a.data(), static_cast<int>(a.size()), n_vars));
}

// ---- bivariate, full domain ----
[[cpp11::register]] SEXP sd_elementwise_mult(SEXP l, SEXP r) { return wrap_expr(new_elementwise_mult(as_expr(l), as_expr(r))); }
[[cpp11::register]] SEXP sd_matmul(SEXP x, SEXP y)           { return wrap_expr(new_matmul(as_expr(x), as_expr(y))); }

// ---- bivariate, restricted domain ----
[[cpp11::register]] SEXP sd_quad_over_lin(SEXP l, SEXP r)          { return wrap_expr(new_quad_over_lin(as_expr(l), as_expr(r))); }
[[cpp11::register]] SEXP sd_rel_entr(SEXP l, SEXP r)               { return wrap_expr(new_rel_entr_vector_args(as_expr(l), as_expr(r))); }
[[cpp11::register]] SEXP sd_rel_entr_first_scalar(SEXP l, SEXP r)  { return wrap_expr(new_rel_entr_first_arg_scalar(as_expr(l), as_expr(r))); }
[[cpp11::register]] SEXP sd_rel_entr_second_scalar(SEXP l, SEXP r) { return wrap_expr(new_rel_entr_second_arg_scalar(as_expr(l), as_expr(r))); }

// ---- non-elementwise, full domain ----
[[cpp11::register]] SEXP sd_prod(SEXP c)           { return wrap_expr(new_prod(as_expr(c))); }
[[cpp11::register]] SEXP sd_prod_axis_zero(SEXP c) { return wrap_expr(new_prod_axis_zero(as_expr(c))); }
[[cpp11::register]] SEXP sd_prod_axis_one(SEXP c)  { return wrap_expr(new_prod_axis_one(as_expr(c))); }

// ---------------------------------------------------------------------------
//  parameters / constants and the atoms that consume them
//  param_id: PARAM_FIXED (-1) for a fixed constant; >= 0 = offset into the theta
//  parameter vector for an updatable parameter. values has length d1*d2
//  (column-major) and is copied by the engine.
// ---------------------------------------------------------------------------
[[cpp11::register]]
SEXP sd_parameter(int d1, int d2, int param_id, int n_vars, SEXP values) {
  return wrap_expr(new_parameter(d1, d2, param_id, n_vars, REAL(values)));
}

// a (a parameter/constant node) combined with child:
[[cpp11::register]] SEXP sd_scalar_mult(SEXP param, SEXP child) { return wrap_expr(new_scalar_mult(as_expr(param), as_expr(child))); }
[[cpp11::register]] SEXP sd_vector_mult(SEXP param, SEXP child) { return wrap_expr(new_vector_mult(as_expr(param), as_expr(child))); }
[[cpp11::register]] SEXP sd_convolve(SEXP param, SEXP child)    { return wrap_expr(new_convolve(as_expr(param), as_expr(child))); }

// quadratic form  x' Q x  (Q square, full, symmetric; CSR components, zero-copy).
[[cpp11::register]]
SEXP sd_quad_form(SEXP child, SEXP Qp, SEXP Qi, SEXP Qx) {
  int n = static_cast<int>(Rf_length(Qp)) - 1;
  CSR_matrix Q = csr_view(Qp, Qi, Qx, n);
  return wrap_expr(new_quad_form(as_expr(child), &Q));  // engine copies Q
}

// constant sparse-matrix products  A @ f(x)  and  f(x) @ A  (A is m x ncol CSR).
[[cpp11::register]]
SEXP sd_left_matmul(SEXP child, SEXP Ap, SEXP Ai, SEXP Ax, int ncol) {
  CSR_matrix A = csr_view(Ap, Ai, Ax, ncol);
  return wrap_expr(new_left_matmul(nullptr, as_expr(child), &A));
}
[[cpp11::register]]
SEXP sd_right_matmul(SEXP child, SEXP Ap, SEXP Ai, SEXP Ax, int ncol) {
  CSR_matrix A = csr_view(Ap, Ai, Ax, ncol);
  return wrap_expr(new_right_matmul(nullptr, as_expr(child), &A));
}

// dense-matrix products. param = NULL with row-major `data` for a constant
// matrix; or a parameter node with empty `data` for a parametric matrix.
[[cpp11::register]]
SEXP sd_left_matmul_dense(SEXP param, SEXP child, int m, int n, SEXP data) {
  const double* dptr = (Rf_length(data) == 0) ? nullptr : REAL(data);
  return wrap_expr(new_left_matmul_dense(as_expr_or_null(param), as_expr(child), m, n, dptr));
}
[[cpp11::register]]
SEXP sd_right_matmul_dense(SEXP param, SEXP child, int m, int n, SEXP data) {
  const double* dptr = (Rf_length(data) == 0) ? nullptr : REAL(data);
  return wrap_expr(new_right_matmul_dense(as_expr_or_null(param), as_expr(child), m, n, dptr));
}

// ---------------------------------------------------------------------------
//  problem construction & evaluation
// ---------------------------------------------------------------------------
[[cpp11::register]]
SEXP sd_problem(SEXP objective, list constraints, bool verbose) {
  int nc = static_cast<int>(constraints.size());
  std::vector<expr*> cons;
  cons.reserve(nc);
  for (int i = 0; i < nc; i++) {
    SEXP s = constraints[i];
    cons.push_back(as_expr(s));
  }
  return wrap_problem(
      new_problem(as_expr(objective), nc ? cons.data() : nullptr, nc, verbose));
}

// Register the problem's updatable parameter nodes (each created by
// sd_parameter with param_id >= 0). The problem keeps weak references; the nodes
// stay alive through the expression DAG and their R handles.
[[cpp11::register]]
void sd_register_params(SEXP prob, list params) {
  problem* p = as_problem(prob);
  int n = static_cast<int>(params.size());
  std::vector<expr*> arr;
  arr.reserve(n);
  for (R_xlen_t k = 0; k < params.size(); k++) { SEXP s = params[k]; arr.push_back(as_expr(s)); }
  problem_register_params(p, n ? arr.data() : nullptr, n);
}

// Update parameter values from the concatenated theta vector (offsets = param_id).
[[cpp11::register]]
void sd_update_params(SEXP prob, SEXP theta) {
  problem_update_params(as_problem(prob), REAL(theta));
}

[[cpp11::register]]
void sd_init_jacobian(SEXP prob) { problem_init_jacobian(as_problem(prob)); }

[[cpp11::register]]
void sd_init_derivatives(SEXP prob) { problem_init_derivatives(as_problem(prob)); }

[[cpp11::register]]
double sd_objective_forward(SEXP prob, SEXP u) {
  return problem_objective_forward(as_problem(prob), REAL(u));
}

[[cpp11::register]]
doubles sd_gradient(SEXP prob) {
  problem* p = as_problem(prob);
  problem_gradient(p);
  writable::doubles g(p->n_vars);
  for (int i = 0; i < p->n_vars; i++) g[i] = p->gradient_values[i];
  return g;
}

[[cpp11::register]]
doubles sd_constraint_forward(SEXP prob, SEXP u) {
  problem* p = as_problem(prob);
  problem_constraint_forward(p, REAL(u));
  writable::doubles cv(p->total_constraint_size);
  for (int k = 0; k < p->total_constraint_size; k++) cv[k] = p->constraint_values[k];
  return cv;
}

// ---------------------------------------------------------------------------
//  sparse derivatives in COO form
//
//  Row/column indices are 0-based (engine convention; the higher-level R layer
//  translates to 1-based). Sparsity is structural and fixed after init; values
//  must be recomputed (objective_forward / constraint_forward populate node
//  values first). The constraint Jacobian COO and the lower-triangular Lagrange
//  Hessian COO mirror CVXPY's diff_engine usage.
// ---------------------------------------------------------------------------
static list coo_sparsity(const COO_matrix* coo) {
  writable::integers rows(coo->nnz), cols(coo->nnz);
  for (int k = 0; k < coo->nnz; k++) { rows[k] = coo->rows[k]; cols[k] = coo->cols[k]; }
  return writable::list({"rows"_nm = rows, "cols"_nm = cols,
                         "nrow"_nm = coo->m, "ncol"_nm = coo->n});
}

// ---- constraint Jacobian ----
[[cpp11::register]]
void sd_init_jacobian_coo(SEXP prob) { problem_init_jacobian_coo(as_problem(prob)); }

[[cpp11::register]]
list sd_jacobian_sparsity(SEXP prob) {
  problem* p = as_problem(prob);
  if (p->jacobian_coo == nullptr) stop("sparsediff: call sd_init_jacobian_coo() first");
  return coo_sparsity(p->jacobian_coo);
}

[[cpp11::register]]
doubles sd_jacobian_values(SEXP prob) {
  problem* p = as_problem(prob);
  if (p->jacobian_coo == nullptr) stop("sparsediff: call sd_init_jacobian_coo() first");
  problem_jacobian(p);  // fills CSR; COO order matches CSR->x order
  int nnz = p->jacobian_coo->nnz;
  writable::doubles vals(nnz);
  for (int k = 0; k < nnz; k++) vals[k] = p->jacobian->x[k];
  return vals;
}

// ---- lower-triangular Lagrange Hessian ----
[[cpp11::register]]
void sd_init_hessian_coo(SEXP prob) {
  problem_init_hessian_coo_lower_triangular(as_problem(prob));
}

[[cpp11::register]]
list sd_hessian_sparsity(SEXP prob) {
  problem* p = as_problem(prob);
  if (p->lagrange_hessian_coo == nullptr) stop("sparsediff: call sd_init_hessian_coo() first");
  return coo_sparsity(p->lagrange_hessian_coo);
}

// obj_w scales the objective Hessian; w (length = total constraint size) scales
// the constraint Hessians: H = obj_w * d2f + sum_i w_i d2g_i.
[[cpp11::register]]
doubles sd_hessian_values(SEXP prob, double obj_w, doubles w) {
  problem* p = as_problem(prob);
  COO_matrix* coo = p->lagrange_hessian_coo;
  if (coo == nullptr) stop("sparsediff: call sd_init_hessian_coo() first");
  std::vector<double> wbuf(w.begin(), w.end());
  problem_hessian(p, obj_w, wbuf.empty() ? nullptr : wbuf.data());
  refresh_lower_triangular_coo(coo, p->lagrange_hessian->x);
  writable::doubles vals(coo->nnz);
  for (int k = 0; k < coo->nnz; k++) vals[k] = coo->x[k];
  return vals;
}
