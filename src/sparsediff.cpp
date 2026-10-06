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
#include <cmath>
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

// Node shape getters (sparsediffpy get_expr_dimensions / get_expr_size).
[[cpp11::register]]
integers sd_get_expr_dimensions(SEXP node) {
  expr* e = as_expr(node);
  writable::integers d({e->d1, e->d2});
  return d;
}
[[cpp11::register]]
int sd_get_expr_size(SEXP node) { return as_expr(node)->size; }

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
// All three relative-entropy atoms require both arguments to be DISTINCT
// variable leaves (CVXPY's canonicalization guarantees it). The engine only
// asserts this (compiled out under NDEBUG); violating it overflows the heap
// while the derivative structure is built (seen under ASAN), so check here.
static void check_rel_entr_args(const char* who, const expr* a, const expr* b) {
  if (a->var_id == NOT_A_VARIABLE || b->var_id == NOT_A_VARIABLE)
    stop("%s: both arguments must be variables (sd_variable)", who);
  if (a->var_id == b->var_id) stop("%s: the two arguments must be different variables", who);
}
// sd_rel_entr dispatches on operand size exactly as sparsediffpy make_rel_entr:
// scalar first argument, scalar second argument, or elementwise.
[[cpp11::register]] SEXP sd_rel_entr(SEXP l, SEXP r) {
  expr* a = as_expr(l);
  expr* b = as_expr(r);
  check_rel_entr_args("sd_rel_entr", a, b);
  if (a->size == 1 && b->size > 1) return wrap_expr(new_rel_entr_first_arg_scalar(a, b));
  if (a->size > 1 && b->size == 1) return wrap_expr(new_rel_entr_second_arg_scalar(a, b));
  return wrap_expr(new_rel_entr_vector_args(a, b));
}
[[cpp11::register]] SEXP sd_rel_entr_first_scalar(SEXP l, SEXP r) {
  expr* a = as_expr(l);
  expr* b = as_expr(r);
  check_rel_entr_args("sd_rel_entr_first_scalar", a, b);
  return wrap_expr(new_rel_entr_first_arg_scalar(a, b));
}
[[cpp11::register]] SEXP sd_rel_entr_second_scalar(SEXP l, SEXP r) {
  expr* a = as_expr(l);
  expr* b = as_expr(r);
  check_rel_entr_args("sd_rel_entr_second_scalar", a, b);
  return wrap_expr(new_rel_entr_second_arg_scalar(a, b));
}

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
  return wrap_expr(new_quad_form_sparse(as_expr(child), &Q));  // engine copies Q
}

// quadratic form  x' Q x  with a DENSE n x n Q, n = length of the vector child.
// Exactly one source: param = NULL with row-major `data` (n*n doubles) for a
// constant Q; or a parameter node (size n*n) with empty `data` for a parametric
// Q, refreshed from the parameter on every forward pass. The engine copies
// `data`. Q must be symmetric (the engine's gradient is 2 Q x): checked here for
// constant data, the caller's contract for a parameter. All validation happens
// before the engine is called, because R builds define NDEBUG (the engine's
// asserts are compiled out) and an engine-side error longjmps out of the call.
[[cpp11::register]]
SEXP sd_quad_form_dense(SEXP param, SEXP child, SEXP data) {
  expr* c = as_expr(child);
  if (c->d1 != 1 && c->d2 != 1)
    stop("sd_quad_form_dense: child must be a vector, not %d x %d", c->d1, c->d2);
  const int n = c->size;
  const bool has_param = (param != R_NilValue);
  const R_xlen_t len = Rf_xlength(data);
  if (has_param == (len > 0))
    stop("sd_quad_form_dense: supply exactly one of `param` and non-empty `data`");
  if (has_param) {
    expr* q = as_expr(param);
    if (q->size != n * n)
      stop("sd_quad_form_dense: parameter has %d entries, expected n*n = %d", q->size, n * n);
    return wrap_expr(new_quad_form_dense(c, n, nullptr, q));
  }
  if (TYPEOF(data) != REALSXP) stop("sd_quad_form_dense: `data` must be a double vector");
  if (len != static_cast<R_xlen_t>(n) * n)
    stop("sd_quad_form_dense: `data` has length %d, expected n*n = %d", static_cast<int>(len), n * n);
  const double* Q = REAL(data);
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      const double a = Q[i * n + j], b = Q[j * n + i];
      if (std::fabs(a - b) > 1e-8 * (1.0 + std::fmax(std::fabs(a), std::fabs(b))))
        stop("sd_quad_form_dense: Q is not symmetric (entries [%d,%d] and [%d,%d] differ)",
             i + 1, j + 1, j + 1, i + 1);
    }
  }
  return wrap_expr(new_quad_form_dense(c, n, Q, nullptr));
}

// Kronecker products  Z = kron(A, B)  with one variable-free operand `param`
// (a parameter, or a constant made with sd_parameter(..., param_id = -1, ...))
// and one variable operand `child`. Mirrors sparsediffpy make_left_kron /
// make_right_kron. left: A = param (p x q), B = child (r x s); right: A = child
// (p x q), B = param (r x s). active_blocks are the 0-based column-major indices
// of the constant operand's nonzero entries (all of them for a parameter); only
// the output rows they cover are materialized. Validated here because the
// engine only asserts (compiled out under NDEBUG).
static SEXP kron_impl(SEXP param, SEXP child, int p, int q, int r, int s,
                      SEXP active_blocks, bool is_left) {
  const char* who = is_left ? "sd_left_kron" : "sd_right_kron";
  if (p < 1 || q < 1 || r < 1 || s < 1) stop("%s: p, q, r, s must be positive", who);
  expr* k = as_expr(param);
  expr* c = as_expr(child);
  const int kd1 = is_left ? p : r, kd2 = is_left ? q : s;  // constant operand
  const int cd1 = is_left ? r : p, cd2 = is_left ? s : q;  // variable operand
  if (k->size != kd1 * kd2)
    stop("%s: `param` has %d entries, expected %d x %d", who, k->size, kd1, kd2);
  if (c->size != cd1 * cd2)
    stop("%s: `child` has %d entries, expected %d x %d", who, c->size, cd1, cd2);
  if (TYPEOF(active_blocks) != INTSXP) stop("%s: `active_blocks` must be an integer vector", who);
  const int n_active = static_cast<int>(Rf_xlength(active_blocks));
  const int* ab = INTEGER(active_blocks);
  for (int b = 0; b < n_active; b++) {
    if (ab[b] == NA_INTEGER || ab[b] < 0 || ab[b] >= kd1 * kd2)
      stop("%s: active_blocks[%d] = %d is outside [0, %d)", who, b + 1, ab[b], kd1 * kd2);
  }
  return wrap_expr(is_left ? new_left_kron(k, c, p, q, r, s, ab, n_active)
                           : new_right_kron(k, c, p, q, r, s, ab, n_active));
}
[[cpp11::register]]
SEXP sd_left_kron(SEXP param, SEXP child, int p, int q, int r, int s, SEXP active_blocks) {
  return kron_impl(param, child, p, q, r, s, active_blocks, true);
}
[[cpp11::register]]
SEXP sd_right_kron(SEXP param, SEXP child, int p, int q, int r, int s, SEXP active_blocks) {
  return kron_impl(param, child, p, q, r, s, active_blocks, false);
}

// constant sparse-matrix products  A @ f(x)  and  f(x) @ A  (A is m x ncol CSR).
// sparsediffpy's sparse form also takes a parameter, but engine 0.6.1 exits on
// a parameter for a sparse matrix, so it is not exposed here; a parametric
// matrix goes through sd_left_matmul_dense / sd_right_matmul_dense.
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
  problem* p = as_problem(prob);
  problem_init_jacobian(p);  // see sd_init_hessian: the Hessian init needs it
  problem_init_hessian_coo_lower_triangular(p);
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

// ---------------------------------------------------------------------------
//  sparse derivatives in CSR form (sparsediffpy problem_jacobian / get_jacobian /
//  problem_init_hessian / problem_hessian / get_hessian). Each returns
//  list(data, indices, indptr, shape): the CSR arrays, 0-based, as in
//  scipy.sparse.csr_array((data, indices, indptr), shape). The Jacobian is the
//  stacked constraint Jacobian; the Hessian is the FULL (both triangles)
//  Lagrangian Hessian obj_w * d2f + sum_i w_i d2g_i.
// ---------------------------------------------------------------------------
static list csr_list(const CSR_matrix* A, int indptr_len) {
  writable::doubles data(A->nnz);
  writable::integers indices(A->nnz), indptr(indptr_len);
  for (int k = 0; k < A->nnz; k++) { data[k] = A->x[k]; indices[k] = A->i[k]; }
  for (int k = 0; k < indptr_len; k++) indptr[k] = A->p[k];
  writable::integers shape({A->m, A->n});
  return writable::list({"data"_nm = data, "indices"_nm = indices,
                         "indptr"_nm = indptr, "shape"_nm = shape});
}

// Evaluate the constraint Jacobian at the point of the last forward pass.
[[cpp11::register]]
list sd_jacobian(SEXP prob) {
  problem* p = as_problem(prob);
  if (p->n_constraints == 0) {
    writable::doubles data(static_cast<R_xlen_t>(0));
    writable::integers indices(static_cast<R_xlen_t>(0)), indptr({0});
    writable::integers shape({0, p->n_vars});
    return writable::list({"data"_nm = data, "indices"_nm = indices,
                           "indptr"_nm = indptr, "shape"_nm = shape});
  }
  if (p->jacobian == nullptr) stop("sparsediff: call sd_init_jacobian() first");
  problem_jacobian(p);
  return csr_list(p->jacobian, p->jacobian->m + 1);
}

// The Jacobian as last evaluated, without re-evaluating it.
[[cpp11::register]]
list sd_get_jacobian(SEXP prob) {
  problem* p = as_problem(prob);
  if (p->jacobian == nullptr) stop("sparsediff: jacobian not initialized - call sd_jacobian() first");
  return csr_list(p->jacobian, p->jacobian->m + 1);
}

// The engine's Hessian init reads the Jacobian structure and segfaults if it
// does not exist yet (problem_init_derivatives orders the two for this reason),
// so initialize the Jacobian first; it is idempotent and purely structural.
[[cpp11::register]]
void sd_init_hessian(SEXP prob) {
  problem* p = as_problem(prob);
  problem_init_jacobian(p);
  problem_init_hessian(p);
}

// Evaluate the Lagrangian Hessian; w has length = total constraint size.
[[cpp11::register]]
list sd_hessian(SEXP prob, double obj_w, doubles w) {
  problem* p = as_problem(prob);
  if (p->lagrange_hessian == nullptr) stop("sparsediff: call sd_init_hessian() first");
  if (static_cast<int>(w.size()) != p->total_constraint_size)
    stop("sd_hessian: `w` has length %d, expected the total constraint size %d",
         static_cast<int>(w.size()), p->total_constraint_size);
  std::vector<double> wbuf(w.begin(), w.end());
  problem_hessian(p, obj_w, wbuf.empty() ? nullptr : wbuf.data());
  return csr_list(p->lagrange_hessian, p->lagrange_hessian->n + 1);
}

// The Hessian as last evaluated, without re-evaluating it.
[[cpp11::register]]
list sd_get_hessian(SEXP prob) {
  problem* p = as_problem(prob);
  if (p->lagrange_hessian == nullptr) stop("sparsediff: hessian not initialized - call sd_hessian() first");
  return csr_list(p->lagrange_hessian, p->lagrange_hessian->n + 1);
}
