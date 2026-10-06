# sparsediff 0.6.1

* The bundled SparseDiffEngine moves from 0.3.0 to 0.6.1, the release that
  'sparsediffpy' 0.6.x and 'CVXPY' 1.9.3 use. From this release on, the package
  version matches the bundled engine version.
* Fixed: a parametric problem re-evaluated after `sd_update_params()` could
  return values computed from the old parameter values (engine fix #107). For a
  modeling layer this meant that re-solving a parametric nonlinear problem
  after changing a parameter returned the previous solution, reported as
  optimal.
* Fixed: selecting entries with repeated indices in `sd_index()` (for example a
  quadratic objective over a symmetric variable) could crash R with a bus error
  when the Jacobian or Hessian was evaluated (engine fix #105).
* New `sd_quad_form_dense()`: the quadratic form with a dense constant or
  parametric matrix. The parametric form follows `sd_update_params()`.
  `sd_quad_form()` keeps its signature and still takes a sparse constant matrix.
* Every 'sparsediffpy' 0.6.1 function now has an `sd_*` counterpart that calls
  the same engine routine with the same arguments; the name map is on the
  package help page (`?sparsediff`). New: `sd_left_kron()` and
  `sd_right_kron()` (Kronecker products, used by the next CVXPY release), the
  compressed-sparse-row derivative accessors `sd_jacobian()`,
  `sd_get_jacobian()`, `sd_init_hessian()`, `sd_hessian()` and
  `sd_get_hessian()`, and `sd_get_expr_dimensions()` / `sd_get_expr_size()`.
  `sd_rel_entr()` now dispatches on operand size like `make_rel_entr`.
* Fixed: `sd_init_hessian_coo()` crashed R when called before the Jacobian was
  initialized; it now initializes the Jacobian first.
* `sd_rel_entr()`, `sd_rel_entr_first_scalar()` and `sd_rel_entr_second_scalar()`
  now check that their arguments are two different variables, which the
  engine requires; other arguments overflowed a heap buffer (found with
  AddressSanitizer).
* The BLAS shim gains `cblas_ddot`, which the 0.6.1 engine calls.
* Fixed: the package failed to compile on Linux systems with glibc 2.36 or
  older (for example Debian 12), because `_GNU_SOURCE` was defined too late for
  `fopencookie()`. It is now set on the compiler command line.
* Documentation: the sparse matrix arguments of `sd_quad_form()` and
  `sd_left_matmul()` / `sd_right_matmul()` are compressed-sparse-row arrays,
  and the dense `data` arguments are row-major. Earlier documentation said
  compressed-sparse-column and column-major. The grouped help pages now have
  usage sections (flagged by current R-devel).
* The package now has a 'testthat' test suite.

# sparsediff 0.4.0

* First CRAN release, bundling SparseDiffEngine 0.3.0.
