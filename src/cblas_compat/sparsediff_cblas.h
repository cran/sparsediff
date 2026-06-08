/* sparsediff_cblas.h
 *
 * Minimal CBLAS declarations for the SparseDiffEngine "RBLAS" build path
 * (engine compiled with -DSPARSEDIFF_USE_RBLAS; see cblas_wrapper.h).
 *
 * The engine is written against the CBLAS C interface, but R guarantees only
 * the *Fortran* BLAS (R's bundled reference BLAS ships no CBLAS symbols). This
 * header declares the small CBLAS subset the engine actually uses; the matching
 * implementations in cblas_shim.c forward to R's Fortran BLAS via F77_CALL().
 * This is the standard "CBLAS over Fortran BLAS" wrapper, restricted to the 4
 * routines SparseDiffEngine calls, so the engine uses the same BLAS as R.
 *
 * Enum values and signatures match the Netlib reference cblas.h so the engine's
 * call sites compile and behave identically.
 */
#ifndef SPARSEDIFF_CBLAS_H
#define SPARSEDIFF_CBLAS_H

#ifdef __cplusplus
extern "C" {
#endif

#define CBLAS_INDEX int

enum CBLAS_ORDER     { CblasRowMajor = 101, CblasColMajor = 102 };
enum CBLAS_TRANSPOSE { CblasNoTrans  = 111, CblasTrans    = 112, CblasConjTrans = 113 };

void cblas_dscal(const int N, const double alpha, double *X, const int incX);

void cblas_dcopy(const int N, const double *X, const int incX,
                 double *Y, const int incY);

void cblas_dgemv(const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA,
                 const int M, const int N, const double alpha,
                 const double *A, const int lda,
                 const double *X, const int incX,
                 const double beta, double *Y, const int incY);

void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB,
                 const int M, const int N, const int K, const double alpha,
                 const double *A, const int lda,
                 const double *B, const int ldb,
                 const double beta, double *C, const int ldc);

#ifdef __cplusplus
}
#endif

#endif /* SPARSEDIFF_CBLAS_H */
