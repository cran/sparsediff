/* cblas_shim.c
 *
 * Implements the 4 CBLAS routines SparseDiffEngine uses (cblas_dscal,
 * cblas_dcopy, cblas_dgemv, cblas_dgemm) by forwarding to R's Fortran BLAS via
 * F77_CALL(). This is the standard CBLAS-over-Fortran-BLAS wrapper (cf. Netlib
 * reference CBLAS), so the engine uses whatever BLAS R is configured with
 * (R's reference BLAS by default), rather than a separate CBLAS implementation.
 *
 * Row-major <-> column-major translation follows the Netlib reference CBLAS:
 * the Fortran BLAS is column-major only, so a row-major operation is computed
 * as its column-major transpose-equivalent (swap operands / flip transpose).
 *
 * USE_FC_LEN_T + FCONE pass the hidden Fortran character-length arguments
 * (required for correctness with modern compilers / LTO; see WRE "Fortran
 * character strings").
 */
#ifndef USE_FC_LEN_T
# define USE_FC_LEN_T
#endif
#include <Rconfig.h>       /* FC_LEN_T */
#include <R_ext/BLAS.h>    /* F77_CALL(dgemm/dgemv/dscal/dcopy), F77_NAME */
#ifndef FCONE
# define FCONE
#endif

#include "sparsediff_cblas.h"

static char trans_char(enum CBLAS_TRANSPOSE t) {
    switch (t) {
        case CblasTrans:     return 'T';
        case CblasConjTrans: return 'C';   /* real matrices: same as 'T' */
        case CblasNoTrans:
        default:             return 'N';
    }
}

void cblas_dscal(const int N, const double alpha, double *X, const int incX) {
    F77_CALL(dscal)(&N, &alpha, X, &incX);
}

void cblas_dcopy(const int N, const double *X, const int incX,
                 double *Y, const int incY) {
    F77_CALL(dcopy)(&N, X, &incX, Y, &incY);
}

void cblas_dgemv(const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA,
                 const int M, const int N, const double alpha,
                 const double *A, const int lda,
                 const double *X, const int incX,
                 const double beta, double *Y, const int incY) {
    if (order == CblasColMajor) {
        char ta = trans_char(TransA);
        F77_CALL(dgemv)(&ta, &M, &N, &alpha, A, &lda, X, &incX, &beta, Y, &incY FCONE);
    } else {
        /* Row-major A (M x N, lda = N) is column-major A^T (N x M, ld = lda).
         * Flip the transpose flag and swap M <-> N. */
        char ta = (TransA == CblasNoTrans) ? 'T' : 'N';
        F77_CALL(dgemv)(&ta, &N, &M, &alpha, A, &lda, X, &incX, &beta, Y, &incY FCONE);
    }
}

void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB,
                 const int M, const int N, const int K, const double alpha,
                 const double *A, const int lda,
                 const double *B, const int ldb,
                 const double beta, double *C, const int ldc) {
    char ta = trans_char(TransA);
    char tb = trans_char(TransB);
    if (Order == CblasColMajor) {
        F77_CALL(dgemm)(&ta, &tb, &M, &N, &K, &alpha, A, &lda, B, &ldb,
                        &beta, C, &ldc FCONE FCONE);
    } else {
        /* Row-major C = op(A) op(B) is computed column-major as
         * C^T = op(B)^T op(A)^T: swap A <-> B, swap M <-> N, swap transA/transB. */
        F77_CALL(dgemm)(&tb, &ta, &N, &M, &K, &alpha, B, &ldb, A, &lda,
                        &beta, C, &ldc FCONE FCONE);
    }
}
