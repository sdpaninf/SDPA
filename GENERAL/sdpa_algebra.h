/*  sdpa_algebra.h

LAPACK+BLAS definitions wrapper

Define macros to mangle the given C identifier (in lower and upper
case), which must not contain underscores, for linking with Fortran.

*/

#ifndef __sdpa_algebra_h__
#define __sdpa_algebra_h__

#define F77_RET_I long int
#define F77_RET_D double

#if defined(__APPLE__) // Dirty...
#define F77_FUNC(name,NAME) name ## _
#endif

#define dtrsm_f77  F77_FUNC (dtrsm, DTRSM)
#define dsyrk_f77  F77_FUNC (dsyrk, DSYRK)
#define dcopy_f77  F77_FUNC (dcopy, DCOPY)
#define daxpy_f77  F77_FUNC (daxpy, DAXPY)
#define dgemm_f77  F77_FUNC (dgemm, DGEMM)
#define dgemv_f77  F77_FUNC (dgemv, DGEMV)
#define dscal_f77  F77_FUNC (dscal, DSCAL)
#define dtrsv_f77  F77_FUNC (dtrsv, DTRSV)
#define dtrmv_f77  F77_FUNC (dtrmv, DTRMV)
#define ddot_f77   F77_FUNC (ddot, DDOT)
#define dtrmm_f77  F77_FUNC (dtrmm, DTRMM)
#define ilaenv_f77 F77_FUNC (ilaenv, ILAENV)
#define dsteqr_f77 F77_FUNC (dsteqr, DSTEQR)
#define dsyev_f77  F77_FUNC (dsyev, DSYEV)
#define dpotrf_f77 F77_FUNC (dpotrf, DPORTRF)


extern "C"
{
// BLAS
  F77_RET_I  dtrsm_f77
      (char* side, char* uplo, char* trans, char* diag,
       SDPA_INT* M, SDPA_INT* N,
       double* alpha,
       double* A, SDPA_INT* lda,
       double* B, SDPA_INT* ldb, SDPA_INT side_len,
       SDPA_INT uplo_len, SDPA_INT trans_len, SDPA_INT diag_len);

  F77_RET_I  dsyrk_f77
      (char* uplo, char* trans, SDPA_INT* N, SDPA_INT* K,
       double* alpha,
       double* A, SDPA_INT* lda,
       double* beta,
       double* C, SDPA_INT* ldc, SDPA_INT uplo_len, SDPA_INT trans_len);

  F77_RET_I  dcopy_f77
      (SDPA_INT* N,
       double* X, SDPA_INT* incX,
       double* Y, SDPA_INT* incY);

  F77_RET_I  daxpy_f77
      (SDPA_INT* N,
       double* alpha,
       double* X, SDPA_INT* incX,
       double* Y, SDPA_INT* incY);

  F77_RET_I  dgemm_f77
      (char* transA, char* transB, SDPA_INT* M, SDPA_INT* N, SDPA_INT* K,
       double* alpha,
       double* A, SDPA_INT* lda,
       double* B, SDPA_INT* ldb,
       double* beta,
       double* C, SDPA_INT* ldc, SDPA_INT transA_len, SDPA_INT transB_len);

  F77_RET_I  dgemv_f77
      (char* trans, SDPA_INT* M, SDPA_INT* N,
       double* alpha,
       double* A, SDPA_INT* lda,
       double* X, SDPA_INT* incX,
       double* beta,
       double* Y, SDPA_INT* incY, SDPA_INT trans_len);

  F77_RET_I  dscal_f77
      (SDPA_INT* N,
       double* alpha,
       double* X, SDPA_INT* incX);

  F77_RET_I  dtrsv_f77
      (char* uplo, char* trans, char* diag, SDPA_INT* N,
       double* A, SDPA_INT* lda,
       double* X, SDPA_INT* incX, SDPA_INT uplo_len,
       SDPA_INT trans_len, SDPA_INT diag_len);

  F77_RET_I  dtrmv_f77
      (char* uplo, char *trans, char* diag, SDPA_INT *N,  
       double *A, SDPA_INT *lda, 
       double *X, SDPA_INT *incX, SDPA_INT uplo_len, SDPA_INT trans_len, SDPA_INT diag_len);

  F77_RET_D  ddot_f77
      (SDPA_INT* N, double* X, SDPA_INT* incX, double* Y, SDPA_INT* incY);

  F77_RET_I  dtrmm_f77
      (char* side, char* uplo, char* trans, char* diag, 
       SDPA_INT* M, SDPA_INT* N,
       double* alpha,
       double* A, SDPA_INT* lda,
       double* B, SDPA_INT* ldb, SDPA_INT side_len, SDPA_INT uplo_len,
       SDPA_INT trans_len, SDPA_INT diag_len);

// LAPACK

  F77_RET_I  ilaenv_f77
      (SDPA_INT *ispec, char *name, char *opts, SDPA_INT *n1, 
	SDPA_INT *n2, SDPA_INT *n3, SDPA_INT *n4, SDPA_INT name_len, SDPA_INT opts_len);

  F77_RET_I  dsteqr_f77
      (char *compz, SDPA_INT *n, double *d, 
	double *e, double *z, SDPA_INT *ldz, double *work, 
	SDPA_INT *info, SDPA_INT compz_len);

  F77_RET_I  dsyev_f77
      (char *jobz, char *uplo, SDPA_INT *n, double *a,
        SDPA_INT *lda, double *w, double *work, SDPA_INT *lwork, 
	SDPA_INT *info, SDPA_INT jobz_len, SDPA_INT uplo_len);

  F77_RET_I  dpotrf_f77
     (char *uplo, SDPA_INT *n, double *a, SDPA_INT *lda,
      SDPA_INT *info, SDPA_INT uplo_len);
}

#endif // __sdpa_algebra_h__
