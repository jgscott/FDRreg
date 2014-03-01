/*
  This file contains the functions that make use of MatrixFrame.
  These include Blas wrappers for matrix arithmetic, LAPACK
  wrappers, and variouis matrix manipulations.

  You probably want to view this file with at least a 160 character
  wide screen.

  The BLAStoC.perl file will help you write BLAS wrappers.
 */

#include "MatrixFrame.h"

//////////////////////////////////////////////////////////////////////
			 // BLAS / LAPACK //
//////////////////////////////////////////////////////////////////////

/*
  See BLAS / LAPACK documentation at netlib.org.  The Fortran source
  code is very regular.  The perl script BLAStoC.perl will take a BLAS
  subroutine/function file and convert it to the necessary C code.
 */

//////////////////////////////////////////////////////////////////////
		    // MatrixFrame BLAS WRAPPER //
//////////////////////////////////////////////////////////////////////

// DOUBLE

void raxpy(int n, double da, double* dx, int incx, double* dy, int incy)
{ daxpy_(&n, &da, dx, &incx, dy, &incy); }

double rdot(int n, double* dx, int incx, double* dy, int incy)
{ return ddot_(&n, dx, &incx, dy, &incy); }

//------------------------------------------------------------------------------

void rgemm(char transa, char transb, int m, int n, int k, double alpha, double* a, int lda, double* b, int ldb, double beta, double* c, int ldc)
{ dgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); }

void rsyr2k(char uplo, char trans, int n, int k, double alpha, double* a, int lda, double* b, int ldb, double beta, double* c, int ldc)
{ dsyr2k_(&uplo, &trans, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); }

void rtrmm(char side, char uplo, char transa, char diag, int m, int n, double alpha, double* a, int lda, double* b, int ldb)
{ dtrmm_(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb); }

void rtrsm(char side, char uplo, char transa, char diag, int m, int n, double alpha, double* a, int lda, double* b, int ldb)
{ dtrsm_(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb); }

// only overwrites uplo portion of c!!!
void rsyrk(char uplo, char trans, int n, int k, double alpha, double* a, int lda, double beta, double* c, int ldc)
{ dsyrk_(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc); }

//------------------------------------------------------------------------------

void rgesv(int n, int nrhs, double* a, int lda, int* ipiv, double* b, int ldb, int& info)
{ dgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info); }

void rposv(char uplo, int n, int nrhs, double* a, int lda, double* b, int ldb, int& info)
{ dposv_(&uplo, &n, &nrhs, a, &lda, b, &ldb, &info); }

void rpotrf(char uplo, int n, double* a, int lda, int& info)
{ dpotrf_(&uplo, &n, a, &lda, &info); }

////////////////////////////////////////////////////////////////////////////////
// FLOAT

#ifndef DISABLE_SINGLE

// It seems that there is a problem with using BLAS Level 1 with float.
void raxpy(int n, float da, float* dx, int incx, float* dy, int incy)
{ 
  saxpy_(&n, &da, dx, &incx, dy, &incy); 
}

// float rdot(int n, float* dx, int incx, float* dy, int incy)
// { 
//   throw std::runtime_error("Error: sdot_ always returns 0.\n");
//   // return sdot_(&n, dx, &incx, dy, &incy); 
// }

//------------------------------------------------------------------------------

void rgemm(char transa, char transb, int m, int n, int k, float alpha, float* a, int lda, float* b, int ldb, float beta, float* c, int ldc)
{ sgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); }

void rsyr2k(char uplo, char trans, int n, int k, float alpha, float* a, int lda, float* b, int ldb, float beta, float* c, int ldc)
{ ssyr2k_(&uplo, &trans, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); }

void rtrmm(char side, char uplo, char transa, char diag, int m, int n, float alpha, float* a, int lda, float* b, int ldb)
{ strmm_(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb); }

void rtrsm(char side, char uplo, char transa, char diag, int m, int n, float alpha, float* a, int lda, float* b, int ldb)
{ strsm_(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb); }

void rsyrk(char uplo, char trans, int n, int k, float alpha, float* a, int lda, float beta, float* c, int ldc)
{ return ssyrk_(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc); }

//------------------------------------------------------------------------------

void rgesv(int n, int nrhs, float* a, int lda, int* ipiv, float* b, int ldb, int& info)
{ sgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info); }

void rposv(char uplo, int n, int nrhs, float* a, int lda, float* b, int ldb, int& info)
{ sposv_(&uplo, &n, &nrhs, a, &lda, b, &ldb, &info); }

void rpotrf(char uplo, int n, float* a, int lda, int& info)
{ spotrf_(&uplo, &n, a, &lda, &info); }

#endif
