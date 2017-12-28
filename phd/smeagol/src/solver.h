#ifndef SOLVER_H_
#define SOLVER_H_
#include "definitions.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "solver.h"
//#include <gsl/gsl_linalg.h>
#include <mkl/mkl.h>
#include <mkl/mkl_blas.h>
#include <mkl/mkl_lapack.h>



extern "C" void dgesvx(char* fact, char* trans, MKL_INT* n, MKL_INT* nrhs,
		double* a, MKL_INT* lda, double* af, MKL_INT* ldaf, MKL_INT* ipiv,
		char* equed, double* r, double* c, double* b, MKL_INT* ldb, double* x,
		MKL_INT* ldx, double* rcond, double* ferr, double* berr, double* work,
		MKL_INT* iwork, MKL_INT* info);

extern "C" void dgesv_(MKL_INT* n, MKL_INT* nrhs, double* a, MKL_INT* lda,
		MKL_INT* ipiv, double* b, MKL_INT* ldb, MKL_INT* info);

//extern "C" int dgesv_(int *n, int *nrhs, double *a, int
//	*lda, int *ipiv, double *b, int *ldb, int *info);

/*
 *
 *  N       (input) INTEGER
 *          The number of linear equations, i.e., the order of the
 *          matrix A.  N >= 0.
 *
 *  NRHS    (input) INTEGER
 *          The number of right hand sides, i.e., the number of columns
 *          of the matrix B.  NRHS >= 0.
 *
 *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
 *          On entry, the N-by-N coefficient matrix A.
 *          On exit, the factors L and U from the factorization
 *          A = P*L*U; the unit diagonal elements of L are not stored.
 *
 *  LDA     (input) INTEGER
 *          The leading dimension of the array A.  LDA >= max(1,N).
 *
 *  IPIV    (output) INTEGER array, dimension (N)
 *          The pivot indices that define the permutation matrix P;
 *          row i of the matrix was interchanged with row IPIV(i).
 *
 *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
 *          On entry, the N-by-NRHS matrix of right hand side matrix B.
 *          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
 *
 *  LDB     (input) INTEGER
 *          The leading dimension of the array B.  LDB >= max(1,N).
 *
 *  INFO    (output) INTEGER
 *          = 0:  successful exit
 *          < 0:  if INFO = -i, the i-th argument had an illegal value
 *          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
 *                has been completed, but the factor U is exactly
 *                singular, so the solution could not be computed. */

int lapack_solver(double * &coeff_matrix, double * &know_vector,
		double*& unknow_vector, long int matrix_dimension);
int gsl_solver(double * coeff_matrix, double * know_vector,
		double*& unknow_vector, long int matrix_dimension);
void mkl_solver(long int matrix_size, double *& matrix, double *&know_vector,
		double * &solution_vecto);

#endif /*SOLVER_H_*/
