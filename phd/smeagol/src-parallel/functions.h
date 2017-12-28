#include <iostream>
#include <vector>
#include <openmpi/mpi.h>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <mkl/mkl.h>
#include <mkl/mkl_scalapack.h>
#include <mkl/mkl_cblas.h>
#include "src/elements.h"
#include "src/smeagol.h"
#include "src/sequential_interface.h"
#include "src/aux_functions.h"
#include "src/definitions.h"

using namespace std;
#define p_of_i(i,bs,p) ( int((i-1)/bs)%p)
#define l_of_i(i,bs,p) (int((i-1)/(p*bs)))
#define x_of_i(i,bs,p) (((i-1)%bs)+1)

const int locR_tag = 1000;
const int myrow_tag = 4001;
const int mycol_tag = 1002;
const int B_loc_tag = 1003;
const int A_loc_tag = 1004;
const int B_loc_tag_2 = 1005;

void find_nps(int np, int &nprow, int & npcol);
//int getIndex(int row, int col,int NCOLS) {return row*NCOLS+col;}
void set_AB(double *&A, double *&B);
void set_rand_AB(int global_N, int seed, double *&A, double *&B);
void print_AB(int size, int size_y, double *A, int rank);
void distribute_2D(int global_N, int locR, int locC, double * A,
		double *& A_loc, int block, int nprow, int npcol, int myrow, int mycol);
void distribute_2D_loc(int global_N, int locR, int locC, double * A,
		double *& A_loc, int block, int nprow, int npcol, int myrow, int mycol);
void distribute_1D(int global_N, int locR, double * B, double *& B_loc,
		int block, int nprow, int npcol, int myrow, int mycol);
void distribute_1D_loc(int global_N, int locR, aem_container& k_vec,
		aem_container& u_vec, aux_functions aquifer, double *& B_loc,
		int block, int nprow, int npcol, int myrow, int mycol);
//void redistribute_1D(int global_N, double * B_loc, double *& B, int block, int nprow, int npcol, int myrow, int mycol);
void redistribute_1D(int locR, double * B_loc, double *& B, int block,
		int nprow, int myrow);
void copy_vector(int global_N, double *B, double *&Bcpy);
void transpose(int locR, int locC, double *A, double *&A_transp);
double frobeninus_check(int global_N, int nrhs, int * descA, int *descB,
		double* A_loc, double *A_loc_cpy, double * B_loc, double *B_loc_cpy,
		double *& R_loc, int ictxt);
int
		check_redistribute(int global_N, double epsilon, double *B,
				double * B_orig);

void build_2D_loc(int global_N, int locR, int locC, aem_container& vec,
		aux_functions aquifer, double *& A_loc, int block, int nprow,
		int npcol, int myrow, int mycol);
void set_aem_elements();

void matrix_distribution(int N, int nprow, int npcol, aem_container& the_list,
	aem_container& the_u_list, aux_functions aquifer_reference,
	double*&A_loc,  int block,int locR, int locC, int myrow, int mycol);
void vector_distribution(int N, int nprow, int npcol,int locR, aem_container the_list,
		aem_container the_u_list, aux_functions aquifer_referecen, int block,
		double *&B_loc, int myrow, int mycol);
void set_results(aux_functions &aquifer, aem_container &uk_list, aem_container &k_list,
		double * unknown_vector, long int size);

extern "C" {
/* BLACS C interface */
void Cblacs_get(int context, int request, int* value);
int Cblacs_gridinit(int* context, char * order, int np_row, int np_col);
void Cblacs_gridinfo(int context, int* np_row, int* np_col, int* my_row,
		int* my_col);
int numroc_(int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);
void Cblacs_gridexit(int ictxt);
void pdgesv_(MKL_INT *n, MKL_INT *nrhs, double *a, MKL_INT *ia, MKL_INT *ja,
		MKL_INT *desca, MKL_INT *ipiv, double *b, MKL_INT *ib, MKL_INT *jb,
		MKL_INT *descb, MKL_INT *info);
void dgesv_(MKL_INT* n, MKL_INT* nrhs, double* a, MKL_INT* lda, MKL_INT* ipiv,
		double* b, MKL_INT* ldb, MKL_INT* info);

void descinit_(int *desc, int *m, int *n, int *mb, int *nb, int *irsrc,
		int *icsrc, int *ictxt, int *lld, int *info);
double pdlamch_(int *ictxt, char *cmach);
double pdlange_(char *norm, int *m, int *n, double *A, int *ia, int *ja,
		int *desca, double *work);

void pdlacpy_(char *uplo, int *m, int *n, double *a, int *ia, int *ja,
		int *desca, double *b, int *ib, int *jb, int *descb);

void pdgemm_(char *TRANSA, char *TRANSB, int * M, int * N, int * K,
		double * ALPHA, double * A, int * IA, int * JA, int * DESCA,
		double * B, int * IB, int * JB, int * DESCB, double * BETA, double * C,
		int * IC, int * JC, int * DESCC);
int indxg2p_(int *indxglob, int *nb, int *iproc, int *isrcproc, int *nprocs);

}

