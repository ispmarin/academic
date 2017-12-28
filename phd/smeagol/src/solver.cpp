#include "solver.h"

int lapack_solver(double * &coeff_matrix, double * &know_vector,
		double*& unknow_vector, long int matrix_dimension) {

//	int INFO=0;
//	int * IPIV;
//	IPIV = new int[matrix_dimension];
//	double * matrix_transp;
//	matrix_transp = (double*) calloc(matrix_dimension * matrix_dimension,
//			sizeof(double));
//
//	//int nhrs = 1;
//
//	/*mat(coeff_matrix,matrix_dimension,0,0) = 1;
//	 mat(coeff_matrix,matrix_dimension,0,1) = 0;
//	 mat(coeff_matrix,matrix_dimension,0,2) = 0;
//
//	 mat(coeff_matrix,matrix_dimension,1,0) = 0;
//	 mat(coeff_matrix,matrix_dimension,1,1) = 1;
//	 mat(coeff_matrix,matrix_dimension,1,2) = 0;
//
//	 mat(coeff_matrix,matrix_dimension,2,0) = 0;
//	 mat(coeff_matrix,matrix_dimension,2,1) = 0;
//	 mat(coeff_matrix,matrix_dimension,2,1) = 1;*/
//
//	unknow_vector[0] = 1;
//	unknow_vector[1] = 2;
//	unknow_vector[2] = 3;
//
//	/*NEEDED!! clapack convention to matrix is different from the row order in C!*/
//	//coeff_matrix = transpose(coeff_matrix, matrix_dimension);
//	//unknow_vector = know_vector;
//
//	//	dgesv_((int*)&matrix_dimension, &nhrs, coeff_matrix,
//	//		(int*)&matrix_dimension, IPIV, unknow_vector,
//	//	(int*)&matrix_dimension, &INFO);
//
//	if (INFO != 0) {
//		printf("Linear system solver PROBLEM %i", INFO);
//		exit(-1);
//	}
//
//	for (int i = 0; i < matrix_dimension; i++) {
//		printf("%f\n", unknow_vector[i]);
//	}
//	return INFO;
	return 0;
}

int gsl_solver(double * coeff_matrix, double * know_vector,
		double*& unknow_vector, long int matrix_dimension) {

	/*mat(coeff_matrix,matrix_dimension,0,0) = 1;
	 mat(coeff_matrix,matrix_dimension,0,1) = 0;
	 mat(coeff_matrix,matrix_dimension,0,2) = 0;

	 mat(coeff_matrix,matrix_dimension,1,0) = 0;
	 mat(coeff_matrix,matrix_dimension,1,1) = 1;
	 mat(coeff_matrix,matrix_dimension,1,2) = 0;

	 mat(coeff_matrix,matrix_dimension,2,0) = 0;
	 mat(coeff_matrix,matrix_dimension,2,1) = 0;
	 mat(coeff_matrix,matrix_dimension,2,1) = 1;

	 unknow_vector[0] = 1;
	 unknow_vector[1] = 2;
	 unknow_vector[2] = 3;*/

	//int i = 0;

	/*int j = 0;

	 for (i = 0; i < matrix_dimension; i++) {
	 for (j = 0; j < matrix_dimension; j++) {
	 printf("%lf\t", mat(coeff_matrix,matrix_dimension,i,j));
	 }
	 printf("\n");
	 }
	 for (j = 0; j < matrix_dimension; j++) {
	 printf("%lf\n", know_vector[j]);
	 }
	 */
//	gsl_matrix_view m = gsl_matrix_view_array(coeff_matrix, matrix_dimension,
//			matrix_dimension);

//	gsl_vector_view b = gsl_vector_view_array(know_vector, matrix_dimension);
//
//	gsl_vector *x = gsl_vector_alloc(matrix_dimension);
//
//	int s;
//
//	gsl_permutation * p = gsl_permutation_alloc(matrix_dimension);
//
//	gsl_linalg_LU_decomp(&m.matrix, p, &s);
//
//	gsl_linalg_LU_solve(&m.matrix, p, &b.vector, x);
//
//	printf("x = \n");
//	gsl_vector_fprintf(stdout, x, "%g");
//
//	for (i = 0; i < matrix_dimension; i++) {
//		unknow_vector[i] = gsl_vector_get(x, i);
//
//	}

	return 0;
}

void mkl_solver(long int matrix_size, double *& matrix, double *&know_vector,
		double * &unknow_vector) {
	//WARNING: the matrix comes from the main program in C style notation. To transverse the matrix, it is needed to use the
	//mat_transp() and NOT mat().

	printf("Using MKL Solver\n");

	int size = (int) matrix_size;
	double * matrix_transp;

	int i = 0, j = 0;

	MKL_INT nrhs = 1;

	int lda = size;
	printf("lda: %i\n", lda);
	int ldb = size;
	int *ipiv;
	MKL_INT info;
	double * work;
	double * solution_vector;

	solution_vector = new double[size];

	work = new double[size * 4];
	ipiv = new int[size];
	matrix_transp = new double[size * size];

	for (i = 0; i < size; i++) {
		solution_vector[i] = know_vector[i];
	}

	printf("Transposing matrix\n");

	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			mat_transp(matrix_transp,size,i,j) = mat(matrix,size,i,j);
		}
	}


	printf("Matrix transposed\n");

	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			printf("%lf\t", mat(matrix_transp,size,i,j));
		}
		printf("\n");
	}
	//iwork = (MKL_INT*)malloc(sizeof(MKL_INT)*size);


	printf("Calling DGESV\n");

	dgesv(&size, &nrhs, matrix_transp, &lda, ipiv, solution_vector, &ldb, &info);

	if(info != 0) {
		printf("INFO returned %i\n",info);
	}

	for (i = 0; i < size; i++) {
		unknow_vector[i] = solution_vector[i];
		printf("x: %f\n", solution_vector[i]);
	}

	delete [] solution_vector ;
	delete [] work ;
	delete [] ipiv ;
	delete [] matrix_transp ;


	//	dgesvx(fact, trans, &n, &nrhs, matrix, &lda, af, &ldaf, ipiv, equed, r, c,
	//			know_vector, &ldb, solution_vector, &ldx, &rcond, ferr, berr, work, iwork, &info);

	//	MKL_INT iter = 0;
	//double * a = matrix;
	//double * b = know_vector;
	//	double * r; //must be allocated - matrix_size
	//	double * c; //must be allocated - matrix_size
	//	double *af; //must be allocated - matrix_size*matrix_size
	//MKL_INT ldaf = size; //must confirm
	//	char *fact = "N";//will be balanced, and the factored form is not supplied
	//	char *trans = "N";//no transpose
	//	MKL_INT n = matrix_size;	//MKL_INT ldx = size; //size of output vector
	//	char *equed = "N";
	//double * x; //must be allocated - matrix_size - output vector =  solution_vector

	//double rcond;
	//double * ferr; //must be allocated - matrix_size - output error factors
	//double * berr; //must be allocated - matrix_size - output error backward
	//r = (double*)malloc(sizeof(double)*size);
	//c = (double*)malloc(sizeof(double)*size);
	//af = (double*)malloc(sizeof(double)*size*size);
	//berr = (double*)malloc(sizeof(double)*size);
	//ferr = (double*)malloc(sizeof(double)*size);
	//MKL_INT * iwork; //must be allocated - matrix_size
	//work = (double*) malloc(sizeof(double) * size * 4);
	//	ipiv = (MKL_INT*) malloc(sizeof(MKL_INT) * size);
	//	matrix_transp = (double *) malloc(sizeof(double) * size);
}
