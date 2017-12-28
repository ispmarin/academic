#include "functions.h"

void find_nps(int np, int &nprow, int & npcol) {

	if (np < 4) {
		cout << "Not enough processes to create grid" << endl;
		MPI::Finalize();
	}
	int min_nprow = 100000;
	int min_npcol = 100000;

	nprow = np;
	npcol = np;

	while (1) {

		npcol--;
		if (np % 2 == 0) {
			if (npcol == 1) {
				nprow--;
				npcol = nprow;
			}
		} else {
			if (npcol == 0) {
				nprow--;
				npcol = nprow;
			}

		}

		if (nprow * npcol == np) {
			min_npcol = npcol;
			if (nprow < min_nprow)
				min_nprow = nprow;
		}

		if (nprow == 1)
			break;

	}

	nprow = min_nprow;
	npcol = min_npcol;

}

void distribute_1D(int global_N, int locR, double * B, double *& B_loc,
		int block, int nprow, int npcol, int myrow, int mycol) {

	int il = 0;
	for (int i = 1; i < global_N + 1; i++) {
		int pi = p_of_i(i,block,nprow);
		int li = l_of_i(i,block,nprow);
		int xi = x_of_i(i,block,nprow);

		il = li * block + xi;
		if ((pi == myrow) && (mycol == 0)) {

			//cout<<"il "<<il-1<<" locR "<<locR<<" myrow "<<myrow<<" mycol "<<mycol<<" B "<<B[i-1]<<endl;
			B_loc[il - 1] = B[i - 1];
		} else if ((pi == myrow) && (mycol != 0)) {
			B_loc[il - 1] = 0;
		}
	}

}

void distribute_1D_loc(int global_N, int locR, aem_container &k_vec,
		aem_container &u_vec, aux_functions aquifer, double *& B_loc,
		int block, int nprow, int npcol, int myrow, int mycol) {

	int li = 0;
	int xi = 0;
	int k = 0;
	for (int i = 0; i < locR; i++) {

		li = i / block;
		xi = i % block;

		k = (li * nprow + myrow) * block + xi;
		B_loc[i] = aquifer.get_k_vector_coeff(k, k_vec, u_vec, aquifer);

	}

}

void redistribute_1D(int locR, double * B_loc, double *& B, int block,
		int nprow, int myrow) {

	int li = 0;
	int xi = 0;

	cout << "locR " << locR << endl;
	for (int i = 0; i < locR; i++) {
		li = i / block;
		xi = i % block;

		int k = (li * nprow + myrow) * block + xi;
		//cout<<"i_loc "<<i<<" li "<<li<<" k "<<k<<" xi "<<xi<<" myrow "<<myrow<<endl;
		B[k] = B_loc[i];

	}

}

int check_redistribute(int global_N, double epsilon, double *B, double * B_orig) {
	int info = 0;
	for (int i = 0; i < global_N; i++) {
		if (B[i] - B_orig[i] >= epsilon) {
			cout << "ERROR! on element " << i << " B " << B[i] << " B_orig "
					<< B_orig[i] << endl;
			info = 1;
			break;
		} else {
			info = 0;
		}

	}
	return info;
}

void copy_vector(int global_N, double *B, double *&Bcpy) {

	for (int i = 0; i < global_N; i++) {
		Bcpy[i] = B[i];
	}
}

void distribute_2D(int global_N, int locR, int locC, double * A,
		double *& A_loc, int block, int nprow, int npcol, int myrow, int mycol) {

	for (int i = 1; i < global_N + 1; i++) {
		for (int j = 1; j < global_N + 1; j++) {

			int pi = p_of_i(i,block,nprow);
			int li = l_of_i(i,block,nprow);
			int xi = x_of_i(i,block,nprow);

			int pj = p_of_i(j,block,npcol);
			int lj = l_of_i(j,block,npcol);
			int xj = x_of_i(j,block,npcol);

			if ((pi == myrow) && (pj == mycol)) {
				int il = li * block + xi;
				int jl = lj * block + xj;
				A_loc[getIndex(il - 1, jl - 1, locC)] = A[getIndex(i - 1,
						j - 1, global_N)];

			}

		}

	}

}

void build_2D_loc(int global_N, int locR, int locC, aem_container& vec,
		aux_functions aquifer, double *& A_loc, int block, int nprow,
		int npcol, int myrow, int mycol) {

	int li = 0;
	int xi = 0;
	int lj = 0;
	int xj = 0;

	for (int i = 0; i < locR; i++) {
		for (int j = 0; j < locC; j++) {

			li = i / block;
			xi = i % block;
			lj = j / block;
			xj = j % block;

			int k = (li * nprow + myrow) * block + xi;
			int l = (lj * npcol + mycol) * block + xj;

			A_loc[getIndex(i, j, locC)] = aquifer.get_matrix_coeff(k, l, vec, aquifer);

		}

	}

}

void distribute_2D_loc(int global_N, int locR, int locC, double * A,
		double *& A_loc, int block, int nprow, int npcol, int myrow, int mycol) {

	int li = 0;
	int xi = 0;
	int lj = 0;
	int xj = 0;

	for (int i = 0; i < locR; i++) {
		for (int j = 0; j < locC; j++) {

			li = i / block;
			xi = i % block;
			lj = j / block;
			xj = j % block;

			int k = (li * nprow + myrow) * block + xi;
			int l = (lj * npcol + mycol) * block + xj;

			A_loc[getIndex(i, j, locC)] = A[getIndex(k, l, global_N)];

		}

	}

}

void transpose(int locR, int locC, double *A, double *&A_transp) {

	if (locR != locC) {
		// cout<<"Remember of changing R and C when rows and colums are different"<<endl;
	}
	for (int i = 0; i < locR; i++) {
		for (int j = 0; j < locC; j++) {
			A_transp[getIndex(j, i, locR)] = A[getIndex(i, j, locC)];
		}
	}
}

void print_AB(int size_row, int size_col, double *A, int choice_rank) {

	int rank = MPI::COMM_WORLD.Get_rank();
	if (choice_rank == -1) {
		cout << " +++++++ rank: " << rank << " locR " << size_row << " locC "
				<< size_col << endl;
		for (int i = 0; i < size_row; i++) {
			for (int j = 0; j < size_col; j++) {
				cout << A[getIndex(i, j, size_col)] << "    ";
			}
			cout << endl;
		}
		cout << "+++++++++++++" << endl;
	} else {
		if (rank == choice_rank) {
			cout << "rank: " << rank << " locR " << size_row << " locC "
					<< size_col << endl;
			for (int i = 0; i < size_row; i++) {
				for (int j = 0; j < size_col; j++) {
					cout << A[getIndex(i, j, size_col)] << "    ";
				}
				cout << endl;
			}
			cout << "+++++++++++++" << endl;

		}
	}
}

double frobeninus_check(int global_N, int nrhs, int * descA, int *descB,
		double* A_loc, double* A_loc_cpy, double * B_loc, double *B_loc_cpy,
		double *& R_loc, int ictxt) {
	int ione = 1;
	double mone = (-1.0e0), pone = (1.0e0);
	double * work = new double [global_N*2];

	pdlacpy_("All", &global_N, &nrhs, B_loc_cpy, &ione, &ione, descB, R_loc,
			&ione, &ione, descB);

	double eps = pdlamch_(&ictxt, "Epsilon");
	pdgemm_("N", "N", &global_N, &nrhs, &global_N, &pone, A_loc_cpy, &ione,
			&ione, descA, B_loc, &ione, &ione, descB, &mone, R_loc, &ione,
			&ione, descB);
	double AnormF = pdlange_("F", &global_N, &global_N, A_loc, &ione, &ione,
			descA, work);
	//double BnormF = pdlange_("F", &global_N, &nrhs, B_loc_cpy, &ione, &ione,
			//descB, work);
	double XnormF = pdlange_("F", &global_N, &nrhs, B_loc, &ione, &ione, descB,
			work);
	double RnormF = pdlange_("F", &global_N, &nrhs, R_loc, &ione, &ione, descB,
			work);
	double residF = RnormF / (AnormF * XnormF * eps * ((double) global_N));

	return residF;

}

void set_rand_AB(int global_N, int seed, double *&A, double *&B) {
	srand(seed);

	for (int i = 0; i < global_N; i++) {
		for (int j = 0; j < global_N; j++) {
			A[getIndex(i, j, global_N)] = ((double) rand())
					/ ((double) RAND_MAX) - 0.5;
		}
		B[i] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
	}

}

void matrix_distribution(int N, int nprow, int npcol, aem_container& the_list,
		aem_container & the_u_list, aux_functions aquifer_reference,
		double*&A_loc, int block, int locR, int locC, int myrow, int mycol) {

	int rank = MPI::COMM_WORLD.Get_rank();
	int izero = 0;

	if (rank == 0) {

		int locR_l = 0;
		int locC_l = 0;
		for (int i = 0; i < nprow; i++) {
			for (int j = 0; j < npcol; j++) {

				int send_to_rank = getIndex(i, j, npcol);

				if (send_to_rank != 0) {
					locR_l = numroc_(&N, &block, &i, &izero, &nprow);
					locC_l = numroc_(&N, &block, &j, &izero, &npcol);
					double * A_loc_2 = new double[locR_l * locC_l];

					build_2D_loc(N, locR_l, locC_l, the_u_list,
							aquifer_reference, A_loc_2, block, nprow, npcol, i,
							j);

					MPI::COMM_WORLD.Send(A_loc_2, locR_l * locC_l, MPI::DOUBLE,
							send_to_rank, A_loc_tag);

					delete[] A_loc_2;
				}
			}
			//rank 0
			build_2D_loc(N, locR, locC, the_u_list, aquifer_reference, A_loc,
					block, nprow, npcol, myrow, mycol);
			print_AB(locR, locC, A_loc, 0);

		}

	} else {

		MPI::COMM_WORLD.Recv(A_loc, locR * locC, MPI::DOUBLE, 0, A_loc_tag);
		print_AB(locR, locC, A_loc, -1);
	}
}

void vector_distribution(int N, int nprow, int npcol, int locR,
		aem_container the_list, aem_container the_u_list,
		aux_functions aquifer_reference, int block, double *&B_loc, int myrow,
		int mycol) {

	int rank = MPI::COMM_WORLD.Get_rank();
	int izero = 0;

	if (rank == 0) {
		cout << "Starting vector distribution" << endl;

		for (int i = 0; i < nprow; i++) {

			int send_to_proc = i * npcol;
			int locR_l = numroc_(&N, &block, &i, &izero, &nprow);
			double * B_loc_2 = new double[locR_l];

			distribute_1D_loc(N, locR_l, the_list, the_u_list,
					aquifer_reference, B_loc_2, block, nprow, npcol, i, 0);

			MPI::COMM_WORLD.Send(B_loc_2, locR_l, MPI::DOUBLE, send_to_proc,
					B_loc_tag_2);
			delete [] B_loc_2;
		}
		distribute_1D_loc(N, locR, the_list, the_u_list, aquifer_reference,
				B_loc, block, nprow, npcol, myrow, mycol);
	} else {
		if (mycol == 0) {
			MPI::COMM_WORLD.Recv(B_loc, locR, MPI::DOUBLE, 0, B_loc_tag_2);
		}
	}
}

void set_results(aux_functions &aquifer, aem_container &uk_list, aem_container &k_list,
		double * unknown_vector, long int size) {

	int i = 0;
	aem_container::aem_element_iterator it;

	for (it = uk_list.begin(); it != uk_list.end(); ++it) {
		it->set_discharge(unknown_vector[i]);
		i++;
	}

	aquifer.set_constant( unknown_vector[size - 1]);

	uk_list.add_vectors(k_list.return_aem_list());

}
