#include "functions.h"

int main(int argc, char ** argv) {

	MPI::Init();

	int rank = MPI::COMM_WORLD.Get_rank();
	int nprocs = MPI::COMM_WORLD.Get_size();

	bool lapack_comp;
	int myrow = 0;
	int mycol = 0;
	int ictxt;
	int nprow = 0, npcol = 0;
	int locR = 0, locC = 0;
	int izero = 0;
	int ione = 1;
	int nrhs = 1;
	int info = 0;
	int locR_vector = 0;
	int N = 0, block = 0;
	double MPI_mat_dist1 = 0, MPI_mat_dist2 = 0;
	double MPI_vec_dist1 = 0, MPI_vec_dist2 = 0;
	double MPI_calc_dist1 = 0, MPI_calc_dist2 = 0;
	double MPI_redist_dist1 = 0, MPI_redist_dist2 = 0;
	double MPI_plot_1 = 0, MPI_plot_2 = 0;
	int root = 0;
	int one = 1;
	block = 2;
	lapack_comp = 0;

	find_nps(nprocs, nprow, npcol);

	Cblacs_get(-1, 0, &ictxt);
	Cblacs_gridinit(&ictxt, "Row", nprow, npcol);
	Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);

	const char * aq_file = "aquifer.dat";
	const char * ref_file = "reference.dat";
	aem_container the_list;
	aem_container the_u_list;
	vector<aux_polygon> my_polygon;
	aux_functions aquifer_reference = open_reference_file(aq_file, ref_file);

	//setting the vector
	if (rank == 0) {



		open_run_file(argv[1], the_list, the_u_list, my_polygon,
				aquifer_reference);

		cout << "Setting the lists" << endl;
		cout << "know_vector size: " << the_list.size() << endl;
		cout << "unknow_vector size: " << the_u_list.size() << endl;

		N = the_u_list.size() +1;//+1 is for the reference point

	}

	MPI::COMM_WORLD.Bcast(&N, one, MPI::INT, root);

	locR = numroc_(&N, &block, &myrow, &izero, &nprow);
	locR_vector = numroc_(&N, &block, &myrow, &izero, &nprow);
	locC = numroc_(&N, &block, &mycol, &izero, &npcol);

	//double * A = new double[N * N];
	//	double * Alapack = new double[N * N];
	//	double * Alapack_transp = new double[N * N];
	double * B = new double[N];
	//	double * Acpy = new double[N * N];
	//	double * Bcpy = new double[N];
	//	double * Blapack = new double[N];

	double * B_loc = new double[locR_vector];
	double * R_loc = new double[locR_vector];
	double * B_loc_cpy = new double[locR_vector];
	double * A_loc = new double[locR * locC];
	double * A_loc_transp = new double[locR * locC];
	double * A_loc_transp_cpy = new double[locR * locC];

	if (lapack_comp) {
		//	memcpy(Alapack, A, N * N * sizeof(double));
		//	memcpy(Blapack, B, N * sizeof(double));
	}

	MPI_mat_dist1 = MPI_Wtime();
	cout << "Matrix starting data: N " << N << " nprow " << nprow << " npcol "
			<< npcol << " block " << block << " locR " << locR << " locC "
			<< locC << " myrow " << myrow << " mycol " << mycol << endl;

	matrix_distribution(N, nprow, npcol, the_list, the_u_list,
			aquifer_reference, A_loc, block, locR, locC, myrow, mycol);
	MPI_mat_dist2 = MPI_Wtime();

	MPI_vec_dist1 = MPI_Wtime();
	vector_distribution(N, nprow, npcol, locR, the_list, the_u_list,
			aquifer_reference, block, B_loc, myrow, mycol);
	print_AB(locR,1,B_loc,2);
	MPI_vec_dist2 = MPI_Wtime();

	//	if (lapack_comp)
	//		distribute_1D(N, locR, B, B_loc_cpy, block, nprow, npcol, myrow, mycol);

	int * ipiv = new int[locC * locR * block];
	int *descA = new int[9];
	int *descB = new int[9];
	int itemp = max(1, locR);

	descinit_(descA, &N, &N, &block, &block, &izero, &izero, &ictxt, &itemp,
			&info);
	descinit_(descB, &N, &nrhs, &block, &block, &izero, &izero, &ictxt, &itemp,
			&info);
	transpose(locR, locC, A_loc, A_loc_transp);

	//	if (lapack_comp)
	//		transpose(locR, locC, A_loc, A_loc_transp_cpy);

	int temp = locR;
	locR = locC;
	locC = temp;

	if (rank == 0)
		cout << "Starting PDGESV" << endl;
	MPI_calc_dist1 = MPI_Wtime();
	pdgesv_(&N, &nrhs, A_loc_transp, &ione, &ione, descA, ipiv, B_loc, &ione,
			&ione, descB, &info);
	MPI_calc_dist2 = MPI_Wtime();

	//setting lapack comparison
	if (lapack_comp) {
		if (rank == 0) {
			cout << "Staring DGESV comparison" << endl;

			//		transpose(N, N, Alapack, Alapack_transp);

			//	dgesv_(&N, &nrhs, Alapack_transp, &N, ipiv, Blapack, &N, &info);

		}
	}
	//cout<<"Calculating Frobenius norm"<<endl;
	double norm = frobeninus_check(N, nrhs, descA, descB, A_loc_transp,
			A_loc_transp_cpy, B_loc, B_loc_cpy, R_loc, ictxt);

	MPI_redist_dist1 = MPI_Wtime();
	if (rank == 0) {
		cout << "Redistributing solution" << endl;
		int locR_recv = 0;
		double * B_reduced = new double[N];

		for (int i = 1; i < nprow; i++) {

			int sending_proc = i * npcol; //receive only from the mycol == 0


			MPI::COMM_WORLD.Recv(&locR_recv, 1, MPI::INT, sending_proc,
					locR_tag);
			double * B_loc_recv = new double[locR_recv];

			MPI::COMM_WORLD.Recv(&B_loc_recv[0], locR_recv, MPI::DOUBLE,
					sending_proc, B_loc_tag);

			redistribute_1D(locR_recv, B_loc_recv, B_reduced, block, nprow, i);

			delete [] B_loc_recv;

		}

		redistribute_1D(locR_vector, B_loc, B_reduced, block, nprow, myrow);
		if (lapack_comp) {
			//if (check_redistribute(N, 0.01, Blapack, B_reduced) == 0)
			//cout << "Lapack/Scalapack/Redistribution check OK" << endl;
		}

		set_results(aquifer_reference, the_u_list, the_list,
				B_reduced, N);

		print_AB(N, 1, B_reduced, 0);

		delete [] B_reduced;
	} else {
		if (mycol == 0) {
			MPI::COMM_WORLD.Send(&locR_vector, 1, MPI::INT, 0, locR_tag);
			MPI::COMM_WORLD.Send(&B_loc[0], locR_vector, MPI::DOUBLE, 0,
					B_loc_tag);
		}
	}

	MPI_redist_dist2 = MPI_Wtime();

	//if(rank == 0) {

	if (norm < 10.0) {
		cout << "Frobenius norm is acceptable: " << norm << " rank " << rank
				<< endl;
	} else {
		cout << "PROBLEM: Frobenius norm is NOT acceptable: " << norm << endl;
	}

	//}

	if(rank ==0 ){
		MPI_plot_1 = MPI_Wtime();
		cout <<"Plotting results"<<endl;
		const char * plot_file = "plot.dat";
		export_to_file plot_iso = open_plot_file(plot_file);
		plot_iso.plot_surface_cond(aquifer_reference, the_u_list, my_polygon);
		MPI_plot_2 = MPI_Wtime();
		cout << "Times: Matriz distribution: " << MPI_mat_dist2 - MPI_mat_dist1
				<< " Vector distribution: " << MPI_vec_dist2 - MPI_vec_dist1
				<< " Computation: " << MPI_calc_dist2 - MPI_calc_dist1
				<< " Result reduction: " << MPI_redist_dist2 - MPI_redist_dist1
				<< " Plotting:  " << MPI_plot_2 -  MPI_plot_1
				<< endl;
	}

	//global delete
	//	delete[] A;
	//	delete[] Alapack;
	//	delete[] Alapack_transp;
	delete[] B;
	//	delete[] Acpy;
	//	delete[] Bcpy;
	//	delete[] Blapack;

	//local delete
	//delete [] B_loc ;
	delete[] R_loc;
	delete[] B_loc_cpy;
	//delete [] A_loc ;
	delete[] A_loc_transp;
	delete[] A_loc_transp_cpy;
	delete[] ipiv;
	//delete work;
	delete[] descA;
	delete[] descB;

	Cblacs_gridexit(ictxt);
	MPI::Finalize();
	return 0;
}

