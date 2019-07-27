#include "sequential_interface.h"

numeric_interface::numeric_interface(long int set_size) {
	size = set_size;
	know_vector = new double[size];
	unknow_vector = new double[size];
	matrix = new double[size * size];
}

numeric_interface::~numeric_interface() {

	delete[] know_vector;
	delete[] unknow_vector;
	delete[] matrix;

}


void numeric_interface::set_k_vector(aux_functions aquifer, aem_container & k_list,
		aem_container & uk_list, double *& k_vector) {

	long int matrix_size = uk_list.size() + 1;
	aem_container::aem_element_iterator it_k;
	aem_container::aem_element_iterator it_uk;

	int i = 0;

	for (it_uk = uk_list.begin(); it_uk != uk_list.end(); ++it_uk) {
		k_vector[i] = 0;//overcautious
		for (it_k = k_list.begin(); it_k != k_list.end(); ++it_k) {
			k_vector[i] -= it_k->get_coeff(it_uk->get_control_point());
		}
		k_vector[i] += aquifer.head_to_pot(it_uk->get_head());
		i++;
	}
	k_vector[matrix_size - 1] = aquifer.head_to_pot(aquifer.get_reference()->get_head());
	//last element of k_vector is reference
	for (i = 0; i < matrix_size; i++) {
		cout << "k  " << k_vector[i] << "    " << i << endl;
	}

}

void numeric_interface::set_matrix(aux_functions aquifer, aem_container & uk_list, double * &matrix) {

	aem_container::aem_element_iterator it_list = uk_list.begin();
	long int matrix_size = uk_list.size() + 1;


	for (int i = 0; i < matrix_size ; i++) {
		for (int j = 0; j < matrix_size ; j++) {
			matrix[getIndex(i, j, matrix_size)] = aquifer.get_matrix_coeff(i, j,
					uk_list,aquifer);
		}
	}

/*
	for (i = 0; i < matrix_size; i++) {
		mat(matrix,matrix_size,i,(matrix_size-1)) = 1;
	}
	//last line
	j = 0;
	int last_line = matrix_size-1;
	for (it_list = uk_list.begin(); it_list != uk_list.end(); ++it_list) {
		matrix[getIndex(last_line,j,matrix_size)] = it_list->get_coeff(
				aquifer.get_reference()->get_position());
		j++;
	}
*/

}

void numeric_interface::set_results(aux_functions &aquifer, aem_container &uk_list, aem_container &k_list,
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


void numeric_interface::solve_system(aux_functions & aquifer_reference,
		aem_container &the_list, aem_container & the_u_list) {

	set_k_vector(aquifer_reference,the_list, the_u_list, know_vector);
	set_matrix(aquifer_reference,the_u_list, matrix);

	//error = gsl_solver(matrix,know_vector,unknow_vector, size);
	mkl_solver(size, matrix, know_vector, unknow_vector);


	set_results(aquifer_reference,the_u_list, the_list, unknow_vector, size);

}

