/*
 * numeric_interface.h
 *
 *  Created on: 19/10/2009
 *      Author: ispmarin
 */

#ifndef NUMERIC_INTERFACE_H_
#define NUMERIC_INTERFACE_H_

#include "aux_functions.h"
#include "definitions.h"
#include "solver.h"

class numeric_interface {
private:
	long int size;
	double * know_vector;
	double * unknow_vector;
	double * matrix;

public:
	numeric_interface(long int size);
	~numeric_interface();

	void solve_system(aux_functions & aquifer_reference, aem_container &the_list,
			aem_container & the_u_list);
	void set_results(aux_functions &aquifer,aem_container &uk_list, aem_container &k_list,
			double * unknown_vector, long int size);
	void set_matrix(aux_functions aquifer, aem_container & uk_list, double * &matrix);
	void set_k_vector(aux_functions aquifer, aem_container & k_list,
			aem_container & uk_list, double *& k_vector);
};

#endif /* NUMERIC_INTERFACE_H_ */
