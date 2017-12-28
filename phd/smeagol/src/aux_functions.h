//
// * aux_functions.h
// *
// *  Created on: 01/06/2009
// *      Author: ispmarin
//

#ifndef AUX_FUNCTIONS_H_
#define AUX_FUNCTIONS_H_

#include <iostream>
#include <vector>
#include <complex>
#include "elements.h"
#include "definitions.h"


//GGL
#include <boost/geometry/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>

using namespace std;
using namespace boost::geometry;

class aux_functions {
private:
	aquifer * aux_ref_aquifer;
	reference * reference_data;
	double constant;

protected:

public:

	aux_functions(double conductivity, double base, double height, complex<
			double> reference_point, double head);
	aux_functions(const aux_functions & aux);
	aux_functions & operator=(const aux_functions & aux_thing);
	aux_functions(){}
	~aux_functions();
	aquifer *get_ref_aquifer() const;
	reference * get_reference() const;
	double get_aquifer_conductivity() const;
	void set_k_vector(aem_container & k_list, aem_container &uk_list,
			double *& know_vector);
	void set_matrix(aem_container & uk_list, double *& matrix);
	double head_to_pot(double head);
	double pot_to_head(double potential, double condutivity, int is_ld);
	void
	set_results(aem_container &uk_list, aem_container & k_list,
			double * unknown_vector, long int size);

	double get_constant(){
		return constant;
	}
	void set_constant(double setconstant) {constant = setconstant;}

	double get_matrix_coeff(int i,int j,aem_container list, aux_functions aquifer);
	double get_k_vector_coeff(int i,  aem_container k_list,
			aem_container u_list, aux_functions aquifer);


};

class aux_polygon {
private:
	double h_conductivity;
	polygon polygon;
public:
	void set_polygon(polygon);
	polygon get_polygon();
	void set_conductivity(double h);
	double get_conductivity();
};




#endif // AUX_FUNCTIONS_H_
