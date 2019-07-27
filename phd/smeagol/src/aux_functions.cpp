//
//* aux_functions.cpp
// *
// *  Created on: 01/06/2009
// *      Author: ispmarin
//

#include "aux_functions.h"
using namespace std;

aux_functions::aux_functions(double conductivity, double base, double height,
		complex<double> reference_point, double head) {
	aux_ref_aquifer = new aquifer(height, base, conductivity);
	reference_data = new reference(reference_point, head);
	constant = 0;
	/*printf(
	 "Aux Functions allocated correctly\n Checking values:\n cond: %lf\n base: %lf\n height: %lf\n ref_x: %lf\t ref_y: %lf\t ref_head: %lf\n",
	 aux_ref_aquifer->get_conduvitity(), aux_ref_aquifer->get_base(),
	 aux_ref_aquifer->get_height(),
	 real(reference_data->get_position()), imag(
	 reference_data->get_position()), reference_data->get_head());*/
}

aux_functions::aux_functions(const aux_functions & aux) {
	constant = 0;
	aux_ref_aquifer = new aquifer(*aux.get_ref_aquifer());
	reference_data = new reference(*aux.get_reference());

}

aux_functions & aux_functions::operator=(const aux_functions & aux) {
	if (this == &aux)
		return *this;

	delete aux_ref_aquifer;
	delete reference_data;
	constant = 0;

	aux_ref_aquifer = new aquifer(*aux.get_ref_aquifer());
	reference_data = new reference(*aux.get_reference());
	return *this;

}

aux_functions::~aux_functions() {

	delete aux_ref_aquifer;
	delete reference_data;
}

aquifer * aux_functions::get_ref_aquifer() const {
	return aux_ref_aquifer;
}

double aux_functions::get_aquifer_conductivity() const {
	return aux_ref_aquifer->get_conduvitity();
}

reference * aux_functions::get_reference() const {
	return reference_data;
}

double aux_functions::head_to_pot(double head) {

	double potential = 0;

	double base = aux_ref_aquifer->get_base();
	double height = aux_ref_aquifer->get_height();
	double h_condutivity = aux_ref_aquifer->get_conduvitity();

	if (head >= (height - base)) {
		potential = h_condutivity * (height - base) * head - 0.5
				* h_condutivity * pow((height - base), 2);
	} else {
		potential = 0.5 * h_condutivity * pow(head, 2);
	}
	return potential;
}

double aux_functions::pot_to_head(double potential, double condutivity,
		int is_ld) {

	double return_head = 0;
	double height = aux_ref_aquifer->get_height();
	double base = aux_ref_aquifer->get_base();
	double h_condutivity = 0;

	if (is_ld) {
		h_condutivity = condutivity;
	} else {
		h_condutivity = aux_ref_aquifer->get_conduvitity();
	}

	double pot_value = 0.5 * h_condutivity * pow((height - base), 2);

	if(potential < 0) {
		return aux_ref_aquifer->get_base();
	}

	if (potential >= pot_value) {
		return_head = (potential + pot_value) / (h_condutivity
				* (height - base));
	} else {
		return_head = sqrt(2 * potential / h_condutivity);
	}

	return return_head;
}

double aux_functions::get_matrix_coeff(int i, int j, aem_container list,
		aux_functions aquifer) {

	int size = list.size();

	if ((i < size) && (j < size)) {

		aem_container::aem_element_iterator it_pos = list.begin() + i;
		aem_container::aem_element_iterator it_elem = list.begin() + j;

		return it_elem->get_coeff(it_pos->get_control_point());

	} else if ((j == size) && (i < size)) {

		aem_container::aem_element_iterator it_elem = list.begin() + i;
		return it_elem->get_coeff(aquifer.get_reference()->get_position());

	} else if (i == size) {

		return 1;

	} else {

		cout << "aux size: " << size << "  i  " << i << "  j  " << j << endl;
		cout << "invalid index! " << endl;
		return -1;

	}
}
double aux_functions::get_k_vector_coeff(int i, aem_container k_list,
		aem_container u_list, aux_functions aquifer) {

	double coeff = 0;
	int size = u_list.size();
	if (i < size) {

		aem_container::aem_element_iterator it_u = u_list.begin() + i;
		aem_container::aem_element_iterator it_k = k_list.begin();

		for (it_k = k_list.begin(); it_k != k_list.end(); ++it_k) {
			coeff -= it_k->get_coeff(it_u->get_control_point());
		}
		coeff += aquifer.head_to_pot(it_u->get_head());

		return coeff;

	} else if (i == size) {

		return aquifer.head_to_pot(aquifer.get_reference()->get_head());

	} else {
		cout << "aux size: " << size << "  i  " << i << endl;
		cout << "invalid index! " << endl;
		return -1;
	}

}



//CLASS Polygon
void aux_polygon::set_polygon(polygon_2d set_polygon) {
	polygon = set_polygon;
}

polygon_2d aux_polygon::get_polygon() {
	return polygon;
}
void aux_polygon::set_conductivity(double h) {
	h_conductivity = h;
}
double aux_polygon::get_conductivity() {
	return h_conductivity;
}

