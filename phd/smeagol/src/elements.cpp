/*
 * elements.cpp
 *
 *  Created on: 27/05/2009
 *      Author: ispmarin
 */

#include "elements.h"
using namespace std;

//++++++++++++++++++++++++++++++++++++++++++++++++++
//AEM_ELEMENT

complex<double> aem_element::coord_transform(complex<double> center, complex<
		double> next, complex<double> point) {
	complex<double> Z;

	Z = (point - 0.5 * (center + next)) / (0.5 * (next - center));
	return Z;

}

//++++++++++++++++++++++++++++++++++++++++++++++++++
//POINT

double well::get_coeff(complex<double> z) {
	double x = 0;
	double y = 0;
	double pos_x = 0;
	double pos_y = 0;
	double tolerance = 1;

	double potential = 0;
	if (abs(z - position) < tolerance) { //bug, must be less than, not equal
		x = real(position) + 1;
		y = imag(position) + 1;
	} else {
		x = real(z);
		y = imag(z);
	}

	pos_x = real(position);
	pos_y = imag(position);
	potential = (1.0 / (4.0 * Pi)) * log((pow((x - pos_x), 2) + pow(
			(y - pos_y), 2)) / radius);

	return potential;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++
//LINE_SINK

double line_sink::get_coeff(complex<double> z) {
	double coeff = 0;
	complex<double> num_dois(2, 0);
	complex<double> num_um(1, 0);
	complex<double> diff(0.01, 0.01);
	complex<double> coeff_complex;

	double L = sqrt(norm(position_end - position));

	if ((position_end == z) || (position == z)) {
		z += diff;//to remove the singularity on the endpoints...
	}

	//complex<double> Z = coord_transform(position, position_end, z);

	//coeff = (L / (4*Pi)) * real( (Z+1.0)*log(Z+1.0) - (Z-1.0)*log(Z-1.0) + 2* log (.5* L ) - 2.0 );

	coeff = (L / (4 * Pi)) * real((2 * log(.5 * L) - 2.0 + 2.0 * (((z
			- position) / (position_end - position)) * log(2.0 * (z - position)
			/ (position_end - position))) - 2.0 * (((z - position_end)
			/ (position_end - position)) * log(2.0 * (z - position_end)
			/ (position_end - position)))));
	//workaround for test purposes: should't the singularities at the tips cancel out?
	if (isnan(coeff)) {
		printf("x: %lf  y: %lf\n", real(z), imag(z));
	}

	return coeff;

}


//++++++++++++++++++++++++++++++++++++++++++++++++++
//LINE DOUBLET

double line_doublet::get_coeff(complex<double> z) {

	complex<double> I(0, 1);
	double coeff = 0;
	double A = 0;

	complex<double> previous = line_doublet::previous->get_control_point();
	complex<double> next = line_doublet::next->get_control_point();

	A = 1.0 / ((ld_condutivity - external_condutivity) / external_condutivity);

	if ((z == previous) || (z == next) || (position == z)) { //I'm a vertice


		if (z == position) { //I'm in myself, matrix element - OR NOT! bug if we are not in the matrix
			coeff = (1.0 / (2.0 * Pi)) * (-1.0) * get_angle();
			return (coeff - A);

		} else if (z == previous) {
			coeff = -real((1.0 / (2.0 * Pi * I)) * (((previous - next) / (next
					- position)) * log((previous - next)
					/ (previous - position))));

		} else if (z == next) {
			coeff = -real((1.0 / (2.0 * Pi * I)) * (((next - position)
					/ (position - previous)) * log((next - position) / (next
					- previous))));
		}
		return coeff;

	} else {//everywhere else


		//coeff = (1.0 / (2.0 * Pi)) * (-1.0) * get_angle();

		coeff = real((1.0 / (2.0 * Pi * I))
		* (((z - previous) / (position - previous)) * log((z - position) / (z - previous))
				+ (( next - z ) / (next- position)) * log((z - next) / (z - position))));
//
//
//		coeff = real((1.0 / (2.0 * Pi * I))
//				* (((z - previous) / (position - previous)) * log(
//						(z - position) / (z - previous)) - ((z - next) / (next
//						- position)) * log((z - next) / (z - position))));

		return coeff;

	}

}

double line_doublet::get_final_potential(complex<double> z) {
	complex<double> I(0, 1);
	double coeff = 0;
	double final = 0;

	//cout << "next discharge " << next->get_discharge() << endl;

	complex<double> previous = line_doublet::previous->get_control_point();
	complex<double> next_c = line_doublet::next->get_control_point();

	complex<double> Z_m = coord_transform(position, next_c, z);

	if ((z == previous) || (z == next_c) || (position == z)) { //I'm a vertice

		if (z == position) {
			//			printf("get_angle: %lf\n\n",(1.0 / (2.0 * Pi)) * (-1.0) *get_angle());
			coeff = -(1.0 / (2.0 * Pi)) * get_angle();

		} else if (z == previous) {
			coeff = -real((1.0 / (2.0 * Pi * I)) * (((previous - next_c)
					/ (next_c - position)) * log((previous - next_c)
					/ (previous - position))));
			//printf("coeff previous: %lf\n\n",coeff);
		} else if (z == next_c) {
			coeff = -real((1.0 / (2.0 * Pi * I)) * (((next_c - position)
					/ (position - previous)) * log((next_c - position)
					/ (next_c - previous))));
			//printf("coeff next: %lf\n\n",coeff);
		}

		final = get_discharge() * coeff;

	} else {//everywhere else


		coeff = real((1.0 / (2.0 * Pi * I))
	* (((z - previous) / (position - previous)) * log((z - position) / (z - previous))
			+ (( next_c - z ) / (next_c- position)) * log((z - next_c) / (z - position))));
	final = get_discharge()*coeff;

	}




//	// percorre todas as line doublets aqui dentro, para saber onde estão as singularidades, e diz que elas somem na soma!
//
//	//cout << "next discharge " << next->get_discharge() << endl;
//
//	complex<double> previous = line_doublet::previous->get_control_point();
//	complex<double> next_c = line_doublet::next->get_control_point();
//
//	complex<double> Z_m = coord_transform(position, next_c, z);
//
//	if ( (z == next_c) || (position == z)) { //I'm a vertice
//		Z_m +=0.01;
//		final = real((1.0 / (2.0 * Pi)) * (get_discharge() * (-0.5
//					* (Z_m - 1.0) * log((Z_m - 1.0) / (Z_m + 1.0)) - 1.0)
//					+ next->get_discharge() * (0.5 * (Z_m + 1.0) * log((Z_m - 1.0)
//							/ (Z_m + 1.0)) + 1.0) ) );
//
//	} else {//everywhere else
//
//
//		final = real((1.0 / (2.0 * Pi)) * (get_discharge() * (-0.5
//				* (Z_m - 1.0) * log((Z_m - 1.0) / (Z_m + 1.0)) - 1.0)
//				+ next->get_discharge() * (0.5 * (Z_m + 1.0) * log((Z_m - 1.0)
//						/ (Z_m + 1.0)) + 1.0) ) );
//
//	}



	return final;
}

double line_doublet::get_angle() {

	//	double A = 1.0 / ((ld_condutivity - external_condutivity)
	//			/ external_condutivity);

	complex<double> previous = line_doublet::previous->get_control_point();
	complex<double> next = line_doublet::next->get_control_point();

	complex<double> vv = next - position;
	complex<double> uu = position - previous;

	double dot_product = real(uu) * real(vv) + imag(uu) * imag(vv);

	double mod = sqrt(real(uu) * real(uu) + imag(uu) * imag(uu)) * sqrt(
			real(vv) * real(vv) + imag(vv) * imag(vv));

	double coeff_ang_novo = (acos(dot_product / mod));

	//printf("angulo entre vetores: %lf\n", coeff_ang_novo);

	return coeff_ang_novo;

}


aem_element* new_clone(const aem_element& a) {
	return a.clone();
}

