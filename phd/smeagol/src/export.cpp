/*
 * export.cpp
 *
 *  Created on: 04/06/2009
 *      Author: ispmarin
 */

#include "export.h"

using namespace std;

void export_to_file::plot_surface_cond(aux_functions &atua,
		aem_container &final_list, vector<aux_polygon> et_lists) {

	ofstream data_exit;
	data_exit.open("head_cond.dat");
	ofstream pot_exit;
	pot_exit.open("resultados.dat");
	ofstream data_exit_head;
	data_exit_head.open("head.dat");
//	ofstream teste;
//	teste.open("discharge.dat");


	int i = 0;
	int j = 0;
	double total;
	double pos_x = 0;
	double pos_y = 0;
	int count = 0;

	//	double max_scalar = -1;
	//	double min_scalar = 1000000;
	//	double media = 0;

	double size_x = x_max - x_min + 1;
	double size_y = y_max - y_min + 1;

	long int total_size_x = size_x / increment;
	long int total_size_y = size_y / increment;
	vector<vector<double> > matrix(total_size_y, vector<double> (total_size_x));

	aem_container::aem_element_iterator it;
	aem_container::aem_element_iterator next = final_list.begin() + 1;
	cout<< "constant in export: "<< atua.get_constant()<<endl;

	vector<aux_polygon>::iterator polygon_it;

	cout<< "Exporting head\n"<<endl;

	pos_x = x_min;
	pos_y = y_min;

	for (i = 0; i < total_size_y; i++) {
		for (j = 0; j < total_size_x; j++) {
			total = 0;

			complex<double> aval_point(pos_x, pos_y);

			for (it = final_list.begin(); it != final_list.end(); ++it) {

				//total += (*it)->get_coeff(aval_point) * (*it)->get_discharge();
				total += it->get_final_potential(aval_point);
				//printf("%i\t%i\t total: %lf\n",i,j,total);

			}

//			for (it = final_list.begin(); it != final_list.end(); ++it) {
//
//				teste << "discharge: " << it->get_discharge() << endl;
//				teste << "next discharge: " << next->get_discharge() << endl;
//				next++;
//				if (next == final_list.end()) {
//					next = final_list.begin();
//				}
//			}

			if (total + atua.get_constant() < 0) {
				if (count < 10) {
					printf(
							"WARNING! NEGATIVE POTENTIAL pot:%lf\t const:%lf %lf\t%lf\n",
							total, atua.get_constant(), pos_x, pos_y);
					count++;
				} else if (count == 10) {
					cout<<"Too many warnings. Suppressing\n"<<endl;
					count++;
				}
			}

			total += atua.get_constant();

			pot_exit << pos_x << "   " << pos_y << "  " << total << endl;

			matrix[i][j] = atua.pot_to_head(total, 0, FALSE);

			data_exit_head << pos_x << "   " << pos_y << "  " << matrix[i][j]
					<< endl;

			for (polygon_it = et_lists.begin(); polygon_it != et_lists.end(); ++polygon_it) {

				if (within(make<point_2d> (pos_x, pos_y),
						polygon_it->get_polygon())) {
					matrix[i][j] = atua.pot_to_head(total,
							(*polygon_it).get_conductivity(), TRUE);
					//					printf("%i\t%i\tConductivity: %lf\thead:%lf\n", i, j,
					//					 (*polygon_it).get_conductivity(),
					//					 matrix[i][j]);

				}

			}

			data_exit << pos_x << "   " << pos_y << "  " << matrix[i][j]
					<< endl;

			pos_y += increment;
		}
		pos_y = 0;
		pos_x += increment;
	}

	cout <<"Exporting heads with correction - Done\n"<<endl;
	matrix.clear();
	data_exit.close();
	pot_exit.close();
	data_exit_head.close();
}

