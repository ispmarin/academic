/*
 * insert_from_file.cpp
 *
 *  Created on: 01/06/2009
 *      Author: ispmarin
 */
#include "insert_from_file.h"
int ld_counter_id=0;

void open_run_file(const char * file, aem_container &the_list,
		aem_container &the_u_list, vector<aux_polygon> &my_polygon,
		aux_functions &aquifer_reference) {

	string elem_run_file = "0";
	string elem_type = "0";
	ifstream run_file;
	run_file.open(file);

	if (!run_file) {
		cerr << "Error! Run file cannot be opened! \n" << endl;
		exit(-1);
	}

	insert_from_file insert_handler;
	while (run_file >> elem_run_file >> elem_type) {
		insert_handler.insert_boost_file(aquifer_reference, elem_run_file,
				elem_type, the_list, the_u_list, my_polygon);
	}

}

export_to_file open_plot_file(const char * file) {
	ifstream plot_file;

	double x_min = 0, y_min = 0, x_max = 0, y_max = 0,increment=0;

	plot_file.open(file);
	plot_file >> x_min >> y_min >> x_max >> y_max >> increment;

	export_to_file plot_iso(x_min,y_min,x_max,y_max, increment);
	return plot_iso;


}
aux_functions open_reference_file(const char * aq_file, const char * ref_file) {
	ifstream aquifer_file;
	ifstream reference_file;
	double height = 0, base = 0, condut = 0, x = 0, y = 0, head = 0;

	aquifer_file.open(aq_file);
	reference_file.open(ref_file);

	reference_file >> x >> y >> head;
	aquifer_file >> condut >> base >> height;

	complex<double> temp(x, y);
	aux_functions temp_aquifer_reference(condut, base, height, temp, head);
	return temp_aquifer_reference;

}


insert_from_file::insert_from_file() {

}

insert_from_file::~insert_from_file() {
	//delete ref_aquifer;
}

int insert_from_file::convert_string_to_int(string type) {
	if (type == "w") {
		return kind_well;
	} else if (type == "ls") {
		return kind_ls;
	} else if (type == "ld") {
		return kind_ld;
	} else if (type == "cf") {
		return kind_cf;
	} else if (type == "aq") {
		return kind_aq;
	} else if (type == "rp") {
		return kind_reference_point;
	} else {
		exit(1);
	}

}

void insert_from_file::insert_boost_file(aux_functions aquifer_reference,
		string elem_run_file, string type, aem_container & the_list,
		aem_container & the_u_list,vector<aux_polygon> &vec_polygon_ld) {
	double x1 = 0;
	double y1 = 0;
	double x2 = 0;
	double y2 = 0;
	double R = 0;
	double head = 0;
	int is_know = 0;
	double discharge = 0;
	double Qx = 0;
	double Qy = 0;
	int int_type = 0;
	double condutivity = 0;
	complex<double> set_previous;
	complex<double> set_next;
	double aq_conductivity = aquifer_reference.get_aquifer_conductivity();


	ifstream data_file;
	data_file.open(elem_run_file.c_str()); //this is c_str(), a f&$ing way to convert a string to char *.
	//const charT* c_str() const >> Returns a pointer to a null-terminated array of characters representing the string's contents.

	//checking for error
	if (!data_file) {
		cerr << "Error! Data file " << data_file << "file cannot be opened! \n"
				<< endl;
		exit(-1);
	}
	vector <line_doublet * >::iterator current;
	vector <line_doublet * >::iterator previous;
	vector <line_doublet * >::iterator next;
	vector <line_doublet * > ld_coords;

	complex<double> previous_coord(-1, -1);
	complex<double> ini;


	vector<point_2d> temp_ld_coords;
	polygon_2d temp_polygon_ld;
	aux_polygon temp_aux_polygon_ld;
	//vector<aux_polygon> vec_polygon_ld;

	int_type = convert_string_to_int(type);//convert the input string to the well_id line_sink_id line_doublet notation

	switch (int_type) {

	case kind_well:
		cout <<"Inserting Wells"<<endl;
		while (data_file >> x1 >> y1 >> R >> head >> discharge >> is_know) {
			complex<double> tp(x1, y1);
			if (is_know) {
				the_list.add_element(new well(is_know, head, discharge, tp, R,kind_well));
			} else {
				the_u_list.add_element(
						new well(is_know, head, discharge, tp, R, kind_well));
			}

		}
		break;
	case kind_ls:
		cout <<"Inserting Line Sinks"<<endl;
		while (data_file >> x1 >> y1 >> x2 >> y2 >> head >> discharge
				>> is_know) {
			complex<double> tp(x1, y1);
			complex<double> tp2(x2, y2);
			if (is_know) {
				the_list.add_element(new line_sink(is_know, head, discharge,
						tp, tp2, kind_ls));
			} else {
				the_u_list.add_element(new line_sink(is_know, head, discharge,
						tp, tp2, kind_ls));
			}

		}
		break;
	case kind_ld:
		cout <<"Inserting Line Doublets"<<endl;
		data_file >> condutivity;

		data_file.clear();
		data_file.seekg(0, ios::beg);
		data_file >> condutivity; //argh! but needed to jump the first line.


		while (data_file >> x1 >> y1) {

			complex<double> tp(x1, y1);


			ld_coords.push_back(new line_doublet(is_know, head, discharge, tp,
											condutivity, aq_conductivity, kind_ld));
			temp_ld_coords.push_back(make<point_2d> (x1, y1));

		}

		//calculates here part of the matrix, maybe adding to a list


		previous = ld_coords.end() - 1;
		current = ld_coords.begin();
		next = ld_coords.begin() + 1;

		for (; current != ld_coords.end(); ++next, ++current, ++previous) {

			if (next == ld_coords.end()) {
				next = ld_coords.begin();
			}
			if (previous == ld_coords.end()) {
				previous = ld_coords.begin();
			}
			(*current)->set_next(*next);
			(*current)->set_previous(*previous);
			the_u_list.add_element(*current);

		}
				assign(temp_polygon_ld, temp_ld_coords);
				correct(temp_polygon_ld);

				cout<< "area "<<area(temp_polygon_ld)<<endl;

				temp_aux_polygon_ld.set_conductivity(condutivity);
				temp_aux_polygon_ld.set_polygon(temp_polygon_ld);
				vec_polygon_ld.push_back(temp_aux_polygon_ld);
				//cout <<"id: "<<ld_counter_id<<endl;

				ld_counter_id++;

		break;

	case kind_cf:

		while (data_file >> Qx >> Qy) {
			complex<double> tp(0, 0);
			the_list.add_element(new constant_flux(0,0,0,tp,kind_cf,Qx,Qy));

		}
		break;

	default:
		cerr << "Element not valid! Wrong number or not implemented" << endl;
		exit(-1);
	}

	cout << "Exiting insert_elements" << endl;


}

