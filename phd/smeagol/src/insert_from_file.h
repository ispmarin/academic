#ifndef INSERT_FROM_FILE_H_
#define INSERT_FROM_FILE_H_

#include <fstream>
#include <iostream>
#include <vector>

#include "definitions.h"
#include "elements.h"
#include "aux_functions.h"
#include "export.h"

using namespace ggl;


export_to_file open_plot_file(const char * file) ;
aux_functions open_reference_file(const char * aq_file, const char * ref_file);
void open_run_file(const char * file, aem_container &the_list,
		aem_container &the_u_list, vector<aux_polygon> &my_polygon,
		aux_functions &aquifer_reference);

class insert_from_file {

private:


public:
	insert_from_file();
	~insert_from_file();
	void insert_file(aux_functions aquifer_reference, string elem_run_file,
			string type,  vector<aem_element> &t_k_list,
			vector<aem_element> &t_uk_list,	vector<aux_polygon> &vec_polygon_ld);

	int convert_string_to_int(string type);

	void insert_boost_file(aux_functions aquifer_reference, string elem_run_file,
			string type, aem_container & the_list, aem_container & the_u_list,vector<aux_polygon> &vec_polygon_ld);

};

#endif /* INSERT_FROM_FILE_H_ */
