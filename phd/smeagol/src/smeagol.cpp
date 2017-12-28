#include "smeagol.h"

int main(int argc, char ** argv) {

	const char * aq_file = "aquifer.dat";
	const char * ref_file = "reference.dat";
	const char * plot_file = "plot.dat";

	aem_container the_list;
	aem_container the_u_list;
	vector<aux_polygon> my_polygon;


	aux_functions aquifer_reference = open_reference_file(aq_file, ref_file);
	open_run_file(argv[1], the_list,the_u_list,my_polygon,aquifer_reference);

	cout << "Setting the lists" << endl;
	cout << "know_vector size: " << the_list.size() << endl;
	cout << "unknow_vector size: " << the_u_list.size() << endl;

	//COMPUTATION
	long int size = (the_u_list.size() + 1);
	numeric_interface calculation_handler(size);
	calculation_handler.solve_system(aquifer_reference, the_list, the_u_list); //parameters passed by reference!!

	//PLOTTING
	export_to_file plot_iso = open_plot_file(plot_file);
	plot_iso.plot_surface_cond(aquifer_reference, the_u_list, my_polygon);

	cout << "Finished Successfully" << endl;

	return 0;
}


