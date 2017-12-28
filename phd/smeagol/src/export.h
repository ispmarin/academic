/*
 * export.h
 *
 *  Created on: 04/06/2009
 *      Author: ispmarin
 */

#ifndef EXPORT_H_
#define EXPORT_H_

#include "elements.h"
#include "aux_functions.h"
#include <fstream>

class export_to_file {
	private:

		double x_min;
		double y_min;
		double x_max;
		double y_max;
		double increment;
	public:
		export_to_file(double set_x_min, double set_y_min, double set_x_max,
				double set_y_max, double set_increment) :
			x_min(set_x_min), y_min(set_y_min), x_max(set_x_max), y_max(
					set_y_max), increment(set_increment) {
		}

		void plot_surface_cond(aux_functions & atua, aem_container &the_list,
				vector<aux_polygon> polygon_list);
};

#endif /* EXPORT_H_ */
