/*
 * elements.h
 *
 *  Created on: 27/05/2009
 *      Author: ispmarin
 */

#ifndef ELEMENTS_H_
#define ELEMENTS_H_

//GGL
#include <boost/geometry/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>

#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <stdio.h>
#include <cstdlib>
#include "definitions.h"

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/assert.hpp>

using namespace boost::geometry;
using namespace std;


//++++++++++++++++++++++++++++++++++++++++++++++++++
//AQUIFER
class aquifer {
private:
	double height;
	double base;
	double condutivity;
protected:
public:
	aquifer(double set_height, double set_base, double set_condutivity) :
		height(set_height), base(set_base), condutivity(set_condutivity) {

	}

	aquifer(const aquifer &aq) :
		height(aq.get_height()), base(aq.get_base()), condutivity(
				aq.get_conduvitity()) {
	}

	~aquifer() {
	}

	double get_height() const {
		return height;
	}

	double get_base() const {
		return base;
	}

	double get_conduvitity() const {
		return condutivity;
	}

};

//++++++++++++++++++++++++++++++++++++++++++++++++++
//REFERENCE
class reference {
private:
	complex<double> position;
	double head;
protected:
public:
	reference(complex<double> set_position, double set_head) :
		position(set_position), head(set_head) {
	}

	reference(const reference & ref) :
		position(ref.get_position()), head(ref.get_head()) {
	}

	~reference() {
	}

	complex<double> get_position() const {
		return position;
	}

	double get_head() const {
		return head;
	}

};

//++++++++++++++++++++++++++++++++++++++++++++++++++
//AEM_ELEMENT
class aem_element: boost::noncopyable {
private:
	int elem_is_know;
	double head;
	double discharge;
	int type;

	virtual aem_element* do_clone() const = 0;

protected:
	complex<double> position;

public:
	aem_element(const aem_element& r) :
		elem_is_know(r.elem_is_know), head(r.head), discharge(r.discharge), type(r.type), position(r.position) {
	}

	explicit aem_element(int set_is_know, double set_head, double set_discharge,
			int ttype, complex<double> set_position=0) :
		elem_is_know(set_is_know), head(set_head), discharge(set_discharge),
				type(ttype), position(set_position) {

	}
	aem_element();

	virtual ~aem_element() {
	}

	int is_know() {
		return elem_is_know;
	}

	double get_discharge() {
		return discharge;
	}

	double get_head() {
		return head;
	}

	int get_type() {
		return type;
	}

	complex<double> coord_transform(complex<double> center,
			complex<double> next, complex<double> point);

	void set_discharge(double set_discharge) {
		discharge = set_discharge;
	}
	;
	virtual complex<double> get_control_point() {
		return position;
	}

	virtual double get_coeff(complex<double> z) {
		cout << "Should not be here!" << endl;
		exit(-1);
		return -1;
	}
	virtual double get_final_potential(complex<double> z) {
		return discharge * get_coeff(z);
	}
	;
	virtual aem_element *get_next() {
		return NULL;
	}

	virtual aem_element *get_prev() {
		return NULL;
	}

	aem_element* clone() const {
		return do_clone();
	}

};

//ALL DERIVED CLASSES FROM AEM_ELEMENT HAS TO BE BELOW THE AEM_ELEMENT DEFINITION


//++++++++++++++++++++++++++++++++++++++++++++++++++
//WELL
class well: public aem_element {
private:
	double radius;

	virtual aem_element* do_clone() const {
		return new well(*this);
	}

protected:
public:
	well(int set_is_know, double set_head, double set_discharge, complex<
			double> set_position, double set_R, int ttype) :
		aem_element(set_is_know, set_head, set_discharge, ttype, set_position),
				radius(set_R) {
		//radius = set_R;
		cout << "Well radius: " << radius << endl;
	}
	well(){}

	~well() {}
	double get_coeff(complex<double> z);

};

//++++++++++++++++++++++++++++++++++++++++++++++++++
//LINE SINK
class line_sink: public aem_element {
private:
	complex<double> position_mid;
	complex<double> position_end;

	virtual aem_element* do_clone() const {
		return new line_sink(*this);
	}
protected:
public:
	line_sink(int set_is_know, double set_head, double set_discharge,
			complex<double> set_position, complex<double> set_position_2, int ttype) :
		aem_element(set_is_know, set_head, set_discharge, ttype, set_position),
				position_end(set_position_2) {

		position_mid = (position + position_end) / 2.0;

	}
	line_sink(){}
	~line_sink(){}
	double get_coeff(complex<double> z);

	complex<double> get_control_point() {return position_mid;}

};

//++++++++++++++++++++++++++++++++++++++++++++++++++
//LINE DOUBLET
class line_doublet: public aem_element {
private:
	aem_element * previous;
	aem_element * next;
	double ld_condutivity;
	double external_condutivity;

	virtual aem_element* do_clone() const {
		return new line_doublet(*this);
	}

protected:
	double get_ld_conductivity() const {
		return ld_condutivity;
	}

	double get_external_condutivity() const {
		return external_condutivity;
	}


public:

	line_doublet(int set_is_know, double set_head,
			double set_discharge, complex<double> set_position,
			double set_ld_condutivity, double set_external_condutivity, int ttype) :
		aem_element(set_is_know, set_head, set_discharge, ttype, set_position),
				ld_condutivity(set_ld_condutivity), external_condutivity(
						set_external_condutivity) {

		cout << "ld:  " << ld_condutivity << "\t ex_ld: " << external_condutivity
				<< "  previous " << previous << "  current  " << position
				<< "  next  " << next << endl;
	}


	~line_doublet(){}

	aem_element * get_previous() const {
		return previous;
	}

	aem_element * get_next() const {
		return next;
	}

	void set_previous(aem_element * set_previous) {
		previous = set_previous;
	}
	void set_next(aem_element * set_next) {
		next = set_next;
	}

	double get_coeff(complex<double> z);
	double get_final_potential(complex<double> z);

	complex<double> get_external_point();
	double get_angle();

};

class constant_flux: public aem_element {
private:
	double Qx;
	double Qy;

	virtual aem_element* do_clone() const {
		return new constant_flux(*this);
	}
public:

	constant_flux(int set_is_know, double set_head, double set_discharge,
			complex<double> set_position,int ttype, double set_Qx, double set_Qy) :
	aem_element(set_is_know, set_head, set_discharge, ttype, set_position), Qx(set_Qx), Qy(set_Qy) {

	}

	double get_coeff(complex<double> z) {
		return -Qx * real(z) - Qy * imag(z);
	}
	double get_final_potential(complex<double> z) {
		return -Qx * real(z) - Qy * imag(z);
	}

};

class aem_container {

private:
	typedef boost::ptr_vector<aem_element> aem_element_list;
	aem_element_list the_list;

public:
	typedef aem_element_list::iterator aem_element_iterator;
	typedef aem_element_list::size_type size_type;
	typedef aem_element_list::auto_type aem_element_transport;

	aem_container() {
	}
	;
	aem_container(aem_element_iterator begin, aem_element_iterator end) :
		the_list(begin, end) {
	}
	;

	aem_element_iterator begin() {
		return the_list.begin();
	}

	aem_element_iterator end() {
		return the_list.end();
	}

	size_type size() const {
		return the_list.size();
	}

	void add_element(aem_element* element) {

		the_list.push_back(element);

	}

	auto_ptr<aem_element_list> return_aem_list() {
		return the_list.release();
	}

	void add_vectors(auto_ptr<aem_element_list> other_list) {
		the_list.transfer(the_list.end(), *other_list);
		BOOST_ASSERT(other_list->empty());
	}

	aem_element_transport remove_aem_element(aem_element_iterator to_remove) {
		if (to_remove == end())
			cout << "Removing the last!" << endl;

		return the_list.release(to_remove);
	}

	double get_coeff_from_index(int i, complex<double> Z) {
		return the_list[i].get_coeff(Z);
	}

	aem_element* clone(const aem_element& a) {
		return a.clone();
	}

};

aem_element* new_clone(const aem_element& a);



#endif /* ELEMENTS_H_ */
