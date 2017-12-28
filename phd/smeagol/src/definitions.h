/*
 * constants.h
 *
 *  Created on: 19/10/2009
 *      Author: ispmarin
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_


#define mat(matriz,coluna,i,j) (matriz[i*coluna+j])
#define mat_transp(matriz,coluna,i,j) (matriz[j*coluna+i])
const double Pi =3.14159265;

const bool TRUE =1;
const bool FALSE = 0;

const int kind_reference_point =0;
const int kind_well = 1;
const int kind_ls = 2;
const int kind_ld =3;
const int kind_aq =4;
const int kind_cf = 5;

inline int getIndex(int row, int col,int NCOLS) {return row*NCOLS+col;}



#endif /* CONSTANTS_H_ */
