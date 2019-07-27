'''
Created on 07/03/2011

@author: ispmarin
'''
import numpy as np

import Crack
import gen_ellipse
import LineDouble
import CircleInhom
import EllipseInhom


def file_len(fname):
    i = 0
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1



def insert_crack(crack_list, filename, n):
    
    if file_len(filename) == 1 : 
        data = np.loadtxt(filename)
        k_int = data[0]
        aperture = data[1]
        z1 = complex(data[2], data[3])
        z2 = complex(data[4], data[5])
        crack_list.append(Crack.Crack(z1, z2, n, k_int, aperture))
    else:
        data = np.loadtxt(filename)
        for field in data:
            k_int = field[0]
            aperture = field[1]
            z1 = complex(field[2], field[3])
            z2 = complex(field[4], field[5])
        
            crack_list.append(Crack.Crack(z1, z2, n, k_int, aperture))

def insert_circle(circle_list, filename, n):
    elem_id = 0
    if file_len(filename) == 1 : 
        data = np.loadtxt(filename)
        k_int = data[0]
        R = data[1]
        center = complex(data[2], data[3])
        circle_list.append(CircleInhom.CircleInhom(elem_id,center, R, k_int, n))
    else:
        data = np.loadtxt(filename)
        for field in data:
            data = np.loadtxt(filename)
            k_int = field[0]
            R = field[1]
            center = complex(field[2], field[3])
            circle_list.append(CircleInhom.CircleInhom(elem_id,center, R, k_int, n))         
            elem_id +=1
            
            
def insert_ellipse(ellipse_list, filename, n):
    elem_id = 0
    if file_len(filename) == 1 : 
        data = np.loadtxt(filename)
        k_int = data[0]
        a = data[1]
        b = data[2]
        center = complex(data[3], data[4])
        angle = data[5]
        ellipse_list.append(EllipseInhom.EllipseInhom(elem_id,a,b, center, angle,n, k_int))
    else:
        data = np.loadtxt(filename)
        for field in data:
            data = np.loadtxt(filename)
            k_int = field[0]
            a = field[1]
            b = field[2]
            center = complex(field[3], field[4])
            angle = field[5]
            ellipse_list.append(EllipseInhom.EllipseInhom(elem_id,a,b,center, angle, n, k_int))        
            elem_id +=1
            
def insert_line_double(line_double_list, filename, elem_id, n, j):
    
    if file_len(filename) == 1 : 
        data = np.loadtxt(filename)
        k_int = data[0]
        far_field = data[1]
        z1 = complex(data[2], data[3])
        z2 = complex(data[4], data[5])
        line_double_list.append(LineDouble.LineDouble(elem_id, z1, z2, n, j, k_int, far_field))
    else:
        data = np.loadtxt(filename)
        for field in data:
            k_int = field[0]
            far_field = field[1]
            z1 = complex(field[2], field[3])
            z2 = complex(field[4], field[5])
        
            line_double_list.append(LineDouble.LineDouble(elem_id, z1, z2, n, j, k_int, far_field))

    
def gen_line_double_ellipse_by_box(line_double_list, elem_id, n, j, x_inf, x_sup, b, subdivisions, k_int,far_field):
    
    ellipse = gen_ellipse.CreateEllipse.by_box(x_inf,x_sup,b,subdivisions)
    ellipse.create_ellipse()
    
    for i in xrange(0,len(ellipse.vertices)-1):
        line_double_list.append(LineDouble.LineDouble(elem_id, ellipse.vertices[i], ellipse.vertices[i+1], n, j, k_int, far_field))
    
    line_double_list.append(LineDouble.LineDouble(elem_id, ellipse.vertices[-1], ellipse.vertices[0], n, j, k_int, far_field))

def gen_line_double_ellipse_by_semiaxes(line_double_list, elem_id, n, j, a, b, center, subdivisions, k_int,far_field):
    
    ellipse = gen_ellipse.CreateEllipse.by_semi_axis(a, b, 0, center, subdivisions)
    ellipse.create_ellipse()
    
    for i in xrange(0,len(ellipse.vertices)-1):
        line_double_list.append(LineDouble.LineDouble(elem_id, ellipse.vertices[i], ellipse.vertices[i+1], n, j, k_int, far_field))
    
    line_double_list.append(LineDouble.LineDouble(elem_id, ellipse.vertices[-1], ellipse.vertices[0], n, j, k_int, far_field)) 
   
