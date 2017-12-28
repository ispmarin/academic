'''
Created on 21/05/2011

@author: ispmarin
'''

import numpy as np


import coord_trans
import aux_functions

elem_distance = 1e-7

class ErrorAnalysis(object):


    def __init__(self,sx1,sx2,sy1,sy2,step,elem_list_1,elem_list_2,u_flow, kind_comparison):
        self.x1 = sx1
        self.x2 = sx2
        self.y1 = sy1
        self.y2 = sy2
        self.step = step
        
        print('Error Analysis set up...')
        
        dim_x = np.linspace(self.x1, self.x2, self.step)
        dim_y = np.linspace(self.y1, self.y2, self.step)
            
        len_x = len(dim_x)
        len_y = len(dim_y)
          
        self.a = np.zeros([len_x, len_y], dtype=complex)
        self.b = np.zeros([len_x, len_y], dtype=complex)
        
        if kind_comparison == 'head':
            for i, y in enumerate(dim_y):
                for j, x in enumerate(dim_x):
                    self.a[i][j] = aux_functions.head_value(complex(x, y), elem_list_1, u_flow)
                    self.b[i][j] = aux_functions.head_value(complex(x, y), elem_list_2, u_flow)
        elif kind_comparison == 'potential':
            for i, y in enumerate(dim_y):
                for j, x in enumerate(dim_x):
                    self.a[i][j] = aux_functions.final_potential(complex(x, y), elem_list_1, u_flow)
                    self.b[i][j] = aux_functions.final_potential(complex(x, y), elem_list_2, u_flow)
                    
        print('done.')

def compare_results_absolute(a, b):
    assert(a.shape == b.shape)
    
    absolute_error = np.abs((a - b))
    print ('Absolute error: ', np.max(absolute_error))
    return absolute_error

def compare_results_relative(a, b):
    assert(a.shape == b.shape)
    
    relative_error = np.abs(((a - b)/b))
    print( 'Relative error: ', np.max(relative_error))
    return relative_error


    
def check_ld_jump(line_double_list, u_flow, aq, tolerance):
    
    print('Checking Line Doublet Jump ')
    #print('Re(Plus) \t Im(Plus) \t Re(Min) \t Im(Min) \t abs(Plus-Min)')
    errors = np.zeros(0)
    for element in line_double_list:
        for pos in np.linspace(-0.9,0.9,20):
            pos_plus  = coord_trans.Z_to_smallz(complex(pos,  2*elem_distance), element.z1, element.z2)
            pos_minus = coord_trans.Z_to_smallz(complex(pos, -2*elem_distance), element.z1, element.z2)
            pot_plus  = aux_functions.final_potential(pos_plus, line_double_list,u_flow)/element.real_k_int
            pot_minus = aux_functions.final_potential(pos_minus, line_double_list, u_flow)/aq.k_ext
            
            if abs(pot_plus.real - pot_minus.real) > tolerance:
                errors = np.append(errors,abs(pot_plus.real - pot_minus.real))
                #print pot_plus.real,'\t', pot_minus.real, '\t',  abs(pot_plus.real - pot_minus.real)
    print 'Max jump error: ', np.max(errors)

def check_circle_jump(circle_list, u_flow, aq, tolerance):

    errors = np.zeros(0)
    print('Checking Circle Jump ')
    for  element in circle_list:
        for theta in np.linspace(0,2*np.pi,100):
            
            pot_plus  = aux_functions.final_potential(coord_trans.c_theta_to_smallz(theta, element.center, element.R + elem_distance), circle_list, u_flow)/aq.k_ext
            pot_minus = aux_functions.final_potential(coord_trans.c_theta_to_smallz(theta, element.center, element.R - elem_distance), circle_list, u_flow)/element.real_k_int
            
            if abs(pot_plus.real - pot_minus.real) > tolerance:
                errors = np.append(errors,abs(pot_plus.real - pot_minus.real))
                #print i, theta,'\t', pot_plus.real/aq.k_ext,'\t',pot_minus.real/element.real_k_int, '\t', abs(pot_plus.real - pot_minus.real)
    print 'Max jump error: ', np.max(errors)

