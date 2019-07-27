'''
Created on 07/03/2011

@author: ispmarin
'''

import math
import numpy as np
import scipy

import Crack
import LineDouble
import CircleInhom
import coord_trans



elem_distance = 1e-7

class BreakIt(Exception): pass

def head_to_pot(head, aq):
    return 0.5 * aq.k_ext * (head ** 2)

def pot_to_head(pot,aq):
    return scipy.sqrt(2 * pot / aq.k_ext)

def get_angle_(z, element):
    epsilon = 0.0001
    if z == element.z1: 
        z = z - complex(epsilon, epsilon)
    
    if z == element.z2:
        z = z - complex(epsilon, epsilon)
    
    return math.atan2(((z - element.z2) / (z - element.z1)).imag, \
    ((z - element.z2) / (z - element.z1)).real)


def get_inside_k(element_list, aq, z):

    angle = 0
    if type(element_list[0]) == CircleInhom.CircleInhom or type(element_list[0]) == CircleInhom.ExactSolution:
        
        for element in element_list:
            Z = coord_trans.c_smallz_to_Z(z, element.center, element.R)
            radius = (Z.real * Z.real + Z.imag * Z.imag)
            if radius <= 1.0:
                return element.real_k_int

    else:
        for element in element_list:
            angle = angle + ((element.real_k_int - aq) / (2.0 * math.pi)) * \
            get_angle_(z, element)
    
    return angle + aq


def potential_value_no_c(z, element_list, u_flow):
    potential = complex(0, 0)
    
    for element in element_list:
        potential = potential + element.get_potential(z)
    
    potential = potential + u_flow.get_potential_no_c(z) 
    
    return potential

def head_value(z, element_list, u_flow, aq):
    potential = complex(0, 0)
    
    for element in element_list:
        potential = potential + element.get_potential(z)
        
    if type(element_list[0]) != CircleInhom.ExactSolution:
        potential = potential + u_flow.get_potential(z)
    
    
    if type(element_list[0]) == Crack.Crack:
        potential = potential / aq.k_ext
    else:
        potential = potential / get_inside_k(element_list, aq.k_ext, z)
        
        
    if potential.real < 0.0:
        potential = 0.0
    
    head = scipy.sqrt(2.0 * potential.real)
    
    return head

def final_potential(z, elem_list, u_flow):
    potential = complex(0, 0)
    
    for element in elem_list:
        potential = potential + element.get_potential(z)
        
    potential = potential + u_flow.get_potential(z) 

    return potential

        
def set_axis(element_list):

    z1 = np.zeros(len(element_list), dtype=complex)
    z2 = np.zeros(len(element_list), dtype=complex)
    
    elem_type = type(element_list[0])
    
    if elem_type == Crack.Crack or elem_type == LineDouble.LineDouble:
        for i, element in enumerate(element_list):
            z1[i] = element.z1
            z2[i] = element.z2
    
    elif elem_type == CircleInhom.CircleInhom or elem_type ==  CircleInhom.ExactSolution:
        for i, element in enumerate(element_list):
            z1[i] = element.center - 1j*element.R - element.R
            z2[i] = element.center + 1j*element.R + element.R
    else:
        print('elem type not recognized:', elem_type)    

    real_axis = np.concatenate((z1.real, z2.real))
    imag_axis = np.concatenate((z1.imag, z2.imag))
    
    max_x = np.amax(real_axis) + 30.0
    min_x = np.amin(real_axis) - 30.0
    max_y = np.amax(imag_axis) + 30.0
    min_y = np.amin(imag_axis) - 30.0
    
    if (max_x - min_x) < (max_y - min_y):
        num_steps = (max_y - min_y) *1.4#/ 2.0
    else:
        num_steps = (max_x - min_x) *1.4#/ 2.0
    
    return min_x, min_y, max_x, max_y, num_steps
