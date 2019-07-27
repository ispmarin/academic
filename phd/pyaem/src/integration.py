# -*- coding: utf-8 -*-
"""
Created on Sun Apr 17 13:13:49 2011

@author: -
"""

import cmath
import math
import numpy as np


import scipy.integrate
import coord_trans
from coord_trans import chi_to_smallz
#import aux_functions


def get_elem_discharge(crack_list, z):
    '''
    Return the discharge of all elements
    '''
    discharge = complex(0.0, 0.0)
    
    for element in crack_list:
        discharge = discharge + element.get_discharge(z)

    return discharge

def func_discharge(theta, n, element_list, being_calculated):
    z = coord_trans.chi_to_smallz_cos(theta, being_calculated.z1, 
                                          being_calculated.z2)
         
    return (get_elem_discharge(element_list, z)* cmath.exp(1j* \
        being_calculated.angle) ).real * math.sin(n * theta) * \
        math.sin(theta)
        
def func_discharge_u_flow(theta, n, u_flow, being_calculated):
         
    return u_flow.get_discharge(being_calculated.angle).real * math.sin(n * theta) * \
        math.sin(theta)
        
        
def func_discharge_m(theta, i, n, element, being_calculated):
    z = coord_trans.chi_to_smallz_cos(theta, being_calculated.z1, 
                                          being_calculated.z2)

        
    return (element.get_discharge_m(z, i) * cmath.exp(1j* \
        being_calculated.angle) ).real * math.sin(n * theta) * \
        math.sin(theta)


def get_elem_potential_ex(crack_list, being_calculated, z):
    '''
    Return the potential of all elements, excluding myself
    '''
    potential = complex(0.0, 0.0)
    
    for element in crack_list:
        if element != being_calculated:
            potential = potential + element.get_potential(z)

    return potential

def get_one_elem_pot_ex(theta, element, being_calculated, n, m, k, M):
    if element != being_calculated:
        z = coord_trans.chi_to_smallz_cos(theta, being_calculated.z1, being_calculated.z2)

        return element.get_potential_m(z, m).real * scipy.cos(theta * n)
    else:
        return 0

def get_elem_potential(crack_list, z):
    '''
    Return the potential of all elements
    '''
    potential = complex(0.0, 0.0)
    
    for element in crack_list:
        potential = potential + element.get_potential(z)

    return potential

def func_potential(theta, n, element_list, being_calculated, u_flow):
    z = coord_trans.chi_to_smallz_cos(theta, being_calculated.z1, 
                                          being_calculated.z2)
                                          
    return (get_elem_potential_ex(element_list,being_calculated, z) + u_flow.get_potential(z)).real \
    * math.cos(float(n) * theta)
    
def func_potential_m(theta, m, n, element, being_calculated):
    
    
    z = coord_trans.chi_to_smallz_cos(theta, being_calculated.z1, 
                                          being_calculated.z2)
                                          
    return (element.get_potential_m(z, m).real) * math.cos(float(n) * theta)
    
    
def func_potential_circle(theta, n, element_list, being_calculated, u_flow):
    z = coord_trans.c_theta_to_smallz(theta, being_calculated.center, 
                                          being_calculated.R)
                        
    return (get_elem_potential_ex(element_list, being_calculated, z) + u_flow.get_potential(z)).real \
    * cmath.exp(-1j*theta*n)
    
def func_potential_circle_m(theta, m, n, element, being_calculated,aq):
    z = coord_trans.c_theta_to_smallz(theta, being_calculated.center, 
                                          being_calculated.R)
    
    if element != being_calculated: #correct
        return element.get_potential_m(z,m).real * scipy.exp(-1j*theta*n)
    else:
        return 0
    
def func_u_flow(theta, n, u_flow, being_calculated):
    z = coord_trans.chi_to_smallz_cos(theta, being_calculated.z1, 
                                          being_calculated.z2)
    
    return u_flow.get_potential(z).real * math.cos(n*theta)

def func_u_flow_no_c(theta, n, u_flow, being_calculated):
    z = coord_trans.chi_to_smallz_cos(theta, being_calculated.z1, 
                                          being_calculated.z2)
    
    return u_flow.get_potential_no_c(z).real * math.cos(n*theta)
    
    #return u_flow.get_potential_no_c(z).real * math.cos(n*theta)
def func_u_flow_circle(theta, n, u_flow, being_calculated, elem_list, aq):
    
    z = coord_trans.c_theta_to_smallz(theta, being_calculated.center, being_calculated.R)
    
    return  u_flow.get_potential(z).real  * scipy.exp(-1j * theta * n) 
    #return  u_flow.get_potential(z).real  * scipy.cos( theta * n)
def func_u_flow_circle_no_c(theta, n, u_flow, being_calculated, elem_list, aq):
    
    z = coord_trans.c_theta_to_smallz(theta, being_calculated.center, being_calculated.R)
    
    return  u_flow.get_potential_no_c(z).real  * scipy.exp(-1j * theta * n) 
   
def func_u_flow_no_c_elem(theta, n,k,M, u_flow, being_calculated):
    #theta = 2 * scipy.pi * k/M
    z = coord_trans.chi_to_smallz_cos(theta, being_calculated.z1, being_calculated.z2)
    
    return  u_flow.get_potential_no_c(z).real  * scipy.cos( theta * n) 

def complex_quadrature(func, a, b,  *args):
    
    def real_func(theta, *args):
        return scipy.real(func(theta,*args))
    def imag_func(theta,*args):
        return scipy.imag(func(theta,*args))
        
    real_integral = scipy.integrate.quad(real_func, a, b, args=(args), limit=100)
    imag_integral = scipy.integrate.quad(imag_func, a, b, args=(args), limit=100)
    
    return real_integral[0] + 1j*imag_integral[0]

def constant_exp(theta, m, n, element, being_calculated,aq):
    
    return scipy.exp(-1j*theta*n)

def constant_cos(theta, m, n, element, being_calculated,aq):
    
    
    return scipy.cos(theta * n)
    #return mpmath.cos(theta * n)

def constant_sin_sin(theta, m, n, element, being_calculated):
    
    
    return scipy.sin(theta*n)*scipy.sin(theta)
    
class Integration(object):
    '''
    Class integration
    '''
    def __init__(self, low_a, upp_b, divs):
        self.lower_lim = low_a
        self.upper_lim = upp_b
        self.divisions = divs
        
        self.step = (self.upper_lim - self.lower_lim) / float(self.divisions)
        
        self.delta_thetaover2 = self.step / 2.0
        
    def test_functional(self,theta,n,b,c):
        return theta*theta 
        
        

        
    def basic_int(self, line_double_list, being_calculated, u_flow, n):
        '''
        Potential integration
        '''
        
        theta = self.delta_thetaover2
        
        integral = complex(0.0, 0.0)
        
        result = 0.0
        
        
        for _ in xrange(self.divisions):
    
            z = coord_trans.chi_to_smallz_cos(theta, being_calculated.z1, 
                                              being_calculated.z2)
        
            for element in line_double_list:
                
                if being_calculated != element:
             
                    integral = integral + element.get_potential(z)  
               
            integral = integral + u_flow.get_potential(z)
            
            result = result + integral.real * math.cos(float(n) * theta) * self.step
    
            theta = theta + self.step
            
            integral = complex(0.0, 0.0)
        
          
        return result
    
    def basic_int_discharge(self, crack_list, being_calculated, n):
        '''
        Discharge integration
        '''
        
        theta = self.delta_thetaover2
        
        integral = complex(0.0, 0.0)
        
        result = 0.0
        
        for _ in xrange(self.divisions):
            
            z = coord_trans.chi_to_smallz_cos(theta, being_calculated.z1, 
                                              being_calculated.z2)
            
            integral = get_elem_discharge(crack_list, z)        
            
            result = result + (integral * cmath.exp(1j* \
            being_calculated.angle) ).real * math.sin(float(n) * theta) * \
            math.sin(theta) * self.step
    
            theta = theta + self.step
            
            integral = complex(0.0, 0.0)
     
        return result
        
    def riemman_rule(self, functional, n, *args ):
        '''
        To calculate correctly the end points, the xrange must be -1, as in the
        rectangle rule.
        '''
        
        theta = self.delta_thetaover2
        
        integral = 0.0

        for _ in xrange( self.divisions ):
  
            integral = integral + functional(theta, n, *args ) 
    
            theta = theta + self.step
          
        
        return integral * self.step
      
    
    def trapezoidal_rule(self, functional, n, *args):
        """Approximate the definite integral of f from a to b by the
        composite trapezoidal rule, using N subintervals
        Note: the trapezoidal must put both ends to -DeltaThetaOver2 to 
        make sense in the context of removing the singular extremities.
        The xrange for number of divisions should not be changed.
        """
        
        a = self.delta_thetaover2
        b = math.pi - self.delta_thetaover2
     
        return (b-a) * ( functional(a, n, *args)/2.0 + \
        functional(b, n, *args)/2.0 + \
        sum([functional(a + (b-a)*k/self.divisions, n, *args) \
        for k in xrange(1,self.divisions )]) ) / self.divisions
        
        
    def rectangle_rule(self,functional, *args):
        '''
        Rectangle rule, exactly equal to the riemman sum if k - i with i = 0. Note that to exclude the final extremity, 
        the xrange must be from 1 to self.divisions-1.
        '''
        a = self.delta_thetaover2
        
        return self.step * np.sum([functional(a + self.step * k ,*args) for k in xrange( self.divisions  )])
    
    def rectangle_rule_m(self,functional, i, n, *args):
        '''
        Rectangle rule, exactly equal to the riemman sum if k - i with i = 0. Note that to exclude the final extremity, 
        the xrange must be from 1 to self.divisions-1.
        '''
        a = self.delta_thetaover2
        
        return self.step * np.sum([functional(a + self.step * k ,i, n, *args) for k in xrange( self.divisions  )])
           
 

        
        

