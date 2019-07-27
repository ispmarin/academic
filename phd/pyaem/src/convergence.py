'''
Created on 24/03/2011

@author: ispmarin
'''

import numpy as np

class convergence_element(object):
    '''
    classdocs
    '''


    def __init__(self, sn, snum_elem, stotal_size, stolerance, elem_list):
        '''
        Constructor
        '''
        
        self.n = sn
        self.num_elem = snum_elem
        self.coeffs= [[1.0]*self.n for _ in xrange(self.num_elem)]
        self.tolerance = stolerance
        
        for i, elem in enumerate(elem_list):
            for j in xrange(self.n):
                self.coeffs[i][j] = elem.coeffs[j]
        
    def check_convergence(self,elem_list):
        flag = 0
        
        for i, elem in enumerate(elem_list):
            for j in xrange(self.n):
                if abs((self.coeffs[i][j] - elem.coeffs[j]) / elem.coeffs[j]) < self.tolerance or -self.tolerance < elem.coeffs[j] < self.tolerance:
                    self.coeffs[i][j] = elem.coeffs[j]
                    flag = flag + 1
                else:
                    self.coeffs[i][j] = elem.coeffs[j]
                    flag = 0
                    
        if flag == self.n*self.num_elem  :
            return 1
        else:
            return 0
        
class convergence_vector(object):
     
    def __init__(self,stotal_size, sabs_tolerance, srel_tolerance, ini_vector):
        '''
        Constructor
        '''
        self.total_size = stotal_size
        self.coeffs= np.zeros(self.total_size)
        self.abs_tolerance = sabs_tolerance
        self.rel_tolerance = srel_tolerance
        
        
        for j in xrange(self.total_size):
                self.coeffs[j] = ini_vector[j]
    
    
    def check_convergence_vector(self,sol_vector):
        flag = 0
        
        for j in xrange(self.total_size):
            if abs((self.coeffs[j] - sol_vector[j]) / sol_vector[j]) < self.rel_tolerance or -self.abs_tolerance < sol_vector[j] < self.abs_tolerance or \
            abs(self.coeffs[j] - sol_vector[j]) < self.abs_tolerance:
                self.coeffs[j] = sol_vector[j]
                flag = flag + 1
            else:
                self.coeffs[j] = sol_vector[j]
                flag = 0
                
        if flag == self.total_size:
            return 1
        else:
            return 0
    