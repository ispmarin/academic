'''
Created on 05/03/2011

@author: ispmarin
'''
import numpy as np
import cmath
import math

from coord_trans import smallz_to_Z
#import aux_functions

class Crack(object):
    '''
    The Crack element
    '''
    def __init__(self, sz1, sz2, sn,sk_int,saperture):
        '''
        Constructor
        '''
        if sz1.real > sz2.real and sz1.imag > sz2.imag:
            temp = sz1
            sz1=sz2
            sz2=temp
        
        self.z1 = sz1
        self.z2 = sz2
        self.n = sn #this is the max n, the sum will be made from 0 to n-1
        self.coeffs = np.zeros(self.n)
        self.real_k_int = sk_int
        self.aperture = saperture
        self.L = cmath.sqrt((self.z2 -self.z1)*(self.z2 - self.z1).conjugate()).real
        self.angle = math.atan2(-(self.z2 - self.z1).imag, -(self.z2 - self.z1).real) + math.pi#cmath.phase(self.z2 - self.z1).real  - math.pi/2.0
        
        #print 'Crack', self.z1, self.z2, 'n', self.n, self.coeffs, 'real_k_int', self.real_k_int,'aperture', self.aperture, 'length', self.L, 'angle', self.angle
    
    def __repr__(self):
        return "Crack(z1=%s, z2=%s, n=%d, k_it=%s, aperture=%s)" % (self.z1, 
                                                   self.z2, 
                                                   self.n,
                                                   self.real_k_int,
                                                   self.aperture)
    
    def get_potential(self, z): 
        Z = smallz_to_Z(z, self.z1, self.z2)
        
        return np.sum(self.coeffs * (Z - np.sqrt(Z - 1.0) * np.sqrt(Z + 1.0)) ** np.arange(self.n))
    
    def get_potential_m(self, z, n): 
        Z = smallz_to_Z(z, self.z1, self.z2)
        
        return np.power((Z - np.sqrt(Z - 1.0) * np.sqrt(Z + 1.0)),n)
    
    def get_discharge(self, z):
        '''
        returns the total discharge of the element. np.sum sums from 0 to n-1
        '''
        Z = smallz_to_Z(z, self.z1, self.z2)
        
#        if abs(Z.real) <= 1.0 + aux_functions.elem_distance and abs(Z.imag) <= aux_functions.elem_distance:
#        #if (Z.real) < 1.0 + aux_functions.elem_distance and (Z.real) > 1.0 - aux_functions.elem_distance and abs(Z.imag) <= aux_functions.elem_distance:
#        #    print 'mylself'
#            return complex (0.0, 0.0);
        
        
        return np.sum(np.arange(self.n) * self.coeffs * (1.0 / (cmath.sqrt(Z - 1.0) * cmath.sqrt(Z + 1.0))) 
                            * ((Z - cmath.sqrt(Z - 1.0) * cmath.sqrt(Z + 1.0)) ** np.arange(self.n))) * (2.0 / (self.z2 - self.z1))
            
        
    def get_discharge_m(self, z, n):

        Z = smallz_to_Z(z, self.z1, self.z2)
        
        return n * (2.0 / (self.z2 - self.z1)) * \
     ((Z - cmath.sqrt(Z - 1.0) * cmath.sqrt(Z + 1.0)) ** n ) / (cmath.sqrt(Z - 1.0) * cmath.sqrt(Z + 1.0)) 

        
