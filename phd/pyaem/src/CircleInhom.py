'''
Created on 10/05/2011

@author: ispmarin
'''

import numpy as np
import scipy

import coord_trans



class CircleInhom(object):
    '''
    classdocs
    '''


    def __init__(self,sid,scenter,sR, sk_int, sn):
        '''
        
        Constructor
        '''
        self.id = sid
        self.center = scenter
        self.R = sR
        self.real_k_int = sk_int
        self.n = sn

        self.coeffs = scipy.zeros(self.n,dtype=complex)
    
    def __repr__(self):
        return "Circle Inhom(id=%s, center_x=%s, center_y=%s, n=%d, R=%s, k_it=%s)" % (self.id, 
                                                                                                         self.center.real,
                                                   self.center.imag,
                                                   self.n,
                                                   self.R,
                                                   self.real_k_int)
        
    def get_potential(self,z):
        Z = coord_trans.c_smallz_to_Z(z,self.center,self.R)
        #print 'circle', Z
        
        radius = (Z.real * Z.real + Z.imag * Z.imag)
        
        if (radius <= 1):
            return self.coeffs[0].real + np.sum(self.coeffs[n]* np.power(Z,n) for n in xrange(1,self.n) )
        else:
            return -np.sum(self.coeffs[n].conjugate()*np.power(1/Z,n) for n in xrange(1,self.n))
        
    def get_potential_m(self,z,n):
        Z = coord_trans.c_smallz_to_Z(z,self.center,self.R)
        
        radius = (Z.real * Z.real + Z.imag * Z.imag)
            
        if (radius <= 1.0):
            return np.power(Z,n)
        else:
            if n == 0:
                return 0
            else:
                return np.power(1/Z,n).real# - np.power(1/Z,n).imag 
        
class ExactSolution(object):
    
    def __init__(self,scenter,sR, sk_int,sQ,aq):
        self.center = scenter
        self.R = sR
        self.real_k_int = sk_int
        self.k_int = aq.k_ext
        self.Q = sQ
        self.k_ext = aq.k_ext
        self.C = 0
        
        #print('Exact Circle Parameters: ', 'center: ',self.center,'R',self.R, 'real_k_int',self.real_k_int,'Q', self.Q,'k_ext', self.k_ext)
        
        Z_R = coord_trans.c_smallz_to_Z(aq.reference_point,self.center,self.R)

        ZZ_R = (Z_R.real) * (Z_R.real) + (Z_R.imag) * (Z_R.imag)

        self.C = self.Q * self.R * (Z_R.real) * (1.0 + ((aq.k_ext - self.real_k_int) / (self.real_k_int + aq.k_ext)) * (1.0 / ZZ_R)) \
                    + 0.5 * aq.k_ext * aq.reference_head * aq.reference_head
        
                    
        print('Exact Constant', self.C)
    
    def get_potential(self,z):
        
        Z = coord_trans.c_smallz_to_Z(z, self.center, self.R)
        
        radius = (Z.real) * (Z.real) + (Z.imag) * (Z.imag)

        if radius <= 1.0:
            return -self.Q * Z * self.R * ((2 * self.real_k_int) / (self.k_ext + self.real_k_int)) + (self.real_k_int / self.k_ext) * self.C
        else: 
            return -self.Q * self.R * (Z + ((self.k_ext - self.real_k_int) / (self.k_ext + self.real_k_int)) * (1.0 / Z)) + self.C
        
    
    
    