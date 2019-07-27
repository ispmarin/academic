'''
Created on 12/05/2011

@author: ispmarin
'''
import coord_trans
import numpy as np
import scipy

class EllipseInhom(object):
    '''
    classdocs
    '''


    def __init__(self,sid,sa,sb,scenter,sangle,sn,sk_int):
        '''
        Constructor
        '''
        self.elem_id = sid
        self.a = sa
        self.b = sb
        self.center = scenter
        self.angle = sangle
        self.n = sn
        self.real_k_int = sk_int
        self.c =  scipy.sqrt(self.a*self.a - self.b*self.b)
        self.coeffs = np.zeros(self.n, dtype=complex)
        self.R = (self.a + self.b)
        
    def __repr__(self):
        return "Ellipse Inhom(id=%s, center_x=%s, center_y=%s, n=%d, a=%s, b=%s, angle=%s, k_it=%s, c=%s, R=%s)" % (self.elem_id, 
                                                   self.center.real,
                                                   self.center.imag,
                                                   self.n,
                                                   self.a,
                                                   self.b,
                                                   self.angle,
                                                   self.real_k_int,
                                                   self.c,
                                                   self.R)
 
 
    def get_potential(self,z):
        
        Z = coord_trans.ellipse_to_unit_circle(z, self.c, self.center, self.R, self.angle)
        
        #print 'ellipse', Z
        radius = (Z.real * Z.real + Z.imag * Z.imag)
        
        if (radius <= 1):
            return self.coeffs[0].real + np.sum(self.coeffs[n]* np.power(Z,n) for n in xrange(1,self.n) )
        else:
            return -np.sum(self.coeffs[n].conjugate()*np.power(1/Z,n) for n in xrange(1,self.n))

