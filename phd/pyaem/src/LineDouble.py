'''
Created on 05/03/2011

@author: ispmarin
'''
import numpy as np
import cmath
import math


from coord_trans import smallz_to_Z
from coord_trans import smallz_to_chi
import aux_functions




class LineDouble(object):
    '''
    The Crack element
    '''
    def __init__(self, sid, sz1, sz2, sn, sj, sk_int, sfar_field):
        '''
        Constructor
        '''
        self.id = sid
        self.z1 = sz1
        self.z2 = sz2
        self.n = sn #this is the max n, the sum will be made from 0 to n-1
        self.j = sj
        self.coeffs = np.zeros(self.n)
        self.real_k_int = sk_int
        self.far_field = sfar_field
        self.L = cmath.sqrt((self.z2 - self.z1) * (self.z2 - self.z1).conjugate()).real
        self.angle = math.atan2(-(self.z2 - self.z1).imag, -(self.z2 - self.z1).real) + math.pi#cmath.phase(self.z2 - self.z1).real  - math.pi/2.0
        self.bj = np.zeros(self.j * (self.n + 1))
        self.an = np.empty([self.n,self.j])
        
        self.set_bj_an()
        
        
        
    def __repr__(self):
        return "Line Double(id=%s, z1=%s, z2=%s, n=%d, j=%s, k_it=%s, far_field=%s, L=%s)" % (self.id, 
                                                                                                         self.z1,
                                                   self.z2,
                                                   self.n,
                                                   self.j,
                                                   self.real_k_int,
                                                   self.far_field,
                                                   self.L)
    
        
    def set_bj_an(self):
        self.bj[0] = 0
        for i in xrange(1, self.j * (self.n + 1), 2):
            self.bj[i] = 1.0 / i
            
        #print self.bj

        #an_temp = np.zeros(self.j)
        for n in xrange(0, self.n): 
            for j in xrange(1, self.j): 

                if j <= n - 1: 
                    self.an[n,j] = 2.0 * (self.bj[n - j] - self.bj[n + j])
#                    an_temp[j] = 2.0 * (self.bj[n - j] - self.bj[n + j])
#                    assert(an_temp[j] == self.an[n,j])
                elif j >= n + 1: 
                    self.an[n,j] = -2.0 * (self.bj[j - n] + self.bj[n + j])
#                    an_temp[j] = -2.0 * (self.bj[j - n] + self.bj[n + j])
#                    assert(an_temp[j] == self.an[n,j])
                elif j == n: 
                    self.an[n,j] = -2.0 * self.bj[2 * j]
#                    an_temp[n,j] = -2.0 * self.bj[2 * j]
#                    assert(an_temp[j] == self.an[n,j])
                else:
                    print("should not be here!")

            
#            self.an.append(an_temp); 
#            an_temp = np.zeros(self.j)
#        
        #print 'an', self.an
        
    def get_far_field_correction(self, chi, n):
        
        #far_field_corr = complex(0.0, 0.0)
        #print chi
        #np.seterr(all='raise')
        return np.sum((np.power(chi, j) + np.power(1.0 / chi, j)) * float(self.bj[n - j]) for j in xrange(1, n))
        #return np.sum(float(self.bj[n - j]) * np.power(chi,-np.arange(1,self.j)))
        
        #for  j in xrange(1, n):

        #    far_field_corr = far_field_corr + (pow(chi, j) + pow(1.0 / chi, j)) * float(self.bj[n - j])
        
        
        #if abs(test - far_field_corr) > 1e-4:
        #    print 'ERROR'
        #    print (pow(chi, j) + pow(1.0 / chi, j))
        #return far_field_corr
    
    def get_far_field(self, chi, n):
        #far_field = complex(0.0, 0.0)
        #np.seterr(all='raise')
        return np.sum(self.an[n,1:self.j] * np.power(1.0/chi,np.arange(1,self.j)))
    
        #return np.sum(self.an[n,j] * pow(chi, -j) for j in xrange(1, self.j))
        #for j in xrange(1, self.j):
        #   far_field = far_field + self.an[n][j] * pow(1.0 / chi, j);
            
        #if abs(test - far_field) > 1e-4:
        #    print 'ERROR'
            
        #return far_field;
    
    def get_F_n(self, chi, n):
#        print (pow(chi, n) + pow(1.0 / chi, n)) \
#                * cmath.log((chi - 1.0) / (chi + 1.0)) + 2.0 * self.bj[n] + 2.0 \
#                * self.get_far_field_correction(chi, n)
        #np.seterr(all='raise')
        return (np.power(chi, n) + np.power(1.0 / chi, n)) \
                * np.log((chi - 1.0) / (chi + 1.0)) + 2.0 * self.bj[n] + 2.0 \
                * self.get_far_field_correction(chi, n)

            
    def get_potential(self, z): 
        Z = smallz_to_Z(z, self.z1, self.z2)
        
        if abs(Z.real) <= 1.0 + aux_functions.elem_distance and abs(Z.imag) <= aux_functions.elem_distance:
            return complex (0.0, 0.0)

        #else: 
        #np.seterr(all='raise')
        chi = smallz_to_chi(z, self.z1, self.z2)
        radius = (chi.real * chi.real + chi.imag * chi.imag)
        
            
        pot = complex(0.0, 0.0)

        if radius > self.far_field:
                
            pot = np.sum(self.coeffs[n] * self.get_far_field(chi, n) for n in xrange(0, self.n))
                #for n in xrange(0, self.n):
                #   pot = pot + self.coeffs[n] * self.get_far_field(chi, n)
            #print 'far field', pot
                
        else:
            pot = np.sum(self.coeffs[n] * self.get_F_n(chi, n) for n in xrange(0, self.n))
                #for n in xrange(0, self.n):
                #    pot = pot + self.coeffs[n] * self.get_F_n(chi, n)
                #print 'coeffs near', self.coeffs[1]
            #print 'near field', pot
        #print "pot ",(1.0/(2.0*math.pi))*pot 
        #print (1.0 / (2.0 * math.pi * 1j)) * pot
        return (1.0 / (2.0 * math.pi * 1j)) * pot


    def get_potential_m(self, z, n): 
        Z = smallz_to_Z(z, self.z1, self.z2)
        
        if abs(Z.real) <= 1.0 + aux_functions.elem_distance and abs(Z.imag) <= aux_functions.elem_distance:
            return complex (0.0, 0.0)

        chi = smallz_to_chi(z, self.z1, self.z2)
        radius = (chi.real * chi.real + chi.imag * chi.imag)
            
        pot = complex(0.0, 0.0)

        if radius > self.far_field:
                
            pot = self.get_far_field(chi, n) 

        else:
            
            pot = self.get_F_n(chi, n)
            
        return (1.0 / (2.0 * math.pi * 1j)) * pot
    
