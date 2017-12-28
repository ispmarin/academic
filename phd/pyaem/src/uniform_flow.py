'''
Created on 07/03/2011

@author: ispmarin
'''
import cmath
import math

class UniformFlowData(object):
    '''
    Uniform flow module. 
    Q0: uniform flow strength
    angle: input angle. Should ALWAYS enter in DEGREES, and converted inside to radians.
    '''


    def __init__(self, sQ0, sbeta):
        '''
        Constructor
        '''
        self.Q0 = sQ0
        self.beta = math.radians(sbeta)
        self.constant = 0.0
        #print ('uniform flow: ', 'Q0', self.Q0, 'beta', self.beta, 'beta angle', sbeta, \
        #'constant', self.constant)  
    
    def get_potential(self, z):
        return -self.Q0 * z * cmath.exp(-1j*self.beta) + self.constant
    
    def get_discharge(self, alpha):
        '''
        Gets the discharge from the uniform flow. The alpha angle should be entered in radians.
        '''
        return self.Q0*cmath.exp(1j*(alpha - self.beta))
    
    def get_potential_no_c(self, z):
        '''
        Returns the potential without the constant.
        '''
        return -self.Q0*z*cmath.exp(-1j*self.beta)
    