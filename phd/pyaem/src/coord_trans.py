'''
Created on 05/03/2011

@author: ispmarin
'''

import scipy
import cmath



def smallz_to_Z(z, z1, z2):
    return (z - 0.5 * (z1 + z2)) / (0.5 * (z2 - z1))

def Z_to_smallz(Z, z1, z2):
    return Z * 0.5 * (z2 - z1) + 0.5 * (z1 + z2)

def chi_to_smallz(chi, z1, z2):
    return 0.5 * (chi + 1.0 / chi) * 0.5 * (z2 - z1) + 0.5 * (z1 + z2)

def chi_to_smallz_cos(chi, z1, z2):
    return scipy.cos(chi) * 0.5 * (z2 - z1) + 0.5 * (z1 + z2)
    #return mpmath.cos(chi) * 0.5 * (z2 - z1) + 0.5 * (z1 + z2)

def smallz_to_chi(z,z1,z2):
    Z = smallz_to_Z(z,z1,z2)
    return Z + scipy.sqrt(Z+1)*scipy.sqrt(Z-1)
    #return Z + mpmath.sqrt(Z+1)*mpmath.sqrt(Z-1)
def c_smallz_to_Z(z,center,R):
    return (z - center)/R

def c_theta_to_smallz(theta,center,R):
    
    return center + R*scipy.exp(1j*theta)

def ellipse_to_circle(z,c,center,angle, R):
    rot = scipy.exp(-1j*angle)
    csi = (z - center) * rot + scipy.sqrt((z - center) * rot - c) * scipy.sqrt(( z - center) * rot + c)
    
    return c_smallz_to_Z(csi, center, R)
    
def unit_circle_to_circle(z,R):
    
    return  R*scipy.exp(1j*cmath.phase(z))

def theta_unit_circle_to_circle(theta,R):
    
    return  R*scipy.exp(1j*theta)


def ellipse_to_unit_circle(z, c, center, R, angle):
    rot = scipy.exp(-1j*angle)
    
    return ((z - center)*rot  + scipy.sqrt((z - center)*rot  - c)*scipy.sqrt((z -center)*rot  + c) ) / R

def unit_circle_to_ellipse(z,c,center,R, angle):
    
    rot = scipy.exp(1j*angle)
    
    z = unit_circle_to_circle(z,R)
    return (0.5*(z + c*c/(z))*rot + center)

def theta_unit_circle_to_ellipse(theta,c,center,R, angle):
    
    rot = scipy.exp(1j*angle)
    
    z = theta_unit_circle_to_circle(theta,R)
    return (0.5*(z + c*c/(z))*rot + center)

