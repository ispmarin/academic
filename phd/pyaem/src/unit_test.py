'''
Created on 05/03/2011

@author: ispmarin
'''
#TODO need more tests for Crack discharge
import unittest
import numpy as np
#import math
from coord_trans import smallz_to_chi 
from coord_trans import chi_to_smallz 
from coord_trans import smallz_to_Z 
from coord_trans import Z_to_smallz 
from Crack import Crack
import uniform_flow
import math
import cmath
import integration


class test_coord_trans(unittest.TestCase):
    
    def setUp(self):
        self.z1 = complex(20, 20)
        self.z2 = complex(20, 40)
        self.z_test1 = complex(20, 30)
        self.z_test2 = complex(50, 50)
        
    def test_smallz_to_Z_and_back(self):
        self.assertEqual(smallz_to_Z(self.z1, self.z1, self.z2), -1)
        self.assertEqual(smallz_to_Z(self.z2, self.z1, self.z2), 1)
        self.assertEqual(Z_to_smallz(-1, self.z1, self.z2), self.z1)
        self.assertEqual(Z_to_smallz(1, self.z1, self.z2) , self.z2)
        pass
    
    def test_chi_to_smallz_and_back(self):
        self.assertEqual(smallz_to_chi(self.z1, self.z1, self.z2), -1)
        self.assertEqual(smallz_to_chi(self.z2, self.z1, self.z2), 1)
        self.assertEqual(chi_to_smallz(1, self.z1, self.z2), self.z2)
        self.assertEqual(chi_to_smallz(-1, self.z1, self.z2), self.z1)
        pass
    
    
    
class test_crack(unittest.TestCase):
    def setUp(self):
        self.z1 = complex(20.0, 20.0)
        self.z2 = complex(20.0, 40.0)
        self.n = 10
        self.coeffs = np.ones(self.n)
        self.real_k_int = 10
        self.aperture = 1
        self.elem_crack = Crack(self.z1, self.z2, self.n,self.real_k_int,self.aperture)
        pass
        
    def test_get_potential(self):
        self.elem_crack.coeffs[0] = 1
        self.assertAlmostEqual(self.elem_crack.get_potential(self.z1), 1,5)
        self.assertAlmostEqual(self.elem_crack.get_potential(self.z2), 1,5)
        
        #checked with Maple
        #F:=(Z-sqrt(Z-1)*sqrt(Z+1))^n;
        #sum(subs(Z=(0.5 - 0.5*((20+I*20)+(20+I*40)))/(0.5*((20+I*40)-(20+I*20))),FF),n=0..9);
        self.elem_crack.coeffs = np.ones(self.n)
        self.assertAlmostEqual(self.elem_crack.get_potential(0.5),complex(.8911069956,-0.6294467672e-1),5)
        pass
    
    def test_get_discharge(self):
        self.elem_crack.coeffs[0] = 1
        #self.assertAlmostEqual(self.elem_crack.get_discharge(self.z1),0,5)
        #self.assertAlmostEqual(self.elem_crack.get_discharge(self.z2),0,5)
        
        #checked with Maple
        #F:=(n*2/((20+I*40)-(20+I*20)))*((Z-sqrt(Z+1)*sqrt(Z-1))^n)/(sqrt(Z-1)*sqrt(Z+1));
        #sum(subs(Z=(0.5 - 0.5*((20+I*20)+(20+I*40)))/(0.5*((20+I*40)-(20+I*20))) ,F),n=0..9);
        self.elem_crack.coeffs = np.ones(self.elem_crack.n)
        self.assertAlmostEqual(self.elem_crack.get_discharge(0.5),complex(0.2790608739e-2,-0.1539150076e-2))
        pass
        #needs more testing
        
class test_uniform_flow(unittest.TestCase):
    def setUp(self):
        self.Q0 = 0.1
        self.beta = 30.0
        self.u_flow = uniform_flow.UniformFlowData(self.Q0,self.beta)
        
    def test_get_potential(self):
        self.assertAlmostEqual(self.u_flow.get_potential(0),0.0,5)
        self.assertAlmostEqual(self.u_flow.get_potential(1.0),-self.Q0*cmath.exp(-1j*math.radians(self.beta)),5)
        pass
        
    def test_get_discharge(self):
        self.assertAlmostEqual(self.u_flow.get_discharge(0), self.Q0*cmath.exp(-1j*math.radians(self.beta)))
        self.assertAlmostEqual(self.u_flow.get_discharge(math.radians(90.0)), self.Q0*cmath.exp(1j*(math.radians(90.0) - math.radians(self.beta))))
        pass
    
class test_integration(unittest.TestCase):
    def setUp(self):
        self.integrator = integration.Integration(0.0,math.pi,450)
        self.n = 1
    def functional(self,theta,n,b,c):
        return math.sin(n*theta)*math.sin(theta)
        

    def test_integrals(self):
        self.assertAlmostEqual(self.integrator.rectangle_rule(self.functional, self.n,0,0), self.integrator.riemman_rule(self.functional, self.n,0,0))
        pass
    def test_integrals_2(self):
        self.assertAlmostEqual(self.integrator.rectangle_rule(self.functional, self.n,0,0), self.integrator.trapezoidal_rule(self.functional, self.n,0,0))
        pass
#    def test_integrals_3(self):
#        self.assertAlmostEqual(self.integrator.trapezoidal_rule(self.functional, self.n,0,0), self.integrator.riemman_rule(self.functional, self.n,0,0))
#        pass
#chi_test = test_coord_trans('test_chi_to_smallz_and_back')
#crack_pot_test = test_crack('test_get_potential')

