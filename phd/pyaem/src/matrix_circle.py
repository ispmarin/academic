'''
Created on 07/03/2011

@author: ispmarin
'''

import math
import time
import numpy as np
import scipy.linalg
#import scipy.sparse.linalg
#import scipy.integrate
#import pprint

import aux_functions
import integration
import IterMatrixSolver


class MatrixSolver():
    '''
    Class to solve the crack problem.
    '''
    
    def __init__(self, sn):
        self.n = sn
        self.max_iterations = 2000
        self.divisions = self.n * 1500
        
        try:
            self.f = open('coeffs_crack.dat', 'w')
        except IOError:
            print ("File does not exist")
        

        print ('Solver', 'n', self.n)
        
        
    def solve_system(self, u_flow, circle_list, aq):
        '''
        Solves the crack system using the matrix.
        '''
#        u_flow.constant = aux_functions.head_to_pot(aq.reference_head) - \
#            (aux_functions.final_potential(aq.reference_point, circle_list, 
#                                           u_flow)).real
        
        print('Initial Constant',u_flow.constant)
        
        #integrator = integration.Integration(0.0, 2*math.pi, self.divisions)
              
        A = scipy.zeros((len(circle_list) * self.n + 1, len(circle_list) * self.n + 1), dtype = complex)
        b = scipy.zeros(len(circle_list) * self.n + 1, dtype = complex)
                     
        matrix_start = time.clock()
        for i, being_calculated in enumerate(circle_list):
            for j, element in enumerate(circle_list):
                for n in xrange( self.n ):
                    for m in xrange( self.n ):
                        if i == j and n  == m :
                            delta_k =  2 * math.pi * aq.k_ext/(being_calculated.real_k_int - aq.k_ext)
                              
                            if n > 0 and m > 0:
                                delta_k =  math.pi*(being_calculated.real_k_int + aq.k_ext) / (being_calculated.real_k_int - aq.k_ext)
                                 
                        else:
                            delta_k = 0

                        A[i * self.n + n  ,j * self.n + m ] = -delta_k +  integration.complex_quadrature(integration.func_potential_circle_m, 0.0, 2 * scipy.pi, m, n, element, being_calculated ,aq) 
                        #print i,j,m,n, integration.complex_quadrature(integration.func_potential_circle_m, 0.0, 2 * scipy.pi, m, n, element, being_calculated ,aq)
        
        
        for j, being_calculated in enumerate(circle_list):
            for n in xrange(self.n):
                A[len(circle_list)*self.n  ,j * self.n + n ] = being_calculated.get_potential_m(aq.reference_point,n)
                A[j * self.n + n, len(circle_list) *self.n ] = integration.complex_quadrature(integration.constant_exp, 0.0, 2 * scipy.pi, n, n, being_calculated, being_calculated ,aq)
          
        A[ len(circle_list)*self.n, len(circle_list)*self.n ] = -1
        
        print A
        
        matrix_stop = time.clock()
        
        
        for i, being_calculated in enumerate(circle_list):
            for n in xrange(self.n):
                b[i * self.n + n ] = -integration.complex_quadrature(integration.func_u_flow_circle_no_c, 0.0, 2 * scipy.pi,  n, u_flow, being_calculated, circle_list, aq )
                #integration.complex_quadrature(integration.func_u_flow_circle, 0.0, 2 * scipy.pi,  n, u_flow, being_calculated, circle_list, aq )
                
                # -scipy.integrate.quad(integration.func_u_flow_circle, 0.0, 2 * scipy.pi, args=( n, u_flow, being_calculated, circle_list, aq ))[0]
        b[len(circle_list)*self.n] = -aux_functions.head_to_pot(aq.reference_head,aq) + \
            (aux_functions.potential_value_no_c(aq.reference_point, circle_list, 
                                           u_flow)).real
        
        print 'b', b 
        solver_start = time.clock()
        
        x = scipy.linalg.solve(A, b)
        
        print('Residual', scipy.linalg.norm(np.dot(A, x) - b)/scipy.linalg.norm(A))
        
        solver_stop = time.clock()
        
        #print 'x', x
        
        for i,element in enumerate(circle_list):
            for n in xrange(self.n):
                element.coeffs[n] = x[i*self.n + n]
            

        u_flow.constant = x[len(circle_list)*self.n].real
#        aux_functions.head_to_pot(aq.reference_head) - \
#            (aux_functions.potential_value_no_c(aq.reference_point, circle_list, 
#                                           u_flow)).real
#        print 'diff',x[len(circle_list)*self.n] #- u_flow.constant 
                                          
        # this is the analytic solution that I set it up after calculating the constant. This is why the thing before was giving wrong results ONLY
        # for a0, and not the rest - a0 can be calculated analytically, and the rest, numerically. Note the values of the constant before and after
        # the calculation for the coeficients.                  
        for i,element in enumerate(circle_list):
            delta_k_0 = (element.real_k_int - aq.k_ext)/aq.k_ext
            element.coeffs[0] = delta_k_0 * (aux_functions.final_potential(element.center, circle_list,u_flow) - element.get_potential(element.center))
            print 'circle',i, element.coeffs
        
        
        print("Final Constant: ", u_flow.constant)
        print( "Matrix Construction Time:", matrix_stop - matrix_start)
        print( "Solver Time:", solver_stop - solver_start) 
        
        #print IterMatrixSolver.Jacobi(A, b, 1e-10, 1e-10, self.max_iterations)
        #print IterMatrixSolver.GaussSeidel(A, b, 1e-10,1e-10, self.max_iterations)
        