'''
Created on 07/03/2011

@author: ispmarin
'''

import math
import time
import numpy as np
import scipy.linalg
import scipy.integrate

import aux_functions
import integration
import IterMatrixSolver


class MatrixSolver():
    '''
    Class to solve the crack problem.
    '''
    
    def __init__(self, sn, smax_iterations):
        self.n = sn
        self.max_iterations = smax_iterations
        self.divisions = self.n * 1500
        
        try:
            self.f = open('coeffs_crack.dat', 'w')
        except IOError:
            print ("File does not exist")
        

        print ('Solver', 'n', self.n, 'max_iterations', self.max_iterations)
        
        
    def solve_system(self, u_flow, crack_list, aq):
        '''
        Solves the crack system using the matrix.
        '''

        #integrator = integration.Integration(0.0, math.pi, self.divisions)
#        
        
        A = np.zeros((len(crack_list) * self.n + 1, len(crack_list) * self.n + 1))
        
        b = np.zeros(len(crack_list) * self.n + 1)
                     
        matrix_start = time.clock()
        for i, being_calculated in enumerate(crack_list):
            for j, element in enumerate(crack_list):
                for n in xrange(1, self.n +1):
                    for m in xrange(1, self.n +1):
                        if i == j and n - 1 == m - 1:
                            delta_k =  math.pi * aq.k_ext / (being_calculated.real_k_int * being_calculated.aperture)
                        else:
                            delta_k = 0

                        A[i * self.n + n - 1 ,j * self.n + m - 1] = delta_k - scipy.integrate.quad(integration.func_discharge_m, 0.0, math.pi, args=( m, n, element, being_calculated),limit=100)[0]
                        
        for i, being_calculated in enumerate(crack_list):
            for n in xrange(self.n): 
                A[len(crack_list) * self.n, i*self.n + n  ] = being_calculated.get_potential_m(aq.reference_point,n).real
                A[i * self.n + n, len(crack_list) * self.n ] = 0         
        
        A[len(crack_list) * self.n,len(crack_list) * self.n] = 1
        
        
        
        matrix_stop = time.clock()
        
        for i, being_calculated in enumerate(crack_list):
            for n in xrange(1,self.n):
                b[i * self.n + n - 1] = scipy.integrate.quad(integration.func_discharge_u_flow, 0, math.pi, args=( n, u_flow, being_calculated))[0]
        
        b[len(crack_list) * self.n] = aux_functions.head_to_pot(aq.reference_head,aq)- \
            (aux_functions.potential_value_no_c(aq.reference_point, crack_list, 
                                           u_flow)).real
        print'VALUE', (aux_functions.potential_value_no_c(aq.reference_point, crack_list, 
                                           u_flow)).real
        #print 'b', b 
        solver_start = time.clock()
        
        x = scipy.linalg.solve(A, b)
        
        solver_stop = time.clock()
        
        #print 'x', x
        for i,element in enumerate(crack_list):
            for n in xrange(1,self.n):
                element.coeffs[n] = (x[i*self.n + (n-1)])
        
        for i, element in enumerate(crack_list):
            print i, element.coeffs
        
        u_flow.constant = aux_functions.head_to_pot(aq.reference_head,aq) - \
            (aux_functions.potential_value_no_c(aq.reference_point, crack_list, 
                                           u_flow)).real
        #u_flow.constant = x[len(crack_list)*self.n]
        
        being_calculated = crack_list[0]
        
        analytic_coeff = (being_calculated.real_k_int * \
                being_calculated.aperture / (aq.k_ext * being_calculated.L \
                + being_calculated.real_k_int * being_calculated.aperture)) * \
                u_flow.get_discharge(being_calculated.angle).real * \
                                                being_calculated.L * 0.5
        
        print('u_flow constant', u_flow.constant)
                                                
        print( "analytic", analytic_coeff)
        print( "Matrix Construction Time:", matrix_stop - matrix_start)
        print( "Solver Time:", solver_stop - solver_start) 
        
        print('Solving Crack Using Jacobi Solver')        
        IterMatrixSolver.Jacobi(A, b,  1e-10, 1e-10, self.max_iterations)
        print('Done.\n')
        print('Solving Crack Using Gauss Seidel Solver')   
        IterMatrixSolver.GaussSeidel(A, b, 1e-10,1e-10, self.max_iterations)
        print('Done.\n')
        
