'''
Created on 07/03/2011

@author: ispmarin
'''

import math
import time
import numpy as np
import scipy.linalg
#import scipy.integrate


import aux_functions
import integration
#import IterMatrixSolver

class LSMSolver():
    '''
    Class to solve the crack problem.
    '''
    
    def __init__(self, sn, sj, sfold, smax_iterations):
        self.n = sn
        self.j = sj
        self.fold = sfold
        self.max_iterations = smax_iterations
        self.divisions = self.n * 1500
        
        try:
            self.f = open('coeffs_crack.dat', 'w')
        except IOError:
            print ("File does not exist")
        

        print ('Solver', 'n', self.n, 'max_iterations', self.max_iterations)
        
        
    def solve_system(self, u_flow, line_double_list, aq):
        '''
        Solves the crack system using the matrix.
        '''
        u_flow.constant = aux_functions.head_to_pot(aq.reference_head) - \
            (aux_functions.potential_value_no_c(aq.reference_point, line_double_list, 
                                           u_flow)).real
                                           
        #integrator = integration.Integration(0.0, math.pi, self.divisions)
        n_size = len(line_double_list) * self.n + 1
        m_size = len(line_double_list) * self.n + 1
        
        A = np.zeros((n_size, m_size))
        b = np.zeros(n_size)
                     
        matrix_start = time.clock()
        for i, element in enumerate(line_double_list):
            for j, being_calculated in enumerate(line_double_list):
                for n in xrange(self.n):
                    for m in xrange(self.n):
                
                        if n == m and i == j:
                            
                            delta_k =  self.n * (being_calculated.real_k_int + aq.k_ext) / (being_calculated.real_k_int - aq.k_ext)
                            
                            if n > 0 and m > 0:
                            
                                delta_k = self.n * (being_calculated.real_k_int + aq.k_ext) / (being_calculated.real_k_int - aq.k_ext)
                        else:
                            
                            delta_k = 0
                    
                        theta = scipy.pi * n /self.n
                        A[i * self.n + n , j * self.n + m] = -delta_k +  integration.get_one_elem_pot_ex(theta, element, being_calculated, n, m, 1, 1)
        
        
        for i, being_calculated in enumerate(line_double_list):
            for n in xrange(self.n):
                theta = scipy.pi * n /self.n
                A[n_size - 1, i*self.n + n] =   integration.constant_cos(theta, 1, n, being_calculated, being_calculated, aq) 
                A[i * self.n + n, m_size - 1] =  being_calculated.get_potential_m(aq.reference_point,n).real
                
        A[n_size -1 , m_size - 1] = 1
        
        print A
        matrix_stop = time.clock()
        
        for i, being_calculated in enumerate(line_double_list):
            for n in xrange(self.n):
                theta = scipy.pi* n/self.n
                b[i * self.n + n ] =  integration.func_u_flow_no_c_elem(theta, n, 1, self.fold, u_flow, being_calculated)
                

        b[m_size - 1] = aux_functions.head_to_pot(aq.reference_head)- \
            (aux_functions.potential_value_no_c(aq.reference_point, line_double_list, 
                                           u_flow)).real
                
        print 'b', b
        solver_start = time.clock()
        
        x = scipy.linalg.solve(A, b)
        solver_stop = time.clock()
        
        #print 'x', x
        for i,element in enumerate(line_double_list):
            for n in xrange(self.n):
                element.coeffs[n] = x[i*self.n + n]
        
            print 'line double',i,  element.coeffs
        
        u_flow.constant = x[len(line_double_list) * self.n]
        
#        u_flow.constant = aux_functions.head_to_pot(aq.reference_head) - \
#            (aux_functions.potential_value_no_c(aq.reference_point, line_double_list, 
#                                           u_flow)).real


        
        print("Final Constant: ", u_flow.constant)
        print( "Matrix Construction Time:", matrix_stop - matrix_start)
        print( "Solver Time:", solver_stop - solver_start) 
        
        
        #IterMatrixSolver.Jacobi(A, b,  1e-10, 1e-10, self.max_iterations)
        #IterMatrixSolver.GaussSeidel(A, b,  1e-14, 1e-14, self.max_iterations)
        
        