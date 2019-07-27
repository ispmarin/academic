'''
Created on 07/03/2011

@author: ispmarin
'''


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
    
    def __init__(self, sn, sfold, smax_iterations):
        self.n = sn

        self.fold = sfold
        self.max_iterations = smax_iterations
        self.divisions = self.n * 1500
        
        try:
            self.f = open('coeffs_crack.dat', 'w')
        except IOError:
            print ("File does not exist")
        

        print ('Solver', 'n', self.n, 'max_iterations', self.max_iterations)
        
        
    def solve_system(self, u_flow, circle_list, aq):
        '''
        Solves the crack system using the matrix.
        '''
        u_flow.constant = aux_functions.head_to_pot(aq.reference_head) - \
            (aux_functions.potential_value_no_c(aq.reference_point, circle_list, 
                                           u_flow)).real
                                           
        #integrator = integration.Integration(0.0, math.pi, self.divisions)
        n_size = len(circle_list) * self.n 
        m_size = len(circle_list) * self.n 
        
        A = np.zeros((n_size, m_size),dtype=complex)
        b = np.zeros(n_size,dtype=complex)
        coeff =0
        matrix_start = time.clock()
        for i, element in enumerate(circle_list):
            for j, being_calculated in enumerate(circle_list):
                for n in xrange(self.n):
                    for m in xrange(self.n):
                
                        if n == m and i == j:
                            
                            delta_k =   2 *  aq.k_ext/(being_calculated.real_k_int - aq.k_ext)
                            
                            if n > 0 and m > 0:
                            
                                delta_k =   (being_calculated.real_k_int + aq.k_ext) / (being_calculated.real_k_int - aq.k_ext)
                        else:
                            
                            delta_k = 0
                            
                        for k in xrange(30):
                            theta = 2 * scipy.pi * k /self.n
                            coeff += (1.0/self.n) *integration.get_one_elem_pot_ex(theta, being_calculated, element, n, m, 1, 1)
                            
                        A[i * self.n + n , j * self.n + m] =  -delta_k + coeff
                        coeff = 0  
        
        
#        for i, being_calculated in enumerate(circle_list):
#            for n in xrange(self.n):
#                theta = 2 * scipy.pi * n /self.n
#                A[n_size - 1, i*self.n + n] =  (1.0/self.n ) * integration.constant_exp(theta, 1, n, being_calculated, being_calculated, aq) 
#                A[i * self.n + n, m_size - 1] =  being_calculated.get_potential_m(aq.reference_point,n).real
#                
#        A[n_size -1 , m_size - 1] = -1
        
        print A
        matrix_stop = time.clock()
        coeff = 0
        
        for i, being_calculated in enumerate(circle_list):
            for n in xrange(self.n):
                
                for k in xrange(30):
                    theta = 2 * scipy.pi * k /self.n
                    coeff += integration.func_u_flow_circle_no_c(theta, n,  u_flow, being_calculated, being_calculated,aq)
                print coeff
                b[i * self.n + n ] = coeff 
                coeff = 0

#        b[m_size - 1] = -aux_functions.head_to_pot(aq.reference_head)+ \
#            (aux_functions.potential_value_no_c(aq.reference_point, circle_list, 
#                                           u_flow)).real
                
        print 'b', b
        solver_start = time.clock()
        
        x = scipy.linalg.solve(A, b)
        solver_stop = time.clock()
        
        #print 'x', x
        for i,element in enumerate(circle_list):
            for n in xrange(self.n):
                element.coeffs[n] = x[i*self.n + n]
        
            print 'circle',i,  element.coeffs
        
        u_flow.constant = x[len(circle_list) * self.n]
        
#        u_flow.constant = aux_functions.head_to_pot(aq.reference_head) - \
#            (aux_functions.potential_value_no_c(aq.reference_point, circle_list, 
#                                           u_flow)).real


        
        print("Final Constant: ", u_flow.constant)
        print( "Matrix Construction Time:", matrix_stop - matrix_start)
        print( "Solver Time:", solver_stop - solver_start) 
        
        
        #IterMatrixSolver.Jacobi(A, b,  1e-10, 1e-10, self.max_iterations)
        #IterMatrixSolver.GaussSeidel(A, b,  1e-14, 1e-14, self.max_iterations)
        
        