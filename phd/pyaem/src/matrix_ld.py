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
    
    def __init__(self, sn, sj, smax_iterations):
        self.n = sn
        self.j = sj
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
        u_flow.constant = aux_functions.head_to_pot(aq.reference_head,aq) - \
            (aux_functions.potential_value_no_c(aq.reference_point, line_double_list, 
                                           u_flow)).real
                                           
        #integrator = integration.Integration(0.0, math.pi, self.divisions)
        
        
        
        
        
        A = np.zeros((len(line_double_list) * self.n +1, len(line_double_list) * self.n +1) )
        b = np.zeros(len(line_double_list) * self.n + 1)
                     
        matrix_start = time.clock()
        for i, being_calculated in enumerate(line_double_list):
            for j, element in enumerate(line_double_list):
                for m in xrange( self.n ):
                    for n in xrange( self.n ):
                        if i == j and n  == m :
                            delta_k = (math.pi/2.0) * ((being_calculated.real_k_int + aq.k_ext) / (being_calculated.real_k_int - aq.k_ext))
                            if n > 0 and m > 0:
                                delta_k = (math.pi/4.0) * ((being_calculated.real_k_int + aq.k_ext) / (being_calculated.real_k_int - aq.k_ext))
                        else:
                            delta_k = 0
                        
                        A[i * self.n + n  ,j * self.n + m ] = -delta_k + scipy.integrate.quad(integration.func_potential_m, 0.0, math.pi, args=( m, n, element, being_calculated),limit=200)[0]
                        #mpmath.fp.quad(lambda x: integration.func_potential_m(x, m, n,element,being_calculated), [0.0, scipy.pi], method='tanh-sinh')
                        #scipy.integrate.quad(integration.func_potential_m, 0.0, math.pi, args=( m, n, element, being_calculated),limit=200)[0]
                
                        
        for i, being_calculated in enumerate(line_double_list):
            for n in xrange(self.n): 
                A[len(line_double_list) * self.n, i*self.n + n  ] = being_calculated.get_potential_m(aq.reference_point,n).real
                A[i * self.n + n, len(line_double_list) * self.n ] = scipy.integrate.quad(integration.constant_cos, 0.0, scipy.pi, args=(n, n, being_calculated, being_calculated ,aq),limit=200)[0]
                #mpmath.fp.mpf(mpmath.quad(lambda x: integration.constant_cos(x,n,n,being_calculated,being_calculated,aq), [0.0, scipy.pi], method='tanh-sinh')) 
                #scipy.integrate.quad(integration.constant_cos, 0.0, scipy.pi, args=(n, n, being_calculated, being_calculated ,aq),limit=200)[0] 
        
        
        A[len(line_double_list) * self.n,len(line_double_list) * self.n] = 1
        
        
        matrix_stop = time.clock()
        
        for i, being_calculated in enumerate(line_double_list):
            for n in xrange(self.n):
                b[i * self.n + n ] = -scipy.integrate.quad(integration.func_u_flow_no_c, 0.0, math.pi, args= (n, u_flow, being_calculated))[0]
                #-mpmath.fp.quad(lambda x: integration.func_u_flow_no_c(x,n,u_flow,being_calculated), [0.0, scipy.pi], method='tanh-sinh') 
                #-scipy.integrate.quad(integration.func_u_flow_no_c, 0.0, math.pi, args= (n, u_flow, being_calculated))[0]


        b[len(line_double_list) * self.n] = aux_functions.head_to_pot(aq.reference_head,aq)- \
            (aux_functions.potential_value_no_c(aq.reference_point, line_double_list, 
                                           u_flow)).real
                
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
        
        
        result, iter_j = IterMatrixSolver.Jacobi(A, b,  1e-7, 1e-7, self.max_iterations)
        result_g, iter_gs = IterMatrixSolver.GaussSeidel(A, b,  1e-7, 1e-7, self.max_iterations)
        print'iteracoes', iter_j, iter_gs
        