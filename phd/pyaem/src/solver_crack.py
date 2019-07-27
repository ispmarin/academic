'''
Created on 07/03/2011

@author: ispmarin
'''

import math
import time
#import scipy.integrate

import convergence
import aux_functions
import integration


class Solver():
    '''
    Class to solve the crack problem.
    '''
    
    def __init__(self, sn, smax_iterations):
        self.n = sn
        self.max_iterations = smax_iterations
        self.divisions = self.n * 500
        
        try:
            self.f = open('coeffs_crack.dat', 'w')
        except IOError:
            print "File does not exist"
        

        print 'Solver', 'n', self.n, 'max_iterations', self.max_iterations
        
        
    def solve_system(self, u_flow, crack_list, aq):
        '''
        Solves the crack system using the iterative scheme.
        '''
        
        iteration = 0
        conv_test = convergence.convergence_element(self.n, len(crack_list), 1e-6, 1e-7,crack_list)
       
        integrator = integration.Integration(0.0, math.pi, self.divisions)
        
        u_flow.constant = aux_functions.head_to_pot(aq.reference_head) - \
            (aux_functions.potential_value_no_c(aq.reference_point, crack_list, 
                                           u_flow)).real
        
        for being_calculated in crack_list:
            being_calculated.coeffs[1] = 0.00001 
#            (being_calculated.real_k_int * \
#            being_calculated.aperture / (aq.k_ext * being_calculated.L \
#            + being_calculated.real_k_int * being_calculated.aperture)) * \
#            u_flow.get_discharge(being_calculated.angle).real * \
#                                            being_calculated.L * 0.5
#                
        
        while iteration < self.max_iterations:
            
            iter_time = time.clock()
                
            self.f.write('%d      ' % iteration)



            for being_calculated in crack_list:
                
                delta_k = (1.0 / math.pi) * (being_calculated.real_k_int * \
                                             being_calculated.aperture / aq.k_ext)
                

                        
                being_calculated.coeffs[1]  = -delta_k * integrator.rectangle_rule(integration.func_discharge, 1, crack_list, being_calculated) 
                #scipy.integrate.quad(integration.func_discharge, 0.0, math.pi, args=( 1, crack_list, being_calculated), limit=100)[0]
                print 'coeff ',1,being_calculated.coeffs[1], integrator.rectangle_rule(integration.func_discharge, 1, crack_list, being_calculated)
                for n in range(2, self.n):
                
                    being_calculated.coeffs[n] = -delta_k * integrator.rectangle_rule(integration.func_discharge, n, crack_list, being_calculated) 
                    #scipy.integrate.quad(integration.func_discharge, 0.0, math.pi, args=( n, crack_list, being_calculated), limit=100)[0]
                    print 'coeff ',n,being_calculated.coeffs[n],integrator.rectangle_rule(integration.func_discharge, n, crack_list, being_calculated)
                    
            if conv_test.check_convergence(crack_list) == 1:
                break
                
            iteration = iteration + 1
            
                 
            self.f.write('\n')
            print iteration, u_flow.constant, 'iteration time: ', \
            time.clock() - iter_time
            
        u_flow.constant = aux_functions.head_to_pot(aq.reference_head) - \
            (aux_functions.potential_value_no_c(aq.reference_point, crack_list, 
                                           u_flow)).real

        for i, element in enumerate(crack_list):
            print 'crack ', i, element.coeffs