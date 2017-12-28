'''
Created on 07/03/2011

@author: ispmarin
'''

#import scipy.integrate
import math
import time

import convergence
import aux_functions
#import output
import integration



class Solver():
    
    def __init__(self, sn, smax_iterations):
        self.n = sn
        self.max_iterations = smax_iterations
        self.divisions = self.n * 10
                
        print( 'Solver', 'n', self.n, 'max_iterations', 
        self.max_iterations)
        
        
        
    def solve_system(self, u_flow, ellipse_list, aq):
        
        try:
            f = open('coeffs_circle.dat', 'w')
        except IOError:
            print "File does not exist"
        
        
        iteration = 0
        u_flow.constant = aux_functions.head_to_pot(aq.reference_head) - \
        (aux_functions.potential_value_no_c(aq.reference_point, ellipse_list, 
                                       u_flow).real)
        
        print 'initial constant', u_flow.constant
        
        conv_test = convergence.convergence(self.n, len(ellipse_list), 1e-5)
        
        #integrator = integration.Integration(0.0, 2*math.pi, self.divisions)
               
        while iteration < self.max_iterations:
            
            iter_time = time.clock()
                
            f.write('%d      ' % iteration)
                        
            for  being_calculated in ellipse_list:
                
                delta_k_0 = (being_calculated.real_k_int - aq.k_ext)/aq.k_ext
                delta_k = (1/math.pi) * ((being_calculated.real_k_int - aq.k_ext)/(aq.k_ext + being_calculated.real_k_int))
                
                
                being_calculated.coeffs[0] = delta_k_0 * (aux_functions.final_potential(being_calculated.center, ellipse_list,u_flow) -
                                                           being_calculated.get_potential(being_calculated.center))
                
                print aux_functions.final_potential(being_calculated.center, ellipse_list,u_flow), being_calculated.get_potential(being_calculated.center)
                
                for n in xrange(1, self.n):
                   
                    being_calculated.coeffs[n] = delta_k * \
                    1#integrator.rectangle_rule(integration.func_potential_ellipse, n, ellipse_list,being_calculated,u_flow)
                    #scipy.integrate.quadrature(integration.func_potential_circle, 0.0, 2*math.pi, args=( n, ellipse_list, being_calculated, u_flow), maxiter=100)[0]
                                                   
            u_flow.constant = aux_functions.head_to_pot(aq.reference_head) - \
            (aux_functions.potential_value_no_c(aq.reference_point, ellipse_list, 
                                       u_flow).real)
            
            if conv_test.check_convergence(ellipse_list) == 1:
                break
            
            iteration = iteration + 1
            f.write('\n')
            print iteration, u_flow.constant, 'iteration time: ', time.clock() - iter_time
       
            
        
        print 'Final constant', u_flow.constant
        
        for i, element in enumerate(ellipse_list):
            print 'ellipse', i, element.coeffs
            

        #aux_functions.check_ld_jump(line_double_list, u_flow,  aq)

