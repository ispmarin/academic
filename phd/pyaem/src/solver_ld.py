'''
Created on 07/03/2011

@author: ispmarin
'''


import math
import time

import convergence
import aux_functions
#import output
import integration



class Solver():
    
    def __init__(self, sn, sj, smax_iterations):
        self.n = sn
        self.j = sj
        self.max_iterations = smax_iterations
        self.divisions = self.n * 10
                
        print( 'Solver', 'n', self.n, 'j', self.j, 'max_iterations', 
        self.max_iterations)
        
        
        
    def solve_system(self, u_flow, line_double_list, aq):
        
        try:
            f = open('coeffs_ld.dat', 'w')
        except IOError:
            print "File does not exist"
        
        
        iteration = 0
        u_flow.constant = aux_functions.head_to_pot(aq.reference_head,aq) - \
        (aux_functions.potential_value_no_c(aq.reference_point, line_double_list, 
                                       u_flow).real)
        
        print 'initial constant', u_flow.constant
        
        conv_test = convergence.convergence_element(self.n, len(line_double_list), 1e-7, 1e-7, line_double_list)
        
        integrator = integration.Integration(0.0, math.pi, self.divisions)
        total_iter_start = time.clock()
        while iteration < self.max_iterations:
            
            iter_time = time.clock()
                
            f.write('%d      ' % iteration)
                        
            for  being_calculated in line_double_list:
              
                
                delta_k = (2.0 / math.pi) * ((being_calculated.real_k_int - aq.k_ext)\
                / (aq.k_ext + being_calculated.real_k_int)) #correct
                
                
                being_calculated.coeffs[0] = delta_k * \
                integrator.rectangle_rule(integration.func_potential, 0,line_double_list,being_calculated,u_flow)
                #scipy.integrate.quad(integration.func_potential, 0.0+1e-3, math.pi-1e-3, args=( 0, line_double_list, being_calculated, u_flow),limit=100)[0]
              

                for n in xrange(1, self.n):

                    being_calculated.coeffs[n] = 2.0 * delta_k * \
                    integrator.rectangle_rule(integration.func_potential, n,line_double_list,being_calculated,u_flow)
                    #scipy.integrate.quad(integration.func_potential, 0.0+1e-3, math.pi-1e-3, args=( n, line_double_list, being_calculated, u_flow),limit=100)[0]

                                                   
                    f.write('%f     ' % (being_calculated.coeffs[n]))
                        
                 
                #print being_calculated.coeffs
            if conv_test.check_convergence(line_double_list) == 1:
                break
            
            iteration = iteration + 1
            f.write('\n')
            print iteration, u_flow.constant, 'iteration time: ', time.clock() - iter_time
       
            u_flow.constant = aux_functions.head_to_pot(aq.reference_head,aq) - \
            (aux_functions.potential_value_no_c(aq.reference_point, line_double_list, 
                                       u_flow).real)
        
        total_iter_stop = time.clock()
        print 'Final constant', u_flow.constant
        
        for i, element in enumerate(line_double_list):
            print 'line double', i, element.coeffs
            

       
        print("Final Constant: ", u_flow.constant)
        print( "Solver Time:", total_iter_stop - total_iter_start)
        

