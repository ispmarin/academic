#! /usr/bin/env python2.6
# -*- coding: utf-8 -*-

'''
Created on 07/03/2011

@author: ispmarin
'''
import file_input
import time
import sys
import numpy as np
from pprint import pprint


import aquifer
import uniform_flow
import solver_crack
import solver_ld
import solver_circle
import solver_ellipse
import output
import aux_functions
import matrix_crack
import matrix_ld
import matrix_circle
import LineDouble
import CircleInhom
import lsm_ld
#import lsm_circle
#import ErrorAnalysis
#import EllipseInhom

def set_u_flow_aq():
    #print('Setting up the system\n')
    
    try:
        open('general_input.dat')
        data = np.loadtxt('general_input.dat')
        
    except IOError:
        sys.exit('ERROR Cannot open general_input.dat')
    
    aquifer_elem = aquifer.AquiferData(data[0],data[1],data[2],data[3],complex(data[4],data[5]))
    u_flow = uniform_flow.UniformFlowData(data[6],data[7])
    
    print('Aquifer and Uniform Flow parameters set\n')
    return u_flow, aquifer_elem

def solve_crack(crack_list,uni_flow, aquifer_elem, kind_solver):
    '''
    Running only the crack
    '''
    start_crack = time.clock()
    try:
        open('input_crack.dat')
        data = np.loadtxt('input_crack.dat')
        
    except IOError:
        sys.exit('ERROR Cannot open input_crack.dat')
    
    if kind_solver == 'iterative':
        print('+++ Solving Crack System with ITERATIVE Solver +++\n')
        crack_solver_handler = solver_crack.Solver(int(data[0]), int(data[1]))
        
    elif kind_solver == 'matrix':
        print('+++ Solving Crack System with MATRIX Solver +++\n')
        crack_solver_handler = matrix_crack.MatrixSolver(int(data[0]), int(data[1]))
        
    else:
        sys.exit('ERROR Solver Strategy not chosen: please use iterative or matrix')
   
    file_input.insert_crack(crack_list, 'crack.dat', crack_solver_handler.n)
        
    #pprint(crack_list)

    crack_solver_handler.solve_system(uni_flow, crack_list, aquifer_elem)
    
    stop_crack = time.clock()
    
    print( 'Crack solver time: ', stop_crack - start_crack,"\n")
    print('Crack Solver Done.')


def solve_circle(circle_list,uni_flow, aquifer_elem, kind_solver):
    '''
    Running only the circle
    '''
    start_circle = time.clock()
    try:
        open('input_crack.dat')
        data = np.loadtxt('input_crack.dat')
        
    except IOError:
        sys.exit('ERROR Cannot open input_crack.dat')

    
    if kind_solver == 'iterative':
        print('\n+++ Solving Circle Inhom System with ITERATIVE Solver +++\n')
        circle_solver_handler = solver_circle.Solver(int(data[0]), int(data[1]))

    elif kind_solver == 'matrix':
        print('\n+++ Solving Circle System with MATRIX Solver +++\n')
        circle_solver_handler = matrix_circle.MatrixSolver(int(data[0]))
        
    else:
        sys.exit('ERROR Solver Strategy not chosen: please use iterative or matrix')
   
    file_input.insert_circle(circle_list, 'circle.dat', circle_solver_handler.n)
        
    #pprint(circle_list)

    circle_solver_handler.solve_system(uni_flow, circle_list, aquifer_elem)
    
    stop_circle = time.clock()
    
    print( 'Circle solver time: ', stop_circle - start_circle,"\n")
    print('Circle Solver Done.')

def solve_dual_circle(circle_list,uni_flow, aquifer_elem, kind_solver):
    '''
    Running only the crack
    '''
    start_circle = time.clock()
    try:
        open('input_crack.dat')
        data = np.loadtxt('input_crack.dat')
        
    except IOError:
        sys.exit('ERROR Cannot open input_crack.dat')

    
   
    circle_solver_handler_iter = solver_circle.Solver(int(data[0]), int(data[1]))
    circle_solver_handler_mat = matrix_circle.MatrixSolver(int(data[0]))
    
    file_input.insert_circle(circle_list, 'circle.dat', circle_solver_handler_iter.n)

    circle_solver_handler_mat.solve_system(uni_flow, circle_list, aquifer_elem)
    circle_solver_handler_iter.solve_system(uni_flow, circle_list, aquifer_elem)
    
    stop_circle = time.clock()
    
    print( 'Circle solver time: ', stop_circle - start_circle,"\n")
  
    
def solve_ellipse(ellipse_list,uni_flow, aquifer_elem, kind_solver):
    '''
    Running only the crack
    '''
    start_circle = time.clock()
    try:
        open('input_crack.dat')
        data = np.loadtxt('input_crack.dat')
        
    except IOError:
        sys.exit('ERROR Cannot open input_crack.dat')

    
    if kind_solver == 'iterative':
        ellipse_solver_handler = solver_ellipse.Solver(int(data[0]), int(data[1]))
        #print('+++ Solving Ellipse Inhom System with ITERATIVE Solver +++\n')
    elif kind_solver == 'matrix':
        sys.exit('ERROR Solver Strategy not implemented yet: please use iterative')
        #crack_solver_handler = matrix_crack.MatrixSolver(int(data[0]), int(data[1]))
        #print('+++ Solving Crack System with MATRIX Solver +++')
    
    else:
        sys.exit('ERROR Solver Strategy not chosen: please use iterative or matrix')
   
    file_input.insert_ellipse(ellipse_list, 'ellipse.dat', ellipse_solver_handler.n)
        
    #pprint(ellipse_list)

    ellipse_solver_handler.solve_system(uni_flow, ellipse_list, aquifer_elem)
    
    stop_circle = time.clock()
    
    print( 'Ellipse solver time: ', stop_circle - start_circle,"\n")


def solve_circle_exact(circle_list, uni_flow,aquifer_elem, exact_vector):
    
    
    exact_vector.append(CircleInhom.ExactSolution(circle_list[0].center, circle_list[0].R, 
                                                  circle_list[0].real_k_int, uni_flow.Q0, aquifer_elem))
    
def solve_circle_ld(circle_list, line_double_list, uni_flow,aquifer_elem, ssubdivisions, kind_solver):
    
    
    
    try:
        open('input_ld.dat')
        data = np.loadtxt('input_ld.dat')
        
    except IOError:
        sys.exit('ERROR Cannot open input_ld.dat')
    
    
    if kind_solver == 'iterative':
        print('\n+++ Solving Line Doublet System with ITERATIVE Solver +++\n')
        solver_ld_handler = solver_ld.Solver(int(data[0]), int(data[1]), 
                                    int(data[2]))

    elif kind_solver == 'matrix':
        print('\n+++ Solving Line Doublet System with MATRIX Solver +++\n')
        solver_ld_handler = matrix_ld.MatrixSolver(int(data[0]), int(data[1]), 
                                        int(data[2]))

    else:
        sys.exit('ERROR Solver Strategy not chosen: please use iterative or matrix')
    print('Solving Exact Circle')
    
    for elem_id, circle in enumerate(circle_list):
        file_input.gen_line_double_ellipse_by_semiaxes(line_double_list, elem_id, solver_ld_handler.n, solver_ld_handler.j, 
                                                   circle.R, circle.R, circle.center, 
                                                   ssubdivisions, circle.real_k_int,1.5)
    
    solver_ld_handler.solve_system(uni_flow, line_double_list, aquifer_elem)
    print('Circle + LD Solver Done.')
    
def solve_circle_exact_ld(circle_list, line_double_list, uni_flow,aquifer_elem, exact_vector, ssubdivisions, kind_solver):
    
    
    
    try:
        open('input_ld.dat')
        data = np.loadtxt('input_ld.dat')
        
    except IOError:
        sys.exit('ERROR Cannot open input_ld.dat')
    
    
    if kind_solver == 'iterative':
        print('\n+++ Solving Line Doublet System with ITERATIVE Solver +++\n')
        solver_ld_handler = solver_ld.Solver(int(data[0]), int(data[1]), 
                                    int(data[2]))

    elif kind_solver == 'matrix':
        print('\n+++ Solving Line Doublet System with MATRIX Solver +++\n')
        solver_ld_handler = matrix_ld.MatrixSolver(int(data[0]), int(data[1]), 
                                        int(data[2]))

    elif kind_solver == 'dual':
        solver_ld_handler_iter = solver_ld.Solver(int(data[0]), int(data[1]), 
                                    int(data[2]))
        solver_ld_handler_mat = matrix_ld.MatrixSolver(int(data[0]), int(data[1]), 
                                        int(data[2]))
    else:
        sys.exit('ERROR Solver Strategy not chosen: please use iterative or matrix')
    
    print('Solving Exact Circle')
    
    exact_vector.append(CircleInhom.ExactSolution(circle_list[0].center, circle_list[0].R, 
                                                  circle_list[0].real_k_int, uni_flow.Q0, aquifer_elem))
    
    if kind_solver == 'dual':
        file_input.gen_line_double_ellipse_by_semiaxes(line_double_list, 0, solver_ld_handler_iter.n, solver_ld_handler_iter.j, 
                                                   circle_list[0].R, circle_list[0].R, circle_list[0].center, 
                                                   ssubdivisions, circle_list[0].real_k_int,1.5)
        
        solver_ld_handler_mat.solve_system(uni_flow, line_double_list, aquifer_elem)
        
        solver_ld_handler_iter.solve_system(uni_flow, line_double_list, aquifer_elem)
    
    else:
        
        file_input.gen_line_double_ellipse_by_semiaxes(line_double_list, 0, solver_ld_handler.n, solver_ld_handler.j, 
                                                   circle_list[0].R, circle_list[0].R, circle_list[0].center, 
                                                   ssubdivisions, circle_list[0].real_k_int,1.5)
    
        solver_ld_handler.solve_system(uni_flow, line_double_list, aquifer_elem)

def solve_ld_crack(crack_list,line_double_list,uni_flow, aquifer_elem, kind_solver):
    '''
    Run the line doublet validation
    '''
    start_line_double = time.clock()
    
    try:
        open('input_ld.dat')
        data = np.loadtxt('input_ld.dat')
        
    except IOError:
        sys.exit('ERROR Cannot open input_ld.dat')
    
    
    if kind_solver == 'iterative':
        print('+++ Solving Line Doublet System with ITERATIVE Solver +++\n')   
        solver_ld_handler = solver_ld.Solver(int(data[0]), int(data[1]), 
                                    int(data[2]))
     
    elif kind_solver == 'matrix':
        print('+++ Solving Line Doublet System with MATRIX Solver +++\n')
        solver_ld_handler = matrix_ld.MatrixSolver(int(data[0]), int(data[1]), 
                                        int(data[2]))

    else:
        sys.exit('ERROR Solver Strategy not chosen: please use iterative or matrix')


    for elem_id, crack in enumerate(crack_list):
        subdivisions = 28      
        far_field = 1.5
        file_input.gen_line_double_ellipse_by_box(line_double_list, elem_id,  
        solver_ld_handler.n, solver_ld_handler.j, crack.z1, crack.z2, 
        crack.aperture, subdivisions, crack.real_k_int, far_field)
    
    
    solver_ld_handler.solve_system(uni_flow, line_double_list, aquifer_elem)
    
    stop_line_double = time.clock()
    
    print( 'Line Doublet solver time: ', stop_line_double - start_line_double,"\n")
    print('Line Doublet Solver Done.')


def solve_ld(input_file, line_double_list,uni_flow, aquifer_elem, file_list, kind_solver):
    '''
    Run the line doublet validation
    '''
    start_line_double = time.clock()
    
    try:
        open(input_file)
        data = np.loadtxt(input_file)
        
    except IOError:
        sys.exit('ERROR Cannot open file')
    
    
    if kind_solver == 'iterative':
        print('+++ Solving Line Doublet System with ITERATIVE Solver +++\n')   
        solver_ld_handler = solver_ld.Solver(int(data[0]), int(data[1]), 
                                    int(data[2]))
     
    elif kind_solver == 'matrix':
        print('+++ Solving Line Doublet System with MATRIX Solver +++\n')
        solver_ld_handler = matrix_ld.MatrixSolver(int(data[0]), int(data[1]), 
                                        int(data[2]))

    else:
        sys.exit('ERROR Solver Strategy not chosen: please use iterative or matrix')

    for elem_id, file in enumerate(file_list):
        
        file_input.insert_line_double(line_double_list, file, elem_id, solver_ld_handler.n, solver_ld_handler.j)
          
      
    solver_ld_handler.solve_system(uni_flow, line_double_list, aquifer_elem)
    
    stop_line_double = time.clock()
    
    print( 'Line Doublet solver time: ', stop_line_double - start_line_double,"\n")
    print('Line Doublet Solver Done.')


def solve_ld_check(line_double_list, uni_flow, aquifer_elem,k_int, kind_solver):
    '''
    Run the line doublet validation
    '''
    start_line_double = time.clock()
    
    try:
        open('input_ld.dat')
        data = np.loadtxt('input_ld.dat')
        
    except IOError:
        sys.exit('ERROR Cannot open input_ld.dat')

    if kind_solver == 'iterative':
        print('+++ Solving Line Doublet System CHECK with ITERATIVE Solver +++\n')
        solver_ld_handler = solver_ld.Solver(int(data[0]), int(data[1]), 
                                    int(data[2]))
        #
    elif kind_solver == 'matrix':
        print('+++ Solving Line Doublet System CHECK with MATRIX Solver +++\n')
        solver_ld_handler = matrix_ld.MatrixSolver(int(data[0]), int(data[1]), 
                                        int(data[2]))
    elif kind_solver == 'lsm':
        solver_ld_handler = lsm_ld.LSMSolver(int(data[0]), int(data[1]), 1, int(data[2]))
    else:
        sys.exit('ERROR Solver Strategy not chosen: please use iterative or matrix')


    line_double_list.append(LineDouble.LineDouble(0,complex(20,20), complex(40,20), solver_ld_handler.n, solver_ld_handler.j, k_int, 1.5 ))
    line_double_list.append(LineDouble.LineDouble(0,complex(40,20), complex(40,40), solver_ld_handler.n, solver_ld_handler.j, k_int, 1.5 ))
    line_double_list.append(LineDouble.LineDouble(0,complex(40,40), complex(20,40), solver_ld_handler.n, solver_ld_handler.j, k_int, 1.5 ))
    line_double_list.append(LineDouble.LineDouble(0,complex(20,40), complex(20,20), solver_ld_handler.n, solver_ld_handler.j, k_int, 1.5 ))
        
    pprint(line_double_list)
    
    solver_ld_handler.solve_system(uni_flow, line_double_list, aquifer_elem)
    
    stop_line_double = time.clock()
    
    print( 'Line Doublet solver time: ', stop_line_double - start_line_double, "\n")

def solve_dual_ld_check(line_double_list, uni_flow, aquifer_elem,k_int, kind_solver):
    '''
    Run the line doublet validation
    '''
    start_line_double = time.clock()
    
    try:
        open('input_ld.dat')
        data = np.loadtxt('input_ld.dat')
        
    except IOError:
        sys.exit('ERROR Cannot open input_ld.dat')

   
    solver_ld_handler_iter = solver_ld.Solver(int(data[0]), int(data[1]), 
                                    int(data[2]))
    solver_ld_handler_mat = matrix_ld.MatrixSolver(int(data[0]), int(data[1]), 
                                        int(data[2]))
    
    line_double_list.append(LineDouble.LineDouble(0,complex(20,20), complex(40,20), solver_ld_handler_iter.n, solver_ld_handler_iter.j, k_int, 1.5 ))
    line_double_list.append(LineDouble.LineDouble(0,complex(40,20), complex(40,40), solver_ld_handler_iter.n, solver_ld_handler_iter.j, k_int, 1.5 ))
    line_double_list.append(LineDouble.LineDouble(0,complex(40,40), complex(20,40), solver_ld_handler_iter.n, solver_ld_handler_iter.j, k_int, 1.5 ))
    line_double_list.append(LineDouble.LineDouble(0,complex(20,40), complex(20,20), solver_ld_handler_iter.n, solver_ld_handler_iter.j, k_int, 1.5 ))
        
    pprint(line_double_list)
    
    solver_ld_handler_mat.solve_system(uni_flow, line_double_list, aquifer_elem)
    solver_ld_handler_iter.solve_system(uni_flow, line_double_list, aquifer_elem)
    
    stop_line_double = time.clock()
    
    print( 'Line Doublet solver time: ', stop_line_double - start_line_double, "\n")
    

def plotting(element_list_1, element_list_2, u_flow, aq, plot_kind, name_plot_1, name_plot_2, kind_solver, save_to_file):
    '''
    Plot both crack and line double
    '''
    rel_tolerance = 1e-7
    abs_tolerance = 1e-7
    start_plot = time.clock()
    
    min_x, min_y, max_x, max_y, num_steps = aux_functions.set_axis(element_list_1)
    
    plot = output.plotter(min_x, min_y, max_x, max_y, num_steps)
    
    if plot_kind == 'head':
        plot.plot_head(element_list_1, u_flow,aq, name_plot_1,kind_solver, save_to_file)
    
    elif plot_kind == 'potential':
        plot.plot_potential(element_list_1, u_flow, name_plot_1,kind_solver, save_to_file)
        
    elif plot_kind == 'compare head':
        plot.plot_potential_matplotlib(element_list_1, element_list_2, u_flow, aq, abs_tolerance,rel_tolerance, 'head',name_plot_1, name_plot_2,kind_solver,save_to_file)
    
    elif plot_kind == 'compare potential':
        plot.plot_potential_matplotlib(element_list_1, element_list_2, u_flow, aq, abs_tolerance,rel_tolerance, 'potential',name_plot_1, name_plot_2, kind_solver,save_to_file)
    
    elif plot_kind == 'check head':
        plot.check_head_line(element_list_1, u_flow, element_list_2[0], element_list_2[1])
    else:
        sys.exit('ERROR Kind of plot not chosen. Please use head, potential, compare head or compare potential')
    
   
    stop_plot = time.clock()
    
    print( 'Plot time: ', stop_plot - start_plot)
  

def crack_alone(UNIFORM_FLOW, AQUIFER, solver_method_crack, solver_method_ld, do_plot):
    
    CRACK_L = []
    LINE_DOUBLE = []
    f = open('timing.dat','a') 
    solve_time = time.clock()
    solve_crack(CRACK_L, UNIFORM_FLOW, AQUIFER, solver_method_crack)
    input_file = 'input_ld.dat'
    #file_list = ['ld.dat','ld2.dat']
    file_list = ['ld.dat']
    #solve_ld(input_file, LINE_DOUBLE, UNIFORM_FLOW, AQUIFER, file_list, solver_method)
    val = 'system ' +solver_method + ' '+ str(time.clock() - solve_time) + '\n'
    f.write(str(val))
    
    if do_plot == 1:
        plot_time = time.clock()
    #plotting(CRACK_L, LINE_DOUBLE,  UNIFORM_FLOW, AQUIFER, 'compare head', 'Fratura', 'Line_Doublet',solver_method, 'yes')
        plotting(CRACK_L, LINE_DOUBLE,  UNIFORM_FLOW,AQUIFER,'head',  'Fratura','nan', solver_method,'yes')
    #plotting(LINE_DOUBLE, LINE_DOUBLE,  UNIFORM_FLOW,AQUIFER,'head', 'Line_Doublet', 'nan',solver_method, 'yes')
        val = 'plot ' +solver_method + ' '+ str(time.clock() - plot_time) + '\n' 
        f.write(str(val))
    
    

def crack_with_ld(UNIFORM_FLOW, AQUIFER, solver_method_crack, solver_method_ld, do_plot):
    
    CRACK_L = []
    LINE_DOUBLE = []
    f = open('timing.dat','a') 
    solve_time = time.clock()
    solve_crack(CRACK_L, UNIFORM_FLOW, AQUIFER, solver_method_crack)
    input_file = 'input_ld.dat'
    #file_list = ['ld.dat','ld2.dat']
    file_list = ['ld.dat']
    solve_ld(input_file, LINE_DOUBLE, UNIFORM_FLOW, AQUIFER, file_list, solver_method)
    val = 'system ' +solver_method + ' '+ str(time.clock() - solve_time) + '\n'
    f.write(str(val))
    
    if do_plot == 1:
        plot_time = time.clock()
        plotting(CRACK_L, LINE_DOUBLE,  UNIFORM_FLOW, AQUIFER, 'compare head', 'Fratura', 'Line_Doublet',solver_method, 'yes')
        plotting(CRACK_L, LINE_DOUBLE,  UNIFORM_FLOW,AQUIFER,'head',  'Fratura','nan', solver_method,'yes')
        plotting(LINE_DOUBLE, LINE_DOUBLE,  UNIFORM_FLOW,AQUIFER,'head', 'Line_Doublet', 'nan',solver_method, 'yes')
        val = 'plot ' +solver_method + ' '+ str(time.clock() - plot_time) + '\n' 
        f.write(str(val))
    
def compare_circle_exact(UNIFORM_FLOW,AQUIFER,solver_method, do_plot):
    
    CIRCLE_L = []
    EXACT_L = []
    f = open('timing.dat','a') 
    solve_time = time.clock()
    solve_circle(CIRCLE_L, UNIFORM_FLOW, AQUIFER, solver_method)
    solve_circle_exact(CIRCLE_L, UNIFORM_FLOW, AQUIFER, EXACT_L)
    val = 'system ' +solver_method + ' '+ str(time.clock() - solve_time) + '\n'
    f.write(str(val))
    
    if do_plot == 1:
        plot_time = time.clock()
        plotting(CIRCLE_L, EXACT_L,  UNIFORM_FLOW, AQUIFER,'compare head', 'Circulo', 'Exato',solver_method, 'yes')
        plotting(CIRCLE_L, EXACT_L,  UNIFORM_FLOW,AQUIFER,'head',  'Circulo','nan', solver_method,'yes')
        plotting(EXACT_L, EXACT_L,  UNIFORM_FLOW,AQUIFER,'head', 'Exato', 'nan',solver_method, 'yes')
        val = 'plot ' +solver_method + ' '+ str(time.clock() - plot_time) + '\n' 
        f.write(str(val))
    
def generate_circle(UNIFORM_FLOW,AQUIFER,solver_method, do_plot):
    
    CIRCLE_L = []
    EXACT_L = []
    f = open('timing.dat','a') 
    solve_time = time.clock()
    solve_circle(CIRCLE_L, UNIFORM_FLOW, AQUIFER, solver_method)
    val = 'system ' +solver_method + ' '+ str(time.clock() - solve_time) + '\n'
    f.write(str(val))
    if do_plot == 1:
        plot_time = time.clock()
        plotting(CIRCLE_L, EXACT_L,  UNIFORM_FLOW, AQUIFER,'head',  'Circulo','nan', solver_method,'yes')
        val = 'plot ' +solver_method + ' '+ str(time.clock() - plot_time) + '\n' 
        f.write(str(val))
    
def compare_circle_exact_ld(UNIFORM_FLOW,AQUIFER, ssubdivisions, solver_method, do_plot):
    
    CIRCLE_L = []
    EXACT_L = []
    LINE_DOUBLE = []
    
    f = open('timing.dat','a') 
    solve_time = time.clock()
    solve_circle(CIRCLE_L, UNIFORM_FLOW, AQUIFER, solver_method)
    solve_circle_exact(CIRCLE_L, UNIFORM_FLOW, AQUIFER, EXACT_L)
    solve_circle_ld(CIRCLE_L, LINE_DOUBLE, UNIFORM_FLOW, AQUIFER, ssubdivisions, solver_method)
    val = 'system ' +solver_method + ' '+ str(time.clock() - solve_time) + '\n'
    f.write(str(val))
    
    if do_plot == 1:
        plot_time = time.clock() 
        plotting(LINE_DOUBLE, EXACT_L,  UNIFORM_FLOW,AQUIFER, 'compare head', 'Line_Doublet', 'Exato', solver_method, 'yes')
        plotting(LINE_DOUBLE, EXACT_L,  UNIFORM_FLOW,AQUIFER,'head',  'Line_Doublet','nan', solver_method,'yes')
        plotting(EXACT_L, EXACT_L,  UNIFORM_FLOW,AQUIFER,'head', 'Exato', 'nan', solver_method, 'yes')
        val = 'plot ' +solver_method + ' '+ str(time.clock() - plot_time) + '\n' 
        f.write(str(val))
    
def generate_ld(UNIFORM_FLOW, AQUFER, solver_method, do_plot):
    LINE_DOUBLE = []
    input_file = 'input_ld.dat'
    file_list = ['ld.dat','ld2.dat']
    #file_list = ['ld.dat']
    f = open('timing.dat','a') 
    solve_time = time.clock()
    
    solve_ld(input_file, LINE_DOUBLE, UNIFORM_FLOW, AQUIFER, file_list, solver_method)
    val = 'system ' +solver_method + ' '+ str(time.clock() - solve_time) + '\n'
    f.write(str(val))
    
    if do_plot == 1:
        plot_time = time.clock() 
        plotting(LINE_DOUBLE, EXACT_L,  UNIFORM_FLOW,AQUIFER,'head',  'Line_Doublet','nan', solver_method,'yes')
        val = 'plot ' +solver_method + ' '+ str(time.clock() - plot_time) + '\n' 
        f.write(str(val))
    
def test_head(UNIFORM_FLOW, AQUIFER, solver_method, do_plot):
    LINE_DOUBLE_1 = []
    LINE_DOUBLE_2 = []
    LINE_DOUBLE_3 = []
    LINE_DOUBLE_4 = []
    
    file_list = ['ld.dat', 'ld2.dat']
    f = open('timing.dat','w') 
    solve_time = time.clock()
    
    input_file = 'input_ld1.dat'
    solve_ld(input_file, LINE_DOUBLE_1, UNIFORM_FLOW, AQUIFER, file_list, solver_method)
    
    input_file = 'input_ld2.dat'
    solve_ld(input_file, LINE_DOUBLE_2, UNIFORM_FLOW, AQUIFER, file_list, solver_method)
    
    input_file = 'input_ld3.dat'
    solve_ld(input_file, LINE_DOUBLE_3, UNIFORM_FLOW, AQUIFER, file_list, solver_method)
    
    input_file = 'input_ld4.dat'
    solve_ld(input_file, LINE_DOUBLE_4, UNIFORM_FLOW, AQUIFER, file_list, solver_method)
    val = 'system ' +solver_method + ' '+ str(time.clock() - solve_time) + '\n'
    f.write(str(val))
    
    LIST_LD = [LINE_DOUBLE_1, LINE_DOUBLE_2,LINE_DOUBLE_3, LINE_DOUBLE_4]
    
    coords = [complex(15,15), complex(85,85)]
    plot_time = time.clock() 
    plotting(LIST_LD, coords,  UNIFORM_FLOW,AQUIFER,'check head', 'Line_Doublet','nan', solver_method,'yes')
    val = 'plot ' +solver_method + ' '+ str(time.clock() - plot_time) + '\n' 
    f.write(str(val))
    
def generate_circle_ld(UNIFORM_FLOW,AQUIFER,ssubdivisions, solver_method, do_plot):
    
    f = open('timing.dat','a') 
    
    CIRCLE_L =[]
    LINE_DOUBLE = []
    EXACT_L = []
    solve_time = time.clock()
    
    solve_circle(CIRCLE_L, UNIFORM_FLOW, AQUIFER, solver_method)
    solve_circle_ld(CIRCLE_L, LINE_DOUBLE, UNIFORM_FLOW, AQUIFER, ssubdivisions, solver_method)
    solve_circle_exact(CIRCLE_L, UNIFORM_FLOW, AQUIFER, EXACT_L)
    
    val = 'system ' +solver_method + ' '+ str(time.clock() - solve_time) + '\n'
    f.write(str(val))
    
    if do_plot == 1:
        plot_time = time.clock()    
        plotting(LINE_DOUBLE, EXACT_L,  UNIFORM_FLOW, AQUIFER,'compare head', 'Line_Doublet', 'Exato', solver_method, 'yes')
        plotting(LINE_DOUBLE, EXACT_L,  UNIFORM_FLOW,AQUIFER,'head',  'Line_Doublet','nan', solver_method,'yes')
    
        val = 'plot ' +solver_method + ' '+ str(time.clock() - plot_time) + '\n' 
        f.write(str(val))
    

def compare_circle_iter_mat(UNIFORM_FLOW,AQUIFER,ssubdivisions, solver_method, do_plot):
    CIRCLE_L =[]
    LINE_DOUBLE = []
    LINE_DOUBLE_2 = []
    EXACT_L = []
    
    f = open('timing.dat','a') 
    solve_time = time.clock()
    solve_circle(CIRCLE_L, UNIFORM_FLOW, AQUIFER, solver_method)
    solve_circle_ld(CIRCLE_L, LINE_DOUBLE, UNIFORM_FLOW, AQUIFER, ssubdivisions, 'matrix')
    solve_circle_ld(CIRCLE_L, LINE_DOUBLE_2, UNIFORM_FLOW, AQUIFER, ssubdivisions, 'iterative')
    solve_circle_exact(CIRCLE_L, UNIFORM_FLOW, AQUIFER, EXACT_L)
    val = 'system ' +solver_method + ' '+ str(time.clock() - solve_time) + '\n'
    f.write(str(val))
    
    if do_plot == 1:
        plot_time = time.clock() 
        plotting(LINE_DOUBLE, LINE_DOUBLE_2,  UNIFORM_FLOW, AQUIFER,'compare head', 'Line_Doublet_Matrix', 'Line_Doublet_Iterative', solver_method, 'yes')
    #plotting(LINE_DOUBLE, EXACT_L,  UNIFORM_FLOW,'head',  'Line Doublet','nan', solver_method,'yes')
        val = 'plot ' +solver_method + ' '+ str(time.clock() - plot_time) + '\n' 
        f.write(str(val))
    
    
def generate_frac(UNIFORM_FLOW, AQUIFER, solver_method, do_plot):
    CRACK_L = []
    LINE_DOUBLE = []
    
    f = open('timing.dat','a') 
    solve_time = time.clock()
    solve_crack( CRACK_L, UNIFORM_FLOW, AQUIFER, solver_method)
    solve_ld_crack(CRACK_L, LINE_DOUBLE,UNIFORM_FLOW, AQUIFER, solver_method)
    val = 'system ' +solver_method + ' '+ str(time.clock() - solve_time) + '\n'
    f.write(str(val))
    
    if do_plot == 1:
        plot_time = time.clock() 
        plotting(CRACK_L, LINE_DOUBLE,  UNIFORM_FLOW,AQUIFER, 'compare head', 'Fratura', 'Line_Doublet', solver_method, 'yes')
        plotting(CRACK_L, LINE_DOUBLE,  UNIFORM_FLOW,AQUIFER,'head',  'Fratura','nan', solver_method,'yes')
        val = 'plot ' +solver_method + ' ' + str(time.clock() - plot_time) + '\n' 
        f.write(str(val))
    
if __name__ == '__main__':

    solver_method = sys.argv[1]
    system = sys.argv[2]
    subdivisions = 20
    do_plot = 1
    
    CRACK_L = []
    CRACK_L2 = []
    LINE_DOUBLE = []
    LINE_DOUBLE2 = []
    CIRCLE_L = []
    CIRCLE_L2 = []
    EXACT_L = []
    EXACT_L2 = []
    ELLIPSE_L = []
    UNIFORM_FLOW, AQUIFER = set_u_flow_aq()
    
    if system == 'circle':
        generate_circle(UNIFORM_FLOW,AQUIFER,solver_method, do_plot)
    elif system == 'compare_circle_exact':
        compare_circle_exact(UNIFORM_FLOW,AQUIFER,solver_method, do_plot)
    elif system == 'compare_circle_exact_ld':
        compare_circle_exact_ld(UNIFORM_FLOW,AQUIFER, subdivisions, solver_method, do_plot)
    elif system == 'iter_mat_compare':
        compare_circle_iter_mat(UNIFORM_FLOW, AQUIFER, subdivisions, solver_method, do_plot)
    elif system == 'crack':
        generate_frac(UNIFORM_FLOW, AQUIFER, solver_method, do_plot)
    elif system == 'only_crack':
        crack_alone(UNIFORM_FLOW, AQUIFER, solver_method, solver_method, do_plot)
    elif system == 'ld':
        generate_ld(UNIFORM_FLOW,AQUIFER,solver_method, do_plot)
    elif system == 'circle_ld':
        generate_circle_ld(UNIFORM_FLOW, AQUIFER, subdivisions, solver_method, do_plot)
    elif system == 'both_crack_ld':
        crack_with_ld(UNIFORM_FLOW, AQUIFER, solver_method, solver_method, do_plot)
    else:
        print "wrong system"
        sys.exit(-1)
    #
    #
    #
    #test_head(UNIFORM_FLOW,AQUIFER,solver_method)
    #
    #
    #generate_frac(UNIFORM_FLOW, AQUIFER, solver_method)
    
    #solve_ellipse(ELLIPSE_L,UNIFORM_FLOW,AQUIFER,solver_method)
    #solve_crack(CRACK_L, UNIFORM_FLOW, AQUIFER, solver_method)
    #solve_ld(CRACK_L, LINE_DOUBLE,UNIFORM_FLOW, AQUIFER, solver_method)
    #solve_ld_check(LINE_DOUBLE, UNIFORM_FLOW, AQUIFER,100, solver_method)
    #solve_ld_check(LINE_DOUBLE2, UNIFORM_FLOW,AQUIFER,10,'iterative')
    #solve_dual_ld_check(LINE_DOUBLE, UNIFORM_FLOW,AQUIFER,10,'iterative')
    #solve_circle(CIRCLE_L, UNIFORM_FLOW, AQUIFER, solver_method)
    #solve_circle(CIRCLE_L2, UNIFORM_FLOW, AQUIFER, 'iterative')
    #solve_circle_exact(CIRCLE_L, UNIFORM_FLOW, AQUIFER, EXACT_L,solver_method)
    #solve_dual_circle(CIRCLE_L, UNIFORM_FLOW, AQUIFER, 'matrix')
    
    #solve_circle(CIRCLE_L, UNIFORM_FLOW, AQUIFER, 'iterative')
    #solve_circle_ld(CIRCLE_L, LINE_DOUBLE, UNIFORM_FLOW, AQUIFER, subdivisions, solver_method)
    #solve_circle_exact_ld(CIRCLE_L, LINE_DOUBLE, UNIFORM_FLOW, AQUIFER, EXACT_L, subdivisions, solver_method)
    #solve_circle_exact_ld(CIRCLE_L, LINE_DOUBLE2, UNIFORM_FLOW, AQUIFER, EXACT_L2, subdivisions, 'iterative')
    #solve_circle_ld(CIRCLE_L, LINE_DOUBLE, UNIFORM_FLOW, AQUIFER, EXACT_L, subdivisions, solver_method)
    
    
    #plotting(CIRCLE_L, CIRCLE_L2 , UNIFORM_FLOW, 'compare head')
    #plotting(LINE_DOUBLE, EXACT_L2 , UNIFORM_FLOW, 'compare head')
    #plotting(LINE_DOUBLE2, EXACT_L2 , UNIFORM_FLOW, 'compare head')
    #plotting(LINE_DOUBLE, CIRCLE_L , UNIFORM_FLOW, 'compare head')
    #plotting(LINE_DOUBLE, EXACT_L,  UNIFORM_FLOW, 'compare head')
    #plotting(LINE_DOUBLE2, EXACT_L2,  UNIFORM_FLOW, 'compare head')
    
    #plotting(CRACK_L, LINE_DOUBLE,  UNIFORM_FLOW, 'compare head')
    #plotting(ELLIPSE_L, CIRCLE_L,  UNIFORM_FLOW, 'head')
    #plotting(LINE_DOUBLE, CIRCLE_L,  UNIFORM_FLOW, 'head')
    
    #x1, y1, x2, y2, step = aux_functions.set_axis(CIRCLE_L)
    #err_analytics = ErrorAnalysis.ErrorAnalysis(x1,y1,x2,y2,step,CIRCLE_L,EXACT_L,UNIFORM_FLOW,'head')
    #x1, y1, x2, y2, step = aux_functions.set_axis(CIRCLE_L)
    #err_analytics = ErrorAnalysis.ErrorAnalysis(x1,y1,x2,y2,step,CIRCLE_L,CIRCLE_L2,UNIFORM_FLOW,'head')
    #ErrorAnalysis.check_circle_jump(CIRCLE_L, UNIFORM_FLOW, AQUIFER, 1e-7)
    #ErrorAnalysis.check_circle_jump(CIRCLE_L2, UNIFORM_FLOW, AQUIFER, 1e-7)
    #ErrorAnalysis.check_ld_jump(LINE_DOUBLE, UNIFORM_FLOW, AQUIFER, 1e-7)
    #ErrorAnalysis.check_circle_jump(EXACT_L, UNIFORM_FLOW, AQUIFER, 1e-7)
    
    #ErrorAnalysis.compare_results_absolute(err_analytics.a, err_analytics.b)
    #ErrorAnalysis.compare_results_relative(err_analytics.a, err_analytics.b)

    #command = """run_it()"""
    #cProfile.runctx( command, globals(), locals(), filename="pyaem.profile" )
    
    