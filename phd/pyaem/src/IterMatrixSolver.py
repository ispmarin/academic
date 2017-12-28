'''
Created on 03/06/2011

@author: ispmarin
'''
import numpy as np
import scipy.linalg

import convergence

def MatrixIteration(Pinv, N, b, dim_A, max_iterations,conv_test):
    
    iteration = 0 
    x_k = np.zeros(dim_A)
    result_vector = np.zeros(dim_A)
    
    while iteration < max_iterations:
        result_vector = np.dot(np.dot(Pinv,N),x_k) + np.dot(Pinv,b) 
        iteration = iteration + 1
            
        if conv_test.check_convergence_vector(result_vector) == 1:
            x_k = result_vector
            break
        x_k = result_vector
            
    return result_vector, iteration
    
def EigenTest(Pinv, N,solver):
    f = open('eigen.dat', 'a')
    eigen = abs(np.max(scipy.linalg.eigvalsh(np.dot(Pinv,N))))
    if eigen >= 1.0:
        print solver, 'Eigenvalue larger than 1', eigen
    print solver, 'Eigenvalue ', eigen
    f.write(str(str(solver) + '   ' + str(eigen)))

def Jacobi(A,b, abs_conv_tol, rel_conv_tol, max_iterations):
    dim_A = len(A)
    M = np.zeros((dim_A, dim_A))
    N = np.zeros((dim_A, dim_A))
        
    for i in xrange(dim_A):
        for j in xrange(dim_A):
            if i == j:
                M[i,j] = A[i,j]
            else:
                N[i,j] = -A[i,j]
                
    
    Pinv = scipy.linalg.inv(M)
    EigenTest(Pinv,N,'Jacobi')
    conv_test = convergence.convergence_vector(dim_A, abs_conv_tol, rel_conv_tol, np.zeros(dim_A))
    
    return MatrixIteration(Pinv, N, b, dim_A, max_iterations,conv_test)   
    
        
def GaussSeidel(A, b, abs_conv_tol,rel_conv_tol, max_iterations):
    dim_A = len(A)
    M = np.zeros((dim_A, dim_A))
    N = np.zeros((dim_A, dim_A))
    
    for i in xrange(dim_A):
        for j in xrange(dim_A):
            if i < j:
                N[i,j] = -A[i,j]
            elif i > j:
                M[i,j] = A[i,j]
            else:
                M[i,j] = A[i,j]
    
    Pinv = scipy.linalg.inv(M)
    EigenTest(Pinv,N,'Gauss Seidel') 
    conv_test = convergence.convergence_vector(dim_A, abs_conv_tol, rel_conv_tol, np.zeros(dim_A))
    
    return MatrixIteration(Pinv, N, b, dim_A, max_iterations,conv_test)
