import numpy as np
import pandas as pd
import numpy.linalg as linal
import sys

max_float_precision = 0.00001

def is_zero(num):
    return (abs(num) < max_float_precision)

def is_LU_possible(A):
    for i in range (0, len(A)):
        leading_minor_A = linal.det(A[0:i + 1, 0:i + 1])

        if (is_zero(leading_minor_A)):
            return False

    return True

def LU_decomposition(A):
    # in case if decomposition is not possible
    if (not is_LU_possible(A)):
        raise -1
    
    A_xshape = A.shape[0]
    A_yshape = A.shape[1]
    
    U = np.zeros((A_xshape, A_yshape))
    L = np.zeros((A_xshape, A_yshape))

    for n in range (0, A_xshape):

        #finding U matrix elements for n's row
        for i in range (n, A_xshape):
            
            #counting (Lnk)*U(ki)
            sum = 0
            for k in range (0, n):
                sum += L[n][k] * U [k][i]

            U[n][i] = A[n][i] - sum

        #finfing L matrix elements for n's column
        for i in range (n, A_xshape):
            
            sum = 0
            for k in range (0, n):
                sum += L[i][k] * U[k][n]

            L[i][n] = (A[i][n] - sum) / U[n][n]

    return L, U


def LU_method(A, f):
    # prepare solution vector
    A_shape = A.shape
    n = A_shape[0]
    solution = np.zeros(n)

    L = 0
    U = 0

    try:
        L, U = LU_decomposition(A)
    except:
        print("LU decomposition can not be performed")
    else:
        y = np.zeros(n)
        #Solving for Ly = f
        for i in range (0, len(y)):

            sum = 0
            for j in range (0, i):
                sum += y[j]*L[i][j]
        
            y[i] = (f[i] - sum) / L[i][i]

        #Solving for Ux = y
        for i in range (len(solution) - 1, -1, -1):
            sum = 0
            for j in range (len(solution) - 1, i, -1):
                sum += solution[j]*U[i][j]

            solution[i] = (y[i] - sum) / U[i][i]

        print (solution)
    return solution


def Gaus_method(A, f):
    # prepare solution vector
    A_shape = A.shape
    n = A_shape[0]
    solution = np.zeros(n)

    # applying gauss elimination
    for j in range(0, n):

        # go for elements to be eliminated  
        for i in range(j + 1, n):
            ratio = A[i][j]/A[j][j]
            
            # substruct line
            for k in range(0, n):
                A[i][k] = A[i][k] - ratio * A[j][k]
            # correct right part vector
            f[i] = f[i] - ratio * f[j]

    # back substitution
    for i in range(n - 1, -1, -1):
        sum = 0
        for j in range(n - 1, i, -1):
            sum += A[i][j]*solution[j]

        solution[i] = (f[i] - sum) / A[i][i]
    
    print(solution)
    
    return solution


A = np.array([[1, 2], [1, 1]])
f = np.array([12, 17])

Gaus_method(A, f)
LU_method(A, f)