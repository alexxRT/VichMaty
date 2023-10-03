import numpy as np
import pandas as pd
import numpy.linalg as linal
import scanf as sc


#here and everywhere next:
# A - matrix of linear sysytem
# f - right part of the system
# typeA = typef = numpy array

# default value, can be specifed on start
default_precision = 5

class MODE:
    PRECISE = 1
    APPROX  = 2


def round_input(A, f, precision):
    A = np.round(A, precision)
    f = np.round(f, precision)

    return A, f

def valid_input(A, f):
    A_lines = len(A)
    A_columns   = len(A[0])

    for i in range (0, A_lines):
        if (len(A[i]) != A_columns):
            return False # A is not a matrix at all
        
    f_lines   = len(f)
    f_columns = np.array(f[0]).size

    if (A_columns != A_lines): # A is not square matrix
        return False

    if (f_columns != 1): # f is not a vector
        return False
    
    if (f_lines != A_lines): # f and A sizes are not coordinated 
        return False
    
    return True



def estimate_precision(A, f, precision):
    mu = 0
    A_inverse = 0

    try:
        A_inverse = linal.inv(A)
    except linal.LinAlgError:
        print ("Matrix A does not have inverse; mu = 1")
        mu = 1

    norm_A = linal.norm(A, 2)
    norm_f = linal.norm(f, 2)
    norm_A_inverse = 0

    if (mu == 0):
        norm_A_inverse = linal.norm(A_inverse, 2)
        mu = norm_A_inverse * norm_A


    delta_A = A - np.round(A, precision)
    delta_f = f - np.round(f, precision)

    norm_delta_A = linal.norm(delta_A, 2)
    norm_delta_f = linal.norm(delta_f, 2)

    resp_delta_A = norm_delta_A / norm_A
    resp_delta_f = norm_delta_f / norm_f


    error_percent  = (mu * (1 / (1 + mu * resp_delta_A)) * (resp_delta_f + resp_delta_A)) * 100

    return error_percent

def my_func(A, f):
    return 0


def calculate_linear (A, f, mode, method, input_precision = default_precision):
    if (not valid_input(A, f)):
        return False
    
    sol_prec = 0
    if (mode == MODE.APPROX):
        sol_prec = estimate_precision(A, f, input_precision)

        print ("Solution error estimation: ", "{:.2f}".format(sol_prec), "%")
        print ("Do u whant to improve/agree error or calculate solution with no round?\n \
               [1] - improve [2] - agree, [3] - no round")
        
        input = 0
        input = sc.scanf("%d")[0]

        while (input != 2 and input != 3):
            print("New precision:", end=" ", flush=True)
            
            input_precision = sc.scanf("%d")[0]
            sol_prec = estimate_precision(A, f, input_precision)

            print ("Solution error estimation: ", "{:.2f}".format(sol_prec), "%")

            print("[1] - improve, [2] - agree [3] - no round")
            input = sc.scanf("%d")[0]

        if (input == 3):
            return method(A, f)
        
        A_round = np.round(A, input_precision)
        f_round = np.round(f, input_precision)

        return method(A_round, f_round) 


    return method(A, f)


A = np.array([[100, 99], [99, 98]])
f = np.array([198.99, 197.01])

calculate_linear(A, f, MODE.APPROX, my_func, 0)